"""Static-HTTP mirroring: session reuse, resume, checksum and content validation.

Everything here talks to static file servers only — plain files with
byte-range support — never to a dynamic cutout renderer. Each worker
thread keeps one `requests.Session` with keep-alive, because at ~100 ms
per TLS handshake a fetch of this size without connection reuse spends
more time handshaking than transferring. Workers are threads, not
processes: the workload is network wait.

Every response is validated on content, not status. The Legacy Survey
estate returns HTTP 200 with short text bodies on error, so a download
only counts when the file starts with the FITS magic and astropy can
open it — or, where the server publishes a per-directory `.sha256sum`
manifest, when the digest matches. Partial downloads resume via HTTP
Range from a `.part` file and are renamed into place only after
validation, so an existing file is always a valid one.
"""

import hashlib
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests
from astropy.io import fits
from tqdm import tqdm

from . import urls

USER_AGENT = "SGA-training-data/1.0 (SteelRail/SGA; GRAF dataset build)"
TIMEOUT = (15, 300)  # connect, read (seconds)
CHUNK = 1 << 20
RETRIES = 5

_local = threading.local()


def get_session():
    """One keep-alive session per thread."""
    session = getattr(_local, "session", None)
    if session is None:
        session = requests.Session()
        session.headers["User-Agent"] = USER_AGENT
        _local.session = session
    return session


def valid_fits(path):
    """True when `path` starts with the FITS (or gzip) magic and every HDU header parses."""
    try:
        with open(path, "rb") as f:
            head = f.read(9)
        if not (head.startswith(b"SIMPLE  =") or head.startswith(b"\x1f\x8b")):
            return False
        with fits.open(path, memmap=True) as hdul:
            for hdu in hdul:
                _ = hdu.header  # forces a seek through the whole HDU list
        return True
    except Exception:
        return False


def sha256_of(path):
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(CHUNK):
            digest.update(chunk)
    return digest.hexdigest()


def fetch_checksums(url, session=None):
    """Parse a published `.sha256sum` manifest into {filename: hexdigest}.

    Returns an empty dict when the manifest does not exist — validation
    then falls back to the FITS content check.
    """
    session = session or get_session()
    try:
        response = session.get(url, timeout=TIMEOUT)
    except requests.RequestException:
        return {}
    if response.status_code != 200 or response.content[:1] == b"<":
        return {}
    checksums = {}
    for line in response.text.splitlines():
        parts = line.split()
        if len(parts) == 2 and len(parts[0]) == 64:
            checksums[parts[1].lstrip("*")] = parts[0]
    return checksums


def _validate(path, checksum):
    if checksum is not None:
        return sha256_of(path) == checksum
    return valid_fits(path)


def fetch(url, dest, checksum=None, session=None):
    """Mirror one file; returns the number of bytes actually transferred.

    Skips work when `dest` already exists (existing files are valid by
    construction — see module docstring). Resumes an interrupted transfer
    from `dest.part` via HTTP Range, retries transient failures with
    exponential backoff, and validates before renaming into place.
    """
    dest = Path(dest)
    if dest.exists():
        return 0
    session = session or get_session()
    dest.parent.mkdir(exist_ok=True, parents=True)
    part = dest.with_suffix(dest.suffix + ".part")

    last_error = None
    for attempt in range(RETRIES):
        if attempt:
            time.sleep(min(2 ** attempt, 60))
        try:
            transferred = _download(session, url, part)
        except requests.RequestException as error:
            last_error = error
            continue
        if _validate(part, checksum):
            part.rename(dest)
            return transferred
        last_error = ValueError(f"validation failed: {url}")
        part.unlink(missing_ok=True)  # do not resume from corrupt bytes
    raise RuntimeError(f"gave up after {RETRIES} attempts: {url}") from last_error


def _download(session, url, part):
    """Stream `url` to `part`, resuming from its current size when possible."""
    headers = {}
    offset = part.stat().st_size if part.exists() else 0
    if offset:
        headers["Range"] = f"bytes={offset}-"
    with session.get(url, headers=headers, stream=True, timeout=TIMEOUT) as response:
        if response.status_code == 416:  # already fully transferred
            return 0
        response.raise_for_status()
        mode = "ab" if offset and response.status_code == 206 else "wb"
        transferred = 0
        with open(part, mode) as f:
            for chunk in response.iter_content(CHUNK):
                f.write(chunk)
                transferred += len(chunk)
    return transferred


def fetch_many(tasks, workers=12, desc="fetching"):
    """Mirror many (url, dest[, checksum]) tasks concurrently.

    Returns (bytes_transferred, failures) where failures is a list of
    (url, exception). Concurrency is bounded by `workers` threads, each
    with its own keep-alive session.
    """
    tasks = [task if len(task) == 3 else (*task, None) for task in tasks]
    transferred = 0
    failures = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(fetch, url, dest, checksum): url
            for url, dest, checksum in tasks
        }
        with tqdm(total=len(futures), desc=desc, unit="file") as progress:
            for future in as_completed(futures):
                try:
                    transferred += future.result()
                except Exception as error:
                    failures.append((futures[future], error))
                progress.update()
                progress.set_postfix_str(f"{transferred / 1e9:.2f} GB")
    return transferred, failures


def brick_dir(root, brickname):
    """Mirror directory of one brick's files under the data root."""
    return Path(root) / "bricks" / brickname[:3] / brickname


def mirror_bricks(bricks, root, workers=12, desc="bricks"):
    """Mirror every coadd file of the given (brickname, hemisphere) pairs.

    Bricks already fully mirrored are skipped without touching the
    network; the rest are verified against the published sha256
    manifests. Returns (bytes_transferred, bricknames_with_a_failed_file)
    so callers decide their own failure policy.
    """
    pending = []
    for brickname, hemisphere in bricks:
        directory = brick_dir(root, brickname)
        files = urls.brick_files(brickname, hemisphere)
        if not all((directory / name).exists() for name in files):
            pending.append((brickname, hemisphere, directory, files))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        manifests = list(pool.map(
            lambda p: fetch_checksums(urls.brick_checksums(p[0], p[1])),
            pending,
        ))
    tasks, owner = [], {}
    for (brickname, _, directory, files), checksums in zip(pending, manifests,
                                                           strict=True):
        for name, url in files.items():
            tasks.append((url, directory / name, checksums.get(name)))
            owner[url] = brickname
    transferred, failures = fetch_many(tasks, workers=workers, desc=desc)
    return transferred, {owner[url] for url, _ in failures}
