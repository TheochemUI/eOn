"""Store path-based con blobs in a readcon-db corpus.

readcon still parses every frame. The corpus is a second copy: one
trajectory per distinct path plus blob, next to the run when ``config.ini``
is at most eight directories up, otherwise beside the con file. A missing
``readcon_db`` install leaves the con file as the only copy.
"""

from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger("eon.concorpus")

_FNV_OFFSET = 14695981039346656037
_FNV_PRIME = 1099511628211
_MASK = (1 << 64) - 1

_unavailable = False
_corpora: dict[str, object] = {}


def traj_id(path: str, text: str) -> int:
    """FNV-1a 64 over the absolute path, a NUL, and the con text."""
    h = _FNV_OFFSET
    for byte in (path + "\0" + text).encode():
        h ^= byte
        h = (h * _FNV_PRIME) & _MASK
    return h or 1


def corpus_dir(con_path: Path) -> Path:
    """Directory of the LMDB corpus for this con path."""
    cur = con_path.resolve().parent
    for _ in range(8):
        if (cur / "config.ini").is_file():
            return cur / "readcon.db"
        parent = cur.parent
        if parent == cur:
            break
        cur = parent
    return con_path.resolve().parent / "readcon.db"


def _corpus(directory: Path):
    global _unavailable
    if _unavailable:
        return None
    key = str(directory)
    cached = _corpora.get(key)
    if cached is not None:
        return cached
    try:
        from readcon_db import ConCorpus
    except ImportError:
        _unavailable = True
        logger.debug("readcon_db is not installed; con files stay authoritative")
        return None
    try:
        db = ConCorpus(key)
    except Exception:
        logger.debug("readcon-db open failed for %s", key, exc_info=True)
        return None
    _corpora[key] = db
    return db


def mirror_con_text(path, text: str) -> None:
    """Insert one con blob. Identical path-and-text pairs are stored once."""
    if not text or not text.strip():
        return
    con_path = Path(path)
    try:
        resolved = str(con_path.resolve())
    except OSError:
        resolved = str(con_path)
    db = _corpus(corpus_dir(con_path))
    if db is None:
        return
    try:
        db.append_trajectory_str(traj_id(resolved, text), text, source=resolved)
    except Exception as exc:
        if "already exists" in str(exc):
            return
        logger.debug("readcon-db mirror skipped for %s: %s", resolved, exc)


def stored_frame_text(path) -> str | None:
    """Frame text for this path's current blob, or None when the corpus has no copy."""
    con_path = Path(path)
    try:
        text = con_path.read_text()
        resolved = str(con_path.resolve())
    except OSError:
        return None
    db = _corpus(corpus_dir(con_path))
    if db is None:
        return None
    try:
        return db.get_frame_text(traj_id(resolved, text), 0)
    except Exception:
        return None


def mirror_con_path(path) -> None:
    """Insert the con file's current bytes."""
    con_path = Path(path)
    try:
        text = con_path.read_text()
    except OSError:
        return
    mirror_con_text(con_path, text)
