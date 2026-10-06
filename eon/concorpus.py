"""Store path-based con blobs in a readcon-db corpus.

readcon still parses every frame. The corpus is a second copy: one
trajectory per distinct path plus blob, next to the run when ``config.ini``
is at most eight directories up, otherwise beside the con file. A missing
``readcon_db`` install leaves the con file as the only copy.

A dynamics movie can also ingest each directory of a simulation tree and
read those frames back with one read-only ``get_frame_texts`` call. The
``.con`` files remain the store.
"""

from __future__ import annotations

import gc
import logging
import os
import tempfile
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


def _loader_error() -> str:
    """dlopen text for a missing libreadcon_db.so, or an empty string."""
    import ctypes

    try:
        ctypes.CDLL("libreadcon_db.so")
    except OSError as exc:
        return str(exc)
    return ""


def _note_library_missing(exc: BaseException) -> None:
    """One warning. Later calls return before this runs."""
    detail = _loader_error() or str(exc)
    logger.warning("libreadcon_db.so failed to load: %s", detail)


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
    except (ImportError, OSError) as exc:
        _unavailable = True
        _note_library_missing(exc)
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


def store_frame_text(path, text: str) -> tuple[int, int] | None:
    """Store one con blob. Returns ``(traj_id, 0)`` for a single frame.

    The corpus takes the con text only. A barrier is not an argument.
    A mode is not an argument either: it rides in the saddle frame as
    the displacements section.
    """
    if not text or not text.strip():
        return None
    con_path = Path(path)
    try:
        resolved = str(con_path.resolve())
    except OSError:
        resolved = str(con_path)
    db = _corpus(corpus_dir(con_path))
    if db is None:
        return None
    key = traj_id(resolved, text)
    try:
        db.append_trajectory_str(key, text, source=resolved)
    except Exception as exc:
        if "already exists" not in str(exc):
            logger.warning("readcon-db store failed for %s: %s", resolved, exc)
            return None
    return key, 0


def load_frame_text(directory: Path, traj: int, frame_idx: int) -> str | None:
    """Frame text for a key in this corpus, or None when it is absent."""
    db = _corpus(Path(directory))
    if db is None:
        return None
    try:
        return db.get_frame_text(int(traj), int(frame_idx))
    except Exception:
        logger.warning(
            "readcon-db get_frame_text failed for %s frame %s:%s",
            directory,
            traj,
            frame_idx,
            exc_info=True,
        )
        return None


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


def directories_with_con(tree: Path) -> list[Path]:
    """Directories under ``tree`` that directly contain ``.con`` or ``.convel``."""
    found: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(tree):
        dirnames.sort()
        if any(name.endswith(".con") or name.endswith(".convel") for name in filenames):
            found.append(Path(dirpath))
    return found


def _index_tree(db, tree: Path) -> dict[str, tuple[int, int]]:
    index: dict[str, tuple[int, int]] = {}
    next_id = 1
    for directory in directories_with_con(tree):
        rows = db.ingest_directory(str(directory), start_traj_id=next_id)
        for traj_id, _nframes, path in rows:
            index[str(Path(path).resolve())] = (int(traj_id), 0)
            next_id = max(next_id, int(traj_id) + 1)
    return index


def frame_texts_for_paths(tree, con_paths, corpus_cls=None) -> list[str]:
    """Ingest ``tree`` and return frame 0 text for each path, in order.

    ``corpus_cls`` defaults to ``readcon_db.ConCorpus``. Tests pass a stand-in
    with the same constructor, ``ingest_directory``, and ``get_frame_texts``.
    """
    if corpus_cls is None:
        from readcon_db import ConCorpus

        corpus_cls = ConCorpus
    tree = Path(tree)
    with tempfile.TemporaryDirectory(prefix="eon-concorpus-") as corpus_dir:
        writer = corpus_cls(corpus_dir)
        index = _index_tree(writer, tree)
        # One LMDB environment per path. Drop the writer before the read-only open.
        close = getattr(writer, "close", None)
        if close is not None:
            close()
        del writer
        gc.collect()
        keys = []
        for path in con_paths:
            resolved = str(Path(path).resolve())
            if resolved not in index:
                raise FileNotFoundError(f"{path} was not ingested from {tree}")
            keys.append(index[resolved])
        reader = corpus_cls(corpus_dir, readonly=True)
        return list(reader.get_frame_texts(keys))
