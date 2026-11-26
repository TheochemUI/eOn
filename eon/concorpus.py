"""Store path-based con blobs in a readcon-db corpus.

readcon still parses every frame. The corpus is a second copy: one
trajectory per distinct path plus blob, next to the run when ``config.ini``
is at most eight directories up, otherwise beside the con file. A missing
``readcon_db`` install leaves the con file as the only copy.

An overwrite of a path starts a generation, a trajectory whose frame 0 is
the written frame, and later appends extend that trajectory. Compressed
con files stay out of the corpus.

A dynamics movie can also ingest each directory of a simulation tree and
read those frames back with one read-only ``get_frame_texts`` call. The
``.con`` files remain the store.
"""

from __future__ import annotations

import atexit
import gc
import logging
import os
import tempfile
import threading
from pathlib import Path

logger = logging.getLogger("eon.concorpus")

_FNV_OFFSET = 14695981039346656037
_FNV_PRIME = 1099511628211
_MASK = (1 << 64) - 1

_GROUP = 64
_unavailable = False
_corpora: dict[str, object] = {}
_pending: dict[str, tuple[int, list[str]]] = {}
_generations: dict[str, int] = {}
_pending_lock = threading.Lock()


def _fnv1a(data: bytes) -> int:
    h = _FNV_OFFSET
    for byte in data:
        h ^= byte
        h = (h * _FNV_PRIME) & _MASK
    return h or 1


def blob_id(path: str, text: str) -> int:
    """FNV-1a 64 over the canonical path, a NUL, and the con text.

    Whole-file blobs are keyed by content. The corpus has no replace, so
    a rewrite of a path stores a new blob under a new id, and a reader
    finds the current one from the file's text. The C++ mirror uses the
    same bytes.
    """
    return _fnv1a((path + "\0" + text).encode())


def generation_id(path: str, first_text: str, generation: int = 0) -> int:
    """Trajectory id of one generation of ``path`` whose frame 0 is ``first_text``.

    Generation 0 is ``blob_id(path, first_text)``, so the blob an overwrite
    stores is the trajectory its appends extend. The corpus cannot drop
    frames, so an overwrite that repeats the first frame of an extended
    generation takes the next one: FNV-1a 64 over the path, a NUL, the
    text, a NUL and the decimal generation. The C++ mirror uses the same
    bytes.
    """
    if generation == 0:
        return blob_id(path, first_text)
    return _fnv1a(f"{path}\0{first_text}\0{generation}".encode())


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


def _resolved(path) -> str:
    con_path = Path(path)
    try:
        return str(con_path.resolve())
    except OSError:
        return str(con_path)


def _extend_batch(resolved: str, tid: int, frames: list[str]) -> None:
    """Extend one batch. A missing symbol or library is a soft skip."""
    if not frames:
        return
    text = "".join(frames)
    db = _corpus(corpus_dir(Path(resolved)))
    if db is None:
        return
    extend_str = getattr(db, "extend_trajectory_str", None)
    try:
        if extend_str is not None:
            extend_str(tid, text, source=resolved)
            return
        extend = getattr(db, "extend_trajectory", None)
        if extend is None:
            logger.debug("readcon-db has no extend entry; corpus batch skipped")
            return
        fd, scratch = tempfile.mkstemp(prefix="eon-corpus-", suffix=".con")
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as fh:
                fh.write(text)
            extend(tid, scratch)
        finally:
            try:
                os.unlink(scratch)
            except OSError:
                pass
    except Exception as exc:
        logger.debug("readcon-db extend skipped for %s: %s", resolved, exc)


def _has_frame(db, tid: int, frame_idx: int) -> bool:
    try:
        db.frame_hash(int(tid), int(frame_idx))
    except Exception:
        return False
    return True


def _stored_frames(db, tid: int, limit: int) -> int:
    """Frames under ``tid``, or ``limit + 1`` when there are more than ``limit``."""
    lo, hi = 0, limit + 1
    while lo < hi:
        mid = (lo + hi) // 2
        if _has_frame(db, tid, mid):
            lo = mid + 1
        else:
            hi = mid
    return lo


def _seed(db, tid: int, text: str, resolved: str) -> bool:
    try:
        db.append_trajectory_str(tid, text, source=resolved)
    except Exception as exc:
        if "already exists" not in str(exc):
            logger.debug("readcon-db seed skipped for %s: %s", resolved, exc)
            return False
    return True


def _flush_path(resolved: str) -> None:
    with _pending_lock:
        slot = _pending.pop(resolved, None)
    if slot is not None:
        _extend_batch(resolved, *slot)


def start_generation(path, text: str) -> int | None:
    """Start the trajectory an overwrite of ``path`` with one frame ``text`` begins.

    Frames still buffered for the previous generation go to it first.
    A generation that holds only this frame is reused. Returns the
    trajectory id, or None when no corpus takes the frame.
    """
    resolved = _resolved(path)
    _flush_path(resolved)
    with _pending_lock:
        _generations.pop(resolved, None)
    if not text or not text.strip():
        return None
    db = _corpus(corpus_dir(Path(resolved)))
    if db is None:
        return None
    generation = 0
    while _has_frame(db, generation_id(resolved, text, generation), 1):
        generation += 1
    tid = generation_id(resolved, text, generation)
    if not _seed(db, tid, text, resolved):
        return None
    with _pending_lock:
        _generations[resolved] = tid
    return tid


def _disk_frame_texts(resolved: str) -> list[str]:
    try:
        import readcon

        frames = readcon.read_con(resolved)
        return [readcon.write_con_string([frame]) for frame in frames]
    except Exception as exc:
        logger.debug("readcon-db resume skipped for %s: %s", resolved, exc)
        return []


def _resume_generation(resolved: str) -> int | None:
    """Bring the generation of an appended file this process did not start up to disk.

    The newest generation for the file's first frame is used unless a
    stored frame disagrees with the file, which starts a new one.
    The just-appended frame is on disk, so it is stored here too.
    """
    db = _corpus(corpus_dir(Path(resolved)))
    if db is None:
        return None
    texts = _disk_frame_texts(resolved)
    if not texts:
        return None
    first = texts[0]
    generation = 0
    while _has_frame(db, generation_id(resolved, first, generation + 1), 0):
        generation += 1
    tid = generation_id(resolved, first, generation)
    stored = _stored_frames(db, tid, len(texts))
    if not _stored_prefix_matches(db, tid, texts, stored):
        tid = generation_id(resolved, first, generation + 1)
        stored = 0
    if stored == 0:
        if not _seed(db, tid, "".join(texts), resolved):
            return None
    else:
        _extend_batch(resolved, tid, texts[stored:])
    with _pending_lock:
        _generations[resolved] = tid
    return tid


def _frame_text(db, tid: int, frame_idx: int) -> str | None:
    try:
        return db.get_frame_text(int(tid), int(frame_idx))
    except Exception:
        return None


_STORED_PRECISION = 6


def _normalized_frame(text: str) -> str | None:
    try:
        import readcon

        return readcon.write_con_string(
            readcon.read_con_string(text), precision=_STORED_PRECISION
        )
    except Exception:
        return None


def _same_frame(stored: str | None, disk: str) -> bool:
    """True when a stored frame and a reread disk frame hold the same values.

    readcon-db 0.1.6 writes extended frames with six decimals, while a
    reread writes round-trip text. Both sides go through that six-decimal
    writer, so a difference below it does not count.
    """
    if stored is None:
        return False
    if stored == disk:
        return True
    normalized = _normalized_frame(stored)
    return normalized is not None and normalized == _normalized_frame(disk)


def _stored_prefix_matches(db, tid: int, texts: list[str], stored: int) -> bool:
    """True when the stored frames are a prefix of ``texts``."""
    if stored > len(texts):
        return False
    for index in range(stored):
        if not _same_frame(_frame_text(db, tid, index), texts[index]):
            return False
    return True


def note_corpus_append(path, frame_text: str) -> None:
    """Buffer one appended frame for the path's generation. Flush at 64 frames."""
    if not frame_text or not frame_text.strip():
        return
    resolved = _resolved(path)
    with _pending_lock:
        tid = _generations.get(resolved)
    if tid is None:
        _resume_generation(resolved)
        return
    ready: tuple[int, list[str]] | None = None
    with _pending_lock:
        slot = _pending.setdefault(resolved, (tid, []))
        slot[1].append(frame_text)
        if len(slot[1]) >= _GROUP:
            ready = slot
            del _pending[resolved]
    if ready is not None:
        _extend_batch(resolved, *ready)


def flush_corpus_appends() -> None:
    """Extend every buffered batch. A crash before this drops at most 63 frames."""
    with _pending_lock:
        pending = list(_pending.items())
        _pending.clear()
    for resolved, (tid, frames) in pending:
        _extend_batch(resolved, tid, frames)


def corpus_buffered_frames(path) -> int:
    resolved = _resolved(path)
    with _pending_lock:
        slot = _pending.get(resolved)
    return len(slot[1]) if slot is not None else 0


def generation_frame_texts(path) -> list[str]:
    """Frame texts of the path's current generation, oldest first.

    Empty when this process holds no generation for the path or no
    corpus is loaded.
    """
    flush_corpus_appends()
    resolved = _resolved(path)
    with _pending_lock:
        tid = _generations.get(resolved)
    if tid is None:
        return []
    db = _corpus(corpus_dir(Path(resolved)))
    if db is None:
        return []
    texts = []
    while (text := _frame_text(db, tid, len(texts))) is not None:
        texts.append(text)
    return texts


atexit.register(flush_corpus_appends)


def mirror_con_text(path, text: str) -> None:
    """Insert one con blob. Identical path-and-text pairs are stored once."""
    if not text or not text.strip():
        return
    resolved = _resolved(path)
    db = _corpus(corpus_dir(Path(resolved)))
    if db is None:
        return
    try:
        db.append_trajectory_str(blob_id(resolved, text), text, source=resolved)
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
    resolved = _resolved(path)
    db = _corpus(corpus_dir(Path(resolved)))
    if db is None:
        return None
    key = blob_id(resolved, text)
    try:
        db.append_trajectory_str(key, text, source=resolved)
    except Exception as exc:
        if "already exists" not in str(exc):
            logger.warning("readcon-db store failed for %s: %s", resolved, exc)
            return None
    return key, 0


def load_frame_text(directory: Path, traj: int, frame_idx: int) -> str | None:
    """Frame text for a key in this corpus, or None when it is absent."""
    flush_corpus_appends()
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
    flush_corpus_appends()
    con_path = Path(path)
    try:
        text = con_path.read_text()
        resolved = str(con_path.resolve())
    except (OSError, UnicodeDecodeError):
        return None
    db = _corpus(corpus_dir(con_path))
    if db is None:
        return None
    try:
        return db.get_frame_text(blob_id(resolved, text), 0)
    except Exception:
        return None


def _path_is_compressed_con(path) -> bool:
    name = Path(path).name.lower()
    return name.endswith(".gz") or name.endswith(".zst")


def mirror_con_path(path) -> None:
    """Insert the con file's current bytes. A .gz or .zst path is skipped."""
    if _path_is_compressed_con(path):
        return
    con_path = Path(path)
    try:
        text = con_path.read_text()
    except (OSError, UnicodeDecodeError):
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
