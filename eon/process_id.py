"""Content-addressed process ids (xxHash).

Process table rows and procdata files are keyed by integer ids. A counter
based on ``len(procs)`` freezes when duplicate ids collapse the dict; a
counter based on ``max(id)+1`` still couples identity to registration
order and external table rewrites.

Instead, seed the id from a stable xxh64 of the process payload (saddle
geometry bytes, barrier tag, direction). Geometric uniqueness remains
``find_repeat``; the hash only assigns free keys that do not depend on
table length. Rare collisions (or an already-occupied pure-content hash)
are resolved by salting until free.
"""

from __future__ import annotations

from typing import Union

import xxhash

BytesLike = Union[bytes, bytearray, memoryview, str]


def _as_bytes(part: BytesLike) -> bytes:
    if isinstance(part, str):
        return part.encode("utf-8")
    if isinstance(part, memoryview):
        return part.tobytes()
    if isinstance(part, bytearray):
        return bytes(part)
    return part


def process_id_from_parts(*parts: BytesLike, salt: int = 0) -> int:
    """xxh64 over *parts* (null-separated), optionally salted.

    Returns a non-negative 63-bit int so it is safe as a dict key and as
    the ``%d`` field in processtable / ``procdata/*_%d.*`` paths.
    """
    h = xxhash.xxh64()
    for part in parts:
        h.update(_as_bytes(part))
        h.update(b"\0")
    if salt:
        h.update(salt.to_bytes(8, "little", signed=False))
    return h.intdigest() & 0x7FFFFFFFFFFFFFFF


def allocate_unique_process_id(
    occupied: dict | set | None,
    *parts: BytesLike,
    max_salts: int = 1024,
) -> int:
    """Return a content-addressed id not already in *occupied*.

    Salt 0 is the pure content hash. If that id is already present (true
    re-registration of the same payload, or an xxh64 collision), salt
    until free. Raises if *max_salts* is exhausted.
    """
    keys = occupied if occupied is not None else ()
    for salt in range(max_salts):
        pid = process_id_from_parts(*parts, salt=salt)
        if pid not in keys:
            return pid
    raise RuntimeError(
        "unable to allocate a free xxh64 process id after %d salts" % max_salts
    )
