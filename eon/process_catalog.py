"""Process catalog for aKMC.

Frames live in the run's readcon-db corpus. The saddle frame carries the
process mode as its readcon displacements section; the line-2 JSON names
that section. The barrier, the prefactor, the same mode components, and
the frame keys live on ``amsel.KdbProcess``. A good saddle is stored when
it is registered. The next search of a matching state refines from the
stored saddle, then a random displacement follows when no suggestion
remains. ``Paths.kdb`` is the directory ``amsel.KdbStore`` opens.
"""

from __future__ import annotations

import hashlib
import json
import logging
from pathlib import Path

import numpy as np

from eon import fileio as io
from eon.concorpus import corpus_dir, store_frame_text
from eon.structure import structure_order

logger = logging.getLogger("kdb")


def pack_frame_key(traj_id: int, frame_idx: int) -> bytes:
    """12-byte readcon-db key: traj_id then frame index, both big-endian."""
    return int(traj_id).to_bytes(8, "big") + int(frame_idx).to_bytes(4, "big")


def unpack_frame_key(blob: bytes) -> tuple[int, int]:
    raw = bytes(blob)
    if len(raw) != 12:
        raise ValueError(f"frame key is {len(raw)} bytes, expected 12")
    return int.from_bytes(raw[:8], "big"), int.from_bytes(raw[8:], "big")


def env_hash(atoms) -> bytes:
    """16-byte hash of element names and positions rounded to 1e-4 angstrom.

    Rounding is finer than ``kdb_dc``, so the same minimum matches and a
    real hop does not.
    """
    positions = np.round(np.asarray(atoms.r, dtype=float), 4)
    names = [str(name) for name in atoms.names]
    payload = json.dumps(
        {"names": names, "r": positions.tolist()},
        separators=(",", ":"),
    )
    return hashlib.blake2s(payload.encode(), digest_size=16).digest()


def _floats(config):
    return float(config.kdb_nf), float(config.kdb_dc), float(config.kdb_mac)


def _open_store(config):
    """Open the catalog at ``config.kdb_path``, or None when amsel is absent."""
    try:
        from amsel import KdbStore
    except ImportError as exc:
        # A missing dependency inside amsel is not "the module was not found".
        if getattr(exc, "name", None) == "amsel":
            logger.error("amsel is not installed; the process catalog is closed")
        else:
            logger.error("process catalog import failed: %s", exc)
        return None
    path = Path(config.kdb_path)
    try:
        path.mkdir(parents=True, exist_ok=True)
        return KdbStore(str(path))
    except Exception:
        logger.exception("amsel.KdbStore failed to open %s", path)
        return None


def _match_dir(config, state) -> Path:
    return Path(config.kdb_scratch_path) / "kdbmatches" / f"state_{int(state.number)}"


def _queried_path(config) -> Path:
    return Path(config.kdb_scratch_path) / "queried"


def _integer_tokens(text: str) -> list[int]:
    numbers = []
    for token in text.split():
        try:
            numbers.append(int(token))
        except ValueError:
            continue
    return numbers


def was_queried(state, config) -> bool:
    path = _queried_path(config)
    if not path.is_file():
        return False
    numbers = _integer_tokens(path.read_text())
    return int(state.number) in numbers


def _mark_queried(state, config) -> None:
    path = _queried_path(config)
    path.parent.mkdir(parents=True, exist_ok=True)
    numbers = _integer_tokens(path.read_text()) if path.is_file() else []
    number = int(state.number)
    if number not in numbers:
        numbers.append(number)
    path.write_text("".join(f"{n}\n" for n in numbers))


def _consumed_path(config, state) -> Path:
    return _match_dir(config, state) / "consumed"


def _consumed(config, state) -> set[int]:
    path = _consumed_path(config, state)
    if not path.is_file():
        return set()
    out = set()
    for line in path.read_text().split():
        try:
            out.add(int(line))
        except ValueError:
            continue
    return out


def _mark_consumed(config, state, index: int) -> None:
    path = _consumed_path(config, state)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        handle.write(f"{int(index)}\n")


def _mode_cosine(mode, reactant, saddle) -> float:
    disp = np.asarray(saddle.r, dtype=float) - np.asarray(reactant.r, dtype=float)
    vec = np.asarray(mode, dtype=float).reshape(-1)
    flat = disp.reshape(-1)
    n = min(vec.size, flat.size)
    if n == 0:
        return 0.0
    vec = vec[:n]
    flat = flat[:n]
    denom = float(np.linalg.norm(vec) * np.linalg.norm(flat))
    if denom == 0.0:
        return 0.0
    return float(np.dot(vec, flat) / denom)


def _max_distance(a, b) -> float:
    ar = np.asarray(a.r, dtype=float)
    br = np.asarray(b.r, dtype=float)
    n = min(len(ar), len(br))
    if n == 0:
        return 0.0
    delta = ar[:n] - br[:n]
    return float(np.max(np.linalg.norm(delta, axis=1)))


def _load_frame(corpus_directory: Path, key: bytes):
    from eon.concorpus import load_frame_text

    traj_id, frame_idx = unpack_frame_key(key)
    text = load_frame_text(corpus_directory, traj_id, frame_idx)
    if not text:
        raise OSError(f"readcon-db has no frame {traj_id}:{frame_idx}")
    return text, traj_id, frame_idx


def _mode_in_structure_order(frame) -> np.ndarray | None:
    """Displacements section in the same row order as ``Structure.r``.

    The file groups atoms by species. ``mode_<id>.dat`` follows ``atom_id``
    order, which is the order :meth:`Structure.from_conframe` restores.
    """
    disp = frame.disp
    if disp is None:
        return None
    vec = np.asarray(disp, dtype=float)
    ids = np.array([atom.atom_id for atom in frame.atoms], dtype=np.uint64)
    return vec[structure_order(ids)]


def _saddle_text_with_mode(saddle_text: str, mode) -> str:
    """Return the saddle con with ``mode`` in the displacements section.

    ``mode`` is Nx3 in ``atom_id`` order, the order of ``mode_<id>.dat``.
    """
    import readcon

    frames = readcon.read_con_string(saddle_text)
    if not frames:
        raise ValueError("saddle con has no frame")
    frame = frames[0]
    vec = np.asarray(mode, dtype=float).reshape(-1, 3)
    atoms = list(frame.atoms)
    if len(atoms) != len(vec):
        raise ValueError(
            f"mode has {len(vec)} rows and the saddle has {len(atoms)} atoms"
        )
    ids = np.array([atom.atom_id for atom in atoms], dtype=np.uint64)
    order = structure_order(ids)
    for struct_i, file_i in enumerate(order):
        atom = atoms[int(file_i)]
        atom.dx = float(vec[struct_i, 0])
        atom.dy = float(vec[struct_i, 1])
        atom.dz = float(vec[struct_i, 2])
    return readcon.write_con_string([frame], 17)


def _suggestion_mode(
    process,
    reactant,
    saddle,
    mac: float,
    saddle_text: str | None = None,
) -> np.ndarray:
    """Direction for a refine.

    The mode in the saddle frame's displacements section is used when its
    cosine with the reactant-to-saddle vector is at least ``mac``. A
    negative cosine flips the mode. Below ``mac`` the direction is that
    vector, so a curved path is still offered. A frame with no
    displacements section uses the mode stored on the catalog row.
    """
    stored = None
    if saddle_text:
        import readcon

        frames = readcon.read_con_string(saddle_text)
        if frames:
            stored = _mode_in_structure_order(frames[0])
    if stored is None:
        stored = np.asarray(list(process.mode), dtype=float).reshape(-1, 3)
    cosine = _mode_cosine(stored, reactant, saddle)
    if abs(cosine) >= float(mac) and float(np.linalg.norm(stored)) > 0.0:
        if cosine < 0.0:
            return -stored
        return stored
    logger.info(
        "kdb mode cosine=%.3f is below mac=%s; using the reactant-to-saddle vector",
        cosine,
        mac,
    )
    return np.asarray(saddle.r, dtype=float) - np.asarray(reactant.r, dtype=float)


def _accepts(process, reactant, corpus_directory: Path, nf: float, dc: float):
    """Return (reactant, saddle) frames when the stored reactant matches."""
    if not process.saddle_frame_key or not process.reactant_frame_key:
        return None
    try:
        reactant_text, _, _ = _load_frame(corpus_directory, process.reactant_frame_key)
        saddle_text, _, _ = _load_frame(corpus_directory, process.saddle_frame_key)
    except Exception:
        logger.exception("readcon-db frame load failed")
        return None
    stored_reactant = io.loadcon(io.StringIO(reactant_text))
    stored_saddle = io.loadcon(io.StringIO(saddle_text))
    allowed = float(dc) * (1.0 + float(nf))
    mismatch = _max_distance(reactant, stored_reactant)
    if mismatch > allowed:
        logger.info(
            "kdb reactant mismatch %.4f A, allowed %.4f A",
            mismatch,
            allowed,
        )
        return None
    return stored_reactant, stored_saddle


def insert(state, process_id, config) -> bool:
    """Store one good process. Returns True when the catalog accepted it.

    Called from process registration. Confidence is not consulted: a
    process with zero repeats is stored.
    """
    nf, dc, mac = _floats(config)
    logger.info(
        "kdb insert nf=%r dc=%r mac=%r path=%s process=%s",
        nf,
        dc,
        mac,
        config.kdb_path,
        process_id,
    )
    row = state.procs[process_id]
    reactant_path = Path(state.proc_reactant_path(process_id))
    saddle_path = Path(state.proc_saddle_path(process_id))
    product_path = Path(state.proc_product_path(process_id))
    reactant_text = reactant_path.read_text()
    saddle_text = saddle_path.read_text()
    product_text = product_path.read_text()
    mode = np.asarray(state.get_process_mode(process_id), dtype=float)
    try:
        saddle_text = _saddle_text_with_mode(saddle_text, mode)
    except Exception:
        logger.exception("saddle frame did not take the mode")
        return False
    # The mode lives on the catalog row. A corpus that rejects the
    # displacements section still stores the saddle geometry.
    saddle_for_corpus = saddle_text
    if store_frame_text(saddle_path, saddle_text) is None:
        saddle_for_corpus = saddle_path.read_text()
    keys = []
    for path, text in (
        (reactant_path, reactant_text),
        (saddle_path, saddle_for_corpus),
        (product_path, product_text),
    ):
        key = store_frame_text(path, text)
        if key is None:
            logger.error(
                "readcon-db did not store %s; the process was not catalogued",
                path,
            )
            return False
        keys.append(pack_frame_key(*key))
    reactant_key, saddle_key, product_key = keys
    product = io.loadcon(str(product_path))
    store = _open_store(config)
    if store is None:
        return False
    # Key by the state reactant. The client's reactant file is the same
    # minimum rewritten, and a 1e-4 A grid would miss it.
    env = env_hash(state.get_reactant())
    try:
        existing = store.lookup(env)
    except Exception:
        logger.exception("amsel.KdbStore lookup failed during insert")
        return False
    for process in existing:
        if bytes(process.saddle_frame_key) == saddle_key and abs(
            float(process.barrier_ev) - float(row["barrier"])
        ) < 1e-8:
            return True
    from amsel import KdbProcess

    record = KdbProcess(
        b"",
        b"",
        float(row["barrier"]),
        float(row["prefactor"]),
        discovery_temperature=float(getattr(config, "main_temperature", 0.0)),
        usage_hint="RefineFirst",
        product_env_hash=env_hash(product),
        reactant_frame_key=reactant_key,
        saddle_frame_key=saddle_key,
        product_frame_key=product_key,
        mode=[float(value) for value in mode.reshape(-1)],
        metadata_json=json.dumps(
            {
                "nf": nf,
                "dc": dc,
                "mac": mac,
                "process_id": int(process_id),
            }
        ),
    )
    try:
        store.insert(env, record)
    except Exception:
        logger.exception("amsel.KdbStore insert failed")
        return False
    logger.info(
        "kdb insert barrier=%.4f eV saddle frame %d:%d into readcon.db",
        float(row["barrier"]),
        *unpack_frame_key(saddle_key),
    )
    return True


def _corpus_directory(state) -> Path:
    reactant = getattr(state, "reactant_path", None)
    if reactant:
        return corpus_dir(Path(reactant))
    return corpus_dir(Path(state.proc_reactant_path(0)))


def _materialize(state, config, store) -> None:
    """Write this state's unused suggestions. Other states' files stay."""
    nf, dc, mac = _floats(config)
    reactant = state.get_reactant()
    try:
        processes = store.lookup(env_hash(reactant))
    except Exception:
        logger.exception("amsel.KdbStore lookup failed")
        raise
    directory = _match_dir(config, state)
    directory.mkdir(parents=True, exist_ok=True)
    corpus_directory = _corpus_directory(state)
    done = _consumed(config, state)
    for index, process in enumerate(processes):
        if index in done:
            continue
        pair = _accepts(process, reactant, corpus_directory, nf, dc)
        if pair is None:
            continue
        stored_reactant, stored_saddle = pair
        saddle_text, traj_id, frame_idx = _load_frame(
            corpus_directory, process.saddle_frame_key
        )
        saddle_path = directory / f"SADDLE_{index}"
        mode_path = directory / f"MODE_{index}"
        barrier_path = directory / f"BARRIER_{index}"
        done_path = directory / f".done_{index}"
        if not saddle_path.is_file():
            saddle_path.write_text(saddle_text)
            mode = _suggestion_mode(
                process, stored_reactant, stored_saddle, mac, saddle_text
            )
            io.save_mode(str(mode_path), mode)
            barrier_path.write_text(f"{float(process.barrier_ev):.10f}\n")
            (directory / f"KEY_{index}").write_text(f"{traj_id} {frame_idx}\n")
            done_path.write_text("")


def query(state, config) -> bool:
    """Look up suggestions for ``state``.

    The state number is appended to ``queried`` only after lookup returns.
    An import failure or an open failure leaves that file unchanged, and
    does not delete another state's ``SADDLE_`` files.
    """
    nf, dc, mac = _floats(config)
    logger.info(
        "kdb query nf=%r dc=%r mac=%r path=%s",
        nf,
        dc,
        mac,
        config.kdb_path,
    )
    store = _open_store(config)
    if store is None:
        return False
    try:
        _materialize(state, config, store)
    except Exception:
        logger.exception("kdb query failed for state %s", getattr(state, "number", "?"))
        return False
    _mark_queried(state, config)
    return True


def make_suggestion(config, state):
    """Return one ``(displacement, mode)`` pair, or ``(None, None)``.

    Pending ``SADDLE_`` files for this state are read first. A process
    stored after the last query is materialized from ``amsel.KdbStore``
    and the saddle frame is loaded from readcon-db. Each process is
    offered once; the caller then uses a random displacement.
    """
    store = _open_store(config)
    if store is None:
        return None, None
    try:
        _materialize(state, config, store)
    except Exception:
        logger.exception("kdb suggestion lookup failed")
        return None, None
    directory = _match_dir(config, state)
    if not directory.is_dir():
        return None, None
    dones = sorted(p for p in directory.glob(".done_*") if p.is_file())
    if not dones:
        return None, None
    number = dones[0].name.split("_", 1)[1]
    saddle_path = directory / f"SADDLE_{number}"
    mode_path = directory / f"MODE_{number}"
    barrier_path = directory / f"BARRIER_{number}"
    key_path = directory / f"KEY_{number}"
    try:
        displacement = io.loadcon(str(saddle_path))
        mode = io.load_mode(str(mode_path))
        barrier = float(barrier_path.read_text().strip())
        traj_id, frame_idx = key_path.read_text().split()
    except (OSError, ValueError):
        logger.exception("kdb suggestion file %s is unreadable", saddle_path)
        return None, None
    for path in (dones[0], saddle_path, mode_path, barrier_path, key_path):
        if path.is_file():
            path.unlink()
    _mark_consumed(config, state, int(number))
    logger.info(
        "KDB suggestion barrier=%.4f eV saddle frame %s:%s from readcon.db",
        barrier,
        traj_id,
        frame_idx,
    )
    return displacement, mode
