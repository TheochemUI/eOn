"""Unit tests for ConfigClass injection and communicator independence."""
from __future__ import annotations

import os
import tempfile
import textwrap
import unittest
from pathlib import Path

import numpy as np

from eon import atoms
from eon import communicator
from eon.config import ConfigClass
from eon.structure import Structure


MINIMAL_INI = textwrap.dedent(
    """
    [Main]
    job = minimization
    temperature = 300
    [Communicator]
    type = local
    """
)


def _write_ini(directory: str, name: str = "config.ini") -> str:
    path = str(Path(directory) / name)
    with open(path, "w") as fh:
        fh.write(MINIMAL_INI)
    return path


def _make_structure(pos, box_scale=10.0):
    s = Structure(1)
    s.r[0] = pos
    s.free[0] = 1.0
    s.names[0] = "H"
    s.mass[0] = 1.0
    s.box = np.eye(3) * box_scale
    return s


class TestCommunicatorIndependence(unittest.TestCase):
    def setUp(self):
        communicator.reset_communicators()

    def tearDown(self):
        communicator.reset_communicators()

    def test_get_communicator_requires_config(self):
        with self.assertRaises(TypeError):
            communicator.get_communicator()  # type: ignore[call-arg]

    def test_two_configs_get_independent_local_comms(self):
        with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
            path_a = _write_ini(d1, "a.ini")
            path_b = _write_ini(d2, "b.ini")
            # Local validates absolute client paths; provide dummy executables.
            client_a = str(Path(d1) / "dummy_client")
            client_b = str(Path(d2) / "dummy_client")
            for c in (client_a, client_b):
                with open(c, "w") as fh:
                    fh.write("#!/bin/sh\n")
                os.chmod(c, 0o755)
            old = os.getcwd()
            try:
                os.chdir(d1)
                cfg_a = ConfigClass()
                cfg_a.init(path_a)
                cfg_a.path_root = d1
                cfg_a.path_scratch = str(Path(d1) / "scratch")
                cfg_a.comm_type = "local"
                cfg_a.comm_local_client = client_a
                cfg_a.comm_local_ncpus = 1
                cfg_a.comm_job_bundle_size = 1

                os.chdir(d2)
                cfg_b = ConfigClass()
                cfg_b.init(path_b)
                cfg_b.path_root = d2
                cfg_b.path_scratch = str(Path(d2) / "scratch")
                cfg_b.comm_type = "local"
                cfg_b.comm_local_client = client_b
                cfg_b.comm_local_ncpus = 2
                cfg_b.comm_job_bundle_size = 1
            finally:
                os.chdir(old)

            comm_a = communicator.get_communicator(cfg_a)
            comm_b = communicator.get_communicator(cfg_b)
            self.assertIsNot(comm_a, comm_b)
            self.assertEqual(comm_a.scratchpath, cfg_a.path_scratch)
            self.assertEqual(comm_b.scratchpath, cfg_b.path_scratch)
            self.assertIs(communicator.get_communicator(cfg_a), comm_a)


class TestAtomsMatchNoGlobalConfig(unittest.TestCase):
    def test_match_uses_explicit_knobs(self):
        a = _make_structure([1.0, 1.0, 1.0])
        b = _make_structure([1.0, 1.0, 1.0])
        self.assertTrue(
            atoms.match(
                a, b, 0.1, 3.0, False,
                check_rotation=False, use_identical=False,
            )
        )
        b2 = _make_structure([5.0, 5.0, 5.0])
        self.assertFalse(
            atoms.match(
                a, b2, 0.1, 3.0, False,
                check_rotation=False, use_identical=False,
            )
        )


class TestModifyConfigUsesPathArg(unittest.TestCase):
    def test_modify_config_reads_given_path(self):
        from eon import fileio as io

        with tempfile.TemporaryDirectory() as d:
            path = _write_ini(d)
            out = io.modify_config(path, [("Main", "temperature", "999")])
            text = out.getvalue()
            self.assertIn("999", text)




class TestAkmcRequiresConfig(unittest.TestCase):
    def test_akmc_without_config_raises(self):
        from eon import akmc
        with self.assertRaises(TypeError):
            akmc.akmc()

    def test_basinhopping_without_config_raises(self):
        from eon import basinhopping
        with self.assertRaises(TypeError):
            basinhopping.basinhopping()

    def test_parallelreplica_without_config_raises(self):
        from eon import parallelreplica
        with self.assertRaises(TypeError):
            parallelreplica.parallelreplica()

if __name__ == "__main__":
    unittest.main()
