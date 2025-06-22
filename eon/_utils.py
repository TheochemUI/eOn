#!/usr/bin/env python3

from pathlib import Path
import sys
import subprocess
import tempfile
import logging
import typing as typ
from dataclasses import dataclass

from eon import fileio as eio
from eon import atoms as eatm
from eon import config as econf


@dataclass
class ScriptConfig:
    """Configuration for running an external atom list script."""

    script_path: Path
    scratch_path: Path
    root_path: Path

    @classmethod
    def from_eon_config(cls, config: econf.ConfigClass) -> typ.Self:
        """
        Factory method to create a ScriptConfig instance from the main EON config.
        """
        # The logic for resolving the script path now lives here,
        # encapsulated within the class itself.
        script_path = Path(config.displace_atom_list_script)
        if not script_path.is_absolute():
            script_path = Path(config.root_path) / script_path

        return cls(
            script_path=script_path,
            scratch_path=Path(config.path_scratch),
            root_path=Path(config.path_root),
        )


def gen_ids_from_con(sconf: ScriptConfig, reactant: eatm.Atoms, logger: logging.Logger):
    script_path = Path(sconf.script_path)
    if not script_path.is_absolute():
        script_path = sconf.root_path / sconf.script_path

    if not script_path.is_file():
        logger.error(f"displace_atom_list_script not found: {script_path}")
        sys.exit(1)

    # Use a secure temporary file to pass the structure to the script
    # The file is automatically deleted when the 'with' block is exited.
    with tempfile.NamedTemporaryFile(
        mode="w+", delete=True, suffix=".con", dir=sconf.scratch_path
    ) as tmpf:
        # Save the displaced structure to the temporary file
        eio.savecon(tmpf, reactant)
        # Ensure all data is written to disk before the script reads it
        tmpf.flush()

        try:
            # Execute the script, passing the temporary filename as an argument
            proc = subprocess.run(
                [sys.executable, str(script_path), tmpf.name],  # Pass filename
                capture_output=True,
                text=True,
                check=True,
                cwd=sconf.root_path,
            )
            atom_list_str = proc.stdout.strip()
            return atom_list_str
        except subprocess.CalledProcessError as e:
            logger.error(f"Error running displace_atom_list_script '{script_path}':")
            logger.error(f"Stderr: {e.stderr.strip()}")
            sys.exit(1)
        except Exception as e:
            logger.error(
                f"An unexpected error occurred while running the displacement script: {e}"
            )
            sys.exit(1)
