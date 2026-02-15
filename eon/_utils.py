from pathlib import Path
import sys
import subprocess
import tempfile
import logging
import typing as typ
from dataclasses import dataclass

from eon import fileio as eio
from eon import config as econf


@dataclass
class ScriptConfig:
    """Configuration for running an external atom list script."""

    script_path: Path
    scratch_path: Path
    root_path: Path

    def __post_init__(self):
        """
        Post-initialization processing to validate and resolve paths.
        This method is automatically called by the dataclass constructor.
        """
        if not self.script_path.is_absolute():
            self.script_path = self.root_path / self.script_path
        if not self.scratch_path.is_absolute():
            self.scratch_path = self.root_path / self.scratch_path

    @classmethod
    def from_eon_config(cls, config: econf.ConfigClass) -> typ.Self:
        """
        Factory method to create a ScriptConfig instance from the main EON config.
        """
        return cls(
            script_path=Path(config.displace_atom_kmc_state_script),
            scratch_path=Path(config.path_scratch),
            root_path=Path(config.path_root),
        )


def gen_ids_from_con(sconf: ScriptConfig, reactant, logger: logging.Logger):
    if not sconf.script_path.is_file():
        logger.error(f"displace_atom_list_script not found: {sconf.script_path}")
        return ""

    # Use a secure temporary file to pass the structure to the script
    # The file is automatically deleted when the 'with' block is exited.
    sconf.scratch_path.mkdir(parents=True, exist_ok=True)
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
                [sys.executable, str(sconf.script_path), tmpf.name],
                capture_output=True,
                text=True,
                check=True,
                cwd=sconf.root_path,
            )
            atom_list_str = proc.stdout.strip()
            return atom_list_str
        except subprocess.CalledProcessError as e:
            logger.error(f"Error running displace_atom_list_script '{sconf.script_path}':")
            logger.error(f"Stderr: {e.stderr.strip()}")
            sys.exit(1)
        except Exception as e:
            logger.error(
                f"An unexpected error occurred while running the displacement script: {e}"
            )
            sys.exit(1)
