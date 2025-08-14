#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NWChem i-PI Socket Driver: A Structured Python Implementation

This script provides a minimal, low-level Python implementation for driving NWChem
calculations using its i-PI-compliant socket interface.

This version has been refactored according to modern Python scripting practices
to enhance readability, robustness, and reusability. It can process multiple
geometries specified on the command line.
"""

# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "numpy==1.26.4",
#   "ase==3.23.0",
#   "click==8.1.8",
#   "rich==14.0.0",
# ]
# ///

import logging
import socket
import struct
import sys
import time
from dataclasses import dataclass, field
from enum import StrEnum

import ase.units as units
import click
import numpy as np
from rich.console import Console
from rich.logging import RichHandler
from rich.table import Table

# --- Setup Structured Logging with Rich ---
logging.basicConfig(
    level="INFO",
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, show_path=False)],
)
log = logging.getLogger("rich")
console = Console()

# --- Protocol and Connection Constants ---
MSGLEN = 12


# --- Data Structures using Enums and Dataclasses ---
class ProtoCommand(StrEnum):
    """Client-side commands sent to the NWChem server."""

    STATUS = "STATUS"
    INIT = "INIT"
    POSDATA = "POSDATA"
    GETFORCE = "GETFORCE"
    EXIT = "EXIT"


class ProtoStatus(StrEnum):
    """Server-side status replies from NWChem."""

    READY = "READY"
    NEEDINIT = "NEEDINIT"
    HAVEDATA = "HAVEDATA"
    FORCEREADY = "FORCEREADY"


@dataclass
class Geometry:
    """A dataclass to hold atomic geometry information."""

    name: str
    nat: int
    positions: np.ndarray
    cell: np.ndarray = field(default_factory=lambda: np.eye(3) * 15.0)

    def __post_init__(self):
        # Convert to numpy arrays if they aren't already
        self.positions = np.asanyarray(self.positions, dtype=float)
        self.cell = np.asanyarray(self.cell, dtype=float)
        # Basic validation
        if self.nat != self.positions.shape[0]:
            raise ValueError("`nat` must match the number of rows in `positions`.")


# --- Pre-defined Geometries for Demonstration ---
GEOMETRIES = {
    "equilibrium": Geometry(
        name="Water Molecule (Equilibrium)",
        nat=3,
        positions=[
            [0.0, 0.0, 0.1173],
            [0.0, 0.7572, -0.4692],
            [0.0, -0.7572, -0.4692],
        ],
    ),
    "stretched": Geometry(
        name="Water Molecule (Stretched)",
        nat=3,
        positions=[[0.0, 0.0, 0.12], [0.0, 0.85, -0.50], [0.0, -0.85, -0.50]],
    ),
}


# =============================================================================
# Communication Helper Functions
# (These functions remain unchanged)
# =============================================================================
def send_header(sock: socket.socket, command: ProtoCommand):
    sock.sendall(command.value.ljust(MSGLEN).encode("ascii"))


def recv_exact(sock: socket.socket, n: int) -> bytes:
    buf = b""
    while len(buf) < n:
        chunk = sock.recv(n - len(buf))
        if not chunk:
            log.critical("Socket closed unexpectedly by peer.")
            raise ConnectionAbortedError("Socket closed unexpectedly by peer")
        buf += chunk
    return buf


def recv_header(sock: socket.socket) -> str:
    return recv_exact(sock, MSGLEN).decode("ascii").strip()


# =============================================================================
# Core Logic Functions
# (The run_calculation function remains unchanged)
# =============================================================================
def run_calculation(conn: socket.socket, geo: Geometry):
    log.info(f"Preparing to run geometry: [bold cyan]{geo.name}[/bold cyan]")
    # 1. Check server status and handle NEEDINIT if necessary.
    send_header(conn, ProtoCommand.STATUS)
    reply = recv_header(conn)
    log.info(f"Initial server status: '{reply}'")
    if reply == ProtoStatus.NEEDINIT:
        log.warning("Server in 'NEEDINIT' state. Sending INIT to reset.")
        send_header(conn, ProtoCommand.INIT)
        conn.sendall(struct.pack("<i", 0))
        conn.sendall(struct.pack("<i", 1))
        conn.sendall(struct.pack("<b", 0))
        send_header(conn, ProtoCommand.STATUS)
        reply = recv_header(conn)
        log.info(f"Status after INIT: '{reply}'")
    if reply != ProtoStatus.READY:
        log.critical(f"Server not ready! Status is '{reply}'. Aborting.")
        raise RuntimeError(f"Server not ready! Status: {reply}")
    # 2. Convert units and prepare data payload.
    cell_bohr = geo.cell / units.Bohr
    invcell_bohr = np.linalg.inv(cell_bohr)
    positions_bohr = geo.positions / units.Bohr
    nat = geo.nat
    assert nat == len(positions_bohr)
    # 3. Send POSDATA.
    send_header(conn, ProtoCommand.POSDATA)
    conn.sendall(struct.pack(f"<{9}d", *cell_bohr.flatten(order="F")))
    conn.sendall(struct.pack(f"<{9}d", *invcell_bohr.flatten(order="F")))
    conn.sendall(struct.pack("<i", nat))
    conn.sendall(struct.pack(f"<{nat * 3}d", *positions_bohr.flatten("F")))
    log.info(f"Sent POSDATA for {nat} atoms.")
    # 4. Poll with STATUS until calculation is done (HAVEDATA).
    log.info("Polling for results...")
    while True:
        send_header(conn, ProtoCommand.STATUS)
        status = recv_header(conn)
        if status == ProtoStatus.HAVEDATA:
            log.info("Calculation complete ('HAVEDATA' received).")
            break
        time.sleep(0.1)
    # 5. Request and receive results.
    send_header(conn, ProtoCommand.GETFORCE)
    hdr = recv_header(conn)
    if hdr != ProtoStatus.FORCEREADY:
        log.critical(f"Expected '{ProtoStatus.FORCEREADY}', but got '{hdr}'.")
        raise RuntimeError(f"Unexpected reply: {hdr}")
    (energy_ha,) = struct.unpack("<d", recv_exact(conn, 8))
    (nat_back,) = struct.unpack("<i", recv_exact(conn, 4))
    assert nat_back == nat, "Atom count mismatch in force data!"
    forces_flat = struct.unpack(f"<{nat * 3}d", recv_exact(conn, nat * 3 * 8))
    _ = struct.unpack("<9d", recv_exact(conn, 9 * 8))
    (extra_len,) = struct.unpack("<i", recv_exact(conn, 4))
    _ = recv_exact(conn, extra_len)
    forces_ha_bohr = np.array(forces_flat).reshape((nat, 3))
    table = Table(
        title=f"Results for {geo.name}",
        title_style="bold green",
        show_header=True,
        header_style="bold magenta",
    )
    table.add_column("Property", justify="right")
    table.add_column("Value", justify="left")
    table.add_row("Energy (Hartree)", f"{energy_ha:.8f}")
    table.add_row("Max Force (Ha/Bohr)", f"{np.max(np.abs(forces_ha_bohr)):.8f}")
    console.print(table)


# =============================================================================
# Command-Line Interface
# =============================================================================


@click.command()
@click.option(
    "--host",
    default="127.0.0.1",
    show_default=True,
    help="Host address for the NWChem server.",
)
@click.option(
    "--port",
    default=9999,
    show_default=True,
    type=int,
    help="Port for the NWChem server.",
)
@click.option(
    "--geometry",
    "geometry_names",  # The variable now stores a tuple of names
    type=click.Choice(GEOMETRIES.keys(), case_sensitive=False),
    multiple=True,  # This allows the option to be used more than once
    required=True,
    help="Name of the geometry to calculate. Can be specified multiple times.",
)
def main(host, port, geometry_names):
    """
    A client driver for running NWChem calculations via an i-PI socket.

    This script connects to a running NWChem instance that has been started
    in socket mode and sends it one or more pre-defined geometries to calculate.
    """
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as server:
        server.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        server.bind((host, port))
        server.listen(1)
        log.info(f"Listening for NWChem on [bold yellow]{host}:{port}[/bold yellow]...")

        try:
            conn, addr = server.accept()
        except KeyboardInterrupt:
            log.warning("\nInterrupted by user. Exiting.")
            sys.exit(0)

        with conn:
            log.info(
                f"Connection established with NWChem: [bold yellow]{addr}[/bold yellow]"
            )
            try:
                # Loop through all geometries provided on the command line
                for name in geometry_names:
                    selected_geo = GEOMETRIES[name]
                    run_calculation(conn, selected_geo)

                log.info(
                    "[bold green]All calculations finished successfully.[/bold green]"
                )
            except (RuntimeError, ConnectionAbortedError, ValueError) as e:
                log.critical(f"A critical error occurred: {e}")
                sys.exit(1)
            finally:
                log.info("Sending EXIT command to NWChem for graceful shutdown.")
                send_header(conn, ProtoCommand.EXIT)


if __name__ == "__main__":
    main()
