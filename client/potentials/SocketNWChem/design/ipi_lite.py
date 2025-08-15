#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NWChem i-PI Socket Driver: A Minimal Python Implementation

This script provides a minimal, low-level Python implementation for driving NWChem
calculations using its i-PI-compliant socket interface.

Design
------
The goal is to interact directly with the NWChem executable as a calculation
server, without relying on high-level libraries like ASE. This offers a clear
view of the underlying communication protocol.

Several notes are relevant about the networking protocol specific to the
implementation in NWChem:

1.  Server Setup (NWChem Side)
    - The NWChem process is launched with an input file that specifies the socket
      driver, e.g., `task scf optimize` and a `driver` block with
      `socket ipi_client <host>:<port>`.
    - The `geometry` block in this initial input file is CRITICAL. NWChem uses
      it at startup to determine the number of atoms and their types,
      allocating memory accordingly. This configuration CANNOT be changed
      during the socket session. All subsequent calculations sent over the
      socket must have the same number of atoms.

2.  Data Serialization (The Binary Format)
    - All multi-byte binary data (doubles, integers) MUST be serialized in
      LITTLE-ENDIAN format (e.g., using '<' in Python's `struct` module).
    - All multi-dimensional arrays (cell matrices, positions) MUST be flattened
      in Fortran 'F' order (`.flatten(order='F')` in NumPy) to match Fortran's
      column-major memory layout.
    - The i-PI protocol mandates the use of ATOMIC UNITS. All lengths
      (positions, cell vectors) must be sent in BOHR, and energies/forces will
      be received in HARTREE and HARTREE/BOHR.

3.  Communication Protocol (The State Machine Handshake)
    - The NWChem server operates as a state machine. The client must respect
      these states by using the `STATUS` command to query the server.
    - The full cycle for a calculation is:
      READY -> (POSDATA sent) -> HAVEDATA -> (GETFORCE sent) -> NEEDINIT
    - To perform a new calculation after the first one, the client must detect
      the `NEEDINIT` state and send an `INIT` command to reset the server
      back to the `READY` state.
    - When all calculations are complete, the client should send `EXIT` to
      terminate the NWChem process gracefully.
"""
import socket
import struct
import time
import numpy as np
import ase.units as units

# --- Protocol and Connection Constants ---
MSGLEN = 12
HOST, PORT = "127.0.0.1", 9999


# =============================================================================
# Communication Helper Functions
# =============================================================================

def send_header(sock: socket.socket, txt: str):
    """Sends a 12-byte, left-justified, ascii-encoded header."""
    sock.sendall(txt.ljust(MSGLEN).encode("ascii"))


def recv_exact(sock: socket.socket, n: int) -> bytes:
    """Receives exactly n bytes, handling short reads."""
    buf = b""
    while len(buf) < n:
        chunk = sock.recv(n - len(buf))
        if not chunk:
            raise ConnectionError("Socket closed unexpectedly by peer")
        buf += chunk
    return buf


def recv_header(sock: socket.socket) -> str:
    """Receives and decodes a 12-byte header."""
    return recv_exact(sock, MSGLEN).decode("ascii").strip()


# =============================================================================
# Main Execution Logic
# =============================================================================

if __name__ == "__main__":
    # Define geometries to run. All coordinates and cell vectors in ANGSTROMS.
    geometries_to_run = [
        {
            "name": "Water Molecule (Equilibrium)",
            "nat": 3,
            "positions": np.array(
                [
                    [0.0, 0.0, 0.1173],
                    [0.0, 0.7572, -0.4692],
                    [0.0, -0.7572, -0.4692],
                ]
            ),
            "cell": np.eye(3) * 15.0,
        },
        {
            "name": "Water Molecule (Stretched)",
            "nat": 3,
            "positions": np.array(
                [[0.0, 0.0, 0.12], [0.0, 0.85, -0.50], [0.0, -0.85, -0.50]]
            ),
            "cell": np.eye(3) * 15.0,
        },
    ]

    # --- Server Setup ---
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as server:
        server.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        server.bind((HOST, PORT))
        server.listen(1)
        print(f"Listening for NWChem on {HOST}:{PORT}")
        conn, addr = server.accept()
        with conn:
            print(f"Connected by NWChem: {addr}")

            # --- Main Calculation Loop ---
            for i, geo in enumerate(geometries_to_run):
                print(f"\n{'='*20} GEOMETRY {i+1}: {geo['name']} {'='*20}")

                # 1. Check server status and handle NEEDINIT if necessary.
                send_header(conn, "STATUS")
                reply = recv_header(conn)
                print(f"Initial status: {reply}")

                if reply == "NEEDINIT":
                    print("Server needs reset. Sending INIT command.")
                    send_header(conn, "INIT")
                    conn.sendall(struct.pack("<i", 0)) # bead index
                    conn.sendall(struct.pack("<i", 1)) # nbytes of extra data
                    conn.sendall(struct.pack("<b", 0)) # a single dummy byte

                    # After sending INIT, re-check status to confirm it's READY.
                    send_header(conn, "STATUS")
                    reply = recv_header(conn)
                    print(f"Status after INIT: {reply}")

                if reply != "READY":
                    raise RuntimeError(f"Server not ready! Status: {reply}")

                # 2. Convert units from Angstrom to Bohr.
                cell_bohr = np.eye(3)
                invcell_bohr = np.eye(3)
                positions_bohr = geo["positions"] / units.Bohr
                nat = geo["nat"]

                # 3. Send POSDATA with data in correct units and memory order.
                send_header(conn, "POSDATA")
                conn.sendall(struct.pack(f"<{9}d", *cell_bohr.flatten(order="F")))
                conn.sendall(struct.pack(f"<{9}d", *invcell_bohr.flatten(order="F")))
                conn.sendall(struct.pack("<i", nat))
                conn.sendall(struct.pack(f"<{nat*3}d", *positions_bohr.flatten("F")))
                print(f"Sent POSDATA for {nat} atoms.")

                # 4. Poll with STATUS until calculation is done.
                print("Polling for results...", end="", flush=True)
                while True:
                    send_header(conn, "STATUS")
                    status = recv_header(conn)
                    if status == "HAVEDATA":
                        print(" Done (HAVEDATA received).")
                        break
                    time.sleep(0.05) # Prevent a busy-wait.

                # 5. Request and receive results.
                send_header(conn, "GETFORCE")
                hdr = recv_header(conn)
                if hdr != "FORCEREADY":
                    raise RuntimeError(f"Expected FORCEREADY, got {hdr}")

                # Unpack the entire payload correctly
                energy_ha, = struct.unpack("<d", recv_exact(conn, 8))
                nat_back, = struct.unpack("<i", recv_exact(conn, 4))
                forces_flat = struct.unpack(f"<{nat_back*3}d", recv_exact(conn, nat_back*3*8))
                virial_flat = struct.unpack("<9d", recv_exact(conn, 9*8))
                extra_len, = struct.unpack("<i", recv_exact(conn, 4))
                extras = recv_exact(conn, extra_len)

                forces_ha_bohr = np.array(forces_flat).reshape((nat_back, 3))

                print(f"  Energy: {energy_ha:.8f} Hartree")
                print( "  Forces (Ha/Bohr):\n", forces_ha_bohr)

            # --- Graceful Shutdown ---
            print("\nAll calculations finished.")
            send_header(conn, "EXIT")
            print("Sent EXIT to NWChem. Connection closed.")
