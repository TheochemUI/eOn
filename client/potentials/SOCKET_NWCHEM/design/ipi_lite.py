import socket
import struct
import time
import numpy as np
import ase.units as units

MSGLEN = 12

def send_header(sock, txt):
    sock.sendall(txt.ljust(MSGLEN).encode("ascii"))

def recv_exact(sock, n):
    buf = b""
    while len(buf) < n:
        chunk = sock.recv(n - len(buf))
        if not chunk:
            raise ConnectionError("Socket closed")
        buf += chunk
    return buf

# Define geometries to run. All coordinates and cell vectors in ANGSTROMS.
geometries_to_run = [
    {  # Water molecule 1
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
    {  # Water molecule 2 (slightly stretched)
        "nat": 3,
        "positions": np.array(
            [[0.0, 0.0, 0.12], [0.0, 0.80, -0.50], [0.0, -0.80, -0.50]]
        ),
        "cell": np.eye(3) * 15.0,
    },
]

HOST, PORT = "127.0.0.1", 9999
with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as server:
    server.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    server.bind((HOST, PORT))
    server.listen(1)
    print(f"Listening on {HOST}:{PORT}")
    conn, addr = server.accept()
    with conn:
        print(f"Connected by {addr}")

        for step, geo in enumerate(geometries_to_run, start=1):
            print(f"\n--- Step {step} ---")

            # 1. Check server status and handle NEEDINIT if necessary.
            send_header(conn, "STATUS")
            reply = recv_exact(conn, MSGLEN).decode().strip()
            print(f"Initial status for step {step}: {reply}")

            if reply == "NEEDINIT":
                print("Server needs reset. Sending INIT.")
                send_header(conn, "INIT")
                conn.sendall(struct.pack("<i", 0))
                conn.sendall(struct.pack("<i", 1))
                conn.sendall(struct.pack("<b", 0))

                # After sending INIT, re-check status to confirm it's READY.
                send_header(conn, "STATUS")
                reply = recv_exact(conn, MSGLEN).decode().strip()
                print(f"Status after INIT: {reply}")

            if reply != "READY":
                print(f"Error: Server is not READY (Status: {reply}). Aborting.")
                break

            # 2. **CRITICAL FIX**: Convert units from Angstrom to Bohr before sending.
            cell_bohr = geo["cell"] / units.Bohr
            invcell_bohr = np.linalg.inv(cell_bohr)
            positions_bohr = geo["positions"] / units.Bohr
            nat = geo["nat"]

            # 3. Send POSDATA with data in correct units and memory order.
            send_header(conn, "POSDATA")
            conn.sendall(struct.pack("<9d", *cell_bohr.flatten(order="F")))
            conn.sendall(struct.pack("<9d", *invcell_bohr.flatten(order="F")))
            conn.sendall(struct.pack("<i", nat))
            conn.sendall(struct.pack(f"<{nat * 3}d", *positions_bohr.flatten(order="F")))
            print(f"Sent POSDATA for {nat} atoms.")

            # 4. Poll until calculation is done.
            while True:
                send_header(conn, "STATUS")
                status = recv_exact(conn, MSGLEN).decode().strip()
                if status == "HAVEDATA":
                    break
                time.sleep(0.01)

            # 5. Request and receive results.
            send_header(conn, "GETFORCE")
            hdr = recv_exact(conn, MSGLEN).decode().strip()
            if hdr != "FORCEREADY":
                print(f"Unexpected header {hdr}, aborting")
                break

            (energy_ha,) = struct.unpack("<d", recv_exact(conn, 8))
            (nat_back,) = struct.unpack("<i", recv_exact(conn, 4))
            forces_flat = struct.unpack(f"<{nat_back*3}d", recv_exact(conn, nat_back*3*8))
            virial_flat = struct.unpack("<9d", recv_exact(conn, 9*8))
            (extra_len,) = struct.unpack("<i", recv_exact(conn, 4))
            extras = recv_exact(conn, extra_len)
            forces_ha_bohr = np.array(forces_flat).reshape((nat_back, 3))

            print(f"Energy (Hartree): {energy_ha}")
            print(f"Forces (Ha/Bohr):\n{forces_ha_bohr}")

        # Send EXIT to close connection cleanly after the loop.
        print("\nAll calculations finished.")
        send_header(conn, "EXIT")
        print("Sent EXIT to NWChem.")
