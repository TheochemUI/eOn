import socket
import struct

import numpy as np

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


HOST, PORT = "127.0.0.1", 9999
with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as server:
    server.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    server.bind((HOST, PORT))
    server.listen(1)
    print(f"Listening on {HOST}:{PORT}")
    conn, addr = server.accept()
    with conn:
        print(f"Connected by {addr}")

        # *** This is the critical change: send the first message ***
        send_header(conn, "STATUS")
        reply = recv_exact(conn, MSGLEN).decode().strip()
        print("Got first reply from NWChem:", reply)

        # Now you can, for example, send POSDATA and then GETFORCE
        if reply == "READY":
            # Example cell & positions
            cell = np.eye(3) * 10.0
            invcell = np.linalg.inv(cell)
            nat = 1
            positions = np.array([[0.71686, 1.24737, -0.16906]])

            send_header(conn, "POSDATA")

            # --- THE FIX: Use order='F' for Fortran array ordering ---
            conn.sendall(struct.pack("<9d", *cell.flatten(order="F")))
            conn.sendall(struct.pack("<9d", *invcell.flatten(order="F")))
            conn.sendall(struct.pack("<i", nat))
            conn.sendall(struct.pack(f"<{nat * 3}d", *positions.flatten(order="F")))

            # Poll until NWChem has HAVEDATA
            while True:
                send_header(conn, "STATUS")
                status = recv_exact(conn, MSGLEN).decode().strip()
                print("Status:", status)
                if status == "HAVEDATA":
                    break

            send_header(conn, "GETFORCE")
            hdr = recv_exact(conn, MSGLEN).decode().strip()
            print("Got header:", hdr)
            if hdr == "FORCEREADY":
                # 1. Read Energy (8-byte double)
                (energy_ha,) = struct.unpack("<d", recv_exact(conn, 8))

                # 2. Read natoms (4-byte integer)
                (nat_back,) = struct.unpack("<i", recv_exact(conn, 4))

                print(f"Energy (Hartree): {energy_ha}")
                print(f"Natoms back: {nat_back}")

                # 3. Read forces (nat * 3 * 8 bytes)
                # Even if nat_back is wrong, we must read what the server sent.
                # The server uses its correct internal `nion=1`.
                forces_flat = struct.unpack(f"<{1 * 3}d", recv_exact(conn, 1 * 3 * 8))
                forces = np.array(forces_flat).reshape((1, 3))

                # 4. MUST READ the virial (9*8 bytes) even though it's zero
                virial_flat = struct.unpack("<9d", recv_exact(conn, 9 * 8))
                virial = np.array(virial_flat).reshape((3, 3), order="F")

                # 5. MUST READ the extra data length (4 bytes) and data (1 byte)
                (extra_len,) = struct.unpack("<i", recv_exact(conn, 4))
                extras = recv_exact(conn, extra_len)

                print(f"Forces (Ha/Bohr):\n{forces}")
                print(f"Virial (Hartree):\n{virial}")
