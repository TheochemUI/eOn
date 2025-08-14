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

        # Initial setup
        cell = np.eye(3) * 10.0
        invcell = np.linalg.inv(cell)
        nat = 1

        # List of positions for multiple geometry steps
        positions_list = [
            np.array([[0.71686, 1.24737, -0.16906]]),
            np.array([[0.72686, 1.24737, -0.16906]]),  # perturbed x coordinate
            np.array([[0.73686, 1.24737, -0.16906]]),
        ]

        for step, positions in enumerate(positions_list, start=1):
            print(f"\n--- Step {step} ---")

            # Send STATUS first to initiate
            send_header(conn, "STATUS")
            reply = recv_exact(conn, MSGLEN).decode().strip()
            if reply == "NEEDINIT":
                print("Server in NEEDINIT state. Sending INIT to reset.")
                send_header(conn, "INIT")
                # The INIT command requires a small dummy payload
                conn.sendall(struct.pack("<i", 0))  # bead index
                conn.sendall(struct.pack("<i", 1))  # nbytes of extra data
                conn.sendall(struct.pack("<b", 0))  # a single dummy byte

                # After INIT, the server must be READY. Let's confirm.
                send_header(conn, "STATUS")
                reply = recv_exact(conn, MSGLEN).decode().strip()
                print(f"Status after INIT: {reply}")

            print(f"Got reply to STATUS: {reply}")
            if reply != "READY":
                print("Unexpected reply, aborting")
                break

            # Send POSDATA with updated positions
            send_header(conn, "POSDATA")
            conn.sendall(struct.pack("<9d", *cell.flatten(order="F")))
            conn.sendall(struct.pack("<9d", *invcell.flatten(order="F")))
            conn.sendall(struct.pack("<i", nat))
            conn.sendall(struct.pack(f"<{nat * 3}d", *positions.flatten(order="F")))
            print(f"Sent POSDATA for {nat} atoms.")

            # Poll for HAVEDATA status indicating calculation finished
            while True:
                send_header(conn, "STATUS")
                status = recv_exact(conn, MSGLEN).decode().strip()
                print("Status:", status)
                if status == "HAVEDATA":
                    break

            # Request forces and energy
            send_header(conn, "GETFORCE")
            hdr = recv_exact(conn, MSGLEN).decode().strip()
            print("Got header:", hdr)
            if hdr != "FORCEREADY":
                print("Unexpected header, aborting")
                break

            (energy_ha,) = struct.unpack("<d", recv_exact(conn, 8))
            (nat_back,) = struct.unpack("<i", recv_exact(conn, 4))
            forces_flat = struct.unpack(
                f"<{nat_back * 3}d", recv_exact(conn, nat_back * 3 * 8)
            )
            virial_flat = struct.unpack("<9d", recv_exact(conn, 9 * 8))
            (extra_len,) = struct.unpack("<i", recv_exact(conn, 4))
            extras = recv_exact(conn, extra_len)

            forces = np.array(forces_flat).reshape((nat_back, 3))
            virial = np.array(virial_flat).reshape((3, 3), order="F")

            print(f"Energy (Hartree): {energy_ha}")
            print(f"Forces (Ha/Bohr):\n{forces}")
            print(f"Virial (Hartree):\n{virial}")
            print(f"Extras bytes length: {len(extras)}")

        # Optionally, send EXIT to close connection cleanly
        send_header(conn, "EXIT")
        print("Sent EXIT to NWChem.")
