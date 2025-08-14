#include "SocketNWChemPot.h"
#include "../../Parameters.h"
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>

// Includes for POSIX sockets
#include <arpa/inet.h>
#include <netdb.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <unistd.h>

// Constructor: Read parameters, connect, and perform initial handshake.
SocketNWChemPot::SocketNWChemPot(std::shared_ptr<Parameters> p)
    : Potential(p),
      is_connected(false) {
  // Read connection parameters from the params file.
  host = p->socket_nwchem_options.host;
  port = p->socket_nwchem_options.port;

  std::cout << "Initializing SocketNWChemPot..." << std::endl;
  connect_to_server();

  // Perform the first handshake to ensure the server is ready.
  char status_buffer[MSG_LEN + 1] = {0};
  send_header("STATUS");
  recv_header(status_buffer);
  if (std::string(status_buffer) != "READY") {
    throw std::runtime_error(
        "NWChem server is not in READY state on connect. Status: " +
        std::string(status_buffer));
  }
  std::cout << "NWChem server is connected and READY." << std::endl;
}

// Destructor: Send EXIT and close the socket.
SocketNWChemPot::~SocketNWChemPot() {
  if (is_connected) {
    std::cout << "Closing connection to NWChem server..." << std::endl;
    try {
      send_header("EXIT");
    } catch (const std::exception &e) {
      // Ignore errors on close, as the server might already be gone.
    }
    close(sock_fd);
  }
}

// The main force call method.
void SocketNWChemPot::force(long N, const double *R, const int *atomicNrs,
                            double *F, double *U, double *variance,
                            const double *box) {
  // --- 1. Check server status and handle NEEDINIT if necessary. ---
  char status_buffer[MSG_LEN + 1] = {0};
  send_header("STATUS");
  recv_header(status_buffer);

  if (std::string(status_buffer) == "NEEDINIT") {
    send_header("INIT");
    // Send dummy INIT payload (bead index, 1 extra byte, a zero byte)
    int32_t init_payload[] = {0, 1};
    char dummy_byte = 0;
    send_exact(&init_payload, sizeof(init_payload));
    send_exact(&dummy_byte, sizeof(dummy_byte));

    // Re-check status to confirm it's now READY.
    send_header("STATUS");
    recv_header(status_buffer);
  }

  if (std::string(status_buffer) != "READY") {
    throw std::runtime_error("NWChem server is not ready! Status: " +
                             std::string(status_buffer));
  }

  // --- 2. Convert units and prepare data for sending. ---
  std::vector<double> pos_bohr(N * 3);
  // The cell is not really used, so it doesn't matter
  std::array<double, 9> cell_bohr;
  std::fill(cell_bohr.begin(), cell_bohr.end(), 0);
  // --- 3. Send POSDATA message ---
  // Even empty things like cell and inverse cell MUST be present
  send_header("POSDATA");
  int32_t nat = N;
  send_exact(cell_bohr.data(), sizeof(cell_bohr));
  send_exact(/*inverse now*/ cell_bohr.data(), sizeof(cell_bohr));
  send_exact(&nat, sizeof(nat));
  send_exact(pos_bohr.data(), pos_bohr.size() * sizeof(double));

  // --- 4. Poll with STATUS until calculation is done ---
  while (true) {
    send_header("STATUS");
    recv_header(status_buffer);
    if (std::string(status_buffer) == "HAVEDATA") {
      break;
    }
    // A small sleep to prevent busy-waiting that consumes 100% CPU.
    usleep(10000); // 10ms
  }

  // --- 5. Request and receive results ---
  send_header("GETFORCE");
  recv_header(status_buffer);
  if (std::string(status_buffer) != "FORCEREADY") {
    throw std::runtime_error("Expected FORCEREADY, got " +
                             std::string(status_buffer));
  }

  // Unpack the results payload.
  double energy_ha;
  int32_t nat_back;
  std::vector<double> forces_ha_bohr(N * 3);
  double virial_ha[9];
  int32_t extra_len;

  recv_exact(&energy_ha, sizeof(energy_ha));
  recv_exact(&nat_back, sizeof(nat_back));
  if (nat_back != N)
    throw std::runtime_error("Atom count mismatch from server!");
  recv_exact(forces_ha_bohr.data(), forces_ha_bohr.size() * sizeof(double));
  recv_exact(&virial_ha, sizeof(virial_ha));
  recv_exact(&extra_len, sizeof(extra_len));
  if (extra_len > 0) {
    std::vector<char> extra_buf(extra_len);
    recv_exact(extra_buf.data(), extra_len);
  }

  // --- 6. Populate output arrays with unit conversions. ---
  // The main code likely expects eV and eV/Angstrom.
  *U = energy_ha * HARTREE_IN_EV;
  for (int i = 0; i < N * 3; ++i) {
    F[i] = forces_ha_bohr[i] * (HARTREE_IN_EV / BOHR_IN_ANGSTROM);
  }
  // Variance is not calculated by NWChem.
  *variance = 0.0;
}

// ===================================================================
// Private Helper Methods for Socket Communication
// ===================================================================

void SocketNWChemPot::connect_to_server() {
  sock_fd = socket(AF_INET, SOCK_STREAM, 0);
  if (sock_fd < 0) {
    throw std::runtime_error("Failed to create socket.");
  }

  sockaddr_in serv_addr;
  serv_addr.sin_family = AF_INET;
  serv_addr.sin_port = htons(port);

  if (inet_pton(AF_INET, host.c_str(), &serv_addr.sin_addr) <= 0) {
    throw std::runtime_error("Invalid address or address not supported.");
  }

  if (connect(sock_fd, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) {
    throw std::runtime_error("Connection to NWChem server failed.");
  }

  is_connected = true;
}

void SocketNWChemPot::send_header(const char *msg) {
  char buffer[MSG_LEN] = {0};
  strncpy(buffer, msg, MSG_LEN);
  send_exact(buffer, MSG_LEN);
}

void SocketNWChemPot::recv_header(char *buffer) {
  recv_exact(buffer, MSG_LEN);
  buffer[MSG_LEN] = '\0'; // Null-terminate
  // Trim trailing whitespace
  for (int i = MSG_LEN - 1; i >= 0 && isspace(buffer[i]); --i) {
    buffer[i] = '\0';
  }
}

void SocketNWChemPot::send_exact(const void *buffer, size_t n_bytes) {
  if (send(sock_fd, buffer, n_bytes, 0) != (ssize_t)n_bytes) {
    throw std::runtime_error("Failed to send all bytes to socket.");
  }
}

void SocketNWChemPot::recv_exact(void *buffer, size_t n_bytes) {
  size_t received = 0;
  while (received < n_bytes) {
    ssize_t n = recv(sock_fd, (char *)buffer + received, n_bytes - received, 0);
    if (n <= 0) {
      throw std::runtime_error(
          "Failed to receive bytes from socket or connection closed.");
    }
    received += n;
  }
}
