#include "SocketNWChemPot.h"

#include <cctype>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include <arpa/inet.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <unistd.h>

// Helper to get element symbols from atomic numbers
using std::this_thread::sleep_for;

namespace {
const char *elementArray[] = {
    "Unknown", "H",  "He", "Li", "Be", "B",    "C",  "N",  "O",  "F",  "Ne",
    "Na",      "Mg", "Al", "Si", "P",  "S",    "Cl", "Ar", "K",  "Ca", "Sc",
    "Ti",      "V",  "Cr", "Mn", "Fe", "Co",   "Ni", "Cu", "Zn", "Ga", "Ge",
    "As",      "Se", "Br", "Kr", "Rb", "Sr",   "Y",  "Zr", "Nb", "Mo", "Tc",
    "Ru",      "Rh", "Pd", "Ag", "Cd", "In",   "Sn", "Sb", "Te", "I",  "Xe",
    "Cs",      "Ba", "La", "Ce", "Pr", "Nd",   "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy",      "Ho", "Er", "Tm", "Yb", "Lu",   "Hf", "Ta", "W",  "Re", "Os",
    "Ir",      "Pt", "Au", "Hg", "Tl", "Pb",   "Bi", "Po", "At", "Rn", "Fr",
    "Ra",      "Ac", "Th", "Pa", "U",  nullptr};
char const *atomicNumber2symbol(int n) { return elementArray[n]; }
} // namespace

SocketNWChemPot::SocketNWChemPot(std::shared_ptr<Parameters> p)
    : Potential(PotType::SocketNWChem, p),
      listen_fd(-1),
      conn_fd(-1),
      is_connected(false) {

  unix_socket_mode = p->socket_nwchem_options.unix_socket_mode;
  nwchem_settings = p->socket_nwchem_options.nwchem_settings;
  mem_in_gb = p->socket_nwchem_options.mem_in_gb;
  make_template_input = p->socket_nwchem_options.make_template_input;

  if (unix_socket_mode) {
    unix_socket_basename = p->socket_nwchem_options.unix_socket_path;
    // NWChem client hardcodes this prefix
    server_address = "/tmp/ipi_" + unix_socket_basename;
    port = -1;
    std::cout << "SocketNWChemPot: Initializing in UNIX mode." << std::endl;
    std::cout << "Listening on socket file: " << server_address << std::endl;
  } else {
    server_address = p->socket_nwchem_options.host;
    port = p->socket_nwchem_options.port;
    std::cout << "SocketNWChemPot: Initializing in TCP mode." << std::endl;
    std::cout << "Listening on: " << server_address << ":" << port << std::endl;
  }

  setup_server();
}

SocketNWChemPot::~SocketNWChemPot() {
  if (is_connected) {
    std::cout << "Closing connection to NWChem client..." << std::endl;
    try {
      send_header("EXIT");
    } catch (...) {
      // Ignore errors during shutdown
    }
  }
  if (conn_fd >= 0)
    ::close(conn_fd);
  if (listen_fd >= 0)
    ::close(listen_fd);
  if (unix_socket_mode) {
    ::unlink(server_address.c_str());
  }
}

// =============================================
// Public Methods
// =============================================

void SocketNWChemPot::write_nwchem_template(
    const std::string &filename, long N,
    const std::vector<std::string> &atom_symbols) {
  std::ofstream outfile(filename);
  if (!outfile.is_open()) {
    throw std::runtime_error("Could not open file to write NWChem template: " +
                             filename);
  }

  outfile << "start nwchem_socket_job\n";
  outfile << "title \"NWChem Server for EON\"\n\n";
  outfile << "memory " << mem_in_gb << " gb\n\n";
  outfile << "geometry units bohr noautosym nocenter noautoz\n";
  // This geometry block is only a template for memory allocation.
  // The atom types and count are what matter.
  for (long i = 0; i < N; ++i) {
    outfile << "  " << atom_symbols[i] << "  0.0 0.0 " << static_cast<double>(i)
            << "\n";
  }
  outfile << "end\n\n";
  outfile << "include " << nwchem_settings << "\n\n";
  outfile << "driver\n";
  if (unix_socket_mode) {
    // For the NWChem input, we provide only the basename. NWChem adds the
    // prefix.
    outfile << "  socket unix " << unix_socket_basename << "\n";
  } else {
    outfile << "  socket ipi_client " << server_address << ":" << port << "\n";
  }
  outfile << "end\n\n";
  outfile << "task scf optimize\n";

  outfile.close();
}

void SocketNWChemPot::force(long N, const double *R, const int *atomicNrs,
                            double *F, double *U, double *variance,
                            const double *box) {
  if (!is_connected) {
    std::vector<std::string> symbols;
    symbols.reserve(N);
    for (long i = 0; i < N; ++i) {
      symbols.emplace_back(atomicNumber2symbol(atomicNrs[i]));
    }
    if (make_template_input) {
      write_nwchem_template("nwchem_socket.nwi", N, symbols);
    }

    std::cout << "Waiting for NWChem client connection..." << std::endl;
    accept_connection();
    std::cout << "NWChem client connected." << std::endl;

    // 1. EON acts as the server: after accepting the NWChem client connection,
    // EON sends "STATUS" to the client to query its status.
    char status_buffer[MSG_LEN + 1] = {0};
    send_header("STATUS");

    // 2. EON (acting as server) then waits for NWChem (the client) to respond
    // with "READY".
    recv_header(status_buffer);
    if (std::string(status_buffer) != "READY") {
      throw std::runtime_error(
          "Handshake failed: NWChem client not READY. It sent: " +
          std::string(status_buffer));
    }
    std::cout << "NWChem server is connected and READY." << std::endl;
  }

  // Check status for this specific force call
  char status_buffer[MSG_LEN + 1] = {0};
  send_header("STATUS");
  recv_header(status_buffer);

  if (std::string(status_buffer) == "NEEDINIT") {
    send_header("INIT");
    // Send dummy INIT payload (bead index, number of bytes in extra string)
    int32_t init_payload[] = {0, 1}; // bead_index=0, nbytes=1
    char dummy_byte = 0;
    send_exact(&init_payload, sizeof(init_payload));
    send_exact(&dummy_byte,
               sizeof(dummy_byte)); // No extra string (just a null terminator)
    send_header("STATUS");
    recv_header(status_buffer);
  }

  if (std::string(status_buffer) != "READY") {
    throw std::runtime_error("NWChem server not ready for new positions!");
  }

  // Convert positions to Bohr
  std::vector<double> pos_bohr(N * 3);
  for (size_t i = 0; i < pos_bohr.size(); ++i) {
    pos_bohr[i] = R[i] / BOHR_IN_ANGSTROM;
  }

  // Per i-PI spec, cell and inverse cell must be sent. NWChem does not use
  // this information for non-periodic calculations, so we send an identity
  // matrix as a safe, non-transforming placeholder to conform to the protocol.
  double invcell_T[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  double cell_T[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};

  send_header("POSDATA");
  int32_t nat = N;
  send_exact(cell_T, sizeof(cell_T));
  send_exact(invcell_T, sizeof(invcell_T));
  send_exact(&nat, sizeof(nat));
  send_exact(pos_bohr.data(), pos_bohr.size() * sizeof(double));

  // Poll for results
  while (true) {
    send_header("STATUS");
    recv_header(status_buffer);
    if (std::string(status_buffer) == "HAVEDATA") {
      break;
    }
    // A small sleep to prevent busy-waiting that consumes 100% CPU.
    sleep_for(10ms);
  }

  // Request and receive results ---
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
    throw std::runtime_error("Atom count mismatch from NWChem");
  recv_exact(forces_ha_bohr.data(), forces_ha_bohr.size() * sizeof(double));
  recv_exact(&virial_ha, sizeof(virial_ha));
  recv_exact(&extra_len, sizeof(extra_len));
  if (extra_len > 0) {
    std::vector<char> extra_buf(extra_len);
    recv_exact(extra_buf.data(), extra_len);
  }

  // Convert results back to EON units (eV and Angstrom)
  *U = energy_ha * HARTREE_IN_EV;
  for (int i = 0; i < N * 3; ++i) {
    F[i] = forces_ha_bohr[i] * (HARTREE_IN_EV / BOHR_IN_ANGSTROM);
  }
  if (variance != nullptr) {
    *variance = 0.0;
  }
}

// =================================================
// Private Helper Methods for Socket Communication
// =================================================

void SocketNWChemPot::setup_server() {
  int domain = unix_socket_mode ? AF_UNIX : AF_INET;
  listen_fd = socket(domain, SOCK_STREAM, 0);
  if (listen_fd < 0) {
    throw std::runtime_error("Failed to create socket.");
  }

  if (unix_socket_mode) {
    ::unlink(server_address.c_str()); // Remove stale socket file if it exists
    sockaddr_un sock_addr{};
    sock_addr.sun_family = AF_UNIX;
    strncpy(sock_addr.sun_path, server_address.c_str(),
            sizeof(sock_addr.sun_path) - 1);

    socklen_t addr_len =
        sizeof(sock_addr.sun_family) + strlen(sock_addr.sun_path);
    if (::bind(listen_fd, (struct sockaddr *)&sock_addr, addr_len) < 0) {
      perror("UNIX bind failed");
      throw std::runtime_error("Failed to bind UNIX socket.");
    }
  } else { // TCP Mode
    int opt = 1;
    setsockopt(listen_fd, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt));
    sockaddr_in sock_addr{};
    sock_addr.sin_family = AF_INET;
    sock_addr.sin_addr.s_addr = inet_addr(server_address.c_str());
    sock_addr.sin_port = htons(port);

    if (::bind(listen_fd, (struct sockaddr *)&sock_addr, sizeof(sock_addr)) <
        0) {
      perror("TCP bind failed");
      throw std::runtime_error("Failed to bind TCP socket.");
    }
  }

  if (::listen(listen_fd, 1) < 0) {
    perror("listen() failed");
    throw std::runtime_error("Socket listen() failed.");
  }
}

void SocketNWChemPot::accept_connection() {
  conn_fd = accept(listen_fd, nullptr, nullptr);
  if (conn_fd < 0) {
    perror("accept() failed");
    throw std::runtime_error("Failed to accept client connection.");
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
  for (int i = MSG_LEN - 1; i >= 0 && isspace((unsigned char)buffer[i]); --i) {
    buffer[i] = '\0';
  }
}

void SocketNWChemPot::send_exact(const void *buffer, size_t n_bytes) {
  size_t sent = 0;
  while (sent < n_bytes) {
    ssize_t n = ::send(conn_fd, (const char *)buffer + sent, n_bytes - sent, 0);
    if (n <= 0) {
      throw std::runtime_error(
          "send_exact failed: connection closed or error.");
    }
    sent += n;
  }
}

void SocketNWChemPot::recv_exact(void *buffer, size_t n_bytes) {
  size_t recvd = 0;
  while (recvd < n_bytes) {
    ssize_t n = ::recv(conn_fd, (char *)buffer + recvd, n_bytes - recvd, 0);
    if (n <= 0) {
      throw std::runtime_error(
          "recv_exact failed: connection closed or error.");
    }
    recvd += n;
  }
}
