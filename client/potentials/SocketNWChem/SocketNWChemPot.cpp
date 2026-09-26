#include "eon/potentials/SocketNWChem/SocketNWChemPot.h"
#include "eon/EonLogger.h"
#include "eon/Parameters.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <format>
#include <fstream>
#include <stdexcept>

#include <arpa/inet.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <unistd.h>

#include <chrono>
#include <cstdint>
#include <readcon-core.hpp>
#include <thread>

SocketNWChemPot::SocketNWChemPot(const eonc::Parameters &p)
    : eonc::Potential(eonc::PotType::SocketNWChem, p) {

  unix_socket_mode = p.socket_nwchem_options().unix_socket_mode;
  nwchem_settings = p.socket_nwchem_options().nwchem_settings;
  mem_in_gb = p.socket_nwchem_options().mem_in_gb;
  make_template_input = p.socket_nwchem_options().make_template_input;

  if (unix_socket_mode) {
    unix_socket_basename = p.socket_nwchem_options().unix_socket_path;
    // NWChem's Fortran i-PI driver truncates the socket name to ~30 chars.
    // The full path is /tmp/ipi_<basename>, so basename must be short.
    server_address = "/tmp/ipi_" + unix_socket_basename;
    if (server_address.size() > 30) {
      EONC_LOG_ERROR(
          "UNIX socket path '{}' is {} characters, past NWChem's ~30 character "
          "limit. Shorten unix_socket_path (currently '{}') to at most {} "
          "characters.",
          server_address, server_address.size(), unix_socket_basename, 30 - 9);
      throw std::runtime_error(
          "unix_socket_path too long for NWChem (max ~21 chars, got " +
          std::to_string(unix_socket_basename.size()) + ")");
    }
    port = -1;
    EONC_LOG_INFO("SocketNWChemPot UNIX socket {}", server_address);
  } else {
    server_address = p.socket_nwchem_options().host;
    port = p.socket_nwchem_options().port;
    EONC_LOG_INFO("SocketNWChemPot TCP {}:{}", server_address, port);
  }

  setup_server();
}

SocketNWChemPot::~SocketNWChemPot() {
  if (is_connected) {
    EONC_LOG_INFO("Closing connection to NWChem client");
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
  outfile << "title \"NWChem Server for eOn\"\n\n";
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
  try {
    forceOnce(N, R, atomicNrs, F, U, variance, box);
    return;
  } catch (const std::runtime_error &) {
    drop_connection();
  }
  forceOnce(N, R, atomicNrs, F, U, variance, box);
}

void SocketNWChemPot::forceOnce(long N, const double *R, const int *atomicNrs,
                                double *F, double *U, double *variance,
                                const double *box) {
  if (!is_connected) {
    std::vector<std::string> symbols;
    symbols.reserve(N);
    for (long i = 0; i < N; ++i) {
      const int z = atomicNrs[i];
      symbols.emplace_back(
          z > 0 ? readcon::z_to_symbol(static_cast<uint64_t>(z)) : "X");
    }
    if (make_template_input) {
      write_nwchem_template("nwchem_socket.nwi", N, symbols);
    }

    EONC_LOG_INFO("Waiting for NWChem client connection");
    accept_connection();
    EONC_LOG_INFO("NWChem client connected");

    // 1. eOn acts as the server: after accepting the NWChem client connection,
    // eOn sends "STATUS" to the client to query its status.
    std::array<char, MSG_LEN + 1> status_buffer{};
    send_header("STATUS");

    // 2. eOn (acting as server) then waits for NWChem (the client) to respond
    // with "READY".
    recv_header(status_buffer.data());
    if (std::string(status_buffer.data()) != "READY") {
      throw std::runtime_error(
          "Handshake failed: NWChem client not READY. It sent: " +
          std::string(status_buffer.data()));
    }
    EONC_LOG_INFO("NWChem server is connected and READY");
  }

  // Check status for this specific force call
  std::array<char, MSG_LEN + 1> status_buffer{};
  send_header("STATUS");
  recv_header(status_buffer.data());

  if (std::string(status_buffer.data()) == "NEEDINIT") {
    send_header("INIT");
    // Send dummy INIT payload (bead index, number of bytes in extra string)
    std::array<int32_t, 2> init_payload{{0, 1}}; // bead_index=0, nbytes=1
    char dummy_byte = 0;
    send_exact(init_payload.data(), init_payload.size() * sizeof(int32_t));
    send_exact(&dummy_byte,
               sizeof(dummy_byte)); // No extra string (just a null terminator)
    send_header("STATUS");
    recv_header(status_buffer.data());
  }

  if (std::string(status_buffer.data()) != "READY") {
    throw std::runtime_error("NWChem server not ready for new positions!");
  }

  // Convert positions to Bohr
  std::vector<double> pos_bohr(N * 3);
  for (size_t i = 0; i < pos_bohr.size(); ++i) {
    pos_bohr[i] = R[i] / BOHR_IN_ANGSTROM;
  }

  // Per i-PI spec, cell and inverse cell must be sent. NWChem does not use
  // this information for non-periodic calculations, so the frame carries an
  // identity cell and its inverse.
  std::array<double, 9> invcell_T{{1, 0, 0, 0, 1, 0, 0, 0, 1}};
  std::array<double, 9> cell_T{{1, 0, 0, 0, 1, 0, 0, 0, 1}};

  send_header("POSDATA");
  int32_t nat = static_cast<int32_t>(N);
  send_exact(cell_T.data(), cell_T.size() * sizeof(double));
  send_exact(invcell_T.data(), invcell_T.size() * sizeof(double));
  send_exact(&nat, sizeof(nat));
  send_exact(pos_bohr.data(), pos_bohr.size() * sizeof(double));

  // Poll for results
  while (true) {
    send_header("STATUS");
    recv_header(status_buffer.data());
    if (std::string(status_buffer.data()) == "HAVEDATA") {
      break;
    }
    // A small sleep to prevent busy-waiting that consumes 100% CPU.
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
  }

  // Request and receive results ---
  send_header("GETFORCE");
  recv_header(status_buffer.data());
  if (std::string(status_buffer.data()) != "FORCEREADY") {
    throw std::runtime_error("Expected FORCEREADY, got " +
                             std::string(status_buffer.data()));
  }

  // Unpack the results payload.
  double energy_ha;
  int32_t nat_back;
  std::vector<double> forces_ha_bohr(N * 3);
  std::array<double, 9> virial_ha{};
  int32_t extra_len;

  recv_exact(&energy_ha, sizeof(energy_ha));
  recv_exact(&nat_back, sizeof(nat_back));
  if (nat_back != N)
    throw std::runtime_error("Atom count mismatch from NWChem");
  recv_exact(forces_ha_bohr.data(), forces_ha_bohr.size() * sizeof(double));
  recv_exact(virial_ha.data(), virial_ha.size() * sizeof(double));
  recv_exact(&extra_len, sizeof(extra_len));
  if (extra_len > 0) {
    std::vector<char> extra_buf(extra_len);
    recv_exact(extra_buf.data(), extra_len);
  }

  // Convert results back to eOn units (eV and Angstrom)
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

namespace {

[[noreturn]] void socket_fail(int fd, const char *what) {
  const int err = errno;
  if (fd >= 0) {
    ::close(fd);
  }
  throw std::runtime_error(std::format("{}: {}", what, std::strerror(err)));
}

} // namespace

void SocketNWChemPot::setup_server() {
  int domain = unix_socket_mode ? AF_UNIX : AF_INET;
  listen_fd = socket(domain, SOCK_STREAM, 0);
  if (listen_fd < 0) {
    throw std::runtime_error(
        std::format("Failed to create socket: {}", std::strerror(errno)));
  }

  if (unix_socket_mode) {
    ::unlink(server_address.c_str()); // Remove stale socket file if it exists
    sockaddr_un sock_addr{};
    sock_addr.sun_family = AF_UNIX;
    std::strncpy(sock_addr.sun_path, server_address.c_str(),
                 sizeof(sock_addr.sun_path) - 1);

    socklen_t addr_len = static_cast<socklen_t>(
        sizeof(sock_addr.sun_family) + std::strlen(sock_addr.sun_path));
    if (::bind(listen_fd, reinterpret_cast<sockaddr *>(&sock_addr), addr_len) <
        0) {
      const int fd = listen_fd;
      listen_fd = -1;
      socket_fail(fd, "Failed to bind UNIX socket");
    }
  } else {
    int opt = 1;
    if (setsockopt(listen_fd, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt)) <
        0) {
      const int fd = listen_fd;
      listen_fd = -1;
      socket_fail(fd, "Failed to set SO_REUSEADDR");
    }
    sockaddr_in sock_addr{};
    sock_addr.sin_family = AF_INET;
    sock_addr.sin_addr.s_addr = inet_addr(server_address.c_str());
    // inet_addr accepts only a dotted IPv4 address. A name fails here instead
    // of being bound as INADDR_NONE.
    if (sock_addr.sin_addr.s_addr == INADDR_NONE) {
      ::close(listen_fd);
      listen_fd = -1;
      throw std::runtime_error(
          "Failed to parse the TCP host as an IPv4 address");
    }
    sock_addr.sin_port = htons(static_cast<uint16_t>(port));

    if (::bind(listen_fd, reinterpret_cast<sockaddr *>(&sock_addr),
               sizeof(sock_addr)) < 0) {
      const int fd = listen_fd;
      listen_fd = -1;
      socket_fail(fd, "Failed to bind TCP socket");
    }
  }

  if (::listen(listen_fd, 1) < 0) {
    const int fd = listen_fd;
    listen_fd = -1;
    socket_fail(fd, "Socket listen() failed");
  }
}

void SocketNWChemPot::accept_connection() {
  conn_fd = ::accept(listen_fd, nullptr, nullptr);
  if (conn_fd < 0) {
    throw std::runtime_error(std::format(
        "Failed to accept client connection: {}", std::strerror(errno)));
  }
  is_connected = true;
}

void SocketNWChemPot::drop_connection() {
  if (conn_fd >= 0) {
    ::close(conn_fd);
    conn_fd = -1;
  }
  is_connected = false;
}

void SocketNWChemPot::send_header(const char *msg) {
  std::array<char, MSG_LEN> buffer{};
  const std::size_t n = std::min(std::strlen(msg), buffer.size());
  std::copy_n(msg, n, buffer.begin());
  send_exact(buffer.data(), buffer.size());
}

void SocketNWChemPot::recv_header(char *buffer) {
  recv_exact(buffer, MSG_LEN);
  buffer[MSG_LEN] = '\0'; // Null-terminate
  // Trim trailing whitespace
  for (int i = MSG_LEN - 1;
       i >= 0 && std::isspace(static_cast<unsigned char>(buffer[i])); --i) {
    buffer[i] = '\0';
  }
}

void SocketNWChemPot::send_exact(const void *buffer, size_t n_bytes) {
  size_t sent = 0;
  while (sent < n_bytes) {
    ssize_t n = ::send(conn_fd, static_cast<const char *>(buffer) + sent,
                       n_bytes - sent, 0);
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
    ssize_t n = ::recv(conn_fd, static_cast<char *>(buffer) + recvd,
                       n_bytes - recvd, 0);
    if (n <= 0) {
      throw std::runtime_error(
          "recv_exact failed: connection closed or error.");
    }
    recvd += n;
  }
}
