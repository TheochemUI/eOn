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

namespace {

const char *elementArray[] = {
    "Unknown", "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na",
    "Mg",      "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",
    "Cr",      "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
    "Kr",      "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
    "Cd",      "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd",      "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf",      "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
    "Po",      "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  NULL};

// guess the atom type from the atomic mass,
std::string mass2atom(double atomicmass) {
  return elementArray[int(atomicmass + .5)];
}

int symbol2atomicNumber(char const *symbol) {
  int i = 0;

  while (elementArray[i] != NULL) {
    if (strcmp(symbol, elementArray[i]) == 0) {
      return i;
    }
    i++;
  }
  // invalid symbol
  return -1;
}

char const *atomicNumber2symbol(int n) { return elementArray[n]; }
} // namespace

void SocketNWChemPot::write_nwchem_template(
    const std::string &filename, long N,
    const std::vector<std::string> &atom_symbols) {
  if (atom_symbols.size() != static_cast<size_t>(N)) {
    throw std::invalid_argument(
        "Number of atoms (N) does not match size of atom_symbols vector.");
  }

  std::cout << "Generating NWChem input template: " << filename << std::endl;
  std::ofstream outfile(filename);

  outfile << "start nwchem_socket_job\n";
  outfile << "title \"NWChem Server for EON\"\n\n";
  outfile << "memory 2 gb\n\n";
  outfile << "geometry units angstroms noautosym nocenter noautoz\n";
  // This geometry block is only a template for memory allocation.
  // The atom types and count are what matter.
  for (long i = 0; i < N; ++i) {
    outfile << "  " << atom_symbols[i] << "  0.0 0.0 " << static_cast<double>(i)
            << "\n";
  }
  outfile << "end\n\n";
  outfile << "basis noprint\n";
  outfile << "  * library 3-21G\n";
  outfile << "end\n\n";
  outfile << "scf\n";
  outfile << "  nopen 0\n";
  outfile << "  thresh 1e-8\n";
  outfile << "  maxiter 2000\n";
  outfile << "end\n\n";
  outfile << "driver\n";
  outfile << "  socket ipi_client " << this->host << ":" << this->port << "\n";
  outfile << "end\n\n";
  outfile << "task scf optimize\n";

  outfile.close();
}

// Constructor: Read parameters, connect, and perform initial handshake.
SocketNWChemPot::SocketNWChemPot(std::shared_ptr<Parameters> p)
    : Potential(PotType::SocketNWChem, p),
      listen_fd(-1),
      conn_fd(-1),
      is_listening(false),
      is_connected(false) {
  // Read connection parameters from the params file.
  host = p->socket_nwchem_options.host;
  port = p->socket_nwchem_options.port;

  std::cout << "Initializing SocketNWChemPot..." << std::endl;
  setup_server();
}

// Destructor: Send EXIT and close the socket.
SocketNWChemPot::~SocketNWChemPot() {
  if (is_connected) {
    std::cout << "Closing connection to NWChem server..." << std::endl;
    try {
      send_header("EXIT");
    } catch (...) {
    }
    close(conn_fd);
  }
  if (is_listening) {
    std::cout << "Closing connection to whatever is listening..." << std::endl;
    close(listen_fd);
  }
}

// The main force call method.
void SocketNWChemPot::force(long N, const double *R, const int *atomicNrs,
                            double *F, double *U, double *variance,
                            const double *box) {
  // On the very first force call, block and wait for NWChem to connect.
  if (!is_connected) {
    std::vector<std::string> symbols;
    symbols.reserve(N);
    for (long i = 0; i < N; ++i) {
      symbols.emplace_back(atomicNumber2symbol(atomicNrs[i]));
    }

    // Write the NWChem input template to disk
    write_nwchem_template("nwchem_socket.nwi", N, symbols);
    accept_connection();

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
  // 1. Create vectors to hold data in atomic units (Bohr).
  std::vector<double> pos_bohr(N * 3);
  std::vector<double> cell_bohr(9);

  // 2. Copy data from input pointers and convert from Angstrom to Bohr.
  for (int i = 0; i < N * 3; ++i) {
    pos_bohr[i] = R[i] / BOHR_IN_ANGSTROM;
  }
  for (int i = 0; i < 9; ++i) {
    cell_bohr[i] = box[i] / BOHR_IN_ANGSTROM;
  }

  // 3. Calculate the inverse cell. A proper 3x3 inverse is needed for
  // non-orthogonal cells.
  std::vector<double> invcell_bohr(9, 0.0);
  if (cell_bohr[0] != 0)
    invcell_bohr[0] = 1.0 / cell_bohr[0];
  if (cell_bohr[4] != 0)
    invcell_bohr[4] = 1.0 / cell_bohr[4];
  if (cell_bohr[8] != 0)
    invcell_bohr[8] = 1.0 / cell_bohr[8];

  // 4. Transpose matrices to Fortran memory order for sending.
  double cell_T[9], invcell_T[9];
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      cell_T[j * 3 + i] = cell_bohr[i * 3 + j];
      invcell_T[j * 3 + i] = invcell_bohr[i * 3 + j];
    }

  // --- 3. Send POSDATA message ---
  // Even empty things like cell and inverse cell MUST be present
  send_header("POSDATA");
  int32_t nat = N;
  send_exact(cell_T, sizeof(cell_T));
  send_exact(invcell_T, sizeof(invcell_T));
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
    F[i] = -1 * forces_ha_bohr[i] * (HARTREE_IN_EV / BOHR_IN_ANGSTROM);
  }
  // Variance is not calculated by NWChem.
  *variance = 0.0;
}

// ===================================================================
// Private Helper Methods for Socket Communication
// ===================================================================

// New private method to set up the listening socket.
void SocketNWChemPot::setup_server() {
  listen_fd = socket(AF_INET, SOCK_STREAM, 0);
  if (listen_fd < 0)
    throw std::runtime_error("Failed to create server socket.");

  // Allow address reuse to prevent "address already in use" errors on restart.
  int opt = 1;
  setsockopt(listen_fd, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt));

  sockaddr_in serv_addr;
  serv_addr.sin_family = AF_INET;
  serv_addr.sin_addr.s_addr = inet_addr(host.c_str());
  serv_addr.sin_port = htons(port);

  if (bind(listen_fd, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) {
    throw std::runtime_error("Failed to bind server socket to port " +
                             std::to_string(port));
  }

  if (listen(listen_fd, 1) < 0) {
    throw std::runtime_error("Failed to listen on server socket.");
  }
  is_listening = true;
  std::cout << "SocketNWChemPot is listening for NWChem on " << host << ":"
            << port << std::endl;
}

// New private method to block and wait for a connection.
void SocketNWChemPot::accept_connection() {
  std::cout << "Waiting for NWChem client to connect..." << std::endl;
  sockaddr_in client_addr;
  socklen_t client_len = sizeof(client_addr);
  conn_fd = accept(listen_fd, (struct sockaddr *)&client_addr, &client_len);

  if (conn_fd < 0) {
    throw std::runtime_error("Failed to accept client connection.");
  }
  is_connected = true;
  char client_ip[INET_ADDRSTRLEN];
  inet_ntop(AF_INET, &client_addr.sin_addr, client_ip, INET_ADDRSTRLEN);
  std::cout << "Accepted connection from NWChem client: " << client_ip
            << std::endl;
}

void SocketNWChemPot::send_header(const char *msg) {
  char buffer[MSG_LEN] = {0};
  strncpy(buffer, msg, MSG_LEN);
  send_exact(buffer, MSG_LEN);
}

void SocketNWChemPot::recv_header(char *buffer) {
  recv_exact(buffer, MSG_LEN);
  buffer[MSG_LEN] = '\0';
  for (int i = MSG_LEN - 1; i >= 0 && isspace(buffer[i]); --i)
    buffer[i] = '\0';
}

void SocketNWChemPot::send_exact(const void *buffer, size_t n_bytes) {
  size_t sent = 0;
  while (sent < n_bytes) {
    ssize_t n = send(conn_fd, (const char *)buffer + sent, n_bytes - sent, 0);
    if (n <= 0) {
      // An error or a closed connection occurred.
      throw std::runtime_error(
          "Failed to send bytes to socket or connection closed.");
    }
    sent += n;
  }
}

void SocketNWChemPot::recv_exact(void *buffer, size_t n_bytes) {
  size_t received = 0;
  while (received < n_bytes) {
    ssize_t n = recv(conn_fd, (char *)buffer + received, n_bytes - received, 0);
    if (n <= 0)
      throw std::runtime_error("Failed to receive bytes or connection closed.");
    received += n;
  }
}
