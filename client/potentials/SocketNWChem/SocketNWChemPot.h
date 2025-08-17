#ifndef SOCKET_NWCHEM_POT_H
#define SOCKET_NWCHEM_POT_H

#include "../../Potential.h"
#include <memory>
#include <string>

/**
 * @brief A potential that gets forces and energy from an NWChem server via a
 * socket.
 *
 * Implements the i-PI socket protocol to communicate with an NWChem executable
 * that has been launched in server mode (e.g., using `task scf optimize` and
 * the `driver socket` directive). This class manages a persistent connection to
 * NWChem. Supports both TCP/IP and UNIX domain sockets.
 */
class SocketNWChemPot : public Potential {
public:
  explicit SocketNWChemPot(std::shared_ptr<Parameters> p);
  ~SocketNWChemPot() override;

  /**
   * @brief Generates an NWChem input file (.nwi) configured to connect to this
   * server.
   * @param filename The path to the output file (e.g., "nwchem_socket.nwi").
   * @param N The number of atoms.
   * @param atom_symbols A vector of atom symbols (e.g., {"O", "H", "H"}).
   */
  void write_nwchem_template(const std::string &filename, long N,
                             const std::vector<std::string> &atom_symbols);

  /**
   * @brief The method called to compute forces and energy.
   *
   * This method implements the full i-PI handshake for a single geometry:
   * 1. Checks the client status, handling the NEEDINIT state if necessary.
   * 2. Sends the new atomic positions and cell (POSDATA).
   * 3. Polls until the calculation is complete (HAVEDATA).
   * 4. Retrieves the results (GETFORCE).
   *
   * @param N The number of atoms.
   * @param R Pointer to the beginning of the position array [x1, y1, z1, ...].
   * @param atomicNrs Pointer to the array of atomic numbers.
   * @param F Pointer to the beginning of the force array to be populated.
   * @param U Pointer to the potential energy value to be set.
   * @param variance Not used by this potential.
   * @param box Pointer to the 3x3 simulation box matrix.
   */
  void force(long N, const double *R, const int *atomicNrs, double *F,
             double *U, double *variance, const double *box) override;

private:
  // --- Private Methods ---
  void setup_server();
  void accept_connection();
  void send_header(const char *msg);
  void recv_header(char *buffer);
  void send_exact(const void *buffer, size_t n_bytes);
  void recv_exact(void *buffer, size_t n_bytes);

  // --- Member Variables ---
  // Socket configuration
  bool unix_socket_mode;
  // For TCP, the host; for UNIX, the full socket path.
  std::string server_address;
  // The original name for the .nwi file.
  std::string unix_socket_basename;
  // File containing actual calculation parameters
  std::string nwchem_settings;
  int port, mem_in_gb;
  // Generate a template input
  bool make_template_input;

  // Socket state
  int listen_fd; // Socket for listening
  int conn_fd;   // Socket for the active connection
  bool is_connected;

  // --- Constants for i-PI protocol and unit conversions ---
  static constexpr int MSG_LEN = 12;
  // From ASE units, values for eV and Angstrom
  static constexpr double BOHR_IN_ANGSTROM = 0.529177210903;
  static constexpr double HARTREE_IN_EV = 27.211386245988;
};

#endif // SOCKET_NWCHEM_POT_H
