#include "catch2/catch_amalgamated.hpp"
#include "eon/Parameters.h"
#include "eon/potentials/SocketNWChem/SocketNWChemPot.h"

#include <atomic>
#include <cctype>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <optional>
#include <sstream>
#include <string>
#include <thread>

#include <arpa/inet.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <unistd.h>

TEST_CASE("SocketNWChem listens and writes a template", "[socketnwchem]") {
  eonc::Parameters params;
  eonc::ParametersLoadAccess::potential_options(params).potential =
      eonc::PotType::SocketNWChem;
  auto &opt = eonc::ParametersLoadAccess::socket_nwchem_options(params);
  opt.unix_socket_mode = true;
  opt.unix_socket_path = "eon_sock_t";
  opt.mem_in_gb = 1;
  opt.nwchem_settings = "theory.nw";

  SocketNWChemPot pot(params);
  const auto path =
      std::filesystem::temp_directory_path() / "eon-nwchem-template.nw";
  pot.write_nwchem_template(path.string(), 2, {"H", "O"});
  std::ifstream in(path);
  std::stringstream text;
  text << in.rdbuf();
  REQUIRE(text.str().find("memory 1 gb") != std::string::npos);
  REQUIRE(text.str().find("geometry units bohr") != std::string::npos);
  std::filesystem::remove(path);
}

namespace {

bool write_all(int fd, const void *data, size_t n) {
  const auto *bytes = static_cast<const char *>(data);
  size_t sent = 0;
  while (sent < n) {
    const ssize_t got = ::send(fd, bytes + sent, n - sent, 0);
    if (got <= 0) {
      return false;
    }
    sent += static_cast<size_t>(got);
  }
  return true;
}

bool read_all(int fd, void *data, size_t n) {
  auto *bytes = static_cast<char *>(data);
  size_t got_n = 0;
  while (got_n < n) {
    const ssize_t got = ::recv(fd, bytes + got_n, n - got_n, 0);
    if (got <= 0) {
      return false;
    }
    got_n += static_cast<size_t>(got);
  }
  return true;
}

bool write_header(int fd, const char *msg) {
  char buffer[12];
  std::memset(buffer, ' ', sizeof(buffer));
  const size_t n = std::min(std::strlen(msg), sizeof(buffer));
  std::memcpy(buffer, msg, n);
  return write_all(fd, buffer, sizeof(buffer));
}

std::string read_header(int fd) {
  char buffer[12];
  if (!read_all(fd, buffer, sizeof(buffer))) {
    return {};
  }
  // The server zero-fills a header. A peer may space-pad it instead.
  std::size_t n = 0;
  while (n < sizeof(buffer) && buffer[n] != '\0') {
    ++n;
  }
  while (n > 0 && std::isspace(static_cast<unsigned char>(buffer[n - 1]))) {
    --n;
  }
  return std::string(buffer, n);
}

int connect_unix(const std::string &path) {
  sockaddr_un addr{};
  addr.sun_family = AF_UNIX;
  std::strncpy(addr.sun_path, path.c_str(), sizeof(addr.sun_path) - 1);
  const socklen_t len = static_cast<socklen_t>(sizeof(addr.sun_family) +
                                               std::strlen(addr.sun_path));
  for (int attempt = 0; attempt < 200; ++attempt) {
    const int fd = ::socket(AF_UNIX, SOCK_STREAM, 0);
    if (fd < 0) {
      return -1;
    }
    if (::connect(fd, reinterpret_cast<sockaddr *>(&addr), len) == 0) {
      return fd;
    }
    ::close(fd);
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
  }
  return -1;
}

void write_force(int fd) {
  const double energy = 1.5;
  const int32_t nat = 2;
  const double forces[6] = {0.2, 0.0, 0.0, 0.0, 0.0, 0.0};
  double virial[9] = {};
  virial[0] = 1.0;
  const int32_t extra = 0;
  write_header(fd, "FORCEREADY");
  write_all(fd, &energy, sizeof(energy));
  write_all(fd, &nat, sizeof(nat));
  write_all(fd, forces, sizeof(forces));
  write_all(fd, virial, sizeof(virial));
  write_all(fd, &extra, sizeof(extra));
}

} // namespace

TEST_CASE("SocketNWChem exchanges one i-PI force", "[socketnwchem]") {
  eonc::Parameters params;
  eonc::ParametersLoadAccess::potential_options(params).potential =
      eonc::PotType::SocketNWChem;
  auto &opt = eonc::ParametersLoadAccess::socket_nwchem_options(params);
  opt.unix_socket_mode = true;
  opt.unix_socket_path = "eonqf";
  opt.mem_in_gb = 1;
  opt.nwchem_settings = "theory.nw";
  opt.make_template_input = true;

  const auto dir =
      std::filesystem::temp_directory_path() /
      ("eon-ipi-" + std::to_string(static_cast<long long>(::getpid())));
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);
  const auto previous = std::filesystem::current_path();
  std::filesystem::current_path(dir);

  std::atomic<bool> client_ok{true};
  std::optional<SocketNWChemPot> pot;
  pot.emplace(params);
  std::thread client([&] {
    const std::string path = "/tmp/ipi_eonqf";
    const int fd = connect_unix(path);
    if (fd < 0) {
      client_ok = false;
      return;
    }
    int phase = 1;
    while (phase < 13) {
      const std::string header = read_header(fd);
      if (header.empty()) {
        client_ok = false;
        break;
      }
      if (header == "EXIT") {
        break;
      }
      if (header == "STATUS" && phase == 1) {
        write_header(fd, "READY");
        phase = 2;
      } else if (header == "STATUS" && phase == 2) {
        write_header(fd, "NEEDINIT");
        phase = 3;
      } else if (header == "INIT" && phase == 3) {
        char extra[9];
        if (!read_all(fd, extra, sizeof(extra))) {
          client_ok = false;
          break;
        }
        phase = 4;
      } else if (header == "STATUS" && phase == 4) {
        write_header(fd, "READY");
        phase = 5;
      } else if (header == "POSDATA" && (phase == 5 || phase == 9)) {
        char payload[196];
        if (!read_all(fd, payload, sizeof(payload))) {
          client_ok = false;
          break;
        }
        phase = phase == 5 ? 6 : 10;
      } else if (header == "STATUS" && (phase == 6 || phase == 10)) {
        write_header(fd, "HAVEDATA");
        phase = phase == 6 ? 7 : 11;
      } else if (header == "GETFORCE" && (phase == 7 || phase == 11)) {
        write_force(fd);
        phase = phase == 7 ? 8 : 12;
      } else if (header == "STATUS" && phase == 8) {
        write_header(fd, "READY");
        phase = 9;
      } else {
        client_ok = false;
        break;
      }
    }
    ::close(fd);
  });

  const double positions[6] = {0.0, 0.0, 0.0, 0.7, 0.0, 0.0};
  const int numbers[2] = {1, 8};
  const double box[9] = {10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0};
  double forces[6] = {};
  double energy = 0.0;
  double variance = -1.0;
  pot->force(2, positions, numbers, forces, &energy, &variance, box);
  pot->force(2, positions, numbers, forces, &energy, &variance, box);
  pot.reset();
  client.join();
  std::filesystem::current_path(previous);

  constexpr double hartree = 27.211386245988;
  constexpr double bohr = 0.529177210903;
  REQUIRE(client_ok);
  REQUIRE(energy == Catch::Approx(1.5 * hartree));
  REQUIRE(forces[0] == Catch::Approx(0.2 * hartree / bohr));
  REQUIRE(variance == 0.0);
  std::ifstream template_in(dir / "nwchem_socket.nwi");
  std::stringstream template_text;
  template_text << template_in.rdbuf();
  REQUIRE(template_text.str().find("socket unix eonqf") != std::string::npos);
  std::filesystem::remove_all(dir);
}

TEST_CASE("SocketNWChem rejects a long socket name and a bad host",
          "[socketnwchem]") {
  eonc::Parameters params;
  eonc::ParametersLoadAccess::potential_options(params).potential =
      eonc::PotType::SocketNWChem;
  auto &opt = eonc::ParametersLoadAccess::socket_nwchem_options(params);
  opt.unix_socket_mode = true;
  opt.unix_socket_path = "this_name_is_far_too_long";
  REQUIRE_THROWS_AS(SocketNWChemPot{params}, std::runtime_error);

  opt.unix_socket_mode = false;
  opt.host = "not-an-ip";
  opt.port = 9;
  REQUIRE_THROWS_AS(SocketNWChemPot{params}, std::runtime_error);

  opt.host = "127.0.0.1";
  opt.port = 40000 + (::getpid() % 20000);
  SocketNWChemPot listening(params);
  REQUIRE(listening.getType() == eonc::PotType::SocketNWChem);
}
