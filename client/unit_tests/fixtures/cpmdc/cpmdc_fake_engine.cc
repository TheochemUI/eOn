#if defined(_WIN32) || defined(_WIN64)
#define RGPOT_CPMDC_BUILD
#endif
#include "rgpot/CPMDPot/cpmd_c_abi.h"
#include "rgpot/abi/Compat.hpp"
#include "rgpot/rpc/Potentials.capnp.h"

#include "cpmd_stress_oracle.hpp"

#include <capnp/message.h>
#include <capnp/serialize.h>
#include <kj/array.h>
#include <kj/exception.h>

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

struct CPMDCSession {
  // Use bytes, not vector<capnp::word>: word is not default-constructible on
  // Apple libc++ (resize/assign fail in CI Xcode 16).
  std::vector<unsigned char> params;
  long long orbital_key = -1;
  bool has_orbital_key = false;
};

namespace {

constexpr CPMDCFeatureEntry kFeatures[] = {
#include "cpmd_feature_table.inc"
};

bool has_flat_message(const void *msg, size_t msg_size_bytes) {
  return msg != nullptr && msg_size_bytes >= sizeof(::capnp::word) &&
         (msg_size_bytes % sizeof(::capnp::word)) == 0;
}

::ForceInput::Reader read_force_input(const void *force_input_capnp,
                                      size_t force_input_capnp_size_bytes,
                                      ::capnp::FlatArrayMessageReader &reader) {
  (void)force_input_capnp;
  (void)force_input_capnp_size_bytes;
  return reader.getRoot<::ForceInput>();
}

bool force_input_ok(const ::ForceInput::Reader &input, size_t *force_count,
                    double *cell_zz) {
  const auto positions = input.getPos();
  const auto atomic_numbers = input.getAtmnrs();
  const auto box = input.getBox();
  if (positions.size() == 0 || positions.size() % 3 != 0)
    return false;
  const size_t n_atoms = positions.size() / 3;
  if (atomic_numbers.size() != n_atoms || box.size() != 9)
    return false;
  if (std::string(input.getLengthUnit().cStr()) != "angstrom" ||
      std::string(input.getEnergyUnit().cStr()) != "eV")
    return false;
  *force_count = positions.size();
  *cell_zz = box[8];
  return std::isfinite(*cell_zz) && *cell_zz > 0.0;
}

std::vector<unsigned char> make_result(size_t force_count, double cell_zz) {
  ::capnp::MallocMessageBuilder msg;
  auto result = msg.initRoot<::PotentialResult>();
  const double hartree = cpmd_stress_oracle::sessionHartree(cell_zz);
  result.setEnergy(hartree * rgpot::units::HARTREE_TO_EV);
  auto forces = result.initForces(static_cast<unsigned int>(force_count));
  for (unsigned int i = 0; i < forces.size(); ++i)
    forces.set(i, 0.011 + 0.001 * static_cast<double>(i));
  auto stress = result.initStress(9);
  for (unsigned int i = 0; i < 9; ++i)
    stress.set(i,
               cpmd_stress_oracle::stressEvPerAngstrom3(static_cast<int>(i)));
  auto words = ::capnp::messageToFlatArray(msg);
  const auto bytes = words.asBytes();
  return std::vector<unsigned char>(bytes.begin(), bytes.end());
}

CPMDCResult fail_result(const char *message) {
  CPMDCResult r;
  r.ok = 0;
  r.energy_h = 0.0;
  std::snprintf(r.message, sizeof(r.message), "%s", message);
  return r;
}

CPMDCResult ok_result(const char *message) {
  CPMDCResult r;
  r.ok = 1;
  r.energy_h = 0.0;
  std::snprintf(r.message, sizeof(r.message), "%s", message);
  return r;
}

// Sessions created through cpmdc_session_create, read back by tests
// through cpmdc_fake_session_create_count.
int g_session_create_count = 0;
// Key named for the last session evaluation; kNoOrbitalKey when none was.
constexpr long long kNoOrbitalKey = -1000000;
long long g_last_orbital_key = kNoOrbitalKey;

} // namespace

extern "C" {

int cpmdc_set_params(const void *params_capnp, size_t params_capnp_size_bytes) {
  return has_flat_message(params_capnp, params_capnp_size_bytes) ? 0 : -1;
}

CPMDCResult cpmdc_energy_gradient(int n_atoms, const double *positions_ang,
                                  const int *atomic_numbers,
                                  const void *params_capnp,
                                  size_t params_capnp_size_bytes,
                                  double *grad_h_bohr) {
  (void)n_atoms;
  (void)positions_ang;
  (void)atomic_numbers;
  (void)params_capnp;
  (void)params_capnp_size_bytes;
  (void)grad_h_bohr;
  return fail_result("fake engine requires session result path");
}

CPMDCResult cpmdc_energy(int n_atoms, const double *positions_ang,
                         const int *atomic_numbers, const void *params_capnp,
                         size_t params_capnp_size_bytes) {
  (void)n_atoms;
  (void)positions_ang;
  (void)atomic_numbers;
  (void)params_capnp;
  (void)params_capnp_size_bytes;
  return fail_result("fake engine requires session result path");
}

CPMDCResult cpmdc_energy_forces(int n_atoms, const double *positions_ang,
                                const int *atomic_numbers,
                                const void *params_capnp,
                                size_t params_capnp_size_bytes,
                                double *forces_h_bohr) {
  (void)n_atoms;
  (void)positions_ang;
  (void)atomic_numbers;
  (void)params_capnp;
  (void)params_capnp_size_bytes;
  (void)forces_h_bohr;
  return fail_result("fake engine requires session result path");
}

CPMDCSession *cpmdc_session_create(const void *params_capnp,
                                   size_t params_capnp_size_bytes) {
  if (!has_flat_message(params_capnp, params_capnp_size_bytes))
    return nullptr;
  ++g_session_create_count;
  auto *session = new CPMDCSession;
  session->params.resize(params_capnp_size_bytes);
  std::memcpy(session->params.data(), params_capnp, params_capnp_size_bytes);
  return session;
}

int cpmdc_session_set_params(CPMDCSession *session, const void *params_capnp,
                             size_t params_capnp_size_bytes) {
  if (session == nullptr ||
      !has_flat_message(params_capnp, params_capnp_size_bytes))
    return -1;
  session->params.resize(params_capnp_size_bytes);
  std::memcpy(session->params.data(), params_capnp, params_capnp_size_bytes);
  return 0;
}

void cpmdc_session_destroy(CPMDCSession *session) { delete session; }

int cpmdc_session_select_orbitals(CPMDCSession *session, long long key) {
  if (session == nullptr)
    return -1;
  session->orbital_key = key;
  session->has_orbital_key = true;
  return 0;
}

CPMDCResult cpmdc_session_energy_gradient(CPMDCSession *session, int n_atoms,
                                          const double *positions_ang,
                                          const int *atomic_numbers,
                                          double *grad_h_bohr) {
  (void)session;
  (void)n_atoms;
  (void)positions_ang;
  (void)atomic_numbers;
  (void)grad_h_bohr;
  return fail_result("fake engine requires session result path");
}

CPMDCResult cpmdc_session_energy(CPMDCSession *session, int n_atoms,
                                 const double *positions_ang,
                                 const int *atomic_numbers) {
  (void)session;
  (void)n_atoms;
  (void)positions_ang;
  (void)atomic_numbers;
  return fail_result("fake engine requires session result path");
}

CPMDCResult cpmdc_session_energy_forces(CPMDCSession *session, int n_atoms,
                                        const double *positions_ang,
                                        const int *atomic_numbers,
                                        double *forces_h_bohr) {
  (void)session;
  (void)n_atoms;
  (void)positions_ang;
  (void)atomic_numbers;
  (void)forces_h_bohr;
  return fail_result("fake engine requires session result path");
}

size_t cpmdc_potential_result_size_for_force_input(
    const void *force_input_capnp, size_t force_input_capnp_size_bytes) {
  if (!has_flat_message(force_input_capnp, force_input_capnp_size_bytes))
    return 0;
  try {
    auto words =
        kj::arrayPtr(static_cast<const ::capnp::word *>(force_input_capnp),
                     force_input_capnp_size_bytes / sizeof(::capnp::word));
    ::capnp::FlatArrayMessageReader reader(words);
    size_t force_count = 0;
    double cell_zz = 0.0;
    if (!force_input_ok(read_force_input(force_input_capnp,
                                         force_input_capnp_size_bytes, reader),
                        &force_count, &cell_zz))
      return 0;
    return make_result(force_count, cell_zz).size();
  } catch (const kj::Exception &) {
    return 0;
  }
}

CPMDCResult cpmdc_session_calculate_forces(CPMDCSession *session,
                                           const void *force_input_capnp,
                                           size_t force_input_capnp_size_bytes,
                                           double *forces_h_bohr,
                                           size_t forces_len) {
  if (session == nullptr ||
      !has_flat_message(force_input_capnp, force_input_capnp_size_bytes) ||
      forces_h_bohr == nullptr)
    return fail_result("invalid fake session arguments");
  try {
    auto words =
        kj::arrayPtr(static_cast<const ::capnp::word *>(force_input_capnp),
                     force_input_capnp_size_bytes / sizeof(::capnp::word));
    ::capnp::FlatArrayMessageReader reader(words);
    size_t force_count = 0;
    double cell_zz = 0.0;
    if (!force_input_ok(read_force_input(force_input_capnp,
                                         force_input_capnp_size_bytes, reader),
                        &force_count, &cell_zz))
      return fail_result("invalid fake ForceInput");
    if (forces_len < force_count)
      return fail_result("force buffer too small");
    for (size_t i = 0; i < force_count; ++i)
      forces_h_bohr[i] = 0.011 + 0.001 * static_cast<double>(i);
    auto result = ok_result("session forces ok");
    result.energy_h = 0.75 + 0.001 * cell_zz;
    return result;
  } catch (const kj::Exception &ex) {
    return fail_result(ex.getDescription().cStr());
  }
}

CPMDCResult cpmdc_session_calculate_result(
    CPMDCSession *session, const void *force_input_capnp,
    size_t force_input_capnp_size_bytes, void *potential_result_capnp,
    size_t potential_result_capnp_capacity_bytes,
    size_t *potential_result_capnp_size_bytes) {
  if (session == nullptr ||
      !has_flat_message(force_input_capnp, force_input_capnp_size_bytes) ||
      potential_result_capnp_size_bytes == nullptr)
    return fail_result("invalid fake session arguments");
  try {
    auto words =
        kj::arrayPtr(static_cast<const ::capnp::word *>(force_input_capnp),
                     force_input_capnp_size_bytes / sizeof(::capnp::word));
    ::capnp::FlatArrayMessageReader reader(words);
    size_t force_count = 0;
    double cell_zz = 0.0;
    if (!force_input_ok(read_force_input(force_input_capnp,
                                         force_input_capnp_size_bytes, reader),
                        &force_count, &cell_zz))
      return fail_result("invalid fake ForceInput");
    const auto result_bytes = make_result(force_count, cell_zz);
    g_last_orbital_key =
        session->has_orbital_key ? session->orbital_key : kNoOrbitalKey;
    const size_t required = result_bytes.size();
    *potential_result_capnp_size_bytes = required;
    if (potential_result_capnp == nullptr ||
        potential_result_capnp_capacity_bytes < required)
      return fail_result("result buffer too small");
    std::memcpy(potential_result_capnp, result_bytes.data(), required);
    return ok_result("session result ok");
  } catch (const kj::Exception &ex) {
    return fail_result(ex.getDescription().cStr());
  }
}

CPMDCResult cpmdc_calculate_result(const void *params_capnp,
                                   size_t params_capnp_size_bytes,
                                   const void *force_input_capnp,
                                   size_t force_input_capnp_size_bytes,
                                   void *potential_result_capnp,
                                   size_t potential_result_capnp_capacity_bytes,
                                   size_t *potential_result_capnp_size_bytes) {
  CPMDCSession *session =
      cpmdc_session_create(params_capnp, params_capnp_size_bytes);
  if (session == nullptr)
    return fail_result("invalid fake params");
  CPMDCResult result = cpmdc_session_calculate_result(
      session, force_input_capnp, force_input_capnp_size_bytes,
      potential_result_capnp, potential_result_capnp_capacity_bytes,
      potential_result_capnp_size_bytes);
  cpmdc_session_destroy(session);
  return result;
}

const char *cpmdc_version(void) { return "cpmdc-fake/0.2"; }

int cpmdc_available(void) { return 1; }

void cpmdc_finalize(void) {}

size_t cpmdc_feature_count(void) {
  return sizeof(kFeatures) / sizeof(kFeatures[0]);
}

const CPMDCFeatureEntry *cpmdc_feature_table(void) { return kFeatures; }

const CPMDCFeatureEntry *cpmdc_feature_find(const char *feature_id) {
  if (feature_id == nullptr)
    return nullptr;
  for (const auto &feature : kFeatures) {
    if (std::strcmp(feature.feature_id, feature_id) == 0)
      return &feature;
  }
  return nullptr;
}

// Minimum-profile symbols (potentials-schema PROFILE.md) so the fake serves
// as the ProfileLoader conformance target. Not all of these appear in
// cpmd_c_abi.h; export them explicitly for MSVC GetProcAddress.

#if defined(_WIN32) || defined(_WIN64)
#define RGPOT_FAKE_ONLY __declspec(dllexport)
#else
#define RGPOT_FAKE_ONLY
#endif

RGPOT_FAKE_ONLY int cpmdc_abi_version(void) { return 0; }

// Test-only: number of cpmdc_session_create calls that returned a
// session. Lets a test assert that a CPMDPot kept or recreated its
// session.
RGPOT_FAKE_ONLY int cpmdc_fake_session_create_count(void) {
  return g_session_create_count;
}

// Test-only: the key named before the last session evaluation, or
// -1000000 when the session never named one.
RGPOT_FAKE_ONLY long long cpmdc_fake_last_orbital_key(void) {
  return g_last_orbital_key;
}

RGPOT_FAKE_ONLY const char *cpmdc_last_error(void) { return ""; }

RGPOT_FAKE_ONLY int cpmdc_configure(const void *config_capnp,
                                    size_t config_capnp_size_bytes) {
  return has_flat_message(config_capnp, config_capnp_size_bytes) ? 0 : -1;
}

RGPOT_FAKE_ONLY CPMDCSession *
cpmdc_session_create_from_config(const void *config_capnp,
                                 size_t config_capnp_size_bytes) {
  if (!has_flat_message(config_capnp, config_capnp_size_bytes))
    return nullptr;
  auto *session = new CPMDCSession();
  const auto *bytes = static_cast<const unsigned char *>(config_capnp);
  session->params.assign(bytes, bytes + config_capnp_size_bytes);
  return session;
}

RGPOT_FAKE_ONLY int cpmdc_session_configure(CPMDCSession *session,
                                            const void *config_capnp,
                                            size_t config_capnp_size_bytes) {
  if (session == nullptr ||
      !has_flat_message(config_capnp, config_capnp_size_bytes))
    return -1;
  const auto *bytes = static_cast<const unsigned char *>(config_capnp);
  session->params.assign(bytes, bytes + config_capnp_size_bytes);
  return 0;
}

RGPOT_FAKE_ONLY CPMDCResult cpmdc_calculate_result_from_config(
    const void *config_capnp, size_t config_capnp_size_bytes,
    const void *force_input_capnp, size_t force_input_capnp_size_bytes,
    void *potential_result_capnp, size_t potential_result_capnp_capacity_bytes,
    size_t *potential_result_capnp_size_bytes) {
  CPMDCSession *session =
      cpmdc_session_create_from_config(config_capnp, config_capnp_size_bytes);
  if (session == nullptr)
    return fail_result("invalid fake config");
  CPMDCResult result = cpmdc_session_calculate_result(
      session, force_input_capnp, force_input_capnp_size_bytes,
      potential_result_capnp, potential_result_capnp_capacity_bytes,
      potential_result_capnp_size_bytes);
  cpmdc_session_destroy(session);
  return result;
}

RGPOT_FAKE_ONLY int
cpmdc_capabilities_result(void *capabilities_capnp,
                          size_t capabilities_capnp_capacity_bytes,
                          size_t *capabilities_capnp_size_bytes) {
  if (capabilities_capnp_size_bytes == nullptr)
    return -1;
  try {
    ::capnp::MallocMessageBuilder msg;
    auto caps = msg.initRoot<::Capabilities>();
    caps.setBackendName("cpmdc");
    caps.setBackendVersion(cpmdc_version());
    caps.setAbiVersion(cpmdc_abi_version());
    caps.setAvailable(cpmdc_available() != 0);
    auto ops = caps.initOperations(3);
    ops.set(0, ::Capabilities::Operation::ENERGY);
    ops.set(1, ::Capabilities::Operation::FORCES);
    ops.set(2, ::Capabilities::Operation::GRADIENT);
    auto kinds = caps.initConfigKinds(1);
    kinds.set(0, "cpmd");
    caps.setSchemaVersion("fake");
    rgpot::abi::fill_compatibility(caps);
    auto words = ::capnp::messageToFlatArray(msg);
    const auto bytes = words.asBytes();
    *capabilities_capnp_size_bytes = bytes.size();
    if (capabilities_capnp == nullptr ||
        capabilities_capnp_capacity_bytes < bytes.size())
      return -1;
    std::memcpy(capabilities_capnp, bytes.begin(), bytes.size());
    return 0;
  } catch (const kj::Exception &) {
    *capabilities_capnp_size_bytes = 0;
    return -1;
  }
}

} // extern "C"
