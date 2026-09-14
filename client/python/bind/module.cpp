/*
** pyeonclient._core — complete nanobind surface for the eOn C++ client.
*/
#include <nanobind/nanobind.h>

// client/python/meson.build reads this from pyproject-pyeonclient.toml.
#ifndef PYEONCLIENT_VERSION
#define PYEONCLIENT_VERSION "0.0.0+unknown"
#endif

namespace eonc::pybind {
namespace nb = nanobind;

void bind_enums(nb::module_ &m);
void bind_parameters(nb::module_ &m);
void bind_potential(nb::module_ &m);
void bind_matter(nb::module_ &m);
void bind_ase(nb::module_ &m);
void bind_jobs(nb::module_ &m);
void bind_neb(nb::module_ &m);
void bind_eigenmode(nb::module_ &m);
void bind_saddle(nb::module_ &m);
void bind_analysis(nb::module_ &m);
void bind_sampling(nb::module_ &m);

} // namespace eonc::pybind

NB_MODULE(_core, m) {
  m.doc() =
      "eOn client core via nanobind: Matter, Parameters, Potential, "
      "NudgedElasticBand, Dimer/ImprovedDimer/Lanczos/Davidson, "
      "MinModeSaddleSearch, Hessian, Prefactor, Dynamics/MC/BH/process_search, "
      "ASE bulk matter_from_ase/matter_to_ase. "
      "Stable ABI (abi3) or free-threaded.";
  m.attr("__version__") = PYEONCLIENT_VERSION;

  m.def(
      "built_with_metatomic",
      []() {
#ifdef WITH_METATOMIC
        return true;
#else
        return false;
#endif
      },
      "True if this wheel/extension was compiled with -Dwith_metatomic=true");
  m.def(
      "built_with_rgpot",
      []() {
#ifdef WITH_RGPOT
        return true;
#else
        return false;
#endif
      },
      "True if compiled with -Dwith_rgpot=true (RGPOT pot / engine dlopen)");

  eonc::pybind::bind_enums(m);
  eonc::pybind::bind_parameters(m);
  eonc::pybind::bind_potential(m);
  eonc::pybind::bind_matter(m);
  eonc::pybind::bind_ase(m);
  eonc::pybind::bind_jobs(m);
  eonc::pybind::bind_neb(m);
  eonc::pybind::bind_eigenmode(m);
  eonc::pybind::bind_saddle(m);
  eonc::pybind::bind_analysis(m);
  eonc::pybind::bind_sampling(m);
}
