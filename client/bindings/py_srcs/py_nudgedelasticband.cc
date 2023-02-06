// clang-format off
#include "py_wrapper.hpp"
#include "../../NudgedElasticBand.h"
// clang-format on
using ListCasterBase = pybind11::detail::list_caster<std::vector<Matter *>, Matter *>;
namespace pybind11 {
namespace detail {
template <>
struct type_caster<std::vector<Matter *>> : ListCasterBase {
    static handle cast(const std::vector<Matter *> &src, return_value_policy, handle parent) {
        return ListCasterBase::cast(src, return_value_policy::reference, parent);
    }
    static handle cast(const std::vector<Matter *> *src, return_value_policy pol, handle parent) {
        return cast(*src, pol, parent);
    }
};
} // namespace detail
} // namespace pybind11
void py_nudgedelasticband(py::module_ &m) {
    py::class_<NudgedElasticBand>(m, "NudgedElasticBand")
        /*
        ** Constructors
        */
        .def(py::init<Matter * /*initialPassed*/,
                      Matter * /*finalPassed*/,
                      Parameters * /*parametersPassed*/>())

        /*
        ** Methods
        */
        .def("compute", &NudgedElasticBand::compute)
        .def("updateForces", &NudgedElasticBand::updateForces)
        .def("convergenceForce", &NudgedElasticBand::convergenceForce)
        .def("findExtrema", &NudgedElasticBand::findExtrema)
        .def("printImageData", &NudgedElasticBand::printImageData)

        /*
        ** Parameters
        */
        .def_readwrite("atoms", &NudgedElasticBand::atoms)
        .def_readwrite("nimages", &NudgedElasticBand::images)
        .def_readwrite("climbingImage", &NudgedElasticBand::climbingImage)
        .def_readwrite("numExtrema", &NudgedElasticBand::numExtrema)
        .def_property(
            "neb_images",
            [](const NudgedElasticBand *neb) -> std::reference_wrapper<std::vector<Matter *>> {
                // TODO: Does this leak?
                // Kanged: https://github.com/pybind/pybind11/issues/637
                auto *ret = new std::vector<Matter *>{neb->image, neb->image + neb->images};
                return *ret;
            },
            [](const NudgedElasticBand *neb, std::vector<Matter *> matvec) -> void {
                size_t iter = 0;
                for (auto &&ref : matvec) {
                    neb->image[iter] = ref;
                    iter++;
                }
            })
        .def_property(
            "neb_tangents",
            [](const NudgedElasticBand *neb) -> std::reference_wrapper<std::vector<AtomMatrix *>> {
                // TODO: Does this leak?
                // Kanged: https://github.com/pybind/pybind11/issues/637
                auto *ret
                    = new std::vector<AtomMatrix *>{neb->tangent, neb->tangent + neb->images};
                return *ret;
            },
            [](const NudgedElasticBand *neb, std::vector<AtomMatrix *> matvec) -> void {
                size_t iter = 0;
                for (auto &&ref : matvec) {
                    neb->tangent[iter] = ref;
                    iter++;
                }
            })
        .def_property(
            "neb_forces",
            [](const NudgedElasticBand *neb) -> std::reference_wrapper<std::vector<AtomMatrix *>> {
                // TODO: Does this leak?
                // Kanged: https://github.com/pybind/pybind11/issues/637
                auto *ret = new std::vector<AtomMatrix *>{
                    neb->projectedForce, neb->projectedForce + neb->images};
                return *ret;
            },
            [](const NudgedElasticBand *neb, std::vector<AtomMatrix *> matvec) -> void {
                size_t iter = 0;
                for (auto &&ref : matvec) {
                    neb->projectedForce[iter] = ref;
                    iter++;
                }
            })
        .def_readwrite("movedAfterForceCall", &NudgedElasticBand::movedAfterForceCall)
        .def_readwrite("extremumEnergy", &NudgedElasticBand::extremumEnergy)
        .def_readwrite("extremumPosition", &NudgedElasticBand::extremumPosition)
        .def_readwrite("extremumCurvature", &NudgedElasticBand::extremumCurvature)
        .def_readwrite("maxEnergyImage", &NudgedElasticBand::maxEnergyImage)

        /*
        ** Python helpers
        */

        .def("__repr__", [](const NudgedElasticBand &a) { return "<NudgedElasticBand object>"; });

    /* Internal Types and Enums */
    py::enum_<NudgedElasticBand::NEBStatus>(m, "nebStatus")
        .value("STATUS_GOOD", NudgedElasticBand::NEBStatus::STATUS_GOOD)
        .value("STATUS_INIT", NudgedElasticBand::NEBStatus::STATUS_INIT)
        .value("STATUS_BAD_MAX_ITERATIONS",
               NudgedElasticBand::NEBStatus::STATUS_BAD_MAX_ITERATIONS);
    // .export_values(); // exports entries into the parent which makes no sense
}
