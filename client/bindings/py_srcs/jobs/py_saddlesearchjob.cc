// clang-format off
#include "../py_wrapper.hpp"
#include "../py_job.hpp"
// Binding
#include "../../../SaddleSearchJob.h"
// clang-format on

template <class SaddleSearchJobBase = SaddleSearchJob>
class PySaddleSearchJob : public PyJob<SaddleSearchJobBase> {
public:
    using PyJob<SaddleSearchJobBase>::PyJob; // Inherit constructor
    // Override pure virtual with non-pure
    std::vector<std::string> run() override {
        PYBIND11_OVERRIDE(std::vector<std::string>, SaddleSearchJobBase, run, );
    };
};

void py_saddlesearchjob(py::module_ &m) {
    py::class_<SaddleSearchJob, Job, PySaddleSearchJob<> >(m, "SaddleSearchJob")
    /*
    ** Constructors
    */
    .def(py::init<Parameters */*params*/>())

    /*
    ** Methods
    */
    .def("run", &SaddleSearchJob::run)

    /*
    ** Python helpers
    */

    .def("__repr__", [](const SaddleSearchJob &a) { return "<SaddleSearchJob object>"; });
}
