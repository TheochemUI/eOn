#include "py_job.hpp"

void py_job(py::module_ &m) {
    py::class_<Job, PyJob<> >(m, "Job")
        .def(py::init())
        /*
        ** Functions
        */
        .def_static("getJob", &Job::getJob)
        .def("run", &Job::run)

        /*
        ** Parameters
        */
        .def_property_readonly_static(
            "PROCESS_SEARCH",
            [](py::object /*self*/) { return Job::PROCESS_SEARCH; })
        .def_property_readonly_static(
            "SADDLE_SEARCH",
            [](py::object /*self*/) { return Job::SADDLE_SEARCH; })
        .def_property_readonly_static(
            "MINIMIZATION", [](py::object /*self*/) { return Job::MINIMIZATION; })
        .def_property_readonly_static(
            "POINT",
            [](py::object /*self*/) { return Job::POINT; })
        .def_property_readonly_static(
            "PARALLEL_REPLICA",
            [](py::object /*self*/) { return Job::PARALLEL_REPLICA; })
        .def_property_readonly_static(
            "SAFE_HYPER",
            [](py::object /*self*/) { return Job::SAFE_HYPER; })
        .def_property_readonly_static(
            "TAD",
            [](py::object /*self*/) { return Job::TAD; })
        .def_property_readonly_static(
            "REPLICA_EXCHANGE",
            [](py::object /*self*/) { return Job::REPLICA_EXCHANGE; })
        .def_property_readonly_static(
            "BASIN_HOPPING",
            [](py::object /*self*/) { return Job::BASIN_HOPPING; })
        .def_property_readonly_static(
            "HESSIAN",
            [](py::object /*self*/) { return Job::HESSIAN; })
        .def_property_readonly_static(
            "FINITE_DIFFERENCE",
            [](py::object /*self*/) { return Job::FINITE_DIFFERENCE; })
        .def_property_readonly_static(
            "NUDGED_ELASTIC_BAND",
            [](py::object /*self*/) { return Job::NUDGED_ELASTIC_BAND; })
        .def_property_readonly_static(
            "DYNAMICS",
            [](py::object /*self*/) { return Job::DYNAMICS; })
        .def_property_readonly_static(
            "PREFACTOR",
            [](py::object /*self*/) { return Job::PREFACTOR; })
        .def_property_readonly_static(
            "GLOBAL_OPTIMIZATION",
            [](py::object /*self*/) { return Job::GLOBAL_OPTIMIZATION; })
        .def_property_readonly_static(
            "STRUCTURE_COMPARISON",
            [](py::object /*self*/) { return Job::STRUCTURE_COMPARISON; })
        .def_property_readonly_static(
            "MONTE_CARLO",
            [](py::object /*self*/) { return Job::MONTE_CARLO; })
        // Unused.
        // .def_property_readonly_static(
        //     "TEST",
        //     [](py::object /*self*/) { return Job::TEST; })

        /*
        ** Python helpers
        */

        .def("__repr__", [](const Job &a) { return "<Job object>"; });
}
