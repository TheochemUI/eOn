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
            [](py::object /*self*/) { return JobStrings::PROCESS_SEARCH; })
        .def_property_readonly_static(
            "SADDLE_SEARCH",
            [](py::object /*self*/) { return JobStrings::SADDLE_SEARCH; })
        .def_property_readonly_static(
            "MINIMIZATION", [](py::object /*self*/) { return JobStrings::MINIMIZATION; })
        .def_property_readonly_static(
            "POINT",
            [](py::object /*self*/) { return JobStrings::POINT; })
        .def_property_readonly_static(
            "PARALLEL_REPLICA",
            [](py::object /*self*/) { return JobStrings::PARALLEL_REPLICA; })
        .def_property_readonly_static(
            "SAFE_HYPER",
            [](py::object /*self*/) { return JobStrings::SAFE_HYPER; })
        .def_property_readonly_static(
            "TAD",
            [](py::object /*self*/) { return JobStrings::TAD; })
        .def_property_readonly_static(
            "REPLICA_EXCHANGE",
            [](py::object /*self*/) { return JobStrings::REPLICA_EXCHANGE; })
        .def_property_readonly_static(
            "BASIN_HOPPING",
            [](py::object /*self*/) { return JobStrings::BASIN_HOPPING; })
        .def_property_readonly_static(
            "HESSIAN",
            [](py::object /*self*/) { return JobStrings::HESSIAN; })
        .def_property_readonly_static(
            "FINITE_DIFFERENCE",
            [](py::object /*self*/) { return JobStrings::FINITE_DIFFERENCE; })
        .def_property_readonly_static(
            "NUDGED_ELASTIC_BAND",
            [](py::object /*self*/) { return JobStrings::NUDGED_ELASTIC_BAND; })
        .def_property_readonly_static(
            "DYNAMICS",
            [](py::object /*self*/) { return JobStrings::DYNAMICS; })
        .def_property_readonly_static(
            "PREFACTOR",
            [](py::object /*self*/) { return JobStrings::PREFACTOR; })
        .def_property_readonly_static(
            "GLOBAL_OPTIMIZATION",
            [](py::object /*self*/) { return JobStrings::GLOBAL_OPTIMIZATION; })
        .def_property_readonly_static(
            "STRUCTURE_COMPARISON",
            [](py::object /*self*/) { return JobStrings::STRUCTURE_COMPARISON; })
        .def_property_readonly_static(
            "MONTE_CARLO",
            [](py::object /*self*/) { return JobStrings::MONTE_CARLO; })
        // Unused.
        // .def_property_readonly_static(
        //     "TEST",
        //     [](py::object /*self*/) { return JobStrings::TEST; })

        /*
        ** Python helpers
        */

        .def("__repr__", [](const Job &a) { return "<Job object>"; });
}
