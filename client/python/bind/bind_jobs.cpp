/*
** Job factory + ClientEON pipeline steps (each step bound separately).
**
** ClientEON sequence (bind each, compose in Python):
**   1. Parameters.load(config.ini)
**   2. make_job / make_potential + Matter
**   3. job.run()  OR  Matter.con2matter → relax → matter2con
**   4. destroy job/potential handles (Python GC / del)
**   5. write_potcall_summary()
**   6. append_results_timing(...)
*/
#include "eon/HelperFunctions.h"
#include "eon/Job.h"
#include "eon/Parameters.h"
#include "eon/PotRegistry.h"

#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <chrono>
#include <format>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace eonc::pybind {
namespace nb = nanobind;

void bind_jobs(nb::module_ &m) {
  nb::class_<eonc::Job>(m, "Job",
                        "Abstract eOn client job (Minimization, NEB, …)")
      .def("get_type", &eonc::Job::getType)
      .def(
          "run",
          [](eonc::Job &self) {
            // Torch autograd in RGPOT metatomic engines refuses the GIL.
            // Release for the whole job (force evals dominate runtime).
            nb::gil_scoped_release release;
            return self.run();
          },
          "Run this job in the *current* working directory. Writes job "
          "artifacts (min.con / neb.dat / results.dat body). Does **not** "
          "write _potcalls.json or timing — call write_potcall_summary and "
          "append_results_timing after destroying the Job.");

  m.def(
      "make_job",
      [](eonc::Parameters &params) {
        auto job =
            eonc::helpers::makeJob(std::make_unique<eonc::Parameters>(params));
        if (!job)
          throw std::runtime_error("make_job: unknown or unsupported job type");
        return std::shared_ptr<eonc::Job>(std::move(job));
      },
      nb::arg("parameters"),
      "Construct a Job from Parameters.main.job (copies Parameters + builds "
      "Potential).");

  // --- ClientEON post-job steps (explicit, not buried in a black box) ---

  m.def(
      "write_potcall_summary",
      [](const std::string &path) {
        PotRegistry::get().write_summary(path);
        return path;
      },
      nb::arg("path") = "_potcalls.json",
      "Write PotRegistry force-call summary (call after Job/Potential are "
      "destroyed so tallies are recorded).");

  m.def(
      "get_process_times",
      []() {
        double real = 0, user = 0, sys = 0;
        eonc::helpers::getTime(&real, &user, &sys);
        return nb::make_tuple(real, user, sys);
      },
      "Return (real, user, system) process times in seconds (same as "
      "ClientEON).");

  m.def(
      "append_results_timing",
      [](const std::string &path, double elapsed_seconds, double user_time,
         double system_time) {
        std::ofstream out(path, std::ios::app);
        if (!out.is_open())
          throw std::runtime_error("append_results_timing: cannot open " +
                                   path);
        // results.dat contract: "<value> <key>" (parse_results /
        // eon_schema.jobs)
        out << std::format("{:.12e} time_seconds\n", elapsed_seconds);
#ifndef _WIN32
        out << std::format("{:.12e} user_time\n", user_time);
        out << std::format("{:.12e} system_time\n", system_time);
#endif
      },
      nb::arg("path") = "results.dat", nb::arg("elapsed_seconds"),
      nb::arg("user_time") = 0.0, nb::arg("system_time") = 0.0,
      "Append ClientEON timing footer to results.dat (value-key lines).");

  m.def(
      "steady_clock_now",
      []() {
        using clock = std::chrono::steady_clock;
        return std::chrono::duration<double>(clock::now().time_since_epoch())
            .count();
      },
      "Monotonic seconds for measuring wall time around job steps.");
}

} // namespace eonc::pybind
