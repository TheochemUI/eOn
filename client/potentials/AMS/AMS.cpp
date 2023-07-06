//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "AMS.h"
#include <iostream>
#include <iterator>

namespace bp = boost::process;

AMS::AMS(std::shared_ptr<Parameters> p)
    : Potential(PotType::AMS, p) {
  // Get the values from the configuration
  // All the parameter values are converted to lowercase in generate_run
  this->engine = p->engine;
  this->forcefield = p->forcefield;
  this->model = p->model;
  this->xc = p->xc;
  this->resources = p->resources;
  this->basis = p->basis;
  this->engine_setup = generate_run(p);
  // Environment
  // TODO: Add more checks for how this can be set
  if (p->amshome.empty() && p->scm_tmpdir.empty() && p->scmlicense.empty() &&
      p->scm_pythondir.empty() && p->amsbin.empty() &&
      p->amsresources.empty()) {
    nativenv = boost::this_process::environment();
  } else {
    nativenv = boost::this_process::environment();
    // Some of these can be derived from the others
    nativenv["AMSHOME"] = p->amshome;
    nativenv["SCM_TMPDIR"] = p->scm_tmpdir;
    nativenv["SCMLICENSE"] = p->scmlicense;
    nativenv["SCM_PYTHONDIR"] = p->scm_pythondir;
    nativenv["AMSBIN"] = p->amsbin;
    nativenv["AMSRESOURCES"] = p->amsresources;
    nativenv["PATH"] += p->amsbin;
  }
  // Do not pass "" in the config files
  // std::cout<<nativenv["PATH"].to_string()<<std::endl;
  amsevals = 0;
  // TODO: Optimize and reuse existing files Currently each Matter will
  // recreate the folders It should instead figure out if results exist and
  // use them
  this->first_run = true;
  // Determine if the engine supports restarts
  if (engine == "MOPAC") {
    this->can_restart = false;
    this->cjob = "amsResults";
    this->pjob = "amsResults";
  } else {
    this->can_restart = true;
    this->cjob = "firstRun";
    this->pjob = "secondRun";
  }
  return;
}

void AMS::cleanMemory(void) { return; }

AMS::~AMS() { cleanMemory(); }

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

void AMS::runAMS() {
  boost::asio::io_context amsRun;
  std::future<std::string> run_out_future, run_err_future;
  std::string runout, runerr;
  chmod("run_AMS.sh", S_IRWXU);
  nativenv["AMS_JOBNAME"] = cjob;
  // assert(validate_order() == true);         // TODO: Debug only
  bp::child c("run_AMS.sh", nativenv,       // set the input
              bp::std_in.close(),           // no input
              bp::std_out > run_out_future, // STDOUT
              bp::std_err > run_err_future, // STDERR
              amsRun);
  amsRun.run();
  runout = run_out_future.get();
  runerr = run_err_future.get();
  if (runerr.find("ERROR") != std::string::npos) {
    bp::spawn("cat myrestart.in");
    bp::spawn("cat run_AMS.sh");
    throw std::runtime_error("\n AMS error while running");
  } else {
    this->amsevals = amsevals + 1;
    // std::cout << "Run completed normally" << std::endl;
    // std::cerr << "NORMAL TERMINATION\n";
  }
}

double AMS::extract_scalar_rkf(std::string key) {
  // The logic here is that the extraction must be asynchronous, as we are not
  // interested in processing this the output line by line, and would prefer to
  // have it all in one place
  std::string execString;
  std::vector<std::string> execDat;
  boost::asio::io_context rkf;
  std::future<std::string> rkf_out_future, rkf_err_future;
  std::string rkfout, rkferr;
  double xval, x;
  std::vector<double> extracted;
  execString =
      fmt::format("dmpkf {jobid:}.results/{engine:}.rkf AMSResults%{key:}",
                  fmt::arg("jobid", this->cjob),
                  fmt::arg("engine", this->engine_lower), fmt::arg("key", key));
  // std::cout << execString << "\n";
  // Extract
  bp::child c(execString, nativenv,         // execute with the environment
              bp::std_in.close(),           // no input
              bp::std_out > rkf_out_future, // STDOUT
              bp::std_err > rkf_err_future, // STDERR
              rkf);

  rkf.run(); // this blocks until the end of the command
  // Populate the strings
  rkfout = rkf_out_future.get();
  rkferr = rkf_err_future.get();
  if (rkferr.find("ERROR") != std::string::npos) {
    throw std::runtime_error(fmt::format(
        "\n AMS error while extracting {}, got:\n {}", key, rkferr));
  } else {
    // std::cout << "Extracting " << key << std::endl;
  }

  execDat = absl::StrSplit(rkfout, '\n');
  // Is energy or some other scalar property
  // First 3 lines are the headers
  // [0] = "AMSResults                      "
  // [1] = "Energy                          "
  // [2] = "         1         1         2"
  // [3] = "   0.135012547958282714E+000"
  // [4] = ""
  if (absl::SimpleAtod(execDat[3], &x)) {
    xval = x * this->energyConversion;
    // std::cout << fmt::format(
    //     "\n Got {:.4f} Hartree from AMS and converted to {:.4f} eV\n", x,
    //     xval);
    return xval;
  } else {
    throw std::runtime_error(
        fmt::format("\n Expected {}, got {} instead", key, execDat[3]));
  }
  // Never reach here
  throw std::runtime_error("Generic AMS dmpkf scalar error \n");
}

std::vector<double> AMS::extract_cartesian_rkf(std::string key) {
  std::string execString;
  std::vector<std::string> execDat, innerDat;
  boost::asio::io_context rkf;
  std::future<std::string> rkf_out_future, rkf_err_future;
  std::string rkfout, rkferr;
  double felem, x;
  std::vector<double> extracted;
  execString =
      fmt::format("dmpkf {jobid:}.results/{engine:}.rkf AMSResults%{key:}",
                  fmt::arg("jobid", this->cjob),
                  fmt::arg("engine", this->engine_lower), fmt::arg("key", key));
  // std::cout << execString << "\n";

  // Extract
  bp::child c(execString, nativenv,         // execute with the environment
              bp::std_in.close(),           // no input
              bp::std_out > rkf_out_future, // STDOUT
              bp::std_err > rkf_err_future, // STDERR
              rkf);

  rkf.run(); // this blocks until the end of the command
  // Populate the strings
  rkfout = rkf_out_future.get();
  rkferr = rkf_err_future.get();

  if (rkferr.find("ERROR") != std::string::npos) {
    throw std::runtime_error(fmt::format(
        "\n AMS error while extracting {}, got:\n {}", key, rkferr));
  } else {
    // std::cout << "Extracting " << key << std::endl;
  }

  execDat = absl::StrSplit(rkfout, '\n');
  // Assume gradients or some other x y z property
  for (int i = 3; i < execDat.size(); i++) {
    // There is one per atom, the number of which is N
    // The first three lines are the header files as before
    std::vector<std::string> strrow = absl::StrSplit(execDat[i], ' ');
    for (auto elem : strrow) {
      // This loop uses an auto variable since some of the elements of strrow
      // are often ' '
      // TODO: Optimize, perhaps with a regex
      // [0] = ""
      // [1] = ""
      // [2] = ""
      // [3] = "0.798885933949825835E-004"
      // [4] = ""
      // [5] = "-0.755400038198822616E-004"
      // [6] = ""
      // [7] = ""
      // [8] = "0.457005740252359742E-004"
      if (!elem.empty() && absl::SimpleAtod(elem, &x)) {
        felem = x * this->forceConversion;
        // std::cout << fmt::format(
        //     "\n Gradient element={grad:.4f} Hartree/Bohr from AMS\n Force "
        //     "element={force:.4f} eV/Angstrom\n",
        //     fmt::arg("grad", x), fmt::arg("force", force));
        extracted.emplace_back(felem);
      }
    }
  }
  // // Debug
  // int counter = 0;
  // for (int a = 0; a < N * 3; a++) {
  //   std::cout << fmt::format("{:.25e} ", forces[a]);
  //   counter++;
  //   if (counter % 3 == 0) {
  //     std::cout << std::endl;
  //   }
  // }
  return extracted;
}

void AMS::updateCoord(long N, const double *R) {
  // The logic used here is that we need to update the PREVIOUS job, since that
  // is the one from which the calculation is to be restarted
  // Only the third line is vaguely complex
  // https://www.scm.com/doc/Scripting/Commandline_Tools/KF_command_line_utilities.html
  std::ofstream updCoord;
  std::string execString, coordDump, newCoord;
  std::vector<std::string> execDat;
  boost::asio::io_context coordio;
  std::future<std::string> err, rdump;
  std::vector<double> gradients;
  // Prep new run
  // Get the previous run's coordinates
  execString = fmt::format("dmpkf {jobid:}.results/ams.rkf Molecule%Coords",
                           fmt::arg("jobid", pjob));
  // std::cout << execString << "\n";
  // Store Coordinates
  // TODO: Simplify this, we only need the first few lines
  bp::child cprog(execString, nativenv, bp::std_in.close(), bp::std_out > rdump,
                  bp::std_err > err, coordio);
  coordio.run();
  execDat = absl::StrSplit(rdump.get(), '\n');
  // Get the first three lines
  int counter = 0;
  for (auto j : execDat) {
    if (counter >= 3) {
      break;
    } else {
      absl::StrAppend(&newCoord, j, "\n");
    }
    counter++;
  }
  // Put the rest of the coordinates
  counter = 0;
  for (int a = 0; a < N * 3; a++) {
    absl::StrAppend(&newCoord, (R[a] * lengthConversion), "  ");
    counter++;
    if (counter % 3 == 0) {
      absl::StrAppend(&newCoord, "\n");
    }
  }
  coordDump = "#!/bin/sh\n udmpkf ";
  absl::StrAppend(&coordDump, pjob, ".results/ams.rkf <<EOF\n", newCoord,
                  "EOF");
  // std::cout << coordDump;
  updCoord.open("updCoord.sh");
  updCoord << coordDump;
  updCoord.close();
  // bp::spawn("chmod +x updCoord.sh");
  chmod("updCoord.sh", S_IRWXU);
  bp::child cuprog("updCoord.sh", nativenv, bp::std_err > bp::null);
  cuprog.wait();
  return;
}

void AMS::switchjob() {
  std::string tmp;
  // std::cout << fmt::format("\nEntered Switch:\nCurrent:{}, Previous:{}\n",
  // cjob, pjob);
  tmp = this->cjob;
  this->cjob = this->pjob;
  this->pjob = tmp;
  // std::cout << fmt::format("\nSwitched\n Current:{}, Previous:{}\n", cjob,
  // pjob);
}

void AMS::write_restart() {
  std::string restart_formatter;
  if (can_restart) {
    restart_formatter = R"(
   EngineRestart {prev:}.results/{engine:}.rkf
   LoadSystem
    File {prev:}.results/ams.rkf
    Section Molecule
   End
  )";
  } else {
    restart_formatter = R"(
   LoadSystem
    File {prev:}.results/ams.rkf
    Section Molecule
   End
  )";
  }
  std::string restart_data =
      fmt::format(restart_formatter, fmt::arg("prev", pjob),
                  fmt::arg("engine", engine_lower));
  restartFrom.open("myrestart.in");
  restartFrom << restart_data;
  restartFrom.close();
  return;
}

void AMS::force(long N, const double *R, const int *atomicNrs, double *F,
                double *U, double *variance, const double *box) {
  variance = nullptr;
  if (not can_restart or first_run) {
    // std::cout << fmt::format("\nCAN_RESTART:{}  FIRST_RUN:{}\n", can_restart,
    // first_run); This is true for all engines with no restart Also if an
    // engine supports being restarted, the first run needs this
    passToSystem(N, R, atomicNrs, box);
    runAMS();
    std::vector<double> frc = extract_cartesian_rkf("Gradients");
    double *ftest = frc.data();
    for (int i = 0; i < N; i++) {
      F[3 * i] = ftest[3 * i];
      F[3 * i + 1] = ftest[3 * i + 1];
      F[3 * i + 2] = ftest[3 * i + 2];
    }
    *U = extract_scalar_rkf("Energy");
    // Update, will still not matter for those without restarts
    if (can_restart) {
      first_run = false;
      // We need the "previous job" to be pre-populated
      // clang-format off
      // std::cout<<fmt::format("\nMoving {} to {} before switching during the first job\n", cjob, pjob);
      const auto copyOptions = std::filesystem::copy_options::overwrite_existing
                             | std::filesystem::copy_options::recursive
                             ;
      // clang-format on
      std::filesystem::copy(fmt::format("./{}.results/", cjob),
                            fmt::format("./{}.results/", pjob), copyOptions);
    }
    return;
  } else {
    // std::cout << fmt::format("\nCAN_RESTART:{}  FIRST_RUN:{}\n", can_restart,
    // first_run);
    smallSys(N, R, atomicNrs, box); // writes run_AMS.sh
    updateCoord(N, R);              // updates coordinates in previous job
    write_restart();                // writes restart file using previous job
    runAMS();
    std::vector<double> frc = extract_cartesian_rkf("Gradients");
    double *ftest = frc.data();
    for (int i = 0; i < N; i++) {
      F[3 * i] = ftest[3 * i];
      F[3 * i + 1] = ftest[3 * i + 1];
      F[3 * i + 2] = ftest[3 * i + 2];
    }
    *U = extract_scalar_rkf("Energy");
    switchjob(); // toggles the jobs
    return;
  }
  // Never reach here
  throw std::runtime_error("Generic AMS force error \n");
}

void AMS::passToSystem(long N, const double *R, const int *atomicNrs,
                       const double *box)
// Creating the standard input file that will be read by the AMS driver
{
  FILE *out;
  out = fopen("run_AMS.sh", "w");

  fprintf(out, "#!/bin/sh\n");
  fmt::print(out, fmt::format("export AMS_JOBNAME={}\n", cjob));
  fprintf(out, "$AMSBIN/ams --delete-old-results <<eor\n");
  fprintf(out, "Task SinglePoint\n");
  fprintf(out, "System\n");
  fprintf(out, " Atoms\n");
  for (int i = 0; i < N; i++) {
    fprintf(out, "  %s\t%.19f\t%.19f\t%.19f\n",
            atomicNumber2symbol(atomicNrs[i]), R[i * 3 + 0], R[i * 3 + 1],
            R[i * 3 + 2]);
  }
  fprintf(out, " End\n");
  if (not model.empty() || not forcefield.empty()) {
    fprintf(out, " Lattice\n");
    for (int i = 0; i < 3; i++) {
      fprintf(out, "  %.19f\t%.19f\t%.19f\n", box[i * 3 + 0], box[i * 3 + 1],
              box[i * 3 + 2]);
    }
    fprintf(out, " End\n");
  }
  fprintf(out, "End\n");
  fmt::print(out, engine_setup);
  fprintf(out, "Properties\n");
  fprintf(out, " Gradients\n");
  fprintf(out, "End\n");
  if (can_restart and not first_run) {
    fprintf(out, "@include myrestart.in\n");
  }
  fprintf(out, "eor");
  fclose(out);
  chmod("run_AMS.sh", S_IRWXU);
  return;
}

void AMS::smallSys(long N, const double *R, const int *atomicNrs,
                   const double *box)
// Creating the truncated input file that will be read by the AMS driver
{
  FILE *out;
  out = fopen("run_AMS.sh", "w");

  fprintf(out, "#!/bin/sh\n");
  fmt::print(out, fmt::format("export AMS_JOBNAME={}\n", cjob));
  fprintf(out, "$AMSBIN/ams --delete-old-results <<eor\n");
  fprintf(out, "Task SinglePoint\n");
  fmt::print(out, engine_setup);
  fprintf(out, "Properties\n");
  fprintf(out, " Gradients\n");
  fprintf(out, "End\n");
  fprintf(out, "@include myrestart.in\n");
  fprintf(out, "eor");
  fclose(out);
  chmod("run_AMS.sh", S_IRWXU);
  return;
}

std::string AMS::generate_run(std::shared_ptr<Parameters> p) {
  std::string engine_block; // Shadows the class variable
  // TODO: Use args everywhere, cleaner logic
  // Ensure capitals and existence
  engine.empty()
      ? throw std::runtime_error("AMS Engine is required \n")
      : std::transform(engine.begin(), engine.end(), engine.begin(), ::toupper);
  // engine is special, it is used as a filename, so we store engine_lower as a
  // lowercase version too
  engine_lower = engine;
  std::transform(engine.begin(), engine.end(), engine_lower.begin(), ::tolower);
  // Prepare the block
  if (engine == "MOPAC") {
    model.empty()
        ? throw std::runtime_error("MOPAC needs a model\n")
        : std::transform(model.begin(), model.end(), model.begin(), ::toupper);
    std::string engine_formatter = R"(
 Engine {engine:}
   Model {model:}
 EndEngine
)";
    engine_block = fmt::format(engine_formatter, fmt::arg("engine", engine),
                               fmt::arg("model", model));
    return engine_block;
  } else if (engine == "ADF" || engine == "BAND") {
    basis.empty()
        ? throw std::runtime_error("ADF/BAND need a basis\n")
        : std::transform(basis.begin(), basis.end(), basis.begin(), ::toupper);
    xc.empty() ? throw std::runtime_error("ADF/BAND need a functional\n")
               : std::transform(xc.begin(), xc.end(), xc.begin(), ::toupper);
    std::string engine_formatter = R"(
   Engine {}
     Basis
       Type {}
     End

     XC
       {}
     End
   EndEngine
  )";
    engine_block = fmt::format(engine_formatter, engine, basis, xc);
    return engine_block;
  } else if (engine == "DFTB") {
    resources.empty() ? throw std::runtime_error("DFTB need resources\n")
                      : std::transform(resources.begin(), resources.end(),
                                       resources.begin(), ::toupper);
    std::string engine_formatter = R"(
   Engine {}
     ResourcesDir {}
   EndEngine
  )";
    engine_block = fmt::format(engine_formatter, engine, resources);
    return engine_block;
  } else if (engine == "reaxff") {
    forcefield.empty() ? throw std::runtime_error("REAXFF needs a forcefield\n")
                       : std::transform(forcefield.begin(), forcefield.end(),
                                        forcefield.begin(), ::toupper);

    std::string engine_formatter = R"(
   Engine {}
     ForceField {}
   EndEngine
  )";
    engine_block = fmt::format(engine_formatter, engine, forcefield);
    return engine_block;
  } else if (engine == "FORCEFIELD") {
    std::string engine_formatter = R"(
   Engine {}
   EndEngine
  )";
    engine_block = fmt::format(engine_formatter, engine);
    return engine_block;
  }

  // Never reach here
  throw std::runtime_error("Generic AMS engine error \n");
}

/*
**
** DEBUGGER
**
** The set of functions below, when activated, include asserts to ensure
*equality with AMS_IO
* These work.
*/

// clang-format off

// void AMS::recieveFromSystem(long N, double *F, double *U) {

//   FILE *in;
//   double junkF;
//   char junkChar[256];
//   double forceX;
//   double forceY;
//   double forceZ;
//   double index;
//   char line[256];

//   in = fopen("ams_output", "r");

//   while (fgets(line, sizeof(line), in)) {

//     if (strcmp(line, "     CALCULATION RESULTS\n") ==
//         0) { // Finding the Energy in the output file

//       fscanf(in, "%s %s %s %lf", junkChar, junkChar, junkChar, U);
//       *U = *U * 27.2114; // Energy in hartree to eV
//     }
//     if (strcmp(line, "  Index   Atom            d/dx            d/dy           "
//                      " d/dz\n") == 0) { // Finding the forces
//       for (int i = 0; i < N; i++) {
//         fscanf(in, "%lf %s %lf %lf %lf", &index, &junkChar, &forceX, &forceY,
//                &forceZ);
//         F[int(i) * 3 + 0] = -forceX;
//         F[int(i) * 3 + 1] = -forceY;
//         F[int(i) * 3 + 2] = -forceZ;
//         // AMS gives gradients, not forces, hence the change.
//       }
//     }
//   }
//   for (int i = 0; i < 3 * N; i++) {
//     F[i] = F[i] * 51.4220862; // Forces from hartree/bohr to eV/Angstrom
//   }

//   fclose(in);
//   return;
// }

// void AMS::force(long N, const double *R, const int *atomicNrs, double *F,
//                 double *U, const double *box, int nImages = 1) {
//   passToSystem(N, R, atomicNrs, box);
//   system("chmod +x run_AMS.sh");
//   system("./run_AMS.sh >> ams_output"); // Run a single point AMS calculation
//                                         // and write the results into ams_output
//   recieveFromSystem(N, F, U);
//   runAMS();
//   double energ = extract_scalar_rkf("Energy"); // Sets energy
//   // assert(fabs(*U-this->energy)<DBL_EPSILON); // Equality == doesn't work well
//   // for floats
//   assert(fabs(*U - energ) < 1e-5); // Equality == doesn't work well for floats
//   // *U = this->energy;
//   // this->forces.clear();        // TODO: Slow!
//   std::vector<double> frc = extract_cartesian_rkf("Gradients"); // Sets forces
//   double *ftest = frc.data();
//   for (int i = 0; i < forces.size(); i++) {
//     assert(fabs(F[i] - ftest[i]) < 1e-5);
//   }
//   *U = energ;
//   F = ftest;
//   return;
// }

// Ordering which works, scopes are correct
// void AMS::force(long N, const double *R, const int *atomicNrs, double *F,
//                 double *U, const double *box, int nImages = 1) {
//   // Somehow broken, even though this is the same as before
//   passToSystem(N, R, atomicNrs, box);
//   runAMS();
//   std::vector<double> frc = extract_cartesian_rkf("Gradients");
//   double *ftest = frc.data();
//   for (int i = 0; i < N; i++) {
//     F[3 * i] = ftest[3 * i];
//     F[3 * i + 1] = ftest[3 * i + 1];
//     F[3 * i + 2] = ftest[3 * i + 2];
//   }
//   *U = extract_scalar_rkf("Energy");
//   return;
// }

// clang-format on

/*
** Debugging Toggles
** These functions can be used to validate the job ordering, by being added to
*runAMS()
*/

// TODO: Only in debug
// std::string AMS::readFile(std::filesystem::path path) {
//   // Kanged: https://stackoverflow.com/a/40903508/1895378
//   std::ifstream f(path, std::ios::in | std::ios::binary);
//   const auto sz = std::filesystem::file_size(path);
//   std::string stres(sz, '\0');
//   f.read(stres.data(), sz);
//   return stres;
// }

// bool AMS::validate_order() {
//   // Validate all inputs
//   if (not first_run) {
//     std::string rfile = readFile("run_AMS.sh");
//     std::string resfile = readFile("myrestart.in");
//     std::string updcoord = readFile("updCoord.sh");
//     // ifstream runfile("run_AMS.sh"), restartfile("myrestart.in"),
//     // updcoord("updCoord.sh"); istream_iterator<string> rfiter(runfile),
//     // refiter(restartfile), uciter(updcoord), eof; vector<string>
//     // rfstore(rfiter, eof), refstore(refiter, eof), ucstore(uciter, eof);
//     // (?<=AMS_JOBNAME=).*$
//     // (?<=udmpkf ).*(?=\.results)
//     // (?<=File ).*(?=\.results)
//     // The logic here is that the current job restarts from the previous one,
//     so
//     // the coordinates of the previous job are updated, the restart
//     references
//     // the previous job, and then finally the current job executes with cjob
//     assert(absl::StrContains(rfile, cjob));
//     assert(absl::StrContains(resfile, pjob));
//     assert(absl::StrContains(updcoord, pjob));
//   }
//   return true;
// }
