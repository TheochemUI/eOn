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
#include <fstream>
#include <iostream>

namespace bp = boost::process;

AMS::AMS(Parameters *p) {
  this->engine_setup = generate_run(p);
  engine = p->engine.c_str();
  forcefield = p->forcefield.c_str();
  model = p->model.c_str();
  xc = p->xc.c_str();
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
  counter = 0;
  // TODO: Optimize and reuse existing files Currently each Matter will
  // recreate the folders It should instead figure out if results exist and
  // use them
  this->first_run = true;
  this->job_one = true;
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
  std::future<std::string> err;
  chmod("run_AMS.sh", S_IRWXU);
  nativenv["AMS_JOBNAME"] = cjob;
  bp::child c("run_AMS.sh", // set the input
              nativenv, bp::std_in.close(),
              bp::std_out > bp::null, // so it can be written without anything
              bp::std_err > err, amsRun);
  amsRun.run();
  // TODO: Rework
  // auto erro = err.get();
  // try {
  //   if (!absl::StrContains(erro, "NORMAL TERMINATION")) {
  //     // throw std::runtime_error("AMS STDERR:\n");
  //   } else {
  //     std::cout << "\nAMS exited with " << erro;
  //   }
  // } catch (const std::exception &e) {
  //   std::cout << e.what();
  //   std::cout << erro;
  //   exit(0);
  // }
}

void AMS::extract_rkf(long N, std::string key) {
  // The logic here is that the extraction must be asynchronous, as we are not
  // interested in processing this the output line by line, and would prefer to
  // have it all in one place
  std::string execString, strEngine(engine);
  std::vector<std::string> execDat, innerDat;
  boost::asio::io_context rkf;
  std::future<std::string> rkf_out_future, rkf_err_future;
  std::string rkfout, rkferr;
  int counter;
  double x;
  std::vector<double> extracted;
  std::transform(strEngine.begin(), strEngine.end(), strEngine.begin(),
                 ::tolower);
  execString =
      fmt::format("dmpkf {jobid:}.results/{engine:}.rkf AMSResults%{key:}",
                  fmt::arg("jobid", cjob), fmt::arg("engine", strEngine),
                  fmt::arg("key", key));
  std::cout << execString << "\n";

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
    std::cout << "Extracting " << key << std::endl;
  }

  execDat = absl::StrSplit(rkfout, '\n');

  if (N == 1) {
    // Is energy or some other scalar property
    // First 3 lines are the headers
    // [0] = "AMSResults                      "
    // [1] = "Energy                          "
    // [2] = "         1         1         2"
    // [3] = "   0.135012547958282714E+000"
    // [4] = ""
    if (absl::SimpleAtod(execDat[3], &x)) {
      this->energy = x * energyConversion;
      std::cout << fmt::format(
          "\n Got {:.4f} Hartree from AMS and converted to {:.4f} eV\n", x,
          energy);
      return; // Early return
    } else {
      throw std::runtime_error(
          fmt::format("\n Expected energy, got {} instead", execDat[3]));
    }
  } else {
    // Assume gradients or some other x y z property
    // There is one per atom, the number of which is N
    // The first three lines are the header files as before
    for (int i = 3; i < N; i++) {
      std::vector<std::string> strrow = absl::StrSplit(execDat[i], ' ');
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
      for (auto elem : strrow) {
        if (absl::SimpleAtod(elem, &x)) {
          double force = x * forceConversion;
          std::cout << fmt::format(
              "\n Gradient element={grad:.4f} Hartree/Bohr from AMS\n Force "
              "element={force:.4f} eV/Angstrom\n",
              fmt::arg("grad", x), fmt::arg("force", force));
          extracted.emplace_back(x);
        }
      }
    }
    this->forces = extracted;
    return; // Early return
  }
  // Never reach here
  throw std::runtime_error("Generic AMS dmpkf error \n");
}

void AMS::updateCoord(long N, const double *R) {
  // TODO: Probably don't need to read in just for the first three lines
  std::ofstream updCoord;
  std::string execString, coordDump, newCoord;
  std::vector<std::string> execDat;
  boost::asio::io_context coordio;
  std::future<std::string> err, rdump;
  int counter;
  std::vector<double> gradients;
  // Prep new run
  absl::StrAppend(&execString, "dmpkf ", cjob,
                  ".results/ams.rkf Molecule%Coords");
  // std::cout << execString << "\n";
  // Store Coordinates
  bp::child cprog(execString, nativenv, bp::std_in.close(), bp::std_out > rdump,
                  bp::std_err > err, coordio);
  coordio.run();
  execDat = absl::StrSplit(rdump.get(), '\n');
  // Get the first three lines
  counter = 0;
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

void AMS::force(long N, const double *R, const int *atomicNrs, double *F,
                double *U, const double *box, int nImages = 1) {
  if (job_one) {
    // True, create first directory
    cjob = "firstRun";
    pjob = "secondRun";
  } else {
    // False, create second directory
    cjob = "secondRun";
    pjob = "firstRun";
  }

  if (first_run) {
    // First run
    restartFrom.open("myrestart.in");
    restartFrom << "\n";
    restartFrom.close();
    restartj.clear();
    passToSystem(N, R, atomicNrs, box);
    runAMS();
  } else {
    updateCoord(N, R);
    runAMS();
  }
  energy = 0;
  extract_rkf(1, "Energy"); // Sets energy
  *U = energy;
  // std::cout<<"Got an energy of "<<*U<<"\n";
  forces.clear();              // TODO: Slow!
  extract_rkf(N, "Gradients"); // Sets forces
  for (int i = 0; i < N; i++) {
    F[3 * i] = forces[3 * i];
    F[3 * i + 1] = forces[3 * i + 1];
    F[3 * i + 2] = forces[3 * i + 2];
  }
  // Toggle job name
  job_one = !job_one;
  // bp::spawn("cat myrestart.in");
  std::string strEngine(engine);
  std::transform(strEngine.begin(), strEngine.end(), strEngine.begin(),
                 ::tolower);
  // TODO: This is a hacky workaround, if there is mopac or an unsupported
  // potential, do not include the restart information
  if (strEngine == "mopac") {
    restartj = "\n";
  } else {
    restartj = "EngineRestart ";
    absl::StrAppend(&restartj, cjob, ".results/", strEngine, ".rkf\n");
  }
  absl::StrAppend(&restartj, "LoadSystem\nFile ", cjob,
                  ".results/ams.rkf\nSection Molecule\nEnd");
  restartFrom.open("myrestart.in");
  restartFrom << restartj;
  restartFrom.close();
  restartj.clear();
  if (first_run) {
    smallSys(N, R, atomicNrs, box);
  }
  // Falsify forever
  first_run = false;
  return;
}

void AMS::passToSystem(long N, const double *R, const int *atomicNrs,
                       const double *box)
// Creating the standard input file that will be read by the AMS driver
{
  FILE *out;
  out = fopen("run_AMS.sh", "w");

  fprintf(out, "#!/bin/sh\n");
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
  if (strlen(model) > 0 || strlen(forcefield)) {
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
  fprintf(out, "@include myrestart.in\n");
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

std::string AMS::generate_run(Parameters *p) {
  // TODO: Use args everywhere, cleaner logic
  std::string engine_block;
  // Get the values from the configuration
  std::string engine(p->engine), forcefield(p->forcefield), model(p->model),
      xc(p->xc), resources(p->resources), basis(p->basis);
  // Ensure capitals and existence
  engine.empty()
      ? throw std::runtime_error("AMS Engine is required \n")
      : std::transform(engine.begin(), engine.end(), engine.begin(), ::toupper);

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
  } else if (engine == "REAXFF") {
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
