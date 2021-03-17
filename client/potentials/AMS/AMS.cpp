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
#include <unistd.h>

namespace bp = boost::process;

AMS::AMS(Parameters *p) {
  engine = p->engine.c_str();
  forcefield = p->forcefield.c_str();
  model = p->model.c_str();
  xc = p->xc.c_str();
  counter = 0;
  // TODO: Optimize and reuse existing files Currently each Matter will
  // recreate the folders It should instead figure out if results exist and
  // use them
  first_run = true;
  job_one = true;
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
  bp::spawn("chmod +x run_AMS.sh");
  bp::child c("run_AMS.sh", // set the input
              bp::env["AMS_JOBNAME"] = jname, bp::std_in.close(),
              bp::std_out > bp::null, // so it can be written without anything
              bp::std_err > err, amsRun);
  amsRun.run();
  auto erro = err.get();
  try {
    if (!absl::StrContains(erro, "NORMAL TERMINATION")) {
      throw std::runtime_error("AMS STDERR:\n");
    } else {
      // std::cout << "\nAMS exited with " << erro;
    }
  } catch (const std::exception &e) {
    std::cout << e.what();
    std::cout << erro;
  }
}

void AMS::extract_rkf(long N, std::string key) {
  std::string execString, strEngine(engine);
  std::vector<std::string> execDat, innerDat;
  boost::asio::io_context rkf;
  std::future<std::string> err, rdump;
  int counter;
  double x;
  std::vector<double> extracted;
  std::transform(strEngine.begin(), strEngine.end(), strEngine.begin(),
                 ::tolower);
  absl::StrAppend(&execString, "dmpkf ", jname, ".results/", strEngine,
                  ".rkf AMSResults%", key);
  // std::cout << execString << "\n";
  // Extract
  bp::child eprog(execString, bp::std_in.close(), bp::std_out > rdump,
                  bp::std_err > err, rkf);

  rkf.run();
  auto erro = err.get();
  try {
    if (!absl::StrContains(erro, "NORMAL TERMINATION")) {
      throw std::runtime_error("AMS STDERR:\n");
    } else {
      // std::cout << "Extract " << key << ": " << erro;
    }
  } catch (const std::exception &e) {
    std::cout << e.what();
    std::cout << erro;
  }
  execDat = absl::StrSplit(rdump.get(), '\n');

  if (N == 1) {
    // Is energy or some other scalar property
    counter = 0;
    for (auto i : execDat) {
      if (counter >= 3) {
        if (absl::SimpleAtod(i, &x)) {
          // std::cout << x << "\n";
          x = x * energyConversion;
          energy = x;
        }
      }
      counter++;
    }
  } else {
    // Assume gradients or some other x y z property
    counter = 0;
    for (auto j : execDat) {
      if (counter >= 3) {
        innerDat = absl::StrSplit(j, ' ');
        for (auto k : innerDat) {
          if (absl::SimpleAtod(k, &x)) {
            x = x * forceConversion;
            // std::cout << x << " ";
            extracted.push_back(x);
          }
        }
        // std::cout << "\n";
      }
      counter++;
    }
    forces = extracted;
  }
  return;
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
  absl::StrAppend(&execString, "dmpkf ", jname,
                  ".results/ams.rkf Molecule%Coords");
  // std::cout << execString << "\n";
  // Store Coordinates
  bp::child cprog(execString, bp::std_in.close(), bp::std_out > rdump,
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
    absl::StrAppend(&newCoord, R[a], "  ");
    counter++;
    if (counter % 3 == 0) {
      absl::StrAppend(&newCoord, "\n");
    }
  }
  coordDump = "#!/bin/sh\n udmpkf ";
  absl::StrAppend(&coordDump, jname, ".results/ams.rkf <<EOF\n", newCoord,
                  "EOF");
  // std::cout << coordDump;
  updCoord.open("updCoord.sh");
  updCoord << coordDump;
  updCoord.close();
  bp::spawn("chmod +x updCoord.sh");
  bp::child cuprog("updCoord.sh", bp::std_err > bp::null);
  cuprog.wait();
  return;
}

void AMS::force(long N, const double *R, const int *atomicNrs, double *F,
                double *U, const double *box, int nImages = 1) {
  if (job_one) {
    // True, create first directory
    jname = "firstRun";
  } else {
    // False, create second directory
    jname = "secondRun";
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
  std::cout<<"Got an energy of "<<*U<<"\n";
  forces.clear(); // TODO: Slow!
  extract_rkf(N, "Gradients"); // Sets forces
  counter = 0;
  for (double f : forces) {
    F[counter] = f;
    counter++;
  }
  // Toggle job name
  job_one = !job_one;
  // bp::spawn("cat myrestart.in");
  restartj = "EngineRestart ";
  absl::StrAppend(&restartj, jname, ".results/reaxff.rkf\n");
  absl::StrAppend(&restartj, "LoadSystem\nFile ", jname,
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
  fprintf(out, "Engine %s\n", engine);
  if (strlen(forcefield) > 0) {
    fprintf(out, "     Forcefield %s\n", forcefield);
  }
  if (strlen(model) > 0) {
    fprintf(out, "     Model %s\n", model);
  }
  if (strlen(xc) > 0) {
    fprintf(out, "xc %s\n");
    fprintf(out, "     hybrid %s\n",
            xc); // basis set not specified (default = DZ)
    fprintf(out, "end\n");
  }
  fprintf(out, "EndEngine\n");
  fprintf(out, "Properties\n");
  fprintf(out, " Gradients\n");
  fprintf(out, "End\n");
  fprintf(out, "@include myrestart.in\n");
  fprintf(out, "eor");
  fclose(out);
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
  fprintf(out, "Engine %s\n", engine);
  if (strlen(forcefield) > 0) {
    fprintf(out, "     Forcefield %s\n", forcefield);
  }
  if (strlen(model) > 0) {
    fprintf(out, "     Model %s\n", model);
  }
  if (strlen(xc) > 0) {
    fprintf(out, "xc %s\n");
    fprintf(out, "     hybrid %s\n",
            xc); // basis set not specified (default = DZ)
    fprintf(out, "end\n");
  }
  fprintf(out, "EndEngine\n");
  fprintf(out, "Properties\n");
  fprintf(out, " Gradients\n");
  fprintf(out, "End\n");
  fprintf(out, "@include myrestart.in\n");
  fprintf(out, "eor");
  fclose(out);
  return;
}
