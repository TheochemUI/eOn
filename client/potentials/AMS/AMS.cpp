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

void AMS::force(long N, const double *R, const int *atomicNrs, double *F,
                double *U, const double *box, int nImages = 1) {
  double forceConversion =
      -51.4220862; // Forces from hartree/bohr to eV/Angstrom
  std::string tmpString, tmpExec;
  passToSystem(N, R, atomicNrs, box);
  bp::spawn("chmod +x run_AMS.sh");
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
  }

  bp::child c("run_AMS.sh", // set the input
              bp::env["AMS_JOBNAME"] = jname, bp::std_in.close(),
              bp::std_out > bp::null, // so it can be written without anything
              bp::std_err > err, ams_run);
  ams_run.run();
  auto erro = err.get();
  try {
    if (!absl::StrContains(erro, "NORMAL TERMINATION")) {
      throw std::runtime_error("AMS STDERR:\n");
    } else {
      std::cout << "\nExited with " << erro;
    }
  } catch (const std::exception &e) {
    std::cout << e.what();
    std::cout << erro;
  }

  absl::StrAppend(&tmpString, "dmpkf ", jname, ".results/reaxff.rkf AMSResults%");
  absl::StrAppend(&tmpExec, tmpString, "Energy");
  std::cout<<tmpExec<<"\n";
  // Extract energy
  bp::child eprog(tmpExec,
                  bp::std_in.close(), bp::std_out > edump, bp::std_err > erre,
                  ams_rkf);

  tmpExec.clear();
  absl::StrAppend(&tmpExec, tmpString, "Gradients");
  std::cout<<tmpExec<<"\n";

  // Extract gradients
  bp::child gprog(tmpExec,
                  bp::std_in.close(), bp::std_out > gdump, bp::std_err > errg,
                  ams_rkf);
  ams_rkf.run();

  energ = absl::StrSplit(edump.get(), '\n');
  grad = absl::StrSplit(gdump.get(), '\n');

  counter = 0;
  for (auto i : energ) {
    if (counter >= 3) {
      if (absl::SimpleAtod(i, &x)) {
        // std::cout << x << "\n";
        energy = x;
      }
    }
    counter++;
  }

  counter = 0;
  for (auto j : grad) {
    if (counter >= 3) {
      tmp = absl::StrSplit(j, ' ');
      for (auto k : tmp) {
        if (absl::SimpleAtod(k, &x)) {
          x = x*forceConversion;
          // std::cout << x << " ";
          gradients.push_back(x);
        }
      }
      // std::cout << "\n";
    }
    counter++;
  }

  // Prep new run
  ams_run.restart();
  // Store Coordinates
  bp::child cprog(
      "dmpkf ams.results/ams.rkf Molecule%Coords",
      bp::std_in.close(),
      bp::std_out > cdump,
      bp::std_err > errc, ams_run);

  ams_run.run();

  std::vector<std::string> coord = absl::StrSplit(cdump.get(), '\n');
  std::string newCoord;

  // Get the first three lines
  counter = 0;
  for (auto j : coord) {
    if (counter >= 3) {
      break;
    } else {
      absl::StrAppend(&newCoord, j, "\n");
    }
    counter++;
  }

  // Put the rest of the coordinates
  counter = 0;
  for (int a = 0; a < N; a++) {
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
  bp::child cuprog("updCoord.sh");
  cuprog.wait();
  // Toggle job name
  job_one = !job_one;
  // Falsify forever
  first_run = false;
  // Reset the io_contexts
  ams_run.restart();
  ams_rkf.restart();
  // bp::spawn("cat myrestart.in");

  restartj = "EngineRestart ";
  absl::StrAppend(&restartj, jname, ".results/reaxff.rkf");
  restartFrom.open("myrestart.in");
  restartFrom << restartj;
  restartFrom.close();
  restartj.clear();
  // std::transform(gradients.begin(), gradients.end(), gradients.begin(),
  //                [&forceConversion](auto &c) { return c * forceConversion; });
  F = gradients.data();
  energy *= 27.2114; // Energy in hartree to eV
  U = &energy;
  std::cout << "energy is " << *U << "\n";
  // for (int m = 0; m < N; m++) {
  //   std::cout << F[m] << " ";
  // }
  return;
}

void AMS::passToSystem(long N, const double *R, const int *atomicNrs,
                       const double *box)
// Creating the standard input file that will be read by the AMS driver
{
  FILE *out;
  out = fopen("run_AMS.sh", "w");

  fprintf(out, "#!/bin/sh\n");
  fprintf(out, "export NSCM_AMSEXTERNAL=$NSCM\n");
  fprintf(out, "export NSCM=1\n");
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
