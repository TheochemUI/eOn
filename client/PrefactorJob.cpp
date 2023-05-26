#include "PrefactorJob.h"
#include "Hessian.h"
#include "Matter.h"
#include "Potential.h"
#include "Prefactor.h"

const char PrefactorJob::PREFACTOR_REACTANT[] = "reactant";
const char PrefactorJob::PREFACTOR_SADDLE[] = "saddle";
const char PrefactorJob::PREFACTOR_PRODUCT[] = "product";

std::vector<std::string> PrefactorJob::run(void) {
  std::vector<std::string> returnFiles;
  VectorXd freqs;

  string reactantFilename("reactant.con");
  string saddleFilename("saddle.con");
  string productFilename("product.con");

  Matter *reactant = new Matter(params);
  Matter *saddle = new Matter(params);
  Matter *product = new Matter(params);

  reactant->con2matter("reactant.con");
  saddle->con2matter("saddle.con");
  product->con2matter("product.con");
  double pref1, pref2;
  Prefactor::getPrefactors(params.get(), reactant, saddle, product, pref1,
                           pref2);
  // printf("pref1: %.3e pref2: %.3e\n", pref1, pref2);

  VectorXi atoms;
  if (params->prefactorAllFreeAtoms) {
    // it is sufficient to pass the configuration
    // for which the frequencies should be determined
    string matterFilename;
    if (params->prefactorConfiguration == PrefactorJob::PREFACTOR_REACTANT) {
      matterFilename = reactantFilename;
    } else if (params->prefactorConfiguration ==
               PrefactorJob::PREFACTOR_SADDLE) {
      matterFilename = saddleFilename;
    } else if (params->prefactorConfiguration ==
               PrefactorJob::PREFACTOR_PRODUCT) {
      matterFilename = productFilename;
    }
    reactant->con2matter(matterFilename);
    saddle->con2matter(matterFilename);
    product->con2matter(matterFilename);

    // account for all free atoms
    atoms = Prefactor::allFreeAtoms(reactant);
  } else {
    reactant->con2matter(reactantFilename);
    saddle->con2matter(saddleFilename);
    product->con2matter(productFilename);

    // determine which atoms moved in the process
    atoms = Prefactor::movedAtoms(params.get(), reactant, saddle, product);
  }
  assert(3 * atoms.rows() > 0);

  // calculate frequencies
  if (params->prefactorConfiguration == PrefactorJob::PREFACTOR_REACTANT) {
    Hessian hessian(params.get(), reactant);
    freqs = hessian.getFreqs(reactant, atoms);
  } else if (params->prefactorConfiguration == PrefactorJob::PREFACTOR_SADDLE) {
    Hessian hessian(params.get(), saddle);
    freqs = hessian.getFreqs(saddle, atoms);
  } else if (params->prefactorConfiguration ==
             PrefactorJob::PREFACTOR_PRODUCT) {
    Hessian hessian(params.get(), product);
    freqs = hessian.getFreqs(product, atoms);
  }

  bool failed = freqs.size() != 3 * atoms.rows();

  FILE *fileResults;
  FILE *fileFreq;

  std::string results_file("results.dat");
  std::string freq_file("freq.dat");

  returnFiles.push_back(results_file);
  returnFiles.push_back(freq_file);

  fileResults = fopen(results_file.c_str(), "wb");
  fileFreq = fopen(freq_file.c_str(), "wb");

  fprintf(fileResults, "%s good\n", failed ? "false" : "true");
  // fprintf(fileResults, "%d force_calls\n", Potential::fcalls);

  if (!failed) {
    for (int i = 0; i < freqs.size(); i++) {
      if (0. < freqs[i]) {
        fprintf(fileFreq, "%f\n", sqrt(freqs[i]) / (2 * M_PI * 10.18e-15));
      } else {
        fprintf(fileFreq, "%f\n", -sqrt(-freqs[i]) / (2 * M_PI * 10.18e-15));
      }
    }
  }

  delete reactant;
  delete product;
  delete saddle;

  return returnFiles;
}
