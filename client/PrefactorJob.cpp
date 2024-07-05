#include "PrefactorJob.h"
#include "Hessian.h"
#include "Matter.h"
#include "Potential.h"
#include "Prefactor.h"

std::vector<std::string> PrefactorJob::run(void) {
  std::vector<std::string> returnFiles;
  VectorType freqs;

  std::string reactantFilename("reactant.con");
  std::string saddleFilename("saddle.con");
  std::string productFilename("product.con");

  auto reactant = std::make_unique<Matter>(pot, params);
  auto saddle = std::make_unique<Matter>(pot, params);
  auto product = std::make_unique<Matter>(pot, params);

  reactant->con2matter("reactant.con");
  saddle->con2matter("saddle.con");
  product->con2matter("product.con");
  double pref1, pref2;
  Prefactor::getPrefactors(params.get(), reactant.get(), saddle.get(),
                           product.get(), pref1, pref2);
  // printf("pref1: %.3e pref2: %.3e\n", pref1, pref2);

  Vector<int> atoms;
  if (params->prefactor.allFreeAtoms) {
    // it is sufficient to pass the configuration
    // for which the frequencies should be determined
    std::string matterFilename;
    // TODO(rg): Also an enum candidate
    if (params->prefactor.configuration == ::Prefactor::TYPE::REACTANT) {
      matterFilename = reactantFilename;
    } else if (params->prefactor.configuration == ::Prefactor::TYPE::SADDLE) {
      matterFilename = saddleFilename;
    } else if (params->prefactor.configuration == ::Prefactor::TYPE::PRODUCT) {
      matterFilename = productFilename;
    }
    reactant->con2matter(matterFilename);
    saddle->con2matter(matterFilename);
    product->con2matter(matterFilename);

    // account for all free atoms
    atoms = Prefactor::allFreeAtoms(reactant.get());
  } else {
    reactant->con2matter(reactantFilename);
    saddle->con2matter(saddleFilename);
    product->con2matter(productFilename);

    // determine which atoms moved in the process
    atoms = Prefactor::movedAtoms(params.get(), reactant.get(), saddle.get(),
                                  product.get());
  }
  assert(3 * atoms.rows() > 0);

  // calculate frequencies
  if (params->prefactor.configuration == ::Prefactor::TYPE::REACTANT) {
    Hessian hessian(params.get(), reactant.get());
    freqs = hessian.getFreqs(reactant.get(), atoms);
  } else if (params->prefactor.configuration == ::Prefactor::TYPE::SADDLE) {
    Hessian hessian(params.get(), saddle.get());
    freqs = hessian.getFreqs(saddle.get(), atoms);
  } else if (params->prefactor.configuration == ::Prefactor::TYPE::PRODUCT) {
    Hessian hessian(params.get(), product.get());
    freqs = hessian.getFreqs(product.get(), atoms);
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

  return returnFiles;
}
