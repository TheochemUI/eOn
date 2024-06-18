#include "CommandLine.h"
#include "Matter.h"
#include "Parameters.h"
#include "Potential.h"
#include "version.h"

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

using namespace std;

void singlePoint(std::unique_ptr<Matter> matter) {
  std::cout << "Energy:         " << std::fixed << std::setprecision(15)
            << matter->getPotentialEnergy() << std::endl;
  std::cout << "(free) Forces:         \n" << matter->getForcesFree() << "\n";
  std::cout << "Max atom force: " << std::scientific << matter->maxForce()
            << std::endl;
}

void minimize(std::unique_ptr<Matter> matter, const string &confileout) {
  matter->relax(false, false);
  if (!confileout.empty()) {
    std::cout << "Saving relaxed structure to " << confileout << std::endl;
  } else {
    std::cout << "No output file specified, not saving" << std::endl;
  }
  matter->matter2con(confileout);
}

void commandLine(int argc, char **argv) {
  bool sflag = false, mflag = false, pflag = false, cflag = false;
  double optConvergedForce = 0.001;
  string potential;
  string confile;
  string optimizer("cg");

  auto params = std::make_shared<Parameters>();

  cxxopts::Options options("eonclient", "The eON client");
  options.add_options()("v,version", "Print version information")(
      "m,minimize", "Minimization of inputConfile saves to outputConfile")(
      "s,single", "Single point energy of inputConfile")(
      "c,compare", "Compare structures of inputConfile to outputConfile")(
      "o,optimizer", "Optimization method",
      cxxopts::value<std::string>()->default_value("cg"))(
      "f,force", "Convergence force",
      cxxopts::value<double>()->default_value("0.001"))(
      "t,tolerance", "Distance tolerance",
      cxxopts::value<double>()->default_value("0.1"))(
      "p,potential", "The potential (e.g. qsc, lj, eam_al)",
      cxxopts::value<std::string>())("h,help", "Print usage");

  try {
    auto result = options.parse(argc, argv);

    if (result.count("help")) {
      std::cout << options.help() << std::endl;
      exit(0);
    }

    if (result.count("version")) {
      std::cout << "eonclient version r" << VERSION << std::endl;
      std::cout << "          compiled " << BUILD_DATE << std::endl;
      exit(0);
    }

    if (result.count("minimize")) {
      mflag = true;
    }

    if (result.count("single")) {
      sflag = true;
    }

    if (result.count("compare")) {
      cflag = true;
    }

    if (result.count("potential")) {
      pflag = true;
      potential = result["potential"].as<std::string>();
    }

    if (result.count("optimizer")) {
      optimizer = result["optimizer"].as<std::string>();
    }

    if (result.count("force")) {
      optConvergedForce = result["force"].as<double>();
    }

    if (result.count("tolerance")) {
      params->distanceDifference = result["tolerance"].as<double>();
    }

    if (sflag && mflag) {
      std::cerr << "Cannot specify both minimization and single point"
                << std::endl;
      exit(2);
    }

    if (!pflag && (sflag || mflag)) {
      std::cerr << "Must specify a potential" << std::endl;
      exit(2);
    }

    if (!cflag) {
      for (auto &ch : potential) {
        ch = tolower(ch);
      }
    }

    auto unmatched = result.unmatched();
    if (unmatched.size() < 1) {
      std::cerr << "At least one non-option argument is required: the con file"
                << std::endl;
      exit(2);
    } else {
      confile = unmatched[0];
    }

    if (!cflag) {
      params->potential = magic_enum::enum_cast<PotType>(
                              potential, magic_enum::case_insensitive)
                              .value_or(PotType::UNKNOWN);
    }

    if (!sflag) {
      params->optMethod = magic_enum::enum_cast<OptType>(
                              optimizer, magic_enum::case_insensitive)
                              .value_or(OptType::CG);
      params->optConvergedForce = optConvergedForce;
    }

    auto pot = helper_functions::makePotential(params);
    auto matter = std::make_unique<Matter>(pot, params);
    auto matter2 = std::make_unique<Matter>(pot, params);
    matter->con2matter(confile);

    string confileout;
    if (unmatched.size() == 2) {
      confileout = unmatched[1];
      if (cflag)
        matter2->con2matter(confileout);
    }

    if (sflag) {
      singlePoint(std::move(matter));
    } else if (mflag) {
      minimize(std::move(matter), confileout);
    } else if (cflag) {
      params->checkRotation = true;
      if (matter->compare(*matter2, true)) {
        std::cout << "Structures match" << std::endl;
      } else {
        std::cout << "Structures do not match" << std::endl;
      }
    }
  } catch (const cxxopts::exceptions::exception &e) {
    std::cerr << "Error parsing options: " << e.what() << std::endl;
    std::cerr << options.help() << std::endl;
    exit(1);
  }
}
