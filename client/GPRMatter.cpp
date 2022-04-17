#include "GPRMatter.h"

GPRMatter::~GPRMatter(){

};

GPRMatter::GPRMatter(Matter initMatter) : truePotMatter{initMatter}  {
    this->eonp = *this->truePotMatter.getParameters();
    // Setup AtomsConf
    auto conf_data = helper_functions::eon_matter_to_frozen_conf_info(&this->truePotMatter,
                                                                      eonp.gprPotActiveRadius);
    this->atmconf = std::get<gpr::AtomsConfiguration>(conf_data);
    this->gpr_parameters.jitter_sigma2 = eonp.gprPotJitterSigma2;
    this->gpr_parameters.sigma2 = eonp.gprPotSigma2;
    this->gprfunc.setParameters(gpr_parameters);
    auto  potparams = helper_functions::eon_parameters_to_gprpot(&eonp);
    for (int i = 0; i < 9; i++) {
        potparams.cell_dimensions.value[i] = this->truePotMatter.getCell()(i);
    }
    this->trainedGPR = std::make_unique<GPRPotential>(&this->eonp);
    gprfunc.initialize(potparams, this->atmconf);
}

void GPRMatter::trainGPR(gpr::Observation obspath){
    this->gprfunc.setHyperparameters(obspath, this->atmconf);
    this->gprfunc.optimize(obspath);
    this->trainedGPR->registerGPRObject(&this->gprfunc);
}

std::pair<double, AtomMatrix> GPRMatter::gpr_energy_forces(Matter& getat){
  int nFreeAtoms = getat.numberOfFreeAtoms();
  auto posdata = getat.getPositionsFree();
  auto celldat = getat.getCell();
  AtomMatrix forces = AtomMatrix::Constant(nFreeAtoms, 3, 0);
  double *pos = posdata.data();
  double *frcs = forces.data();
  double *bx = celldat.data();
  double energy{0};
  double *erg = &energy;
  this->trainedGPR->force(nFreeAtoms, pos, nullptr, frcs, erg, bx, 1);
  return std::make_pair(energy, forces);
}

std::pair<double, AtomMatrix> GPRMatter::true_free_energy_forces(Matter& getat){
  return std::make_pair(getat.getPotentialEnergy(), getat.getForcesFree());
}

bool GPRMatter::isCloseTo(Matter& testat, double eps){
    // TODO: Use more convergence methods
    double convval {((gpr_energy_forces(testat)).second
    - testat.getForcesFree()).norm()};
    return (convval < eps);
}
