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
    this->curpath = obspath;
    this->gprfunc.setHyperparameters(this->curpath, this->atmconf);
    this->gprfunc.optimize(this->curpath);
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

void GPRMatter::retrainGPR(gpr::Observation& newpath){
    this->gprfunc.setHyperparameters(newpath, this->atmconf);
    this->gprfunc.optimize(newpath);
    this->trainedGPR->registerGPRObject(&this->gprfunc);
}

GPRobj::~GPRobj(){

};

GPRobj::GPRobj(Matter initMatter, Parameters eonp): eonp{eonp}, trainedGPR{&eonp} {
    // Setup AtomsConf
    auto conf_data = helper_functions::eon_matter_to_frozen_conf_info(&initMatter,
                                                                      eonp.gprPotActiveRadius);
    this->atmconf = std::get<gpr::AtomsConfiguration>(conf_data);
    gpr::GPRSetup gpr_parameters;
    gpr_parameters.jitter_sigma2 = eonp.gprPotJitterSigma2;
    gpr_parameters.sigma2 = eonp.gprPotSigma2;
    this->gprfunc.setParameters(gpr_parameters);
    auto  potparams = helper_functions::eon_parameters_to_gprpot(&eonp);
    for (int i = 0; i < 9; i++) {
        potparams.cell_dimensions.value[i] = initMatter.getCell()(i);
    }
    gprfunc.initialize(potparams, this->atmconf);
    // TODO: Figure out if the first point should be added to curpath
}

void GPRobj::trainGPR(std::vector<Matter>& initialPoints){
    Potential* potter = Potential::getPotential(&this->eonp);
    this->curpath.clear();
    this->curpath = GPRobj::prepobs(initialPoints, potter);
    this->gprfunc.setHyperparameters(this->curpath, this->atmconf);
    this->gprfunc.optimize(this->curpath);
}

void GPRobj::retrainGPR(std::vector<Matter>& newPoints){
    Potential* potter = Potential::getPotential(&this->eonp);
    this->curpath = GPRobj::prepobs(newPoints, potter);
    this->gprfunc.setHyperparameters(this->curpath, this->atmconf, false);
    this->gprfunc.optimize(this->curpath);
}

gpr::Observation GPRobj::prepobs(std::vector<Matter>& matvec, Potential* pot){
    size_t nfree(matvec.front().numberOfFreeAtoms());
    gpr::Observation obspath{this->curpath}, obs;
    for (auto& mat : matvec){
        obs.clear();
        obs.R.resize(1, nfree * 3);
        obs.G.resize(1, nfree * 3);
        obs.E.resize(1);
        auto pe_forces = helper_functions::energy_and_forces(&mat, pot);
        obs.E.set(std::get<double>(pe_forces));
        auto forces = std::get<AtomMatrix>(pe_forces);
        AtomMatrix freepos = mat.getPositionsFree();
        for (size_t idx{0}; idx < nfree * 3; ++idx){
            obs.R[idx] = freepos.data()[idx];
            obs.G[idx] = -1 * forces.data()[idx];
        }
        obspath.append(obs);
    }
    return obspath;
}
