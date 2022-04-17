#include "GPRMatter.h"

GPRMatter::~GPRMatter(){

};

GPRMatter::GPRMatter(Matter initMatter, std::shared_ptr<GPRobj> gpf) : truePotMatter{initMatter}, gprobj{gpf}  {
}

std::pair<double, AtomMatrix> GPRMatter::gpr_energy_forces(){
  int nFreeAtoms = truePotMatter.numberOfFreeAtoms();
  auto posdata = truePotMatter.getPositionsFree();
  auto celldat = truePotMatter.getCell();
  AtomMatrix forces = AtomMatrix::Constant(nFreeAtoms, 3, 0);
  double *pos = posdata.data();
  double *frcs = forces.data();
  double *bx = celldat.data();
  double energy{0};
  double *erg = &energy;
  auto trainedGPR = this->gprobj->yieldGPRPot();
  trainedGPR.force(nFreeAtoms, pos, nullptr, frcs, erg, bx, 1);
  return std::make_pair(energy, forces);
}

std::pair<double, AtomMatrix> GPRMatter::true_free_energy_forces(){
  return std::make_pair(truePotMatter.getPotentialEnergy(),
                        truePotMatter.getForcesFree());
}

bool GPRMatter::isCloseToTrue(double eps){
    // TODO: Use more convergence methods
    double convval {((gpr_energy_forces()).second
    - truePotMatter.getForcesFree()).norm()};
    return (convval < eps);
}

void GPRMatter::updateMatter(const Matter otherMatter){
    this->truePotMatter = otherMatter;
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

GPRPotential GPRobj::yieldGPRPot(){
    GPRPotential gppot{&eonp};
    gppot.registerGPRObject(&gprfunc);
    return gppot;
}
