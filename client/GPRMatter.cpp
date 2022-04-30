#include "GPRMatter.h"
#include <cmath>

#include "NudgedElasticBand.h"

GPRMatter::~GPRMatter(){

};

GPRMatter::GPRMatter(Matter& initMatter, std::shared_ptr<GPRobj> gpf) :
   truePotMatter{initMatter},
   trueForcesFree{initMatter.getForcesFree()},
   truePotEnergy{initMatter.getPotentialEnergy()},
   trueForces{initMatter.getForces()},
   gprobj{gpf} {
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
  // TODO: Caching
  return std::make_pair(truePotMatter.getPotentialEnergy(),
                        truePotMatter.getForcesFree());
}

std::pair<double, AtomMatrix> GPRMatter::true_energy_forces(){
  // TODO: Caching
  return std::make_pair(truePotMatter.getPotentialEnergy(),
                        truePotMatter.getForces());
}

bool GPRMatter::areForcesCloseToTrue(double eps){
    // TODO: Use more convergence methods
    auto gprforces = (gpr_energy_forces()).second;
    auto trueforces = (true_free_energy_forces()).second;
    double convval {(gprforces - trueforces).norm()};
    // std::cout<<"Force norm diff "<<convval<<std::endl;
    return (convval < eps);
}

bool GPRMatter::isForceMaxElementLow(double eps){
  auto gprforces = (gpr_energy_forces()).second;
  std::vector<double> gpfvec(gprforces.size());
  Eigen::Map<AtomMatrix>(gpfvec.data(),
                         gprforces.rows(),
                         gprforces.cols()) = gprforces;
  auto itr_res = std::max_element(gpfvec.begin(),
                                  gpfvec.end(),
                                  helper_functions::abs_compare<double>);
  auto max_force_elem = std::abs(*itr_res);
  // std::cout<<"Max force element is "<<max_force_elem<<"\n";
  return (max_force_elem < eps);
}

bool GPRMatter::areEnergiesCloseToTrue(double eps){
    // TODO: Use more convergence methods
    double gpr_energy = std::abs((gpr_energy_forces()).first);
    // std::cout<<"GPR Energy "<<gpr_energy<<std::endl;
    double true_energy = std::abs((true_free_energy_forces()).first);
    // std::cout<<"True Energy "<<true_energy<<std::endl;
    double convval {std::abs(true_energy-gpr_energy)};
    // std::cout<<"Energy diff "<<convval<<std::endl;
    return (convval < eps);
}

void GPRMatter::updateMatter(Matter& otherMatter){
    this->truePotMatter = otherMatter;
    this->trueForcesFree = otherMatter.getForcesFree();
    this->trueForces = otherMatter.getForces();
    this->truePotEnergy = otherMatter.getPotentialEnergy();
}

GPRobj::~GPRobj(){

};

GPRobj::GPRobj(Matter initMatter, Parameters eonp): eonp{eonp} {
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

bool helper_functions::maybeUpdateGPRobj(NudgedElasticBand& neb, std::shared_ptr<GPRobj>& gpf){
  size_t totalImages(neb.images+2);
  bool result{false};
  for (size_t idx{0}; idx < totalImages; idx++){
    GPRMatter gp_testp(*neb.image[idx], gpf);
    if (gp_testp.areEnergiesCloseToTrue() && gp_testp.areForcesCloseToTrue()){
      result = true;
    } else{
      result = false;
    }
  }
  if (result){
    std::vector<Matter> ia2;
    for (size_t idx{1}; idx < totalImages - 1; idx++){
      ia2.push_back(*neb.image[idx]);
    }
    std::cout<<"Retraining";
    gpf->retrainGPR(ia2);
  }
  return result;
}

bool helper_functions::max_true_energy(GPRMatter& matOne, GPRMatter& matTwo){
      auto pef_one = matOne.true_energy_forces();
      auto pef_two = matTwo.true_energy_forces();
      auto eOne = std::get<double>(pef_one);
      auto eTwo = std::get<double>(pef_two);
      return (eOne < eTwo);
}

bool helper_functions::max_gpr_energy(GPRMatter& matOne, GPRMatter& matTwo){
      auto pef_one = matOne.gpr_energy_forces();
      auto pef_two = matTwo.gpr_energy_forces();
      auto eOne = std::get<double>(pef_one);
      auto eTwo = std::get<double>(pef_two);
      return (eOne < eTwo);
}

bool helper_functions::true_force_norm_converged(GPRMatter& mat, double eps){
  auto pef = mat.true_free_energy_forces();
  auto force_norm = std::get<AtomMatrix>(pef).norm();
  std::cout<<"\n truenorm: "<<force_norm<<" vs eps: "<<eps<<std::endl;
  return (force_norm < eps);
}

std::vector<GPRMatter> helper_functions::prepGPRMatterVec(std::vector<Matter>& newPath, std::shared_ptr<GPRobj>& gpf){
  std::vector<GPRMatter> gpmvec;
  std::transform(newPath.begin(),
                 newPath.end(),
                 std::back_inserter(gpmvec),
                 [&](Matter mat)->GPRMatter{
                   GPRMatter tmpmatter(mat, gpf);
                   return tmpmatter;
                 });
  return gpmvec;
}

void helper_functions::peSliceSurface(const std::pair<double, double> xrange,
                                      const std::pair<double, double> yrange,
                                      const double zlvl,
                                      const std::pair<size_t, size_t> elemXY,
                                      const Matter baseObject,
                                      const std::shared_ptr<GPRobj> gpf,
                                      const size_t gpr_num){
  auto fname = fmt::format("gpr_surface_{}.dat", gpr_num);
  auto oname = fmt::output_file(fname);
  std::vector<double> xvals(elemXY.first), yvals(elemXY.second);
  xvals.front() = xrange.first;
  xvals.back() = xrange.second;
  yvals.front() = yrange.first;
  yvals.back() = yrange.second;
  double xdelta = (xrange.second - xrange.first) / elemXY.first;
  double ydelta = (yrange.second - yrange.first) / elemXY.second;
  double xtmp{xrange.first}, ytmp{yrange.first};
  std::generate(xvals.begin(), xvals.end(), [&]{ xtmp = xtmp + xdelta; return xtmp; });
  std::generate(yvals.begin(), yvals.end(), [&]{ ytmp = ytmp + ydelta; return ytmp; });
  // Temporaries
  Matter mtmp = baseObject;
  oname.print(fmt::format("{:^15} {:^15} {:^15} {:^15} {:^15}\n",
                          "xval", "yval", "zval", "gpr_energy", "true_energy"));
  for (size_t xtick{0}; xtick < elemXY.first; xtick++){
    for (size_t ytick{0}; ytick < elemXY.second; ytick++){
      Vector3d pos{xvals[xtick], yvals[ytick], zlvl};
      mtmp.setPositionsFreeV(pos);
      GPRMatter gpmat(mtmp, gpf);
      auto fmstring = fmt::format("{:<10.8e} {:<10.8e} {:<10.8e} {:<10.8e} {:<10.8e}\n",
                                  pos[0], pos[1], pos[2],
                                  std::get<double>(gpmat.gpr_energy_forces()),
                                  mtmp.getPotentialEnergy());
      oname.print(fmstring);
    }
  }
}

std::vector<Matter> helper_functions::getMatterVector(const std::vector<GPRMatter>& gpmatvec){
  std::vector<Matter> matvec;
  std::transform(gpmatvec.begin(),
                 gpmatvec.end(),
                 std::back_inserter(matvec),
                 [&](GPRMatter gpmat)->Matter{
                   return gpmat.truePotMatter;
                 });
  return matvec;
}
