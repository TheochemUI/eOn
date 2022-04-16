#include "GPRMatter.h"

GPRMatter::~GPRMatter(){

};

GPRMatter::GPRMatter(Matter initMatter) : truePotMatter{initMatter} {
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
    gprfunc.initialize(potparams, this->atmconf);
}
