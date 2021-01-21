// An interface to the GPDimer library

#include "HelperFunctions.h"
#include "AtomicGPDimer.h"
#include "Log.h"
#include <cmath>
#include <cassert>

using namespace helper_functions;

const char AtomicGPDimer::OPT_SCG[] = "scg";
const char AtomicGPDimer::OPT_LBFGS[] = "lbfgs";

AtomicGPDimer::AtomicGPDimer(Matter *matter, Parameters *params)
{
    parameters    = params;
    x0            = new Matter(parameters);
    x1            = new Matter(parameters);
    *x0           = *matter;
    *x1           = *matter;
    tau.resize(3*matter->numberOfAtoms());
    tau.setZero();
    totalForceCalls = 0;

    if(parameters->dimerOptMethod == OPT_SCG){
        init_cg = true;
    }

    atomic_dimer = new atmd::AtomicDimer;

    typedef gpr::Field<double> FieldDbl;

    std::map<std::string, gpr::Coord*> dict_coord;
    std::map<std::string, FieldDbl*> dict_field;

    io::FileManager fm;

    aux::ProblemSetUp problem_setup;
    gpr::Coord R_all_init;  // no initial data points
    FieldDbl E_all_init;    // no initial data
    gpr::Coord G_all_init;  // no initial data
    gpr::Coord R_init;
    FieldDbl E_init;
    gpr::Coord G_init;
    gpr::Coord orient_init, orient_init_tmp, orient_start;
    gpr::Coord R_sp, R;
    ConfInfo conf_info;
    ConfInfo conf_info_inactive;
    FieldDbl conf;
    Params parameters;
    Index_t D = 0;
    Index_t N_mov = 0;
    Observation init_observations;
    Observation init_middle_point;
    Index_t number_of_H2 = 2;
    Index_t number_of_Cu = 216;
    gpr::Field<Index_t> conf_atom;
    AtomsConfiguration atoms_config;

    // Read input parameters
    fm.readInputFile("input/input.dat", parameters);

    // Read the initial orientation from the input file
    dict_coord["orient_init"] = &orient_init_tmp;
    dict_coord["orient_start"] = &orient_start;
    fm.readDataFile("input/orients_CuH2.dat", dict_coord);
    dict_coord.clear();

    // Read the initial middle point from the input file
    dict_coord["R_sp"] = &R_sp;
    dict_coord["R"] = &R;
    fm.readDataFile("input/MEP_CuH2.dat", dict_coord);
    dict_coord.clear();

    // Read the initial configuration data from the input file
    dict_field["conf"] = &conf;
    fm.readDataFile("input/CuH2-idpp8.dat", dict_field);
    dict_field.clear();

    // Set up a structure with configuration of all atoms
    atoms_config.assignFromField(conf);

    // Define the initial middle point of the dimer
    // Initial middle point not observed
    R_init.resize(1, R_sp.getNj());
    double dist_sp = parameters.dist_sp.value[parameters.i_dist.value];
    for(Index_t n = 0; n < R_sp.getNj(); ++n) {
        R_init[n] = R_sp[n] + dist_sp * orient_start(parameters.i_run.value, n);
    }
    E_init.clear();
    G_init.clear();

    // Define the initial orientation of the dimer (unit vector along the
    // direction of the dimer)
    orient_init.append(orient_init_tmp, parameters.i_run.value);

    // Misc
    D = R_init.getNj(); // number of dimensions
    N_mov = D / 3; // number of moving atoms

    // 'conf_info' is a structure including information about the configurations
    // necessary for the GP model
    // 'conf_info.conf_fro': coordinates of active frozen atoms (N_fro x 3)
    // In the beginning, list none of the frozen atoms as active
    conf_info.conf_fro.clear();
    // 'conf_info.atomtype_mov': atomtype indices for moving atoms (1 x N_mov)
    conf_info.atomtype_mov.resize(1, N_mov);
    conf_info.atomtype_mov.set(0); // H atoms
    // 'conf_info.atomtype_fro': atomtype indices for active frozen atoms (1 x N_fro)
    conf_info.atomtype_fro.clear();
    // The atomtypes must be indexed as 1,2,...,n_at.
    Index_t n_at = 2; // number of atomtypes (including also the types of inactive frozen atoms)
    // 'conf_info.pairtype': pairtype indices for pairs of atomtypes (n_at x n_at)
    // Active pairtypes are indexed as 1,2,...,n_pt. Inactive pairtypes are given index 0.
    conf_info.pairtype.resize(n_at, n_at);
    conf_info.pairtype.set(EMPTY);
    // 'conf_info.n_pt': number of active pairtypes
    conf_info.n_pt = 0;
    // Set pairtype indices for moving+moving atom pairs (and update number of active pairtypes):
    problem_setup.set_pairtype_mov(conf_info.atomtype_mov, conf_info.n_pt, conf_info.pairtype);

    conf_atom.resize(1, number_of_Cu + number_of_H2);
    for(Index_t i = 0; i < conf.getNi(); ++i) {
        conf_atom[i] = conf(i, 4);
    }

    // 'conf_info_inactive' is a structure including information about inactive frozen atoms:
    // 'conf_info_inactive.conf_ifro': coordinates of inactive frozen atoms (N_ifro x 3)
    // In the beginning, list all frozen atoms as inactive:
    conf_info_inactive.conf_fro.resize(1, 3 * number_of_Cu);
    for(Index_t i = 0; i < number_of_Cu; ++i) {
        conf_info_inactive.conf_fro.set(0, i, {conf(i, 0), conf(i, 1), conf(i, 2)});
    }

    // 'conf_info_inactive.atomtype_ifro': atomtype indices for inactive frozen atoms (1 x N_ifro)
    conf_info_inactive.atomtype_fro.resize(1, conf_info_inactive.conf_fro.getNumPoints());
    conf_info_inactive.atomtype_fro.set(1); // Cu atoms

    // Activate frozen atoms within activation distance:
    problem_setup.activateFrozenAtoms(R_init, parameters.actdist_fro.value,
                                      conf_info, conf_info_inactive);

    init_observations.clear();
    init_middle_point.clear();
    init_middle_point.R = R_init;

    atomic_dimer->initialize(parameters, init_observations, init_middle_point,
                            orient_init, conf_info, conf_info_inactive, atoms_config);
}

AtomicGPDimer::~AtomicGPDimer()
{
    delete x0;
    delete x1;
    delete atomic_dimer;
}

void AtomicGPDimer::compute(Matter *matter, AtomMatrix initialDirectionAtomMatrix)
{
    atomic_dimer->execute();
}

double AtomicGPDimer::getEigenvalue()
{
    atomic_dimer->getFinalCurvature();
}

AtomMatrix AtomicGPDimer::getEigenvector()
{
    // return atomic_dimer->getFinalOrientation();
    return MatrixXd::Map(tau.data(), x0->numberOfAtoms(), 3);
}
