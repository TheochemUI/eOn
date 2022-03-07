#ifndef STATICDATA_H
#define STATICDATA_H

#pragma once

#include <string>
#include <vector>
using namespace std::string_literals; // For ""s

namespace elements {
// These are all null terminated strings so each element can be used as char* with blah.data()
    static const std::vector<std::string> symbols {"Unknown"s, "H"s,  "He"s, "Li"s, "Be"s, "B"s,  "C"s,  "N"s,  "O"s,  "F"s,  "Ne"s, "Na"s,
       "Mg"s,      "Al"s, "Si"s, "P"s,  "S"s,  "Cl"s, "Ar"s, "K"s,  "Ca"s, "Sc"s, "Ti"s, "V"s,
       "Cr"s,      "Mn"s, "Fe"s, "Co"s, "Ni"s, "Cu"s, "Zn"s, "Ga"s, "Ge"s, "As"s, "Se"s, "Br"s,
       "Kr"s,      "Rb"s, "Sr"s, "Y"s,  "Zr"s, "Nb"s, "Mo"s, "Tc"s, "Ru"s, "Rh"s, "Pd"s, "Ag"s,
       "Cd"s,      "In"s, "Sn"s, "Sb"s, "Te"s, "I"s,  "Xe"s, "Cs"s, "Ba"s, "La"s, "Ce"s, "Pr"s,
       "Nd"s,      "Pm"s, "Sm"s, "Eu"s, "Gd"s, "Tb"s, "Dy"s, "Ho"s, "Er"s, "Tm"s, "Yb"s, "Lu"s,
       "Hf"s,      "Ta"s, "W"s,  "Re"s, "Os"s, "Ir"s, "Pt"s, "Au"s, "Hg"s, "Tl"s, "Pb"s, "Bi"s,
       "Po"s,      "At"s, "Rn"s, "Fr"s, "Ra"s, "Ac"s, "Th"s, "Pa"s, "U"s};
}

namespace StatCode{
}
// TODO: Make a strings:: namespace

namespace JobStrings {
	//! Runs a \ref ProcessSearchJob "Process Search" job 
    static const std::string PROCESS_SEARCH {"process_search"s};
	//! Runs a \ref SaddleSearchJob "Saddle Search" job
    static const std::string SADDLE_SEARCH             { "saddle_search"s };
	//! ref MinimizationJob
    static const std::string MINIMIZATION              { "minimization"s };
        //! ref PointJob
    static const std::string POINT                     { "point"s };
        //! ref ParallelReplicaJob
    static const std::string PARALLEL_REPLICA          { "parallel_replica"s };
	//! ref ReplicaExchangeJob
    static const std::string REPLICA_EXCHANGE          { "replica_exchange"s };
	//! ref BasinHopingJob
    static const std::string BASIN_HOPPING             { "basin_hopping"s };
	//! ref HessianJob
    static const std::string HESSIAN                   { "hessian"s };
	//! ref FiniteDifferenceJob
    static const std::string FINITE_DIFFERENCE         { "finite_difference"s };
	//! ref NudgedElasticBandJob
    static const std::string NUDGED_ELASTIC_BAND       { "nudged_elastic_band"s };
	//! ref DynamicsJob
    static const std::string DYNAMICS                  { "molecular_dynamics"s };
        //! ref SafeHyperdynamicsJob
    static const std::string SAFE_HYPER                { "safe_hyper"s };
	//! ref TADJob
    static const std::string TAD                       { "tad"s };
	//! ref PrefactorJob
    static const std::string PREFACTOR                 { "prefactor"s };
	//! ref GlobalOptimizationJob
    static const std::string GLOBAL_OPTIMIZATION       { "global_optimization"s };
	//! ref StructureComparisonJob
    static const std::string STRUCTURE_COMPARISON      { "structure_comparison"s };
	//! ref MonteCarloJob
    static const std::string MONTE_CARLO               { "monte_carlo"s };
	// //! ref TestJob
    //     static const std::string TEST;
}

namespace PrefactorStrings {
    static const std::string RATE_HTST             { "htst"s };
    static const std::string RATE_QQHTST           { "qqhtst"s };
    static const std::string FILTER_CUTOFF         { "cutoff"s };
    static const std::string FILTER_FRACTION       { "fraction"s };

}

namespace PrefactorJobStrings {
    static const std::string PREFACTOR_REACTANT { "reactant"s };
    static const std::string PREFACTOR_SADDLE { "saddle"s };
    static const std::string PREFACTOR_PRODUCT { "product"s };

}

namespace PotentialStrings {

static const std::string POT_EMT          { "emt"s };
static const std::string POT_EXT         { "ext_pot"s };
static const std::string POT_LJ          { "lj"s };
static const std::string POT_LJCLUSTER   { "lj_cluster"s };
static const std::string POT_MORSE_PT    { "morse_pt"s };
static const std::string POT_NEW         { "new_pot"s };

#ifdef IMD_POT
static const std::string POT_IMD         { "imd"s };
#endif

#ifdef WITH_GPRD
static const std::string POT_GPR         { "gpr_pot"s };
#endif

#ifdef WITH_WATER
static const std::string POT_TIP4P       { "tip4p"s };
static const std::string POT_TIP4P_PT    { "tip4p_pt"s };
#ifdef WITH_FORTRAN
static const std::string POT_TIP4P_H     { "tip4p_h"s };
#endif
static const std::string POT_SPCE        { "spce"s };
#endif

#ifdef WITH_FORTRAN
static const std::string POT_EAM_AL      { "eam_al"s };
static const std::string POT_EDIP        { "edip"s };
static const std::string POT_FEHE        { "fehe"s };
static const std::string POT_LENOSKY_SI  { "lenosky_si"s };
static const std::string POT_SW_SI       { "sw_si"s };
static const std::string POT_TERSOFF_SI  { "tersoff_si"s };
#endif

#ifdef LAMMPS_POT
static const std::string POT_LAMMPS      { "lammps"s };
#endif

#ifdef EONMPI
static const std::string POT_MPI         { "mpi"s };
#endif

#ifdef WITH_PYTHON
 #ifdef PYAMFF_POT
static const std::string POT_PYAMFF      { "pyamff"s };
 #endif
static const std::string POT_QSC         { "qsc"s };
#endif

#ifdef WITH_AMS
static const std::string POT_AMS         { "ams"s };
static const std::string POT_AMS_IO      { "ams_io"s };
#endif

#ifdef WITH_VASP
static const std::string POT_VASP        { "vasp"s };
#endif

}

namespace EpiCentersStrings {
    static const std::string DISP_LOAD               { "load"s };
    static const std::string DISP_NOT_FCC_OR_HCP     { "not_fcc_hcp_coordinated"s };
    static const std::string DISP_MIN_COORDINATED    { "least_coordinated"s };
    static const std::string DISP_LAST_ATOM          { "last_atom"s };
    static const std::string DISP_RANDOM             { "random"s };

}

namespace LowestEigenmodeStrings {
static const std::string MINMODE_DIMER {  "dimer"s };
static const std::string MINMODE_GPRDIMER {  "gprdimer"s };
static const std::string MINMODE_LANCZOS {"lanczos"s };
}

namespace ImprovedDimerStrings {
const std::string OPT_SD { "sd"s };
const std::string OPT_CG { "cg"s };
const std::string OPT_LBFGS { "lbfgs"s };
}

namespace DynamicsStrings {
static const std::string ANDERSEN { "andersen"s };
static const std::string NOSE_HOOVER { "nose_hoover"s };
static const std::string LANGEVIN { "langevin"s };
static const std::string NONE { "none"s };
}

namespace Hyperdynamics {
static const std::string NONE { "none"s };
static const std::string BOND_BOOST { "bond_boost"s };
}
#endif /* STATICDATA_H */
