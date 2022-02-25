// clang-format off
#include "py_wrapper.hpp"
// clang-format on

void py_matter(py::module_ &m) {

    py::class_<MatterPrivateData>(
        m, "MatterPrivateData", py::dynamic_attr()) // dynamic incurs a penalty
        .def(py::init())
        /*
        ** Functions
        */

        /*
        ** Parameters
        */
        .def_readwrite("parameters", &MatterPrivateData::parameters, "Parameters")
        .def_readwrite("nAtoms", &MatterPrivateData::nAtoms, "Number of atoms")
        .def_readwrite("positions", &MatterPrivateData::positions, "Positions")
        .def_readwrite("velocities", &MatterPrivateData::velocities, "Velocities")
        .def_readwrite("forces", &MatterPrivateData::forces, "Forces")
        // TODO: BondBoost isn't implemented yet
        // .def_readwrite("biasPotential", &MatterPrivateData::biasPotential, "Bias potential")
        .def_readwrite("masses", &MatterPrivateData::masses, "Masses")
        .def_readwrite("atomicNrs", &MatterPrivateData::atomicNrs, "Atomic numbers")
        .def_readwrite("isFixed",
                       &MatterPrivateData::isFixed,
                       "An array of bool, false for movable atom, true for fixed")
        .def_readwrite("cell", &MatterPrivateData::cell, "Cell")
        .def_readwrite("cellInverse", &MatterPrivateData::cellInverse, "Cell Inverse")
        .def_readwrite("potentialEnergy", &MatterPrivateData::potentialEnergy, "Potential Energy")

        /*
        ** Python helpers
        */

        .def("__repr__", [](const MatterPrivateData &a) { return "<MatterPrivateData object>"; });

    py::class_<Matter>(m, "Matter", py::dynamic_attr()) // dynamic incurs a penalty
        /*
        ** Constructors
        */
        .def(py::init<Parameters *>())
        .def(py::init<Parameters *, long double>())
        .def(py::init<const Matter &>())
        /*
        ** Operators
        */
        // operator= copy constructor not bound yet
        // .def(py::self = py::self)
        /*
        ** Methods
        */
        .def("compare",
             &Matter::compare,
             "Comparison",
             py::arg("matter"),
             py::arg("indistinguishable"))
        .def("distanceTo",
             &Matter::distanceTo,
             "The distance to the given matter object",
             py::arg("matter"))
        .def("perAtomNorm",
             &Matter::perAtomNorm,
             "The maximum distance between two atoms in the Matter objects",
             py::arg("matter"))
        // TODO: See why this is unresolved
        // .def("setPotential", &Matter::setPotential, "Sets the potential function")
        .def("resize", &Matter::resize, "Set or reset the number of atoms", py::arg("nAtoms"))
        .def("numberOfAtoms", &Matter::numberOfAtoms, "Return the number of atoms")
        .def("getCell", &Matter::getCell, "Returns the cell dimensions")
        .def("setCell", &Matter::setCell, "Sets the cell dimensions", py::arg("newCell"))
        .def("getPosition",
             &Matter::getPosition,
             "Return the position of an atom along one of the axis",
             py::arg("atom"),
             "axis"_a)
        .def("setPosition",
             &Matter::setPosition,
             "Set the position of atom along axis to position",
             py::arg("atom"),
             "axis"_a,
             "position"_a) //
        .def("setVelocity",
             &Matter::setVelocity,
             "Set the velocity of atom along axis to velocity",
             py::arg("atom"),
             "axis"_a,
             "velocity"_a)
        .def("relax",
             &Matter::relax,
             "Relax structure",
             py::arg("quiet"),
             "writeMovie"_a,
             "checkpoint"_a,
             "prefixMovie"_a,
             "prefixCheckpoint"_a)
        .def("pbc", &Matter::pbc, "Periodic boundary conditions", py::arg("diff"))
        .def("pbcV", &Matter::pbcV, "PBC as vector", py::arg("diff"))
        // Getters and Setters
        .def("getPositions", &Matter::getPositions, "Return coordinates of atoms in array")
        .def("getPositionsV",
             &Matter::getPositionsV,
             "Return coordinates of atoms in array (as vector)")
        .def("getPositionsFree",
             &Matter::getPositionsFree,
             "Return coordinates of free atoms in array")
        .def("getPositionsFreeV",
             &Matter::getPositionsFreeV,
             "Return coordinates of free atoms in array (as vector)")
        .def("setPositions",
             &Matter::setPositions,
             "Update Matter with the new positions of the atoms in the array",
             py::arg("AtomMatrix_pos"))
        .def("setPositionsV",
             &Matter::setPositionsV,
             "Update Matter with the new positions of the atoms in the array (as vector)",
             py::arg("VectorXd_pos"))
        .def("setPositionsFree",
             &Matter::setPositionsFree,
             "Update Matter with the new positions of the free atoms given in array",
             py::arg("AtomMatrix_pos"))
        .def("setPositionsFreeV",
             &Matter::setPositionsFreeV,
             "Update Matter with the new positions of the free atoms given in array (as vector)",
             py::arg("VectorXd_pos"))
        .def("getVelocities", &Matter::getVelocities, "Returns velocities")
        .def("setVelocities", &Matter::setVelocities, "Sets velocities", py::arg("AtomMatrix_v"))
        .def("setBiasForces",
             &Matter::setBiasForces,
             "Sets the bias forces",
             py::arg("AtomMatrix_bf"))
        // TODO: BondBoost isn't implemented yet
        // .def("setBiasPotential",
        //      &Matter::setBiasPotential,
        //      "Sets a BondBoost potential",
        //      py::arg("BondBoost"))
        .def("setForces", &Matter::setForces, "Sets the forces", py::arg("AtomMatrix_f"))
        .def("getAccelerations", &Matter::getAccelerations, "Returns accelerations")
        .def("getForces", &Matter::getForces, "Returns forces applied on all atoms")
        .def("getBiasForces", &Matter::getBiasForces, "Returns bias forces applied on all atoms")
        .def("getForcesV", &Matter::getForcesV, "Returns forces applied on all atoms (as vector)")
        .def("getForcesFree", &Matter::getForcesFree, "Returns free forces applied on all atoms")
        .def("getForcesFreeV",
             &Matter::getForcesFreeV,
             "Returns free forces applied on all atoms (as vector)")
        .def("getMass", &Matter::getMass, "Return the mass of the atom specified", py::arg("atom"))
        .def("setMass", &Matter::setMass, "Set the mass of an atom", py::arg("atom"), "mass"_a)
        .def("setMasses",
             &Matter::setMasses,
             "Set the mass of an atom",
             py::arg("VectorXd_massesIn"))
        .def("getAtomicNr",
             &Matter::getAtomicNr,
             "Return the atomic number of the atom specified",
             py::arg("atom"))
        .def("setAtomicNr",
             &Matter::setAtomicNr,
             "Set the atomic number of an atom",
             py::arg("atom"),
             "atomicNr"_a)
        .def("getFixed",
             &Matter::getFixed,
             "Return true if the atom is fixed, false if it is movable",
             py::arg("atom"))
        .def("setFixed",
             &Matter::setFixed,
             "Set the atom to fixed (true) or movable (false)",
             py::arg("atom"),
             "isFixed"_a)
        .def("getPotentialEnergy", &Matter::getPotentialEnergy, "Return the potential energy")
        .def("getKineticEnergy", &Matter::getKineticEnergy, "Return the kinetic energy")
        .def("getMechanicalEnergy",
             &Matter::getMechanicalEnergy,
             "Return the mechanical energy (i.e. kinetic plus potential energy)")
        // TODO: Look into this
        // .def("distance",
        //      py::overload_cast<long int /*index1*/, long int /*index2*/>(&Matter::distance),
        //      "Return the distance between two atoms in same configuration",
        //      py::arg("index1"),
        //      "index2"_a)
        // .def("distance",
        //      py::overload_cast<const Matter &, long int>(&Matter::distance),
        //      "Returns the distance between instances of the same atom in different
        //      configurations", py::arg("matter"), "index"_a)
        .def("pdistance",
             &Matter::pdistance,
             "Returns the distance between two atoms along an axis",
             py::arg("index1"),
             "index2"_a,
             "axis"_a)
        .def("numberOfFreeAtoms",
             &Matter::numberOfFreeAtoms,
             "Return the number of free (or movable) atoms")
        .def("numberOfFixedAtoms", &Matter::numberOfFixedAtoms, "Return the number of fixed atoms")
        .def("getForceCalls", &Matter::getForceCalls, "Return the number of force calls so far")
        .def("resetForceCalls", &Matter::resetForceCalls, "Zeros the value of force calls")
        .def("maxForce", &Matter::maxForce, "Returns the maximum force")
        // I/O functions
        .def("con2matter",
             py::overload_cast<std::string /*filename*/>(&Matter::con2matter),
             "Read con file into Matter, return true if successful",
             py::arg("filename"))
        .def("con2matter",
             py::overload_cast<std::ifstream & /*file*/>(&Matter::con2matter),
             "Read con file and load data into Matter, return true if successful",
             py::arg("file"))
        .def("convel2matter",
             py::overload_cast<std::string /*filename*/>(&Matter::convel2matter),
             "Read con file with both coordinates and velocities into Matter",
             py::arg("filename"))
        .def("convel2matter",
             py::overload_cast<std::ifstream & /*file*/>(&Matter::convel2matter),
             "Read con file with both coordinates and velocities into Matter",
             py::arg("file"))
        // BUG: These raise segfaults
        .def("matter2con",
             py::overload_cast<std::string /*filename*/, bool /*append*/>(&Matter::matter2con),
             "Write con file from data in Matter",
             py::arg("filename"),
             "append"_a = false)
        .def("matter2con",
             py::overload_cast<std::ofstream & /*file*/>(&Matter::matter2con),
             "Write con file from data in Matter",
             py::arg("file"))
        .def("matter2convel",
             py::overload_cast<std::string /*filename*/>(&Matter::matter2convel),
             "Write con file with both coordinates and velocities in Matter",
             py::arg("filename"))
        .def("matter2convel",
             py::overload_cast<std::ofstream & /*file*/>(&Matter::matter2convel),
             "Write con file with both coordinates and velocities in Matter",
             py::arg("file"))
        .def("matter2xyz",
             py::overload_cast<std::string /*filename*/, bool /*append*/>(&Matter::matter2xyz),
             "Write xyz file from data in Matter",
             py::arg("filename"),
             "append"_a = false)
        // More getters and setters
        .def("getFree", &Matter::getFree, "Get free")
        .def("getFreeV", &Matter::getFreeV, "Get free (as vector)")
        .def("getMasses", &Matter::getMasses, "Get masses")

        /*
        ** Parameters
        */

        /*
        ** Python helpers
        */

        .def("__repr__", [](const Matter &a) { return "<Matter object>"; });
}
