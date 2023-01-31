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
        .def("setPotential", &Matter::setPotential, "Sets the potential function")
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
        // Properties
        // The V variants are just ravel() variants
        .def_property("positions", &Matter::getPositions, &Matter::setPositions)
        .def_property("free_positions", &Matter::getPositionsFree, &Matter::setPositionsFree)
        .def_property("velocities",
                      &Matter::getVelocities,
                      &Matter::setVelocities) // Always sets free velocities
        .def_property(
            "forces",
            &Matter::getForces,
            &Matter::setForces) // Always sets free forces, that is with 0's in the right places
                                // free_forces are basically ndarray[~np.all(m1.forces==0, axis=1)]
        // Optional additional read only properties
        .def_property_readonly("free_forces", &Matter::getForcesFree)
        .def_property_readonly("pot_energy", &Matter::getPotentialEnergy)
        // Getters and Setters
        .def("getAccelerations", &Matter::getAccelerations, "Returns accelerations") // no setter
        .def("getBiasForces", &Matter::getBiasForces, "Returns bias forces applied on all atoms")
        .def("getMass", &Matter::getMass, "Return the mass of the atom specified", py::arg("atom"))
        .def("setMass", &Matter::setMass, "Set the mass of an atom", py::arg("atom"), "mass"_a)
        .def("setMasses",
             &Matter::setMasses,
             "Set the masses of all atoms",
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
        .def("distance",
             py::overload_cast<long int /*index1*/, long int /*index2*/>(&Matter::distance,
                                                                         py::const_),
             "Return the distance between two atoms in same configuration")
        .def("distance",
             py::overload_cast<const Matter &, long int>(&Matter::distance, py::const_),
             "Returns the distance between instances of the same atom in different configurations",
             py::arg("matter"),
             "index"_a)
        .def("pdistance",
             &Matter::pdistance,
             "Returns the distance between two atoms along an axis",
             py::arg("index1"),
             "index2"_a,
             "axis"_a)
        .def_property_readonly("numberOfFreeAtoms", &Matter::numberOfFreeAtoms)
        .def_property_readonly("numberOfFixedAtoms", &Matter::numberOfFixedAtoms)
        .def_property_readonly("force_calls", &Matter::getForceCalls)
        .def("resetForceCalls", &Matter::resetForceCalls, "Zeros the value of force calls")
        .def("maxForce", &Matter::maxForce, "Returns the maximum force")
        // I/O functions
        .def("con2matter",
             py::overload_cast<std::string /*filename*/>(&Matter::con2matter),
             "Read con file into Matter, return true if successful",
             py::arg("filename"))
        .def("con2matter",
             py::overload_cast<FILE* /*file*/>(&Matter::con2matter),
             "Read con file and load data into Matter, return true if successful",
             py::arg("file"))
        .def("convel2matter",
             py::overload_cast<std::string /*filename*/>(&Matter::convel2matter),
             "Read con file with both coordinates and velocities into Matter",
             py::arg("filename"))
        .def("convel2matter",
             py::overload_cast<FILE*  /*file*/>(&Matter::convel2matter),
             "Read con file with both coordinates and velocities into Matter",
             py::arg("file"))
        // BUG: These raise segfaults
        .def("matter2con",
             py::overload_cast<std::string /*filename*/, bool /*append*/>(&Matter::matter2con),
             "Write con file from data in Matter",
             py::arg("filename"),
             "append"_a = false)
        .def("matter2con",
             py::overload_cast<FILE* /*file*/>(&Matter::matter2con),
             "Write con file from data in Matter",
             py::arg("file"))
        .def("matter2convel",
             py::overload_cast<std::string /*filename*/>(&Matter::matter2convel),
             "Write con file with both coordinates and velocities in Matter",
             py::arg("filename"))
        .def("matter2convel",
             py::overload_cast<FILE* /*file*/>(&Matter::matter2convel),
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
        .def("setBiasForces",
             &Matter::setBiasForces,
             "Sets the bias forces",
             py::arg("AtomMatrix_bf"))
        // TODO: BondBoost isn't implemented yet

        /*
        ** Parameters
        */

        /*
        ** Python helpers
        */

        .def("__repr__", [](const Matter &a) { return "<Matter object>"; });
}
