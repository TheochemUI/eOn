import logging
import re
from math import cos, sin

import ase
import numpy
from ase.vibrations import Vibrations
from xtb.ase.calculator import XTB

from eon import atoms
from eon import fileio as io
from eon.config import ConfigClass  # Typing
from eon.config import config as EON_CONFIG

logger = logging.getLogger("displace")

class DisplacementManager:
    def __init__(self, reactant, moved_atoms, config: ConfigClass = EON_CONFIG):
        self.config = config
        self.reactant = reactant

        # Initialize all possible displacement objects if their weight is > 0
        if self.config.displace_random_weight > 0:
            self.random = Random(
                self.reactant,
                self.config.disp_magnitude,
                self.config.disp_radius,
                hole_epicenters=moved_atoms,
                config=self.config,
            )
        if self.config.displace_under_coordinated_weight > 0:
            self.under = Undercoordinated(
                self.reactant,
                self.config.disp_max_coord,
                self.config.disp_magnitude,
                self.config.disp_radius,
                hole_epicenters=moved_atoms,
                cutoff=self.config.comp_neighbor_cutoff,
                use_covalent=self.config.comp_use_covalent,
                covalent_scale=self.config.comp_covalent_scale,
                config=self.config,
            )
        if self.config.displace_least_coordinated_weight > 0:
            self.least = Leastcoordinated(
                self.reactant,
                self.config.disp_magnitude,
                self.config.disp_radius,
                hole_epicenters=moved_atoms,
                cutoff=self.config.comp_neighbor_cutoff,
                use_covalent=self.config.comp_use_covalent,
                covalent_scale=self.config.comp_covalent_scale,
                config=self.config,
            )
        if self.config.displace_listed_atom_weight > 0:
            self.listed_atoms = ListedAtoms(
                self.reactant,
                self.config.disp_magnitude,
                self.config.disp_radius,
                hole_epicenters=moved_atoms,
                cutoff=self.config.comp_neighbor_cutoff,
                use_covalent=self.config.comp_use_covalent,
                covalent_scale=self.config.comp_covalent_scale,
                displace_all=self.config.displace_all_listed,
                config=self.config,
            )
        if self.config.displace_listed_type_weight > 0:
            self.listed_types = ListedTypes(
                self.reactant,
                self.config.disp_magnitude,
                self.config.disp_radius,
                hole_epicenters=moved_atoms,
                cutoff=self.config.comp_neighbor_cutoff,
                use_covalent=self.config.comp_use_covalent,
                covalent_scale=self.config.comp_covalent_scale,
                config=self.config,
            )
        if self.config.displace_not_FCC_HCP_weight > 0:
            self.not_FCC_HCP = NotFCCorHCP(
                self.reactant,
                self.config.disp_magnitude,
                self.config.disp_radius,
                hole_epicenters=moved_atoms,
                cutoff=self.config.comp_neighbor_cutoff,
                use_covalent=self.config.comp_use_covalent,
                covalent_scale=self.config.comp_covalent_scale,
                config=self.config,
            )
        if self.config.displace_not_TCP_BCC_weight > 0:
            self.not_TCP_BCC = NotTCPorBCC(
                self.reactant,
                self.config.disp_magnitude,
                self.config.disp_radius,
                hole_epicenters=moved_atoms,
                cutoff=self.config.comp_neighbor_cutoff,
                use_covalent=self.config.comp_use_covalent,
                covalent_scale=self.config.comp_covalent_scale,
                config=self.config,
            )
        if self.config.displace_not_TCP_weight > 0:
            self.not_TCP = NotTCP(
                self.reactant,
                self.config.disp_magnitude,
                self.config.disp_radius,
                hole_epicenters=moved_atoms,
                cutoff=self.config.comp_neighbor_cutoff,
                use_covalent=self.config.comp_use_covalent,
                covalent_scale=self.config.comp_covalent_scale,
                config=self.config,
            )
        if self.config.displace_water_weight > 0:
            self.water = Water(
                self.reactant,
                self.config.stdev_translation,
                self.config.stdev_rotation,
                self.config.molecule_list,
                self.config.disp_at_random,
            )
        if self.config.displace_softest_mode_weight > 0:
            self.softest = SoftestMode(
                self.reactant,
                self.config.disp_magnitude,
                self.config.disp_radius,
                hole_epicenters=moved_atoms,
                config=self.config,
            )

        # Calculate total weight from all possible displacement types
        total = (
            self.config.displace_random_weight
            + self.config.displace_listed_atom_weight
            + self.config.displace_listed_type_weight
            + self.config.displace_under_coordinated_weight
            + self.config.displace_least_coordinated_weight
            + self.config.displace_not_FCC_HCP_weight
            + self.config.displace_not_TCP_BCC_weight
            + self.config.displace_not_TCP_weight
            + self.config.displace_water_weight
            + self.config.displace_softest_mode_weight
        )

        # Build the cumulative probability list
        if total == 0.0:
            # If no fractions are defined, do 100% random displacements.
            self.random = Random(
                self.reactant,
                self.config.disp_magnitude,
                self.config.disp_radius,
                hole_epicenters=moved_atoms,
                config=self.config,
            )
            self.plist = [1.0]
        else:
            self.plist = []
            cumulative_prob = 0.0
            cumulative_prob += self.config.displace_random_weight / total
            self.plist.append(cumulative_prob)
            cumulative_prob += self.config.displace_listed_atom_weight / total
            self.plist.append(cumulative_prob)
            cumulative_prob += self.config.displace_listed_type_weight / total
            self.plist.append(cumulative_prob)
            cumulative_prob += self.config.displace_under_coordinated_weight / total
            self.plist.append(cumulative_prob)
            cumulative_prob += self.config.displace_least_coordinated_weight / total
            self.plist.append(cumulative_prob)
            cumulative_prob += self.config.displace_not_FCC_HCP_weight / total
            self.plist.append(cumulative_prob)
            cumulative_prob += self.config.displace_not_TCP_BCC_weight / total
            self.plist.append(cumulative_prob)
            cumulative_prob += self.config.displace_not_TCP_weight / total
            self.plist.append(cumulative_prob)
            cumulative_prob += self.config.displace_water_weight / total
            self.plist.append(cumulative_prob)
            cumulative_prob += self.config.displace_softest_mode_weight / total
            self.plist.append(cumulative_prob)

    def make_displacement(self):
        disp_types = [
            "random",
            "listed_atoms",
            "listed_types",
            "under",
            "least",
            "not_FCC_HCP",
            "not_TCP_BCC",
            "not_TCP",
            "water",
            "softest_mode",
        ]

        # Add a special case for the 100% random fallback
        if len(self.plist) == 1:
            logger.debug("Made random displacement (fallback)")
            return self.random.make_displacement()

        r = numpy.random.random_sample()
        i = 0
        while i < len(self.plist) - 1 and self.plist[i] < r:
            i += 1
        disp_type = disp_types[i]

        if disp_type == "random":
            logger.debug("Made random displacement")
            return self.random.make_displacement()
        elif disp_type == "listed_atoms":
            logger.debug("Made listed atoms displacement")
            return self.listed_atoms.make_displacement()
        elif disp_type == "listed_types":
            logger.debug("Made listed atom types displacement")
            return self.listed_types.make_displacement()
        elif disp_type == "under":
            logger.debug("Made under-coordinated displacement")
            return self.under.make_displacement()
        elif disp_type == "least":
            logger.debug("Made least-coordinated displacement")
            return self.least.make_displacement()
        elif disp_type == "not_FCC_HCP":
            logger.debug("Made not-FCC-or-HCP displacement")
            return self.not_FCC_HCP.make_displacement()
        elif disp_type == "not_TCP_BCC":
            logger.debug("Made not-TCP-or-BCC displacement")
            return self.not_TCP_BCC.make_displacement()
        elif disp_type == "not_TCP":
            logger.debug("Made not-TCP displacement")
            return self.not_TCP.make_displacement()
        elif disp_type == "water":
            logger.debug("Made water displacement")
            return self.water.make_displacement()
        elif disp_type == "softest_mode":
            logger.debug("Made softest mode displacement")
            return self.softest.make_displacement()

        raise DisplaceError("Failed to select a valid displacement method.")


class NotImplementedError(Exception):
    pass


class DisplaceError(Exception):
    pass


class Displace:
    def __init__(
        self,
        reactant,
        std_dev,
        radius,
        hole_epicenters,
        config: ConfigClass = EON_CONFIG,
    ):
        """Reactant is an Atoms object. std_dev is the standard deviation
        of the normal distribution used to create the random displacements.
        radius is the distance to neighbors that will also be displaced.
        """
        self.config = config
        self.reactant = reactant
        self.std_dev = std_dev
        self.radius = radius
        self.hole_epicenters = hole_epicenters

        # ## mike w.
        self.void_bias_fraction = 0.2  # self.config.void_bias_fraction

        # temporary numpy array of same size as self.reactant.r
        self.temp_array = numpy.zeros(self.reactant.r.shape)

        self.neighbors_list = None

    def make_displacement(self):
        """Writes the displacement_passed.con and mode_passed.dat to path."""
        raise NotImplementedError

    def get_displacement(self, atom_index):
        """Returns a displacement to be added to self.reactant.r"""
        if self.neighbors_list is None:
            self.neighbors_list = atoms.neighbor_list(
                self.reactant, self.radius, self.config.comp_brute_neighbors
            )
        displacement_norm = 0.0
        displacement = numpy.zeros(self.reactant.r.shape)
        if hasattr(atom_index, "__getitem__"):
            logger.debug("Displacement epicenters: ", atom_index)
            neighbors = [
                self.neighbors_list[i]
                for i in range(len(self.neighbors_list))
                if i in atom_index
            ]
            # flatten
            neighbors = sum(neighbors, [])
            neighbors = numpy.array(list(set(neighbors)), dtype=int)

            displaced_atoms = numpy.append(atom_index, neighbors)
        else:
            logger.debug("Displacement epicenter: %d" % atom_index)
            displaced_atoms = [atom_index] + self.neighbors_list[atom_index]

        # ensures that the total displacement vector exceeds a given length, but the current default is zero
        while displacement_norm <= self.config.disp_min_norm:
            displacement = numpy.zeros(self.reactant.r.shape)
            for i in range(len(displaced_atoms)):
                # don't displace frozen atoms
                if not self.reactant.free[displaced_atoms[i]]:
                    continue
                # Displace one of the free atoms by a gaussian distributed
                # random number with a standard deviation of self.std_dev.
                displacement[displaced_atoms[i]] = numpy.random.normal(
                    scale=self.std_dev, size=3
                )
                if self.config.displace_1d:
                    displacement[displaced_atoms[i]] = displacement[
                        displaced_atoms[i]
                    ] * [1, 0, 0]
            displacement_norm = numpy.linalg.norm(displacement)

        ### Mike W.
        ## this is a fraction of the total displacement magnitude so a small that
        ## a hard coded value will not break things later. It's essentially
        ## negligible at this size.
        if self.void_bias_fraction > 1e-6:
            self.neighbor_list_vectors = atoms.neighbor_list_vectors(
                self.reactant, self.radius, self.config.comp_brute_neighbors
            )

            #            print config.random_mode
            #            print sorted(displaced_atoms)
            #            print numpy.linalg.norm(displacement)

            #            for atom_index in displaced_atoms:
            #                print (atom_index)
            #                print (self.neighbors_list[atom_index])
            #                dist_list = []
            #                for vec in self.neighbor_list_vectors[atom_index]:
            #                    dist_list.append(numpy.linalg.norm(vec))
            #                print (dist_list)

            ## treats the nearest neighbors as repulsive, since I keep finding
            ## interstitials
            pseudoelectrostatic_force = numpy.zeros(self.reactant.r.shape)
            for atom_index in displaced_atoms:
                for vec in self.neighbor_list_vectors[atom_index]:
                    mag = (
                        numpy.linalg.norm(vec) + 1e-6
                    )  # I just want to prevent NaNs in perfectly symmetric situations
                    pseudoelectrostatic_force[atom_index] += -vec / (mag**3)
            # now we norm it for mixing
            void_vec = pseudoelectrostatic_force / numpy.linalg.norm(
                pseudoelectrostatic_force
            )

            displacement = (
                self.void_bias_fraction * displacement_norm * void_vec
                + (1 - self.void_bias_fraction) * displacement
            )

        displacement_atoms = self.reactant.copy()
        displacement_atoms.r += displacement

        if self.config.random_mode:
            displacement *= 0.0
            for i in range(len(displaced_atoms)):
                if not self.reactant.free[displaced_atoms[i]]:
                    continue
                displacement[displaced_atoms[i]] = numpy.random.normal(
                    scale=self.std_dev, size=3
                )
                if self.config.displace_1d:
                    displacement[displaced_atoms[i]] = displacement[
                        displaced_atoms[i]
                    ] * [1, 0, 0]

        displacement /= numpy.linalg.norm(displacement)
        return displacement_atoms, displacement / numpy.linalg.norm(displacement)

    def filter_epicenters(self, epicenters):
        """Returns the epicenters that lie in the hole defined by Displace.hole_epicenters.
        If Displace.hole_epicenters is None, all of the epicenters are accepted."""
        if self.hole_epicenters is None:
            return epicenters
        new_epicenters = []
        for e in epicenters:
            if e in self.hole_epicenters:
                new_epicenters.append(e)

        # GH: added to prevent case in which there are no displaceable atoms in the active region
        if len(new_epicenters) == 0:
            logger.warning(
                "No displaceable atoms found in the active region around atoms that moved in the previous transition;"
            )
            logger.warning("  reverting to the full list of any displaceable atoms.")
            return epicenters

        return new_epicenters


class Undercoordinated(Displace):
    def __init__(
        self,
        reactant,
        max_coordination,
        std_dev=0.05,
        radius=5.0,
        hole_epicenters=None,
        cutoff=3.3,
        use_covalent=False,
        covalent_scale=1.3,
        config: ConfigClass = EON_CONFIG,
    ):
        Displace.__init__(self, reactant, std_dev, radius, hole_epicenters, config)

        self.max_coordination = max_coordination
        self.undercoordinated_atoms = []

        self.coordination_distance = cutoff
        self.initialized = False

    def init(self):
        self.initialized = True
        # cns is an array of the coordination numbers for each of the atoms.
        cns = atoms.coordination_numbers(self.reactant, self.coordination_distance)

        # Only allow displacements for atoms <= the maximum coordination and that are free.
        self.undercoordinated_atoms = [
            i
            for i in range(len(cns))
            if cns[i] <= self.max_coordination and self.reactant.free[i] == 1
        ]

        self.undercoordinated_atoms = self.filter_epicenters(
            self.undercoordinated_atoms
        )

        if len(self.undercoordinated_atoms) == 0:
            errmsg = "No free atoms have a coordination of %i or less"
            errmsg = errmsg % self.max_coordination
            raise DisplaceError(errmsg)

    def make_displacement(self):
        """Select an undercoordinated atom and displace all atoms in a radius about it."""
        # TODO: We should make sure that the amount of I/O to disk is what we think it should be:
        #       about 100 kB or so per make_displacement().
        if not self.initialized:
            self.init()
        epicenter = self.undercoordinated_atoms[
            numpy.random.randint(len(self.undercoordinated_atoms))
        ]
        return self.get_displacement(epicenter)


class Leastcoordinated(Displace):
    def __init__(
        self,
        reactant,
        std_dev=0.05,
        radius=5.0,
        hole_epicenters=None,
        cutoff=3.3,
        use_covalent=False,
        covalent_scale=1.3,
        config: ConfigClass = EON_CONFIG,
    ):
        Displace.__init__(self, reactant, std_dev, radius, hole_epicenters, config)

        self.leastcoordinated_atoms = []

        self.coordination_distance = cutoff

        self.leastcoordinated_atoms = atoms.least_coordinated(
            self.reactant, self.coordination_distance
        )
        self.leastcoordinated_atoms = [
            i for i in self.leastcoordinated_atoms if self.reactant.free[i] == 1
        ]

        self.leastcoordinated_atoms = self.filter_epicenters(
            self.leastcoordinated_atoms
        )

        if len(self.leastcoordinated_atoms) == 0:
            errmsg = "The least coordinated atoms are all frozen"
            raise DisplaceError(errmsg)

    def make_displacement(self):
        """Select an undercoordinated atom and displace all atoms in a radius about it."""
        epicenter = self.leastcoordinated_atoms[
            numpy.random.randint(len(self.leastcoordinated_atoms))
        ]
        return self.get_displacement(epicenter)


class ListedAtoms(Displace):
    def __init__(
        self,
        reactant,
        std_dev=0.05,
        radius=5.0,
        hole_epicenters=None,
        cutoff=3.3,
        use_covalent=False,
        covalent_scale=1.3,
        displace_all=False,
        config: ConfigClass = EON_CONFIG,
    ):
        Displace.__init__(self, reactant, std_dev, radius, hole_epicenters, config)

        self.displace_all = displace_all
        # each item in this list is the index of a free atom
        self.listed_atoms = [
            i for i in self.config.disp_listed_atoms if self.reactant.free[i]
        ]
        # print "self.listed_atoms:"
        # print self.listed_atoms
        self.listed_atoms = self.filter_epicenters(self.listed_atoms)
        logger.debug("Listed atoms: %s", self.listed_atoms)

        # print self.listed_atoms
        if len(self.listed_atoms) == 0:
            # raise DisplaceError("Listed atoms are all frozen")
            raise DisplaceError("Listed atoms are all frozen")

    def make_displacement(self):
        """Select a listed atom and displace all atoms in a radius about it."""
        # chose a random atom from the supplied list
        if self.displace_all:
            epicenter = self.listed_atoms
        else:
            epicenter = self.listed_atoms[numpy.random.randint(len(self.listed_atoms))]
        logger.debug("Listed atom displacement epicenters: %s", epicenter)
        return self.get_displacement(epicenter)


class ListedTypes(Displace):
    def __init__(
        self,
        reactant,
        std_dev=0.05,
        radius=5.0,
        hole_epicenters=None,
        cutoff=3.3,
        use_covalent=False,
        covalent_scale=1.3,
        displace_all=False,
        config: ConfigClass = EON_CONFIG,
    ):
        Displace.__init__(self, reactant, std_dev, radius, hole_epicenters, config)

        #        print self.config.disp_listed_types

        self.displace_all = displace_all
        # each item in this list is the index of a free atom
        self.listed_atoms = [
            i
            for i in range(len(self.reactant.free))
            if (self.reactant.free[i] == 1)
            and (self.reactant.names[i] in self.config.disp_listed_types)
        ]

        self.listed_atoms = self.filter_epicenters(self.listed_atoms)

        if len(self.listed_atoms) == 0:
            raise DisplaceError(
                "Listed atom types are all frozen or not found in reactant"
            )

    #        print self.listed_atoms

    def make_displacement(self):
        """Select a listed atom and displace all atoms in a radius about it."""
        # chose a random atom from the supplied list
        if self.displace_all:
            epicenter = self.listed_atoms
        else:
            epicenter = self.listed_atoms[numpy.random.randint(len(self.listed_atoms))]
        return self.get_displacement(epicenter)


class Random(Displace):
    def __init__(
        self,
        reactant,
        std_dev=0.05,
        radius=5.0,
        hole_epicenters=None,
        config: ConfigClass = EON_CONFIG,
    ):
        Displace.__init__(self, reactant, std_dev, radius, hole_epicenters, config)

        # each item in this list is the index of a free atom
        self.free_atoms = [
            i for i in range(len(self.reactant.free)) if self.reactant.free[i]
        ]

        self.free_atoms = self.filter_epicenters(self.free_atoms)

        if len(self.free_atoms) == 0:
            raise DisplaceError("There are no free atoms in the reactant")

    def make_displacement(self):
        """Select a random atom and displace all atoms in a radius about it."""
        # chose a random atom
        epicenter = self.free_atoms[numpy.random.randint(len(self.free_atoms))]
        return self.get_displacement(epicenter)


class NotFCCorHCP(Displace):
    def __init__(
        self,
        reactant,
        std_dev=0.05,
        radius=5.0,
        hole_epicenters=None,
        cutoff=3.3,
        use_covalent=False,
        covalent_scale=1.3,
        config: ConfigClass = EON_CONFIG,
    ):
        Displace.__init__(self, reactant, std_dev, radius, hole_epicenters, config)

        self.not_HCP_or_FCC_atoms = []

        self.coordination_distance = cutoff

        self.not_HCP_or_FCC_atoms = atoms.not_HCP_or_FCC(
            self.reactant, self.coordination_distance
        )

        self.not_HCP_or_FCC_atoms = [
            i for i in self.not_HCP_or_FCC_atoms if self.reactant.free[i] == 1
        ]

        self.not_HCP_or_FCC_atoms = self.filter_epicenters(self.not_HCP_or_FCC_atoms)

        if len(self.not_HCP_or_FCC_atoms) == 0:
            errmsg = "The atoms without FCC or HCP coordination are all frozen"
            raise DisplaceError(errmsg)

    def make_displacement(self):
        """Select an atom without HCP or FCC coordination and displace all atoms in a radius about it."""
        epicenter = self.not_HCP_or_FCC_atoms[
            numpy.random.randint(len(self.not_HCP_or_FCC_atoms))
        ]
        return self.get_displacement(epicenter)


class NotTCPorBCC(Displace):
    def __init__(
        self,
        reactant,
        std_dev=0.05,
        radius=5.0,
        hole_epicenters=None,
        cutoff=3.3,
        use_covalent=False,
        covalent_scale=1.3,
        config: ConfigClass = EON_CONFIG,
    ):
        Displace.__init__(self, reactant, std_dev, radius, hole_epicenters, config)

        self.not_TCP_or_BCC_atoms = []

        self.coordination_distance = cutoff

        self.not_TCP_or_BCC_atoms = atoms.not_TCP_or_BCC(
            self.reactant, self.coordination_distance
        )

        self.not_TCP_or_BCC_atoms = [
            i for i in self.not_TCP_or_BCC_atoms if self.reactant.free[i] == 1
        ]

        self.not_TCP_or_BCC_atoms = self.filter_epicenters(self.not_TCP_or_BCC_atoms)

        if len(self.not_TCP_or_BCC_atoms) == 0:
            errmsg = "The atoms without BCC or TCP coordination are all frozen."
            raise DisplaceError(errmsg)

    def make_displacement(self):
        """Select an atom without TCP or BCC coordination and displace all atoms in a radius about it."""
        epicenter = self.not_TCP_or_BCC_atoms[
            numpy.random.randint(len(self.not_TCP_or_BCC_atoms))
        ]
        return self.get_displacement(epicenter)


class NotTCP(Displace):
    def __init__(
        self,
        reactant,
        std_dev=0.05,
        radius=5.0,
        hole_epicenters=None,
        cutoff=3.3,
        use_covalent=False,
        covalent_scale=1.3,
        config: ConfigClass = EON_CONFIG,
    ):
        Displace.__init__(self, reactant, std_dev, radius, hole_epicenters, config)

        self.not_TCP_atoms = []

        self.coordination_distance = cutoff

        self.not_TCP_atoms = atoms.not_TCP(self.reactant, self.coordination_distance)

        self.not_TCP_atoms = [
            i for i in self.not_TCP_atoms if self.reactant.free[i] == 1
        ]

        self.not_TCP_atoms = self.filter_epicenters(self.not_TCP_atoms)

        if len(self.not_TCP_atoms) == 0:
            errmsg = "The atoms without TCP coordination are all frozen"
            raise DisplaceError(errmsg)

    def make_displacement(self):
        """Select an atom without HCP or FCC coordination and displace all atoms in a radius about it."""
        epicenter = self.not_TCP_atoms[numpy.random.randint(len(self.not_TCP_atoms))]
        return self.get_displacement(epicenter)


# XXX(rg): Why doesn't this actually form a child class of Displace? No initialization..
class Water(Displace):
    """Displace molecules of water without streatching them."""

    def __init__(
        self, reactant, stdev_translation, stdev_rotation, molecule_list=[], random=0
    ):
        (
            """reactant: structure to be displaced\n"""
            """stdev_translation: translational standard deviation (Angstrom)\n"""
            """stdev_rotation: rotational standard deviation (radian)"""
            """molecules: list of indices of the molecules to displace or None to displace all the molecules"""
            """random: if 0 displace all molecules in molecule_list, if 'random > 0' picked up at random in 'molecule_list' a number of molecules equal to the number soted in 'random' and displace only these"""
        )
        assert isinstance(molecule_list, list)
        self.reactant = reactant
        self.stdev_translation = stdev_translation
        self.stdev_rotation = stdev_rotation
        for i, name in enumerate(reactant.names):
            if not re.search("^H", name):
                break
        # For water assume that all the hydrogen are listed first, then all the oxygen
        self.n_water = i / 2
        if len(molecule_list) == 0:
            molecule_list = list(range(self.n_water))
        self.molecule_list = molecule_list
        self.random = random

    def make_displacement(self):
        """Returns Atom object containing displaced structure and an array containing the displacement."""
        return self.get_displacement()

    def get_displacement(self):
        """Returns Atom object containing displaced structure and an array containing the displacement."""
        free = self.reactant.free
        displaced_atoms = self.reactant.copy()
        if self.random > 0:
            n = len(self.molecule_list)
            molecule_list = list()
            for i in range(self.random):
                i = int(numpy.random.uniform(0, n))
                molecule_list.append(self.molecule_list[i])
        else:
            molecule_list = self.molecule_list
        for i in molecule_list:
            # print 'displacing', i
            h1 = i * 2
            h2 = i * 2 + 1
            o = i + self.n_water * 2
            # don't displace if any of the three atoms is fixed
            # if not (free[h1] and free[h2] and free[o]):
            #    continue
            # Displace one of the free atoms by a gaussian distributed
            # random number with a standard deviation of self.std_dev.
            disp = numpy.random.normal(scale=self.stdev_translation, size=3)
            displaced_atoms.r[h1] += disp
            displaced_atoms.r[h2] += disp
            displaced_atoms.r[o] += disp
            rh1 = displaced_atoms.r[h1]
            rh2 = displaced_atoms.r[h2]
            ro = displaced_atoms.r[o]
            disp = numpy.random.normal(scale=self.stdev_rotation, size=3)
            rh1, rh2, ro = self.rotate_water(rh1, rh2, ro, disp[0], disp[1], disp[2])
            displaced_atoms.r[h1] = rh1
            displaced_atoms.r[h2] = rh2
            displaced_atoms.r[o] = ro
        displacement = displaced_atoms.r - self.reactant.r
        return displaced_atoms, displacement

    ## Rotate a molecule of water.
    # Rotate around the centre of gravity of the molecule.
    # @param[in] "hydrogen1, hydrogen2, oxygen" numpy.array: coordinates of atoms.
    # @param[in] "psi, theta, phi" float: Angle in @em radians of rotation around <em> x, y, z </em> axes.
    # @param[in] "hydrogenMass, oxygenMass" float: masses of the hydrogen and oxygen atoms.
    # @return "hydrogen1, hydrogen2, oxygen" coordinates of atoms after rotation
    # The rotations uses the @em x, y, z convention (pitch-roll-yaw). The fisrt rotation to take place is around z, then y, then x.
    # Equations and notations are from: http://mathworld.wolfram.com/EulerAngles.html .
    @staticmethod
    def rotate_water(
        hydrogen1,
        hydrogen2,
        oxygen,
        psi,
        theta,
        phi,
        hydrogen_mass=1.0,
        oxygen_mass=16.0,
    ):
        G = (hydrogen_mass * (hydrogen1 + hydrogen2) + oxygen_mass * oxygen) / (
            hydrogen_mass * 2.0 + oxygen_mass
        )
        rot = numpy.array(
            [
                [cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta)],
                [
                    sin(psi) * sin(theta) * cos(phi) - cos(psi) * sin(phi),
                    sin(psi) * sin(theta) * sin(phi) + cos(psi) * cos(phi),
                    cos(theta) * sin(psi),
                ],
                [
                    cos(psi) * sin(theta) * cos(phi) + sin(psi) * sin(phi),
                    cos(psi) * sin(theta) * sin(phi) - sin(psi) * cos(phi),
                    cos(theta) * cos(psi),
                ],
            ]
        )
        rh1 = numpy.tensordot(rot, (hydrogen1 - G), 1) + G
        rh2 = numpy.tensordot(rot, (hydrogen2 - G), 1) + G
        ro = numpy.tensordot(rot, (oxygen - G), 1) + G
        return rh1, rh2, ro


if __name__ == "__main__":
    import sys
    import time

    if len(sys.argv) < 3:
        print("%s: reactant.con outpath" % sys.argv[0])
        sys.exit(1)

    reactant = io.loadcon(sys.argv[1])
    d = Random(reactant, 0.05, 5.0)
    # d = Undercoordinated(reactant, 11, 0.05, 5.0)
    t0 = time.time()
    ntimes = 1000
    for i in range(ntimes):
        d.make_displacement()
    dt = time.time() - t0
    outf = io.savecon(sys.argv[2], reactant)
    print("%.2f displacements per second" % (float(ntimes) / dt))


class SoftestMode(Displace):
    """
    Displaces atoms along the softest vibrational mode, typically corresponding
    to the imaginary frequency of a transition state.
    """

    def __init__(
        self,
        reactant,
        magnitude,
        radius,
        hole_epicenters=None,
        config: ConfigClass = EON_CONFIG,
    ):
        super().__init__(reactant, magnitude, radius, hole_epicenters, config)
        self.softest_mode_vector = None
        self._calculate_softest_mode()

    def _calculate_softest_mode(self):
        """
        Performs a robust vibrational analysis using XTB to find the softest mode,
        correctly handling electronic structure and boundary conditions.
        """
        logger.info(
            "Calculating vibrational modes with XTB to find the softest mode..."
        )

        atoms_for_vib = ase.Atoms(
            symbols=self.reactant.names,
            positions=self.reactant.r,
            cell=self.reactant.box,
            pbc=False,
        )

        charge = getattr(self.config, "charge", 0)
        multiplicity = getattr(self.config, "multiplicity", None)

        if multiplicity is None:
            total_electrons = atoms_for_vib.get_atomic_numbers().sum() - charge
            if total_electrons % 2 == 1:
                multiplicity = 2
                logger.info(
                    "Guessed multiplicity=2 (doublet) from odd number of electrons."
                )
            else:
                multiplicity = 1
                logger.info(
                    "Guessed multiplicity=1 (singlet) from even number of electrons."
                )

        uhf_electrons = multiplicity - 1
        atoms_for_vib.calc = XTB(method="GFN2-xTB", charge=charge, uhf=uhf_electrons)

        try:
            vib = Vibrations(atoms_for_vib, name="vib_temp_softest_mode")
            vib.run()

            all_frequencies = vib.get_frequencies()
            freq = all_frequencies[0]
            mode_vector = vib.get_mode(0)

            logger.info(f"Found softest mode with frequency: {freq.real:.4f} cm^-1")

            if freq.real > 0:
                logger.warning(
                    "The softest mode has a real frequency. "
                    "This structure may be a minimum, not a saddle point."
                )

            self.softest_mode_vector = mode_vector / numpy.linalg.norm(mode_vector)
            vib.clean()

        except Exception as e:
            logger.error(f"XTB vibrational analysis failed: {e}")
            raise DisplaceError(
                "Could not calculate the softest mode for displacement."
            )

    def make_displacement(self):
        """
        Generates the displaced Atoms object and the normalized displacement vector.
        """
        if self.softest_mode_vector is None:
            raise DisplaceError("Softest mode vector has not been calculated.")

        displacement = self.std_dev * self.softest_mode_vector

        for i, is_free in enumerate(self.reactant.free):
            if not is_free:
                displacement[i] = [0.0, 0.0, 0.0]

        current_norm = numpy.linalg.norm(displacement)
        if current_norm > 1e-9:
            displacement *= self.std_dev / current_norm

        displaced_atoms = self.reactant.copy()
        displaced_atoms.r += displacement

        norm = numpy.linalg.norm(displacement)
        if norm < 1e-9:
            return displaced_atoms, displacement
        final_normalized_displacement = displacement / norm
        return displaced_atoms, final_normalized_displacement
