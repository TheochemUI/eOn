#!/usr/bin/env python
import ase
import numpy as np
import tsase
from random import random
from sys import argv


atoms = ase.Atoms(argv[1])
atoms.set_cell((30.0,30.0,30.0))
atoms.set_positions(np.random.random((len(atoms),3))*10.0)
atoms.center(30.0)

tsase.io.write_con('matter1.con', atoms)

atoms.positions += random()*3.0
atoms.rattle(stdev=.002, seed=int(random()*2**30))
atoms.rotate_euler(center='COM', phi=random()*np.pi*2, theta=random()*np.pi, psi=random()*np.pi/2)

tsase.io.write_con('matter2.con', atoms)
