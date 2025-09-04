#!/usr/bin/env python
import ase
import numpy as np
import tsase
from random import random
from sys import argv

cell_size=5.0

atoms = ase.Atoms(argv[1])
atoms.set_cell((cell_size,cell_size,cell_size))
atoms.set_positions(np.random.random((len(atoms),3))*10.0)

tsase.io.write_con('matter1.con', atoms)

atoms.positions += random()*3.0
atoms.rattle(stdev=.915, seed=int(random()*2**30))
atoms.rotate_euler(center='COM', phi=random()*np.pi*2, theta=random()*np.pi, psi=random()*np.pi/2)

tsase.io.write_con('matter2.con', atoms)
