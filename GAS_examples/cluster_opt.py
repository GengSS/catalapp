from ase.calculators.vasp.vasp2 import Vasp2
from ase.optimize import BFGS
from ase.io.trajectory import Trajectory
from ase.io import read, write
import sys
import os
import numpy as np
from ase.constraints import FixAtoms
from ase.calculators.singlepoint import SinglePointCalculator

# first, create the input file from DIMER or IRC
if os.path.isfile("sd_optimized.traj"):
    # if we have IRC result, we choose IRC results
    a0 = read("sd_optimized.traj")
else:
    a0 = read("dimer_optimized.traj")
ch4_indexes = []
for c in a0.constraints:
    if isinstance(c, FixAtoms):
        continue
    for i in c.get_indices():
        if i not in ch4_indexes:
            ch4_indexes.append(i)
a0.set_constraint()
del a0[ch4_indexes]
write("cluster_initial.traj",a0)

calc = Vasp2(encut=400, ispin=2, ediff=1.0e-5, nelm=120, nelmin=5, xc='pbe',
             kpts=(1, 1, 1), gamma=True, prec="N", algo="N", ismear=0, sigma=0.1,
             npar=8, lreal="Auto", lcharg=False, lwave=True, directory='calculator')


atoms = read("cluster_initial.traj")
atoms.set_calculator(calc)
dyn = BFGS(atoms, trajectory="cluster_opt.traj",logfile="cluster_opt.log")
dyn.run(fmax=0.01)
write("cluster_optimized.traj", atoms)
