from ase.dimer import DimerControl, MinModeAtoms, MinModeTranslate
from ase.dimer import normalize, read_eigenmode
from ase.calculators.vasp.vasp2 import Vasp2
from ase.io import read, write
import numpy as np
from ase.constraints import FixBondLengths
import os
from ase.calculators.singlepoint import SinglePointCalculator


current_trajectory = "dimer_optimization.traj"
current_logfile = "dimer_optimization.log"
current_optimized = "dimer_optimized.traj"
current_control = "dimer_control.log"
current_eigenmode = "dimer_eigenmodes.dat"

atoms = read("vasp_optimized.traj")

# find out which C-H is breaking
symbols = atoms.get_chemical_symbols()
carbon_index = symbols.index("C")
hydrogen_index = None

# remove the constraint of FIxBondLength on C-H
new_constraints = []
for c in atoms.constraints:
    if isinstance(c, FixBondLengths):
        pairs = c.pairs.ravel()
        assert len(pairs) == 2
        if pairs[0] == carbon_index:
            hydrogen_index = pairs[1]
        elif pairs[1] == carbon_index:
            hydrogen_index = pairs[0]
    else:
        new_constraints.append(c)

# remove all FixBondLengths
atoms.set_constraint(new_constraints)

# start DIMER
initial_mode_method = "displacement"
displacement_method = 'vector'

# set up initial guess for DIMER
natoms = atoms.get_number_of_atoms()
displacement_vector = np.zeros(shape=(natoms, 3))
pos = atoms.get_positions()
_delta_x = normalize(pos[carbon_index] - pos[hydrogen_index]) * 0.01
displacement_vector[hydrogen_index] = _delta_x
mask = [False] * natoms
mask[hydrogen_index] = True

calc = Vasp2(encut=400, ispin=2, ediff=1.0e-7, nelm=120, nelmin=5, xc='pbe',
             kpts=(1, 1, 1), gamma=True, prec="N", algo="N", ismear=0, sigma=0.1,
             npar=8, lreal="Auto", lcharg=False, lwave=True, directory='calculator')


atoms.set_calculator(calc)

# Set up the dimer
d_control = DimerControl(initial_eigenmode_method='displacement',
                         displacement_method='vector',
                         logfile='dimer_control.log',
                         dimer_separation=0.008, max_num_rot=2, mask=mask,
                         eigenmode_logfile='dimer_eigenmodes.dat')

d_atoms = MinModeAtoms(atoms, d_control)
d_atoms.displace(displacement_vector=displacement_vector)


dim_rlx = MinModeTranslate(d_atoms, trajectory="dimer_opt.traj", logfile="dimer_opt.log")
dim_rlx.run(fmax=0.03)

# remember to add in the constraint
original_atoms = read("vasp_optimized.traj")
energy = atoms.get_potential_energy()
forces = atoms.get_forces()
a = atoms.copy()
a.set_constraint(original_atoms.constraints)
a.set_calculator(SinglePointCalculator(a, energy=energy, forces=forces))
write("dimer_optimized.traj", a)
