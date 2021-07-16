from ase.io import read, write
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator

# define the optimization target, which is the free energy of the TS
def get_the_optimization_target(a=None):
    e = a.get_potential_energy()
    symbols = a.get_chemical_symbols()
    if symbols.count('C') == 1:
        mu = {'H': -4.067} # chemical potential of H
        e += 4.0*-4.067
        e -= -24.031 #This is the electronic energy of Methane
    else:
        mu = {'H': -4.067}
    for ii, s in enumerate(a.get_chemical_symbols()):
        e = e - mu.get(s, 0.0)
    return e

cluster = read("cluster_optimized.traj")
t1 = get_the_optimization_target(cluster)

ts = read("dimer_optimized.traj")
t2 = get_the_optimization_target(ts)

# E is a artifical target, that BH will optimize (minimize)
E = max([t1, t2])

atoms = ts.copy()
atoms.set_calculator(SinglePointCalculator(atoms, energy=E, forces=np.zeros_like(atoms.get_positions())))


# BH code only cares about optimized.traj
write("optimized.traj",atoms)
