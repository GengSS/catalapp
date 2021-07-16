from ase.io import read, write
import sys
import os
from ase.optimize.optimize import Optimizer
from ase.dimer import normalize, read_eigenmode
from ase.constraints import FixBondLengths
import numpy as np
from ase.calculators.vasp.vasp2 import Vasp2
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io.trajectory import Trajectory

class StopTag(Exception):
    pass

#steepest descent
class SteepestDescent(Optimizer):
    def __init__(self, atoms, restart=None, logfile='-', trajectory=None,
                 maxstep=None, dt=None, master=None):
        """Parameters:

        atoms: Atoms object
            The Atoms object to relax.

        restart: string
            Pickle file used to store hessian matrix. If set, file with
            such a name will be searched and hessian matrix stored will
            be used, if the file exists.

        trajectory: string
            Pickle file used to store trajectory of atomic movement.

        maxstep: float
            Used to set the maximum distance an atom can move per
            iteration (default value is 0.2 Angstroms).

        logfile: string
            Text file used to write summary information.

        master: boolean
            Defaults to None, which causes only rank 0 to save files.  If
            set to true,  this rank will save files.
        """
        Optimizer.__init__(self, atoms, restart, logfile, trajectory, master)

        if dt is not None:
            self.dt = dt
        else:
            self.dt = 0.2
        if maxstep is not None:
            self.maxstep = maxstep
        else:
            self.maxstep = 0.05

    def initialize(self):
        self.v = None

    def read(self):
        self.v, self.dt = self.load()

    def step(self, f=None):
        atoms = self.atoms

        if f is None:
            f = atoms.get_forces()

        self.v = 0.5 * self.dt * f

        r = atoms.get_positions()
 
        disp = self.dt*self.v 

        max_disp = np.sqrt(np.power(disp,2).sum(axis=1).max())
        if max_disp > self.maxstep:
            disp = disp/max_disp*self.maxstep
        atoms.set_positions(r + disp)
        self.dump((self.v, self.dt))


# read the dimer results
atoms_dimer = read("dimer_optimized.traj")
eigen_mode = read_eigenmode("dimer_eigenmodes.dat")

atoms = atoms_dimer.copy()
symbols = atoms.get_chemical_symbols()
carbon_index = symbols.index("C")
hydrogen_index = None

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

# remove the constraints
atoms.set_constraint(new_constraints)

s = np.sqrt(np.power(eigen_mode, 2).sum(axis=1))
disp = 0.05*eigen_mode/s.max()
a1 = atoms.copy()
a2 = atoms.copy()
pos0 = atoms.get_positions()
a1.set_positions(pos0+disp)
a2.set_positions(pos0-disp)

if a1.get_distance(carbon_index,hydrogen_index,mic=True) >  a2.get_distance(carbon_index,hydrogen_index,mic=True):
    initial = a2
else:
    initial = a1

write("sd_initial.traj",initial)

calc = Vasp2(encut=400, ispin=2, ediff=1.0e-5, nelm=120, nelmin=5, xc='pbe',
             kpts=(1, 1, 1), gamma=True, prec="N", algo="N", ismear=0, sigma=0.1,
             npar=8, lreal="Auto", lcharg=False, lwave=True, directory='calculator')

def dch(a=initial):
    d = initial.get_distance(carbon_index,hydrogen_index,mic=True)
    if d < 1.25:
        raise StopTag("{:.3f}".format(d))
    else:
        sys.stderr.write("{:.3f}\n".format(d))

initial.set_calculator(calc)
dyn = SteepestDescent(initial,trajectory="sd_opt.traj",logfile="sd_opt.log")
# we want SteepestDescent runs at least 200 steps or CH4 formed,
# to leave the TS region, where fmax is smalll
dyn.attach(dch,interval=1)
try:
    dyn.run(fmax=0.0000001, steps=200)
except StopTag:
    pass

trajectory = Trajectory('sd_opt.traj', mode='a', atoms=initial, master=True)
dyn = SteepestDescent(initial,trajectory=trajectory,logfile="sd_opt.log")
dyn.run(fmax=0.05)
energy = initial.get_potential_energy()
forces = initial.get_forces()
a = initial.copy()
a.set_constraint(atoms_dimer.constraints)
a.set_calculator(SinglePointCalculator(a, energy=energy, forces=forces))
write("sd_optimized.traj", a)
