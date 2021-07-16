from ase.calculators.vasp.vasp2 import Vasp2
from ase.optimize import BFGS
from ase.io import read, write


initial_structure = "input.traj"
calc = Vasp2(encut=400, ispin=2, ediff=1.0e-5, nelm=120, nelmin=5, xc='pbe',
             kpts=(1, 1, 1), gamma=True, prec="N", algo="N", ismear=0, sigma=0.1,
             npar=8, lreal="Auto", lcharg=False, lwave=True, directory='calculator')

atoms = read(initial_structure)
if "initial_magmoms" in atoms.arrays.keys():
    del atoms.arrays['initial_magmoms']

atoms.set_calculator(calc)
opt = BFGS(atoms, trajectory="vasp_opt.traj",logfile="vasp_opt.log")
opt.run(fmax=0.01)
write("vasp_optimized.traj", atoms)
