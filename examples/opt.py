from ase.io import read,write
from ase.calculators.eam import EAM
from ase.optimize import BFGS


calc=EAM(potential='NiAlH_jea.eam.alloy')

atoms=read("input.traj") # input is always input.traj
atoms.set_calculator(calc)

dyn=BFGS(atoms,logfile="opt.log")
dyn.run(fmax=0.05)

write("optimized.traj",atoms) # mandatory
