reference_no_irc/
No IRC, run in sequence 
python vasp_opt.py
python run_dimer.py
python cluster_opt.py
python write_final.py


reference_with_irc
with IRC, run in sequence 
python vasp_opt.py
python run_dimer.py
python irc_steepest_descent.py
python cluster_opt.py
python write_final.py


input.traj: the input structure, with a CH4 and cluster
           The CH4 contains ase.constraints
            FixBondLengths is used for the breaking C-H
            Hookean is used for the other three C-H bonds
vasp_opt.py: running constrained optimization, which reads "input.traj", generates "vasp_optimized.traj"
run_dimer.py: running DIMER, which reads "vasp_optimized.traj", generates "dimer_optimized.traj"
irc_steepest_descent.py: running IRC (steepest gradient), find the CH4 adsorbe state, generates  "sd_optimized.traj"
cluster_opt.py:  remove CH4 from "sd_optimzied.traj" (if found) otherwise from "dimer_optimized.traj", run local optimization and produces "cluster_optimized.traj"
write_final.py: compare free energies of "dimer_optimized.traj" and "cluster_optimized.traj", writes "optimized.traj"

