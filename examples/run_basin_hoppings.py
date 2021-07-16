import sys
import os
path=os.path.realpath("../gcbh/")
sys.path.insert(0,path)

from ase.io import read
from gcbasin3_for_PtH import GrandCanonicalBasinHopping
from pygcga import mutation_atoms,add_molecule_on_cluster,remove_one_adsorbate

if not os.path.isfile("NiAlH_jea.eam.alloy"):
    print("You should download file 'NiAlH_jea.eam.alloy' into this folder\n from https://www.ctcms.nist.gov/potentials/system/Al before runing this example")
    exit()

filescopied=['opt.py','NiAlH_jea.eam.alloy'] # files required to complete an optimization

atoms=read("Ni17.xyz")
bh_run=GrandCanonicalBasinHopping(atoms=atoms, bash_script="optimize.sh", files_to_copied=filescopied,
                                  restart=True, chemical_potential="chemical_potentials.dat")
bh_run.add_modifier(add_molecule_on_cluster, name="add", outer_sphere=False, distribution=1,
                    molecule="H", metal=["Ni"])
#bh_run.add_modifier(add_molecule_on_cluster, name="add_H2", outer_sphere=False, distribution=1,
#                    molecule="H2", metal=["Al"],weight=0.7,anchor_atom='H')
bh_run.add_modifier(remove_one_adsorbate, name='remove', molecule="H")
bh_run.add_modifier(mutation_atoms, name='mutation', weight=3.0)
bh_run.run(20)
