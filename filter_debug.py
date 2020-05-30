import sys
import pyrosetta
from pyrosetta import rosetta

import hydrogen_bonds
import metrics

# Run
#  ls OrigPDB_andStats/*.pdb | awk '{print "python filter_debug.py "$0}' | sed 's/\.pdb//' | sh

# no .pdb
pdb_file = sys.argv[1]

#uuid = sys.argv[2]

pyrosetta.init()
print("Current file..." + pdb_file)

#from pyrosetta.toolbox import cleanATOM
#cleanATOM(pdb_file+".pdb")
#pose = pyrosetta.pose_from_pdb(pdb_file+".clean.pdb")
pose = pyrosetta.pose_from_pdb(pdb_file+".pdb")

print(pose)
#print("> UUID: %s" % (uuid))
print("> Sequence: %s" % (pose.sequence()))


# Create a PyRosetta score function using:
from pyrosetta.teaching import *
sfxn = get_score_function(True)
print("> Total Energy: %f" % sfxn(pose))


# Hydrophobicity 
hydro_sasa = metrics.calc_buried_np_SASA(pose, residue_types=['A', 'F', 'I', 'L', 'M', 'V', 'W', 'Y']) / pose.size()
#print("> Hydrophobicity: %3.2f" % hydro_sasa)

# Solvent Accessible Surface Area
rsd_sasa = pyrosetta.rosetta.utility.vector1_double()
rsd_hydrophobic_sasa = pyrosetta.rosetta.utility.vector1_double()
rosetta.core.scoring.calc_per_res_hydrophobic_sasa(pose, rsd_sasa, rsd_hydrophobic_sasa, 1.4) 
#print("> Solvent Accessible Surface Area: %3.2f" % sum(rsd_sasa))
#print("> Hydrophobic SASA: %3.2f" % sum(rsd_hydrophobic_sasa))

# Polarity Hydrogen bond donors and acceptors
buried_unsat_acceptors, buried_unsat_donors = hydrogen_bonds.get_buhs_for_each_atom(pose)
#print("> HBond Acceptors: %d" % len(buried_unsat_acceptors))
#print(buried_unsat_acceptors)
#print("HBond Acceptors: %s" % buried_unsat_acceptors)
#print("> HBond Donors: %d" % len(buried_unsat_donors))
#print(buried_unsat_donors)
#print("HBond Donors: %3.2f" % buried_unsat_donors)

#atom_sasa = hydrogen_bonds.get_atom_sasa(pose)
#print(atom_sasa)
OUTPUT = ">> "
OUTPUT += "%s | " % (pose.sequence())
OUTPUT += "hydro | %3.2f | " % hydro_sasa
OUTPUT += "SASA | %3.2f | " % sum(rsd_sasa)
#OUTPUT += "atomSASA | %3.2f | " % sum(atom_sasa)
OUTPUT += "hydroSASA | %3.2f | " % sum(rsd_hydrophobic_sasa)
OUTPUT += "hbond-acc | %d | " % len(buried_unsat_acceptors)
OUTPUT += "hbond-donor | %d |" % len(buried_unsat_donors)
OUTPUT += "total-energy | %d \n" % sfxn(pose)
#print(OUTPUT)

outfile=pdb_file+".txt"
with open(outfile, "w") as outfile:
            outfile.write(OUTPUT)


#

