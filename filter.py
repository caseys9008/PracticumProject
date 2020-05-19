import sys
import pyrosetta
from pyrosetta import rosetta

import hydrogen_bonds
import metrics

pdb_file = sys.argv[1]
uuid = sys.argv[2]

pyrosetta.init()
pose = pyrosetta.pose_from_pdb(pdb_file)
print(pose)
print("> UUID: %s" % (uuid))
print("> Sequence: %s" % (pose.sequence()))

# Create a PyRosetta score function using:
from pyrosetta.teaching import *
sfxn = get_score_function(True)
print("> Total Energy: %f" % sfxn(pose))


# Hydrophobicity 
hydro_sasa = metrics.calc_buried_np_SASA(pose, residue_types=['A', 'F', 'I', 'L', 'M', 'V', 'W', 'Y']) / pose.size()
print("> Hydrophobicity: %3.2f" % hydro_sasa)

# Solvent Accessible Surface Area
rsd_sasa = pyrosetta.rosetta.utility.vector1_double()
rsd_hydrophobic_sasa = pyrosetta.rosetta.utility.vector1_double()
rosetta.core.scoring.calc_per_res_hydrophobic_sasa(pose, rsd_sasa, rsd_hydrophobic_sasa, 1.4) 
print("> Solvent Accessible Surface Area: %3.2f" % sum(rsd_sasa))
print("> Hydrophobic SASA: %3.2f" % sum(rsd_hydrophobic_sasa))

# Polarity Hydrogen bond donors and acceptors
buried_unsat_acceptors, buried_unsat_donors = hydrogen_bonds.get_buhs_for_each_atom(pose)
print("> HBond Acceptors: %d" % len(buried_unsat_acceptors))
#print(buried_unsat_acceptors)
#print("HBond Acceptors: %s" % buried_unsat_acceptors)
print("> HBond Donors: %d" % len(buried_unsat_donors))
#print(buried_unsat_donors)
#print("HBond Donors: %3.2f" % buried_unsat_donors)

#atom_sasa = hydrogen_bonds.get_atom_sasa(pose)
#print(atom_sasa)
