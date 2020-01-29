#Ex 3.1 BioPython
# Find all residues in trypsin (code 2PTC from pdb) that have more than two close contacts to the inhibitor.

from Bio import *
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser

parser = PDBParser()
structure=parser.get_structure("2PTC", "/home/lpp/BIOINFORMATICS/sb2019/week5/2ptc.pdb")

# Check if two residues are close
def is_close(res1, res2):
    """Return True if two residues have more than two
    contact that are less than 3.5 angstroms apart"""
    count = 0
    for atom1 in res1:
        for atom2 in res2:
                if abs(atom1 - atom2) < 3.5:
                    count += 1
                    if count > 2:
                        return True
    return False

chain_E = structure[0]['E']
chain_I = structure[0]['I']
close_residues = []

# Compare all enzyme and inhibitor residues
print("\nComparing chain E to chan I")
for residue_E in chain_E:
    for residue_I in chain_I:
        if is_close(residue_E, residue_I):
            close_residues.append((residue_E, residue_I))

# Output on the screen
print("%s pairs found" % len(close_residues))
for res1, res2 in close_residues:
    print("%3s-%3d is close to %3s-%3d" % (res1.get_resname(), res1.get_id()[1], res2.get_resname(), res2.get_id()[1]))




        
            
                
                
