#Ex 3.4 BioPython
# Output all atoms within a sphere of 10 Angstrom of the center of trypsin to a separate PDB file. 
# Use the CA atoms to calculate the center.

from Bio import *
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser

# Define a class based on Bio.PDB.Select to output only atoms close to the center of mass 
class CenterSelect(Select):
    def __init__(self, center):
        self.center = center                             # Set center of molecule (will be determined by com)

    def accept_atom(self, atom):
        """Accept atoms close to center"""
        diff = self.center - atom.get_vector()           # Diff between com and the atom position (vector)
        dist = diff.norm()                               # .norm for length
        if dist < 10:
            return True                                  # It is the same of using 1 and 0 (used for filtering with PDBIO)
        else:
            return False

# Load structure
p = PDBParser()
s = p.get_structure("2PTC", "/home/lpp/BIOINFORMATICS/sb2019/week5/2ptc.pdb")
enzyme = s[0]["E"]                                                # structure 0 / chain E

# Find CA center-of-mass                                          # center-of-mass = sum(CA vectors) / (n. CA vectors)
n = 0                                                             # CA counter
atom_sum = Vector(0., 0., 0.)                                     # Bio.PDB.Vector to calculate center-of-mass (com)

for res in enzyme:
    if res.has_id("CA"):                                          # Check if res is sane (not water or missing CA)
        atom_sum += res["CA"].get_vector()                        # Add the CA vector to the atom_sum
        n += 1                                                    # Add 1 to the counter

com = atom_sum ** (1/n)                                           # ** multiplies Vector with Vector or scalar                               

print(com)                                                        # For debugging

## Output with selection
io = PDBIO()                                                      # Create PDBIO object
io.set_structure(s)                                               # Set the structure
outfile = '/home/lpp/BIOINFORMATICS/sb2019/week5/2ptc-center.pdb' # Filename for saving
 
# Create an object of class CenterSelect defined previously and pass the calculated center-of-mass
select = CenterSelect(com)

# Save structure using CenterSelect filter
io.save(outfile, select)
print("Printed %s" % outfile)