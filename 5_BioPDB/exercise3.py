#Ex 3.3 BioPython
# Output Trypsine to a separate PDB file.

from Bio import *
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser

# Define a class based on Bio.PDB.Select to output chain E only
class E_Select(Select):
    def accept_chain(self, chain):
        if chain.get_id()=='E':
            return True
        else:
            return False

# Load structure
p = PDBParser()
s = p.get_structure("trypsin", "2ptc.pdb")

# Create PDBIO object
io = PDBIO()

# Set the structure
io.set_structure(s)

# Filename to save to
outfile = 'out.pdb'

# Create an object of class E_Select defined previously
select = E_Select()

# Save the structure using the E_Select object as filter
io.save(outfile, select)
print("Printed %s" % outfile)

