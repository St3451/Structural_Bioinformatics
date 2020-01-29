#!/usr/bin/env python3

import Bio.PDB as PDB

class CenterSelect(PDB.Select):
    def __init__(self, center):
        # Set center of molecule
        self.center = center

    def accept_atom(self, atom):
        """Accept atoms close to center"""
        diff = self.center - atom.get_vector()        
        dist = diff.norm()                            
        if dist < 10:
            return 1                            
        else:
            return 0


# Create parser
p = PDB.PDBParser()

# Get trypsin
s = p.get_structure("2PTC", "/home/lpp/BIOINFORMATICS/sb2019/week5/2ptc.pdb")
enzyme = s[0]["E"]

## Find CA center-of-mass
# CA counter
n = 0 

# Bio.PDB.Vector to calculate center-of-mass
atom_sum = PDB.Vector(0., 0., 0.)
for res in enzyme:
    # Check if res is sane (not water or missing CA)
    if res.has_id("CA"):
        # In-place elementwise addition
        atom_sum += res["CA"].get_vector()
        n += 1
# Divide by the number of atoms to get center-of-mass
# Note: ** multiplies Vector with Vector or scalar
com = atom_sum ** (1/n)

# For debugging
print(com)

## Create PDBIO object
io = PDB.PDBIO()

# Set the structure
io.set_structure(s)

# Filename for saving
outfile = '/home/lpp/BIOINFORMATICS/sb2019/week5/2ptc-center.pdb'

# Create CenterSelect object and pass the calculated center-of-mass
select = CenterSelect(com)

# Save structure using CenterSelect filter
io.save(outfile, select)

print("Printed %s" % outfile)