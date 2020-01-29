#!/usr/bin/env python3

import Bio.PDB as PDB

# define a class based on Bio.PDB.Select to output chain E only
class TrypsinSelect(PDB.Select):
    def accept_chain(self, chain):
        if chain.get_id() == "E":
            return True
        else:
            return False


# Load structure
p = PDB.PDBParser()
s = p.get_structure("2PTC", "2PTC.pdb")

# Create PDBIO object
io = PDB.PDBIO()

# Set the structure
io.set_structure(s)

# Filename to save to
outfile = 'out.pdb'

# Create an object of class TrypsinSelect defined previously
select = TrypsinSelect()

# Save the structure using the TrypsinSelect object as filter
io.save(outfile, select)

print("Printed %s" % outfile)
