#!/usr/bin/env python3

from Bio.PDB import PDBParser

# Create parser object
p = PDBParser()
# Get structure
s = p.get_structure("2PTC", "2PTC.pdb")

def is_close(res1, res2):
    """ Returns True if the two residues have
    more than two contacts that are less than 3.5 angstroms apart"""
    close_contacts = 0
    
    # Compare all atoms between the two residues
    for atom1 in res1:
        for atom2 in res2:
            if atom1 - atom2 < 3.5:     # Calculate distance using Atom class minus operator overload
                close_contacts += 1
                # Return True if more than two atom pairs were close
                if close_contacts > 2:
                    return True

    # Finished comparing atoms and there were fewer than 3, so return False
    return False

# List to store found residue pairs
close_pairs = []

# Compare all enzyme residue - inhibitor residue pairs
print("Comparing chain E to chain I")
for res1 in s[0]["E"]:
    for res2 in s[0]["I"]:
        if is_close(res1, res2):
            close_pairs.append((res1, res2))

# Output
print("Found %s close pairs" % len(close_pairs))
for r1, r2 in close_pairs:
    # Output format
    print("%3s%-3d is close to %3s%-3d" %
          (r1.get_resname(), r1.get_id()[1], r2.get_resname(), r2.get_id()[1]))
