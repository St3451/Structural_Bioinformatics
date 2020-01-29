#!/usr/bin/env python3

import Bio.PDB as PDB
from Bio.PDB import Vector, calc_dihedral 

# Create parser object
p = PDB.PDBParser()
# Get structure
s = p.get_structure("2PTC", "2PTC.pdb")

#####################################
# First using the Polypeptide class #
#####################################

def get_phi_psi(structure):
    """
    Calculate phi,psi dihedral angles and return lists.
    Uses the polypeptide class. A polypeptide is a normal (index zero-based)
    list of residues that are covalently bonded to the next and previous
    list element. No waters or ions.
    """
    # Create a list of  polypeptide objects
    ppb = PDB.PPBuilder()
    pp_list = ppb.build_peptides(structure)

    # Get phi and psi angles
    phi_list = []
    psi_list = []
    # Iterate over polypeptide molecules
    for pp in pp_list:
        # Calculate phi and psi angles and unpack list and tuple
        for phi, psi in pp.get_phi_psi_list():
            # put them in the lists
            phi_list.append(phi)
            psi_list.append(psi)

    return phi_list, psi_list

# Print (phi, psi) of residue 2
phi_list, psi_list = get_phi_psi(s)
print("Phi angles:")
print(phi_list[1])
print("Psi angles:")
print(psi_list[1])

##############################
# Second, own implementation #
##############################

def calc_phi_psi(res1, res2, res3):
    """Return the list of phi/psi dihedral angles."""
    n = res2['N'].get_vector()
    ca = res2['CA'].get_vector()
    c = res2['C'].get_vector()
    # Phi
    cp = res1['C'].get_vector()
    phi = calc_dihedral(cp, n, ca, c)
    # Psi
    nn = res3['N'].get_vector()
    psi = calc_dihedral(n, ca, c, nn)
    # Return phi, psi tuple
    return (phi, psi)

# Print (phi, psi) of residue 2
# Get all residues
r=list(s.get_residues())
print(calc_phi_psi(r[0], r[1], r[2]))


