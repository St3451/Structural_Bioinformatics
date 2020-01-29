#Ex 3.2 BioPython
# Print (phi, psi) angles for both proteins, using the Polypeptide class and using your own code (use calc_dihedral).

from Bio import *
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser

parser = PDBParser()
structure=parser.get_structure("Trypsin", "2ptc.pdb")

# Print phi and psi angles for both proteins using Polypeptide class
ppb = PPBuilder()
pp_list = ppb.build_peptides(structure)                                 # Create a list of polypeptide objects

for pp in pp_list:                                              # Just to check what we just created
    print(pp)
    for residue in pp:
        print(residue)

def find_phi_psi(structure):
    """Calculate phi,psi dihedral angles and return lists.
    Uses the polypeptide class. A polypeptide is a normal (index zero-based)
    list of residues that are covalently bonded to the next and previous
    list element. No waters or ions."""
    phi_list = []
    psi_list = []

    for pp in pp_list:
        for phi, psi in pp.get_phi_psi_list():                  # Calculate phi and psi angles and unpack list and tuple
            phi_list.append(phi)
            psi_list.append(psi)

    return phi_list, psi_list

phi_list, psi_list = find_phi_psi(s)
print("Phi angle:")
print(phi_list[1])
print("Psi angle:")
print(psi_list[1])

# Print phi and psi by own implementation
def calc_psi_phi(res1, res2, res3):                   # we need a function that take three residues (r-1, r, r+1):
    n = res2['N'].get_vector()                        #               r-1 ;       r      ; r+1                
    ca = res2['CA'].get_vector()                      #                C  ; N, Calpha, C ;  N
    c = res2['C'].get_vector()                        # from r-1, r   ( 4 atoms for phi  )          
    # Calculate phi                                   # from r, r+1       (  4 atoms for psi )            
    cp = res1['C'].get_vector()                       # cp = carbon previous residue (r-1)
    phi = calc_dihedral(cp, n, ca, c)
    # Calculate psi
    nn = res3['N'].get_vector()                       # nn = nitrogen next residue (r-1)
    psi = calc_dihedral(n, ca, c, nn)

    return (phi, psi)

res_list = list(s.get_residues())
print(calc_psi_phi(res_list[0], res_list[1], res_list[2]))





