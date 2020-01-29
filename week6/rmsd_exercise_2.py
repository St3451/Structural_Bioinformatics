###### PART 2 - Work with PDB protein

from Bio import *
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser

from numpy import *
from numpy.linalg import *

p = PDBParser()
s = p.get_structure("1LCD", "/home/lpp/BIOINFORMATICS/sb2019/week6/1lcd.pdb")

chainA_x = s[0]["A"]
chainA_y = s[1]["A"]


def extract_c_alpha(chain):
    """Take a chain and return a matrix of C alpha, where rows are x,y,z and col are the vectors of C alpha"""
    c_alpha = [[],[],[]]
    residues = []
    for residue in chain:
        if is_aa(residue):
            if residue.has_id("CA"):                                 # Check if it is residue (probably redundant)
                residues.append(residue)
                ca = residue["CA"].get_vector()                      # Check if res is sane (not water or missing CA)
                ca_x = ca[0]
                ca_y = ca[1]
                ca_z = ca[2]
                c_alpha[0].append(ca_x)
                c_alpha[1].append(ca_y)
                c_alpha[2].append(ca_z)
    return (asmatrix(c_alpha))


x = extract_c_alpha(chainA_x)
y = extract_c_alpha(chainA_y)
print(x.shape, y.shape)    # 3 row, 51 col


# Calculate center of mass (we need to center x and y so move them to their center of mass)
center_of_mass_x = x.sum(1) / x.shape[1]
center_of_mass_y = y.sum(1) / y.shape[1]

# Center
x_centered = x - center_of_mass_x
y_centered = y - center_of_mass_y

r = y_centered * x_centered.T 
v, s, wt = svd(r)                
U = wt.T * v.T 

# Check for reflection  
z = diag([1,1,-1])
if round(det(U)) == -1:
    U = wt.T * z * v.T
    print(round(det(U)) == 1, "Roto-reflection fixed")

# RMSD 1
len_x = norm(x_centered)
len_y = norm(y_centered)
E_0 = len_x**2 + len_y**2
L_U = sum(s)

if round(det(U)) == 1:                    # if the U matrix is a rotation we calculate the RMSD
    RMSD_1 = sqrt( (E_0 - 2 * L_U) / x.shape[1])

print(RMSD_1)

# RMSD 2
RMSD_2 = sqrt(sum((norm(x_centered - U * y_centered, axis = 0))**2) / x.shape[1]) # with axis = 0 output the sum of each col
print(RMSD_2) 

print(U)

