###### PART 2 - Work with PDB protein

# Import the various modules needed
from Bio import *
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser

from numpy import *
from numpy.linalg import *

# Load the PDB structure and select the models (0 and 1) and chain A
p = PDBParser()
s = p.get_structure("1LCD", "/home/lpp/BIOINFORMATICS/sb2019/week6/1lcd.pdb")
chainA_x = s[0]["A"]
chainA_y = s[1]["A"]

# Define a function that extract the C alpha from a chain
def extract_c_alpha(chain):
    """Take a chain as input and return a matrix of C alpha vectors, where rows are 
    x,y,z axis and columns are the vectors of C alpha in the three dimensional space"""
    c_alpha = [[],[],[]]
    # Iterate through all the residues in the chain
    for residue in chain:
        # Check if it is a residue (probably redundant)
        if is_aa(residue):
            # Check if res is sane (not water or missing CA)
            if residue.has_id("CA"):      
                # Get and store the coordinate of the C alpha                           
                ca = residue["CA"].get_vector()                      
                ca_x = ca[0]
                ca_y = ca[1]
                ca_z = ca[2]
                c_alpha[0].append(ca_x)
                c_alpha[1].append(ca_y)
                c_alpha[2].append(ca_z)
    # return the set of coordinates of the C alpha of the chain as 3 by n matrix, where n are the number of residues
    return (asmatrix(c_alpha))

# Define a function that superimpose two chains by applying a rotation matrix U that minimize the RMSD
def calc_RMSD(chain_a, chain_b):

    # Call the previously define function to extract the C alpha
    x, y = extract_c_alpha(chain_a), extract_c_alpha(chain_b)
    print("x, y matrices shape check: ", x.shape, y.shape) # should be 3 row, n col

    # Calculate center of mass (we need to centre x and y, so move them to their center of mass)
    center_of_mass_x = x.sum(1) / x.shape[1]
    center_of_mass_y = y.sum(1) / y.shape[1]

    # Center
    x_centered = x - center_of_mass_x
    y_centered = y - center_of_mass_y

    # Calculate R, the correlation matrix of x and y
    r = y_centered * x_centered.T 

    # Apply the SVD to R in order to find the rotation matrix U to apply to Y in order to minimize the RMSD
    v, s, wt = svd(r)                
    U = wt.T * v.T 

    # Check for reflection  
    z = diag([1,1,-1])
    if round(det(U)) == -1:
        U = wt.T * z * v.T
        print(round(det(U)) == 1, "Roto-reflection fixed")

    # Use the closed formula to calculate the RMSD
    len_x = norm(x_centered)
    len_y = norm(y_centered)
    E_0 = len_x**2 + len_y**2
    L_U = sum(s)
    if round(det(U)) == 1:
        RMSD_1 = sqrt( (E_0 - 2 * L_U) / x.shape[1])

    # Use the rotated coordinates to calculate the RMSD
    RMSD_2 = sqrt(sum((norm(x_centered - U * y_centered, axis = 0))**2) / x.shape[1]) # with axis = 0 output the sum of each col
    
    # Return the RMSD calculated by the two different methods and the rotation matrix U
    return(RMSD_1, RMSD_2, U)

print("RMSD by formula and by rotated coordinates: ", calc_RMSD(chainA_x, chainA_y)[0], calc_RMSD(chainA_x, chainA_y)[1])
print("Rotation matrix:\n", calc_RMSD(chainA_x, chainA_y)[2])