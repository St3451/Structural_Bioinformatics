##### PART 1 - Work with toy proteins coordinate

from Bio import *
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser

from numpy import *
from numpy.linalg import *

x = matrix([[ 18.92238689,  9.18841188,  8.70764463,  9.38130981,  8.53057997],      # x
        [ 1.12391951,  0.8707568 ,  1.01214183,  0.59383894,  0.65155349],           # y 
        [ 0.46106398,  0.62858099, -0.02625641,  0.35264203,  0.53670857]], 'f')     # z

y = matrix([[ 1.68739355,  1.38774297,  2.1959675 ,  1.51248281,  1.70793414],
        [ 8.99726755,  8.73213223,  8.86804272,  8.31722197,  8.9924607 ],
        [ 1.1668153 ,  1.1135669 ,  1.02279055,  1.06534992,  0.54881902]], 'f')

print(x.shape) # Check the shape of the popypeptide, should return a 3 by n matrix where n equal the number of residues

# Define a function that superimpose two chains by applying a rotation matrix U that minimize the RMSD
def calc_RMSD(chain_a, chain_b):
    """Take two sets of vectors representing the coordinates of the alpha carbons and return a tuple
     of three elements containing the two RMSD calculated by two methods and the rotation matrix U"""

    # Call the previously define function to extract the C alpha
    x, y = (chain_a), (chain_b)

    # Calculate center of mass (we need to centre x and y, so move them to their center of mass)
    center_of_mass_x = x.sum(1) / x.shape[1]
    center_of_mass_y = y.sum(1) / y.shape[1]
    print("center of mass XXXXX", center_of_mass_x)

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
    if round(det(U)) == -1:                       # here I also should change the third value of the s vector
        U = wt.T * z * v.T
        print(round(det(U)) == 1, "Roto-reflection fixed")

    # Use the closed formula to calculate the RMSD (purely using linear algebra)
    len_x = norm(x_centered)                                   # with x.A we can use x as array (element wise moltiplication) 
    len_y = norm(y_centered)
    E_0 = len_x**2 + len_y**2
    L_U = sum(s)
    if round(det(U)) == 1:
        RMSD_1 = sqrt( (E_0 - 2 * L_U) / x.shape[1])

    # Use the rotated coordinates to calculate the RMSD
    RMSD_2 = sqrt(sum((norm(x_centered - U * y_centered, axis = 0))**2) / x.shape[1]) # with axis = 0 output the sum of each col
    
    # Return the RMSD calculated by the two different methods and the rotation matrix U
    return(RMSD_1, RMSD_2, U)

part1 = calc_RMSD(x, y)
print("\n\nRMSD by formula and by rotated coordinates:\n", part1[0],";", part1[1],"\n")
print("Rotation matrix:\n", calc_RMSD(x, y)[2],"\n\n")


###### PART 2 - Work with PDB protein

# Load the PDB structure and select the models (0 and 1) and chain A
p = PDBParser(QUIET = True)
s = p.get_structure("1LCD", "/home/lpp/BIOINFORMATICS/sb2019/week6/1lcd.pdb")
chain_Ax = s[0]["A"]
chain_Ay = s[1]["A"]

# Define a function that extract the C alpha from a chain
def extract_c_alpha(chain):
    """Take a chain as input and return a matrix of alpha C coordinates, where rows are 
    x,y,z axis and columns are the vectors of alpha C in the three dimensional space"""
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

part2 = calc_RMSD(extract_c_alpha(chain_Ax), extract_c_alpha(chain_Ay))

print("\n\nRMSD by formula and by rotated coordinates:\n", part2[0], ";", part2[1], "\n")
print("Rotation matrix:\n", part2[2])