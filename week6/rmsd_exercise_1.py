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

print(x.shape)

###### PART 1 - Work with toy proteins coordinate

#### Calculate RMSD with formula (we will not use U, because we are not rotating any protein, we are just using close formula)

# Calculate center of mass (we need to center x and y so move them to their center of mass)
center_of_mass_x = x.sum(1) / x.shape[1]
center_of_mass_y = y.sum(1) / y.shape[1]

# Center
x_centered = x - center_of_mass_x
y_centered = y - center_of_mass_y

# we need to multiply the two matrix r = y * x.t 
r = y_centered * x_centered.T 

# Singular Value Decomposition, in order to find U (rotational matrix to apply to y)
v, s, wt = svd(r)                # if Q = ortogonal -> Q * Q.T = Q.T * Q = I0
print(v, "v")                    # v and wt are two ortogonal matrices, wt.T * v.t = U = identity matrix that dont change the value but only rotate
print(s, "s")                    # s gives only a vector of the diagonals
print(wt, "wt")

U = wt.T * v.T 

# Z = diag([1,1,-1])

#Check for reflection
        #if round(det(U)) ==  -1 :
                # U = wt.T * Z * v.T

# Rotate Y by applying U
        # Y_rotated = U * y

## Calculate RMSD = sqrt(E_O - 2 * L_U)

# E_0
norm_x = norm(x_centered)
norm_y = norm(y_centered)
E_0 = norm_x**2 + norm_y**2

# L_U         
L_U = sum(s)                              # s = diagonal elements = singular values

# RMSD
if round(det(U)) == 1:                    # if the U matrix is a rotation we calculate the RMSD
        RMSD_1 = sqrt( (E_0 - 2 * L_U) / x.shape[1])
# elif round(det(U)) == -1:               # if the U matrix is a reflection we change U multiply by Z
print(RMSD_1)


#### Calculate RMSD with rotation formula                            #x.shape[0] = rows, x.shape[1] = col
RMSD_2 = sqrt(sum((norm(x_centered - U * y_centered, axis = 0))**2) / x.shape[1]) # with axis = 0 output the sum of each col
print(RMSD_2)                                                                     # shape[0] = rows
                                                                                  # anyway I don't need to specify it
print(U)