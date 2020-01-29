

from Bio import *
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser

from numpy import *
from numpy.linalg import *

#### Center a chain

x = matrix([[ 1,  2,  3, 5],      # x
            [ 2,  4,  3, 2],           # y 
            [ 1,  5,  2, 1]], 'f')     # z

y = matrix([[ 3,  1,  4, 4],
            [ 4,  1,  2, 2],
            [ 3,  2,  2, 3]], 'f')

print(x.shape) # Check the shape of the popypeptide, should return a 3 by n matrix where n equal the number of residues


center_of_mass_x = x.sum(1) / x.shape[1]  # center of mass is a 1 by 3 matrix = [sum of x of all residues / n. of residues..],[y],[z] 
center_of_mass_y = y.sum(1) / y.shape[1]

print(x.shape[1])
print(x.sum(1))
print(center_of_mass_x)

x_centered = x - center_of_mass_x
y_centered = y - center_of_mass_y

print(x_centered)

lista = [1,2,3,4,5]
for n in range(len(lista)):
    print(lista[n])