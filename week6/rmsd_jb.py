import numpy as np
from numpy import * # for matrix, sum, sqrt, diag, array,...
from numpy.linalg import svd, det # singular value decomposition, determinant

x=matrix([[ 18.92238689,  9.18841188,  8.70764463,  9.38130981,  8.53057997],
        [ 1.12391951,  0.8707568 ,  1.01214183,  0.59383894,  0.65155349],
        [ 0.46106398,  0.62858099, -0.02625641,  0.35264203,  0.53670857]], 'f') #floating point
y=matrix([[ 1.68739355,  1.38774297,  2.1959675 ,  1.51248281,  1.70793414],
        [ 8.99726755,  8.73213223,  8.86804272,  8.31722197,  8.9924607 ],
        [ 1.1668153 ,  1.1135669 ,  1.02279055,  1.06534992,  0.54881902]], 'f')

# # X,Y are {3,N} matrices
# Move X, Y to center of mass

# 1. To center the 3 by N matrix x to its center of mass:

# # Calculate center of mass of x and y
center_of_mass_x = x.sum(1)/x.shape[1] #array.shape() gives the dimensions of it, [1] is the number of column (python starts at 0), here = 3
center_of_mass_y = y.sum(1)/y.shape[1] #array.sum(1) sum the vectors of the array by row, to get the vector of the center of mass.

# # Center x and y 
centered_x = x-center_of_mass_x
centered_y = y-center_of_mass_y

# R=YX t
x_t = centered_x.T
R = centered_y*x_t

# # Singular Value Decomposition
# V, S, W t =SVD(R)
v, s, w_t = svd(R)
print(v, "V")
print(s, "s") #gives only a vector of the diagonals
print(w_t, "w_t")

# U=WV t
w = w_t.T
v_t = v.T
u = w*v_t

#DO THE DOT PRODUCT FOR THE NORM OF X AGAINST ITSELF -> TRANSFORM INTO ARRAY
norm_x = np.linalg.norm(centered_x)
norm_y = np.linalg.norm(centered_y)
E_0 = norm_x**2 + norm_y**2

if round(det(u)) == 1:
    RMSD = sqrt( (E_0 - 2*sum(s)) / x.shape[1]) #number of columns
elif round(det(u)) == -1:
    "calculate for the roto-reflection"

print(RMSD)


# norm_x = sqrt(sum(centered_x.A * centered_x.A)) #this works!
# norm_y = sqrt(sum(centered_y.A * centered_y.A))


#NOT LOOPING BECOZ SUPER COMPUTATIONALLY EXPENSIVE
# for i in range(x.shape[1]):
#     norm_x = sqrt((x[0, i]**2) + (x[1, i]**2) + (x[2, i]**2))
#     norm_y = sqrt((y[0, i]**2) + (y[1, i]**2) + (y[2, i]**2))
#     E_0 += norm_x**2 + norm_y**2
# print(E_0)
# print(E_U)
# # Check for reflection
# Z=diag(1,1,­1)
# if det(U)==­1:
# U=WZV t
# Rotate Y by applying U
# Calculate RMSD (either in real space or by formula)