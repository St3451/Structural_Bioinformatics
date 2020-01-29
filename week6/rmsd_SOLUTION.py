from numpy import *
from numpy.linalg import svd, det
# Test data
# Coordinate vectors along columns 
a=matrix([[ 18.92238689,  9.18841188,  8.70764463,  9.38130981, 8.53057997],
         [ 1.12391951,  0.8707568 ,  1.01214183,  0.59383894, 0.65155349],
         [ 0.46106398,  0.62858099, -0.02625641,  0.35264203, 0.53670857]], 'f')
b=matrix([[ 1.68739355,  1.38774297,  2.1959675 , 1.51248281,  1.70793414],
         [ 8.99726755,  8.73213223,  8.86804272, 8.31722197,  8.9924607 ],
         [ 1.1668153 ,  1.1135669 ,  1.02279055, 1.06534992,  0.54881902]], 'f')

print(a.sum(0))
print(a.shape)
# Code
def center(m):
   # Returns centered m
   # Calculate center of mass of x
   center_of_mass_m=m.sum(1)/m.shape[1]
   # Center m
   centered_m=m-center_of_mass_m
   return centered_m
def sup(x, y):
    # Nr of atoms
    N=x.shape[1]
    ########################
    # Center x and y
    x=center(x)
    y=center(y)
    # correlation matrix
    r=y*x.T
    # SVD of correlation matrix
    v, s, wt=svd(r)
    w=wt.T
    vt=v.T
    # Rotation matrix
    u=w*vt
    # Check for roto-reflection
    if det(u)<0:
        z=diag((1,1,-1))
        u=w*z*vt
        s[2]=-s[2]
    # Calculate RMSD
    e0=sum(x.A*x.A+y.A*y.A)
    print("E0 ", e0)
    rmsd=sqrt((e0-2*sum(s))/N)
    print('RMSD (svd) ', rmsd)

    print('RMSD (svd) ', rmsd)
    print("u estimated ")
    print(u)
    # Calculate RMSD from the coordinates 
    d=x-u*y
    d=d.A*d.A
    rmsd=sqrt(sum(d)/N)
    print('RMSD (real) ', rmsd)
    
if __name__=="__main__":
    # Run some stuff
    # 1. Test data
    sup(a,b)
    # 2. Protein
    # Parse PDB file
    from Bio.PDB import *
    p=PDBParser()
    s=p.get_structure("XXX", "1lcd.pdb")
    # Chain A in model 0
    chain1=s[0]["A"]
    # Chain A in model 1
    chain2=s[1]["A"]
    def get_coordinates(chain):
        # Get the CA coordinates of a cahin and return 3xn matrix
        coords=[]
        for res in chain:
            try:
                # Extract CA coordinate
                a=res["CA"]
                c=a.get_coord()
                coords.append(c)
            except:
                # No CA atom - skip
                pass
        # Turn coordinate list in 3xn numpy matrix
        coords=matrix(coords) # nx3
        coords=coords.T       # 3xn
        return coords
    # Get the 3xn coordinate matrices for the chains
    a=get_coordinates(chain1)
    b=get_coordinates(chain2)
    # Superimpose and done
    sup(a,b)