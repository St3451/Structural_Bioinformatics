from Bio import *
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser
from numpy import *
from numpy.linalg import *
import matplotlib.pyplot as plt

def extract_structures(directory):
    """
    Return a dictionary with entries for all protein files in the directory. 
    Each entry contain the parsed structure of the protein file.
    """
    p = PDBParser(QUIET = True) # Quiet ignore warnings 
    list_structures = []
    # Add to a dictionary the structure of each protein file and the name of the protein as key
    for filename in os.listdir(directory):
        print("Parsing %s..." % filename)
        # Try to parse the structure
        try:                                                                
            s = p.get_structure(filename, os.path.join(directory, filename))  
            list_structures.append(s)                                      
        except:
            print(filename, "can't be parsed")
    return(list_structures)

def get_centered_sidechain(res):
    """
    Returns the centered coordinates of the side chain  
    as a 3 by n numpy matrix.
    """
#    print("\n\n")
    coords = []
    # Select the atoms that I don't want to consider for the RMSD calculation
    atoms_to_avoid = ["O","2HA","1HA","3HB","3HG1","HH","HH2","HZ3","HZ2","HE3","HG1","HZ","HD2","HD1","3HE","3HZ","2HZ","1HZ","2HE","1HE","3HD2","HB","3HD1","2HD1","1HD1","3HG2","2HG2","1HG2","2HG1","1HG1","OXT","N","C","HA","1HB","2HB","1HG","2HG","1HD","2HD","H","HE","HE1","HE2","1HH1","2HH1","1HH2","2HH2","1HD2","2HD2","HG","1HE2","2HE2"]

    # Add coordinates excluding the atoms to avoid
    for atom in res.get_atoms():
        #print(atom.get_id())
        if atom.get_id() not in atoms_to_avoid:       
#            print(atom.get_id())
            coords.append(atom.get_coord())    
    # Turn coordinate list in 3 by n numpy matrix
    coords = matrix(coords) 
    coords = coords.T
    # Calculate side chain center of mass
    center = coords.sum(1)/coords.shape[1]
    # Calculate centered side chain coordinates
    centered_coords = coords - center   
    # Return center
#    print(centered_coords)
    return (centered_coords)

def calc_RMSD(x, y):
    """
    Superimpose two residues, represented by two 3 by n numpy
    matrices, and return their minimum RMSD.
    """
    # Nr of atoms
    N = x.shape[1]
    # Correlation matrix
    r = y*x.T
    # SVD of correlation matrix
    v,s,wt = svd(r)
    w = wt.T
    vt = v.T
    # Rotation matrix
    u = w*vt
    # Check for roto-reflection
    if det(u) < 0:
        z = diag((1,1,-1))
        u = w*z*vt
        s[2] = -s[2]
    # Calculate RMSD
    e0 = sum(x.A*x.A+y.A*y.A)
#    print("E0 ", e0)
    rmsd = sqrt((e0-2*sum(s))/N)
#   print('RMSD (svd) ', rmsd)
#    print("u estimated ")
#    print(u)
    # Calculate RMSD from the coordinates 
#    d = x-u*y
#    d = d.A*d.A
#    rmsd = sqrt(sum(d)/N)
#    print('RMSD (real) ', rmsd, "\n")
    return(rmsd)

def extract_RMSD_list(aa, list_structures):
    """
    Given a list of structures and an amino acid name, extract 1000 pairs
    of the selected amino acids, randomly and from two different proteins.
    """
    print("Computing %s RMSD scores..." % aa)
    n = 0
#    cycle = 0
    rmsd_list = []
    # Repeat the while loop until 1000 residue pairs, with equal number of atoms, are selected
    while n < 1000:
#        print("\n", cycle, "\n")
#        cycle += 1
#        print("\n", n, "\n")
        # Select two random different proteins
        i,j = random.choice(len(list_structures), size = 2, replace = False)   
        s1 = list_structures[i]
        s2 = list_structures[j]    
        # Extract all the selected amino acid from the two selected proteins
        list_res1 = []
        for res in s1[0].get_residues(): 
            if res.get_resname() == aa:
                list_res1.append(res)
        list_res2 = []
        for res in s2[0].get_residues(): 
            if res.get_resname() == aa:
                list_res2.append(res)    
        # Check if the amino acid lists are not empty
        if len(list_res1) != 0 and len(list_res2) != 0:
            # Select randomly one amino acid from each of the two proteins
            i1 = random.choice(len(list_res1))
            i2 = random.choice(len(list_res2))
            res1 = list_res1[i1]
            res2 = list_res2[i2]
            # Get the centered coordinates of the two side_chains
            sd1 = get_centered_sidechain(res1)
            sd2 = get_centered_sidechain(res2)
            # Check if they have both the same number of atoms
            if sd1.shape == sd2.shape:
                # Clculate RMSD and add it to a list
                rmsd = calc_RMSD(sd1, sd2)
                rmsd_list.append(rmsd)
                # Add 1 to the counter
                n += 1 
            # Discharge the pair if the number of atoms are not the same
            else:
                pass
        # If the selected amino acid is not found, 
        else:
            pass
#            print("%s not found, picking up other two proteins.." % aa)

        #except:
        #    pass
#    print(len(rmsd_list))
    return (rmsd_list)

def save_histogram(rmsd_list, aa):
    """
    Make histogram and save as png file.
    """
    # The histogram of the data
    n, bins, patches = plt.hist(rmsd_list, 50, density=True)
    # Labels
    plt.xlabel('RMSD')
    plt.ylabel('Probability')
    plt.title('Histogram of RMSD for %s' % aa)
    # x and y limits of plot
    plt.xlim(0, 1.4)
    plt.ylim(0, 40)
    plt.grid(True)
    # Save
    plt.savefig("/home/lpp/BIOINFORMATICS/sb2019/ESAME/PROTEIN/rmsd_hist_"+aa+".png")
    # Clear canvas for next plot
    plt.clf()

structures_list = extract_structures("/home/lpp/BIOINFORMATICS/sb2019/ESAME/PROTEIN/top500H")

aa_names = ["ARG","ASN","ASP","CYS","GLN","GLU","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]
aa_names3 = []

rmsd_dict = {}
for aa in aa_names:
    rmsd_dict[aa] = extract_RMSD_list(aa, structures_list)

for element in rmsd_dict:
    save_histogram(rmsd_dict[element], element)
    print("mean %s:" % element, mean(rmsd_dict[element]))





