#### The implementation of the RMSD algorithm (and most of the rest of the code) is 
#### based on weekly exercises solutions provided by our Structural Bioinformatics
#### professor Thomas Hamelryck (Associate professor, Computational and RNA Biology; 
#### University of Copenhagen).

from Bio import *
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser
from numpy import *
from numpy.linalg import *
import matplotlib.pyplot as plt

def extract_structures(directory):
    """
    Return a list which elements correspond to the parsed 
    protein structures of PDB files in a given directory.
    """
    # Create parser object
    p = PDBParser(QUIET = True)
    list_structures = []
    # Iterate through all the files in a given directory
    for filename in os.listdir(directory):
        print("Parsing %s..." % filename)
        # Try to parse the structure
        try:       
            # Get structure object s                                                         
            s = p.get_structure(filename, os.path.join(directory, filename))  
            # Append the structure object to a list
            list_structures.append(s)                                      
        except:
            print(filename, "can't be parsed")      
    return(list_structures)

def get_centered_sidechain(res):
    """
    Returns the centered coordinates of the  
    side chain as a 3 by n numpy matrix.
    """
    coords = []
    # Select the atoms that I do not want to consider for the RMSD calculation
    atoms_to_avoid = ["O","2HA","1HA","3HB","3HG1","HH","HH2","HZ3","HZ2",
                     "HE3","HG1","HZ","HD2","HD1","3HE","3HZ","2HZ","1HZ",
                     "2HE","1HE","3HD2","HB","3HD1","2HD1","1HD1","3HG2",   
                     "2HG2","1HG2","2HG1","1HG1","OXT","N","C","HA","1HB",  
                     "2HB","1HG","2HG","1HD","2HD","H","HE","HE1","HE2",
                     "1HH1","2HH1","1HH2","2HH2","1HD2","2HD2","HG","1HE2",
                     "2HE2"]
    # Add coordinates excluding the atoms to avoid
    for atom in res.get_atoms():
        if atom.get_id() not in atoms_to_avoid:  
            coords.append(atom.get_coord())    
    # Turn coordinate list in 3 by n numpy matrix
    coords = matrix(coords) 
    coords = coords.T
    # Calculate side chain center of mass
    center = coords.sum(1)/coords.shape[1]
    # Calculate centered side chain coordinates
    centered_coords = coords - center   
    # Return centered coordinates
    return (centered_coords)

def calc_RMSD(x, y):
    """
    Superimpose two residues, represented by two 3 by
    n numpy matrices, and return their minimum RMSD.
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
    # Calculate RMSD from the coordinates 
    d = x-u*y
    d = d.A*d.A
    rmsd = sqrt(sum(d)/N)                       
    return(rmsd)

def extract_RMSD_list(aa, list_structures):
    """
    Given a list of structures and an amino acid name, extract 1000 pairs
    of the selected amino acids, randomly and from two different proteins.
    """
    print("Computing %s RMSD scores..." % aa)
    n = 0
    rmsd_list = []
    # Repeat the while loop until 1000 residue pairs with equal number of atoms are selected
    while n < 1000:
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
        # Check if both proteins have the selected amino acid
        if len(list_res1) != 0 and len(list_res2) != 0:
            # Select randomly one amino acid from each of the two proteins
            i1 = random.choice(len(list_res1))
            i2 = random.choice(len(list_res2))
            res1 = list_res1[i1]
            res2 = list_res2[i2]
            # Get the centered coordinates of the two side_chains
            sc1 = get_centered_sidechain(res1)
            sc2 = get_centered_sidechain(res2)
            # Check if the side chains have the same number of atoms
            if sc1.shape == sc2.shape:
                # Calculate RMSD and append it to the RMSD list
                rmsd = calc_RMSD(sc1, sc2)
                rmsd_list.append(rmsd)
                # Add 1 to the counter, the loop stop when 1000 pairs are extracted
                n += 1 
            # Discharge the pair if the number of atoms are not the same
            else:
                pass
        # Repeat the cycle if the selected amino acid is not present in both proteins
        else:
            pass
    return (rmsd_list)

def save_histogram(rmsd_list, aa):
    """
    Make histogram and save as png file.
    """
    # The histogram of the data
    plt.hist(rmsd_list, 50, density=True)
    # Labels
    plt.xlabel('RMSD')
    plt.ylabel('Probability')
    plt.title('Histogram of RMSD for %s' % aa)
    # x and y limits of plot
    plt.xlim(0, 1.4)
    plt.ylim(0, 40)
    plt.grid(True)
    # Save
    plt.savefig("rmsd_hist_"+aa+".png")
    # Clear canvas for next plot
    plt.clf()

# Create a list containing the amino acids name that we want to investigate
aa_names = ["ARG","ASN","ASP","CYS","GLN","GLU","HIS","ILE","LEU",
            "LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]
# Extract the protein structures from PDB files in a given directory
list_structures = extract_structures("top500H")
# Create a dictionary in which each entry contains a list of 1000 RMSD from pairs of a specific amino acid
rmsd_dict = {}
for aa in aa_names:
    rmsd_dict[aa] = extract_RMSD_list(aa, list_structures)
# Create and save the plot of the RMSD scores distribution for each amino acid
for element in rmsd_dict:
    save_histogram(rmsd_dict[element], element)