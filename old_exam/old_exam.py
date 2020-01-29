from Bio import *
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser

from numpy import *
from numpy.linalg import *


######## 1)

def extract_protein_structures(directory):
    """ Take a folder as input and return a dictionary with entries for all the protein files in the folder. 
    Each entry contain the parsed structure of the protein"""
    p = PDBParser(QUIET = True) # Quiet ignore warnings 
    s_protein_dict = {}
    for filename in os.listdir(directory): # Add to a dictionary the structure of each protein file and the name of the protein as key
        try:                                                                
            s = p.get_structure(filename, os.path.join(directory, filename))  # create the structure
            s_protein_dict[filename] = s                                   # add the structure to the dictionary
        except:
            print(filename, " can't be parsed")
    return(s_protein_dict)

dict_structures = extract_protein_structures("/home/lpp/BIOINFORMATICS/sb2019/old_exam/top100H")

# I want a dictionary where each entry correspond to a list of segments each with 9 aa from a certain protein
def extract_protein_segments(n, dictionary_structures):
    ppbuilder = PPBuilder()
    dict_segments = {}
    for protein in dictionary_structures:                  # Iterate trough all stuctures extracted from the directory
            structure = dictionary_structures[protein]
            pp_list = ppbuilder.build_peptides(structure) # create a list of polypeptide for each structure
            list_segments = []
            for pp in pp_list:
                for i in range(len(pp)-(n-1)):    # ignore the last 8 amino acids
                    if is_aa(pp[i]) and is_aa(pp[i+1]) and is_aa(pp[i+2]) and is_aa(pp[i+3]) and is_aa(pp[i+4]) and is_aa(pp[i+5]) and is_aa(pp[i+6]) and is_aa(pp[i+7]) and is_aa(pp[i+8]):
                        sublist = pp[i:i+n]
                        list_segments.append(sublist)
            dict_segments[protein] = list_segments
    return(dict_segments)

dict_segments = extract_protein_segments(9, dict_structures)
#print(dict_segments["1aacH"][0][4]) # in this way I select the entry of the 1aacH protein, the first segment, the fifth residue


######## 2)

#### Part a) For Ala, Gly, Pro, Phe, Asp, Arg and Leu, select 500 random segment pairs that have that amino acid type as central residue
def exclude_keys(dictionary, keys):
    """Filters a dict by excluding certain keys."""
    key_set = set(dictionary.keys())
    key_set.remove(keys)
    return [key for key in key_set]


def select_500_gly(dictionary_segments):
    import random
    gly_list_of_pairs = []
    n = 0                      # n = number of pairs added to the list
    while n < 500:
        protein1 = random.choice(list(dictionary_segments.keys())) # select the first protein from which we can extract the segments
        tmp_keys = exclude_keys(dictionary_segments, protein1) # create a list of keys where we removed the key used for the first protein
        protein2 = random.choice(tmp_keys) # select the second protein

        len1 = len(dict_segments[protein1])        
        len2 = len(dict_segments[protein2])
        
        gly_pair = []
        debug = 0 # set the debug value = 0, if at the end of the for loop it didnt change, it means I added 1 pair 

        m = 0                                         
        debug_count1 = 0
        while m == 0:                                  # select a random segment until the fifth residue is a glycine (m=1)
            debug_count1 += 1
            index1 = random.randint(0, len1-1)         # the max index of a list is equal to len(list) -1
            segment1 = dict_segments[protein1][index1] 
                                                        
            if segment1[4].get_resname() == "GLY":      # check if the fifth residue is glycine
                gly_pair.append(segment1)               # append the first segment to the gly_pair list
                m = 1
            if debug_count1 == 5000:                    # if after 5000 try it can't find a res with glycine at fifth position, select an other protein
                debug = 1
                break

        m = 0
        debug_count2 = 0 
        while m == 0:                                   # do the same for the second segment of the pair
            debug_count2 += 1
            index2 = random.randint(0, len2-1)
            segment2 = dict_segments[protein2][index2]

            if segment2[4].get_resname() == "GLY":
                gly_pair.append(segment2)
                m = 1

            if debug_count2 == 5000:
                debug = 1
                break
        if debug == 0:                
            n += 1
            gly_list_of_pairs.append(gly_pair)          # append it to the glycine list of pairs 
            print(n, " Glycine pairs added")
                                                  # here I could procede in the same way adding other while loops for each amino acid list I want to create                  
    return(gly_list_of_pairs)

gly500 = select_500_gly(dict_segments)

#### Part b) Calculate RMSD

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

# Define a function that extract all atoms from a chain
def extract_atoms(chain):
    """Take a chain as input and return a matrix of atoms coordinates, where rows are 
    x,y,z axis and columns are the vectors of the atoms in the three dimensional space"""
    atoms = [[],[],[]]
    # Iterate through all the residues in the chain
    for residue in chain:
        # Check if it is a residue
        if is_aa(residue):
            # Get and store the coordinate of the C alpha  
            for atom in residue:                   
                atm = atom.get_vector()                      
                atm_x = atm[0]
                atm_y = atm[1]
                atm_z = atm[2]
                atoms[0].append(atm_x)
                atoms[1].append(atm_y)
                atoms[2].append(atm_z)
    # return the set of coordinates of the C alpha of the chain as 3 by n matrix, where n are the number of residues
    return (asmatrix(atoms))

# Define a function that superimpose two chains by applying a rotation matrix U that minimize the RMSD
def calc_RMSD(fragment_a, fragment_b):
    """Take two sets of vectors representing the coordinates of the alpha carbons and return a tuple
     of three elements containing the two RMSD calculated by two methods and the rotation matrix U"""

    x, y = extract_c_alpha(fragment_a), extract_c_alpha(fragment_b)

    # We want to center each fragments on the C-alpha atom of its central amino acid
    x_central_ca_x = extract_c_alpha(fragment_a)[0,4] # ca_x_coordinate
    x_central_ca_y = extract_c_alpha(fragment_a)[1,4] # ca_y_coordinate
    x_central_ca_z = extract_c_alpha(fragment_a)[2,4] # ca_z_coordinate
    x_central_ca = matrix([[x_central_ca_x], [x_central_ca_y], [x_central_ca_z]])

    y_central_ca_x = extract_c_alpha(fragment_b)[0,4]
    y_central_ca_y = extract_c_alpha(fragment_a)[1,4]
    y_central_ca_z = extract_c_alpha(fragment_a)[2,4]
    y_central_ca = matrix([[y_central_ca_x], [y_central_ca_y], [y_central_ca_z]])
    
    # Center
    x_centered = x - x_central_ca
    y_centered = y - y_central_ca

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





######## 3)

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.ticker import PercentFormatter

def calc_RMSD_pairs(pairs_list):
    """Take a list of pairs of protein segments and 
    calculate the RMSD between each pair"""
    list_RMSD = []
    for n in range(len(pairs_list)):
        RMSD = calc_RMSD(pairs_list[n][0], pairs_list[n][1])
        list_RMSD.append(RMSD[0])
    return(list_RMSD)

gly_rmsd_scores = calc_RMSD_pairs(gly500)
print(calc_RMSD_pairs(gly500))

x = array(gly_rmsd_scores)
plt.hist(x, bins= 20, rwidth=0.85, color='#0504aa')
plt.grid(axis='y', alpha=0.75)
plt.xlabel('RMSD');
plt.ylabel("Count")
plt.title("RMSD of 500 polypeptide segments with Glycine as central residue")
plt.show()

# try to calculate the RMSD by not all atoms but only N, CA, C, so I just need to change the extract atoms in extract ncac atoms