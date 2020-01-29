from Bio import *
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser

import matplotlib.pyplot as plt
import numpy as np

# Center (average of the cooridnates of all CA of the protein)

# We expect that the mean of PHE will be shorter than ASP because it is an hydrophilic aa

def calc_center(structure):
    """
    Returns center of structure s as vector.
    """
    ca_sum = Vector(0., 0., 0.)
    n = 0
    for residue in structure[0].get_residues():
        if is_aa(residue):
            try:
                # extract CA coordinate
                if residue.has_id("CA"):
                    ca_sum += residue["CA"].get_vector()
                    n += 1
            except:
                print(residue, " skipped, no CA atom")  
    center = ca_sum ** (1/n)  # cross product multiply a vector for a number
    return(center)

def calc_distance(structure, aa_name, center):
    """
    For given center vector and structure s, return list of 
    distances for amino acid type aa_name.
    """
    rlist = [] 
    assert(is_aa(aa_name)) # Make sure aa_name is amino acid
    for residue in structure[0].get_residues():
        if residue.get_resname() == aa_name:
            try:
                # Calculate distance
                ca = residue["CA"]
                vector = ca.get_vector()
                diff = vector - center
                r = diff.norm()
                rlist.append(r)
            except:
                pass
    return(rlist)

def save_histogram(rlist, aa):
    """
    Make histogram and save as png file.
    """
    # The histogram of the data
    n, bins, patches = plt.hist(rlist, 50, density=True, 
        facecolor='g', alpha=0.75)
    # Labels
    plt.xlabel('Distance')
    plt.ylabel('Probability')
    plt.title('Histogram of distances to center for %s' % aa)
    # x and y limits of plot
    plt.xlim(0, 40)
    plt.ylim(0, 0.08)
    plt.grid(True)
    # Save
    plt.savefig("/home/lpp/BIOINFORMATICS/sb2019/week7/dist_hist_"+aa+".png")
    # Clear canvas for next plot
    plt.clf()

# what is after that is not executed if I import this script as module
if __name__=="__main__":  
    import glob

    # Dictionary of global distance lists for given amino acids
    aa_names=["ALA", "GLY", "PRO", "HIS", "ASP", "PHE"]
    rlists={}
    for aa in aa_names:
        rlists[aa]=[]

    p = PDBParser(QUIET = True) # Quiet ignore warnings 
    s_protein_list = []
    directory = "/home/lpp/BIOINFORMATICS/sb2019/week7/top500H"
    for filename in os.listdir(directory): # Add to a dictionary the structure of each protein file and the name of the protein as key
        print("parsing ", filename, "...")
        try:
            # Get structure                                                                
            s = p.get_structure(filename, os.path.join(directory, filename))  # create the structure
            # Calculate center
            center = calc_center(s)
            # Calculate distance lists and append to global distance lists
            calc_distance(s, "ALA", center)
            for aa in aa_names:
                rlists[aa] += calc_distance(s, aa, center)
        except:
            print(filename, " can't be parsed")
    
    # At this point I have a dictionary with a list of distances for each amino acid
    for aa in aa_names:
        # Print count, mean, standard deviation
        print("\n")
        l = len(rlists[aa])
        print("Count for %s: %i" % (aa, l))
        m = np.mean(rlists[aa])
        print("Mean for %s: %.2f" % (aa, m))
        sd = np.std(rlists[aa])
        print("Std. dev.  for %s: %.2f" % (aa, sd))
        # Save histogram
        save_histogram(rlists[aa], aa)