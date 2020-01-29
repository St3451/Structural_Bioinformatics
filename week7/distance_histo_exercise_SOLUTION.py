import Bio.PDB as bio 
import numpy as np
import matplotlib.pyplot as plt
"""
This scripts plots histograms of the distances of the C-alpha atoms of
various amino acids to the center of their proteins.
"""
def get_center(s):
    """
    Returns center of structure s as vector.
    """
    v=bio.Vector(0,0,0)
    n=0.0 # nr of atoms counter
    for atom in s.get_atoms():
        v=v+atom.get_vector()
        n+=1
    # Average
    c=v**(1/n)
    # Return center
    return c

def get_distances(s, aa_name, center):
    """
    For given center vector and structure s, return list of 
    distances for amino acid type aa_name.
    """
    rlist=[] # List of distances
    assert(bio.is_aa(aa_name)) # Make sure aa_name is amino acid
    for r in s.get_residues():
        if r.get_resname()==aa_name:
            try:
                ca=r["CA"]                     
                # Calculate distance
                v=ca.get_vector()
                diff=v-center
                r=diff.norm()
                rlist.append(r)
            except:
                # Missing C-alpha
                pass
    return rlist

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
    plt.savefig("dist_hist_"+aa+".png")
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

    # Loop over files in directory
    for fname in glob.glob("top500H/*"):
        print("Parsing "+fname+"...")
        p=bio.PDBParser(QUIET=True) # QUIET suppresses warnings
        try:
            # Get structure object s
            s=p.get_structure("", fname)
            # Get center for s
            center=get_center(s)
            # Get distance lists and append to global distance lists
            for aa in aa_names:
                rlists[aa]+=get_distances(s, aa, center)
        except:
            # Something went wrong pasring PDB file - skip
            print("\tERROR in "+fname)
            
    for aa in aa_names:
        # Print count, mean, standard deviation
        print("\n")
        l=len(rlists[aa])
        print("Count for %s: %i" % (aa, l))
        m=np.mean(rlists[aa])
        print("Mean for %s: %.2f" % (aa, m))
        sd=np.std(rlists[aa])
        print("Std. dev.  for %s: %.2f" % (aa, sd))
        # Save histogram
        save_histogram(rlists[aa], aa)