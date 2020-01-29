# Parse PDB file
from Bio import *
from Bio.PDB import *
import Bio.PDB as bio
from Bio.PDB.PDBParser import PDBParser
import numpy as np
from numpy import *
from numpy.linalg import svd, det
import matplotlib.pyplot as plt
import glob
import random

# Loop over files in directory
parsed_list = []
for fname in glob.glob("top500H/*"):
    print("Parsing "+fname+"...")
    p = bio.PDBParser(QUIET=True) # QUIET suppresses warnings
    try:
        # Get structure object s and append it to a list
        s = p.get_structure("", fname)
        parsed_list.append(s)
    except:
        # Something went wrong parsing PDB file - skip
        print("\tERROR in "+fname) 

def center(matrix):
    '''calculate protein center of mass and centers the atoms'''
   # Calculate center of mass of protein side chain
    center_of_mass_coords = matrix.sum(1) / matrix.shape[1]
   # Center the matrix with the atoms coordinates
    centered_coords = matrix - center_of_mass_coords
    return centered_coords

def sup(sc1, sc2):
    '''superimpose the 2 side chains , given the centered matrix, and return the rmsd score'''
    # Nr of atoms
    N = sc1.shape[1]
    # Center side_chain_1 and side_chain_2
    sc1 = center(sc1)
    sc2 = center(sc2)
    # correlation matrix
    r = sc2 * sc1.T 
    # SVD of correlation matrix
    v, s, wt = svd(r)
    w = wt.T
    vt = v.T
    # Rotation matrix
    u = w * vt
    # Check for roto-reflection
    if det(u) < 0:
        z = diag((1, 1, -1))
        u = w * z * vt
        s[2] = -s[2]
    # Calculate RMSD from the coordinates 
    d = sc1 - u * sc2
    d = d.A*d.A
    rmsd = sqrt(sum(d)/N) 
    return(rmsd)

def get_coordinates(residue):
    '''get coordinates of the allowed atoms from a residue'''
    coordinates = []
    # atoms to exclude
    excluded_atoms = ['C', 'N', 'O', 'H', '1HB', '2HB', '1HG', '2HG', '1HD', '2HD', '1HH1', '2HH1', '1HH2', '2HH2', 'HE', '2HD2', 'HA', 'HG', '1HE2', '2HE2', 'HD2', 'HE1', 'HE2', '1HG1', '2HG1', '1HG2', '2HG2', '3HG2', '1HD1', '2HD1', '3HD1', '1HD2', '3HD2', '1HE', '2HE', '1HZ', '2HZ', '3HZ', '3HE', 'HD1', 'HZ', 'HB', 'HG1', 'HE3', 'HZ2', 'HZ3', 'HH2', 'HH', 'OH', '3HG1']
    # iterate through all the atoms in the residue
    for atom in residue.get_atoms():
        id_a = atom.get_id()
        # select the correct atoms and get their coordinates
        if id_a not in excluded_atoms:
            c = atom.get_coord() 
            coordinates.append(c)
            # turn list of coordinates into 3 x n numpy matrix
            coordinates_matr = matrix(coordinates)
            coordinates_matr = coordinates_matr.T
    return(coordinates_matr)

def save_histogram(distr, aa):
    '''make histogram, one for each amino acid, reporting the distribution of its rmsd and save the plot as png file'''
    # The histogram of the data
    plt.hist(distr, 50, density=True, color = 'r')

    # Labels
    plt.xlabel('RMSD')
    plt.ylabel('Probability')
    plt.title('Histogram of RMSD for %s side chain' % aa)
    # x and y limits of plot
    plt.xlim(0, 1.5)
    plt.ylim(0, 40)
    plt.grid(True)
    # Save
    plt.savefig("rmsd_hist_"+aa+".png")
    # Clear canvas for next plot
    plt.clf()

# list the 18 amino acids of interest
amino_acids = ['ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
for res_name in amino_acids:
    print("\nGetting RMSD distribution of", res_name, "...")
    # create a list to store the indexes of the proteins where the desired amino acid hasn't been found after 2000 tentatives
    no_residue_index = []
    count_0 = 0
    rmsd_distr = []
    # continue the cycle until we get 1000 pairs of the same residue from different proteins
    while count_0 < 1000:
        count_ch = 0
        # choose 2 proteins randomly from the list of parsed structures
        while count_ch == 0:
            ch1 = np.random.choice(len(parsed_list), size = 1)
            ch2 = np.random.choice(len(parsed_list), size = 1)
            # check that the 2 chosen proteins are not the same and that the indexes of the chosen proteins don't correspond to those stored in the list no_residue_index
            if ch1 != ch2 and (ch1 not in no_residue_index and ch2 not in no_residue_index):
                count_ch += 1
        structure1 = parsed_list[int(ch1)]
        structure2 = parsed_list[int(ch2)]
        list_res1 = list(structure1.get_residues())       # list of residues in protein 1
        list_res2 = list(structure2.get_residues())       # list of residues in protein 2

        # try to find the same residue in 2 different proteins, superimpose the side chains, return its rmsd score
        # append the scores into a list 
        try:
            count_1 = 0
            # when debug1 reach 2000, the desired amino acid has not been found in the protein
            debug1 = 0
            # continue until I get the first residue
            while count_1 == 0:
                debug1 += 1
                # if after 2000 attempt I don't get the desired residue
                if debug1 == 2000:
                    # append the index of the protein that does not contain the desired amino acid in no_residue_index
                    no_residue_index.append(ch1)
                    # cause an error that will break the try block and will restard the while cycle
                    discard == True
                # choose randomly one residue from the list
                v1 = int(np.random.choice(len(list_res1), size = 1))
                res1 = list_res1[v1]
                # if the chosen residue is the same as the one in the amino acid list, it will get its coordinates
                if res1.get_resname() == res_name:
                    c_r_1 = get_coordinates(res1)
                    count_1 += 1
            # repeat the same procedure for the second residue
            count_2 = 0
            debug2 = 0
            # continue until we get the second residue
            while count_2 == 0:
                debug2 += 1
                if debug2 == 2000:
                    no_residue_index.append(ch2)
                    discard == True 
                v2 = int(np.random.choice(len(list_res2), 1))
                res2 = list_res2[v2]
                if res2.get_resname() == res_name:
                    c_r_2 = get_coordinates(res2)
                    count_2 += 1 
            # check if the number of atoms in the 2 residues is the same
            if c_r_1.shape == c_r_2.shape:
                # superimpose the two residues
                sp = sup(c_r_1, c_r_2)
                # store the obtained score into a list
                rmsd_distr.append(sp)
                count_0 = count_0 + 1
                # report the number of scores obtained (once every 100 scores)
                if count_0 % 100 == 0:
                    print("Obtained", count_0, "scores")
        # if no amino acid is found, the residue is just saved in the variable and the while cycle restarts                      
        except:
            not_print = (res_name, "not found")
    # check how many proteins for each amino acid are analyzed before getting 1000 pairs       
    print(len(no_residue_index), "proteins marked without any or with few", res_name)
    # save the histogram
    save_histogram(rmsd_distr, res_name)