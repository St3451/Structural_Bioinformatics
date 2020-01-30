# Final exam

## A) Protein practical 

Your task is to examine the variability in terms of RMSD of the side chains of the 18 amino acids excluding Gly and Ala. Gly and Ala are excluded because they lack degrees of freedom in their side chains due to their small size.   
• As protein data base, use the top500 collection of high quality protein structures.  
• Use Bio.PDB to implement the script.  
    * Disregard any structures that cannot be parsed by the Bio.PDB parser, but ignore warnings.  
• For each of the 18 amino acids (Gly and Ala excluded), select 1000 pairs randomly sampled from the protein data set with replacement.  
    * The two amino acids in each pair should come from different proteins.   
• Superimpose the side chain atoms using the optimal RMSD algorithm.  
    * Side chain atoms are here defined as C-alpha, C-beta and anything attached beyond the C-beta. Main chain atoms (N, C and O) and hydrogens are excluded.
    * Make a well-justified decision on how you are going to center the atoms before applying the optimal RMSD superposition.  
        * Make a histogram of the RMSD distribution for each of the 18 amino acids. Make sure all histograms use the same scale on the x- and y-axis.  
    * Discuss and interpret the results.  

## 1. Introduction  
The task of this exercise is to investigate the distributions of RMSD scores between side chains pairs of 18 different amino acids. The two side chains from each pair should come from different proteins, this is a way to compare the variability of the amino acid structures and their differences. Gly and Ala are excluded from the analysis since their side chains are too small, and without enough degrees of freedom, to present a significative difference in terms of structural variability.  

## 2. Materials and methods  
For the exercise we used the Top500 database of PDB files, available from [Richardson Lab Web Site](http://kinemage.biochem.duke.edu/subindex.php#database). Richardson and colleagues, from Duke University, used this data for their Ramachandran and rotamer studies. This is a selection of 500 files from the Protein Data Bank (PDB) that are high resolution (1.8 Å or better), low homology, and high quality. The PDB format provides a standard representation for macromolecular structure data derived from X-ray diffraction and NMR studies. 

The programming language we used to perform the analysis is Python 3. In addition we used NumPy package to do operations with vectors and matrices, Matplotlib to plot the histograms and Bio.Python to work with the PDB files. In particular we used a Bio.Python module called Bio.PDB, the module has been developed by Thomas Hamelryck and focuses on working with crystal structures of biological macromolecules. It contains a parser for PDB files that makes the atomic information available in an easy-to-use but powerful data structure.

### 2.1 Implementation
For my implementation I used five functions that I will not completely report here to avoid redundancy. The first function extract the protein structures from a given directory. The second function extract the atoms coordinates of the side chain of a given residue, calculate the side chain center of mass and return the centered set of coordinates. The third function superimpose two side chains, represented by two 3 by N NumPy matrices, and return their minimum RMSD. The fourth function extract 1000 side chains pairs (of a certain amino acid) randomly sampled from the protein data set, and calculate their RMSD. The last function is used to make and save the plots of the RMSD distributions.

#### Parsing the structure
I start my implementation by parsing all the structure contained in the Top500H directory, try and except are used in order to ignore the structure that cannot be parsed by the PDBParser. I use the miscellaneous operating system interfaces (OS) module to access the PDB files contained in the directory.

```python
p = PDBParser(QUIET = True)
    list_structures = []
for filename in os.listdir(directory):
try:                                                            
    s = p.get_structure(filename, os.path.join(directory, filename))  
    list_structures.append(s)                                      
except:
    print(filename, "can't be parsed")
```

#### Extract side chain pairs
I select randomly two protein structures using random.choice() from NumPy package. Then I extract all the desired amino acids from each protein and, if both proteins contain the selected amino acid, I choose randomly one residue from each protein.

```python
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
```

Then I obtain the centered coordinates of their side chains (method described in the next subsection) and, if the two side chains have the same number of atoms, I compute the RMSD score.

```python
    # Get the centered coordinates of the two side_chains
    sc1 = get_centered_sidechain(res1)
    sc2 = get_centered_sidechain(res2)
    # Check if the side chains have the same number of atoms
    if sc1.shape == sc2.shape:
        # Calculate RMSD and append it to the RMSD list
        rmsd = calc_RMSD(sc1, sc2)
        rmsd_list.append(rmsd)
```
#### Optimal RMSD superposition
The implementation of the RMSD algorithm (and most of the rest of the code) is based on weekly exercises solutions provided by our Structural Bioinformatics professor Thomas Hamelryck. 
In order to measure the structural similarity between side chain pairs, we used the root-mean-square deviation (RMSD) of atomic positions which is simply the square root of the distance between all atoms divided by their number.
We want to apply a U rotation matrix to y, until the RMSD is minimized.



## B. RNA practical

The overall accuracy of RNA secondary structure prediction can be improved if for some stretches of the RNA sequence the structure is known. Knowledge about such substructures could come from either experimental or computational studies. Here, we will extend the Nussinov algorithm to consider a single known substructure at a given location. In addition, we require a minimum loop length of 3.

Use your implementation to predict the structure of the sequence:
```
Sequence: GGGGGUAUAGCUCAGGGGUAGAGCAUUUGACUGCAGAUCAAGAGGUCCCUGGUUCAAAUCCAGGUGCCCCCU
```
for which the following substructure is already known (from position 26 to 42):
```
Constraint: .........................(((((xxxxxxx)))))..............................
```