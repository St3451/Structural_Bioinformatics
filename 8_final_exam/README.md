# Final exam
The exam is composed of three parts: protein theory, protein practice and rna. The theory part (which can be found in `protein_theory.pdf`) is about an interesting comparison between two revolutionary protein structure prediction methods, DeepMind's AlphaFold and AlQuaraishi's end-to-end prediction. The code of the practical exercises can be found in `protein.py` and `rna.py`.

------------------------------------------------------------------------------------------------------------------------------------------

## A) Protein practical part

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

### 2.1. Implementation
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

In the exercise the centers of mass of the two sets of vectors used for the RMSD calculation are not at their origin, so I centered the atoms before applying the optimal RMSD superposition. Since the task was to compare the structural similarities between side chains of the same amino acid, I calculated the center of mass (COM) by adding all the coordinates of the side chain atoms to a vector, including the alpha carbon and excluding all hydrogen. Then I divided that vector for the number of atoms (N).  

Finally I centered the atoms by subtracting the center of mass to each coordinate vector in the set. Once I obtained the centered coordinates, I wanted to find the rotation matrix U that minimize the RMSD score, so I applied the singular value decomposition (SVD) to the correlation matrix R.  

Sometimes the rotation matrix U that minimize the RMSD is aroto-inversion, that will superimpose a mirror image. To avoid that we have to multiply the components of the rotation matrix U for Z = diag(1,1,−1), and we also change the sign to the third element of the diagonal matrix S (even that this is not necessary in the calculation of the RMSD from the coordinates). Than I applyed the rotation matrix U to y and I finally calculated the RMSD from the two set of coordinates.  

## 3. Results
The final results of the exercise are the plots of the RMSD scores distributions of the 18 different amino acids. It is possible to observe that, with some exceptions (ASN and ASP), they show three main distributions shapes.

![](pictures/rmsd_hist_CYS.png)
![](pictures/rmsd_hist_SER.png)
![](pictures/rmsd_hist_THR.png)
![](pictures/rmsd_hist_VAL.png)
![](pictures/rmsd_hist_HIS.png)
![](pictures/rmsd_hist_ILE.png)
![](pictures/rmsd_hist_LEU.png)
![](pictures/rmsd_hist_PHE.png)
![](pictures/rmsd_hist_PRO.png)
![](pictures/rmsd_hist_TRP.png)
![](pictures/rmsd_hist_TYR.png)
![](pictures/rmsd_hist_ARG.png)
![](pictures/rmsd_hist_GLN.png)
![](pictures/rmsd_hist_GLU.png)
![](pictures/rmsd_hist_LYS.png)
![](pictures/rmsd_hist_MET.png)
![](pictures/rmsd_hist_ASN.png)
![](pictures/rmsd_hist_ASP.png)

## 4. Conclusions
It is possible to observe that the distributions of RMSD scores exhibit certain patterns depending on the flexibility of the side chain, their size and other characteristics.  

Serine, Threonine, Valine and Cystein have the smallest side chain and they have a more compact structure. As expected they have the lowest RMSD score and exhibit the lowest variability in the score, meaning that they have the lowest ability to move in the three-dimensional space.  
![Image4.1](pictures/small_aa.png)

Tryptophan, Phenylalanine, Tyrosine, Isoleucine, Leucine, Proline and Histidine are the amino acids that present a characteristic pattern in their structural variability. In particular, they show a clear bimodal distribution of their RMSD scores, that may indicate that they can be found in two main different spatial arrangements. Probably the most reasonable explanation of this phenomenon is that, since most of these amino acids are present in β-sheets, the two peaks may be related to the fact that the side chains alternately point in opposite directions from the sheet.  
![Image4.2](pictures/medium_aa5.png)

Arginine, Lysine, Glutamate, Glutamine and Methionine are the amino acids that present the largest RMSD scores and have the largest variability in their distributions. With the exception of Methionine they are polar amino acids, they all have long side chains and, as it is possible to observe, due to the nature of their structure they have less constraints than others amino acids and they are therefore more flexible and have more freedom of movement.  
![Image4.3](pictures/large_aa2.png)

------------------------------------------------------------------------------------------------------------------------------------------

## B) RNA part

The overall accuracy of RNA secondary structure prediction can be improved if for some stretches of the RNA sequence the structure is known. Knowledge about such substructures could come from either experimental or computational studies. Here, we will extend the Nussinov algorithm to consider a single known substructure at a given location. In addition, we require a minimum loop length of 3.

1. Explain how this constraint can be implemented (hint: consider the initialization of the dynamic programming matrix). Then implement this constraint folding in your choice of Nussinov implementation (your own or one of those already available).

2. Use your implementation to predict the structure of the sequence:
```
Sequence: GGGGGUAUAGCUCAGGGGUAGAGCAUUUGACUGCAGAUCAAGAGGUCCCUGGUUCAAAUCCAGGUGCCCCCU
```
for which the following substructure is already known (from position 26 to 42):
```
Constraint: .........................(((((xxxxxxx)))))..............................
```

3. Predict the structure for the full sequence (without the folding constraint); again set the minimal loop length to 3. Then, compute the base pair distance between the two dot-bracket strings you have obtained (with/without folding constraint). Discuss the difference in the structure and provide a sketch of the two structures. Does one of the structures resembles a known structure?

4. Annotate the RNA sequence by a method of your choice. Which structure prediction (constraint/unconstraint) is compatible with the annotation?

5. Run the RNAfold webserver without and with constraint. Compare the foldings and describe your observations.

## 1. Introduction
The Nussinov algorithm is historically the first attempts at RNA secondary structure prediction. It is a dynamic programming algorithm that, given a RNA sequence, recursively finds the secondary structure that maximize the number of base pairs. In the Nussinov algorithm we first initialize a scoring matrix, we decide a minimum loop size and then we start filling each cell by iterating through the diagonals of the matrix. We can give a score to a cell 𝑖,𝑗 by looking at the three neighbor cells and the possible branching structures (bifurcation). We can give a score to a cell 𝑖,𝑗 by looking at the three neighbor cells and the possible branching structures (bifurcation). If we consider 𝐸(𝑖,𝑗) as the maximum number of base pairs for the subsequence 𝑥\[𝑖;𝑗\] , than 

![](pictures/nussi.png)

where  𝑠(𝑖,𝑗)  can be just  1  in case of base pair or can have different scores (3,2,1) depending on the specific bases that form a pair in that cell. After the matrix is filled, the optimal RNA secondary structure is obtained by backtracking from the top right cell.
The Nussinov algorithm is a really simplified version of RNA folding and it has several limitations. It has a computational complexity 𝑂(𝑁^3) and cannot consider pseudoknots, otherwise the computational complexity could reach 𝑂(𝑁^6). An other problem is the ambiguity because often the same structure can be procuded in several ways and the structures with the maximum number of base pairs are often not unique.  

Δ𝐺𝑙𝑜𝑜𝑝(𝑛)=Δ𝐺size (𝑛)+Δ𝐺sequence +Δ𝐺special,  

in this way this method provide a much better prediction and an unambiguous solution to the folding problem using free energy minimization [3].