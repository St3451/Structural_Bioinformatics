# Final exam

## 1. Protein practical 

Your task is to examine the variability in terms of RMSD of the side chains of the 18 amino acids excluding Gly and Ala. Gly and Ala are excluded because they lack degrees of freedom in their side chains due to their small size.   
• As protein data base, use the top500 collection of high quality protein structures.  
• Use Bio.PDB to implement the script. 
     Disregard any structures that cannot be parsed by the Bio.PDB parser, but ignore warnings.  
• For each of the 18 amino acids (Gly and Ala excluded), select 1000 pairs randomly sampled from the protein data set with replacement.
     The two amino acids in each pair should come from dierent proteins. 
• Superimpose the side chain atoms using the optimal RMSD algorithm. 
     Side chain atoms are here dened as C-alpha, C-beta and anything attached beyond the C-beta. Main chain atoms (N, C and O) are excluded, exclude hydrogens. 
     Make a well-justified decision on how you are going to center the atoms before applying the optimal RMSD superposition. 
         Make a histogram of the RMSD distribution for each of the 18 amino acids. Make sure all histograms use the same scale on the x- and y-axis. 
     Discuss and interpret the results. 
```
```
## 2. RNA practical

The overall accuracy of RNA secondary structure prediction can be improved if for some stretches of the RNA sequence the structure is known. Knowledge about such substructures could come from either experimental or computational studies. Here, we will extend the Nussinov algorithm to consider a single known substructure at a given location. In addition, we require a minimum loop length of 3.

Use your implementation to predict the structure of the sequence:
```
Sequence:   GGGGGUAUAGCUCAGGGGUAGAGCAUUUGACUGCAGAUCAAGAGGUCCCUGGUUCAAAUCCAGGUGCCCCCU
```
for which the following substructure is already known (from position 26 to 42):
```
Constraint: .........................(((((xxxxxxx)))))..............................
```