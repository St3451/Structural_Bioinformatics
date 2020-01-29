from Bio import *
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser

from numpy import *
from numpy.linalg import *

p = PDBParser()
#s = p.get_structure("1aacH", "/home/lpp/BIOINFORMATICS/sb2019/old_exam/top100H/1aacH")

## Get a list of all protein files in the directory
import glob
proteins_list = glob.glob("/home/lpp/BIOINFORMATICS/sb2019/old_exam/top100H/*")

## Use reg ex to get the structure of each protein
import re
pattern = re.compile("(\/home\/lpp\/BIOINFORMATICS\/sb2019\/old_exam\/top100H\/)(.*)") 

# Example
test = ("/home/lpp/BIOINFORMATICS/sb2019/old_exam/top100H/1aacH")
match = pattern.search(test)
print(match.group())                                                                    # entire matching string                                           
print(match.group(2))                                                                   # second group
s2 = p.get_structure(match.group(2), match.group())

# Extract with reg ex
def extract_protein_structures(protein_files):
    import glob
    import re
    p = PDBParser()
    # Get a list of all protein files in the directory
    proteins_list = glob.glob(protein_files)     
    s_protein_dict = {}
    # Set the pattern for the name of the protein (group(2))
    pattern = re.compile("(\/home\/lpp\/BIOINFORMATICS\/sb2019\/old_exam\/top100H\/)(.*)")  

    # Add to a dictionary the structure of each protein file and the name of the protein as key
    for file in proteins_list:
        match = pattern.search(file)
        if match:
            protein_name = match.group(2)
            protein_path = match.group()
            s = p.get_structure(protein_name, protein_path)
            s_protein_dict[protein_name] = s
    return(s_protein_dict)

dict_structures = extract_protein_structures("/home/lpp/BIOINFORMATICS/sb2019/old_exam/top100H/*")    
print(dict_structures)

# extract with reg ex
def extract_4protein_structures(protein_files_path):
    import glob
    import re
    p = PDBParser(QUIET = True)
    files_list = glob.glob(protein_files_path)     
    s_protein_dict = {}
    pattern = re.compile("(\/home\/lpp\/BIOINFORMATICS\/sb2019\/old_exam\/top100H\/)(.*)")  
  
    for file in files_list:
        match = pattern.search(file)
        if match:
            protein_name = match.group(2)
            protein_path = match.group()
            try:     
                s = p.get_structure(protein_name, protein_path)
                s_protein_dict[protein_name] = s
            except:
                print(protein_name, " can't be parsed")
    return(s_protein_dict)

# Extract with os
def extract_protein_structures(folder):
    p = PDBParser(QUIET = True) # Quiet ignore warnings   
    s_protein_dict = {}

    # Add to a dictionary the structure of each protein file and the name of the protein as key
    for filename in os.listdir(folder):
        try:                                                                
            s = p.get_structure(filename, os.path.join(folder, filename))
            s_protein_dict[filename] = s
        except:
            print(filename, " can't be parsed")
    return(s_protein_dict)

dict_structures = extract_protein_structures("/home/lpp/BIOINFORMATICS/sb2019/old_exam/top100H")  

# I want a dictionary where each entry correspond to a list of 9 aa fragments from a certain protein
ppbuilder = PPBuilder()
dict_fragments = {}
for protein in dict_structures:
    structure = dict_structures[protein]
    pp_list = ppbuilder.build_peptides(structure)      # create a list of polypeptides for each structure
    list_fragments = []                   # create a list to store all the 9aa fragments from one protein
    for pp in pp_list:
        for i in range(len(pp)-(9-1)):    # ignore the last 8 amino acids
            if is_aa(pp[i]) and is_aa(pp[i+1]) and is_aa(pp[i+2]) and is_aa(pp[i+3]) and is_aa(pp[i+4]) and is_aa(pp[i+5]) and is_aa(pp[i+6]) and is_aa(pp[i+7]) and is_aa(pp[i+8]):
                sublist = pp[i:i+9]
                list_fragments.append(sublist)
    dict_fragments[protein] = list_fragments



#### Center a chain

x = matrix([[ 1,  2,  3, 5],      # x
            [ 2,  4,  3, 2],           # y 
            [ 1,  5,  2, 1]], 'f')     # z

y = matrix([[ 3,  1,  4, 4],
            [ 4,  1,  2, 2],
            [ 3,  2,  2, 3]], 'f')

print(x.shape) # Check the shape of the popypeptide, should return a 3 by n matrix where n equal the number of residues


center_of_mass_x = x.sum(1) / x.shape[1]  # center of mass is a 1 by 3 matrix = [sum of x of all residues / n. of residues..],[y],[z] 
center_of_mass_y = y.sum(1) / y.shape[1]

print(x.shape[1])
print(x.sum(1))
print(center_of_mass_x)

x_centered = x - center_of_mass_x
y_centered = y - center_of_mass_y

print(x_centered)



def column(matrix, i):
    """Take a matrix and an index, extract the i'th column from the matrix"""
    return [row[i] for row in matrix]