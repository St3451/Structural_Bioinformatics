from Bio import *
from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser
parser = PDBParser()
s = parser.get_structure("Trypsin", "/home/lpp/BIOINFORMATICS/sb2019/week5/2ptc.pdb")

### How to access data
for model in s:
    for chain in model:
        for residue in chain:
            for atom in residue:
                print (atom)

model = s[0]
chain_E = model['E']
chain_I = model['I']
residue_E = chain_E[100]
atoms_E = residue_E['CA']

atom = s[0]['E'][100]['CA']

for residue in model.get_residues():
    print (residue)
for residues in chain_E.get_residues():
    print(residues)
for atoms in chain_E.get_atoms():
    print(atoms)
for atom in s.get_atoms():
    print (atom)

chains_mzero = model.get_chains()
atoms_E = chain_E.get_atoms()
atoms_E_index = model['E'].get_atoms()
if atoms_E != atoms_E_index:
    print("They are not the same")

res_list = Selection.unfold_entities(s, 'R')                 # Get all residues from a structure
atom_list = Selection.unfold_entities(chain_E, 'A')          # Get all atoms from a chain
residue_list = Selection.unfold_entities(atom_list, 'R')     # etc
chain_list = Selection.unfold_entities(atom_list, 'C')


### How do I extract information from an Atom object and a residue object

a.get_name()       # atom name (spaces stripped, e.g. 'CA')
a.get_id()         # id (equals atom name)
a.get_coord()      # atomic coordinates
a.get_vector()     # atomic coordinates as Vector object
a.get_bfactor()    # isotropic B factor
a.get_occupancy()  # occupancy
a.get_altloc()     # alternative location specifier
a.get_sigatm()     # std. dev. of atomic parameters
a.get_siguij()     # std. dev. of anisotropic B factor
a.get_anisou()     # anisotropic B factor
a.get_fullname()   # atom name (with spaces, e.g. '.CA.')

r.get_resname()    # return the residue name (eg. 'GLY')
r.is_disordered()  # 1 if the residue has disordered atoms
r.get_segid()      # return the SEGID
r.has_id(name)     # test if a residue has a certain atom


### Measure

# Measure distance
ca1 = residue1['CA']
ca2 = residue2['CA']
distance = ca1 - ca2

# Measure angles
vector1 = atom1.get_vector()
vector2 = atom2.get_vector()
vector3 = atom3.get_vector()
angle = calc_angle(vector1, vector2, vector3)

# Measure torsion angles
vector1 = atom1.get_vector()
vector2 = atom2.get_vector()
vector3 = atom3.get_vector()
vector4 = atom4.get_vector()
angle = calc_dihedral(vector1, vector2, vector3, vector4)

### Polypeptide Class
ppb = PPBuilder()                              # Using C-N
for pp in ppb.build_peptides(s):
    print(pp.get_sequence())

ppb = CaPPBuilder()                            # Using CA-CA
for pp in ppb.build_peptides(s):
    print(pp.get_sequence())

### Get sequence
seq = polypeptide.get_sequence()

is_aa(residue)