#!/usr/bin/env python3

#************************************************How to run this script:************************************************
# python3 secondary_structure_SNP.py
# or
# chmod +x secondary_structure_SNP.py
# ./secondary_structure_SNP.py

#***********************************************************************************************************************

# RNA structures for exercise 1
RNA_EX1_STRUCT_WT = '.....(((((.(((..(((((((((....)))))))))..))))))))..'
RNA_EX1_STRUCT_MUT = '......(.((((((..(((((((((....)))))))))..)))))).)..'

# RNA structures for exercise 2
RNA_EX2_STRUCT_WT = '(((((((((((((............))))))........)))))))....'
RNA_EX2_STRUCT_MUT = '(((((((..((((((..........))))))........)))))))....'


#***********************************************************************************************************************

# Function that given two strings computes the Hamming distance between them.
# Hamming distance = number of positions in the data structures at which the elements differ.

def compute_hamming_distance(structure_wt: str, structure_mut: str):
    h_dist = 0
    for i in range(0, len(structure_wt)):
        if structure_wt[i] != structure_mut[i]:
            h_dist += 1
    return h_dist


# Function that given two dot-brackets strings representing RNA structures computes the base pair distance between them.
# Base pair distance = number of base pairs that are present in the first structure and not in the second one PLUS
# number of base pairs that are present in the second structure and not in the first one.
# NOTE: we use two stacks to store the dot-bracket structures. Brackets that are open at the same position on
# the two stacks and are closed at the same time in both stacks are conserved base pairs.
# If an open bracket is closed in one stack and not in the other at the same iteration, then that base pair is not
# conserved. Dots are included in the stack but skipped when a bracket is being closed (to make the bulge/loop..).
def compute_base_pair_distance(str_wt: str, str_mut: str):
    diff = 0
    stack_wt = []
    stack_mut = []
    for i in range(len(str_wt)):
        if str_wt[i] != ')':
            stack_wt.append(str_wt[i])
        else:
            elem_wt = stack_wt.pop()
            while elem_wt == '.':
                elem_wt = stack_wt.pop()
        if str_mut[i] != ')':
            stack_mut.append(str_mut[i])
        else:
            elem_mut = stack_mut.pop()
            while elem_mut == '.':
                elem_mut = stack_mut.pop()
        if str_mut[i] == ')' and len(stack_mut) != len(stack_wt):
            diff += 1
        if str_wt[i] == ')' and len(stack_mut) != len(stack_wt):
            diff += 1
    return diff


#***********************************************************************************************************************

hamming_distance_ex_1 = compute_hamming_distance(structure_wt=RNA_EX1_STRUCT_WT, structure_mut=RNA_EX1_STRUCT_MUT)
hamming_distance_ex_2 = compute_hamming_distance(structure_wt=RNA_EX2_STRUCT_WT, structure_mut=RNA_EX2_STRUCT_MUT)
base_pair_distance_ex1 = compute_base_pair_distance(str_wt=RNA_EX1_STRUCT_WT, str_mut=RNA_EX1_STRUCT_MUT)
base_pair_distance_ex2 = compute_base_pair_distance(str_wt=RNA_EX2_STRUCT_WT, str_mut=RNA_EX2_STRUCT_MUT)

print('##############################')
print('Results exercise 1:')
print('Hamming distance = %i' % hamming_distance_ex_1)
print('Base pair distance = %i' % base_pair_distance_ex1)
print('******************************')
print('Results exercise 2:')
print('Hamming distance = %i' % hamming_distance_ex_2)
print('Base pair distance = %i' % base_pair_distance_ex2)
print('##############################')

#***********************************************************************************************************************