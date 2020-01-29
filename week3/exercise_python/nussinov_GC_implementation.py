#! /usr/bin/env python3
import argparse as ap
import pandas as pd
import numpy as np
from Bio import SeqIO

# The program takes in input a fasta file. For each sequence in the file, it computes the Nussinov Matrix and the backtracking.
# It includes bifurcations and minimum loop size, that can be set up with the variable "MIN_SIZE_LOOP"
# Copyright: 2018 Giulia I. Corsi giulia@rth.dk
#**********************How to run the program:********************#
# python3 nussinov_iterative_with_min_loop.py -i rna_seq.fasta
# or
# chmod +x nussinov_iterative_with_min_loop.py
# ./nussinov_iterative_with_min_loop.py -i rna_seq.fasta
# where rna_seq.fasta is the fasta file containing your sequence
#*****************************************************************#

MIN_SIZE_LOOP = 1

# Here we initialize the matrix
def initialize_matrix(RNA_seq: str):
    nucleotides_list = [nt for nt in RNA_seq]
    df = pd.DataFrame(index=nucleotides_list, columns=nucleotides_list,
                      data=np.triu([[-1] * len(RNA_seq)] * len(RNA_seq), +MIN_SIZE_LOOP + 1))
    return df

# This function returns 1 if the two nucleotides in "tuple nucleotides" can base pair, 0 otherwise
def base_pair_score(tuple_nucleotides: tuple):
    if 'A' in tuple_nucleotides and 'U' in tuple_nucleotides:
        return 1
    if 'G' in tuple_nucleotides and ('U' in tuple_nucleotides or 'C' in tuple_nucleotides):
        return 1
    else:
        return 0

# This function computes the nussinov score, with bifurcations, in one cell (whose position is defined by row and col)
def max_bp(df: pd.DataFrame, row: int, col: int):
    matrix = df.values
    max_bp_no_bifurcation = max(matrix[row][col - 1],
                                matrix[row + 1][col],
                                matrix[row + 1][col - 1] + base_pair_score(
                                    tuple_nucleotides=(df.columns[row], df.columns[col])))
    # List containing all the scores that can be obtained with bifurcations
    bifurcation_bp_scores = [matrix[row, fork_point] + matrix[fork_point + 1, col] for fork_point in
                             range(row + 1, col - 1)]
    if len(bifurcation_bp_scores) > 0:
        return max(max_bp_no_bifurcation, max(bifurcation_bp_scores))
    else:
        return max_bp_no_bifurcation

# Given a sequence of nucleotides, this function produces a matrix and updates its values by computing max base pair scores
def nussinov(seq: str):
    df = initialize_matrix(seq)
    print('*****Initialized matrix*****')
    print(df)
    len_seq = len(df)
    # Initialize matrix: upper triangle is set to -1 (which means the cell has not been scanned yet).
    # In case MIN_SIZE_LOOP is != 0 the upper triangle filled with -1 is moved up-right
    matrix = df.values
    # At each iteration we move to a different diagonal by increasing "offset"
    print('*****Filling the matrix*****')
    for offset in range(MIN_SIZE_LOOP + 1, len_seq):
        # At each iteration we move on the diagonal, scanning all its cells
        for col in range(offset, len_seq):
            row = col - offset
            # Here we are updating the cell value
            matrix[row][col] = max_bp(df=df, row=row, col=col)
    return df


def extract_dot_brackets_notation(row: int, col: int, dot_brackets_seq : list):
    matrix = df_seq_nussinov.values
    while row < col + MIN_SIZE_LOOP and matrix[row, col]!=0:
        # look down: can we come from there?
        if matrix[row, col] == matrix[row + 1, col]:
            row += 1
        # look left: can we come from there?
        elif matrix[row, col] == matrix[row, col - 1]:
            col -= 1
        # look down left: can we come from there?
        elif matrix[row, col] == matrix[row + 1, col - 1] + base_pair_score(
                tuple_nucleotides=(df_seq_nussinov.columns[row], df_seq_nussinov.columns[col])):
            dot_brackets_seq[row] = '('
            dot_brackets_seq[col] = ')'
            row += 1
            col -= 1
        else:
            # none of the previous cases was accepted..is it a bifurcation?
            for k in range(row + 1, col - 1):
                if matrix[row, col] == matrix[row, k] + matrix[k + 1, col]:
                    dot_brackets_seq = extract_dot_brackets_notation(row=row, col=k, dot_brackets_seq = dot_brackets_seq)
                    dot_brackets_seq = extract_dot_brackets_notation(row=k + 1, col=col, dot_brackets_seq=dot_brackets_seq)
                    break
            else:
                # it was not a bifurcation, and no match was found for the previous cases either.
                # Then there must be an error in the matrix,
                raise RuntimeError('There is an error in the nussinov algorithm implementation. '
                                   'Cell in position row=%i col=%i has a value that does '
                                   'not match any case in the Nussinov equation'
                                   % (row, col))
            break
    return dot_brackets_seq


# if __name__ == '__main__':
#     # ArgumentParser module is used to read inputs
#     parser = ap.ArgumentParser()
#     parser.add_argument('-i', '--input', help='Input RNA sequences in fasta format',
#                         type=str)
#     args = parser.parse_args()
#     # The SeqIO module of Biopython is used to read all sequences in the fasta file
#     with open(args.input, 'r') as RNA_seq:
#         RNA_sequences = SeqIO.parse(RNA_seq, 'fasta')
#         for RNA_seq in RNA_sequences:
#             print('Sequence name:\n%s' % RNA_seq.description)
#             sequence = str(RNA_seq.seq)
#             if 'T' in sequence:
#                 print('The DNA sequence is transformed to RNA')
#                 sequence = sequence.replace('T','U')
sequence = "GGGGGUAUAGCUCAGGGGUAGAGCAUUUGACUGCAGAUCAAGAGGUCCCUGGUUCAAAUCCAGGUGCCCCCU"
df_seq_nussinov = nussinov(seq=sequence)
print('Nussinov Matrix:')
print(df_seq_nussinov)
dot_brackets_seq = extract_dot_brackets_notation(row=0, col=len(RNA_seq) - 1, dot_brackets_seq = list('.' * len(RNA_seq)))
print('Sequence:\n%s' % RNA_seq.seq)
print('A possible structure:')
print(''.join(dot_brackets_seq),'\n')
