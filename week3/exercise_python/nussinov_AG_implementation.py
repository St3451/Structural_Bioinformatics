#! /usr/bin/env python3

# lib for parsing fasta files
from Bio import SeqIO
# helper to make command line tools
import argparse
# an excellent array library
# SCIentific PYthon
import scipy as sp


def is_complementary(seq, a, b):
    """Check for base-pairing. Returns boolean as number."""
    # list of valid base pairings
    valid = {'AU', 'CG', 'GU'}
    # alphabetically sorted input
    pair = [seq[a], seq[b]]
    pair.sort()
    pair = ''.join(pair)
    # check
    return int(pair in valid)

assert is_complementary('ACGU', 0, 3) == 1
assert is_complementary('ACGU', 1, 3) == 0

def score(seq, mat, i, j):
    # compute simple cases
    extendLeft = mat[i + 1, j]
    extendRight = mat[i, j - 1]
    pair = mat[i + 1, j - 1] + is_complementary(seq, i, j)
    # get maximum for simple cases
    maxSimple = max(extendLeft, extendRight, pair)
    # check biforcations
    fork = []
    for k in range(i + 1, j - 1):
        fork.append(mat[i, k] + mat[k + 1, j])
    # if there was no biforcation, stick to the simple cases
    if len(fork) == 0:
        return maxSimple
    else:
        maxFork = max(fork)
        return max(maxSimple, maxFork)

def nussinov(seq):
    """Run the nussinov algorithm."""
    n = len(seq)
    # Create an empty matrix
    # this would resemble the 'java/c-like' init without preset values,
    # but a short negative value is easier to look at and to compare
    #mat = sp.empty((n, n), dtype=int)
    mat = -1 * sp.ones((n, n))
    # we could also just init everything with zeros
    #mat = sp.zeros((n, n))
    # we could also use as list of lists
    # then: indexing is slightly different
    #       mat[i, j] -> mat[i][j]
    #mat = [[-1 for i in range(n)] for j in range(n)]

    # Set diagonal to zero
    for i in range(n):
        mat[i, i] = 0
        if i > 0:
            mat[i, i - 1] = 0
    # Compute scores
    for offset in range(1, n):
        for j in range(offset, n):
            i = j - offset
            mat[i][j] = score(seq, mat, i, j)
    # Get the traceback
    # -1 because zero-based indexing
    dots = list('.' * n)
    visitNext = [(0, n - 1)]
    while len(visitNext) > 0:
        i, j = visitNext.pop()
        if i >= j: 
            # reached diagonal
            continue
        # Check if maximum came from a simple case in this order of
        if mat[i, j] ==  mat[i + 1, j]:
            # extend left
            visitNext.append((i + 1, j))
        elif mat[i, j] ==  mat[i, j - 1]:
            # extend right
            visitNext.append((i, j - 1))
        elif mat[i, j] ==  mat[i + 1, j - 1] + is_complementary(seq, i, j):
            # base pair
            dots[i] = '('
            dots[j] = ')'
            visitNext.append((i + 1, j - 1))
        else:
            # check biforcation
            for k in range(i + 1, j - 1):
                if mat[i, j] == mat[i, k] + mat[k + 1, j]:
                    visitNext.append((i, k))
                    visitNext.append((k + 1, j))
                    break
            else:
                raise RuntimeError("Did not find an origin for {}, {}".format(i, j))

    # convert dot list to string
    dotStr = ''.join(dots)
    return mat, dotStr


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=("Compute 2nd structure for "
        "the sequences in the input"))
    parser.add_argument('fasta', help='sequence file in fasta format',
        type=str)
    args = parser.parse_args()

    with open(args.fasta, 'r') as h:
        sequences = SeqIO.parse(h, 'fasta')
        for seq in sequences:
            assert(set(seq.seq) <= {'A', 'C', 'G', 'U'})
            mat, dots = nussinov(seq.seq)
            print(seq.description)
            print(mat)
            print(seq.seq)
            print(dots)
            print("Score: {}".format(mat[0, -1]))

