#!/usr/bin/python
# Copyright Stefan E Seemann <seemann@rth.dk>, Alexander Junge <ajunge@rth.dk>

# the sequence to be folded
seq = 'AAACUUUCCCAGGG'
#seq = 'CGUAAAGC'

# sequence length
l = len(seq)

# minimum loop size
h = 1

# basepair scores - contains allowed base pairs and allows to score base pairs differently, 
# e.g., according to their energy, although this is not done so far
scores = {'AU': 1, 'UA': 1, 'GU': 1, 'UG': 1, 'GC': 1, 'CG': 1}

# initialize matrix
m = [[0 for i in range(l)] for j in range(l)]

#fill scoring matrix
for j0 in range(h+1, l): # the diagonal of the matrix to loop over
	for i in range(0, l-j0): # the entry on the diagonal to fill
		j= i + j0
		#rule 1) i,j paired
		if seq[i] + seq[j] in scores:
			m[i][j] = m[i+1][j-1] + scores[seq[i] + seq[j]]
		#rule 2) i unpaired
		if m[i+1][j] > m[i][j]:
			m[i][j] = m[i+1][j]
		#rule 3) j unpaired 
		if m[i][j-1] > m[i][j]:
			m[i][j] = m[i][j-1]
		#rule 4) bifurcation k
		for k in range(i+1+h, j-1-h):
			if m[i][k] + m[k+1][j] > m[i][j]:
				m[i][j] = m[i][k] + m[k+1][j]
		#print i,j,seq[i],seq[j],m[i][j]
	
#print out original scoring matrix
print "RAW         ",
for i in range(0,l):
	print (" %2d " % i),
print ""
print "RAW         ",
for i in range(0,l):
	print ("  %s " % seq[i]),
print ""
for i in range(0,l):
	print "RAW  %2d  %s  " % (i, seq[i]),
	for j in range(0,l):
		print (" %2d " % m[i][j]),
	print ""
print ""


#backtracking - simple implementation of pseudo-code given in lecture
structure = ['.' for i in range(l)]
stack = []
stack.append((0,l-1))
while len(stack) > 0:
	top = stack.pop(),
	i = top[0][0]
	j = top[0][1]
	if i >= j:
		continue
	elif m[i+1][j] == m[i][j]:
		stack.append((i+1,j))
	elif m[i][j-1] == m[i][j]:
		stack.append((i,j-1))
	elif seq[i] + seq[j] in scores and m[i+1][j-1] + scores[seq[i]+seq[j]] == m[i][j]:
		#record basepair i,j
		#print "BP",i,j,":",seq[i],seq[j],"\n",
		structure[i] = "("
		structure[j] = ")"
		stack.append((i+1,j-1))
	else:
		for k in range(i+1+h, j-1-h):
			if m[i][k]+m[k+1][j] == m[i][j]:
				stack.append((k+1,j))
				stack.append((i,k))
				break

#output
print seq, l, "\n", ''.join(structure), "\n", "Score:", m[0][l-1]
