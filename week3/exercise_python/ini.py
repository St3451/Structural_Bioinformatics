s = "GUAUCCGA"

##Test1

# Create the 8X8 matrix (list of the list)
table =[[],[],[],[],[],[],[],[]]
k = 0
for n in table:
    if k > len(s) - 1:
        break
    for i in range(0, 8):
        table[k].append(-100)
    k += 1

# Print the matrix
k = 0
for element in table:
    print(table[k])
    k += 1

print("\n"*2)



##Test2

# Create a dictionary
d_table = {}

# Iterate through the diagonals
for n in range(0, len(s)):                                         # n should be the size, it moves the diagonal to j+1 after the first cycle
    j = 0 + n
    i = 0
    for e in range(0, len(s)):                                     # iterate over the diagonal 
        bp = 0
        if j > len(s) - 1:
            break
        else:
            if n > 2:                                              # gives the score
                if (s[i] == "A" and s[j] == "U") or (s[i] == "U" and s[j] == "A") or (s[i] == "G" and s[j] == "C") or (s[i] == "C" and s[j] == "G") or (s[i] == "G" and s[j] == "U") or (s[i] == "U" and s[j] == "G"):
                    bp += 1
            
            d_table[i,j] = bp                                      # assign the bp score to the cell
            i += 1 
            j += 1



# Print the dictionary
g = 1
r = len(s)
for element in d_table:
    print(element)
    print(d_table[element])
    g += 1
    if g > r:
        print("\n")
        g = 1
        r = r - 1

