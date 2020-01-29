# Assignment 3 - Measures to detect the eﬀect of SNPs on RNA secondary structure

The effect of SNPs on RNA secondary structure can be predicted by comparing the structures of wild-type and mutant (with SNP) RNA. The structures being compared may be either optimal (MFE) structure or ensemble structure (obtained from partition function).
While comparing the optimal structure between wild-type and mutant, the structures can be considered as two distinct strings and the diﬀerence at the pure string level be used to measure how divergent they are. This strategy is employed by the RNAmute program (Churkin and Barash, 2006):
1. Hamming distance - The number of position at which the corresponding symbols are diﬀerent.
2. base-pair distance - The total number of base pairs that are diﬀerent between two structures

## Exercise 1

Compute both the Hamming and base pair distance of the following pairs of sequence:  

1. __Pair 1__ 
```
WT   GCGGGCCCCGC ((((...)))) 
MUT  ACGGGCCCCGC .(((...))).
```
2. __Pair 2__
```
WT   CAAUCCCGGCUGCGUCCCAGUUGGAUUUAUCCAGCUGGUUCGUGCUGGUU .....(((((.(((..(((((((((....)))))))))..)))))))).. 
MUT  CAAUCCCGGCUGCGUCCCAGUUGGAUUUAUCCAGCUGGUUCGUGGUGGUU ......(.((((((..(((((((((....)))))))))..)))))).)..
```

3. __Pair 3__
```
WT   AGCGGGGGAGACAUAUAUCACAGCCUGUCUCGUGCCCGACCCCGCUGGUU .....(((((.(((..(((((((((....)))))))))..)))))))).. 
MUT  AGCGGGGGAGAGAUAUAUCACAGCCUGUCUCGUGCCCGACCCCGCUGGUU (((((((..((((((..........))))))........)))))))....
```

4. __Pair 3__ 
```python
WT   (((((..((((((((........(((((......)))))........)))))(((((...))))))))...)))))...((((((.((((((....)))))).).)))))..((((((...................))))))...((((((((((((.(((((((....))))))))))..((((((.....(((.((((((((.....))))))))....))).....))))))....))))))).))..
MUT  (((((..((((((((........(((((......)))))........)))))(((((...))))))))...)))))...((((((.((((((....)))))).).)))))..((((((...................))))))...(((((((((..(((((((..((((((...........))))))....))))))).....((((((....))))))...((......))......))))))).))..

```

### Implementation

* __Hamming distance__
```python
def calculate_hamming_distance(seq1_db, seq2_db):
    """
    Take two dot bracket sequences and 
    return their hamming distance
    """
    hamming_difference = 0
    for index in range(0, len(seq1_db)):
        if seq1_db[index] != seq2_db[index]:
            hamming_difference += 1
    return hamming_difference
```

* __Base-pair distance__
```python
def calculate_bp_distance(seq1_db, seq2_db):
    """
    Take two dot bracket sequences and 
    return their base pairs distance
    """
    lst_opening_1 = []
    lst_couples_1 = []
    # Iterate through the symbols of the first dot bracket sequence
    for index in range(len(seq1_db)):
        # If I have an open "(" append the position of the "(" to the lst_opening
        if seq1_db[index] == "(":
            lst_opening_1.append(index)               
        # if I have a closed ")", remove (pop) the last element of lst_opening and create a bp_coordinate with:
        elif seq1_db[index] == ")":                          # - element removed from lst_opening (position opening bracket)
            elem = lst_opening_1.pop()                       # - current position (position closing bracket)                              
            pair = (elem, index)
            lst_couples_1.append(pair)
    # Repeat for the second dot bracket sequence
    lst_opening_2 = []
    lst_couples_2 = []
    for index in range(len(seq2_db)):
        if seq2_db[index] == "(":
            lst_opening_2.append(index)                      
        elif seq2_db[index] == ")":                  
            elem = lst_opening_2.pop()                                                              
            pair = (elem, index)
            lst_couples_2.append(pair)
    # Check if there are differences between base pairs
    set1 = set(lst_couples_1)                         
    set2 = set(lst_couples_2)
    diff = len(set1.difference(set2)) + len(set2.difference(set1))
    return diff
```

* __Assign the pairs and print the output__
```python
wt1 = "((((...))))"
mut1= ".(((...)))."
wt2 = ".....(((((.(((..(((((((((....)))))))))..)))))))).."
mut2= "......(.((((((..(((((((((....)))))))))..)))))).).."
wt3 = "(((((((((((((............))))))........)))))))...."
mut3= "(((((((..((((((..........))))))........)))))))...."
wt4 = "(((((..((((((((........(((((......)))))........)))))(((((...))))))))...)))))...((((((.((((((....)))))).).)))))..((((((...................))))))...((((((((((((.(((((((....))))))))))..((((((.....(((.((((((((.....))))))))....))).....))))))....))))))).)).."
mut5= "(((((..((((((((........(((((......)))))........)))))(((((...))))))))...)))))...((((((.((((((....)))))).).)))))..((((((...................))))))...(((((((((..(((((((..((((((...........))))))....))))))).....((((((....))))))...((......))......))))))).)).."

print("\nPair 1")
print("Hamming difference: " + str(calculate_hamming_distance(wt1, mut1)))
print("BP distance: " + str(calculate_bp_distance(wt1, mut1)))
print("\nPair 2")
print("Hamming difference: " + str(calculate_hamming_distance(wt2, mut2)))
print("BP distance: " + str(calculate_bp_distance(wt2, mut2)))
print("\nPair 3")
print("Hamming difference: " + str(calculate_hamming_distance(wt3, mut3)))
print("BP distance: " + str(calculate_bp_distance(wt3, mut3)))
print("\nPair 4")
print("Hamming difference: " + str(calculate_hamming_distance(wt4, mut4)))
print("BP distance: " + str(calculate_bp_distance(wt4, mut4)))
```

## Exercise 2

Using the base pair probability information of the affected region in `dot_plot_matrix.txt` file (obtained from [RNAsnp webserver](https://rth.dk/resources/rnasnp/)), find the base pair in the wild type sequence with the highest probability and compare it to the probability of the same base pair position in the mutant version.

```python
def file_to_matrix(text_file):
    """
    Take a text tile and return 
    a list of list (matrix)
    """
    file_opened = open(text_file, "r")
    lst_matrix = []
    for line in file_opened.readlines():
        line = line.rstrip()
        lst_matrix.append(line.split())
    return lst_matrix

def find_max_top_matrix_triangle(matrice):
    """
    Take a matrix (list of list) as input, find the 
    maximum of its top right triangle and the value 
    of the same cell in the botton left triangle
    """
    WT_max = 0
    WT_max_cell = ()
    for row, elem in enumerate(matrice):           # row is the position and element is the row_list that contain the values
        for col, value in enumerate(elem):         # the position in the row_list are the columns. With enumerate() the first element is the index and the second is the element itself
            if row < col:                          # This condition ensure to be in the top right triangle of the matrix
                if float(value) > WT_max:
                    WT_max = float(value)
                    WT_max_cell = (row, col)
    MT_relat_value = float(lst_matrix[WT_max_cell[1]][WT_max_cell[0]])   # I just invert the coordinate (col, row) instead of (row, col)
    return WT_max_cell, WT_max, MT_relat_value

lst_matrix = file_to_matrix("dot_plot_matrix.txt")

WT_max_cell = find_max_top_matrix_triangle(lst_matrix)[0]
WT_max = find_max_top_matrix_triangle(lst_matrix)[1]
MT_relat_value = find_max_top_matrix_triangle(lst_matrix)[2] 

print("\nWT max cell (+1, +1): (" + str(WT_max_cell[0]) + ", " + str(WT_max_cell[1]) + ")")
print("WT max value: " + str(WT_max))
print("MT relative position value: " + str(MT_relat_value) + "\n") 
```

## Exercise 3 
In the local region with altered RNA structure find the number of base pairs with pair probabilities higher then 0.5 and 0.8. Do you see a consistent pattern?

```python
def find_cell_with_higher_prob(matrice):
    """
    Take a bp prob. matrix as an input and return 
    a list of all the values of the cells with a 
    probability higher than threeshold
    """
    wt_08plus = []
    wt_05plus = []
    mt_08plus = []
    mt_05plus = []

    for row, elem in enumerate(matrice):           # row is the position and element is the row_list that contain the values
        for col, value in enumerate(elem):         # the position in the row_list are the columns
            if row < col:                          # Check in wild type triangle (top right triangle)
                if float(value) > 0.5:
                    wt_05plus.append(float(value))
                    if float(value) > 0.8:
                        wt_08plus.append(float(value))
            if row > col:                          # Check in mutant triangle
                if float(value) > 0.5:
                    mt_05plus.append(float(value))
                    if float(value) > 0.8:
                        mt_08plus.append(float(value))
    return wt_08plus, wt_05plus, mt_08plus, mt_05plus

wt_08plus = find_cell_with_higher_prob(lst_matrix)[0]
wt_05plus = find_cell_with_higher_prob(lst_matrix)[1]
mt_08plus = find_cell_with_higher_prob(lst_matrix)[2]
mt_05plus = find_cell_with_higher_prob(lst_matrix)[3]

ratio_MT_WT08 = len(mt_08plus)/len(wt_08plus)
ratio_MT_WT05 = len(mt_05plus)/len(wt_05plus)

print("WT higher than 0.8: " + str(len(wt_08plus)))
print("WT higher than 0.5: " + str(len(wt_05plus)))
print("Mutant higher than 0.8: " + str(len(mt_08plus)))
print("Mutant higher than 0.5: " + str(len(mt_05plus)))
print("Ratio MT/WT bp with prob higher than 0.8: " + str(ratio_MT_WT08))
print("Ratio MT/WT bp with prob higher than 0.5: " + str(ratio_MT_WT05))
```