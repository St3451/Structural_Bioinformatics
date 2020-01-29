
# Hamming distance
def calculate_hamming_distance(seq1_db, seq2_db):
    """Take two dot bracket sequences and return their hamming distance"""
    hamming_difference = 0
    for index in range(0, len(seq1_db)):
        if seq1_db[index] != seq2_db[index]:
            hamming_difference += 1
    return hamming_difference


# Base-pair distance

#WT
def calculate_bp_distance(seq1_db, seq2_db):
    """Take two dot bracket sequences and return their base pairs distance"""
    lst_opening_1 = []
    lst_couples_1 = []
    for index in range(len(seq1_db)):
        # If I have an open "(" append the position of the "(" to the lst_opening
        if seq1_db[index] == "(":
            lst_opening_1.append(index)               
        # if I have a closed ")", remove (pop) the last element of lst_opening and create a bp_coordinate with:
        elif seq1_db[index] == ")":                   # - element removed from lst_opening (position opening bracket)
            elem = lst_opening_1.pop()                # - current position (position closing bracket)                               
            pair = (elem, index)
            lst_couples_1.append(pair)

    lst_opening_2 = []
    lst_couples_2 = []
    for index in range(len(seq2_db)):
        if seq2_db[index] == "(":
            lst_opening_2.append(index)               
            
        elif seq2_db[index] == ")":                  
            elem = lst_opening_2.pop()               
                                                      
            pair = (elem, index)
            lst_couples_2.append(pair)

    set1 = set(lst_couples_1)                         # Check if there are differences between bp
    set2 = set(lst_couples_2)
    diff = len(set1.difference(set2)) + len(set2.difference(set1))
    return diff

## Set the sequences
seq1 = "(((((((.((((((....).))))..)).)).(((((((....))))(((((.(...).))))))))))))."
seq2 = "((((((.(.(((((....).)))))(((((.......))))))(..((((((.(...).))))).)))))))"

## Execute 
print("\nNormal:     " + seq1 + "\nConstraint: " + seq2)
print("Hamming distance: " + str(calculate_hamming_distance(seq1, seq2)))
print("BP distance: " + str(calculate_bp_distance(seq1, seq2)))