#### My implementation of the Nussinov algorithm is based on Giulia Corsi (PhD fellow Animal Genetics, 
#### Bioinformatics and Breeding; University of Copenhagen) implementation and on line guides I received
#### in week 3 lecture exercises of Structural Bioinformatics.

## 1) Implement the constraint folding in the Nussinov algorithm implementation

# Create a list(matrix) to store our scores
def create_matrix(seq):
    """
    Take a sequence and return 
    a matrix initialize with 0
    """
    lst = []
    for row in range(len(seq)):
        lst.append([])
        for col in range(len(seq)):
            lst[row].append(0)
    return lst

# Function that print a matrix
def print_matrix(matrix:list):
    """
    Take a nested list and print it as a matrix 
    (only works with small sequenceses)
    """
    # If the sequence is longer than 42 letters, than the matrix start to overlap. Anyway it can be helpful for debugging.
    if len(matrix) < 43:
        for row in matrix:
            print(row)

# Score cell function
def score_cell(seq, constraint, score_matrix, i, j, minimum_loop_size = 2):
    """
    Take a DNA seq, a score_matrix and a min_loop_size and 
    calculate the score of a cell(i,j) of the matrix
    """
    score_list = []
    bp_score = 0
    # Check if there is a base pair in the constraint sequence, in that case force to have base pair
    if (constraint[i] == "(" and constraint[j] == ")") or (constraint[j] == "(" and constraint[i] == ")"):
        if (seq[i] == "A" and seq[j] == "U") or (seq[i] == "U" and seq[j] == "A"):
            bp_score = 2
        elif (seq[i] == "G" and seq[j] == "C") or (seq[i] == "C" and seq[j] == "G"):
            bp_score = 3
        elif (seq[i] == "G" and seq[j] == "U") or (seq[i] == "U" and seq[j] == "G"):
            bp_score = 1
        # Append the score from the diagonal (match)
        score_list.append(score_matrix[i+1][j-1] + bp_score)   
    # Check if the base pair is forbidden in the constraint sequence, in that case force to avoid base pair
    elif (constraint[i] == "x"):
        # Append score from the left
        score_list.append(score_matrix[i][j-1])                 
        # Append score from the bottom
        score_list.append(score_matrix[i+1][j])                 
    # Check if there are not constraints, in that case allow both base pairing (diagonal) and not base pairing (left or bottom)
    elif (constraint[i] == "." and constraint[j] == "."):
        if (seq[i] == "A" and seq[j] == "U") or (seq[i] == "U" and seq[j] == "A"):
            bp_score = 2
        elif (seq[i] == "G" and seq[j] == "C") or (seq[i] == "C" and seq[j] == "G"):
            bp_score = 3
        elif (seq[i] == "G" and seq[j] == "U") or (seq[i] == "U" and seq[j] == "G"):
            bp_score = 1
        # Append the score from the diagonal (match)
        score_list.append(score_matrix[i+1][j-1] + bp_score)    
        # Append the score from the left
        score_list.append(score_matrix[i][j-1])                
        # Append the score from the bottom
        score_list.append(score_matrix[i+1][j])
    # Compute and append the bifurcation scores
    k_scores = []                                          
    for k in range(i,j - minimum_loop_size):
        score = score_matrix[i][k] + score_matrix[k+1][j]
        k_scores.append(score)
    # Append to the score list the max of the bifurcation scores
    score_list.append(max(k_scores))
    # Score the cell with the max of the four scores (diagonal + basepair score, left, bottom and bifurcation)
    score_matrix[i][j] = max(score_list)                                

# Create a list for the dot bracket notation sequence
def create_lst_db(seq):
    """
    Take a seq and initiate a list to 
    store dot_bracket notation symbols
    """
    lst_db_ini = []
    for letter in seq:
        lst_db_ini.append(".")
    # Return an initialized list of dots for the dot bracket notation sequence
    return lst_db_ini

# Iterate through the diagonals               
def dinamic_programming_folding(seq, con, matrix_lst, min_loop_size, constrain = False):
    if constrain == False:
        con = create_lst_db(seq)
    for n in range(min_loop_size + 1, len(seq)):                      
        for j in range(n, len(seq)):                             
            i = j - n                                             
            score_cell(seq, con, matrix_lst, i, j, min_loop_size)  

# Backtracking
def db_build(row, col, dot_brackets_lst, min_loop_size, seq, score_matrix):   # the order of how we decide the conditions will determine which structure generate (same pairs)
    """
    Take starting row, starting col, a list_dot_brackets, min_loop_size, 
    a DNA seq and its score_matrix. Update the list of dot_brackets by 
    traceback using base pair weighted score
    """
    # Add stop criteria  
    while row < col + min_loop_size and score_matrix[row][col]!=0:
        # Check if the score of the cell is obtained from the diagonal (match)
        if (score_matrix[row][col] == score_matrix[row+1][col-1] + 2 and ((seq[row] == "A" and seq[col] == "U") or (seq[row] == "U" and seq[col] == "A"))) or (score_matrix[row][col] == score_matrix[row+1][col-1] + 3 and ((seq[row] == "G" and seq[col] == "C") or (seq[row] == "C" and seq[col] == "G"))) or (score_matrix[row][col] == score_matrix[row+1][col-1] + 1 and ((seq[row] == "G" and seq[col] == "U") or (seq[row] == "U" and seq[col] == "G"))):  
            dot_brackets_lst[col] = ")"
            dot_brackets_lst[row] = "("
            col -= 1
            row += 1
        # Check if the score of the cell is obtained from the bottom
        elif score_matrix[row][col] == score_matrix[row+1][col]:              
            row += 1
        # Check if the score of the cell is obtained from the left
        elif score_matrix[row][col] == score_matrix[row][col-1]:             
            col -= 1                                                         
        else:
            # Check for possible bifurcations
            for k in range(row, col - min_loop_size):                         
                if score_matrix[row][col] == score_matrix[row][k] + score_matrix[k+1][col]:
                    db_build(row, k, dot_brackets_lst, min_loop_size, seq, score_matrix)
                    db_build(k+1, col, dot_brackets_lst, min_loop_size, seq, score_matrix)
                    break
            else:
                raise RuntimeError("\n ERRORE not match with anything in Cell: (" + str(row) + ", " + str(col) + ")")
            break

# Set the sequence
seq = "GGGGGUAUAGCUCAGGGGUAGAGCAUUUGACUGCAGAUCAAGAGGUCCCUGGUUCAAAUCCAGGUGCCCCCU"
# Set the constraint if needed
con = ".........................(((((xxxxxxx))))).............................."
# Set the minimum loop size
ml_size = 3                                              
# Generate the initialized matrix
matrix1 = create_matrix(seq)    
matrix2 = create_matrix(seq)                          
# Iterate through the diagonals of the matrix and score the cells
dinamic_programming_folding(seq, con, matrix1, ml_size, constrain = False)
dinamic_programming_folding(seq, con, matrix2, ml_size, constrain = True)   
# Generate an initialized list of dot brackets simbols                               
lst_db1 = create_lst_db(seq)     
lst_db2 = create_lst_db(seq)                                 
# Backtracking function to update dot brackets notation list
row = 0
col = len(seq) - 1
db_build(row, col, lst_db1, ml_size, seq, matrix1) 
db_build(row, col, lst_db2, ml_size, seq, matrix2) 
# Print minimum loop size
print("\nMin loop size: " + str(ml_size))
# Print sequence
print("Seq:            " + seq)
# Print the dot brackets as a string
normal_str = "".join(lst_db1)
constraint_str = "".join(lst_db2)
print("Sdb normal:     " + normal_str)   
print("Sdb constraint: " + constraint_str)

## 2) Calculate base pair distance between normal and constraint structure

# Base-pair distance
def calculate_bp_distance(seq1_db, seq2_db):
    """Take two dot bracket sequences and return their base pairs distance"""
    lst_opening_1 = []
    lst_couples_1 = []
    # Iterate through the symbols of the first dot bracket sequence
    for index in range(len(seq1_db)):
        # If I have an open "(" append the position of the "(" to the lst_opening
        if seq1_db[index] == "(":
            lst_opening_1.append(index)               
        # if I have a closed ")", remove (pop) the last element of lst_opening and create a bp_coordinate with:
        elif seq1_db[index] == ")":                   # - element removed from lst_opening (position opening bracket)
            elem = lst_opening_1.pop()                # - current position (position closing bracket)                               
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

# Execute 
print("BP distance: ", calculate_bp_distance(normal_str, constraint_str), "\n")