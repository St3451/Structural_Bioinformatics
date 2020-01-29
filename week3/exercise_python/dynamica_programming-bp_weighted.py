#### Nussinov algo 3

s = "GUAUCCGA"
s2 = "AAACUUUCCCAGGG"
s4 = "UCUUAcAGUGGcaugugaCCGUUUAAGG"


s3 = "GUAUCCGA"
s2 = "AGCACACAGGC"
s2 = "GGGGGUAUAGCUCAGGGGUAGAGCAUUUGACUGCAGAUCAAGAGGUCCCUGGUUCAAAUCCAGGUGCCCCCU"
s2 = "GGGGGUAUAGCUCAGGGGUAGAGCA.................AGGUCCCUGGUUCAAAUCCAGGUGCCCCCU"
co = ".........................(((((xxxxxxx))))).............................."

















s4 = "GGGGGUAUAGCUCAGGGGUAGAGCAUUUGACUGCAGAUCAAGAGGUCCCUGGUUCAAAUCCAGGUGCCCCCU"
MY = "((((((.(.(((((....).)))))(((((.......))))))(..((((((.(...).))))).)))))))"
NI = "((((((.(.(((((....).)))))(((((.......))))))(..((((((.(...).))))).)))))))"
ii = "012345678901234567890123456789012345678901234567890123456789012345678901"












# Create a list(matrix) to store our scores
def create_matrix(seq):
    """Take a seq and return a matrix initialize with 0"""
    lst = []
    for row in range(len(seq)):
        lst.append([])
        for col in range(len(seq)):
            lst[row].append(0)
    return lst

# Function that print the matrix
def print_matrix(matrix:list):
    """Take a nested list and print it as a matrix"""
    for row in matrix:
        print(row)

# Score cell function
def score_cell(seq, score_matrix, i, j, minimum_loop_size = 2):
    """Take a DNA seq, a score_matrix and a min_loop_size and calculate the score of a cell(i,j) of the matrix"""
    score_list = []
    print("\nCell: (" + str(i) + ", " + str(j) + ")")
    print("BP: " + seq[i] + seq[j])
    bp_score = 0
    if (seq[i] == "A" and seq[j] == "U") or (seq[i] == "U" and seq[j] == "A"):
        bp_score = 2
    elif (seq[i] == "G" and seq[j] == "C") or (seq[i] == "C" and seq[j] == "G"):
        bp_score = 3
    elif (seq[i] == "G" and seq[j] == "U") or (seq[i] == "U" and seq[j] == "G"):
        bp_score = 1
    score_list.append(score_matrix[i+1][j-1] + bp_score)    # diagonal
    score_list.append(score_matrix[i][j-1])                 # from left
    score_list.append(score_matrix[i+1][j])                 # from bottom
    
    k_scores = []                                           # bifurcation
    for k in range(i,j - minimum_loop_size):
        score = score_matrix[i][k] + score_matrix[k+1][j]
        k_scores.append(score)
    score_list.append(max(k_scores))

    score_matrix[i][j] = max(score_list)              # score the cell with the max of the 4 options     
    print("MAX: " + str(max(score_list)))  
    print(k_scores)
    print(score_list)       
         


# Iterate through the diagonals               
print("Seq: " + s2)

def dinamic_programming_folding(seq, matrix_lst, min_loop_size):
    for n in range(m_l_size + 1, len(seq)):                      
        for j in range(n, len(seq)):                             
            i = j - n                                             
            score_cell(seq, matrix_lst, i, j, min_loop_size)                    


# Backtracking
def create_lst_db(seq):
    """Take a seq and initiate a list to store dot_bracket notation symbols"""
    lst_db_ini = []
    for letter in seq:
        lst_db_ini.append(".")
    return lst_db_ini

def db_build(row, col, dot_brackets_lst, min_loop_size, seq, score_matrix):   # the order of how we decide the conditions will determine which structure generate (same pairs)
    """Take starting row, starting col, a list_dot_brackets, min_loop_size, a DNA seq and its score_matrix.
    Update the list of dot_brackets by traceback using base pair weighted score"""
    print("\nStarting traceback:")
    while row < col + min_loop_size and score_matrix[row][col]!=0:
                                                                              # look diagonal (this time for different bp_score depending on the pair)
        if (score_matrix[row][col] == score_matrix[row+1][col-1] + 2 and ((seq[row] == "A" and seq[col] == "U") or (seq[row] == "U" and seq[col] == "A"))) or (score_matrix[row][col] == score_matrix[row+1][col-1] + 3 and ((seq[row] == "G" and seq[col] == "C") or (seq[row] == "C" and seq[col] == "G"))) or (score_matrix[row][col] == score_matrix[row+1][col-1] + 1 and ((seq[row] == "G" and seq[col] == "U") or (seq[row] == "U" and seq[col] == "G"))):  
            print("Cell: (" + str(row) + ", " + str(col) + "), generated from diagonal (****MATCH****)")
            dot_brackets_lst[col] = ")"
            dot_brackets_lst[row] = "("
            col -= 1
            row += 1
        elif score_matrix[row][col] == score_matrix[row+1][col]:              # look bottom
            print("Cell: (" + str(row) + ", " + str(col) + "), generated from bottom")
            row += 1
        elif score_matrix[row][col] == score_matrix[row][col-1]:              # look left
            print("Cell: (" + str(row) + ", " + str(col) + "), generated from left")
            col -= 1                                                         
        else:
            print(">>>>> Try bifurcation, enter else state <<<<<")
            for k in range(row, col - min_loop_size):                         # look possible bifurcations
                if score_matrix[row][col] == score_matrix[row][k] + score_matrix[k+1][col]:
                    print("Cell: (" + str(row) + ", " + str(col) + "), generated from Cells: (" + str(row) + ", " + str(k) + ") + (" + str(k + 1) + ", " + str(col) + ") (****BIFURCATION****)")
                    db_build(row, k, dot_brackets_lst, min_loop_size, seq, score_matrix)
                    db_build(k+1, col, dot_brackets_lst, min_loop_size, seq, score_matrix)
                    break
            else:
                raise RuntimeError("\n ERRORE NOT MATCH WITH ANYTHING in Cell: (" + str(row) + ", " + str(col) + ")")
            break


m_l_size = 3                                                # set min loop size
lista = create_matrix(s2)                                   # generate an initialized (with 0) matrix list
dinamic_programming_folding(s2, lista, m_l_size)            # iterate through the diagonals and score the cells
print("\n")
print_matrix(matrix=lista)                                  # print the updated matrix with scores
lst_db  = create_lst_db(s2)                                 # generate an initialized (with .) list of dot brackets symbols 
db_build(0, (len(s2) - 1), lst_db, m_l_size, s2, lista)     # backtracking function to update dot brackets list
print("\n")
print("Min loop size: " + str(m_l_size))
print("Seq: " + s2)
print("Sdb: " + "".join(lst_db))                            # convert the dot brackets list into a string  
lst_index =[]
for i in range(0, len(lst_db)):
    lst_index.append(str(i % 10))

print("Ind: " + "".join(lst_index))
