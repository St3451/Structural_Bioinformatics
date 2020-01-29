# Implementation of the Nussinov Algorithm

The Nussinov algorithm for simple RNA folding by computation of the maximum number of base pairs employ a dynamic progrmaming algorithm. Here, this algorithm will be implemented.

## Background

1. Before imlementing the Nussinov algorithm work a scheme on how you iterate your way through the matrix for the version without bifurcation.
2. Now, implement a simple version computing the maximum number of base pairs still ignoring the bifurcation term.
3. Implement backtrack.
4. Repeat the steps above now including bifurcation.

## Implementation

This is a non-weighted implementation of the Nussinov algorithm, the final version (weighted score function that can consider specific constraint) can be found in the `final_exam` folder.

```python
#### My implementation of the Nussinov algorithm is based on Giulia Corsi (PhD fellow Animal Genetics, 
#### Bioinformatics and Breeding; University of Copenhagen) implementation and on line guides I received
#### in week 3 lecture exercises of Structural Bioinformatics. 
```

* __Create a list of list (matrix) to store our scores__
```python
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
```

# __Function that print a matrix__
```python
def print_matrix(matrix:list):
    """
    Take a nested list and print it as a matrix 
    (only works with small sequenceses)
    """
    # If the sequence is longer than 42 letters, than the matrix start to overlap. Anyway it can be helpful for debugging.
    if len(matrix) < 43:
        for row in matrix:
            print(row)
```

##### Score cell function
```python
def score_cell(seq, score_matrix, i, j, minimum_loop_size = 2):
    """
    Take a DNA seq, a constraint sequence, a score_matrix and a min_loop_size and 
    calculate the score of a cell(i,j) of the matrix
    """    
    score_list = []
    bp = 0
    # Check if there is a base pair in the constraint sequence, in that case force to have base pair    
    if (seq[i] == "A" and seq[j] == "U") or (seq[i] == "U" and seq[j] == "A") or (seq[i] == "G" and seq[j] == "C") or (seq[i] == "C" and seq[j] == "G") or (seq[i] == "G" and seq[j] == "U") or (seq[i] == "U" and seq[j] == "G"):
        bp += 1
        # Append the score from the diagonal (match)
        score_list.append(score_matrix[i+1][j-1] + bp)   
    # Append score from the left
    score_list.append(score_matrix[i][j-1])
    # Append score from the bottom
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
```

##### Iterate through the diagonals              
```python
def dinamic_programming_folding(seq, matrix_lst, min_loop_size):
    """
    Iterate through the diagonals of the 
    matrix and assign a score to each cell
    """
    for n in range(min_loop_size + 1, len(seq)):                      
        for j in range(n, len(seq)):                             
            i = j - n                                             
            score_cell(seq, matrix_lst, i, j, min_loop_size)                                 
```

##### Create a list for the dot bracket notation sequence
```python
def create_lst_db(seq):
    """
    Take a sequence and initialize a list to 
    store dot_bracket notation symbols
    """
    lst_db_ini = []
    for letter in seq:
        lst_db_ini.append(".")
    return lst_db_ini
```

##### Backtracking
```python
def db_build(row, col, dot_brackets_lst, min_loop_size, seq, score_matrix):   # the order of how we decide the conditions will determine which structure will be generated (same base pairs)
    """
    Take starting row, starting col, a list_dot_brackets, min_loop_size, 
    a DNA seq and its score_matrix. Update the list of dot_brackets by 
    traceback using base pair score
    """
    # Add stop criteria  
    while row < col + min_loop_size and score_matrix[row][col]!=0:           
        # Check if the score of the cell is obtained from the diagonal (match)                                                                    
        if score_matrix[row][col] == score_matrix[row+1][col-1] + 1 and (seq[row] == "A" and seq[col] == "U") or (seq[row] == "U" and seq[col] == "A") or (seq[row] == "G" and seq[col] == "C") or (seq[row] == "C" and seq[col] == "G") or (seq[row] == "G" and seq[col] == "U") or (seq[row] == "U" and seq[col] == "G"):  
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
        # Check for possible bifurcations                                                    
        else:
            for k in range(row, col - min_loop_size):
                if score_matrix[row][col] == score_matrix[row][k] + score_matrix[k+1][col]:
                    db_build(row, k, dot_brackets_lst, min_loop_size, seq, score_matrix)     # recursion (i,k)
                    db_build(k+1, col, dot_brackets_lst, min_loop_size, seq, score_matrix)   #           (k+1,j)
                    break
            else:
                print("\n ERROR, not match with anything in cell: (" + str(row) + ", " + str(col) + ")")
            break
```

##### Execution
```python
# Set the sequence
sequence = "AAACUUUCCCAGGG"  
# Set the minimum loop size                                           
m_l_size = 1                                                           
# Generate the initialized matrix
lista = create_matrix(sequence) 
# Iterate through the diagonals of the matrix and score the cells                                        
dinamic_programming_folding(sequence, lista, m_l_size)                  
# Print the updated matrix
print_matrix(matrix=lista)
# Generate an initialized list of dot brackets simbols                                                                             
lst_db  = create_lst_db(sequence)                                       
# Backtracking function to update dot brackets notation list
db_build(0, (len(sequence) - 1), lst_db, m_l_size, sequence, lista)     
# Print minimum loop size
print("\nMin loop size: " + str(m_l_size))
# Print sequence
print("Seq: " + sequence)
# Print the dot brackets as a string
print("Sdb: " + "".join(lst_db))                                        
```