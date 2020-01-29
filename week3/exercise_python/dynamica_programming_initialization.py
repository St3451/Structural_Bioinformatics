#### Nussinov algo 1

s = "GUAUCCGA"
s3 = "AAACUUUCCCAGGG"
s2 = "UCUUACAGUGGCAUGUGACCGUUUAAGG"
s2= " AAACUUUCCC"


##Test2

# Create a list

def create_matrix(seq):
    lst = []
    for row in range(len(seq)):
        lst.append([])
        for col in range(len(seq)):
            lst[row].append(0)
    return lst

lista = create_matrix(s2)

# Function that print the matrix
def print_matrix(matrix:list):
    for row in matrix:
        print(row)

# Scor cell function
def score_cell(seq, minimum_loop_size = 2):
    score_list = []
    bp = 0
    print("\nCell: (" + str(i) + ", " + str(j) + ")")
    print("seq[i]: " + seq[i])
    print("seq[j] " + seq[j])
    if (seq[i] == "A" and seq[j] == "U") or (seq[i] == "U" and seq[j] == "A") or (seq[i] == "G" and seq[j] == "C") or (seq[i] == "C" and seq[j] == "G") or (seq[i] == "G" and seq[j] == "U") or (seq[i] == "U" and seq[j] == "G"):
        bp += 1
    score_list.append(lista[i+1][j-1] + bp)    # diagonal
    score_list.append(lista[i][j-1])           # from left
    score_list.append(lista[i+1][j])           # from bottom
    
    k_scores = []                              # bifurcation
    for k in range(i,j - minimum_loop_size):
        score = lista[i][k] + lista[k+1][j]
        k_scores.append(score)
    score_list.append(max(k_scores))

    lista[i][j] = max(score_list)              # append the max of the 4 options      
    print(k_scores)
    print(score_list)         
         


# Iterate through the diagonals
print(s2)
for n in range(2+1, len(s2)):           # 2 is min loop size, 1 is never access first diagonal                     # n should be the size, it moves the diagonal to j+1 after the first cycle
    for j in range(n, len(s2)):                                     # iterate over the diagonal 
        i = j - n
        score_cell(s2, 2)                 

print("\n")
print_matrix(matrix=lista)
print("\n")


# Backtracking
#lst_db = create_matrix(s)

def create_lst_db(seq):               # create list to store dot bracket notation symbols
    lst_db_ini = []
    for letter in seq:
        lst_db_ini.append(".")
    return lst_db_ini

lst_db  = create_lst_db(s2)

def db_build(row, col, dot_brackets_seq, min_size_loop, seq):
    print(row, col)
    print
    while row < col + min_size_loop and lista[row][col]!=0:
        print(lista[row][col])
        if lista[row][col] == lista[row+1][col]:                # look bottom
            row += 1
            print("down\n")
            print(row, col)
        elif lista[row][col] == lista[row][col-1]:              # look left
            col -= 1
            print("left\n")
            print(row, col)                                     # look diagonal 
        elif lista[row][col] == lista[row+1][col-1] + 1 and (seq[row] == "A" and seq[col] == "U") or (seq[row] == "U" and seq[col] == "A") or (seq[row] == "G" and seq[col] == "C") or (seq[row] == "C" and seq[col] == "G") or (seq[row] == "G" and seq[col] == "U") or (seq[row] == "U" and seq[col] == "G"):  
            print("***MATCH**** " + str(row) + " " + str(col))
            lst_db[col] = ")"
            lst_db[row] = "("
            col -= 1
            row += 1
            print("diagonal\n")
            print(row, col)
        else:
            for k in range(row, col - min_size_loop):
                if lista[row][col] == lista[row][k] + lista[k+1][col]:
                    db_build(row, k, dot_brackets_seq, min_size_loop)
                    db_build(k+1, col, dot_brackets_seq, min_size_loop)
                    break
                else:
                    print("\n ERRORE NOT MATCH WITH ANYTHING")
                    print(row, col)
                    break
                    #raise RuntimeError('Error:'
                     #                  'Cell in position row=%i col=%i has a value that does '
                     #                  'not match any case in the Nussinov equation'
                     #                   % (row, col))
            break






db_build(0, (len(s2) - 1), lst_db, 2, s2)                    # call the function(row, col, list_db, min_loop_size)
print("\n")
print("".join(lst_db))