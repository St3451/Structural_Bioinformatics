#3.Ex 
# Using the base pair probability information of the affected region find the base pair in the wild type sequence with the
# highest probability and compare it to the probability of the same base pair position in the mutant version.

def file_to_matrix(text_file):
    """Take a text tile and return a list of list (matrix)"""
    file_opened = open(text_file, "r")
    lst_matrix = []
    for line in file_opened.readlines():
        line = line.rstrip()
        lst_matrix.append(line.split())
    return lst_matrix


lst_matrix = file_to_matrix("/home/lpp/BIOINFORMATICS/sb2019/week4/exercise_distance_SNP_RNA/dot_plot_matrix.txt")

def find_max_top_matrix_triangle(matrice):
    """Take a matrix (list of list) as input, find the maximum 
    of its top right triangle and the value of the same cell in the botton
    left triangle"""
    WT_max = 0
    WT_max_cell = ()
    for row, elem in enumerate(matrice):           # row is the position and element is the row_list that contain the values
        for col, value in enumerate(elem):         # the position in the row_list are the columns
            if row < col:
                if float(value) > WT_max:
                    WT_max = float(value)
                    WT_max_cell = (row, col)
    return WT_max_cell, WT_max, float(lst_matrix[WT_max_cell[1]][WT_max_cell[0]])

WT_max_cell = find_max_top_matrix_triangle(lst_matrix)[0]
WT_max = find_max_top_matrix_triangle(lst_matrix)[1]
MT_relat_value = find_max_top_matrix_triangle(lst_matrix)[2] 
print("\nWT max cell (+1, +1): (" + str(WT_max_cell[0]) + ", " + str(WT_max_cell[1]) + ")")
print("WT max value: " + str(WT_max))
print("MT relative position value: " + str(MT_relat_value) + "\n")  # I just invert the coordinate (col, row) instead of (row, col)


#4.Ex
# In the local region with altered RNA structure find the number of base pairs with pair probabilities 
# higher then 0.5 and 0.8. Do you see a consistent pattern?

def find_cell_with_higher_prob(matrice):
    """Take a matrix as an input and return a list of all the values 
    of the cells with a probability higher than threeshold"""
    wt_08plus = []
    wt_05plus = []
    mt_08plus = []
    mt_05plus = []

    for row, elem in enumerate(matrice):           # row is the position and element is the row_list that contain the values
        for col, value in enumerate(elem):         # the position in the row_list are the columns
            if row < col:
                if float(value) > 0.5:
                    wt_05plus.append(float(value))
                    if float(value) > 0.8:
                        wt_08plus.append(float(value))
            if row > col:
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
