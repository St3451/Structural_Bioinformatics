# Assignment 3 - Measures to detect the eﬀect of SNPs on RNA secondary structure

The effect of SNPs on RNA secondary structure can be predicted by comparing the structures of wild-type and mutant (with SNP) RNA. The structures being compared may be either optimal (MFE) structure or ensemble structure (obtained from partition function).
While comparing the optimal structure between wild-type and mutant, the structures can be considered as two distinct strings and the diﬀerence at the pure string level be used to measure how divergent they are. This strategy is employed by the RNAmute program (Churkin and Barash, 2006):
* Hamming distance - The number of position at which the corresponding symbols are diﬀerent.
* base-pair distance - The total number of base pairs that are diﬀerent between two structures

Compute both the Hamming and base pair distance of the following pairs of sequence:  

__Pair 1__ 
```
WT   GCGGGCCCCGC ((((...)))) 
MUT  ACGGGCCCCGC .(((...))).
```
__Pair 2__
```
WT   CAAUCCCGGCUGCGUCCCAGUUGGAUUUAUCCAGCUGGUUCGUGCUGGUU .....(((((.(((..(((((((((....)))))))))..)))))))).. 
MUT  CAAUCCCGGCUGCGUCCCAGUUGGAUUUAUCCAGCUGGUUCGUGGUGGUU ......(.((((((..(((((((((....)))))))))..)))))).)..
```

__Pair 3__
```
WT   AGCGGGGGAGACAUAUAUCACAGCCUGUCUCGUGCCCGACCCCGCUGGUU .....(((((.(((..(((((((((....)))))))))..)))))))).. 
MUT  AGCGGGGGAGA*G*AUAUAUCACAGCCUGUCUCGUGCCCGACCCCGCUGGUU (((((((..((((((..........))))))........)))))))....
```
![#f03c15](https://placehold.it/15/f03c15/000000?text=ASDASDSAD) `#f03c15`

__Pair 3__ 
```python
WT   (((((..((((((((........(((((......)))))........)))))(((((...))))))))...)))))...((((((.((((((....)))))).).)))))..((((((...................))))))...((((((((((((.(((((((....))))))))))..((((((.....(((.((((((((.....))))))))....))).....))))))....))))))).))..
MUT  (((((..((((((((........(((((......)))))........)))))(((((...))))))))...)))))...((((((.((((((....)))))).).)))))..((((((...................))))))...(((((((((..(((((((..((((((...........))))))....))))))).....((((((....))))))...((......))......))))))).))..

```

## Implementation

##### Hamming distance
```
def calculate_hamming_distance(seq1_db, seq2_db):
    """Take two dot bracket sequences and return their hamming distance"""
    hamming_difference = 0
    for index in range(0, len(seq1_db)):
        if seq1_db[index] != seq2_db[index]:
            hamming_difference += 1
    return hamming_difference
```

##### Base-pair distance
```
def calculate_bp_distance(seq1_db, seq2_db):
    """Take two dot bracket sequences and return their base pairs distance"""
    lst_opening_1 = []                               
    lst_couples_1 = []
    for index in range(len(seq1_db)):
        if seq1_db[index] == "(":                     # A stack allows us to add or remove elements only from the top of it (last in first out):
            lst_opening_1.append(index)               # 1) if I have an open "(" append the position of the "(" to the lst_opening

        elif seq1_db[index] == ")":                   # 2) if I have a closed ")", remove (pop) the last element of lst_opening (position of "(" )
            elem = lst_opening_1.pop()                #    and create a bp_coordinate with: element removed from lst_opening (position opening bracket)
                                                      #                                     current position (position closing bracket)
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

    set1 = set(lst_couples_1)                         # set() can be used to check if there are differences between two list (two sets now) 
    set2 = set(lst_couples_2)                         # using .difference method
    diff = len(set1.difference(set2)) + len(set2.difference(set1))    # set1.difference(set2) output a set with elements present in
                                                                      # set1 but not in set2
    return diff
```

##### Pair 1
```
wt0 = "((((...))))"
mut0= ".(((...)))."
```

## 1
```
wt1 = ".....(((((.(((..(((((((((....)))))))))..)))))))).."
mut1= "......(.((((((..(((((((((....)))))))))..)))))).).."
```

wt2 = "(((((((((((((............))))))........)))))))...."
mut2= "(((((((..((((((..........))))))........)))))))...."

## 2
wt3 = "(((((..((((((((........(((((......)))))........)))))(((((...))))))))...)))))...((((((.((((((....)))))).).)))))..((((((...................))))))...((((((((((((.(((((((....))))))))))..((((((.....(((.((((((((.....))))))))....))).....))))))....))))))).)).."
mut3= "(((((..((((((((........(((((......)))))........)))))(((((...))))))))...)))))...((((((.((((((....)))))).).)))))..((((((...................))))))...(((((((((..(((((((..((((((...........))))))....))))))).....((((((....))))))...((......))......))))))).)).."

## Execute 
print("\n0.Ex\n>>>> WT sequence: " + wt0 + "\n>>>> MT sequence: " + mut0)
print("Hamming difference: " + str(calculate_hamming_distance(wt0, mut0)))
print("BP distance: " + str(calculate_bp_distance(wt0, mut0)))

print("\n1.Ex Case1\n>>>> WT sequence: " + wt1 + "\n>>>> MT sequence: " + mut1)
print("Hamming difference: " + str(calculate_hamming_distance(wt1, mut1)))
print("BP distance: " + str(calculate_bp_distance(wt1, mut1)))

print("\n1.Ex Case2\n>>>> WT sequence: " + wt2 + "\n>>>> MT sequence: " + mut2)
print("Hamming difference: " + str(calculate_hamming_distance(wt2, mut2)))
print("BP distance: " + str(calculate_bp_distance(wt2, mut2)))

print("\n2.Ex Case2\n>>>> WT sequence: " + wt3 + "\n>>>> MT sequence: " + mut3)
print("Hamming difference: " + str(calculate_hamming_distance(wt3, mut3)))
print("BP distance: " + str(calculate_bp_distance(wt3, mut3)))