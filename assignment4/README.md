# Assignment 3 - Measures to detect the eﬀect of SNPs on RNA secondary structure

The effect of SNPs on RNA secondary structure can be predicted by comparing the structures of wild-type and mutant (with SNP) RNA. The structures being compared may be either optimal (MFE) structure or ensemble structure (obtained from partition function).
While comparing the optimal structure between wild-type and mutant, the structures can be considered as two distinct strings and the diﬀerence at the pure string level be used to measure how divergent they are. This strategy is employed by the RNAmute program (Churkin and Barash, 2006):
* Hamming distance - The number of position at which the corresponding symbols are diﬀerent.
* base-pair distance - The total number of base pairs that are diﬀerent between two structures

Compute both the Hamming and base pair distance of the following pairs of sequence:
1. 
```
WT   GCGGGCCCCGC ((((...)))) 
MUT  ACGGGCCCCGC .(((...))).
```
2. 
```
WT   CAAUCCCGGCUGCGUCCCAGUUGGAUUUAUCCAGCUGGUUCGUGCUGGUU .....(((((.(((..(((((((((....)))))))))..)))))))).. 
MUT  CAAUCCCGGCUGCGUCCCAGUUGGAUUUAUCCAGCUGGUUCGUGGUGGUU ......(.((((((..(((((((((....)))))))))..)))))).)..
```
3. 
```python
WT   (((((..((((((((........(((((......)))))........)))))(((((...))))))))...)))))...((((((.((((((....)))))).).)))))..((((((...................))))))...((((((((((((.(((((((....))))))))))..((((((.....(((.((((((((.....))))))))....))).....))))))....))))))).))..
MUT  (((((..((((((((........(((((......)))))........)))))(((((...))))))))...)))))...((((((.((((((....)))))).).)))))..((((((...................))))))...(((((((((..(((((((..((((((...........))))))....))))))).....((((((....))))))...((......))......))))))).))..

```
