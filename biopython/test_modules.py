

from Bio.Seq import Seq
my_seq = Seq("AGTACACTGGT")
print(my_seq)
print(my_seq.complement())
print(my_seq.reverse_complement())

from Bio import SeqIO
for seq_record in SeqIO.parse("/home/lpp/sb2019/biopython/ls_orchid.gbk", "genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))