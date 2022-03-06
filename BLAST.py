from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

rat_BTG1 = SeqIO.read("BTG1_rat.fasta", format = "fasta") 
result_handle = NCBIWWW.qblast("blastn", "nt", rat_BTG1.seq)

blast_record = NCBIXML.read(result_handle)

E_VALUE_THRESH = 0.04
for alignment in blast_record.alignments:
  for hsp in alignment.hsps:
    if hsp.expect < E_VALUE_THRESH:
      print("***Alignment***")
      print("sequence:", alignment.title)
      print("length:", alignment.length)
      print("e value:", hsp.expect)