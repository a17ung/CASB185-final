from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.Align import substitution_matrices
import matplotlib.pyplot as plt
# for seq_record in SeqIO.parse("Downloads/BTG1_rat.fasta", "fasta"):
#   print(seq_record.id)
#   print(repr(seq_record.seq))
#   print(len(seq_record))
# record = SeqIO.read("Downloads/BTG1_rat.fasta", "fasta")
# print(record.id)
# print(record.description)
# #print(record.seq)
# print(record.translate())

# rat = SeqIO.read("Rat_ST6.fasta", "fasta")
# human = SeqIO.read("Human_ST6.fasta", "fasta")
# alignments = pairwise2.align.localxx(rat.seq, human.seq, score_only = True)
# print(alignments)

#read in sequence info (Fasta format)
#separate by space (NOT COMMA) (ie. "BTG1_rat.fasta ST6_rat.fasta")
def getSeq():
  r_seq = input("Enter rat sequences: ").split()
  h_seq = input("Enter human sequences: ").split()
  genes = input("Enter gene names: ").split()
  return r_seq, h_seq, genes

#convert sequences to SeqRecord objects
def convertSeq(r, h):
  rat = list()
  human = list()

  for i in range(len(r)):
    rat.append(SeqIO.read(r[i], "fasta"))
    human.append(SeqIO.read(h[i], "fasta"))
  return rat, human

#calculate alignment scores
def align(r, h):
  alignments = list()
  for i in range(len(r)):
    alignments.append(pairwise2.align.localxx(r[i].seq, h[i].seq, score_only = True))
  return alignments

#VARIABLES
rat_seq, human_seq, gene_names = getSeq()
rat, human = convertSeq(rat_seq, human_seq)
alignments = align(rat, human)
print(alignments)

#GRAPHS
left = list(range(len(alignments))) #x-coordinates
height = alignments #y-coordinates

tick_label = gene_names 

#plot graph
plt.bar(left, height, tick_label = tick_label, width = 0.8)
plt.xlabel('gene')
plt.ylabel('alignment score')
plt.title('Alignment Scores of Genes')
plt.show()
