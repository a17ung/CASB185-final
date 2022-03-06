from Bio import SeqIO
from Bio.Seq import Seq
import pylab

window = 7
rat_record = SeqIO.read("BTG1_rat.fasta", "fasta")
human_record = SeqIO.read("BTG1_human.fasta", "fasta")

rat_seq = rat_record.seq.upper()
hum_seq = human_record.seq.upper()

dict_rat = {}
dict_human = {}

for(seq, section_dict) in [ 
  (rat_seq, dict_rat),
  (hum_seq, dict_human),
]:
  for i in range(len(seq) - window):
    section = seq[i : i + window]
    try:
      section_dict[section].append(i)
    except KeyError:
      section_dict[section] = [i]
matches = set(dict_rat).intersection(dict_human)
print("%i unique matches" % len(matches))

x = []
y = []
for section in matches:
  for i in dict_rat[section]:
    for j in dict_human[section]:
      x.append(i)
      y.append(j)

pylab.gray()
pylab.scatter(x, y)
pylab.xlim(0, len(rat_record) - window)
pylab.ylim(0, len(human_record) - window)
pylab.xlabel("%s (length %i bp)" % (rat_record.id, len(rat_record)))
pylab.ylabel("%s (length %i bp)" % (human_record.id, len(human_record)))
pylab.title("Dot plot using window size %i\n(allowing no mis-matches)" % window)
pylab.show()

# data = [
#   [
#     (rat_seq[i: i + window] != hum_seq[j: j + window])
#     for j in range(len(rat_seq) - window)
#   ]
#   for i in range(len(hum_seq) - window)
# ]

# pylab.gray()
# pylab.imshow(data)
# pylab.xlabel("%s (length %i bp)" % (rat_record.id, len(rat_record)))
# pylab.ylabel("%s (length %i bp)" % (human_record.id, len(human_record)))
# pylab.title("Dot plot using window size %i\n(allowing no mismatches)"% window)
# pylab.show()