import collections
from DNA_Reverse_Compliment import revCompIterative


import numpy as np
from OrfFinderComparison import orfComparison

GetORFInput = open('./GetORF_K-12.txt', mode='rb')

GetORFOutput = open('./GetORF_Metrics.csv', mode='wb')


Genome_Input = open('/home/nick/Git/Open_Reading_Frame_Comparison/K-12.txt', mode='rb')


Alignments = open('/home/nick/Git/Open_Reading_Frame_Comparison/K-12.gff', mode='rb')

Genome = ""

for line in Genome_Input:
    line = line.replace("\n","")
    if ">" not in line:
        Genome += str(line)


##############################################

GetORFORFs = collections.OrderedDict()
GetORFORFsFiltered = collections.OrderedDict()
Genes = collections.OrderedDict()

GetORFLengths = []
############################
count = 0
for line in Alignments:
    if "Chromosome	ena	CDS" in line:
        count += 1
        Start = int(line.split('\t')[3])

        Stop = int(line.split('\t')[4])
        Gene = str(Start)+','+str(Stop)
        Genes.update({count:Gene})

##################################
Start_Codons = collections.OrderedDict()
Prev_Stop = 0
OL = False
Genome_Size = len(Genome)
Genome_rev = revCompIterative(Genome)

Seen = []
for line in GetORFInput:
    if ">Chromosome" in line :


        Temp = line.split('_')[1]
        loci = Temp.partition('[')[-1].rpartition(']')[0]
        po1 = int(loci.split('-')[0])
        po2 = int(loci.split('-')[1])
        # if po1 in Seen:
        #     print "sss"
        # if po2 in Seen:
        #     print "ddd"
        if po2 > po1:
            po = str(po1) + ',' + str(po2)
            Start_Codons.update({po:"+"})
            GetORFORFs.update({po1: po2})
        elif po1 > po2:
            r_start = Genome_Size - po2
            po = str(r_start) + ',' + str(po2)
            Start_Codons.update({po: "-"})
            GetORFORFs.update({po2: po1})
        Seen.append(po1)
        Seen.append(po2)










Output,Start_Precision,Stop_Precision = orfComparison(Genes,GetORFORFs,Start_Codons,Genome)

GetORFOutput.write(str(Output)+'\n')
GetORFOutput.write(str(Start_Precision)+'\n')
GetORFOutput.write(str(Stop_Precision)+'\n')



























