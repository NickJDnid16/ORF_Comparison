import collections
from DNA_Reverse_Compliment import revCompIterative


import numpy as np
from OrfFinderComparison import orfComparison
import linecache
orfMInput = open('./orfM_K-12.txt', mode='rb')

orfMOutput = open('./orfM_Metrics.csv', mode='wb')


Genome_Input = open('/home/nick/Git/Open_Reading_Frame_Comparison/K-12.txt', mode='rb')


Alignments = open('/home/nick/Git/Open_Reading_Frame_Comparison/K-12.gff', mode='rb')

Genome = ""

for line in Genome_Input:
    line = line.replace("\n","")
    if ">" not in line:
        Genome += str(line)


##############################################
#orfMORFsTemps = collections.OrderedDict()
orfMORFs = collections.OrderedDict()
orfMFiltered = collections.OrderedDict()
Genes = collections.OrderedDict()

orfMLengths = []
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
count = 0

for line in orfMInput:
    if ">Chromosome" in line :

        count += 1
        Start = int(line.split('_')[1])
        AA = linecache.getline('orfM_K-12.txt',count+1).rstrip('\n')
        Stop = (len(AA)*3) + Start
        count += 1
        if Start < Prev_Stop and OL == False:
            print "Overlapping"
            OL = True


        reverse = '4', '5', '6'
        forward = '1', '2', '3'
        if line.split('_')[2] in forward:
            Strand = '+'
            po = str(Start)+','+str(Stop)
            Start_Codons.update({po:Strand})

        elif line.split('_')[2] in reverse:
            Strand = '-'
            r_start = Genome_Size - Stop
            po = str(r_start)+','+str(Stop)
            Start_Codons.update({po:Strand})

        Prev_Stop = Stop

        orfMORFs.update({Start:Stop})


########################################################### Remove nested and overlapping
# print "Started"
# for T_Pos, T_Strand in orfMORFsTemps.items(): #Remove shorter versions of same gene
#     reverse = '4', '5', '6'
#     forward = '1', '2', '3'
#     if T_Strand in forward:
#         T_Start = int(T_Pos.split(',')[0])
#         T_Stop = int(T_Pos.split(',')[1])
#         po = str(T_Start) + ',' + str(T_Stop)
#         Strand = '+'
#
#     elif T_Strand in reverse:
#         T_Start = int(T_Pos.split(',')[0])
#         T_Stop = int(T_Pos.split(',')[1])
#
#         po = str(T_Start) + ',' + str(T_Stop)
#         Strand = '-'
#
#
#     if len(orfMORFs) == 0:
#         orfMORFs.update({T_Start: T_Stop})
#         Start_Codons.update({po:Strand})
#
#
#     for G_Start, G_Stop in orfMORFs.items():
#         if T_Start <= G_Start and T_Stop > G_Start:
#             orfMORFs.update({T_Start:T_Stop})
#
#             Start_Codons.update({po: Strand})
#             # if len(orfMORFs) > 0:
#             #     del orfMORFs[G_Start]
#
#         elif T_Start >= G_Start and T_Start <= G_Stop:
#             orfMORFs.update({T_Start:T_Stop})
#
#             Start_Codons.update({po: Strand})
#             # if len(orfMORFs) > 0:
#             #     del orfMORFs[G_Start]
#
#
#
#         else:
#             print"Derped"
#             print T_Start
#
#
# print "Done"
Output,Start_Precision,Stop_Precision = orfComparison(Genes,orfMORFs,Start_Codons,Genome)


orfMOutput.write(str(Output)+'\n')
orfMOutput.write(str(Start_Precision)+'\n')
orfMOutput.write(str(Stop_Precision)+'\n')































