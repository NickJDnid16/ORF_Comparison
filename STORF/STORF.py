import collections
from DNA_Reverse_Compliment import revCompIterative
from DNA_Reverse_Compliment import revCompIterative
from OrfFinderComparison_STORF_FP import orfComparison
import sys


STORFInput = open('/home/nick/Git/Open_Reading_Frame_Comparison/STORF/Myco/Myco_100_TAA_TAG.txt', mode='rb')

STORFOutput = open('./STORF_Metrics_Meh.csv', mode='wb')


Genome_Input = open('/home/nick/Git/Open_Reading_Frame_Comparison/Myco.txt', mode='rb')


Alignments = open('/home/nick/Git/Open_Reading_Frame_Comparison/Myco.gff', mode='rb')

Genome = ""

for line in Genome_Input:
    line = line.replace("\n","")
    if ">" not in line:
        Genome += str(line)


##############################################

STORFORFs = collections.OrderedDict()
STORFORFsFiltered = collections.OrderedDict()
Genes = collections.OrderedDict()

STORFLengths = []
############################
count = 0
for line in Alignments:
    if "Chromosome	ena	gene" in line:
        count += 1
        Start = int(line.split('\t')[3])
        #Start = Start-3
        Stop = int(line.split('\t')[4])
        Gene = str(Start)+','+str(Stop)
        Genes.update({count:Gene})

##################################
Start_Codons = collections.OrderedDict()
Prev_Stop = 0
OL = False
Genome_Size = len(Genome)
Genome_rev = revCompIterative(Genome)




#for line in STORFInput:
#     if "+" in line or '-' in line :
#         Start = int(line.split()[0])
#         Stop = int(line.split()[2])
#         if Start < Prev_Stop and OL == False:
#             print "Overlapping"
#             OL = True
#         Prev_Stop = Stop
#         reverse = '-'
#         forward = '+'
#         if '-' in line: #Reverse Compliment starts and stops to confirm to our deffinition
#             #Switched to match Sense Strand
#             r_stop = Genome_Size - Stop
#
#             po = str(r_stop)+','+str(Stop)
#             Start_Codons.update({po:reverse})
#         elif '+' in line:
#             po = str(Start)+','+str(Stop)
#             Start_Codons.update({po:forward})
#
#
#         STORFORFs.update({Start:Stop})
#
# Output, Start_Precision, Stop_Precision = orfComparison(Genes, STORFORFs, Start_Codons, Genome)

for line in STORFInput:
    if "+" in line or "-" in line:
        line = line.split(',')[1]
        Start = int(line.split()[0])
        Stop = int(line.split()[2])
        if Start < Prev_Stop and OL == False:
            print "Overlapping"
            OL = True
        Prev_Stop = Stop
        reverse = '-'
        forward = '+'

        if '+' in line:
            po = str(Start) + ',' + str(Stop)
            Start_Codons.update({po: forward})

        elif '-' in line: #Reverse Compliment starts and stops to confirm to our deffinition
         #Switched to match Sense Strand
            r_start = Genome_Size - Stop



            po = str(Start) + ',' + str(Stop)
            Start_Codons.update({po:reverse})

        STORFORFs.update({po:Stop})



Output, Start_Precision, Stop_Precision = orfComparison(Genes, STORFORFs, Start_Codons, Genome)









STORFOutput.write(str(Output) + '\n')
STORFOutput.write(str(Start_Precision) + '\n')
STORFOutput.write(str(Stop_Precision) + '\n')
















