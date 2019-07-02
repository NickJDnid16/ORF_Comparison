import collections
from DNA_Reverse_Compliment import revCompIterative

from OrfFinderComparison import orfComparison

ProkkaInput = open('./PROKKA_K-12.gff', mode='rb')

ProkkaOutput = open('./Prokka_Metrics.csv', mode='wb')

Genome_Input = open('/home/nick/Git/Open_Reading_Frame_Comparison/K-12.txt', mode='rb')


Alignments = open('/home/nick/Git/Open_Reading_Frame_Comparison/K-12.gff', mode='rb')

Genome = ""

for line in Genome_Input:
    line = line.replace("\n","")
    if ">" not in line:
        Genome += str(line)


##########################################


Genes = collections.OrderedDict()
ProkkaORFs = collections.OrderedDict()
ProkkaORFsFiltered = collections.OrderedDict()


ProkkaLengths = []
count = 0

for line in Alignments:
    if "Chromosome	ena	gene" in line:
        count += 1
        Start = int(line.split('\t')[3])

        Stop = int(line.split('\t')[4])
        Gene = str(Start)+','+str(Stop)
        Genes.update({count:Gene})



#########################################################
Start_Codons = collections.OrderedDict()
Stop_Codons = collections.OrderedDict()
Prev_Stop = 0
OL = False
Genome_Size = len(Genome)
Genome_rev = revCompIterative(Genome)

for line in ProkkaInput:
    if "Chromosome	prodigal:2.6	CDS" in line :
        Start = int(line.split()[3])
        Stop = int(line.split()[4])
        if Start < Prev_Stop and OL == False:
            print "Overlapping"
            OL = True
        Prev_Stop = Stop
        reverse = '-'
        forward = '+'
        if '-' in line.split()[6]: #Reverse Compliment starts and stops to confirm to our deffinition
            #Switched to match Sense Strand
            r_start = Genome_Size - Stop
            r_stop = Genome_Size - Start
            po = str(r_start) + ',' + str(Stop)
            Start_Codons.update({po: reverse})
            po = str(Start) + ',' + str(r_stop)
            Stop_Codons.update({po: reverse})
        elif '+' in line.split()[6]:
            po = str(Start)+','+str(Stop)
            Stop_Codons.update({po: forward})
            Start_Codons.update({po:forward})

        ProkkaORFs.update({Start:Stop})

Output, Start_Precision, Stop_Precision = orfComparison(Genes, ProkkaORFs, Start_Codons,Stop_Codons, Genome)

ProkkaOutput.write(str(Output)+'\n')
ProkkaOutput.write(str(Start_Precision)+'\n')
ProkkaOutput.write(str(Stop_Precision)+'\n')