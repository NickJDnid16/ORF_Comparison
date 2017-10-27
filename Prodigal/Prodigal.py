import collections
from DNA_Reverse_Compliment import revCompIterative
from DNA_Reverse_Compliment import revCompIterative
from OrfFinderComparison import orfComparison

ProdigalInput = open('./Prodigal_K-12.txt', mode='rb')

ProdigalOutput = open('./Prodigal_Metrics.csv', mode='wb')


Genome_Input = open('/home/nick/Git/Open_Reading_Frame_Comparison/K-12.txt', mode='rb')


Alignments = open('/home/nick/Git/Open_Reading_Frame_Comparison/K-12.gff', mode='rb')

Genome = ""

for line in Genome_Input:
    line = line.replace("\n","")
    if ">" not in line:
        Genome += str(line)


##############################################

ProdigalORFs = collections.OrderedDict()
ProdigalORFsFiltered = collections.OrderedDict()
Genes = collections.OrderedDict()

ProdigalLengths = []
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
Stop_Codons = collections.OrderedDict()
Prev_Stop = 0
OL = False
Genome_Size = len(Genome)
Genome_rev = revCompIterative(Genome)

for line in ProdigalInput:
    if "CDS" in line :
        line = line.replace("<","")
        if 'complement' in line: #Reverse Compliment starts and stops to confirm to our deffinition
            #Switched to match Sense Strand

            loci = line.split()[1]
            loci = loci.partition('(')[-1].rpartition(')')[0]
            Start = int(loci.split('..')[0])
            Stop = int(loci.split('..')[1])
            r_start = Genome_Size - Stop
            r_stop = Genome_Size - Start
            po = str(r_start)+','+str(Stop)
            Start_Codons.update({po:'-'})
            po = str(Start) + ',' + str(r_stop)
            Stop_Codons.update({po: '-'})
            ProdigalORFs.update({Start: Stop})
        else:
            loci = line.split()[1]
            Start = int(loci.split('..')[0])
            Stop = int(loci.split('..')[1])
            po = str(Start)+','+str(Stop)

            Start_Codons.update({po:'+'})
            Stop_Codons.update({po: '+'})

            ProdigalORFs.update({Start:Stop})

Output,Start_Precision,Stop_Precision = orfComparison(Genes,ProdigalORFs,Start_Codons,Stop_Codons,Genome)

ProdigalOutput.write(str(Output) + '\n')
ProdigalOutput.write(str(Start_Precision) + '\n')
ProdigalOutput.write(str(Stop_Precision) + '\n')
















