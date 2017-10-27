import collections
from DNA_Reverse_Compliment import revCompIterative
from OrfFinderComparison import orfComparison


Input = open('./Easy_Gene_Staph_Pseudo.txt', mode='rb')

EasyGeneOutput = open('./Easy_Gene_Metrics.csv', mode='wb')

Genome_Input = open('/home/nick/Git/Open_Reading_Frame_Comparison/Pseudo.txt', mode='rb')

Alignments = open('/home/nick/Git/Open_Reading_Frame_Comparison/Pseudo.gff', mode='rb')

Genome = ""

for line in Genome_Input:
    line = line.replace("\n","")
    if ">" not in line:
        Genome += str(line)


##############################################


Genes = collections.OrderedDict()
EasyGeneORFs = collections.OrderedDict()
Lengths = []
############################
count = 0
for line in Alignments:
    if "Chromosome	ena	gene" in line :
        count += 1
        Start = int(line.split('\t')[3])

        Stop = int(line.split('\t')[4])
        Gene = str(Start) + ',' + str(Stop)
        Genes.update({count: Gene})
Genome_Size = len(Genome)
Genome_rev = revCompIterative(Genome)

##################################
Start_Codons = collections.OrderedDict()
Prev_Stop = 0
OL = False
Genome_Size = len(Genome)
Genome_rev = revCompIterative(Genome)

for line in Input:
    if "Chromosome	SAM03	CDS" in line or "Chromosome	EC03	CDS" in line:
        Start = int(line.split('\t')[3])
        Stop = int(line.split('\t')[4])
        if Start < Prev_Stop and OL == False:
            print "Overlapping"
            OL = True
        Prev_Stop = Stop
        reverse = '-'
        forward = '+'
        if '-' in line.split('\t')[6]: #Reverse Compliment starts and stops to confirm to our deffinition
            #Switched to match Sense Strand
            r_start = Genome_Size - Stop

            po = str(r_start)+','+str(Stop)
            Start_Codons.update({po:reverse})
        elif '+' in line.split('\t')[6]:
            po = str(Start)+','+str(Stop)

            Start_Codons.update({po:forward})

        EasyGeneORFs.update({Start:Stop})


Output,Start_Precision,Stop_Precision = orfComparison(Genes,EasyGeneORFs,Start_Codons,Genome)


EasyGeneOutput.write(str(Output)+'\n')
EasyGeneOutput.write(str(Start_Precision)+'\n')
EasyGeneOutput.write(str(Stop_Precision)+'\n')

















