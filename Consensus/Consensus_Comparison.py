import collections
from DNA_Reverse_Compliment import revCompIterative

from OrfFinderComparison_STORF_FULL_NT import orfComparison

#ConsensusInput = open('./Direct_Consensus_Comparison/Consensus_ORFs.txt', mode='rb')
ConsensusInput = open('./Majority_Data/Myco_Majority.txt', mode='rb')

ConsensusOutput = open('./Direct_Consensus_Comparison/Consensus_Metrics_Majority.csv', mode='wb')

Genome_Input = open('/home/nick/Git/Open_Reading_Frame_Comparison/Myco.txt', mode='rb')


Alignments = open('/home/nick/Git/Open_Reading_Frame_Comparison/Myco.gff', mode='rb')

Genome = ""

for line in Genome_Input:
    line = line.replace("\n","")
    if ">" not in line:
        Genome += str(line)


##########################################


Genes = collections.OrderedDict()
ConsensusORFs = collections.OrderedDict()
ConsensusORFsFiltered = collections.OrderedDict()


ConsensusLengths = []
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
Prev_Stop = 0
OL = False
Genome_Size = len(Genome)
Genome_rev = revCompIterative(Genome)

for line in ConsensusInput:

    Start = int(line.split(',')[0])
    Stop = int(line.split(',')[1])

    reverse = '-'
    forward = '+'
    if '-' in line.split(',')[2]: #Reverse Compliment starts and stops to confirm to our deffinition
        #Switched to match Sense Strand
        r_start = Genome_Size - Stop

        po = str(r_start)+','+str(Stop)
        Start_Codons.update({po:reverse})
    elif '+' in line.split(',')[2]:
        po = str(Start)+','+str(Stop)
        Start_Codons.update({po:forward})


    ConsensusORFs.update({Start:Stop})

Output, Start_Precision, Stop_Precision = orfComparison(Genes, ConsensusORFs, Start_Codons, Genome)

ConsensusOutput.write(str(Output)+'\n')
ConsensusOutput.write(str(Start_Precision)+'\n')
ConsensusOutput.write(str(Stop_Precision)+'\n')
