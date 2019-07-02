import collections
from DNA_Reverse_Compliment import revCompIterative
from DNA_Reverse_Compliment import revCompIterative
from OrfFinderComparison import orfComparison

def prodigal(Genome):

ProdigalORFs = collections.OrderedDict()
ProdigalORFsFiltered = collections.OrderedDict()

ProdigalLengths = []
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
















