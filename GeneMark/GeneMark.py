import collections
from DNA_Reverse_Compliment import revCompIterative
from OrfFinderComparison import orfComparison

GeneMarkInput = open('./GeneMark_Staph_K-12.txt', mode='rb')

GeneMarkOutput = open('./GeneMark_Metrics.csv', mode='wb')


Genome_Input = open('/home/nick/Git/Open_Reading_Frame_Comparison/K-12.txt', mode='rb')


Alignments = open('/home/nick/Git/Open_Reading_Frame_Comparison/K-12.gff', mode='rb')

Genome = ""

for line in Genome_Input:
    line = line.replace("\n","")
    if ">" not in line:
        Genome += str(line)


##############################################

GeneMarkORFsTemp = collections.OrderedDict()
GeneMarkORFs = collections.OrderedDict()
GeneMarkORFsFiltered = collections.OrderedDict()
Genes = collections.OrderedDict()

GeneMark_Lengths = []
############################
count = 0
for line in Alignments:
    if "Chromosome	ena	gene" in line:
        count += 1
        Start = int(line.split('\t')[3])
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
Num = 0
for line in GeneMarkInput:
    if "fr " in line and "0." in line:
        Start = int(line.split()[0])
        Stop = int(line.split()[1])
        if Start < Prev_Stop and OL == False:
            print "Overlapping"
            OL = True
        Prev_Stop = Stop
        Pos = str(Start)+','+str(Stop)
        GeneMarkORFsTemp.update({Pos:line.split()[2]})
        Num = Num+1
print "Num"
print Num

########################################################### Remove nested and overlapping
GeneMarkORFs = collections.OrderedDict()
SeenStops = []
count = 0





for T_Pos, T_Strand in GeneMarkORFsTemp.items(): #Remove shorter versions of same gene
    count +1
    if 'direct' in T_Strand:
        T_Start = int(T_Pos.split(',')[0])
        T_Stop = int(T_Pos.split(',')[1])
        po = str(T_Start) + ',' + str(T_Stop)
        po_s = po
        Strand = '+'
    elif 'complement' in T_Strand:
        T_Start = int(T_Pos.split(',')[0])
        T_Stop = int(T_Pos.split(',')[1])

        r_start = Genome_Size - T_Stop
        r_stop = Genome_Size - Start
        po = str(r_start) + ',' + str(T_Stop)
        po_s = str(Start) + ',' + str(r_stop)
        Strand = '-'
    if len(GeneMarkORFs) == 0:
        GeneMarkORFs.update({T_Start: T_Stop})
        Start_Codons.update({po:Strand})

    if T_Stop not in SeenStops:
        SeenStops.append(T_Stop)
        GeneMarkORFs.update({T_Start:T_Stop})
        Start_Codons.update({po: Strand})
        Stop_Codons.update({po_s: Strand})
    if count == 1000:
        print "1000"

print str(Start_Codons)
print str(GeneMarkORFs)
Output,Start_Precision,Stop_Precision = orfComparison(Genes,GeneMarkORFs,Start_Codons,Stop_Codons,Genome)


GeneMarkOutput.write(str(Output) + '\n')
GeneMarkOutput.write(str(Start_Precision) + '\n')
GeneMarkOutput.write(str(Stop_Precision) + '\n')


##################################################










