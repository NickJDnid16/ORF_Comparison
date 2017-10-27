import collections
from DNA_Reverse_Compliment import revCompIterative
from DNA_Reverse_Compliment import revCompIterative
from OrfFinderComparison import orfComparison

GLIMMERInput = open('./GLIMMER_K-12.txt', mode='rb')

GLIMMEROutput = open('./GLIMMER_Metrics.csv', mode='wb')


Genome_Input = open('/home/nick/Git/Open_Reading_Frame_Comparison/K-12.txt', mode='rb')


Alignments = open('/home/nick/Git/Open_Reading_Frame_Comparison/K-12.gff', mode='rb')

Genome = ""

for line in Genome_Input:
    line = line.replace("\n","")
    if ">" not in line:
        Genome += str(line)


##############################################

GLIMMERORFs = collections.OrderedDict()
GLIMMERORFsFiltered = collections.OrderedDict()
Genes = collections.OrderedDict()

GLIMMERLengths = []
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
Starts = []
Stops = []
for line in GLIMMERInput:
    if "orf0" in line :
        Start = int(line.split()[1])
        Stop = int(line.split()[2])
        if Start not in Starts:
            Starts.append(Start)
        else:
            print "STARTS"
        if Stop not in Stops:
            Stops.append(Stop)
        else:
            print "STOPS"
        if Start < Prev_Stop and OL == False:
            print "Overlapping"
            OL = True
        Prev_Stop = Stop
        reverse = '-'
        forward = '+'
        if '-' in line.split()[3]: #Reverse Compliment starts and stops to confirm to our deffinition
            #Switched to match Sense Strand
            r_start = Genome_Size - Stop

            r_stop = Genome_Size - Start
            po = str(r_stop)+','+str(Start)
            Start_Codons.update({po:reverse})
            GLIMMERORFs.update({Stop: Start})
            po = str(Stop) + ',' + str(r_start)
            Stop_Codons.update({po: reverse})
        elif '+' in line.split()[3]:
            po = str(Start)+','+str(Stop)

            Start_Codons.update({po:forward})
            Stop_Codons.update({po:forward})

            GLIMMERORFs.update({Start:Stop})

Output, Start_Precision, Stop_Precision = orfComparison(Genes, GLIMMERORFs, Start_Codons,Stop_Codons, Genome)

GLIMMEROutput.write(str(Output) + '\n')
GLIMMEROutput.write(str(Start_Precision) + '\n')
GLIMMEROutput.write(str(Stop_Precision) + '\n')
















