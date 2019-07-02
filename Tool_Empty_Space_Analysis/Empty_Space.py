import collections
from DNA_Reverse_Compliment import revCompIterative
from OrfFinderComparison_STORF_NT_Empty import orfComparison


ORF_NUM = 0

Genome_Input = open('/home/nick/Git/Open_Reading_Frame_Comparison/K-12.txt', mode='rb')

Alignments = open('/home/nick/Git/Open_Reading_Frame_Comparison/K-12.gff', mode='rb')

Genome = ""

for line in Genome_Input:
    line = line.replace("\n","")
    if ">" not in line:
        Genome += str(line)

Genome_Size = len(Genome)
Genome_rev = revCompIterative(Genome)

Combi_ORFs = collections.OrderedDict()

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


######################################################## EasyGene
Input = open('./Easy_Gene_K-12_K-12.txt', mode='rb')
Prev_Stop = 0
OL = False
EasyGeneOutput = open('./Easy_Gene_Metrics.csv', mode='wb')


EasyGeneORFs = collections.OrderedDict()

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
            r_stop = Genome_Size - Start

            po = str(r_start)+','+str(Stop)

            po = str(Start) + ',' + str(r_stop)

        elif '+' in line.split('\t')[6]:
            po = str(Start)+','+str(Stop)



        Combi_ORFs.update({ORF_NUM:str(str(Start) + ',' + str(Stop))})
        ORF_NUM += 1

 ######################################### Augustus
AugustusORFs = collections.OrderedDict()
Prev_Stop = 0
OL = False
AugustusInput = open('./Augustus_K-12.txt', mode='rb')
for line in AugustusInput:
    if "Chromosome	AUGUSTUS	CDS" in line :
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
            r_stop = Genome_Size - Start
            po = str(r_start)+','+str(Stop)

            po = str(Start) + ',' + str(r_stop)

        elif '+' in line.split('\t')[6]:
            po = str(Start)+','+str(Stop)

        Combi_ORFs.update({ORF_NUM:str(str(Start) + ',' + str(Stop))})
        ORF_NUM += 1


############################################## FragGeneScan
FragGeneScanInput = open('./FragGeneScan_K-12.txt', mode='rb')
FragGeneScanORFs = collections.OrderedDict()
Prev_Stop = 0
OL = False
for line in FragGeneScanInput:
    if "Chromosome	FGS	CDS" in line :
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

            po = str(r_start)+','+str(Stop)

        elif '+' in line.split()[6]:
            po = str(Start)+','+str(Stop)

        Combi_ORFs.update({ORF_NUM:str(str(Start) + ',' + str(Stop))})
        ORF_NUM += 1
############################################## GeneMark
GeneMarkInput = open('./GeneMark_K-12.txt', mode='rb')
GeneMarkORFsTemp = collections.OrderedDict()
GeneMarkORFs = collections.OrderedDict()
Prev_Stop = 0
OL = False
GeneMarkORFsFiltered = collections.OrderedDict()
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

#Remove nested and overlapping
GeneMarkORFs = collections.OrderedDict()
SeenStops = []
Prev_Stop = 0
OL = False
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
        po = str(r_start) + ',' + str(T_Stop)
        Strand = '-'
    if len(GeneMarkORFs) == 0:
        Combi_ORFs.update({ORF_NUM:str(T_Start) + ',' + str(T_Stop)})
        ORF_NUM += 1


    if T_Stop not in SeenStops:
        SeenStops.append(T_Stop)
        Combi_ORFs.update({ORF_NUM:str(T_Start) + ',' + str(T_Stop)})
        ORF_NUM += 1

    if count == 1000:
        print "1000"

#######################################GeneMarkS
GeneMark_S_Input = open('./GeneMark_S_K-12.txt', mode='rb')
GeneMark_S_ORFsTemp = collections.OrderedDict()
GeneMark_S_ORFs = collections.OrderedDict()
Prev_Stop = 0
OL = False
for line in GeneMark_S_Input:
    if "Chromosome" in line:
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

        elif '+' in line.split('\t')[6]:
            po = str(Start)+','+str(Stop)

        Combi_ORFs.update({ORF_NUM:str(Start) + ',' + str(Stop)})
        ORF_NUM += 1

##############################################GeneMarkHA
GeneMark_HA_Input = open('./GeneMark_HA_K-12.txt', mode='rb')

OL = False
GeneMark_HA_ORFsTemp = collections.OrderedDict()
GeneMark_HA_ORFs = collections.OrderedDict()
GeneMark_HA_ORFsFiltered = collections.OrderedDict()
Prev_Stop = 0
for line in GeneMark_HA_Input:
    if "Chromosome" in line:
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

        elif '+' in line.split('\t')[6]:
            po = str(Start)+','+str(Stop)

        Combi_ORFs.update({ORF_NUM:str(Start) + ',' + str(Stop)})
        ORF_NUM += 1


##################################GeneMarkHMM
GeneMark_HMMInput = open('./GeneMark_HMM-K-12.txt', mode='rb')
Prev_Stop = 0
OL = False
GeneMark_HMMORFs = collections.OrderedDict()
GeneMark_HMMORFsFiltered = collections.OrderedDict()
for line in GeneMark_HMMInput:
    if "Chromosome" in line :
        Start = int(line.split('\t')[3])
        Stop = int(line.split('\t')[4])
        if Start < Prev_Stop and OL == False:
            print "Overlapping"
            OL = True
        Prev_Stop = Stop
        reverse = '-'
        forward = '+'
        if '-' in line.split('\t')[6]:  # Reverse Compliment starts and stops to confirm to our deffinition
            # Switched to match Sense Strand
            r_start = Genome_Size - Stop
            r_stop = Genome_Size - Start
            po = str(r_start) + ',' + str(Stop)

        elif '+' in line.split('\t')[6]:
            po = str(Start) + ',' + str(Stop)

        Combi_ORFs.update({ORF_NUM:str(Start) + ',' + str(Stop)})
        ORF_NUM += 1

########################################~GLIMMER
GLIMMERInput = open('./GLIMMER_K-12.txt', mode='rb')
Prev_Stop = 0
OL = False
GLIMMERORFs = collections.OrderedDict()
GLIMMERORFsFiltered = collections.OrderedDict()
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
            r_stop = Genome_Size - Start

            po = str(r_stop)+','+str(Start)

            Combi_ORFs.update({ORF_NUM: str(Stop) + ',' + str(Start)})
            ORF_NUM += 1
        elif '+' in line.split()[3]:
            po = str(Start)+','+str(Stop)

            Combi_ORFs.update({ORF_NUM: str(Start) + ',' + str(Stop)})
            ORF_NUM += 1


####################################MetaGeneAnotator

MetaGeneAnnotatorInput = open('./MetaGeneAnnotator_K-12.txt', mode='rb')
MetaGeneAnnotator_ORFs = collections.OrderedDict()
MetaGeneAnnotator_ORFsFiltered = collections.OrderedDict()
Prev_Stop = 0
OL = False
for line in MetaGeneAnnotatorInput:
    if "gene_" in line:
        Start = int(line.split('\t')[1])
        Stop = int(line.split('\t')[2])
        if Start < Prev_Stop and OL == False:
            print "Overlapping"
            OL = True
        Prev_Stop = Stop
        reverse = '-'
        forward = '+'
        if '-' in line.split('\t')[3]: #Reverse Compliment starts and stops to confirm to our deffinition
            #Switched to match Sense Strand
            r_start = Genome_Size - Stop

            po = str(r_start)+','+str(Stop)

        elif '+' in line.split('\t')[3]:
            po = str(Start)+','+str(Stop)

        Combi_ORFs.update({ORF_NUM:str(Start) + ',' + str(Stop)})
        ORF_NUM += 1

################################MetaGeneMark

MetaGeneMarkInput = open('./MetaGeneMark_K-12.txt', mode='rb')
MetaGeneMarkORFs = collections.OrderedDict()
MetaGeneMarkORFsFiltered = collections.OrderedDict()
Prev_Stop = 0
OL = False
for line in MetaGeneMarkInput:
    if "Chromosome" in line:
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

        elif '+' in line.split('\t')[6]:
            po = str(Start)+','+str(Stop)

        Combi_ORFs.update({ORF_NUM:str(Start) + ',' + str(Stop)})
        ORF_NUM += 1

#######################################prodigal
ProdigalInput = open('./Prodigal_K-12.txt', mode='rb')

ProdigalORFs = collections.OrderedDict()
ProdigalORFsFiltered = collections.OrderedDict()
Prev_Stop = 0
OL = False
for line in ProdigalInput:
    if "CDS" in line :
        line = line.replace("<","")
        if 'complement' in line: #Reverse Compliment starts and stops to confirm to our deffinition
            #Switched to match Sense Strand

            loci = line.split()[1]
            loci = loci.partition('(')[-1].rpartition(')')[0]
            Start = int(loci.split('..')[0])
            Stop = int(loci.split('..')[1])
            r_stop = Genome_Size - Stop
            po = str(r_stop)+','+str(Stop)

            Combi_ORFs.update({ORF_NUM: str(Start) + ',' + str(Stop)})
            ORF_NUM += 1
        else:
            loci = line.split()[1]
            Start = int(loci.split('..')[0])
            Stop = int(loci.split('..')[1])
            po = str(Start)+','+str(Stop)

            Combi_ORFs.update({ORF_NUM: str(Start) + ',' + str(Stop)})
            ORF_NUM += 1


#####################################Prokka
ProkkaInput = open('./PROKKA_K-12.gff', mode='rb')
ProkkaORFs = collections.OrderedDict()
ProkkaORFsFiltered = collections.OrderedDict()
Prev_Stop = 0
OL = False
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

            po = str(r_start)+','+str(Stop)

        elif '+' in line.split()[6]:
            po = str(Start)+','+str(Stop)

        Combi_ORFs.update({ORF_NUM:str(Start) + ',' + str(Stop)})
        ORF_NUM += 1

#########################################~TransDecoder
TransDecoderInput = open('./TransDecoder_K-12.txt', mode='rb')
TransDecoderORFs = collections.OrderedDict()
TransDecoderORFsFiltered = collections.OrderedDict()
Prev_Stop = 0
OL = False
for line in TransDecoderInput:
    if "Chromosome	transdecoder	CDS" in line :
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

        elif '+' in line.split('\t')[6]:
            po = str(Start)+','+str(Stop)



        Combi_ORFs.update({ORF_NUM:str(str(Start) + ',' + str(Stop))})
        ORF_NUM += 1

Output = orfComparison(Combi_ORFs, Genes, Genome)
print Output
Out = open('./All_FP.txt', mode='wb')
for i, x in Combi_ORFs:
    start = x.split(',')[0]
    stop = x.split(',')[1]
    Out.write(start+" + "+stop)






















































