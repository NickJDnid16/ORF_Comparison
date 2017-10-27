import collections




ConsensusORFs = collections.OrderedDict()
NON_ConsensusORFs = collections.OrderedDict()
PROKKA_ORFS = collections.OrderedDict()
GLIMMER_ORFs = collections.OrderedDict()
GeneMark_S_ORFs = collections.OrderedDict()
GeneMark_HA_ORFs = collections.OrderedDict()
MetaGeneMark_ORFs = collections.OrderedDict()
MetaGeneAnnotator_ORFs = collections.OrderedDict()
FragGeneScan_ORFs = collections.OrderedDict()
Prodigal_ORFs = collections.OrderedDict()


PROKKA_Input = open('./Input_Data/Myco/PROKKA_Myco.gff', mode='rb')
Keyi = 0
for line in PROKKA_Input:
    if "Chromosome	Prodigal:2.6	CDS" in line :
        Start = int(line.split()[3])
        Stop = int(line.split()[4])
        if '-' in line.split()[6]:
            po = str(Start)+','+str(Stop)+ ',-'
            PROKKA_ORFS.update({Keyi:po})
        elif '+' in line.split()[6]:
            po = str(Start)+','+str(Stop)+ ',+'
            PROKKA_ORFS.update({Keyi:po})
        Keyi +=1



GLIMMER_Input = open('./Input_Data/Myco/GLIMMER_Myco.txt', mode='rb')

for line in GLIMMER_Input:
    if "orf0" in line :
        Start = int(line.split()[1])
        Stop = int(line.split()[2])
        if '-' in line.split()[3]:
            po = str(Stop)+','+str(Start)+ ',-'
            GLIMMER_ORFs.update({Keyi:po})
        elif '+' in line.split()[3]:
            po = str(Start)+','+str(Stop)+ ',+'
            GLIMMER_ORFs.update({Keyi:po})
        Keyi += 1


GeneMark_S_Input = open('./Input_Data/Myco/GeneMark_S_Myco.txt', mode='rb')

for line in GeneMark_S_Input:
    if "Chromosome" in line:
        Start = int(line.split('\t')[3])
        Stop = int(line.split('\t')[4])
        if '-' in line.split('\t')[6]:
            po = str(Start)+','+str(Stop) + ',-'
            GeneMark_S_ORFs.update({Keyi:po})
        elif '+' in line.split('\t')[6]:
            po = str(Start)+','+str(Stop) + ',+'
            GeneMark_S_ORFs.update({Keyi:po})
        Keyi += 1




GeneMark_HA_Input = open('./Input_Data/Myco/GeneMark_HA_Myco.txt', mode='rb')
for line in GeneMark_HA_Input:
    if "Chromosome" in line:
        Start = int(line.split('\t')[3])
        Stop = int(line.split('\t')[4])
        if '-' in line.split('\t')[6]:
            po = str(Start) + ',' + str(Stop) + ',-'
            GeneMark_HA_ORFs.update({Keyi:po})
        elif '+' in line.split('\t')[6]:
            po = str(Start)+','+str(Stop) + ',+'
            GeneMark_HA_ORFs.update({Keyi:po})
        Keyi += 1


MetaGeneMark_Input = open('./Input_Data/Myco/MetaGeneMark_Myco.txt', mode='rb')
for line in MetaGeneMark_Input:
    if "Chromosome" in line:
        Start = int(line.split('\t')[3])
        Stop = int(line.split('\t')[4])
        if '-' in line.split('\t')[6]:
            po = str(Start) + ',' + str(Stop) + ',-'
            MetaGeneMark_ORFs.update({Keyi:po})
        elif '+' in line.split('\t')[6]:
            po = str(Start) + ',' + str(Stop) + ',+'
            MetaGeneMark_ORFs.update({Keyi:po})
        Keyi +=1


MetaGeneAnnotator_Input = open('./Input_Data/Myco/MetaGeneAnnotator_Myco.txt', mode='rb')
for line in MetaGeneAnnotator_Input:
    if "gene_" in line:
        Start = int(line.split('\t')[1])
        Stop = int(line.split('\t')[2])
        if '-' in line.split('\t')[3]:
            po = str(Start) + ',' + str(Stop) + ',-'
            MetaGeneAnnotator_ORFs.update({Keyi:po})
        elif '+' in line.split('\t')[3]:
            po = str(Start) + ',' + str(Stop) + ',+'
            MetaGeneAnnotator_ORFs.update({Keyi:po})
        Keyi +=1


Prodigal_Input = open('./Input_Data/Myco/Prodigal_Myco.txt', mode='rb')
for line in Prodigal_Input:
    if "CDS" in line:
        line = line.replace("<", "")
        line = line.replace(">", "")
        if 'complement' in line:  # Reverse Compliment starts and stops to confirm to our deffinition
            # Switched to match Sense Strand
            loci = line.split()[1]
            loci = loci.partition('(')[-1].rpartition(')')[0]
            Start = int(loci.split('..')[0])
            Stop = int(loci.split('..')[1])
            po = str(Start) + ',' + str(Stop) + ',-'
            Prodigal_ORFs.update({Keyi:po})
        else:
            loci = line.split()[1]
            Start = int(loci.split('..')[0])
            Stop = int(loci.split('..')[1])
            po = str(Start) + ',' + str(Stop) + ',+'
            Prodigal_ORFs.update({Keyi:po})
        Keyi += 1

FragGeneScan_Input = open('./Input_Data/Myco/FragGeneScan_Myco.txt', mode='rb')
for line in FragGeneScan_Input:
    if "Chromosome	FGS	CDS" in line:
        Start = int(line.split()[3])
        Stop = int(line.split()[4])
        if '-' in line.split()[6]:  # Reverse Compliment starts and stops to confirm to our deffinition
            po = str(Start) + ',' + str(Stop) + ',-'
            FragGeneScan_ORFs.update({Keyi: po})
        elif '+' in line.split()[6]:
            po = str(Start) + ',' + str(Stop) + ',+'
            FragGeneScan_ORFs.update({Keyi: po})
        Keyi +=1







#########################################

# d  = PROKKA_ORFS.viewvalues() & Prodigal_ORFs.viewvalues()
#
# print d
PROKKA_Set = set(PROKKA_ORFS.values())
PRODIGAL_Set = set(Prodigal_ORFs.values())
GLIMMER_Set = set(GLIMMER_ORFs.values())
GeneMark_S_Set = set(GeneMark_S_ORFs.values())
GeneMark_HA_Set = set(GeneMark_HA_ORFs.values())
MetaGeneAnnotator_Set = set(MetaGeneAnnotator_ORFs.values())
MetaGeneMark_Set = set(MetaGeneMark_ORFs.values())
FragGeneScan_Set = set(FragGeneScan_ORFs.values())
#intersection = list(PROKKA_Set & PRODIGAL_Set & GLIMMER_Set & GeneMark_S_Set & GeneMark_HA_Set & MetaGeneAnnotator_Set & MetaGeneMark_Set & FragGeneScan_Set)
internot = list(PROKKA_Set ^ PRODIGAL_Set ^ GLIMMER_Set ^ GeneMark_S_Set ^ GeneMark_HA_Set ^ MetaGeneAnnotator_Set ^ MetaGeneMark_Set ^ FragGeneScan_Set)


intersection = list(PROKKA_Set & PRODIGAL_Set  )

Consensus_Out = open('./Direct_Consensus_Comparison/Consensus_ORFs.txt', mode='wb')
NON_Consensus_Out = open('./Non_Consensus_ORFs.txt', mode='wb')

for po in intersection:
    Consensus_Out.write(po+'\n')
for po in internot:
    NON_Consensus_Out.write(po+'\n')

print len(intersection)
print len(internot)




























#for key, po in intersection:

Near_Intersection = []

############################ Number of non-consensus ORFS per tool - And Near Consensus
# NP_PROKKA = 0
# I_PROKKA = 0
# Near_PROKKA = 0
# Filtered_PROKKA = []
# for po in PROKKA_ORFS.itervalues():
#     if po in intersection:
#         I_PROKKA +=1
#     elif po not in intersection:
#         I_PROKKA +=1
#
#         Start = int(po.split(',')[0])
#         Stop = int(po.split(',')[1])
#         Di = str(po.split(',')[2])
#         for n_po in internot:
#             n_Start = int(n_po.split(',')[0])
#             n_Stop = int(n_po.split(',')[1])
#             n_Di = str(n_po.split(',')[2])
#             if Di in n_Di:
#                 Start_Diff = abs(Start-n_Start)
#                 Stop_Diff= abs(Stop-n_Stop)
#                 if Start_Diff <= 10 and Stop_Diff <= 10:
#                     Near_PROKKA =+ 1
#                     Near_Intersection.append(po)




#
#
# #print I_PROKKA
# print NP_PROKKA
# ###########################
# NP_Prodigal = 0
# I_Prodigal = 0
# Filtered_Prodigal = []
# for po in Prodigal_ORFs.itervalues():
#     if po not in intersection:
#         NP_Prodigal +=1
#     elif po in intersection:
#         I_Prodigal +=1
# #print I_Prodigal
# print NP_Prodigal
# ###########################
# NP_GLIMMER = 0
# I_GLIMMER = 0
# Filtered_GLIMMER = []
# for po in GLIMMER_ORFs.itervalues():
#     if po not in intersection:
#         NP_GLIMMER +=1
#     elif po in intersection:
#         I_GLIMMER +=1
# #print I_GLIMMER
# print NP_GLIMMER
# ###########################
# NP_GeneMark_S = 0
# I_GeneMark_S = 0
# Filtered_GeneMark_S = []
# for po in GeneMark_S_ORFs.itervalues():
#     if po not in intersection:
#         NP_GeneMark_S +=1
#     elif po in intersection:
#         I_GeneMark_S +=1
# #print I_GeneMark_S
# print NP_GeneMark_S
# ###########################
# NP_GeneMark_HA = 0
# I_GeneMark_HA = 0
# Filtered_GeneMark_HA = []
# for po in GeneMark_HA_ORFs.itervalues():
#     if po not in intersection:
#         NP_GeneMark_HA +=1
#     elif po in intersection:
#         I_GeneMark_HA +=1
# #print I_GeneMark_HA
# print NP_GeneMark_HA
# ##########################
# NP_MetaGeneAnnotator = 0
# I_MetaGeneAnnotator = 0
# Filtered_MetaGeneAnnotator = []
# for po in MetaGeneAnnotator_ORFs.itervalues():
#     if po not in intersection:
#         NP_MetaGeneAnnotator +=1
#     elif po in intersection:
#         I_MetaGeneAnnotator +=1
# #print I_MetaGeneAnnotator
# print NP_MetaGeneAnnotator
# ##########################
# NP_MetaGeneMark = 0
# I_MetaGeneMark = 0
# Filtered_MetaGeneMark = []
# for po in MetaGeneMark_ORFs.itervalues():
#     if po not in intersection:
#         NP_MetaGeneMark +=1
#     elif po in intersection:
#         I_MetaGeneMark +=1
# #print I_MetaGeneMark
# print NP_MetaGeneMark
# ##########################
# NP_FragGeneScan = 0
# I_FragGeneScan = 0
# Filtered_FragGeneScan = []
# for po in FragGeneScan_ORFs.itervalues():
#     if po not in intersection:
#         NP_FragGeneScan +=1
#     elif po in intersection:
#         I_FragGeneScan +=1
# #print I_FragGeneScan
# print NP_FragGeneScan
# ##########################
#





