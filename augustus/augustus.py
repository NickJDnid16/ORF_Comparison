import collections
from DNA_Reverse_Compliment import revCompIterative
from OrfFinderComparison import orfComparison

def augustus(input_to_analyse,Genome,Genes):
##############################################

    AugustusORFs = collections.OrderedDict()




    ##################################
    Start_Codons = collections.OrderedDict()
    Stop_Codons = collections.OrderedDict()
    Prev_Stop = 0
    OL = False
    Genome_Size = len(Genome)
    Genome_rev = revCompIterative(Genome)


    with open(input_to_analyse,'r') as augustus_input:
        for line in augustus_input:
            if "Chromosome	AUGUSTUS	CDS" in line :
                Start = int(line.split('\t')[3])
                Stop = int(line.split('\t')[4])
                if Start < Prev_Stop and OL == False:
                    print ("Overlapping")
                    OL = True
                Prev_Stop = Stop
                reverse = '-'
                forward = '+'
                if '-' in line.split('\t')[6]: #Reverse Compliment starts and stops to confirm to our deffinition
                    #Switched to match Sense Strand
                    r_start = Genome_Size - Stop
                    r_stop = Genome_Size - Start
                    po = str(r_start)+','+str(Stop)
                    if po  in Start_Codons.keys():
                        print("Here")
                    Start_Codons.update({po:reverse})
                    po = str(Start) + ',' + str(r_stop)
                    if po  in Stop_Codons.keys():
                        print("Here")
                    Stop_Codons.update({po: reverse})

                elif '+' in line.split('\t')[6]:
                    po = str(Start)+','+str(Stop)
                    if po  in Start_Codons.keys():
                        print("Here")
                    Start_Codons.update({po:forward})
                    if po  in Stop_Codons.keys():
                        print("Here")
                    Stop_Codons.update({po: forward})

                AugustusORFs.update({Start:Stop})



    print (len(Stop_Codons))
    print (len(Start_Codons))

    outputDescription,output,Start_Precision,Stop_Precision = orfComparison(Genes,AugustusORFs,Start_Codons,Stop_Codons,Genome)

    return outputDescription,output,Start_Precision,Stop_Precision
































