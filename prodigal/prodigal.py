import collections
from OrfFinderComparison import orfComparison

def prodigal(input_to_analyse,Genome,Genes):

    ProdigalORFs = collections.OrderedDict()
    Start_Codons = collections.OrderedDict()
    Stop_Codons = collections.OrderedDict()
    Genome_Size = len(Genome)

    with open(input_to_analyse,'r') as prodigal_input:
        for line in prodigal_input:
            if "CDS" in line :
                line = line.replace("<","")
                if 'complement' in line:
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

    return Output,Start_Precision,Stop_Precision
















