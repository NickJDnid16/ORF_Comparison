import collections
import click
import csv

@click.command()
@click.option('--tool', default='GFF', help='Which tool/format would you like to analyse?')
@click.option('--tool_model', default='', help='Model of tool if needed.')
@click.option('--input_to_analyse', help='Location of file to compare.')
@click.option('--genome_to_compare', help='Which genome to analyse?')


def comparator(tool,tool_model,input_to_analyse,genome_to_compare):


    Genome = ""
    with open('genomes/'+genome_to_compare+'.txt', 'r') as genome:
        for line in genome:
            line = line.replace("\n","")
            if ">" not in line:
                Genome += str(line)


##############################################
    Genes = collections.OrderedDict()

    count = 0
    with open('genomes/'+genome_to_compare+'.gff','r') as genome_gff:
        for line in genome_gff:
            if "Chromosome	ena	CDS" in line:
                count += 1
                Start = int(line.split('\t')[3])
                Stop = int(line.split('\t')[4])
                Gene = str(Start)+','+str(Stop)
                Genes.update({count:Gene})
    print(len(Genes))
    #
# ##################################
    from importlib import import_module
    tool_compare = import_module(tool+'.'+tool)
    tool_compare = getattr(tool_compare,tool)


    outputDescription,output, start_precision, stop_precision = tool_compare(input_to_analyse,Genome,Genes)

    with open(tool+'/'+genome_to_compare+'_'+tool+'_'+tool_model+'.csv', 'w') as out_file:
        tool_out = csv.writer(out_file, delimiter=',')
        tool_out.writerow(outputDescription)
        tool_out.writerow(output)


















if __name__ == "__main__":
    comparator()








