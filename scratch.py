import matplotlib.pyplot as plt
import networkx as nx

orthologs = 'rawdata/Orthologs-Acidaminococcaceae/Orthologues_0_Acidaminococcaceae_Acidaminococcus_P/0_Acidaminococcaceae_Acidaminococcus_P__v__0_Acidaminococcaceae_Acidaminococcus-P-GCA_000425045.1.tsv'


def readOrthologs(ortholog_file):
    graph = nx.DiGraph()
    with open(ortholog_file) as file:
        next(file)
        line = file.readline()
        while line:
            values = line.strip().split('\t')
            values[1] = values[1].split('|')[0]
            node = '-'.join(values[1:])
            print(node)
            line = file.readline()
    return 0


readOrthologs(orthologs)
