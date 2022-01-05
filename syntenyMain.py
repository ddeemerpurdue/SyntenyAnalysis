import glob
import os
import sys

import argparse

import syntenyLibrary as SL
import syntenyTracker as ST


def makeDir(path):
    if not os.path.exists(path):
        os.makedirs(path)
    return 0


def writeSyntenyOutput(output, bin_a, bin_b, synteny):
    synteny_dic, a_not, b_not = synteny
    annotations = ['Transfer_RNAs', 'COG20_CATEGORY', 'COG20_PATHWAY',
                   'COG20_FUNCTION', 'KEGG_Module', 'KEGG_Class', 'KOfam',
                   'FigFams', 'Pfam', 'TIGRFAM']
    header = '\t'.join(annotations)
    header += f'\tContig\txUpStreamDownstream\txStartStop\t{bin_a}\t'
    header += f'{bin_b}\tyStartStop\tyUpstreamDownstream\tContig\t'
    annotations.reverse()
    header += '\t'.join(annotations) + '\n'
    with open(output, 'w') as out:
        out.write(f"{header}\n")
        for seed in synteny_dic:
            for Gene, Ortho in zip(synteny_dic[seed][0], synteny_dic[seed][1]):
                writeline = SL.formatWriteLine(Gene=Gene, Ortholog=Ortho)
                out.write(writeline)
            split = '\t'.join(['X'] * 28) + '\n'
            out.write(split)
        out.write('\n')
        # Writing non-orthologs
        for seed in a_not:
            Gene = a_not[seed]
            writeline = SL.formatWriteLine(Gene=Gene, Ortholog=False)
            out.write(writeline)
        for seed in b_not:
            Ortholog = b_not[seed]
            writeline = SL.formatWriteLine(Gene=False, Ortholog=Ortholog)
            out.write(writeline)


def parseOrthologDirectory(directory):
    for ortholog_file in glob.glob(f"{directory}/rawdata/Orthologs/*.tsv"):
        values = ortholog_file.split('__')
        genome_a = f"{os.path.basename(values[0])}.gff3"
        genome_a_gff = os.path.join(f"{directory}/rawdata/GFF3/", genome_a)
        genome_b = f"{values[2].strip('.tsv')}.gff3"
        genome_b_gff = os.path.join(f"{directory}/rawdata/GFF3/", genome_b)

        makeDir(f"{directory}/output")

        uniq_name = os.path.basename(ortholog_file).replace('.tsv', '.txt')
        output = f"{directory}/output/{uniq_name}"
        print(
            f'Running synteny analyis using:\n{genome_a_gff}\n{genome_a_gff}\n{ortholog_file}\n')
        synteny = ST.traverseSynteny(genome_a_gff, genome_b_gff, ortholog_file)
        writeSyntenyOutput(output, genome_a, genome_b, synteny)


parseOrthologDirectory(sys.argv[1])
