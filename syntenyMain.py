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
    # annotations = ['FigFams', 'FigFams-ACC', 'FigFams-eval', 'KEGG_Module', 'KEGG_Module-ACC',
    #                'KEGG_Module-eval', 'COG20_PATHWAY', 'COG20_PATHWAY-ACC', 'COG20_PATHWAY-eval',
    #                'TIGRFAM', 'TIGRFAM-ACC', 'TIGRFAM-eval', 'KOfam', 'KOfam-ACC', 'KOfam-eval', 'KEGG_Class',
    #                'KEGG_Class-ACC', 'KEGG_Class-eval', 'Transfer_RNAs', 'Transfer_RNAs-ACC', 'Transfer_RNAs-eval',
    #                'COG20_FUNCTION', 'COG20_FUNCTION-ACC', 'COG20_FUNCTION-eval', 'Pfam', 'Pfam-ACC', 'Pfam-eval',
    #                'COG20_CATEGORY', 'COG20_CATEGORY-ACC', 'COG20_CATEGORY-eval']
    annotations = ['ID', 'product']
    header = '\t'.join(annotations)
    header += f'\tContig\txUpStreamDownstream\txStartStop\t{bin_a}\t'
    header += f'{bin_b}\tyStartStop\tyUpstreamDownstream\tContig\t'
    annotations.reverse()
    header += '\t'.join(annotations) + '\n'
    print(header)
    with open(output, 'w') as out:
        out.write(f"{header}")
        for seed in synteny_dic:
            for Gene, Ortho in zip(synteny_dic[seed][0], synteny_dic[seed][1]):
                writeline = SL.formatWriteLine(Gene=Gene, Ortholog=Ortho)
                out.write(writeline)
            split = '\t'.join(['x'] * 12) + '\n'
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


def parseOrthologDirectory(orthos, gff3):
    print('Starting analysis...')
    for cnt, ortholog_file in enumerate(glob.glob(f"{orthos}/*.tsv")):
        if cnt % 250 == 0:
            print(f'Processed {cnt} samples')
        ortholog_file_base = os.path.basename(ortholog_file)
        values = os.path.splitext(ortholog_file_base)[0].split('__v__')
        genome_a = f"{values[0]}.gff"
        genome_a_gff = os.path.join(gff3, genome_a)
        genome_b = f"{values[1]}.gff"
        genome_b_gff = os.path.join(gff3, genome_b)
        print(f'OrthoF: {ortholog_file}')
        print(f'Gff3: {genome_a_gff}, {genome_b_gff}')

        makeDir(f"Output/")

        uniq_name = os.path.basename(ortholog_file).replace('.tsv', '.txt')
        output = f"Output/{uniq_name}"

        # print(
        #    f'Running synteny analyis using:\n{genome_a_gff}\n{genome_a_gff}\n{ortholog_file}\n')
        synteny = ST.traverseSynteny(genome_a_gff, genome_b_gff, ortholog_file)
        print(synteny[0])
        print(f'Genomes:\n{genome_a}: {genome_b}')
        writeSyntenyOutput(output, genome_a, genome_b, synteny)


parseOrthologDirectory(sys.argv[1], sys.argv[2])
