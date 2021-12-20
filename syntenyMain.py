import glob
import os

import argparse

import syntenyLibrary as SL
import syntenyTracker as ST


def getBins(files):
    mybins = []
    for file in files:
        bin_name = os.path.basename(file)
        if any([string in bin_name for string in ['GCA', 'GCF']]):
            pass
        else:
            mybins.append(bin_name.strip('Orthologues_'))
    return mybins


def grabBinB(ortholog_file):
    bin_b = os.path.basename(ortholog_file)
    return bin_b.split('__v__')[1].strip('.tsv')


def grabGffB(bin_b, family):
    if any([string in bin_b for string in ['GCA', 'GCF']]) and '-' in bin_b:
        gffB_file = glob.glob(f'rawdata/GFF-{family}/{bin_b}*.gff')
        assert len(gffB_file) == 1, 'Invalid gffB_file pattern matching!'
        gffB_file = gffB_file[0]
        Gb_flag = gffB_file
    else:
        gffB_file = f'rawdata/GFF-Internal/{bin_b}.gff'
        Gb_flag = False
    return gffB_file, Gb_flag


def makeDir(path):
    if not os.path.exists(path):
        os.makedirs(path)
    return 0


def writeSyntenyOutput(output, bin_a, bin_b, synteny_dic, summary):
    with open(output, 'w') as out:
        out.write(f"{bin_a}\t{bin_b}\n")  # Update header
        for seed in synteny_dic:
            for gene_a, gene_b, loc_a, loc_b in zip(synteny_dic[seed][0], synteny_dic[seed][1], synteny_dic[seed][2], synteny_dic[seed][3]):
                if gene_a in ['-', '+']:
                    vals = '\t'.join(['', '', '', '', '', ''])
                else:
                    vals = summary[gene_a]
                    vals = '\t'.join(vals)
                out.write(f"{gene_a}\t{gene_b}\t{loc_a}\t{loc_b}\t{vals}\n")
            out.write(f"X\tX\n")


def writeNoOrthologs(output,)


def main(family, summary_file):
    files = glob.glob(f'rawdata/Orthologs-{family}/*')
    mybins = getBins(files)

    print(f'Starting master analysis:')

    for bin_a in mybins:
        print(f"Analyzing data for bin: {bin_a}")
        orth_prefix = f'rawdata/Orthologs-{family}/Orthologues_{bin_a}/'
        print(f"\tOrtholog prefix: {orth_prefix}")
        for ortholog_file in glob.glob(f"{orth_prefix}/*"):
            print(f"\t\tMining ortho file: {os.path.basename(ortholog_file)}")
            bin_b = grabBinB(ortholog_file)
            print(f'\t\tAnalyzing {bin_a} vs {bin_b}')

            gffA_file = f'rawdata/GFF-Internal/{bin_a}.gff'
            gffB_file, Gb_flag = grabGffB(bin_b, family)

            print(f"\t\t\tGff a: {gffA_file}\n\t\t\tGff b: {gffB_file}")

            makeDir(f"Output/{family}")

            uniq_name = os.path.basename(ortholog_file)
            output = f'Output/{family}/{uniq_name}'
            print(f"\t\t\t\tWriting analysis results to file: {output}\n\n")

            summary = SL.mineSummaryFile(summary_file, bin_a)
            synteny_dic, no_orthos = ST.traverseSynteny(summary_file, gffA_file, gffB_file,
                                                        ortholog_file, Gb=Gb_flag)

            uniq_name = os.path.basename(ortholog_file)
            output = f'Output/{family}/{uniq_name}'
            writeSyntenyOutput(output, bin_a, bin_b, synteny_dic, summary)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Parser")
    parser.add_argument("-f", "--Family", help="Family",
                        required=True)
    parser.add_argument("-s", "--Summary", help="Summary file",
                        required=True)
    argument = parser.parse_args()
    main(argument.Family, argument.Summary)
