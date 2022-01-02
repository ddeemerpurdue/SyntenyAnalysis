'''
'''
import sys

import syntenyLibrary as SL


def combineSummaryAndGFF(summary_file, genome, genecall_file, output):
    annotations = ['KeggAccession', 'KeggDescription', 'CogFxAccession',
                   'CogFxDescription', 'CogPathway', 'COG20Category']
    summary = SL.mineSummaryFile(summary_file, genome)
    genecalls = SL.mineGeneCalls(genecall_file, 'all')

    with open(output, 'w') as out:
        for gene in genecalls:
            summary_entry = summary[gene]
            genecall_entry = genecalls[gene]
            attr = f'contig={genecall_entry.contig}'
            for annotation, entry in zip(annotations, summary_entry):
                attr = f"{attr};{annotation}={entry.replace(' ', '_')}"
            attr = attr + '\n'
            contig = genecall_entry.contig
            source = 'Anvio'
            feature = 'gene'
            start = str(genecall_entry.start)
            end = str(genecall_entry.stop)
            score = '.'  # Can incorporate this
            strand = genecall_entry.strand
            frame = '0'
            full_line = '\t'.join(
                [contig, source, feature, start, end, score, strand, frame, attr])
            out.write(full_line)

    return 0


combineSummaryAndGFF('rawdata/All-PS-Summary.txt', '0_Acidaminococcaceae_Acidaminococcus_P',
                     'rawdata/GFF-Internal/0_Acidaminococcaceae_Acidaminococcus_P.gff',
                     'tmp.gff')
