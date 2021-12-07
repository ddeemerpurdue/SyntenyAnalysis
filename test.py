'''

'''
from collections import namedtuple
from collections import OrderedDict
import sys

import syntenyLibrary as SL
import syntenyClasses as SC


def createAnnotationNamedTuple():
    return namedtuple("Annotations", "cog, kegg, tigrfam, pfam, figfam, kofam")


def grabAnnotationFields(fields):
    annotation_indices = [12, 14, 16, 18, 20, 28]
    return [fields[v] for v in annotation_indices]


def mineSummaryFile(summary_file, genome_name):
    count = 0
    genome_summary = {}
    with open(summary_file) as file:
        line = file.readline()  # Skip header
        line = file.readline()
        while line:
            CurrentEntry = SC.SummaryFileEntry(line)
            if CurrentEntry.genome == genome_name:
                count += 1
                genome_summary[CurrentEntry.gene] = CurrentEntry.annotations
            line = file.readline().strip()

    return genome_summary


# sfile = 'rawdata/Hot_Lakes_2021_UCC_ALORS-Clustering-WithAllReference_gene_clusters_summary.txt'
# a = mineSummaryFile(sfile, 'Algoriphagus__UCCA')

def mineOrthologFile(ortholog_file):
    # TODO - Add redundancy for multi-multi-hits
    ortholog_summary = {}
    with open(ortholog_file) as file:
        line = file.readline().strip()
        genome_a, genome_b = line.split('\t')[1:]
        line = file.readline().strip()
        while line:
            CurrentOrthologEntry = SC.OrthologFileEntry(line)
            for a_ortholog in CurrentOrthologEntry.a_orthologs:
                ortholog_summary[a_ortholog] = CurrentOrthologEntry.b_ortholog
            line = file.readline().strip()
    return ortholog_summary, genome_a, genome_b


# ortholog_file = 'rawdata/Algoriphagus__UCCA__v__Algoriphagus__UCCL.tsv'

# mineOrthologFile(ortholog_file)


def mineGeneCalls(genecall_file, valid_gene_ids):
    genecall_summary = {}
    Upstream, Current, Downstream = SL.initNamedTuples()
    with open(genecall_file) as file:
        current_line = file.readline()  # Skip header
        current_line = file.readline().strip()
        upstream_line = None
        UpstreamGeneEntry = SC.GeneCallEntry(None)
        while current_line:
            CurrentGeneEntry = SC.GeneCallEntry(current_line)

            if CurrentGeneEntry.gene in valid_gene_ids:
                if not SL.sameContigAsUpstream(CurrentGeneEntry.contig, UpstreamGeneEntry.contig):
                    CurrentGeneEntry.setUpstreamEntry(None)
                    UpstreamGeneEntry.setDownstreamEntry(None)
                else:
                    CurrentGeneEntry.setUpstreamEntry(upstream_line)
                    UpstreamGeneEntry.setDownstreamEntry(current_line)

                if UpstreamGeneEntry.gene not in ['-', '+']:
                    genecall_summary[UpstreamGeneEntry.gene] = UpstreamGeneEntry
                UpstreamGeneEntry = CurrentGeneEntry
                current_line = upstream_line

            current_line = file.readline().strip()

        UpstreamGeneEntry.setDownstreamEntry(None)
        genecall_summary[UpstreamGeneEntry]

    return genecall_summary
