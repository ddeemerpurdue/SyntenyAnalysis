'''
'''
from collections import namedtuple
import syntenyTracker as ST


def grabValidGeneIds(summary_file, genome_a, genome_b):
    valid_gene_ids_a = list(ST.mineSummaryFile(summary_file, genome_a).keys())
    valid_gene_ids_b = list(ST.mineSummaryFile(summary_file, genome_b).keys())
    return valid_gene_ids_a, valid_gene_ids_b


class Synteny:

    def __init__(self, gff_file, valid_gene_ids):
        self.gff_file = gff_file
        self.valid_gene_ids = valid_gene_ids

    def test(self):
        if self.gff_file.endswith('.txt'):
            return True
        else:
            return True

    @staticmethod
    def initNamedTuples():
        Upstream = namedtuple("Upstream", "gene, contig, start, stop")
        Current = namedtuple("Current", "gene, contig, start, stop")
        Downstream = namedtuple("Downstream", "gene, contig, start, stop")
        return Upstream, Current, Downstream

    @staticmethod
    def setStartingValues():
        return [True, '', '', '', '']

    @staticmethod
    def sameContigAsPrevious(current, previous):
        if current == previous:
            return True
        return False

    @staticmethod
    def setCurrentValues(fields):
        return [fields[val] for val in [0, 1, 2, 3]]

    def readGenecallFile(self):
        genecall_summary = {}
        Upstream, Current, Downstream = self.initNamedTuples()
        with open(self.gff_file) as file:
            line = file.readline()
            line = file.readline().strip()
            first, prev_contig, prev_gene, prev_start, prev_stop = self.setStartingValues()
            while line:
                fields = line.split('\t')
                current_gene, current_contig, start, stop = self.setCurrentValues(
                    fields)
                if current_gene in self.valid_gene_ids:
                    if not self.sameContigAsPrevious(current_contig, prev_contig):
                        upstream_gene, upstream_contig, prev_start, prev_stop = [
                            '-'] * 4
                        prev_downstream_gene = '+'
                    else:
                        upstream_gene = prev_gene
                        upstream_contig = prev_contig
                        prev_downstream_gene = current_gene

                    # Add information to the dictionary
                    genecall_summary.setdefault(current_gene, [])
                    genecall_summary[current_gene].append(
                        Upstream(upstream_gene, upstream_contig, prev_start, prev_stop))
                    genecall_summary[current_gene].append(
                        Current(current_gene, current_contig, start, stop))
                    if first:
                        first = False
                    else:
                        if not self.sameContigAsPrevious(current_contig, prev_contig):
                            genecall_summary[prev_gene].append(
                                Downstream(prev_downstream_gene, '+', '+', '+'))
                        else:
                            genecall_summary[prev_gene].append(Downstream(
                                prev_downstream_gene, current_contig, start, stop))

                    prev_gene, prev_contig, prev_start, prev_stop = current_gene, current_contig, start, stop
                else:
                    pass
                line = file.readline().strip()
            genecall_summary[prev_gene].append(Downstream('+', '+', '+', '+'))

        return genecall_summary


summary_file = 'rawdata/Hot_Lakes_2021_UCC_ALORS-Clustering-WithAllReference_gene_clusters_summary.txt'
gffA = 'rawdata/A-genecalls'
valid_gene_ids = grabValidGeneIds(
    summary_file, 'Algoriphagus__UCCA', 'Algoriphagus__UCCL')

A_Synteny = Synteny(gffA, valid_gene_ids)
















# SL.appendBothIgnore(ignore, switch, next_.gene, current_ortholog.gene)
# SL.appendValues(values, switch, next_.gene, current_ortholog.gene, seed_direction)



# SL.recordSyntenyInStone(synteny_dic, seed[1], values)
# break
