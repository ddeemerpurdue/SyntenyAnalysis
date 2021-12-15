'''
'''
import sys


class SummaryFileEntry:
    def __init__(self, line):
        self.values = line.strip().split('\t')
        self.genome = self.values[3]
        self.gene = self.values[4]
        self.annotation_fields = [12, 13, 14, 15, 16, 19]
        if self.testAnnotationFieldsValid():
            self.annotations = self.parseLine()
        else:
            self.annotations = ['No_annotations']

    def testAnnotationFieldsValid(self):
        if len(self.values) > self.annotation_fields[-1]:
            return True
        return False

    def assertValidLineToParse(self):
        assert len(self.values) > 12, 'Not enough fields! in the summary file!'
        return 0

    def parseLine(self):
        return [self.values[v] for v in self.annotation_fields]


class OrthologFileEntry:
    def __init__(self, line):
        self.values = line.strip().split('\t')
        assert len(self.values) == 3, 'Invalid number of fields'
        self.a_orthologs = self.getAGenes()
        self.b_ortholog = self.getBGenes()

    def getAGenes(self):
        a_genes = [gene.split('|')[0].strip()
                   for gene in self.values[1].split(',')]
        return a_genes

    def getBGenes(self):
        b_genes = self.values[2].split(',')[0].split('|')[0].strip()
        return b_genes


class GeneCallEntry:
    def __init__(self, line):
        self.line = line
        if line == '+':
            self.gene, self.contig, self.start, self.stop = ['+'] * 4
        elif line == '-':
            self.gene, self.contig, self.start, self.stop = ['-'] * 4
        elif line:
            self.values = line.strip().split('\t')
            assert len(
                self.values) in [10, 11, 12], 'Invalid number of fields in GeneCallEntry'
            self.entry_list = self.values[0:4]
            self.gene, self.contig, self.start, self.stop = self.entry_list
            self.testValues()
        else:
            self.gene, self.contig, self.start, self.stop = [None] * 4

    def testValues(self):
        try:
            self.start, self.stop = int(self.start), int(self.stop)
        except ValueError:
            print(self.line)
            raise ValueError(
                f'One of start/stop values is not an integer: {self.start}:{self.stop}')
        if not self.gene.strip():
            raise ValueError(f'Line {self.values} contains no gene!')
        if not self.contig.strip():
            raise ValueError(f'Line {self.values} contains no contig!')
        return 0

    def setUpstreamEntry(self, line):
        if line is None:
            self.upstream_entry_list = ['-'] * 4
            self.upstream_gene, self.upstream_contig = ['-', '-']
            self.upstream_start, self.upstream_stop = ['-', '-']
        else:
            self.upstream_values = line.strip().split('\t')
            assert len(
                self.upstream_values) in [10, 11, 12], 'Invalid number of fields in GeneCallEntry'
            self.upstream_entry_list = self.upstream_values[0:4]
            self.upstream_gene, self.upstream_contig, self.upstream_start, self.upstream_stop = self.upstream_entry_list

    def setDownstreamEntry(self, line):
        if line is None:
            self.downstream_entry_list = ['+'] * 4
            self.downstream_gene, self.downstream_contig = ['+', '+']
            self.downstream_start, self.downstream_stop = ['+', '+']
        else:
            self.downstream_values = line.strip().split('\t')
            assert len(
                self.downstream_values) in [10, 11, 12], 'Invalid number of fields in GeneCallEntry'
            self.downstream_entry_list = self.downstream_values[0:4]
            self.downstream_gene, self.downstream_contig, self.downstream_start, self.downstream_stop = self.downstream_entry_list

    def sameContigAsPrevious(self):
        if self.contig == self.upstream_contig:
            return True
        else:
            return False
