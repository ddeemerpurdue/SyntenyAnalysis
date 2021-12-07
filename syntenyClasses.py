'''
'''


class SummaryFileEntry:
    def __init__(self, line):
        self.values = line.strip().split('\t')
        self.annotation_fields = [12, 14, 16, 18, 20, 28]
        self.annotations = self.parseLine()
        self.genome = self.values[3]
        self.gene = self.values[4]

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
        return [gene for gene in self.values[1].split(',')]

    def getBGenes(self):
        return self.values[2].split(',')[0]


class GeneCallEntry:
    def __init__(self, line):
        self.values = line.strip().split('\t')
        assert len(self.values) == 10, 'Invalid number of fields in GeneCallEntry'
        self.entry_list = self.values[0:4]
        self.gene, self.contig, self.start, self.stop = self.entry_list

    def setUpstreamEntry(self, line):
        if line is None:
            self.upstream_entry_list = ['-'] * 4
            self.upstream_gene, self.upstream_contig = ['-', '-']
            self.upstream_start, self.upstream_stop = ['-', '-']
        else:
            self.upstream_values = line.stip().split('\t')
            assert len(
                self.values) == 10, 'Invalid number of fields in GeneCallEntry'
            self.upstream_entry_list = self.values[0:4]
            self.upstream_gene, self.upstream_contig, self.upstream_start, self.upstream_stop = self.entry_list

    def setDownstreamEntry(self, line):
        if line is None:
            self.downstream_entry_list = ['+'] * 4
            self.downstream_gene, self.downstream_contig = ['+', '+']
            self.downstream_start, self.downstream_stop = ['+', '+']
        else:
            self.downstream_values = line.stip().split('\t')
            assert len(
                self.values) == 10, 'Invalid number of fields in GeneCallEntry'
            self.downstream_entry_list = self.values[0:4]
            self.downstream_gene, self.downstream_contig, self.downstream_start, self.downstream_stop = self.entry_list

    def sameContigAsPrevious(self):
        if self.contig == self.upstream_gene:
            return True
        else:
            return False
