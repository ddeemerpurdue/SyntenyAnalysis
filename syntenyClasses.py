'''
'''


class SummaryFileEntry:
    def __init__(self, line):
        self.values = line.strip().split('\t')
        self.cluster = self.values[1]
        self.genome = self.values[3]
        self.gene = self.values[4]
        self.annotation_fields = [12, 13, 14, 15, 16, 19]
        if self.testAnnotationFieldsValid():
            self.annotations = self.parseLine()
        else:
            self.annotations = self.parseIncompleteLine()

    def testAnnotationFieldsValid(self):
        if len(self.values) > self.annotation_fields[-1]:
            return True
        return False

    def assertValidLineToParse(self):
        assert len(self.values) > 12, 'Not enough fields! in the summary file!'
        return 0

    def parseLine(self):
        return [self.values[v] for v in self.annotation_fields]

    def parseIncompleteLine(self):
        annotations = []
        for field in self.annotation_fields:
            try:
                annotations.append(self.values[field])
            except IndexError:
                annotations.append('No_annotation')
        return annotations


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
                self.values) == 9, 'Invalid number of fields in GeneCallEntry'
            self.contig = self.values[0]
            self.type = self.values[2]
            self.start = self.values[3]
            self.stop = self.values[4]
            self.strand = self.values[6]
            self.attributes = self.values[8]
            self.gene = self.setSelfGene(self.attributes)
            self.testValues()
        else:
            self.gene, self.contig, self.start, self.stop = [None] * 4

    def setSelfGene(self, attributes):
        values = attributes.split(';')
        for value_pair in values:
            values = value_pair.split('=')
            if values[0] == 'Locus':
                return values[1]
        raise ValueError(
            f'Locus field not found in attribute line: {attributes}')

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
                self.upstream_values) == 9, 'Invalid number of fields in GeneCallEntry'
            self.upstream_contig = self.upstream_values[0]
            self.upstream_type = self.upstream_values[2]
            self.upstream_start = self.upstream_values[3]
            self.upstream_stop = self.upstream_values[4]
            self.upstream_strand = self.upstream_values[6]
            self.upstream_attributes = self.upstream_values[8]
            self.upstream_gene = self.setSelfGene(self.upstream_attributes)

    def setDownstreamEntry(self, line):
        if line is None:
            self.downstream_entry_list = ['+'] * 4
            self.downstream_gene, self.downstream_contig = ['+', '+']
            self.downstream_start, self.downstream_stop = ['+', '+']
        else:
            self.downstream_values = line.strip().split('\t')
            assert len(
                self.downstream_values) == 9, 'Invalid number of fields in GeneCallEntry'
            self.downstream_contig = self.downstream_values[0]
            self.downstream_type = self.downstream_values[2]
            self.downstream_start = self.downstream_values[3]
            self.downstream_stop = self.downstream_values[4]
            self.downstream_strand = self.downstream_values[6]
            self.downstream_attributes = self.downstream_values[8]
            self.downstream_gene = self.setSelfGene(self.downstream_attributes)

    def sameContigAsPrevious(self):
        if self.contig == self.upstream_contig:
            return True
        else:
            return False


class SeedDirection:
    def __init__(self):
        self.seed_direction = 'Forward'
        self.direction_x = 'Downstream'
        self.direction_y = 'Upstream'

    def changeSeedAndXDirection(self):
        self.seed_direction = 'Reverse'
        self.direction_x = 'Upstream'
        return 0

    def restartFromSeed(self):
        if self.seed_direction == 'Forward':
            self.changeSeedAndXDirection()
            return True
        return False

    def XtoYandYtoNeither(self):
        self.direction_x = self.direction_y
        self.direction_y = 'Neither'


class CurrentRun:
    def __init__(self):
        pass
