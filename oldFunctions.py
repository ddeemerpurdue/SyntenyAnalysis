# 1
def mineSummaryFile(summary_file, genome_name):
    gene_to_annotation = {}
    count = 0
    with open(summary_file) as file:
        line = file.readline()
        while line:
            CurrentEntry = SC.SummaryFileEntry(line)
            if CurrentEntry.genome == genome_name:
                gene_to_annotation[CurrentEntry.gene] = CurrentEntry.annotations
                count += 1
            line = file.readline().strip()
    if count > 0:
        return gene_to_annotation
    else:
        raise ValueError(
            f'Genome: {genome_name} was not found in the summary file.')


# 2
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


# 3
def getGenomeNames(ortholog_file, override):
    if override:
        assert len(override) == 2, 'Invalid override format'
        return override[0], override[1]
    else:
        with open(ortholog_file) as f:
            genomes = f.readline().split('\t')[1:]
            assert len(
                genomes) == 2, 'Invalid field number in ortholog file when assigning genomes.'
            return genomes[0], genomes[1]


# 4
def grabValidGeneIds(gene_to_annotation, gff=False):
    # TODO Can remove genetoannotation
    assertError = 'Invalid number of fields in GeneCallEntry when grabbing valid gene ids'
    if not gff:
        valid_gene_ids = list(gene_to_annotation.keys())
    else:
        valid_gene_ids = []
        with open(gff) as file:
            line = file.readline()
            if 'GCA' in gff or 'GCF' in os.path.basename(gff):
                if '-' in os.path.basename(gff):  # TODO
                    line = file.readline()
            assert len(line.split('\t')) in [10, 11, 12], assertError
            while line:
                valid_gene_ids.append(line.split('\t')[0])
                line = file.readline()
    valid_gene_ids.sort()
    return valid_gene_ids


# 5
def testSummaryGeneIDsMatchOrthologFile(gene_to_annotation, orthologs_a):
    summary_geneIDs = list(gene_to_annotation.keys())
    ortholog_geneIDs = list(orthologs_a.keys())
    if allXElementsInY(ortholog_geneIDs, summary_geneIDs):
        return True
    else:
        invald_ids = diffXElementsNotInY(ortholog_geneIDs, summary_geneIDs)
        invald_ids_print = '\n'.join(invald_ids)
        raise IndexError(
            f"Could not match the following ortholog gene ids:\n{invald_ids_print}")


# 6
def allXElementsInY(x, y):
    if all(val in y for val in x):
        return True
    else:
        return False


# 7
def diffXElementsNotInY(x, y):
    x = [str(v) for v in x]
    y = [str(v) for v in y]
    return list(set(x) - set(y))


# 8
def testValidGeneIDsLongerEqualSynteny(valid_gene_ids, synteny):
    if len(valid_gene_ids) >= len(synteny):
        return True
    else:
        invald_ids = diffXElementsNotInY(synteny, valid_gene_ids)
        invald_ids_print = '\n'.join([str(v) for v in invald_ids])
        print(invald_ids_print)
        raise ValueError('More genes in synteny file than valid gene ids.')


# 9
def geneInValidGeneIDs(GeneCallEntry, valid_gene_ids):
    # TODO: Currently the traverseSynteny cannot have
    # valid gene ids == 'all'
    if valid_gene_ids == 'all' or GeneCallEntry.gene in valid_gene_ids:
        return True
    else:
        return False


# 10
def geneIsValid(gene, ignore, ortholog_dic, no_orthos):
    if geneInIgnore(gene, ignore):
        return False
    if geneHasOrtholog(gene, ortholog_dic, no_orthos):
        return True
    return False


# 11
def addCurrentGeneToList(GeneCallEntry, list_):
    testGeneCallEntry(GeneCallEntry)
    list_.append(GeneCallEntry.gene)
    return 0


# 12
def moveStream(synteny_dic, next_gene):
    return synteny_dic[next_gene]


# 13
def moveBothStream(synteny, CurrentGene, CurrentOrtholog, switch):
    testGeneCallEntry()
    x_synteny = moveStream(synteny[switch], CurrentGene.gene)
    y_synteny = moveStream(synteny[switchFlag(switch)], CurrentOrtholog.gene)
    return [x_synteny, y_synteny]


# 14
def appendMiscValues(loc_info, switch, Gene, Ortholog, SeedDirection):
    testGeneCallEntry(Gene)
    up_down = f"{Gene.upstream_gene}-{Gene.downstream_gene}"
    start_stop = f"{Gene.start}-{Gene.stop}"
    full_misc = f"{up_down}\t{start_stop}"

    up_down_orth = f"{Ortholog.upstream_gene}-{Ortholog.downstream_gene}"
    start_stop_orth = f"{Ortholog.start}-{Ortholog.stop}"
    full_misc_orth = f"{up_down_orth}\t{start_stop_orth}"

    if SeedDirection.direction_x == 'Downstream':
        loc_info[switch].append(full_misc)
        loc_info[switchFlag(switch)].append(full_misc_orth)
    elif SeedDirection.direction_x == 'Upstream':
        loc_info[switch].insert(0, full_misc)
        loc_info[switchFlag(switch)].insert(0, full_misc_orth)
    return 0


# 15
def testGeneInSynteny(gene, synteny):
    if gene in synteny:
        return True
    else:
        raise ValueError(f'Gene {gene} is not in synteny dictionary.')


# 16
''' Below, last part of return function for syntenyTracker
c = 0
all_written_x = set()
all_written_y = set()
for seed in synteny_dic.keys():
    c += len(synteny_dic[seed][0])
    for gene in synteny_dic[seed][0]:
        all_written_x.add(gene.gene)
    for gene in synteny_dic[seed][0]:
        all_written_y.add(gene.gene)
print(f"Syntenic loci found: {c}")
no_orthos.sort()
x_not_written = set(X_Genecalls.keys()) - all_written_x
y_not_written = set(Y_Genecalls.keys()) - all_written_y
print(f"Not written for a: {len(x_not_written)}")
print(f"Not written for b: {len(y_not_written)}")

a_not = {}
b_not = {}
for gene in x_not_written:
    a_not[gene] = X_Genecalls[gene]
for gene in y_not_written:
    b_not[gene] = Y_Genecalls[gene]

return [synteny_dic, a_not, b_not]
'''
