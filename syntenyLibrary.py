'''

'''
from collections import namedtuple
import os
import sys

import syntenyClasses as SC


'''
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
'''


def mineOrthologFile(ortholog_file):
    # TODO - Add redundancy for multi-multi-hits
    orthologs = namedtuple("Orthologs", "X, Y")
    ortholog_summary_F = {}
    with open(ortholog_file) as file:
        next(file)
        line = file.readline().strip()
        while line:
            CurrentOrthologEntry = SC.OrthologFileEntry(line)
            for a_ortholog in CurrentOrthologEntry.a_orthologs:
                ortholog_summary_F[a_ortholog] = CurrentOrthologEntry.b_ortholog
            line = file.readline().strip()
    ortholog_summary_R = {v: k for k, v in ortholog_summary_F.items()}
    orthologs.X, orthologs.Y = ortholog_summary_F, ortholog_summary_R
    return orthologs


def mineGff3File(gff3_file):
    gene_loci = {}
    with open(gff3_file) as file:  # TODO: Below
        # if 'GCA' in gff3_file or 'GCF' in os.path.basename(gff3_file):
        #     if '-' in os.path.basename(gff3_file):
        #         next(file)  # Skip header :-)
        next(file)
        current_line = file.readline().strip()
        # print(current_line)
        upstream_line = None
        UpstreamGeneEntry = SC.GeneCallEntry('-')
        while current_line:
            CurrentGeneEntry = SC.GeneCallEntry(current_line)

            # if not geneInValidGeneIDs(CurrentGeneEntry, valid_gene_ids):
            #     upstream_line = current_line
            #     current_line = file.readline().strip()
            #     continue

            if sameContigAsUpstream(CurrentGeneEntry.contig, UpstreamGeneEntry.contig):
                CurrentGeneEntry.setUpstreamEntry(upstream_line)
                UpstreamGeneEntry.setDownstreamEntry(current_line)
            else:
                CurrentGeneEntry.setUpstreamEntry(None)
                UpstreamGeneEntry.setDownstreamEntry(None)

            if UpstreamGeneEntry.gene not in ['-', '+', 'x', None]:
                gene_loci[UpstreamGeneEntry.gene] = UpstreamGeneEntry

            UpstreamGeneEntry = CurrentGeneEntry
            upstream_line = current_line
            current_line = file.readline().strip()

        UpstreamGeneEntry.setDownstreamEntry(None)

        gene_loci[UpstreamGeneEntry.gene] = UpstreamGeneEntry

    return gene_loci


def printNDicValues(dict_, N):
    for cnt, (k, v) in enumerate(dict_.items()):
        if cnt < N:
            print(f"{k} - {v}")
        else:
            break


def sameContigAsUpstream(current, previous):
    if current == previous:
        return True
    return False


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


'''
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
'''


def geneInValidGeneIDs(GeneCallEntry, valid_gene_ids):
    # TODO: Currently the traverseSynteny cannot have
    # valid gene ids == 'all'
    if valid_gene_ids == 'all' or GeneCallEntry.gene in valid_gene_ids:
        return True
    else:
        return False


def allXElementsInY(x, y):
    if all(val in y for val in x):
        return True
    else:
        return False


def diffXElementsNotInY(x, y):
    x = [str(v) for v in x]
    y = [str(v) for v in y]
    return list(set(x) - set(y))


'''
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
'''


'''
def testValidGeneIDsLongerEqualSynteny(valid_gene_ids, synteny):
    if len(valid_gene_ids) >= len(synteny):
        return True
    else:
        invald_ids = diffXElementsNotInY(synteny, valid_gene_ids)
        invald_ids_print = '\n'.join([str(v) for v in invald_ids])
        print(invald_ids_print)
        raise ValueError('More genes in synteny file than valid gene ids.')
'''


# ~~~
# BELOW - WITHIN THE MAIN FOR LOOP
# ~~~


def geneInIgnore(gene, ignore):
    if gene in ignore:
        return True
    else:
        return False


def geneHasOrtholog(gene, ortholog_dic, switch=0):
    assert isinstance(gene, str), 'Gene was not passed as a string!'
    if switch == 0 and gene in ortholog_dic.X:
        return True
    elif switch == 1 and gene in ortholog_dic.Y:
        return True
    else:
        return False


def geneIsValid(gene, ignore, ortholog_dic, no_orthos):
    if geneInIgnore(gene, ignore):
        return False
    if geneHasOrtholog(gene, ortholog_dic, no_orthos):
        return True
    return False


def getOrtholog(orthologs, gene):
    assert isinstance(gene, str), 'Gene was not passed as a string!'
    try:
        return orthologs[gene]
    except ValueError:
        sys.exit(f'Gene: {gene} does not have an ortholog!')


def setOrthologGeneCall(GeneCallEntry, orthologs, _gene_calls_, switch):
    if switch == 0:
        current_ortholog = getOrtholog(orthologs.X, GeneCallEntry.gene)
    elif switch == 1:
        current_ortholog = getOrtholog(orthologs.Y, GeneCallEntry.gene)

    try:
        return _gene_calls_[switchFlag(switch)][current_ortholog]
    except KeyError:
        raise KeyError(
            f"Current ortholog:-{current_ortholog}- is not in _gene_calls_ {switchFlag(switch)}")


def appendIgnore(ignore, GeneCallEntry):
    # assert isinstance(gene, str), 'Gene was not passed as a string!'
    testGeneCallEntry(GeneCallEntry)
    if GeneCallEntry in ignore:
        return 0
    ignore.add(GeneCallEntry)
    return 0


def appendBothIgnore(ignore, switch, Gene, Ortholog):
    # TODO: modify tetGeneCallEntry to have variable vars passed
    # and loop to unpack and test
    testGeneCallEntry(Gene)
    testGeneCallEntry(Ortholog)
    appendIgnore(ignore[switch], Gene)
    appendIgnore(ignore[switchFlag(switch)], Ortholog)
    return 0


def moveStream(synteny_dic, next_gene):
    return synteny_dic[next_gene]


# Dont need below?
def moveBothStream(synteny, CurrentGene, CurrentOrtholog, switch):
    testGeneCallEntry()
    x_synteny = moveStream(synteny[switch], CurrentGene.gene)
    y_synteny = moveStream(synteny[switchFlag(switch)], CurrentOrtholog.gene)
    return [x_synteny, y_synteny]


def testGeneCallEntry(GeneCallEntry):
    assert type(GeneCallEntry) == SC.GeneCallEntry, 'GeneCallEntry not passed.'
    return 0


def directionNotNeither(direction):
    if direction == 'Neither':
        return False
    else:
        return True


def geneIsEndFlag(GeneCallEntry):
    testGeneCallEntry(GeneCallEntry)
    if GeneCallEntry.gene in ['+', '-']:
        return True


def assignDirectionEndFlag(direction):
    assert isinstance(
        direction, str), f'Invalid type for direction: {direction}'
    if direction == 'Downstream':
        end_flag = '+'
    elif direction == 'Upstream':
        end_flag = '-'
    elif direction == 'Neither':
        end_flag = 'x'
    elif direction is None:
        return
    else:
        assert ValueError(f"Direction: {direction} is invalid.")
    return end_flag


def endOfContig(GeneCallEntry, direction):
    testGeneCallEntry(GeneCallEntry)
    end_flag = assignDirectionEndFlag(direction)

    if GeneCallEntry.upstream_gene == '-':
        if end_flag in ['-']:
            return True
    elif GeneCallEntry.downstream_gene == '+':
        if end_flag in ['+']:
            return True
    else:
        return False


def assignDirection(GeneCallEntry, CurrentSeedDirection):
    testGeneCallEntry(GeneCallEntry)
    if GeneCallEntry.upstream_gene == '-':
        CurrentSeedDirection.direction_y = 'Downstream'
        return 0
    elif GeneCallEntry.downstream_gene == '+':
        CurrentSeedDirection.direction_y = 'Upstream'
        return 0
    else:
        raise ValueError(
            f"{GeneCallEntry.gene} determined as end but is not the case.")


def setNextXGene(GeneCallEntry, direction_x, synteny_x, ignore):
    testGeneCallEntry(GeneCallEntry)

    if direction_x == 'Downstream':
        if GeneCallEntry.downstream_gene == '+':
            return SC.GeneCallEntry('+')
        elif GeneCallEntry.downstream_gene not in ignore:
            # TODO - need to catch when a value isn't in a synteny dictionary
            return synteny_x[GeneCallEntry.downstream_gene]
        else:
            return SC.GeneCallEntry(None)
    elif direction_x == 'Upstream':
        if GeneCallEntry.upstream_gene == '-':
            return SC.GeneCallEntry('-')
        elif GeneCallEntry.upstream_gene not in ignore:
            return synteny_x[GeneCallEntry.upstream_gene]
        else:
            return SC.GeneCallEntry(None)
    elif direction_x == 'Neither':
        raise ValueError(f'Direction {direction_x} should not be Neither.')
    else:
        raise ValueError(f'Direction {direction_x} improperly set')


def onlyOneGeneOnContig(CurrentGene):
    assert type(CurrentGene) == SC.GeneCallEntry, 'GeneCallEntry not passed.'
    if CurrentGene.upstream_gene == '-' and CurrentGene.downstream_gene == '+':
        return True
    else:
        return False


def testGeneInSynteny(gene, synteny):
    if gene in synteny:
        return True
    else:
        raise ValueError(f'Gene {gene} is not in synteny dictionary.')


def setNextXGeneIfDirectionIsNeither(GeneCallEntry, synteny_x, ignore):
    testGeneCallEntry(GeneCallEntry)

    if GeneCallEntry.upstream_gene == '-':
        if GeneCallEntry.downstream_gene not in ignore:
            return synteny_x[GeneCallEntry.downstream_gene], 'Downstream'
        else:
            return SC.GeneCallEntry(None), 'Neither'

    elif GeneCallEntry.downstream_gene == '+':
        if GeneCallEntry.upstream_gene not in ignore:
            return synteny_x[GeneCallEntry.upstream_gene], 'Upstream'
        else:
            return SC.GeneCallEntry(None), 'Neither'
    else:
        return SC.GeneCallEntry(None), 'Neither'


def switchFlag(flag):
    if flag == 0:
        return 1
    elif flag == 1:
        return 0
    else:
        raise ValueError


'''
def addCurrentGeneToList(GeneCallEntry, list_):
    testGeneCallEntry(GeneCallEntry)
    list_.append(GeneCallEntry.gene)
    return 0
'''


def appendValues(values, switch, NextGene, NextOrtholog, SeedDirection):
    testGeneCallEntry(NextGene)
    if SeedDirection.direction_x == 'Downstream':
        values[switch].append(NextGene)
        values[switchFlag(switch)].append(NextOrtholog)
    elif SeedDirection.direction_x == 'Upstream':
        values[switch].insert(0, NextGene)
        values[switchFlag(switch)].insert(0, NextOrtholog)
    return 0


def processAttributeLine(Gene):
    testGeneCallEntry(Gene)

    attributes = Gene.attributes
    values = attributes.split(';')
    annotation_dic = {}
    for value_pair in values:
        try:
            annotation, description = value_pair.split('=', 1)
            annotation_dic[annotation] = description
        except ValueError:
            raise ValueError(
                f'Too many values to unpack for: {value_pair} in \n{Gene.gene}')

    return annotation_dic


def formatAttributeLineForWriting(Gene, ortho=False):
    annotations = ['Transfer_RNAs', 'COG20_CATEGORY', 'COG20_PATHWAY',
                   'COG20_FUNCTION', 'KEGG_Module', 'KEGG_Class', 'KOfam',
                   'FigFams', 'Pfam', 'TIGRFAM']
    if ortho:
        annotations.reverse()
    annotation_dic = processAttributeLine(Gene)
    formatted_line = ''
    for annotation in annotations:
        try:
            description = annotation_dic[annotation]
            # evalue = annotation_dic[f"{annotation}-eval"]
            formatted_line += f'{description}\t'
        except KeyError:
            raise IndexError(
                f'Annotation {annotation} not found in Gene {Gene.gene}')
    return ''.join(formatted_line.rsplit('\t', 1))  # Replace last \t with ''


def formatWriteLine(Gene=False, Ortholog=False):
    if Gene:
        testGeneCallEntry(Gene)
    if Ortholog:
        testGeneCallEntry(Ortholog)

    if Gene:
        gene = Gene.gene
        contig = Gene.contig
        up_down = f"{Gene.upstream_gene}-{Gene.downstream_gene}"
        start_stop = f"{Gene.start}-{Gene.stop}"
        annotations = formatAttributeLineForWriting(Gene)
        gene_line = '\t'.join([annotations, contig, up_down, start_stop, gene])

    if Ortholog:
        ortholog = Ortholog.gene
        contig_ortho = Ortholog.contig
        up_down_orth = f"{Ortholog.upstream_gene}-{Ortholog.downstream_gene}"
        start_stop_orth = f"{Ortholog.start}-{Ortholog.stop}"
        annotations_orth = formatAttributeLineForWriting(Ortholog, ortho=True)
        ortho_line = '\t'.join(
            [ortholog, start_stop_orth, up_down_orth, contig_ortho, annotations_orth])

    if (Gene and Ortholog):
        return '\t'.join([gene_line, ortho_line]) + '\n'
    elif Gene:
        ortho_line = '\t'.join([''] * 14)
        return '\t'.join([gene_line, ortho_line]) + '\n'
    elif Ortholog:
        gene_line = '\t'.join([''] * 14)
        return '\t'.join([gene_line, ortho_line]) + '\n'
    else:
        raise ValueError('Neither Gene nor Ortholog were provided!')


'''
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
# REFORMAT FOR WRITING ONLY, NOT LOOP-RECORDING!!!
'''


def recordSyntenyInStone(synteny_dic, Seed, values):
    testGeneCallEntry(Seed)
    synteny_dic[Seed.gene] = [values[0], values[1]]


def testNext(direction_y, upstr_gene, downstr_gene):
    if direction_y == 'Downstream':
        if upstr_gene == '-':
            return True
    elif direction_y == 'Upstream':
        if downstr_gene == '-':
            return True
    else:
        if (upstr_gene == '-' or downstr_gene == '-'):
            return True
    return False


def logEvent(file, message):
    with open(file, 'a') as out:
        out.write(f"{message}\n")
    return 0
