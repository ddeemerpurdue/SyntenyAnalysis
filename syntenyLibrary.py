'''
'''
from collections import namedtuple
import os
import sys

import syntenyClasses as SC


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


def mineGeneCalls(genecall_file, valid_gene_ids):
    gene_loci = {}
    with open(genecall_file) as file:  # TODO: Below
        if 'GCA' in genecall_file or 'GCF' in os.path.basename(genecall_file):
            if '-' in os.path.basename(genecall_file):
                next(file)  # Skip header :-)

        current_line = file.readline().strip()
        upstream_line = None
        UpstreamGeneEntry = SC.GeneCallEntry('-')
        while current_line:
            CurrentGeneEntry = SC.GeneCallEntry(current_line)

            if not geneInValidGeneIDs(CurrentGeneEntry, valid_gene_ids):
                upstream_line = current_line
                current_line = file.readline().strip()
                continue

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


def grabValidGeneIds(gene_to_annotation, gff=False):
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


def testValidGeneIDsLongerEqualSynteny(valid_gene_ids, synteny):
    if len(valid_gene_ids) >= len(synteny):
        return True
    else:
        invald_ids = diffXElementsNotInY(synteny, valid_gene_ids)
        invald_ids_print = '\n'.join([str(v) for v in invald_ids])
        print(invald_ids_print)
        raise ValueError('More genes in synteny file than valid gene ids.')


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

    return orthologs[gene]


def setCurrentOrthologsSynteny(GeneCallEntry, orthologs, synteny, switch):
    if switch == 0:
        current_ortholog = getOrtholog(orthologs.X, GeneCallEntry.gene)
    elif switch == 1:
        current_ortholog = getOrtholog(orthologs.Y, GeneCallEntry.gene)

    try:
        return synteny[switchFlag(switch)][current_ortholog]
    except KeyError:
        raise KeyError(
            f"Current ortholog:-{current_ortholog}- is not in synteny {switchFlag(switch)}")


def appendIgnore(ignore, gene):
    assert isinstance(gene, str), 'Gene was not passed as a string!'
    if gene in ignore:
        return 0
    ignore.add(gene)
    return 0


def appendBothIgnore(ignore, switch, gene, ortholog):
    appendIgnore(ignore[switch], gene)
    appendIgnore(ignore[switchFlag(switch)], ortholog)
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


def addCurrentGeneToList(GeneCallEntry, list_):
    testGeneCallEntry(GeneCallEntry)
    list_.append(GeneCallEntry.gene)
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


def appendValues(values, switch, NextGene, NextOrtholog, direct):
    testGeneCallEntry(NextGene)
    if direct == 'Forward':
        values[switch].append(NextGene.gene)
        values[switchFlag(switch)].append(NextOrtholog.gene)
    elif direct == 'Reverse':
        values[switch].insert(0, NextGene.gene)
        values[switchFlag(switch)].insert(0, NextOrtholog.gene)
    return 0


def appendMiscValues(loc_info, switch, NextGene, NextOrtholog, direct):
    testGeneCallEntry(NextGene)
    if direct == 'Forward':
        loc_info[switch].append(
            f"{NextGene.upstream_gene}-{NextGene.downstream_gene}\t{NextGene.start}-{NextGene.stop}")
        loc_info[switchFlag(switch)].append(
            f"{NextOrtholog.upstream_gene}-{NextOrtholog.downstream_gene}\t{NextOrtholog.start}-{NextOrtholog.stop}")
    elif direct == 'Reverse':
        loc_info[switch].insert(
            0, f"{NextGene.upstream_gene}-{NextGene.downstream_gene}\t{NextGene.start}-{NextGene.stop}")
        loc_info[switchFlag(switch)].insert(
            0, f"{NextOrtholog.upstream_gene}-{NextOrtholog.downstream_gene}\t{NextOrtholog.start}-{NextOrtholog.stop}")
    return 0


def recordSyntenyInStone(synteny_dic, Seed, values, loc_info):
    testGeneCallEntry(Seed)
    synteny_dic[Seed.gene] = [values[0], values[1], loc_info[0], loc_info[1]]


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
