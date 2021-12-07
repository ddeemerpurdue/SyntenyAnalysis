import sys


def sameContigAsUpstream(current, previous):
    if current == previous:
        return True
    return False


def geneInIgnore(gene, ignore):
    if gene in ignore:
        return True
    else:
        return False


def geneHasOrtholog(gene, ortholog_dic, no_orthos):
    if gene in ortholog_dic:
        return True
    no_orthos.append(gene)
    return False


def geneIsValid(gene, ignore, ortholog_dic, no_orthos):
    if geneInIgnore(gene, ignore):
        return False
    if geneHasOrtholog(gene, ortholog_dic, no_orthos):
        return True
    return False


def getOrtholog(orthologs, gene):
    return orthologs[gene]


def addCurrentGeneToList(synteny, list_):
    return list_.append(synteny[1].gene)


def directionNotNeither(direction):
    if direction == 'Neither':
        return False
    else:
        return True


def setNextXGene(direction_x, x_synteny, ignore):
    if direction_x == 'Downstream':
        if x_synteny[2].gene not in ignore:
            return x_synteny[2]
        else:
            return None
    elif direction_x == 'Upstream':
        if x_synteny[0].gene not in ignore:
            return x_synteny[0]
        else:
            return None
    elif direction_x == 'Neither':
        return None
    else:
        print('Direction improperly set')
        sys.exit()


def setNextXGeneIfDirectionIsNeither(x_synteny, ignore):
    # Test A: if one ends, and B if in ignore
    if x_synteny[0].gene == '-':
        if x_synteny[2].gene not in ignore:
            return x_synteny[2], 'Downstream'
        else:
            return None, 'Neither'
    elif x_synteny[2].gene == '+':
        if x_synteny[0].gene not in ignore:
            return x_synteny[0], 'Upstream'
        else:
            return None, 'Neither'
    else:
        return None, 'Neither'


def onlyOneGeneOnContig(x_synteny):
    if x_synteny[0].gene == '-' and x_synteny[2].gene == '+':
        return True
    else:
        return False


def setCurrentOrthologsSynteny(orthologs, switch, gene, synteny):
    current_ortholog = getOrtholog(orthologs[switch], gene)
    Os_synteny = synteny[switchFlag(switch)][current_ortholog]
    return Os_synteny


def moveStream(synteny_dic, next_gene):
    return synteny_dic[next_gene]


def moveBothStream(synteny, gene, ortholog, switch):
    x_synteny = moveStream(synteny[switch], gene)
    y_synteny = moveStream(synteny[switchFlag(switch)], ortholog)
    return [x_synteny, y_synteny]


def switchFlag(flag):
    if flag == 0:
        return 1
    elif flag == 1:
        return 0
    else:
        raise ValueError


def appendIgnore(ignore, gene):
    ignore.add(gene)
    return 0


def appendBothIgnore(ignore, switch, gene, ortholog):
    appendIgnore(ignore[switch], gene)
    appendIgnore(ignore[switchFlag(switch)], ortholog)
    return 0


def appendValues(values, switch, gene, match_gene, direct):
    if direct == 'Forward':
        values[switch].append(gene)
        values[switchFlag(switch)].append(match_gene)
    elif direct == 'Reverse':
        values[switch].insert(0, gene)
        values[switchFlag(switch)].insert(0, match_gene)
    return 0


def recordSyntenyInStone(synteny_dic, seed, values):
    synteny_dic[seed.gene] = [values[0], values[1]]


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
