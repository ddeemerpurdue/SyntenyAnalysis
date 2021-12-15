'''
'''
from collections import namedtuple
import os

import syntenyLibrary as SL


def createAnnotationNamedTuple():
    return namedtuple("Annotations", "kegg, cog20fx")


def grabAnnotationFields(fields):
    annotation_indices = [12, 14]
    return [fields[v] for v in annotation_indices]


def mineSummaryFile(summary_file, genome_name):
    count = 0
    allgenes = []
    Annotations = createAnnotationNamedTuple()
    genome_summary = {}
    genome_index = 3
    gene_callers_index = 4
    with open(summary_file) as file:
        line = file.readline()  # Skip header
        line = file.readline().strip()
        while line:
            fields = line.split('\t')
            current_genome = fields[genome_index]
            if current_genome == genome_name:
                count += 1
                current_gene = fields[gene_callers_index]
                try:
                    annotations = grabAnnotationFields(fields)
                except IndexError:
                    fields.extend(['', '', '', '', '', '', '', '', ''])
                    annotations = grabAnnotationFields(fields)
                allgenes.append(current_gene)
                genome_summary[current_gene] = Annotations(*annotations)

            line = file.readline().strip()
    if count == 0:
        raise ValueError(
            f"Genome {genome_name} was not found in summary file!")

    return genome_summary


def mineOrthologFile(ortholog_file):
    # TODO - Add redundancy for multi-multi-hits
    ortholog_summary = {}
    with open(ortholog_file) as file:
        line = file.readline()
        genome_a, genome_b = line.strip().split('\t')[1:]
        line = file.readline().strip()
        while line:
            ortho_a, ortho_b = line.split('\t')[1:]
            if len(ortho_b.split(',')) > 1:  # Always go with first for B
                ortho_b = ortho_b.split(',')[0].split('|')[0]
            if len(ortho_a.split(',')) > 1:
                for ortho_x in ortho_a.split(','):
                    ortho_x = ortho_x.strip().split('|')[0]
                    ortholog_summary[ortho_x] = ortho_b
            else:
                ortholog_summary[ortho_a.split('|')[0]] = ortho_b.split('|')[0]
            line = file.readline().strip()
    return ortholog_summary, genome_a, genome_b


def initNamedTuples():
    Upstream = namedtuple("Upstream", "gene, contig, start, stop")
    Current = namedtuple("Current", "gene, contig, start, stop")
    Downstream = namedtuple("Downstream", "gene, contig, start, stop")
    return Upstream, Current, Downstream


def setStartingValues():
    return [True, '', '', '', '']


def setCurrentValues(fields):
    return [fields[val] for val in [0, 1, 2, 3]]


def sameContigAsPrevious(current, previous):
    if current == previous:
        return True
    return False


def mineGeneCalls(genecall_file, valid_gene_ids):
    genecall_summary = {}
    Upstream, Current, Downstream = initNamedTuples()
    with open(genecall_file) as file:
        line = file.readline()
        line = file.readline().strip()
        first, prev_contig, prev_gene, prev_start, prev_stop = setStartingValues()
        while line:
            fields = line.split('\t')
            current_gene, current_contig, start, stop = setCurrentValues(
                fields)
            if current_gene in valid_gene_ids:
                if not sameContigAsPrevious(current_contig, prev_contig):
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
                    if not sameContigAsPrevious(current_contig, prev_contig):
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


def grabValidGeneIds(summary_file, genome_a, genome_b, gff_b):
    valid_gene_ids_a = list(mineSummaryFile(summary_file, genome_a).keys())
    # valid_gene_ids_b = list(mineSummaryFile(summary_file, genome_b).keys())
    valid_gene_ids_b = []
    with open(gff_b) as file:
        line = file.readline()
        line = file.readline()
        while line:
            valid_gene_ids_b.append(line.split('\t')[0])
            line = file.readline()
    return valid_gene_ids_a, valid_gene_ids_b


def traverseSynteny(gffA, gffB, orthologs, summary_file):
    synteny_dic = {}

    a_orthologs, genome_a, genome_b = mineOrthologFile(orthologs)
    print('a_orthologs')
    SL.printNDicValues(a_orthologs, 10)
    genome_a = 'Bin_111_Acidaminococcaceae_Phascolarctobacterium_SM'

    b_orthologs = {v: k for k, v in a_orthologs.items()}
    print('b_orthologs')
    SL.printNDicValues(b_orthologs, 10)
    orthologs = [a_orthologs, b_orthologs]

    valid_gene_ids_a, valid_gene_ids_b = grabValidGeneIds(
        summary_file, genome_a, genome_b, gffB)
    print(valid_gene_ids_a[0:10])
    print(valid_gene_ids_b[0:10])

    print('Analyzing A_synteny')
    a_synteny = mineGeneCalls(gffA, valid_gene_ids_a)
    print('a_synteny')
    SL.printNDicValues(a_synteny, 10)
    print('Analyzing B_synteny')
    b_synteny = mineGeneCalls(gffB, valid_gene_ids_b)
    print('b_synteny')
    SL.printNDicValues(b_synteny, 10)
    synteny = [a_synteny, b_synteny]

    logfile = 'Log-Debugger.txt'

    direction_x = 'Downstream'
    direction_y = 'Neither'
    ignore = [set(), set()]
    no_orthos = []
    a_genes = list(a_synteny.keys())
    b_genes = list(b_synteny.keys())
    both_genes = [a_genes, b_genes]

    ignore_count = 0
    switch = 0
    top_count = 0
    append_ = False
    break_ = False

    for current_gene in both_genes[switch]:
        switch = 0
        top_count += 1
        seed_direction = 'Forward'
        direction_x = 'Downstream'
        direction_y = 'Upstream'

        if SL.geneInIgnore(current_gene, ignore[switch]):
            ignore_count += 1
        elif not SL.geneHasOrtholog(current_gene, orthologs[switch], no_orthos):
            pass
        else:
            values = [[], []]

            current_ortholog = SL.getOrtholog(orthologs[switch], current_gene)
            SL.appendBothIgnore(ignore, switch, current_gene, current_ortholog)

            x_synteny, y_synteny = SL.moveBothStream(
                synteny, current_gene, current_ortholog, switch)

            seed = x_synteny
            seed_ortholog = y_synteny

            SL.addCurrentGeneToList(
                synteny[switch][current_gene], values[switch])
            SL.addCurrentGeneToList(
                synteny[SL.switchFlag(switch)][current_ortholog], values[SL.switchFlag(switch)])

            while True:
                if SL.directionNotNeither(direction_x):
                    next_ = SL.setNextXGene(
                        direction_x, x_synteny, ignore[switch])
                else:
                    if SL.onlyOneGeneOnContig(x_synteny):
                        append_ = True
                        message = f"{x_synteny[1]}+{y_synteny[1]} broke after switching due to immediate contig end."
                        SL.logEvent(logfile, message)
                        switch = 0
                        if seed_direction == 'Forward':
                            seed_direction = 'Reverse'
                            direction_x = 'Upstream'
                            x_synteny, y_synteny = SL.moveBothStream(
                                synteny, seed[1].gene, seed_ortholog[1].gene, switch)
                        else:
                            break_ = True
                    else:
                        next_, direction_x = SL.setNextXGeneIfDirectionIsNeither(
                            x_synteny, ignore[switch])

                if next_ is None:
                    break_ = True
                    append_ = False
                elif SL.geneHasOrtholog(next_.gene, orthologs[switch], no_orthos):
                    nexts_ortholog = SL.setCurrentOrthologsSynteny(
                        orthologs, switch, next_.gene, synteny)
                    Os_upstream, current_ortholog, Os_downstream = nexts_ortholog

                    if current_ortholog.gene == y_synteny[2].gene:
                        direction_y = 'Downstream'
                        append_ = True

                        x_synteny, y_synteny = SL.moveBothStream(
                            synteny, next_.gene, current_ortholog.gene, switch)

                    elif current_ortholog.gene == y_synteny[0].gene:
                        direction_y = 'Upstream'
                        append_ = True

                        x_synteny, y_synteny = SL.moveBothStream(
                            synteny, next_.gene, current_ortholog.gene, switch)

                    else:
                        if ((y_synteny[0].gene == '-' and direction_y in ['Downstream', 'Neither']) or
                                (y_synteny[2].gene == '+' and direction_y in ['Upstream', 'Neither'])):
                            if Os_upstream.gene == '-':
                                direction_y = 'Downstream'
                                append_ = True

                                x_synteny, y_synteny = SL.moveBothStream(
                                    synteny, next_.gene, current_ortholog.gene, switch)

                            elif Os_downstream.gene == '+':
                                direction_y = 'Upstream'
                                append_ = True

                                x_synteny, y_synteny = SL.moveBothStream(
                                    synteny, next_.gene, current_ortholog.gene, switch)

                            else:
                                message = "After switching, neither up/downstr of the next ortholog was the end of a contig!"
                                message = message + \
                                    f" {x_synteny[1]}+{y_synteny[1]}"
                                SL.logEvent(logfile, message)
                                switch = 0

                                if seed_direction == 'Forward':
                                    seed_direction = 'Reverse'
                                    direction_x = 'Upstream'
                                    x_synteny, y_synteny = SL.moveBothStream(
                                        synteny, seed[1].gene, seed_ortholog[1].gene, switch)
                                else:
                                    break_ = True
                        else:
                            message = 'Current lagging is not at start/end and next_.genes ortholog doesnt match upstr_/dwnstr_match'
                            message = message + \
                                f" {x_synteny[1]}+{y_synteny[1]}"
                            SL.logEvent(logfile, message)
                            switch = 0

                            if seed_direction == 'Forward':
                                seed_direction = 'Reverse'
                                direction_x = 'Upstream'
                                x_synteny, y_synteny = SL.moveBothStream(
                                    synteny, seed[1].gene, seed_ortholog[1].gene, switch)
                            else:
                                break_ = True

                elif next_.gene in ['+', '-']:
                    x_synteny, y_synteny = y_synteny, x_synteny
                    direction_x = direction_y
                    direction_y = 'Neither'

                    if direction_x == 'Downstream':
                        if x_synteny[2].gene == '+':
                            message = f'Contigs end at same spot: {x_synteny[1]}+{y_synteny[1]}'
                            SL.logEvent(logfile, message)
                            switch = 0
                            if seed_direction == 'Forward':
                                seed_direction = 'Reverse'
                                direction_x = 'Upstream'
                                x_synteny, y_synteny = SL.moveBothStream(
                                    synteny, seed[1].gene, seed_ortholog[1].gene, switch)
                            else:
                                break_ = True
                        else:
                            switch = SL.switchFlag(switch)
                    elif direction_x == 'Upstream':
                        if x_synteny[0].gene == '-':
                            message = f'Contigs end at same spot: {x_synteny[1]}+{y_synteny[1]}'
                            SL.logEvent(logfile, message)
                            switch = 0
                            if seed_direction == 'Forward':
                                seed_direction = 'Reverse'
                                direction_x = 'Upstream'

                                x_synteny, y_synteny = SL.moveBothStream(
                                    synteny, seed[1].gene, seed_ortholog[1].gene, switch)
                            else:
                                break_ = True
                        else:
                            switch = SL.switchFlag(switch)
                    else:
                        switch = SL.switchFlag(switch)
                else:
                    message = f'No ortholog for {x_synteny[1]}'
                    SL.logEvent(logfile, message)
                    switch = 0
                    if seed_direction == 'Forward':
                        seed_direction = 'Reverse'
                        direction_x = 'Upstream'
                        x_synteny, y_synteny = SL.moveBothStream(
                            synteny, seed[1].gene, seed_ortholog[1].gene, switch)
                    else:
                        break_ = True
                if append_:
                    SL.appendBothIgnore(
                        ignore, switch, next_.gene, current_ortholog.gene)
                    SL.appendValues(values, switch, next_.gene,
                                    current_ortholog.gene, seed_direction)
                if break_:
                    SL.recordSyntenyInStone(synteny_dic, seed[1], values)
                    append_ = False
                    break_ = False
                    break
                append_ = False

    print(f"Ignore1 :{len(ignore[0])} Ignore2: {len(ignore[1])}")
    print(
        f"Counts:\nTop count: {top_count}\nIgnore: {ignore_count}")
    print(f"Set of no_orthos: {len(set(no_orthos))}")
    print(f"No orthos: {no_orthos}")
    print(f"{set(no_orthos)}")
    return synteny_dic


if __name__ == '__main__':
    summary_file = 'rawdata/AllPS-Genomes-Clean-21Sep21_gene_clusters_summary.NOAA.txt'
    genome = 'Bin_111_Acidaminococcaceae_Phascolarctobacterium_SM'

    orth_prefix = 'rawdata/Orthologs-Acidaminococcaceae/Orthologues_111_Acidaminococcaceae_Phascolarctobacterium_SM/'
    orth_file = '111_Acidaminococcaceae_Phascolarctobacterium_SM__v__111_Acidaminococcaceae_Phascolarctobacterium-SM-GCA_014648015.1.tsv'
    ortholog_file = os.path.join(orth_prefix, orth_file)

    gffA = 'rawdata/GFF-Internal/SupernatantGeneCalls.txt'
    gffL = 'rawdata/GFF-Acidaminococcaceae/111_Acidaminococcaceae_Phascolarctobacterium-S*-GCA_014648015.1.gff'

    summary = mineSummaryFile(summary_file, genome)
    SL.printNDicValues(summary, 10)

    a = traverseSynteny(gffA, gffL, ortholog_file, summary_file)

    output = f'Output/test.txt'
    with open(output, 'w') as out:
        out.write(f"Algoriphagus__UCCA\tAlgoriphagus_marincola\n")
        for seed in a:
            for gene_a, gene_b in zip(a[seed][0], a[seed][1]):
                if gene_a in ['-', '+']:
                    vals = '\t'.join(['', '', '', '', '', ''])
                else:
                    vals = summary[gene_a]
                    vals = '\t'.join(vals)
                out.write(f"{gene_a}\t{gene_b}\t{vals}\n")
            out.write(f"X\tX\n")
