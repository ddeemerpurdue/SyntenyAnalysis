'''
'''
import syntenyLibrary as SL


def traverseSynteny(summary_file, gffA_file, gffB_file, ortholog_file,
                    override=False, Ga=False, Gb=False):
    synteny_dic = {}

    # Get the genome names for mining summary file
    genome_a, genome_b = SL.getGenomeNames(ortholog_file, override)
    # gene_to_annotation just provides us with annotations
    # Below has a built-in check that genome_a is indeed in the summary file
    _gene_to_annotation_ = SL.mineSummaryFile(summary_file, genome_a)

    # Interal checks
    _orthologs_ = SL.mineOrthologFile(ortholog_file)
    # Checks we can match up these ortholog ids to the summary file
    SL.testSummaryGeneIDsMatchOrthologFile(_gene_to_annotation_, _orthologs_.X)

    # Have the option to grab the validGeneIDs from summary_file or GFF
    valid_gene_ids_a = SL.grabValidGeneIds(_gene_to_annotation_, gffA_file)
    valid_gene_ids_b = SL.grabValidGeneIds(_gene_to_annotation_, gffB_file)

    print('Analyzing A_synteny')
    print(f"GFF_a file: {gffA_file}")
    a_synteny = SL.mineGeneCalls(gffA_file, valid_gene_ids_a)

    print('Analyzing B_synteny')
    print(f"GFF_b file: {gffB_file}")
    b_synteny = SL.mineGeneCalls(gffB_file, valid_gene_ids_b)
    synteny = [a_synteny, b_synteny]

    SL.testValidGeneIDsLongerEqualSynteny(valid_gene_ids_a, a_synteny)
    SL.testValidGeneIDsLongerEqualSynteny(valid_gene_ids_b, b_synteny)

    direction_x = 'Downstream'
    direction_y = 'Neither'
    ignore = [set(), set()]
    no_orthos = []

    ignore_count = 0
    switch = 0
    top_count = 0
    append_ = False
    break_ = False

    # Go through each valid_gene_id
    for gene in valid_gene_ids_a:
        switch = 0
        top_count += 1
        seed_direction = 'Forward'
        direction_x = 'Downstream'
        direction_y = 'Upstream'

        # Turn ignore into a namedtuple?
        if SL.geneInIgnore(gene, ignore[switch]):
            ignore_count += 1
            continue
        # Split the below function to add to no_orthos
        elif not SL.geneHasOrtholog(gene, _orthologs_):
            no_orthos.append(gene)
            continue

        # Below, don't like how the order changes from 'geneHasOrtholog'
        CurrentGene = a_synteny[gene]
        ortholog = SL.getOrtholog(_orthologs_.X, CurrentGene.gene)  # IC
        CurrentOrtholog = b_synteny[ortholog]

        SL.appendIgnore(ignore[0], CurrentGene.gene)  # IC
        SL.appendIgnore(ignore[1], CurrentOrtholog.gene)  # IC

        # Maybe a class or a named tuple?
        values = [[], []]
        debug = [[], []]
        # Add the current gene and current ortholog to values
        SL.addCurrentGeneToList(CurrentGene, values[0])  # Internal check
        SL.addCurrentGeneToList(CurrentOrtholog, values[1])

        Seed = CurrentGene
        SeedOrtholog = CurrentOrtholog

        while True:
            if SL.directionNotNeither(direction_x):  # TODO - IC
                NextGene = SL.setNextXGene(
                    CurrentGene, direction_x, synteny[switch], ignore[switch])
            else:
                if SL.onlyOneGeneOnContig(CurrentGene):
                    append_ = True
                    switch = 0  # Dont think I need this
                    if seed_direction == 'Forward':
                        seed_direction = 'Reverse'
                        direction_x = 'Upstream'
                        CurrentGene, CurrentOrtholog = Seed, SeedOrtholog
                    else:
                        break_ = True
                else:
                    # BELOW - SWITCH IS MESSED UP - 1 WHEN SHOULD BE 0
                    NextGene, direction_x = SL.setNextXGeneIfDirectionIsNeither(
                        CurrentGene, synteny[switch], ignore[switch])

            # TODO - Change this break - add the stuff after while broken?
            # SUB 1A - NEXT GENE IS NONE
            if NextGene.gene is None:
                break_ = True
                append_ = False

            # SUB 1B - NEXT GENE HAS AN ORTHOLOG
            elif SL.geneHasOrtholog(NextGene.gene, _orthologs_, switch=switch):
                NextOrtholog = SL.setCurrentOrthologsSynteny(
                    NextGene, _orthologs_, synteny, switch)

                # Below is making a run without switching contigs/directions
                # TODO - simplify
                if NextOrtholog.gene == CurrentOrtholog.downstream_gene:
                    direction_y = 'Downstream'
                    append_ = True
                    CurrentGene, CurrentOrtholog = NextGene, NextOrtholog

                elif NextOrtholog.gene == CurrentOrtholog.upstream_gene:
                    direction_y = 'Upstream'
                    append_ = True
                    CurrentGene, CurrentOrtholog = NextGene, NextOrtholog

                else:
                    # Below, means let's try switching contigs or directions
                    if SL.endOfContig(CurrentOrtholog, direction_y):
                        if SL.endOfContig(NextOrtholog, direction='Neither'):
                            direction_y, append_ = SL.assignDirection(
                                NextOrtholog)
                            CurrentGene, CurrentOrtholog = NextGene, NextOrtholog
                            append_ = True
                        else:
                            switch = 0
                            # Below, we want to flip and search seed in opposite
                            if seed_direction == 'Forward':
                                seed_direction = 'Reverse'
                                direction_x = 'Upstream'
                                CurrentGene, CurrentOrtholog = Seed, SeedOrtholog
                            else:
                                break_ = True
                    else:
                        switch = 0
                        # Below, we want to flip and search seed in opposite
                        if seed_direction == 'Forward':
                            seed_direction = 'Reverse'
                            direction_x = 'Upstream'
                            CurrentGene, CurrentOrtholog = Seed, SeedOrtholog
                        else:
                            break_ = True

            # SUB 1C - NEXT GENE IS ACTUALLY A +/-
            elif SL.geneIsEndFlag(NextGene):
                CurrentGene, CurrentOrtholog = CurrentOrtholog, CurrentGene
                direction_x = direction_y
                direction_y = 'Neither'

                if SL.endOfContig(CurrentGene, direction_x) or SL.onlyOneGeneOnContig(CurrentGene):
                    switch = 0
                    if seed_direction == 'Forward':
                        seed_direction = 'Reverse'
                        direction_x = 'Upstream'
                        CurrentGene, CurrentOrtholog = Seed, SeedOrtholog
                    else:
                        break_ = True
                else:
                    switch = SL.switchFlag(switch)

            else:
                switch = 0
                if seed_direction == 'Forward':
                    seed_direction = 'Reverse'
                    direction_x = 'Upstream'
                    CurrentGene, CurrentOrtholog = Seed, SeedOrtholog
                else:
                    break_ = True

            if append_:
                SL.appendBothIgnore(
                    ignore, switch, NextGene.gene, NextOrtholog.gene)
                SL.appendValues(values, switch, NextGene,
                                NextOrtholog, seed_direction)
                # SL.appendMiscValues(debug, switch, CurrentGene,
                #                     CurrentOrtholog, seed_direction)
            if break_:
                SL.recordSyntenyInStone(synteny_dic, Seed, values)
                append_ = False
                break_ = False
                break
            append_ = False

    print(f"Ignore1 :{len(ignore[0])} Ignore2: {len(ignore[1])}")
    print(
        f"Counts:\nTop count: {top_count}\nIgnore: {ignore_count}")
    print(f"Set of no_orthos: {len(set(no_orthos))}")
    # print(f"No orthos: {no_orthos}")
    # print(f"{set(no_orthos)}")
    return synteny_dic
