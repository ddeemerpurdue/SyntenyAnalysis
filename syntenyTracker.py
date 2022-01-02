import syntenyClasses as SC
import syntenyLibrary as SL


def traverseSynteny(summary_file, gffA_file, gffB_file, ortholog_file,
                    override=False):
    synteny_dic = {}

    genome_a, genome_b = SL.getGenomeNames(ortholog_file, override)

    _gene_to_annotation_ = SL.mineSummaryFile(summary_file, genome_a)
    _orthologs_ = SL.mineOrthologFile(ortholog_file)

    SL.testSummaryGeneIDsMatchOrthologFile(_gene_to_annotation_, _orthologs_.X)

    # Have the option to grab the validGeneIDs from summary_file or GFF
    valid_gene_ids_a = SL.grabValidGeneIds(_gene_to_annotation_, gffA_file)
    valid_gene_ids_b = SL.grabValidGeneIds(_gene_to_annotation_, gffB_file)

    a_synteny = SL.mineGeneCalls(gffA_file, valid_gene_ids_a)
    b_synteny = SL.mineGeneCalls(gffB_file, valid_gene_ids_b)
    synteny = [a_synteny, b_synteny]

    SL.testValidGeneIDsLongerEqualSynteny(valid_gene_ids_a, a_synteny)
    SL.testValidGeneIDsLongerEqualSynteny(valid_gene_ids_b, b_synteny)

    ignore = [set(), set()]
    no_orthos = []

    ignore_count = 0
    append_ = False
    break_ = False

    # Go through each valid_gene_id
    for top_count, gene in enumerate(valid_gene_ids_a):
        switch = 0
        CurrentSeedDirection = SC.SeedDirection()

        if SL.geneInIgnore(gene, ignore[switch]):
            ignore_count += 1
            continue
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
        loc_info = [[], []]
        # Add the current gene and current ortholog to values
        SL.addCurrentGeneToList(CurrentGene, values[0])  # Internal check
        SL.addCurrentGeneToList(CurrentOrtholog, values[1])
        # SL.addCurrentGeneToList(CurrentGene, loc_info[0])  # Fx - rearrange
        # SL.addCurrentGeneToList(CurrentOrtholog, loc_info[1])
        SL.appendMiscValues(loc_info, switch, CurrentGene,
                            CurrentOrtholog, CurrentSeedDirection.direction_x)

        Seed = CurrentGene
        SeedOrtholog = CurrentOrtholog

        while True:
            if SL.directionNotNeither(CurrentSeedDirection.direction_x):
                NextGene = SL.setNextXGene(
                    CurrentGene, CurrentSeedDirection.direction_x, synteny[switch], ignore[switch])
            else:
                if SL.onlyOneGeneOnContig(CurrentGene):
                    append_ = True
                    switch = 0  # Dont think I need this
                    if CurrentSeedDirection.restartFromSeed():
                        CurrentGene, CurrentOrtholog = Seed, SeedOrtholog
                    else:
                        break_ = True
                else:
                    NextGene, CurrentSeedDirection.direction_x = SL.setNextXGeneIfDirectionIsNeither(
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

                if NextOrtholog.gene == CurrentOrtholog.downstream_gene:
                    CurrentSeedDirection.direction_y = 'Downstream'
                    append_ = True
                    CurrentGene, CurrentOrtholog = NextGene, NextOrtholog

                elif NextOrtholog.gene == CurrentOrtholog.upstream_gene:
                    CurrentSeedDirection.direction_y = 'Upstream'
                    append_ = True
                    CurrentGene, CurrentOrtholog = NextGene, NextOrtholog

                else:
                    if (SL.endOfContig(CurrentOrtholog, CurrentSeedDirection.direction_y) and
                            SL.endOfContig(NextOrtholog, direction='Neither')):
                        SL.assignDirection(NextOrtholog)
                        CurrentGene, CurrentOrtholog = NextGene, NextOrtholog
                        append_ = True
                    else:
                        switch = 0
                        if CurrentSeedDirection.restartFromSeed():
                            CurrentGene, CurrentOrtholog = Seed, SeedOrtholog
                        else:
                            break_ = True

            # SUB 1C - NEXT GENE IS ACTUALLY A +/-
            elif SL.geneIsEndFlag(NextGene):
                CurrentGene, CurrentOrtholog = CurrentOrtholog, CurrentGene
                CurrentSeedDirection.XtoYandYtoNeither()

                if SL.endOfContig(CurrentGene, CurrentSeedDirection.direction_x) or SL.onlyOneGeneOnContig(CurrentGene):
                    switch = 0
                    if CurrentSeedDirection.restartFromSeed():
                        CurrentGene, CurrentOrtholog = Seed, SeedOrtholog
                    else:
                        break_ = True
                else:
                    switch = SL.switchFlag(switch)

            else:
                switch = 0
                if CurrentSeedDirection.restartFromSeed():
                    CurrentGene, CurrentOrtholog = Seed, SeedOrtholog
                else:
                    break_ = True

            if append_:
                SL.appendBothIgnore(
                    ignore, switch, NextGene.gene, NextOrtholog.gene)
                SL.appendValues(values, switch, NextGene,
                                NextOrtholog, CurrentSeedDirection.direction_x)
                SL.appendMiscValues(loc_info, switch, NextGene,
                                    NextOrtholog, CurrentSeedDirection.direction_x)
            if break_:
                SL.recordSyntenyInStone(synteny_dic, Seed, values, loc_info)
                append_ = False
                break_ = False
                break
            append_ = False

    print(f"Ignore1 :{len(ignore[0])} Ignore2: {len(ignore[1])}")
    print(f"Ignore1 :{len(set(ignore[0]))} Ignore2: {len(set(ignore[1]))}")
    print(
        f"Counts:\nTop count: {top_count}\nIgnore: {ignore_count}")
    print(f"Set of no_orthos: {len(set(no_orthos))}")
    a_not_written = set(valid_gene_ids_a) - set(ignore[0])
    b_not_written = set(valid_gene_ids_b) - set(ignore[1])
    print(f"Not written for a: {len(a_not_written)}")
    print(f"Not written for b: {len(b_not_written)}")

    return synteny_dic, [a_not_written, b_not_written]
