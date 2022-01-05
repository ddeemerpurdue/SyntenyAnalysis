import syntenyClasses as SC
import syntenyLibrary as SL


def traverseSynteny(gffA_file, gffB_file, ortholog_file, override=False):
    synteny_dic = {}

    genome_a, genome_b = SL.getGenomeNames(ortholog_file, override)

    _orthologs_ = SL.mineOrthologFile(ortholog_file)

    A_Genecalls = SL.mineGff3File(gffA_file)
    B_Genecalls = SL.mineGff3File(gffB_file)
    _gene_calls_ = [A_Genecalls, B_Genecalls]

    ignore = [set(), set()]
    no_orthos = []

    ignore_count = 0
    append_ = False

    runs = 0
    top_ignore = 0
    mid_top_ignore = 0
    body_ignore = 0

    for top_count, gene in enumerate(A_Genecalls):
        switch = 0
        CurrentSeedDirection = SC.SeedDirection()
        CurrentGene = A_Genecalls[gene]

        if SL.geneInIgnore(CurrentGene, ignore[switch]):
            top_ignore += 1
            ignore_count += 1
            continue
        elif not SL.geneHasOrtholog(gene, _orthologs_):
            no_orthos.append(gene)
            continue

        CurrentOrtholog = SL.setOrthologGeneCall(
            CurrentGene, _orthologs_, _gene_calls_, switch)

        SL.appendBothIgnore(ignore, switch, CurrentGene, CurrentOrtholog)  # IC

        mid_top_ignore += 1

        # Maybe a class or a named tuple?
        values = [[], []]
        #  loc_info = [[], []]  # TODO: REMOVE - DO NOT NEED!!!

        SL.appendValues(values, switch, CurrentGene,
                        CurrentOrtholog, CurrentSeedDirection)

        # SL.appendMiscValues(loc_info, switch, CurrentGene,
        #                     CurrentOrtholog, CurrentSeedDirection)

        Seed = CurrentGene
        SeedOrtholog = CurrentOrtholog

        while True:
            if SL.directionNotNeither(CurrentSeedDirection.direction_x):
                NextGene = SL.setNextXGene(
                    CurrentGene, CurrentSeedDirection.direction_x, _gene_calls_[switch], ignore[switch])
            else:
                if SL.onlyOneGeneOnContig(CurrentGene):
                    append_ = True
                    switch = 0  # Dont think I need this
                    if CurrentSeedDirection.restartFromSeed():
                        CurrentGene, CurrentOrtholog = Seed, SeedOrtholog
                    else:
                        break
                else:
                    NextGene, CurrentSeedDirection.direction_x = SL.setNextXGeneIfDirectionIsNeither(
                        CurrentGene, _gene_calls_[switch], ignore[switch])

            # SUB 1A - NEXT GENE IS NONE
            if NextGene.gene is None:
                break

            # SUB 1B - NEXT GENE HAS AN ORTHOLOG
            elif SL.geneHasOrtholog(NextGene.gene, _orthologs_, switch=switch):
                NextOrtholog = SL.setOrthologGeneCall(
                    NextGene, _orthologs_, _gene_calls_, switch)

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
                            break

            # SUB 1C - NEXT GENE IS ACTUALLY A +/-
            elif SL.geneIsEndFlag(NextGene):
                CurrentGene, CurrentOrtholog = CurrentOrtholog, CurrentGene
                CurrentSeedDirection.XtoYandYtoNeither()

                if SL.endOfContig(CurrentGene, CurrentSeedDirection.direction_x) or SL.onlyOneGeneOnContig(CurrentGene):
                    switch = 0
                    if CurrentSeedDirection.restartFromSeed():
                        CurrentGene, CurrentOrtholog = Seed, SeedOrtholog
                    else:
                        break
                else:
                    switch = SL.switchFlag(switch)

            else:
                # if CurrentSeedDirection.direction_x == 'Downstream':
                #     no_orthos.append(gene)
                switch = 0
                if CurrentSeedDirection.restartFromSeed():
                    CurrentGene, CurrentOrtholog = Seed, SeedOrtholog
                else:
                    break

            if append_:
                SL.appendBothIgnore(
                    ignore, switch, NextGene, NextOrtholog)
                SL.appendValues(values, switch, CurrentGene,
                                CurrentOrtholog, CurrentSeedDirection)
                #  SL.appendMiscValues(loc_info, switch, NextGene,
                #                    NextOrtholog, CurrentSeedDirection)
                append_ = False
        # Out of while loop
        runs += 1
        SL.recordSyntenyInStone(synteny_dic, Seed, values)
        append_ = False

    print(f"Total genes in A: {len(A_Genecalls.keys())}")
    print(f"Total genes in B: {len(B_Genecalls.keys())}")
    print(f"Runs: {runs}")
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
    a_not_written = set(A_Genecalls.keys()) - all_written_x
    b_not_written = set(B_Genecalls.keys()) - all_written_y
    print(f"Not written for a: {len(a_not_written)}")
    print(f"Not written for b: {len(b_not_written)}")

    a_not = {}
    b_not = {}
    for gene in a_not_written:
        a_not[gene] = A_Genecalls[gene]
    for gene in b_not_written:
        b_not[gene] = B_Genecalls[gene]
    # print(a_not)

    return [synteny_dic, a_not, b_not]
