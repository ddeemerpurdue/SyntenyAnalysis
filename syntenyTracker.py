import syntenyClasses as SC
import syntenyLibrary as SL


def traverseSynteny(gffA_file, gffB_file, ortholog_file, override=False):
    synteny_dic = {}

    genome_a, genome_b = SL.getGenomeNames(ortholog_file, override)

    _orthologs_ = SL.mineOrthologFile(ortholog_file)
    print(_orthologs_)

    A_Genecalls = SL.mineGff3File(gffA_file, type_field="gene")
    B_Genecalls = SL.mineGff3File(gffB_file, type_field="gene")
    _gene_calls_ = [A_Genecalls, B_Genecalls]
    # tmp = _gene_calls_[0].keys()
    # tmp = sorted(tmp)
    # print('A_Genecalls')
    # print(tmp)
    # tmp = _gene_calls_[1].keys()
    # tmp = sorted(tmp)
    # print('B_Genecalls')
    # print(tmp)

    ignore = [set(), set()]
    no_orthos = [set(), set()]
    already_analyzed_a = []
    already_analyzed_b = []

    ignore_count = 0
    append_ = False

    runs = 0

    for top_count, gene in enumerate(A_Genecalls):
        switch = 0
        CurrentSeedDirection = SC.SeedDirection()
        CurrentGene = A_Genecalls[gene]

        if SL.geneInIgnore(CurrentGene, ignore[switch]):  # already analyzed
            ignore_count += 1
            continue
        if not SL.geneHasOrtholog(gene, _orthologs_):
            no_orthos[switch].add(CurrentGene.gene)
            continue

        CurrentOrtholog = SL.setOrthologGeneCall(
            CurrentGene, _orthologs_, _gene_calls_, switch)
        SL.appendBothIgnore(ignore, switch, CurrentGene, CurrentOrtholog)
        values = [[], []]

        SL.appendValues(values, switch, CurrentGene,
                        CurrentOrtholog, CurrentSeedDirection)

        Seed = CurrentGene
        SeedOrtholog = CurrentOrtholog

        while True:
            # SETTING THE NEXT GENE
            if SL.onlyOneGeneOnContig(CurrentGene):
                append_ = True
                if CurrentSeedDirection.restartFromSeed():
                    CurrentGene, CurrentOrtholog = Seed, SeedOrtholog
                    continue
                else:
                    break
            elif SL.directionXNotNeither(CurrentSeedDirection):
                NextGene = SL.setNextXGene(
                    CurrentGene, CurrentSeedDirection, _gene_calls_[switch], ignore[switch])
            else:
                NextGene, CurrentSeedDirection.direction_x = SL.setNextXGeneIfDirectionIsNeither(
                    CurrentGene, _gene_calls_[switch], ignore[switch])

            if NextGene.gene is None:
                break

            # NextGene is [+, -]
            if SL.geneIsEndFlag(NextGene):
                CurrentGene, CurrentOrtholog = CurrentOrtholog, CurrentGene
                CurrentSeedDirection.XtoYandYtoNeither()

                if SL.endOfContig(CurrentGene, CurrentSeedDirection.direction_x):
                    if CurrentSeedDirection.restartFromSeed():
                        CurrentGene, CurrentOrtholog = Seed, SeedOrtholog
                        switch = 0
                        continue
                    else:
                        break
                else:
                    switch = SL.switchFlag(switch)
                    continue

            # NextGene has no ortholog
            if not SL.geneHasOrtholog(NextGene.gene, _orthologs_, switch=switch):
                no_orthos[switch].add(NextGene.gene)
                if CurrentSeedDirection.restartFromSeed():
                    CurrentGene, CurrentOrtholog = Seed, SeedOrtholog
                    switch = 0
                    continue
                else:
                    break

            NextOrtholog = SL.setOrthologGeneCall(
                NextGene, _orthologs_, _gene_calls_, switch)
            # print(NextGene.gene)
            # print(f'NextOrtholog.gene={NextOrtholog.gene}')
            # print(f'CurrentOrtholog.gene={CurrentOrtholog.gene}')
            # print(f'Dn NextOrtho={NextOrtholog.downstream_gene}')
            # print(f'Dn CurrOtho={CurrentOrtholog.downstream_gene}')
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
                    SL.assignDirection(NextOrtholog, CurrentSeedDirection)
                    CurrentGene, CurrentOrtholog = NextGene, NextOrtholog
                    append_ = True
                elif NextOrtholog.gene is None:
                    break
                else:
                    if CurrentSeedDirection.restartFromSeed():
                        switch = 0
                        CurrentGene, CurrentOrtholog = Seed, SeedOrtholog
                    else:
                        break

            if append_:
                SL.appendBothIgnore(
                    ignore, switch, NextGene, NextOrtholog)
                SL.appendValues(values, switch, CurrentGene,
                                CurrentOrtholog, CurrentSeedDirection)
                append_ = False
        # Out of while loop
        runs += 1
        SL.recordSyntenyInStone(synteny_dic, Seed, values)
        append_ = False

    c = 0
    all_written_x = set()
    all_written_y = set()
    for seed in synteny_dic.keys():
        c += len(synteny_dic[seed][0])
        for gene in synteny_dic[seed][0]:
            all_written_x.add(gene.gene)
        for gene in synteny_dic[seed][1]:
            all_written_y.add(gene.gene)
    # print(f"Syntenic loci found: {c}")
    a_not_written = set(A_Genecalls.keys()) - all_written_x
    a_not_written.update(no_orthos[0])
    b_not_written = set(B_Genecalls.keys()) - all_written_y
    b_not_written.update(no_orthos[1])
    # print(f"Not written for a: {len(a_not_written)}")
    # print(f"Not written for b: {len(b_not_written)}")

    a_not = {}
    b_not = {}
    for gene in a_not_written:
        a_not[gene] = A_Genecalls[gene]
    for gene in b_not_written:
        b_not[gene] = B_Genecalls[gene]

    return [synteny_dic, a_not, b_not]
