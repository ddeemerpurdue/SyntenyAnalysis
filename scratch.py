import os
import sys

# summary = 'rawdata/All-PS-Summary.txt'
# gff = 'rawdata/GFF-Internal/SupernatantGeneCalls.txt'

# mydic = {}

# with open(summary) as file:
#     line = file.readline()
#     line = file.readline()
#     while line:
#         bin_ = line.split('\t')[3]
#         if bin_.endswith(('S', 'SM')):
#             gene_id = line.split('\t')[4]
#             mydic[gene_id] = bin_
#         line = file.readline()

# with open(gff) as file:
#     line = file.readline()
#     line = file.readline()
#     while line:
#         gene_id = line.split('\t')[0]
#         if gene_id in mydic:
#             output = f"rawdata/GFF-Internal/{mydic[gene_id]}.gff"
#             with open(output, 'a') as out:
#                 out.write(line)
#         line = file.readline()

# # val = 'dane_greg_deemer_allen'
# # newval = val[0:-1] + 'X'
# # print(newval)
# # import os
# # if os.path
# # # if __name__ == '__main__':
# # #     summary_file = 'rawdata/All-PS-Summary.txt'

# # #     orth_prefix = 'rawdata/Orthologs-Acidaminococcaceae/Orthologues_111_Acidaminococcaceae_Phascolarctobacterium_SM/'
# # #     orth_file = '111_Acidaminococcaceae_Phascolarctobacterium_SM__v__111_Acidaminococcaceae_Phascolarctobacterium-SM-GCA_014648015.1.tsv'
# # #     ortholog_file = os.path.join(orth_prefix, orth_file)

# # #     gffA_file = 'rawdata/GFF-Internal/111_Acidaminococcaceae_Phascolarctobacterium_SM.gff '
# # #     gffB_file = 'rawdata/GFF-Acidaminococcaceae/111_Acidaminococcaceae_Phascolarctobacterium-S*-GCA_014648015.1.gff'

# # #     summary = mineSummaryFile(
# # #         summary_file, '111_Acidaminococcaceae_Phascolarctobacterium_SM')
# # #     synteny_dic = traverseSynteny(summary_file, gffA_file, gffB_file,
# # #                                   ortholog_file, Gb=gffB_file)

# # #     output = f'Output/test.txt'
# # #     with open(output, 'w') as out:
# # #         out.write(f"GenomeA\tGenomeB\n")
# # #         for seed in synteny_dic:
# # #             for gene_a, gene_b in zip(synteny_dic[seed][0], synteny_dic[seed][1]):
# # #                 if gene_a in ['-', '+']:
# # #                     vals = '\t'.join(['', '', '', '', '', ''])
# # #                 else:
# # #                     vals = summary[gene_a]
# # #                     vals = '\t'.join(vals)
# # #                 out.write(f"{gene_a}\t{gene_b}\t{vals}\n")
# # #             out.write(f"X\tX\n")


mystring = 'dane_gregory'


# import glob
#     # family = 'Acidaminococcaceae'
#     family = 'Desulfovibrionaceae'
#     files = glob.glob(f'rawdata/Orthologs-{family}/*')
#     mybins = []
#     for file in files:
#         bin_name = os.path.basename(file)
#         if any([string in bin_name for string in ['GCA', 'GCF']]):
#             pass
#         else:
#             mybins.append(bin_name.strip('Orthologues_'))
#     summary_file = 'rawdata/All-PS-Summary.txt'

#     print('Starting master analysis...')
#     for bin_ in mybins:
#         print(f"Analyzing data for bin: {bin_}")
#         orth_prefix = f'rawdata/Orthologs-{family}/Orthologues_{bin_}/'
#         print(f"\tOrtholog prefix: {orth_prefix}")
#         for ortholog_file in glob.glob(f"{orth_prefix}/*"):
#             print(f"\t\tMining ortho file: {os.path.basename(ortholog_file)}")
#             bin_b = os.path.basename(ortholog_file)
#             bin_b = bin_b.split('__v__')[1].strip('.tsv')
#             print(f'\t\tAnalyzing {bin_} vs {bin_b}')

#             gffA_file = f'rawdata/GFF-Internal/{bin_}.gff'
#             if any([string in bin_b for string in ['GCA', 'GCF']]):
#                 gffB_file = glob.glob(f'rawdata/GFF-{family}/{bin_b}*.gff')
#                 assert len(
#                     gffB_file) == 1, 'Invalid gffB_file pattern matching!'
#                 gffB_file = gffB_file[0]
#                 Gb_flag = gffB_file
#             else:
#                 gffB_file = f'rawdata/GFF-Internal/{bin_b}.gff'
#                 Gb_flag = False
#             print(f"\t\t\tGff a: {gffA_file}")
#             print(f"\t\t\tGff b: {gffB_file}")

#             if os.path.exists(f"Output/{family}"):
#                 pass
#             else:
#                 os.makedirs(f"Output/{family}")

#             uniq_name = os.path.basename(ortholog_file)
#             output = f'Output/{family}/{uniq_name}'
#             print(f"\t\t\t\tWriting analysis results to file: {output}\n\n")

#             summary = mineSummaryFile(summary_file, bin_)
#             synteny_dic = traverseSynteny(summary_file, gffA_file, gffB_file,
#                                           ortholog_file, Gb=Gb_flag)
#             with open(output, 'w') as out:
#                 out.write(f"{bin_}\t{bin_b}\n")
#                 for seed in synteny_dic:
#                     for gene_a, gene_b in zip(synteny_dic[seed][0], synteny_dic[seed][1]):
#                         if gene_a in ['-', '+']:
#                             vals = '\t'.join(['', '', '', '', '', ''])
#                         else:
#                             vals = summary[gene_a]
#                             vals = '\t'.join(vals)
#                         out.write(f"{gene_a}\t{gene_b}\t{vals}\n")
#                     out.write(f"X\tX\n")
