# annotations = ['COG20_CATEGORY', 'COG20_CAT_ACC', 'COG20_PATHWAY',
#                'COG20_PATHWAY_ACC', 'COG20_FUNCTION', 'COG20_FUNCTION_ACC',
#                'KEGG_Module', 'KEGG_Module_ACC', 'KEGG_Class', 'KEGG_Class_ACC',
#                'KOfam', 'KOfam_ACC', 'FigFams', 'FigFams_ACC', 'Pfam', 'Pfam_ACC'
#                'TIGRFAM', 'TIGRFAM_ACC']

import glob
import os
import sys


# def zipDirectories(mypath):
#     family = mypath.split('/')[0]
#     print(family)
#     for item in glob.glob(mypath + '*/*'):
#         print(item)
#         if os.path.isdir(item):
#             save_directory = '/'.join(item.split('/')[0:-1])
#             print(f'Save directory: {save_directory}')
#             basename = os.path.basename(item)
#             command = f'tar -czf {save_directory}/{family}_{basename}.tar.gz {item}'
#             print(f'Command: {command}')
#             os.system(command)


# def deleteDirectories(mypath):
#     for item in glob.glob(mypath + '*/*'):
#         if os.path.isdir(item):
#             print(f'Deleting: {item}')
#             delete_command = f'rm -r {item}'
#             os.system(delete_command)


# if __name__ == '__main__':
#     zipDirectories(sys.argv[1])
#     deleteDirectories(sys.argv[1])
print('hello', 'dane', end='')
