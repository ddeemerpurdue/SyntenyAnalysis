import glob
import os

mypath = '../Anvio-MAGs/GFF3-Final/'

for file in glob.glob(f"{mypath}*"):
    print(file)
    family = file.split('-')[-1]
    print(family)
    newlocation = f'{family}/rawdata/GFF3/'
    command = f'cp {file}/*.gff3 {newlocation}'
    print(command)
    #os.system(command)
