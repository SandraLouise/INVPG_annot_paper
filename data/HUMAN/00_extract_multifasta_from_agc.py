from os import system, remove
from os.path import exists

# Extracting sequences names
if not exists('agc-1.1_x64-linux'):
    system('curl -L https://github.com/refresh-bio/agc/releases/download/v1.1/agc-1.1_x64-linux.tar.gz|tar -zxvf - agc-1.1_x64-linux/agc')
system('agc-1.1_x64-linux/agc listset HPRC-yr1.agc > temp.txt')

#
paths: list = list()
with open('temp.txt', 'r', encoding='utf-8') as reader:
    for line in reader:
        paths.append(line.strip()[:-2])

# We get the genomes from the .tsv file
chain: str = 'NA19434	HG00096	HG00171	HG00268	HG00512	HG00513	HG00514	HG00731	HG00732	HG00733	HG00864	HG01114	HG01352	HG01505	HG01573	HG01596	HG02011	HG02018	HG02059	HG02106	HG02492	HG02587	HG02818	HG03009	HG03065	HG03125	HG03371	HG03486	HG03683	HG03732	HG04217	NA12329	NA12878	NA18534	NA18939	NA19036	NA19238	NA19239	NA19240	NA19650	NA19983	NA20509	NA20847	NA24385'
genomes: list = chain.split()

# We extract intersection of genomes
for ind in set(genomes).intersection(set(paths)):
    for hap in ['1', '2']:
        if not exists(f'HPP_fastas/{ind}.{hap}.fa'):
            system(
                f'agc-1.1_x64-linux/agc getset HPRC-yr1.agc {ind}.{hap} > {ind}.{hap}.fa'
            )

# We clean temp files
remove('temp.txt')
