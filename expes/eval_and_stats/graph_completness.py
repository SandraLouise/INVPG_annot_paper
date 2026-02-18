from Bio import SeqIO
from sys import argv
from os import system,remove
from os.path import exists
from pathlib import Path

_,gfa_file_path,fasta_file_path,mode = argv


# mode must be either 'cactus', 'mc' or 'pggb'


tempname: str = "temp.fa"

if not exists('rs-pancat-paths'):
    system("wget https://github.com/dubbsieg/rs-pancat-path/releases/download/0.1.1/rs-pancat-paths")
    system("chmod +x rs-pancat-paths")

# Compute path reconstruction
system(f"./rs-pancat-paths {gfa_file_path} -C > {tempname}")

print(f"----- {Path(gfa_file_path).stem} -----")
# Loading reference fasta
with open(fasta_file_path) as fa:
    # Loading graph fasta
    status:bool = False
    for fa_record in SeqIO.parse(fa, "fasta"):
        if mode == 'cactus' or mode == 'mc':
            head,particle,tail = fa_record.id.split('#')
            header = head + '#' + particle + '#' + head + '.' + particle + ('_' if tail[0].isupper() else '.') + tail + '#0'
        elif mode == 'pggb':
            header = fa_record.id
        else:
            pass
        found:bool = False
        with open(tempname) as gfa:
            for gfa_record in SeqIO.parse(gfa, "fasta"):
                if header == gfa_record.id:
                    found = True
                    if fa_record.seq.lower() != gfa_record.seq.lower():
                        status = True
                        print(f"Record {fa_record.id} not identical ")
                    break
        if not found:
            status = True
            print(f"Record {fa_record.id} not found in graph")
    if status:
        print("Graph is incomplete")
    else:
        print("Graph is complete")

#Remove temporary file
if exists(tempname):
    remove(tempname)