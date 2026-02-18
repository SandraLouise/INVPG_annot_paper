#!/usr/bin/env python3
from BubbleGun.Graph import Graph as BubbleGraph
from BubbleGun.find_bubbles import find_bubbles
from gfagraphs import Graph as GfaGraph, Orientation
from sys import argv
from os import mkdir, system, path


def revcomp(string: str, compl: dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}) -> str:
    """Tries to compute the reverse complement of a sequence

    Args:
        string (str): original character set
        compl (dict, optional): dict of correspondances. Defaults to {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}.

    Raises:
        IndexError: Happens if revcomp encounters a char that is not in the dict

    Returns:
        str: the reverse-complemented string
    """
    try:
        return ''.join([compl[s] for s in string][::-1])
    except IndexError as exc:
        raise IndexError(
            "Complementarity does not include all chars in sequence.") from exc


"""
This script aims to output a VCF-like file from the output of BubbleGun.

Steps:
1. Read the input GFA file with BubbleGun
2. Find the bubbles in the graph using BubbleGun
3. From node lists, extract the subgraphs representing the bubbles with odgi
4. Find in each subgraph the paths that represent the bubble
5. Select the subgraphs (bubbles) that contains the reference path
6. For each selected bubble, output a VCF line with the information of the bubble
"""

# argv[1] = input GFA file
# argv[2] = reference path

# STEP 1
graph: BubbleGraph = BubbleGraph(argv[1])
reference: str = argv[2]
temp_folder: str = 'temp/'

if not path.exists(temp_folder):
    mkdir(temp_folder)

# STEP 2
find_bubbles(graph)

# STEP 3
for i, val in enumerate(graph.bubbles.values()):
    output_bgf: str = f"{temp_folder}/bubbles_{i}.bgf"
    with open(output_bgf, 'a', encoding='utf-8') as gfa_writer:
        gfa_writer.write(
            '\n'.join(val.list_bubble())
        )
    # STEP 4
    output_og: str = f"{temp_folder}/odgiex_{i}.og"
    system(
        f"odgi extract -i {argv[1]} -o {output_og} --node-list {output_bgf}"
    )

# STEP 5
print('##fileformat=VCFv4.2')
print('##source=BubbleGun')
print('##reference=' + reference)
print(f'##total_bubbles={len(graph.bubbles)}')
print("#ref_name\tstart_pos\tx\tref_seq\talt_seq\tx\tx\tinfo_field\tx\tx")
for i in range(len(graph.bubbles)):
    output_og: str = f"{temp_folder}/odgiex_{i}.og"
    output_gfa: str = f"{temp_folder}/gfaview_{i}.gfa"
    system(
        f"odgi view -g -i {output_og} > {output_gfa}"
    )
    gfa_graph: GfaGraph = GfaGraph(
        gfa_file=output_gfa,
        with_sequence=True,
    )

    # STEP 6
    if any(bool_vector := [x.startswith(reference) for x in list(gfa_graph.paths.keys())]):
        reference_name: str = list(gfa_graph.paths.keys())[
            bool_vector.index(True)]
        reference_path: str = ''.join(
            [f'>{x}' if ori == Orientation.FORWARD else f'<{x}' for x,
                ori in gfa_graph.paths[reference_name]['path']]
        )
        reference_allele: str = ''.join(
            [gfa_graph.segments[x]['seq'] if ori == Orientation.FORWARD else revcomp(gfa_graph.segments[x]['seq']) for x,
                ori in gfa_graph.paths[reference_name]['path']]
        )
        alternates_paths: str = ','.join(
            [
                ''.join(
                    [f'>{x}' if ori == Orientation.FORWARD else f'<{x}' for x,
                        ori in gfa_graph.paths[path]['path']]
                ) for path in gfa_graph.paths.keys() if path != reference_name
            ]
        )
        alternates_alleles: str = ','.join(
            [
                ''.join(
                    [gfa_graph.segments[x]['seq'] if ori == Orientation.FORWARD else revcomp(gfa_graph.segments[x]['seq']) for x,
                        ori in gfa_graph.paths[path]['path']]
                ) for path in gfa_graph.paths.keys() if path != reference_name
            ]
        )
        start_pos: str = reference_name.split(':')[1].split('-')[0]

        print(
            reference_name,
            start_pos,
            '.',
            reference_allele,
            alternates_alleles,
            '.',
            '.',
            'AT='+reference_path+','+alternates_paths,
            '.',
            '.',
            sep='\t',
        )

# Clean up
system(f"rm -r {temp_folder}")
