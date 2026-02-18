#!/usr/bin/env python3
import re
from collections import defaultdict
from sys import argv

"""
This script detects one-node inversions in a gfa file by comparing 2 paths, the output (std out) is a bed file with inversion coordinates (similar format as inv-pg-annot)

2 input arguments: gfa_file and chromosome_name
Usage: python detect_one_node_inversions_in_gfa_paths.py mgc.gfa scaffold_6 > mgc_inversion_1node.bed

Note: limited to 2 paths
Note: looks for this signal x+,y+,z+ in one path and x+,y-,z+ in the other path ; (or can be x-,y-,z- / x-,y+,z- ; or even x+,y+,x- / x+,y-,x-) : the requirement = first and third nodes with identical id and strand in both paths and middle node same id but strand different ; 
Note: the position returned corresponds to the starting position of y in path 1 (first path encountered in the gfa file P lines)  
Note: does not look for signals x+,y+,z+ / z-,y+,x-  which should be valid 1node-inversions too (inversion inside another inversion)...

Steps:
1. Read the input GFA file and extract paths as sequences of nodes and a dictionary of node sizes 
2. Index all nodes in both paths + all minus nodes
3. Find common nodes that appear as '-' in either path
4. Check node triplets around these nodes
5. Computes the positions of inversions and output in bed format
"""

def get_gfa_paths(file_path):
    """
    Extracts and returns a list of path lines from a GFA file + a dict with node sizes

    :param file_path: Path to the GFA file.
    :return: a tuple of a list of path lines as strings and a dictionary of node sizes
    """
    path_sequences = []
    node_sizes = defaultdict(int)
    
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("S"):
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    node_id = parts[1]
                    sequence = parts[2]
                    node_sizes[node_id] = len(sequence)
            if line.startswith('P'):
                parts = line.strip().split('\t')  # Split by tab
                if len(parts) >= 3:  # Ensure correct format
                    segments = parts[2].split(',')  # Extract sequence and split by comma
                    parsed_segments = [re.match(r'(\w+)([+-])', seg).groups() for seg in segments]
                    path_sequences.append(parsed_segments)   
    return path_sequences, node_sizes

def index_nodes_and_minus_nodes(path):
    """
    Indexes all positions of each node in the path.

    :param path: List of [name, orientation] segments.
    :return: Dictionary mapping node names to their indices.
    """
    node_positions = defaultdict(list)
    minus_nodes = defaultdict(list)

    for i, (node, orientation) in enumerate(path):
        node_positions[node].append(i)
        if orientation == '-':
            minus_nodes[node].append(i)
    return node_positions, minus_nodes

def get_triplet(path, index):
    """
    Extracts a triplet centered around a given index.

    :param path: List of [name, orientation] segments.
    :param index: Center node index.
    :return: Tuple of (prev, curr, next) or None if out of bounds.
    """
    if 0 < index < len(path) - 1:
        return (tuple(path[index - 1]), tuple(path[index]), tuple(path[index + 1]))
    return None


#Arguments
# argv[1] = input GFA file
# argv[2] = reference name to put in output bed

graph_file = argv[1]
ref_name = argv[2]

# Step1: Read the input GFA file and extract paths as sequences of nodes and a dictionary of node sizes 
path_sequences, node_sizes = get_gfa_paths(graph_file)

path1 = path_sequences[0]
path2 = path_sequences[1]


# Step 2: Index all nodes in both paths + all minus nodes
node_positions_1, minus_nodes_1 = index_nodes_and_minus_nodes(path1)
node_positions_2, minus_nodes_2 = index_nodes_and_minus_nodes(path2)
    
# Step 3: Find common nodes that appear as '-' in either path
possible_inversion_nodes = set(minus_nodes_1.keys()) ^ set(minus_nodes_2.keys())
#print(len(possible_inversion_nodes))

    
# Step 4: Check node triplets around these nodes
inversions = []
for node in possible_inversion_nodes:
    positions1 = node_positions_1[node]
    positions2 = node_positions_2[node]

    for pos1 in positions1:
        triplet1 = get_triplet(path1, pos1)
        if not triplet1:
            continue  # Skip if triplet is incomplete

        for pos2 in positions2:
            triplet2 = get_triplet(path2, pos2)
            if not triplet2:
                continue

            # Compare triplets: first & last must match, middle must flip
            if (triplet1[0] == triplet2[0] and  # First node same (name and orientation)
                triplet1[2] == triplet2[2] and  # Last node same
                triplet1[1][0] == triplet2[1][0] and  # Middle node name same
                triplet1[1][1] != triplet2[1][1]):  # Middle node orientation flips
                
                inversions.append([pos1,triplet1])  # Store from path1
    
# Step 5: computes the positions of inversions
#bed = []
for position1, triplet in inversions:
    begin = 0
    for index in range(position1):
        begin += node_sizes[path1[index][0]]
    inversion_size = node_sizes[triplet[1][0]]
    #bed.append([begin, inversion_size])
    print(f"{ref_name}\t{begin}\t{begin+inversion_size}\tINV:path:1.0,{triplet[0][0]}:{triplet[2][0]}")

