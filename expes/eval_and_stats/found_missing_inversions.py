#! /bin/python3

import sys

input_files = sys.argv[1:]

true_inv = input_files[0]
all_intersect = input_files[1]
recall_intersect = input_files[2]

# List all true inversions
inv_status = dict()
inv_to_screen = list()
with open(true_inv, "r") as file:
    for line in file:
        inv_id = line.split("\t")[1]
        inv_status[inv_id] = None
        inv_to_screen.append(inv_id)

# List found inversions (passing recall requirements)
with open(recall_intersect, "r") as file:
    for line in file:
        status = [None, None]
        parsed_line = line.rstrip().split("\t")
        inv_id = parsed_line[2]

        if parsed_line[-1] != ".":
            status[0] = "Precise"
            if inv_id in inv_to_screen:
                inv_to_screen.remove(inv_id)

            if "path" in parsed_line[-1] and "aln" in parsed_line[-1]:
                status[1] = "Mixed"
            elif "path" in parsed_line[-1]:
                status[1] = "Path-explicit"
            else:
                status[1] = "Alignment-rescued"
            
        inv_status[inv_id] = status

# List imprecise inversion bubbles
with open(all_intersect, "r") as file:
    for line in file:
        status = [None, None]
        parsed_line = line.rstrip().split("\t")
        inv_id = parsed_line[1]
        overlap = int(parsed_line[-1])

        if inv_id in inv_to_screen and overlap > 0:
            status[0] = "Imprecise"
            if inv_id in inv_to_screen:
                inv_to_screen.remove(inv_id)

            if "path" in parsed_line[-1] and "aln" in parsed_line[-1]:
                status[1] = "Mixed"
            elif "path" in parsed_line[-1]:
                status[1] = "Path-explicit"
            else:
                status[1] = "Alignment-rescued"
            
            inv_status[inv_id] = status

# Label missing inversions
for inv_id in inv_to_screen:
    status = ["Unannotated", "Unannotated"]
    inv_status[inv_id] = status

# Output results
print('\t'.join(["chr", "start", "end", "size", "quality", "type"]))

with open(true_inv, "r") as file:
    for line in file:
        parsed_line = line.rstrip().split("\t")
        inv_id = parsed_line[1]
        start = int(inv_id)
        end = int(parsed_line[2])

        parsed_line.append(str(end - start + 1))
        parsed_line.extend(inv_status[inv_id])

        print('\t'.join(parsed_line))