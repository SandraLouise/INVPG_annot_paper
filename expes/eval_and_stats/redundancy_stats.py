#! /bin/python3

import sys

intersect = sys.argv[1]
recall_intersect = sys.argv[2]

# Set indexes
CHROM_ID = 0
INV_START = 1
INV_END = 2

with open(intersect, "r") as file:
    for line in file:
        if int(line.strip('\n').split("\t")[-1]) > 0:
            parsedline = line.strip('\n').split("\t")
            field_count = len(parsedline)
            for i in range(INV_END+1, field_count):
                if parsedline[i] == parsedline[CHROM_ID]:
                    BUBBLE_OFFSET = i
            break

BUBBLE_START = BUBBLE_OFFSET + 1
BUBBLE_END = BUBBLE_OFFSET + 2
OVERLAP = field_count - 1

# Save all_intersect in dict
d_intersect = dict()

with open(intersect, "r") as file:
    for line in file:

        parsed_line = line.rstrip().split("\t")

        if parsed_line[BUBBLE_START] == "-1":
            continue

        if parsed_line[INV_START] not in d_intersect.keys():
            d_intersect[parsed_line[INV_START]] = list()
        
        d_intersect[parsed_line[INV_START]].append(parsed_line)

# Save recall_intersect in list (save the bubble start)
s_recall_intersect = set()

with open(recall_intersect, "r") as file:
    for line in file:

        parsed_line = line.rstrip().split("\t")

        # Ignore inversions without recall-quality bubble
        if parsed_line[BUBBLE_START] == "-1":
            continue
        
        # Save bubble start in recall list
        s_recall_intersect.add(parsed_line[BUBBLE_START])

# Divide annotated bubbles
nonRedundant_TP = set()
redundant = set()
imprecise = set()

for inv in d_intersect.keys():

    ## Redundant bubbles
    if len(d_intersect[inv]) > 1:

        for entry in d_intersect[inv]:
            redundant.add(entry[BUBBLE_START])

            # CORRECTION: redundant can also be imprecise
            ## Redundant imprecise
            if entry[BUBBLE_START] not in s_recall_intersect:
                imprecise.add(entry[BUBBLE_START])

    elif len(d_intersect[inv]) == 1:
        entry = d_intersect[inv][0]

        ## Precise AND non redundant (entry[BUBBLE_START] should be in s_recall_intersect)
        if entry[BUBBLE_START] in s_recall_intersect:
            nonRedundant_TP.add(entry[BUBBLE_START])

        ## Non redundant imprecise
        else:
            imprecise.add(entry[BUBBLE_START])

# Output stats

nonRedundant_TP = sorted(list(nonRedundant_TP))
redundant = sorted(list(redundant))
imprecise = sorted(list(imprecise))

for entry in nonRedundant_TP:

    if entry in redundant:
        redundant.remove(entry)
    if entry in imprecise:
        imprecise.remove(entry)

print(f"NonRedundant_bubbles\t{str(len(nonRedundant_TP))}")
print(f"Redundant_bubbles\t{str(len(redundant))}")
print(f"Imprecise_bubbles\t{str(len(imprecise))}")