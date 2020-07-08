#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 16:49:47 2016

Seg Stitch

@author: dpenaherrera
"""
import sys

# ==============================================================================
# WHERE ALL THE MAGIC HAPPENS
# Grab all .seg files in path and convert them to .bed files


def seg_extend(file):
    with open(file, 'r') as f:
        seglist = f.readlines()

    # Parse begins
    all_lines = []

    for line in seglist:
        all_lines.append(line.rstrip().split('\t'))

    del seglist
    header = all_lines.pop(0)

    # Set these to integer
    #       STARTEND
    #       NUM_POINTS_COPY_RATIO
    #       NUM_POINTS_ALLELE_FRACTION
    #   and these to float
    #       LOG2_COPY_RATIO_POSTERIOR_10
    #       LOG2_COPY_RATIO_POSTERIOR_50
    #       LOG2_COPY_RATIO_POSTERIOR_90
    #       MINOR_ALLELE_FRACTION_POSTERIOR_10
    #       MINOR_ALLELE_FRACTION_POSTERIOR_50
    #       MINOR_ALLELE_FRACTION_POSTERIOR_90
    for idx in range(len(all_lines)):
        all_lines[idx][1:5] = [int(x) for x in all_lines[idx][1:5]]
        all_lines[idx][5:] = [float(x) for x in all_lines[idx][5:]]

# -------FINDING SEGMENTS THAT CROSS THE CENTROMERE-----------------------------
    for chrom in chromes:
        i = 1
        while i in range(len(all_lines)):
            line = all_lines[i]
            if line[0] == chrom:
                if (int(line[1]) < centromeres[chrom][0]) and (int(line[2]) > centromeres[chrom][1]):
                    newline = [line[0], int(centromeres[chrom][1]), line[2], line[3], line[4], line[5],
                               line[6], line[7], line[8], line[9], line[10]]
                    all_lines.insert(i+1, newline)
                    all_lines[i][2] = centromeres[chrom][0]
                    del newline
                    i += 1
                else:
                    i += 1
            else:
                i += 1

# ---STITCHING SEMGENTS---------------------------------------------------------
    for k in range(len(all_lines)-1):
        for chrom in chromes:
            if (all_lines[k][0] == chrom) and (all_lines[k+1][0] == chrom):
                if (all_lines[k][2] != centromeres[chrom][0]) and (all_lines[k+1][1] != centromeres[chrom][1]):
                    all_lines[k][2] = round((int(all_lines[k][2])+int(all_lines[k+1][1]))/2)
                    all_lines[k+1][1] = all_lines[k][2]
                elif (all_lines[k][2] == centromeres[chrom][0]) and (all_lines[k+1][1] == centromeres[chrom][1]):
                    continue

# ---WRITE THE OUTPUTS TO FILE--------------------------------------------------
    for idx in range(len(all_lines)):
        all_lines[idx][:] = [str(x) for x in all_lines[idx][:]]

    with open(file, 'w') as f:
        f.write("\t".join(elem for elem in header) + "\n")

        for line in all_lines:
            f.write("\t".join(elem for elem in line)+"\n")


# ==============================================================================
if __name__ == '__main__':
    centromeres = {}

    with open(sys.argv[1]) as centro:
        for interval in centro:
            (key, start, stop) = interval.rstrip().split('\t')
            centromeres[key] = [int(start), int(stop)]

    chromes = list(centromeres.keys())

    seg_extend(sys.argv[2])
