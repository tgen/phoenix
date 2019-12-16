#!/usr/bin/env python3

# MYC blacklist region chr8:129571461-129571635 (b37)
# MYC blacklist region chr8:128559215-128559389 (b38)

#################################################################
# Configure the enviroment
import pandas as pd
import argparse
from multiprocessing import Pool

#################################################################

parser = argparse.ArgumentParser(description='Process pairoscope discordant reads to make myeloma Ig translocation calls.')

parser.add_argument('-i', '--input_file',
                    required=True,
                    metavar='File.tsv',
                    dest="input_file",
                    help='Discordant_read_table')
parser.add_argument('-s', '--specimen',
                    required=True,
                    help='Specimen Name, must match discordant read table')
parser.add_argument('-o', '--output_file',
                    required=True,
                    metavar='File.tsv',
                    help='Output filename')
parser.add_argument('-w', '--window',
                    default=2000,
                    type=int,
                    metavar='INT',
                    help='Genomic window size to query (recommend 2.5x insert size)')
parser.add_argument('-m', '--window_min',
                    default=100,
                    type=int,
                    metavar='INT',
                    help='Required discordant read island size (recommend 0.1x insert size)')
parser.add_argument('-c', '--call_requirement',
                    default=3,
                    type=int,
                    metavar='INT',
                    help='Required number of discordant read pairs meeting requirements to define a call (recommend >= 3)')

args = parser.parse_args()

#################################################################
# Capture Run Parameters

window = args.window
call_requirement = args.call_requirement
window_min = args.window_min
sample = args.specimen

#################################################################
# Read in the discordant reads table

dreads = pd.read_csv(args.input_file, sep="\t")

#################################################################
# Define Functions

def call_translocations(sample, gene, gene_chr, order, window_start=0, window_end=300000000):
    # initialize variables
    gene_window = 0
    gene_maxwindowcount = 0
    gene_nextlargestwindowcount = 0
    gene_windowend = 0
    gene_maxwindowwidth = 0
    gene_maxwindowlocation = 0
    gene_nextwindowwidth = 0
    gene_nextwindowlocation = 0
    gene_ighsource = 0

    # Debug
    print('')
    print('------------------------------')
    print(sample)
    print(gene)
    print(gene_chr)
    print(order)

    # Determine which immunoglobulin loci is the most likely to contain a translocation
    # Watch ChrA and ChrB as they are in numeric order so dependign on Tx and igH region partner
    # and target can switch columns
    if order == "standard":
        print('1 - In Standard Loop')
        if gene != 'MYC':
            print('1 - In MYC Sub-loop')
            gene_igh = dreads[(dreads.Specimen == sample) & (dreads.ChrA == gene_chr) & (dreads.ChrB == "chr14")
                              & (dreads.PositionA >= window_start) & (dreads.PositionA <= window_end)]

            count_gene_igh = len(gene_igh.index)

            gene_igk = dreads[(dreads.Specimen == sample) & (dreads.ChrA == "chr2") & (dreads.ChrB == gene_chr)
                              & (dreads.PositionB >= window_start) & (dreads.PositionB <= window_end)]
            count_gene_igk = len(gene_igk.index)
            gene_igl = dreads[(dreads.Specimen == sample) & (dreads.ChrA == gene_chr) & (dreads.ChrB == "chr22")
                              & (dreads.PositionA >= window_start) & (dreads.PositionA <= window_end)]
            count_gene_igl = len(gene_igl.index)
        else:
            print('1 - In Non-MYC Sub-loop')
            gene_igh = dreads[(dreads.Specimen == sample) & (dreads.ChrA == gene_chr) & (dreads.ChrB == "chr14")
                              & (dreads.PositionA >= window_start) & (dreads.PositionA <= window_end)
                              & ((dreads.PositionA <= 128559215) | (dreads.PositionA >= 128559389))]
            count_gene_igh = len(gene_igh.index)
            gene_igk = dreads[(dreads.Specimen == sample) & (dreads.ChrA == "chr2") & (dreads.ChrB == gene_chr)
                              & (dreads.PositionB >= window_start) & (dreads.PositionB <= window_end)
                              & ((dreads.PositionB <= 128559215) | (dreads.PositionB >= 128559389))]
            count_gene_igk = len(gene_igk.index)
            gene_igl = dreads[(dreads.Specimen == sample) & (dreads.ChrA == gene_chr) & (dreads.ChrB == "chr22")
                              & (dreads.PositionA >= window_start) & (dreads.PositionA <= window_end)
                              & ((dreads.PositionA <= 128559215) | (dreads.PositionA >= 128559389))]
            count_gene_igl = len(gene_igl.index)
        print(count_gene_igh)
        print(count_gene_igk)
        print(count_gene_igl)
        if count_gene_igh >= count_gene_igk and count_gene_igh >= count_gene_igl:
            print('2-In IgH Loop')
            # Suspect igH Translocation as igH count greater then igK and igL
            # if all three have zeros or the same number we default to the igH)
            table_gene = gene_igh
            count_gene = len(table_gene.index)
            gene_ighsource = 1
        elif count_gene_igk > count_gene_igh and count_gene_igk > count_gene_igl:
            print('2-In IgK Loop')
            # Suspect igK Translocation
            # WARNiNG THE CHROMOSOME ORDER BECOMES UNEXPECTED, HOW TO FiX?
            # Create Flipped Table for the kappa counts (	ChrA	PositionA	ChrB	PositionB)
            # Update the column headers
            new_cols = ['Speciment', 'ChrB', 'PositionB', 'ChrA', 'PositionA']
            gene_igk.columns = new_cols
            # Create new GENE table with the relabelled columns in the expected order
            cols = ['Speciment', 'ChrA', 'PositionA', 'ChrB', 'PositionB']
            table_gene = gene_igk[cols]
            count_gene = len(table_gene.index)
            gene_ighsource = 2
        elif count_gene_igl > count_gene_igh and count_gene_igl > count_gene_igk:
            print('2-In IgL Loop')
            # Suspect igL Translocation
            table_gene = gene_igl
            count_gene = len(table_gene.index)
            gene_ighsource = 3
        else:
            # Error capture NOT AN EXPECTED EVENT
            # Could be 0 in all three potentially, then what?
            # What if two but not all three have the same number of counts?
            print('2 - Argh - WHAT CAUSES THiS TO HAPPEN')
            if count_gene_igh < 3 and count_gene_igk < 3 and count_gene_igl < 3:
                print('DEFAULTiNG TO iGH COUNTS BECAUSE ALL LESS THAN 3')
                table_gene = gene_igh
                count_gene = len(table_gene.index)
                gene_ighsource = 1
            elif count_gene_igh < count_gene_igk and count_gene_igh < count_gene_igl and count_gene_igk == count_gene_igl:
                print('DEFAULTiNG TO iGH COUNTS BECAUSE IgK and IgL ARE EQUAL EVEN THOUGH LESS')
                table_gene = gene_igh
                count_gene = len(table_gene.index)
                gene_ighsource = 1

    elif order == "reverse":
        print('1 - In Reverse Loop')
        gene_igh = dreads[(dreads.Specimen == sample) & (dreads.ChrA == "chr14") & (dreads.ChrB == gene_chr)
                          & (dreads.PositionB >= window_start) & (dreads.PositionB <= window_end)]
        count_gene_igh = len(gene_igh.index)
        gene_igk = dreads[(dreads.Specimen == sample) & (dreads.ChrA == "chr2") & (dreads.ChrB == gene_chr)
                          & (dreads.PositionB >= window_start) & (dreads.PositionB <= window_end)]
        count_gene_igk = len(gene_igk.index)
        gene_igl = dreads[(dreads.Specimen == sample) & (dreads.ChrA == gene_chr) & (dreads.ChrB == "chr22")
                          & (dreads.PositionA >= window_start) & (dreads.PositionA <= window_end)]
        count_gene_igl = len(gene_igl.index)

        print(count_gene_igh)
        print(count_gene_igk)
        print(count_gene_igl)

        if count_gene_igh >= count_gene_igk and count_gene_igh >= count_gene_igl:
            print('2-In IgH Loop')
            # Suspect igH Translocation as igH count greater then igK and igL
            # if all three have zeros or the same number we default to the igH
            # Update the column headers
            new_cols = ['Speciment', 'ChrB', 'PositionB', 'ChrA', 'PositionA']
            gene_igh.columns = new_cols
            # Create new GENE table with the relabelled columns in the expected order
            cols = ['Speciment', 'ChrA', 'PositionA', 'ChrB', 'PositionB']
            table_gene = gene_igh[cols]
            count_gene = len(table_gene.index)
            gene_ighsource = 1
        elif count_gene_igk > count_gene_igh and count_gene_igk > count_gene_igl:
            print('2-In IgK Loop')
            # Suspect igK Translocation
            # WARNiNG THE CHROMOSOME ORDER BECOMES UNEXPECTED, HOW TO FiX?
            # Create Flipped Table for the kappa counts (	ChrA	PositionA	ChrB	PositionB)
            # Update the column headers
            new_cols = ['Speciment', 'ChrB', 'PositionB', 'ChrA', 'PositionA']
            gene_igk.columns = new_cols
            # Create new GENE table with the relabelled columns in the expected order
            cols = ['Speciment', 'ChrA', 'PositionA', 'ChrB', 'PositionB']
            table_gene = gene_igk[cols]
            count_gene = len(table_gene.index)
            gene_ighsource = 2
        elif count_gene_igl > count_gene_igh and count_gene_igl > count_gene_igk:
            print('2-In IgL Loop')
            # Suspect igL Translocation
            table_gene = gene_igl
            count_gene = len(table_gene.index)
            gene_ighsource = 3
        else:
            # Error capture NOT AN EXPECTED EVENT
            # Could be 0 in all three potentially, then what?
            # What if two but not all three have the same number of counts?
            print('2 - Argh - WHAT CAUSES THiS TO HAPPEN')
            if count_gene_igh < 3 and count_gene_igk < 3 and count_gene_igl < 3:
                print('DEFAULTiNG TO iGH COUNTS BECAUSE ALL LESS THAN 3')
                # Update the column headers
                new_cols = ['Speciment', 'ChrB', 'PositionB', 'ChrA', 'PositionA']
                gene_igh.columns = new_cols
                # Create new GENE table with the relabelled columns in the expected order
                cols = ['Speciment', 'ChrA', 'PositionA', 'ChrB', 'PositionB']
                table_gene = gene_igh[cols]
                count_gene = len(table_gene.index)
                gene_ighsource = 1
            elif count_gene_igh < count_gene_igk and count_gene_igh < count_gene_igl and count_gene_igk == count_gene_igl:
                print('DEFAULTiNG TO iGH COUNTS BECAUSE IgK and IgL ARE EQUAL EVEN THOUGH LESS')
                # Update the column headers
                new_cols = ['Speciment', 'ChrB', 'PositionB', 'ChrA', 'PositionA']
                gene_igh.columns = new_cols
                # Create new GENE table with the relabelled columns in the expected order
                cols = ['Speciment', 'ChrA', 'PositionA', 'ChrB', 'PositionB']
                table_gene = gene_igh[cols]
                count_gene = len(table_gene.index)
                gene_ighsource = 1
    # Do Stuff
    if count_gene > 1:
        # Sort the table by PositionA
        gene_sorted = table_gene.sort_values(['PositionA'], ascending=True, inplace=False, kind='quicksort', na_position='last')
        gene_top = gene_sorted.head(1)
        gene_bottom = gene_sorted.tail(1)
        gene_top_position = gene_top.iat[0, 2]
        gene_bottom_position = gene_bottom.iat[0, 2]
        gene_window = gene_bottom_position - gene_top_position
        # There are likely some true breakpoints with two breakpoint clusters and these create large windows
        # To account for this issue we need to count discordant reads per window
        for row in gene_sorted.index:
            # Get PositionA value, calculate +/- window, capture rows within the +/- window
            position_current = gene_sorted.at[row, 'PositionA']
            position_negative = (position_current - window)
            position_positive = (position_current + window)
            position_table = gene_sorted[
                (gene_sorted.PositionA > position_negative) & (gene_sorted.PositionA < position_positive)]
            position_count = len(position_table.index)
            # Create Tracking Variables
            # if the current window cound is greater than existing Max, update Max and Next accordingly
            if position_count > gene_maxwindowcount:
                gene_nextlargestwindowcount = gene_maxwindowcount
                gene_maxwindowcount = position_count
                # NEW FEATURE START
                gene_nextwindowwidth = gene_maxwindowwidth
                gene_nextwindowlocation = gene_maxwindowlocation
                gene_max_window_top = position_table.head(1)
                gene_max_window_bottom = position_table.tail(1)
                gene_max_window_top_position = gene_max_window_top.iat[0, 2]
                gene_max_window_bottom_position = gene_max_window_bottom.iat[0, 2]
                gene_maxwindowwidth = gene_max_window_bottom_position - gene_max_window_top_position
                gene_maxwindowlocation = (gene_max_window_bottom_position + gene_max_window_top_position) / 2
                # NEW FEATURE END
                # if the two bundle counts will be greater than total count reset the NextLargest to 0
                if (gene_nextlargestwindowcount + gene_maxwindowcount) > count_gene:
                    gene_nextlargestwindowcount = 0
                    # NEW FEATURE START
                    gene_nextwindowwidth = 0
                    gene_nextwindowlocation = 0
                    # NEW FEATURE END
                gene_windowend = position_positive
            # if the second bundle in order is greater than the current next but less than the max then capture current as a next
            elif position_count > gene_nextlargestwindowcount and position_count <= gene_maxwindowcount and position_negative >= gene_windowend:
                gene_nextlargestwindowcount = position_count
                # NEW FEATURE START
                # Reusing the Max calculation variables BUT only to calculate the Next window and location (KiSS)
                gene_max_window_top = position_table.head(1)
                gene_max_window_bottom = position_table.tail(1)
                gene_max_window_top_position = gene_max_window_top.iat[0, 2]
                gene_max_window_bottom_position = gene_max_window_bottom.iat[0, 2]
                gene_nextwindowwidth = gene_max_window_bottom_position - gene_max_window_top_position
                gene_nextwindowlocation = (gene_max_window_bottom_position + gene_max_window_top_position) / 2
                # NEW FEATURE END

    # Make result calls
    if gene_maxwindowcount >= call_requirement and gene_maxwindowwidth >= window_min and gene_nextlargestwindowcount >= call_requirement and gene_nextwindowwidth >= window_min:
        gene_call = 1
        gene_bundle_count = 2
    elif gene_maxwindowcount >= call_requirement and gene_maxwindowwidth >= window_min and (
            gene_nextlargestwindowcount < call_requirement or gene_nextwindowwidth < window_min):
        gene_call = 1
        gene_bundle_count = 1
        # remove the width and location calculations from the 1sies and 2sies (below call requirement)
        gene_nextwindowwidth = 0
        gene_nextwindowlocation = 0
    else:
        gene_call = 0
        gene_bundle_count = 0
        gene_maxwindowwidth = 0
        # remove the width and location calculations from the 1sies and 2sies (below call requirement)
        gene_maxwindowlocation = 0
        gene_nextwindowwidth = 0
        gene_nextwindowlocation = 0
    print('Window Count = ' + str(count_gene))
    return (count_gene, gene_window, gene_maxwindowcount, gene_maxwindowwidth, gene_maxwindowlocation,
            gene_nextlargestwindowcount, gene_nextwindowwidth, gene_nextwindowlocation, gene_call, gene_ighsource,
            gene_bundle_count)


#################################################################
# Make Column Names for Translocation Table
list_of_genes = ['NSD2', 'CCND3', 'MYC', 'MAFA', 'CCND1', 'CCND2', 'MAF', 'MAFB']
list_of_features = ['COUNT', 'WINDOW', 'MAXWINDOWCOUNT', 'MAXWINDOWWIDTH', 'MAXWINDOWLOCATION',
                    'NEXTWINDOWCOUNT', 'NEXTWINDOWWIDTH', 'NEXTWINDOWLOCATION', 'CALL', 'IGSOURCE', 'BUNDLECOUNT']

column_names = ['SAMPLE']
# Make list of column_names:
for gene in list_of_genes:
    for feature in list_of_features:
        header = '_'.join([gene, feature])
        column_names.append(header)


#################################################################
# Function to call all functions for a sample

def make_sample_calls(sample):
    sample_id = sample

    # MMSET
    mmset_call = call_translocations(gene="NSD2", gene_chr="chr4", order="standard", sample=sample_id)

    # CCND3
    ccnd3_call = call_translocations(gene="CCND3", gene_chr="chr6", order="standard", sample=sample_id)

    # MYC
    myc_call = call_translocations(gene="MYC", gene_chr="chr8", order="standard", window_start=124987758,
                                   window_end=129487754, sample=sample_id)

    # MAFA
    mafa_call = call_translocations(gene="MAFA", gene_chr="chr8", order="standard", window_start=142918584,
                                    window_end=143925832, sample=sample_id)

    # CCND1
    ccnd1_call = call_translocations(gene="CCND1", gene_chr="chr11", order="standard", sample=sample_id)

    # CCND2
    ccnd2_call = call_translocations(gene="CCND2", gene_chr="chr12", order="standard", sample=sample_id)

    # MAF
    maf_call = call_translocations(gene="MAF", gene_chr="chr16", order="reverse", sample=sample_id)

    # MAFB
    mafb_call = call_translocations(gene="MAFB", gene_chr="chr20", order="reverse", sample=sample_id)

    results = ((sample,), mmset_call, ccnd3_call, myc_call, mafa_call, ccnd1_call, ccnd2_call, maf_call, mafb_call)
    return (sum(results, ()))


#################################################################
# Create pool for multithreading

pool = Pool(processes=2)

jobs = []
job = pool.apply_async(make_sample_calls, [sample])
jobs.append(job)

results = []
for job in jobs:
    print(job)
    result = job.get()
    print(result)
    results.append(result)


translocationsTable = pd.DataFrame(results, columns=column_names)

# Save results
translocationsTable.to_csv(args.output_file, sep="\t", index=False, na_rep="NaN")

#################################################################

# Message all done
print('ALL PROCESSES COMPLETE')
