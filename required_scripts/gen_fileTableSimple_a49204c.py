#!/usr/bin/env python3

# MYC blacklist region chr8:129571461-129571635 (b37)
# MYC blacklist region chr8:128559215-128559389 (b38)

#################################################################
# Configure the enviroment
import pandas as pd
import argparse
import numpy as np
import re
from multiprocessing import Pool

#################################################################

parser = argparse.ArgumentParser(
  description='Process pairoscope discordant reads to make myeloma Ig translocation calls.')

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

table_headers = ['Specimen', 'ChrA', 'PositionA', 'alt', 'svtype', 'PRAll',
                 'SRAll']  # , 'RO', 'RUT1','RUT2','RDISTDISC1','RDISTDISC2','RCDIS1','RCDIS2','DRNOISE1','DRNOISE1','EVALRCNOISE','CLMQ1','CLMQ2'];
IgChrList = ['chr2', 'chr22', 'chr14']
dreads = pd.read_csv(args.input_file, sep="\t", header=0, names=table_headers)
# set chrB and posnB
dreads['ChrB'] = ''
dreads['PositionB'] = 0
dreads['PR_T'] = 0
dreads['SR_T'] = 0
add_cols = ['ChrB', 'PositionB', 'PR_T', 'SR_T'];
tot_cols = table_headers
tot_cols.append(add_cols)
tot_cols[1] = 'Speciment'

for row in dreads.index:
  vals = []
  alt = dreads.at[row, 'alt']
  if (alt.find(']') != -1):
    print(alt)
    vals = re.split('\:|\]', alt)
    # print("x0" + vals[0])
    print("x1" + vals[1])
    print("x2" + vals[2])
    dreads.at[row, 'ChrB'] = vals[1]
    dreads.at[row, 'PositionB'] = vals[2]
  elif (alt.find('[') != -1):
    print("y" + alt)
    # print("y1" + vals[1])
    vals = re.split('\:|\[', alt)
    print("y0" + vals[1])
    print("y1" + vals[2])
    dreads.at[row, 'ChrB'] = vals[1]
    dreads.at[row, 'PositionB'] = vals[2]

  if (dreads.at[row, "ChrA"] in IgChrList):
    print("swap")
    # swap chrA. chrB
    temp_chr = dreads.at[row, 'ChrA']
    temp_pos = dreads.at[row, 'PositionA']

    dreads.at[row, "ChrA"] = dreads.at[row, 'ChrB']
    dreads.at[row, "ChrB"] = temp_chr

    dreads.at[row, "PositionA"] = dreads.at[row, 'PositionB']
    dreads.at[row, 'PositionB'] = temp_pos

  print(dreads.at[row, 'PRAll'])
  print(dreads.at[row, 'SRAll'])
  PR = dreads.at[row, 'PRAll']
  SR = dreads.at[row, 'SRAll']
  if (PR.find(',') != -1):
    print("PR is " + str(PR))
    PRvals = re.split(',', PR)
    dreads.at[row, 'PR_T'] = PRvals[1]

  if (SR.find(',') != -1):
    print("SR is " + str(SR))
    SRvals = re.split(',', SR)
    dreads.at[row, 'SR_T'] = SRvals[1]


# dreads=format_basic(dreads)
#################################################################
# Define Functions
# def format_basic()


def call_translocations(sample, gene, gene_chr, order, window_start=0, window_end=300000000):
  # initialize variables
  # "CCND3 Count\tCCND3 PR\tCCND3 SR\tCCND3 Posn\tCCND3 Partner Posn\tCCND3 Call\t";
  # gene_window = 0
  # gene_maxwindowcount = 0
  # gene_nextlargestwindowcount = 0
  # gene_windowend = 0
  # gene_maxwindowwidth = 0
  # gene_maxwindowlocation = 0
  # gene_nextwindowwidth = 0
  # gene_nextwindowlocation = 0
  gene_igChr = ""
  gene_igPosn = 0
  gene_PR = 0
  gene_SR = 0
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
    print("IG counts: ")
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
      gene_igChr = "chr14"
      # gene_igPosn = gene_igh.at['PositionB']

    elif count_gene_igk > count_gene_igh and count_gene_igk > count_gene_igl:
      print('2-In IgK Loop')
      # Suspect igK Translocation
      # WARNiNG THE CHROMOSOME ORDER BECOMES UNEXPECTED, HOW TO FiX?
      # Create Flipped Table for the kappa counts (	ChrA	PositionA	ChrB	PositionB)
      # Update the column headers
      print(gene_igh)
      print("************")
      # 'Specimen', 'ChrA', 'PositionA', 'alt', 'svtype','PRAll','SRAll
      new_cols = ['Speciment', 'ChrA', 'PositionA', 'alt', 'svtype', 'PRAll', 'SRAll', 'ChrB', 'PositionB', 'PR_T',
                  'SR_T']
      gene_igk.columns = new_cols
      # Create new GENE table with the relabelled columns in the expected order
      cols = ['Speciment', 'ChrA', 'PositionA', 'alt', 'svtype', 'PRAll', 'SRAll', 'ChrB', 'PositionB', 'PR_T', 'SR_T']
      table_gene = gene_igk[cols]
      count_gene = len(table_gene.index)
      gene_ighsource = 2
      gene_igChr = "chr2"
      # gene_igPosn = gene_igk.loc[0,8]
    elif count_gene_igl > count_gene_igh and count_gene_igl > count_gene_igk:
      print('2-In IgL Loop')
      # Suspect igL Translocation
      table_gene = gene_igl
      count_gene = len(table_gene.index)
      gene_ighsource = 3
      gene_igChr = "chr22"
      # gene_igPosn = gene_igl.loc[0,8]
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
        gene_igChr = "chr14"
        # gene_igPosn = "XX"
      elif count_gene_igh < count_gene_igk and count_gene_igh < count_gene_igl and count_gene_igk == count_gene_igl:
        print('DEFAULTiNG TO iGH COUNTS BECAUSE IgK and IgL ARE EQUAL EVEN THOUGH LESS')
        table_gene = gene_igh
        count_gene = len(table_gene.index)
        gene_ighsource = 1
        gene_igChr = "chr14"
        # gene_igPosn = "XX"

  elif order == "reverse":
    print('1 - In Reverse Loop')
    gene_igh = dreads[(dreads.Specimen == sample) & (dreads.ChrB == "chr14") & (dreads.ChrA == gene_chr)
                      & (dreads.PositionA >= window_start) & (dreads.PositionA <= window_end)]
    count_gene_igh = len(gene_igh.index)
    gene_igk = dreads[(dreads.Specimen == sample) & (dreads.ChrB == "chr2") & (dreads.ChrA == gene_chr)
                      & (dreads.PositionA >= window_start) & (dreads.PositionA <= window_end)]
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
      new_cols = ['Speciment', 'ChrA', 'PositionA', 'alt', 'svtype', 'PRAll', 'SRAll', 'ChrB', 'PositionB', 'PR_T',
                  'SR_T']
      #            [Specimen,     ChrA, PositionA, alt, svtype, PRAll, SRAll, chrB, PositionB, PR_T, SR_T, ChrB]
      print(gene_igh)
      gene_igh.columns = new_cols
      # Create new GENE table with the relabelled columns in the expected order
      cols = ['Speciment', 'ChrA', 'PositionA', 'alt', 'svtype', 'PRAll', 'SRAll', 'ChrB', 'PositionB', 'PR_T', 'SR_T']
      table_gene = gene_igh[cols]
      count_gene = len(table_gene.index)
      gene_ighsource = 1
      gene_igChr = "chr14"
      # gene_igPosn = gene_igh.loc[0,8]
    elif count_gene_igk > count_gene_igh and count_gene_igk > count_gene_igl:
      print('2-In IgK Loop')
      # Suspect igK Translocation
      # WARNiNG THE CHROMOSOME ORDER BECOMES UNEXPECTED, HOW TO FiX?
      # Create Flipped Table for the kappa counts (	ChrA	PositionA	ChrB	PositionB)
      # Update the column headers
      new_cols = ['Speciment', 'ChrA', 'PositionA', 'alt', 'svtype', 'PRAll', 'SRAll', 'ChrB', 'PositionB', 'PR_T',
                  'SR_T']
      gene_igk.columns = new_cols
      # Create new GENE table with the relabelled columns in the expected order
      cols = ['Speciment', 'ChrA', 'PositionA', 'alt', 'svtype', 'PRAll', 'SRAll', 'ChrB', 'PositionB', 'PR_T', 'SR_T']
      table_gene = gene_igk[cols]
      count_gene = len(table_gene.index)
      gene_ighsource = 2
      gene_igChr = "chr2"
      # gene_igPosn = gene_igk.loc[0,8]
    elif count_gene_igl > count_gene_igh and count_gene_igl > count_gene_igk:
      print('2-In IgL Loop')
      # Suspect igL Translocation
      table_gene = gene_igl
      count_gene = len(table_gene.index)
      gene_ighsource = 3
      gene_igChr = "chr22"
      # gene_igPosn = gene_igk.loc[0,8]
    else:
      # Error capture NOT AN EXPECTED EVENT
      # Could be 0 in all three potentially, then what?
      # What if two but not all three have the same number of counts?
      print('2 - Argh - WHAT CAUSES THiS TO HAPPEN')
      if count_gene_igh < 3 and count_gene_igk < 3 and count_gene_igl < 3:
        print('DEFAULTiNG TO iGH COUNTS BECAUSE ALL LESS THAN 3')
        # Update the column headers
        new_cols = ['Speciment', 'ChrA', 'PositionA', 'alt', 'svtype', 'PRAll', 'SRAll', 'ChrB', 'PositionB', 'PR_T',
                    'SR_T']
        gene_igh.columns = new_cols
        # Create new GENE table with the relabelled columns in the expected order
        cols = ['Speciment', 'ChrA', 'PositionA', 'alt', 'svtype', 'PRAll', 'SRAll', 'ChrB', 'PositionB', 'PR_T',
                'SR_T']
        table_gene = gene_igh[cols]
        count_gene = len(table_gene.index)
        gene_ighsource = 1
        gene_igChr = "chr14"
        gene_igPosn = "XX"  # gene_igh.loc[0,8]
      elif count_gene_igh < count_gene_igk and count_gene_igh < count_gene_igl and count_gene_igk == count_gene_igl:
        print('DEFAULTiNG TO iGH COUNTS BECAUSE IgK and IgL ARE EQUAL EVEN THOUGH LESS')
        # Update the column headers
        new_cols = ['Speciment', 'ChrA', 'PositionA', 'alt', 'svtype', 'PRAll', 'SRAll', 'ChrB', 'PositionB', 'PR_T',
                    'SR_T']
        gene_igh.columns = new_cols
        # Create new GENE table with the relabelled columns in the expected order
        cols = ['Speciment', 'ChrA', 'PositionA', 'alt', 'svtype', 'PRAll', 'SRAll', 'ChrB', 'PositionB', 'PR_T',
                'SR_T']
        table_gene = gene_igh[cols]
        count_gene = len(table_gene.index)
        gene_ighsource = 1
        gene_igChr = "chr14"
        gene_igPosn = "XX"  # gene_igk.loc[0,8]
  # gene_PR = 0
  # gene_SR = 0
  print("Table  Gene is ....")
  print(table_gene)
  if (count_gene > 0):
    gene_call = 1
    # Using hard coded index values is not good, update this later
    # gene_igChr = table_gene.iat[0,7]
    gene_igPos = table_gene.iat[0, 8]
    gene_BP = table_gene.iat[0, 2]
    gene_PR = table_gene.iat[0, 9]
    gene_SR = table_gene.iat[0, 10]
  else:
    gene_call = 0
    gene_igChr = "0"
    gene_igPos = 0
    gene_BP = 0
    gene_PR = 0
    gene_SR = 0
    gene_ighsource = 0

  return (gene_BP, gene_igPos, gene_igChr, gene_PR, gene_SR, gene_call, gene_ighsource)


#  return (count_gene, gene_window, gene_maxwindowcount, gene_maxwindowwidth, gene_maxwindowlocation,
#          gene_nextlargestwindowcount, gene_nextwindowwidth, gene_nextwindowlocation, gene_call, gene_ighsource,
#          gene_bundle_count)


#################################################################
# Make Column Names for Translocation Table
list_of_genes = ['NSD2', 'CCND3', 'MYC', 'MAFA', 'CCND1', 'CCND2', 'MAF', 'MAFB']
list_of_features = ['GENE_BP', 'IG_BP', 'IG_CHR', 'PR_T', 'SR_T', 'CALL', 'IGSOURCE']
# list_of_features = ['COUNT', 'WINDOW', 'MAXWINDOWCOUNT', 'MAXWINDOWWIDTH', 'MAXWINDOWLOCATION',
#                    'NEXTWINDOWCOUNT', 'NEXTWINDOWWIDTH', 'NEXTWINDOWLOCATION', 'CALL', 'IGSOURCE', 'BUNDLECOUNT']

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
  mmset_call = call_translocations(gene="NSD2", gene_chr="chr4", order="standard", window_start=1798273,
                                   window_end=1998273, sample=sample_id)
  # CCND3
  ccnd3_call = call_translocations(gene="CCND3", gene_chr="chr6", order="standard", window_start=41632262,
                                   window_end=42332262, sample=sample_id)

  # MYC
  myc_call = call_translocations(gene="MYC", gene_chr="chr8", order="standard", window_start=124987758,
                                 window_end=129487754, sample=sample_id)

  # MAFA
  mafa_call = call_translocations(gene="MAFA", gene_chr="chr8", order="standard", window_start=142918584,
                                  window_end=143925832, sample=sample_id)

  # CCND1
  ccnd1_call = call_translocations(gene="CCND1", gene_chr="chr11", order="standard", window_start=68732532,
                                   window_end=69685232, sample=sample_id)

  # CCND2
  ccnd2_call = call_translocations(gene="CCND2", gene_chr="chr12", order="standard", window_start=3690834,
                                   window_end=4690834, sample=sample_id)

  # MAF
  maf_call = call_translocations(gene="MAF", gene_chr="chr16", order="reverse", window_start=78096103,
                                 window_end=79866103, sample=sample_id)

  # MAFB
  mafb_call = call_translocations(gene="MAFB", gene_chr="chr20", order="reverse", window_start=39671358,
                                  window_end=40971360, sample=sample_id)

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
