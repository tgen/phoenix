# Usage
# python Process_Discordant_SAM_File.py <Collapsed_TophatFusion_File> <Discordant_Reads.sam>

# requires python 3.4 with pandas, numpy, pysam packages

print("Importing Packages")
# Configure Enviroment
import pandas as pd
import numpy as np
import pysam
import sys
import os.path
import re
from subprocess import call

# Variables
# call_threshold = 3
# minimum_window_count = 2
window_size = int(sys.argv[3])
min_reads = 5

# tophat = int(sys.argv[4])
minimum_count_threshold = 0.75

print("Defining Functions")


# Function to process SAM file from SAMBLASTER DISCORDANT EXPORT
def bam_to_df(bam, chr=None, start=None, stop=None, file_name=None):
  file_name = os.path.basename(file_name)
  file_name = file_name.replace('\_R1\_Trinity\_sorted\.bam', '')
  print('Updated file name is ' + file_name)
  index = 0
  seq = []
  name = []
  r1_chr = []
  r1_pos = []
  r1_is_read1 = []
  r1_is_reversed = []
  r1_cigar = []
  r1_mapq = []
  r2_is_reversed = []
  r2_chr = []
  r2_pos = []
  frag_length = []
  r2_cigar = []
  r2_mapq = []
  pair_orientation = []
  for read in bam.fetch(chr, start, stop):
    # index = index + 1
    seq.append(read.query_sequence)
    # name.append(file_name+":"+read.query_name) #+"_"+str(index))
    name.append(read.query_name)
    r1_chr.append(read.reference_name)
    r1_pos.append(read.reference_start)
    r1_is_read1.append(read.is_read1)
    r1_is_reversed.append(read.is_reverse)
    r1_cigar.append(read.cigarstring)
    r1_mapq.append(read.mapping_quality)
    r2_is_reversed.append(read.mate_is_reverse)
    r2_chr.append(read.next_reference_name)
    r2_pos.append(read.next_reference_start)
    frag_length.append(read.template_length)
    if read.has_tag('MC') == True:
      r2_cigar.append(read.get_tag('MC'))
    else:
      r2_cigar.append('NA')
    if read.has_tag('MQ') == True:
      r2_mapq.append(read.get_tag('MQ'))
    else:
      r2_mapq.append(255)
    if read.is_reverse is False and read.mate_is_reverse is True:
      pair_orientation.append("FR")
    elif read.is_reverse is True and read.mate_is_reverse is False:
      pair_orientation.append("RF")
    elif read.is_reverse is False and read.mate_is_reverse is False:
      pair_orientation.append("FF")
    elif read.is_reverse is True and read.mate_is_reverse is True:
      pair_orientation.append("RR")
    else:
      pair_orientation.append("Error")
    # print(str(r1_cigar) + " and  " + str(name) + " mapq "+str(r1_mapq))
  return pd.DataFrame({'seq': seq,
                       'name': name,
                       'r1_pos': r1_pos,
                       'r1_chr': r1_chr,
                       'r1_is_read1': r1_is_read1,
                       'r1_is_reversed': r1_is_reversed,
                       'r1_cigar': r1_cigar,
                       'r1_mapq': r1_mapq,
                       'r2_is_reversed': r2_is_reversed,
                       'r2_chr': r2_chr,
                       'r2_pos': r2_pos,
                       'frag_length': frag_length,
                       'r2_cigar': r2_cigar,
                       'r2_mapq': r2_mapq,
                       'pair_orientation': pair_orientation})


# function to get junction break
def get_junctionBreak(cigarStr):
  # split cigar by M,S,H
  Mval = 0
  if (cigarStr != ''):
    tmp_cigar = cigarStr
    tmp_cigar = tmp_cigar.replace('H', 'M')
    tmp_cigar = tmp_cigar.replace('S', 'M')
    tmp_cigar = tmp_cigar.replace('I', 'M')
    tmp_cigar = tmp_cigar.replace('D', 'M')
    # print( tmp_cigar)
    cigar_list = tmp_cigar.split('M')
    for values in cigar_list:
      if (values != '' and cigarStr.find(str(values) + 'M') != -1):
        if (int(values) > Mval):
          Mval = int(values)
  return Mval


# function to get junction break
def get_contigLength(cigarStr):
  con_len = 0
  if (cigarStr != ''):
    tmp_cigar = cigarStr
    tmp_cigar = tmp_cigar.replace('H', 'M')
    tmp_cigar = tmp_cigar.replace('S', 'M')
    tmp_cigar = tmp_cigar.replace('I', 'M')
    tmp_cigar = tmp_cigar.replace('D', 'M')
    cigar_list = tmp_cigar.split('M')
    for values in cigar_list:
      if (values != ''):
        con_len = con_len + int(values)
  return con_len


def get_longest_seqforcontig(contig_read_table_tmp):
  longest_seq = ""
  for row in contig_read_table_tmp.index:
    next_seq = contig_read_table_tmp.at[row, 'seq']
    if (len(next_seq) > len(longest_seq)):
      longest_seq = next_seq
  return longest_seq


########################################
#
# Only doing a complement not reverse
#
########################################
def reverseComplement(seq):
  # Reverse sequence string
  rseq = seq  # [::-1]

  rseq = rseq.replace('A', 'B')
  rseq = rseq.replace('T', 'A')
  rseq = rseq.replace('B', 'T')

  rseq = rseq.replace('C', 'B')
  rseq = rseq.replace('G', 'C')
  rseq = rseq.replace('B', 'G')

  return rseq


# Function to check contigs
def check_contigs(contig_table, fastq_path, reads_sam_path, igregions, window_size=200, min_reads=5, contig_perc=0.1):
  list_names = []

  # find unique contigs
  for row in contig_table.index:
    r1_cigar = contig_table.at[row, 'r1_cigar']
    name = contig_table.at[row, 'name']
    if (name not in list_names):
      list_names.append(name)

  # loop through to unique contigs to generate dataframe
  names2 = []
  for names in list_names:
    # extract by name
    print(names)

    # add data to
    names2.append(names)
  final_table = pd.DataFrame({'name': names2
                              })

  # Process each contig
  for currIndex in final_table.index:
    names = final_table.at[currIndex, 'name']
    print("\nNow Processing" + names)
    contig_table_by_read = contig_table[(contig_table.name == names)]
    # if contigs aligns to multiple locations
    count = len(contig_table_by_read.index)
    if (count > 1):  # only need Txs
      # initialize vars for loop
      loop_var = 1
      # track mininum mval
      min_mVal = 100000  # set to a large val
      # track variables for IG and Gene
      IG_found = ""
      Gene_found = ""

      for contig_row in contig_table_by_read.index:
        contig_chr = contig_table_by_read.at[contig_row, 'r1_chr']
        contig_pos = contig_table_by_read.at[contig_row, 'r1_pos']
        # check if gene is in Tx region
        # @Sara Need to read updated bedfile as agreed by JJK
        mygene = isKnownTx(contig_chr, int(contig_pos), igregions)

        if (mygene != "" and contig_chr in ['chr2', 'chr14', 'chr22'] and IG_found == ""):
          IG_found = mygene
        if (mygene != "" and contig_chr in ['chr4', 'chr6', 'chr8', 'chr11', 'chr12', 'chr16',
                                            'chr20'] and Gene_found == ""):
          Gene_found = mygene

        # print("contig_chr="+contig_chr +"contig_pos " + str(contig_pos) + mygene)

        mVal = get_junctionBreak(contig_table_by_read.at[contig_row, 'r1_cigar'])
        if (mVal <= min_mVal):
          min_mVal = mVal
          print("Value lower than thr " + str(mVal))
        juncbreak = mVal

        contig_len = get_contigLength(contig_table_by_read.at[contig_row, 'r1_cigar'])

        # since bwa hardclips supplementary reads find the longest seq for contig
        contig = get_longest_seqforcontig(contig_table_by_read)  # contig_table_by_read.at[contig_row,'seq']

        # extract fastq name from contig name
        namesSplit = names.split(':')

        fastq = namesSplit[0]
        fastq = fastq + '_R1_001.fastq.gz'

        ffastq = fastq_path + "/" + fastq
        reads_sam = reads_sam_path + "/" + namesSplit[0] + "/" + namesSplit[0] + "_ReadstoContigs.bam"
        # print("\nfastq name is " + fastq)

        # extract 25mers left and right of contig
        jright = juncbreak + 25
        jleft = juncbreak - 25

        # extract chromosome coordinates 25mers around contig
        location = contig_table_by_read.at[contig_row, 'name'] + ":" + str(jleft) + "-" + str(jright)

        rr = contig[jleft:jright:1]
        rr_rev = reverseComplement(rr)
        # get junction count  read 1 and read2
        r1_count = getReadsatjunction(rr, contig, ffastq) + getReadsatjunction(rr_rev, contig, ffastq)
        ffastq_2 = ffastq
        ffastq_2 = ffastq_2.replace("_R1", "_R2")

        r2_count = getReadsatjunction(rr, contig, ffastq_2) + getReadsatjunction(rr_rev, contig, ffastq_2)

        # get unique fragments at junction
        frag_count = getFragsatJunction_samtools(location, contig, reads_sam)
        print("read count is " + str(r1_count))

        # define df columns
        CIGAR = "cigar_" + str(loop_var)
        GENE = "Gene_" + str(loop_var)
        LEN = "length_" + str(loop_var)
        PERC = "percent_of_contig_at_Gene_" + str(loop_var)
        R1_C = "R1_reads_at_junc_" + str(loop_var)
        R2_C = "R2_reads_at_junc_" + str(loop_var)
        FRAG_C = "Fragments_at_junc_" + str(loop_var)
        POS_START = "pos_" + str(loop_var) + "_start"
        POS_END = "pos_" + str(loop_var) + "_end"
        REV = "r1_is_reversed" + str(loop_var)

        if (mVal >= window_size):

          final_table.at[[currIndex], CIGAR] = contig_table_by_read.at[contig_row, 'r1_cigar']
          final_table.at[[currIndex], GENE] = mygene
          final_table.at[[currIndex], LEN] = mVal
          final_table.at[[currIndex], PERC] = mVal * 100 / contig_len
          final_table.at[[currIndex], R1_C] = r1_count
          final_table.at[[currIndex], R2_C] = r2_count
          final_table.at[[currIndex], FRAG_C] = r2_count
          final_table.at[[currIndex], POS_START] = contig_table_by_read.at[contig_row, 'r1_pos']
          final_table.at[[currIndex], POS_END] = contig_table_by_read.at[contig_row, 'r1_pos'] + mVal
          final_table.at[[currIndex], REV] = contig_table_by_read.at[contig_row, 'r1_is_reversed']

          # increment counter

          loop_var = loop_var + 1
          print("\nIG " + IG_found + "gene " + Gene_found)
          if (count > 1 and IG_found != "" and Gene_found != ""):
            final_table.at[[currIndex], 'IgTxCalled'] = 1
  return final_table


###############################################
# Function to check if given  location (chr:pos) is in list of IgTx regions
###############################################
def isKnownTx(qchr, qpos, listofRegions):
  txGene = ''
  chr_region = listofRegions[(listofRegions.chr == qchr)]
  for row in chr_region.index:
    if ((chr_region.at[row, 'start'] <= qpos) and (chr_region.at[row, 'stop'] >= qpos)):
      if (qchr == 'chr14'):
        txGene = "IGH"
      if (qchr == 'chr2'):
        txGene = "IGK"
      if (qchr == 'chr22'):
        txGene = "IGL"
      if (qchr == 'chr4'):
        txGene = "NSD2"
      if (qchr == 'chr11'):
        txGene = "CCND1"
      if (qchr == 'chr12'):
        txGene = "CCND2"
      if (qchr == 'chr6'):
        txGene = "CCND3"
      if (qchr == 'chr16'):
        txGene = "MAF"
      if (qchr == 'chr8' and qpos >= 142918584 and qpos <= 143925832):
        txGene = "MAFA"
      if (qchr == 'chr20'):
        txGene = "MAFB"
      if (qchr == 'chr8' and qpos >= 124987758 and qpos <= 129487754):
        txGene = "MYC"
  return txGene


###############################################
#  Simple approach to find reads at junction
#  Need to be replaced with samtools in next
#  update
###############################################
def getReadsatjunction(region, contig, fastq):
  r1_count = 0
  # grep region1 from fastq
  print(fastq)
  print("\n Now checking for this sequence in fastq \n")
  print(region)
  # @Bryce Should be changed to a tmp location
  tmp_file = sys.argv[5] + "/r1counts.txt"
  status = call("zcat " + fastq + " | grep " + region + " | wc -l > " + tmp_file, shell=True)
  if status < 0:
    print("### Cat Command Failed....now exiting!!")
    sys.exit(-1)
  else:
    with open(tmp_file, "r") as myfile:
      data = myfile.readlines()
      r1_count = int(data[0])
  return r1_count


###############################################
#  Simple approach to find reads at junction
#  Need to be replaced with samtools in next
#   update. NOT CURRENTLY USED
###############################################
def getAllMappedReadsatJunction(region, contig, readsam):
  r1_count = 0
  r2_count = 0
  # grep region1 from fastq
  # print(fastq)
  print("\n Now checking for this sequence in fastq \n" + contig)
  print(region)

  region_rev = reverseComplement(region)

  # @Bryce Should be changed to a tmp location
  tmp_file = sys.argv[5] + "/r1counts.txt"
  # status = call("/home/snasser/toolkit_snasser/IgTxAssembly/findinfile.sh "+ fastq +" "+region+" "+tmp_file, shell=True)
  print("egrep \"" + region + "|" + region_rev + "\" " + readsam + " | awk '{ print $1 }' | sort | uniq | wc -l ")
  status = call(
    "egrep \"" + region + "|" + region_rev + "\" " + readsam + " | awk '{ print $1 }' | sort | uniq | wc -l  > " + tmp_file,
    shell=True)
  if status < 0:
    print("### Cat Command Failed....now exiting!!")
    sys.exit(-1)
  else:
    with open(tmp_file, "r") as myfile:
      data = myfile.readlines()
      r1_count = int(data[0])
  return r1_count


###############################################
#  Simple approach to find reads at junction
#  Need to be replaced with samtools in next
#   update
###############################################
def getFragsatJunction_samtools(location, contig, bam):
  r1_count = 0
  r2_count = 0
  # grep region1 from fastq
  # print(fastq)
  print("\n Now checking for this sequence in fastq \n")
  # print(region)
  # @Bryce Should be changed to a tmp location
  tmp_file = sys.argv[5] + "/r1counts.txt"
  # "zcat " + fastq + " | grep " + region + " > " + tmp_file, shell=True

  # region_rev = reverseComplement(region)

  status = call(
    "samtools view " + bam + " | grep " + location + " | awk '{ print $1 }' | sort | uniq | wc -l > " + tmp_file,
    shell=True)
  if status < 0:
    print("### Cat Command Failed....now exiting!!")
    sys.exit(-1)
  else:
    with open(tmp_file, "r") as myfile:
      data = myfile.readlines()
      r1_count = int(data[0])
  return r1_count


# END FUNCTION DEFINITIONS

print("Importing Data")

print("Reading SAM File")

# Create python object for the SAM file
assm_sam_file = sys.argv[1]

# get fastq
fastq_str = sys.argv[2]
fastq_path = fastq_str
reads_sam_file_path = sys.argv[5]
fastq_list = []  # fastq_str.split(',')

# get Ig/parter bed file
ig_headers = ['chr', 'start', 'stop', 'name']
igregions = pd.read_csv(sys.argv[4], sep="\t", header=None, names=ig_headers)

# get outpath
out_path = sys.argv[5]

results_full_table = pd.DataFrame()
index_sam = 0
table_hrd = []
# for assm_sam_file in bam_list:
if (assm_sam_file != ""):
  samfile = pysam.AlignmentFile(assm_sam_file, "r")

  print("Creating SAM File Table")
  # Call Function to convert SAM file into pandas dataframe

  table = bam_to_df(samfile, file_name=assm_sam_file)
  table_hrd = list(table.columns)
  results_table = check_contigs(table, fastq_path, reads_sam_file_path, igregions, window_size, min_reads)

# filter results table
filt_table = results_table[
  results_table.Gene_1.isin(['IGH', 'IGL', 'NSD2', 'CCND1', 'CCND2', 'CCND3', 'IGK', 'MYC', 'MAFA', 'MAFB', 'MAF'])]
print("\nWriting outfile ****************")

out_file = out_path + "/ContigResults.txt"
out_file_filt = out_path + "/FilteredContigResults.txt"
results_table.to_csv(out_file, sep="\t", index=False, na_rep=0, float_format='%.0f')
filt_table.to_csv(out_file_filt, sep="\t", index=False, na_rep=0, float_format='%.0f')
print("Test Done")
