# Usage
#
# (c) TGEN 2020
# requires python 3.0 with pandas, numpy, pysam packages
#
#

print("Importing Packages")
# Configure Enviroment
import pandas as pd
import numpy as np
import pysam
import sys
import os.path
import re
import argparse
from subprocess import call
import itertools

# Variables

minimum_count_threshold = 0.75

parser = argparse.ArgumentParser(description='Assemble DLE fastqs and call myeloma Ig translocation calls.')

parser.add_argument('-i', '--input_bam',
                    required=True,
                    metavar='File.bam',
                    dest="input_file",
                    help='Merged bam of all regional contigs from DLE')
parser.add_argument('-s', '--specimen',
                    required=True,
                    help='Specimen Name, must match the bam')
parser.add_argument('-o', '--output_path',
                    required=True,
                    metavar='File.tsv',
                    help='Output path')
parser.add_argument('-c', '--contig',
                    default=200,
                    type=int,
                    metavar='INT',
                    help='Minimum contig overlap length at each breakend ')
parser.add_argument('-r', '--reads',
                    default=5,
                    type=int,
                    metavar='INT',
                    help='Minimum fragments at junction')
parser.add_argument('-w', '--window',
                    default=25,
                    type=int,
                    metavar='INT',
                    help='Window length in basepairs to check for fragments at junction ')
parser.add_argument('-f', '--fastq_dir',
                    default='Path',
                    required=True,
                    metavar='/scratch/',
                    help='Path to original fastqs')
parser.add_argument('-b', '--bam_path',
                    default='Path',
                    required=True,
                    metavar='/scratch/',
                    help='Path to bams for each region with fastq mapped to bam')
parser.add_argument('-p', '--pairoscope_bed_file',
                    default='BED',
                    required=True,
                    metavar='regions.bed',
                    help='BED file with MM Ig Tx regions ')

args = parser.parse_args()

###########################################
# args to variables

window_size = args.contig
bp_win_at_junc = args.window
min_reads = args.reads
sample_name = args.specimen
out_path = args.output_path
assm_sam_file = args.input_file
bed_file = args.pairoscope_bed_file
reads_sam_file_path = args.bam_path
fastq_path = args.fastq_dir

# define constants
ig_dict = {
  "IGH": "1",
  "IGK": "2",
  "IGL": "3"
}

print("Defining Functions")


# Function to process SAM file from SAMBLASTER DISCORDANT EXPORT
def bam_to_df(bam, chr=None, start=None, stop=None, file_name=None):
  file_name = os.path.basename(file_name)
  file_name = file_name.replace('\_R1\_Trinity\_sorted\.bam', '')
  # print('Updated file name is ' + file_name)
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
    seq.append(read.query_sequence)
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


###############################
# function to get junction break
###############################
def get_junctionBreak(cigarStr):
  # split cigar by M,S,H
  Mval = 0
  if (cigarStr != ''):
    tmp_cigar = cigarStr
    tmp_cigar = tmp_cigar.replace('H', 'M')
    tmp_cigar = tmp_cigar.replace('S', 'M')
    tmp_cigar = tmp_cigar.replace('I', 'M')
    tmp_cigar = tmp_cigar.replace('D', 'M')

    total_len = 0
    cigar_list = tmp_cigar.split('M')
    for values in cigar_list:
      if (values != '' and cigarStr.find(str(values) + 'M') != -1):
        if (int(values) > Mval):
          Mval = int(values)
  return Mval


##################################
# function to get junction break
##################################
def get_junctionBreakNew(cigarStr):
  Mval = 0
  Sval = 0
  Hval = 0
  Junc = 0
  if (cigarStr != ''):
    tmp_cigar = cigarStr
    tmp_cigar = tmp_cigar.replace('H', 'M')
    tmp_cigar = tmp_cigar.replace('S', 'M')
    tmp_cigar = tmp_cigar.replace('I', 'M')
    tmp_cigar = tmp_cigar.replace('D', 'M')

    total_len = 0
    cigar_list = tmp_cigar.split('M')
    for values in cigar_list:
      if (values != '' and cigarStr.find(str(values) + 'M') != -1):
        if (int(values) > Mval):
          Mval = int(values)
      # Generall only one S and/or H, keep first S or H
      if (values != '' and cigarStr.find(str(values) + 'S') != -1 and Sval == 0):
        Sval = int(values)
      if (values != '' and cigarStr.find(str(values) + 'H') != -1 and Hval == 0):
        Hval = int(values)

    # if clipping is sbefore M  then Mval =Mval -length
    sindex = cigarStr.find('S')
    hindex = cigarStr.find('H')
    mindex = cigarStr.find('M')

    if (sindex > 0 and hindex > 0):
      if (sindex < mindex and sindex < hindex):
        Junc = Sval
      elif (hindex < mindex and hindex < sindex):
        Junc = Hval
    elif (sindex > 0 and sindex < mindex):  # mval after s
      Junc = Sval
    elif (hindex > 0 and hindex < mindex):
      Junc = Hval
    elif (mindex < sindex):
      Junc = Mval
    elif (mindex < hindex):
      Junc = Mval
  return Junc


##########################################
#
##########################################
def get_junctionBreakNewMPos(cigarStr):
  Mval = 0
  Sval = 0
  Hval = 0
  Junc = 0
  if (cigarStr != ''):
    tmp_cigar = cigarStr
    tmp_cigar = tmp_cigar.replace('H', 'M')
    tmp_cigar = tmp_cigar.replace('S', 'M')
    tmp_cigar = tmp_cigar.replace('I', 'M')
    tmp_cigar = tmp_cigar.replace('D', 'M')

    total_len = 0
    final_len = 0
    cigar_list = tmp_cigar.split('M')
    for values in cigar_list:
      if (values != ''):
        if (cigarStr.find(str(values) + 'M') != -1):
          if (int(values) > Mval):
            Mval = int(values)
        final_len = final_len + int(values)
    # get posn
    for values in cigar_list:
      if (values != ''):
        total_len = total_len + int(values)
        if ((cigarStr.find(str(values) + 'M') != -1) and (int(values) == int(Mval))):
          if (total_len == final_len):
            total_len = total_len - Mval  # if M was the last matched go to prev junction
          break

    Junc = total_len
  return Junc


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


########################################
#
########################################
def get_longest_seqforcontig(contig_read_table_tmp, cur_row):
  longest_seq = ""
  cur_row_rev = contig_read_table_tmp.at[cur_row, 'r1_is_reversed']

  for row in contig_read_table_tmp.index:
    next_seq = contig_read_table_tmp.at[row, 'seq']
    if (len(next_seq) >= len(longest_seq)):
      longest_seq = next_seq
      rev_chr_str = contig_read_table_tmp.at[row, 'r1_chr'] + 'r'
      longest_rev = contig_read_table_tmp.at[row, 'r1_is_reversed']
      if (contig_read_table_tmp.at[row, 'r1_is_reversed']):  # and not cur_row_rev):
        longest_seq = reverseComplement(next_seq, 1)
        print("longest reversed")
  return longest_seq


########################################
#
########################################
def orient_match_with_longest(contig_read_table_tmp, cur_row):
  longest_seq = ""
  cur_row_rev = contig_read_table_tmp.at[cur_row, 'r1_is_reversed']
  longest_seq_orient = cur_row_rev

  for row in contig_read_table_tmp.index:
    next_seq = contig_read_table_tmp.at[row, 'seq']
    if (len(next_seq) >= len(longest_seq)):
      longest_seq = next_seq
      longest_seq_orient = contig_read_table_tmp.at[row, 'r1_is_reversed']

  if (longest_seq_orient == cur_row_rev):
    return 1

  return 0


########################################
#
# Only doing a complement not reverse
#
########################################
def reverseComplement(seq, rflag):
  # Reverse sequence string
  if (rflag == 1):
    rseq = seq[::-1]
  else:
    rseq = seq

  rseq = rseq.replace('A', 'B')
  rseq = rseq.replace('T', 'A')
  rseq = rseq.replace('B', 'T')

  rseq = rseq.replace('C', 'B')
  rseq = rseq.replace('G', 'C')
  rseq = rseq.replace('B', 'G')

  return rseq


########################################
# Function to check contigs
########################################
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
  final_table = pd.DataFrame({'name': names2,
                              'Gene_1': "",
                              'IgTxCalled': 0})

  # Process each contig
  for currIndex in final_table.index:
    names = final_table.at[currIndex, 'name']
    print("\nNow Processing " + names)
    contig_table_by_read = contig_table[(contig_table.name == names)]
    # if contigs aligns to multiple locations
    count = len(contig_table_by_read.index)
    if (count > 1):  # only need Txs
      print("In test loop " + str(count))
      # initialize vars for loop
      loop_var = 1
      # track mininum mval
      min_mVal = 100000  # set to a large val
      # track variables for IG and Gene
      IG_found = ""
      Gene_found = ""
      strand = ""
      gene_len = 0
      IG_len = 0
      for contig_row in contig_table_by_read.index:
        contig_chr = contig_table_by_read.at[contig_row, 'r1_chr']
        contig_name = contig_table_by_read.at[contig_row, 'name']
        r_contig = "_" + contig_chr + "r"
        f_contig = "_" + contig_chr + "f"
        print("\n  " + contig_chr + "===" + r_contig + " ^^^^******* " + f_contig + "=" + name)
        if (f_contig in contig_name):
          strand = "Pos"
        if (r_contig in contig_name):
          strand = "Neg"
        contig_pos = contig_table_by_read.at[contig_row, 'r1_pos']
        # check if gene is in Tx region
        # @Sara Need to read updated bedfile as agreed by JJK

        mygene = isKnownTx(contig_chr, int(contig_pos), igregions)

        if (mygene != "" and contig_chr in ['chr2', 'chr14', 'chr22'] and IG_found == ""):
          IG_found = mygene
        if (mygene != "" and contig_chr in ['chr4', 'chr6', 'chr8', 'chr11', 'chr12', 'chr16',
                                            'chr20'] and Gene_found == ""):
          Gene_found = mygene

        # junctions in contig
        mVal = get_junctionBreak(contig_table_by_read.at[contig_row, 'r1_cigar'])
        juncbreak = get_junctionBreakNewMPos(contig_table_by_read.at[contig_row, 'r1_cigar'])
        if (mVal <= min_mVal):
          min_mVal = mVal

        # fine IG and Gene length
        if (Gene_found == mygene and Gene_found != ""):
          gene_len = mVal
        if (IG_found == mygene and IG_found != ""):
          IG_len = mVal

        # find longest seq for contig from bam to address  hard clipping
        contig_len = get_contigLength(contig_table_by_read.at[contig_row, 'r1_cigar'])

        # since bwa hardclips supplementary reads find the longest seq for contig
        longest_contig = get_longest_seqforcontig(contig_table_by_read,
                                                  contig_row)  # contig_table_by_read.at[contig_row,'seq']
        if (len(longest_contig) > len(contig_table_by_read.at[contig_row, 'seq'])):
          # ddi orientation switch?
          same_orient_as_longest = orient_match_with_longest(contig_table_by_read, contig_row)
          contig = longest_contig
        else:
          same_orient_as_longest = 1
          contig = contig_table_by_read.at[contig_row, 'seq']

        # updating contigs for reverse

        # switch only if orientation of longest os different than current
        # case 1 query contig is reversed by longest is not, since longest contig func reverses dont do anything
        if ((contig_table_by_read.at[contig_row, 'r1_is_reversed']) and (
            same_orient_as_longest == 0)):  # (rev_chr_str in contig_table_by_read.at[contig_row,'name'])):
          print("reverse contig junc" + str(juncbreak))
        # if longest was reversed and query is not
        elif ((not contig_table_by_read.at[contig_row, 'r1_is_reversed']) and (same_orient_as_longest == 0)):
          contig = reverseComplement(contig, 1)

        namesSplit = names.split(':')

        fastq = namesSplit[0]
        fastq = fastq + '_R1_001.fastq.gz'

        ffastq = fastq_path + "/" + fastq
        reads_sam = reads_sam_path + "/" + namesSplit[0] + "/" + namesSplit[0] + "_ReadstoContigs.bam"

        # extract 25mers left and right of contig
        jright = juncbreak + int(bp_win_at_junc)
        jleft = juncbreak - int(bp_win_at_junc)

        # extract chromosome coordinates 25mers around contig
        location = contig_table_by_read.at[contig_row, 'name'] + ":" + str(jleft) + "-" + str(jright)

        rr = contig[jleft:jright:1]
        rr_rev = reverseComplement(rr, 0)

        # get junction count  read 1 and read2
        r1_count = getReadsatjunction(rr, contig, ffastq) + getReadsatjunction(rr_rev, contig, ffastq)
        ffastq_2 = ffastq
        ffastq_2 = ffastq_2.replace("_R1", "_R2")

        r2_count = getReadsatjunction(rr, contig, ffastq_2) + getReadsatjunction(rr_rev, contig, ffastq_2)

        # get unique fragments at junction
        frag_count = getFragsatJunction_samtools(location, contig, reads_sam)

        # define df columns
        CIGAR = "cigar_" + str(loop_var)
        GENE = "Gene_" + str(loop_var)
        LEN = "Contig_length_" + str(loop_var)
        JUNC = "Contig_BP_" + str(loop_var)
        PERC = "percent_of_contig_at_Gene_" + str(loop_var)
        R1_C = "R1_reads_at_junc_" + str(loop_var)
        R2_C = "R2_reads_at_junc_" + str(loop_var)
        FRAG_C = "fragments_at_junc_" + str(loop_var)
        POS_START = "pos_" + str(loop_var) + "_start"
        POS_END = "pos_" + str(loop_var) + "_end"
        REV = "Contig_reversed" + str(loop_var)
        CHR = 'chr_' + str(loop_var)
        MVAL = 'aligned_length_' + str(loop_var)
        STRAND = 'strand_' + str(loop_var)
        SEQ = 'seq'

        if (loop_var == 1):
          final_table.at[[currIndex], SEQ] = contig

        print("IG_found =" + IG_found + " and Gene_found " + Gene_found)
        # if( mygene !=''):
        if (mVal >= window_size and mygene != ''):
          final_table.at[[currIndex], CIGAR] = contig_table_by_read.at[contig_row, 'r1_cigar']
          final_table.at[[currIndex], GENE] = mygene
          final_table.at[[currIndex], LEN] = len(contig)
          final_table.at[[currIndex], PERC] = mVal * 100 / contig_len
          final_table.at[[currIndex], MVAL] = mVal
          final_table.at[[currIndex], R1_C] = r1_count
          final_table.at[[currIndex], R2_C] = r2_count
          final_table.at[[currIndex], FRAG_C] = frag_count
          final_table.at[[currIndex], POS_START] = contig_table_by_read.at[contig_row, 'r1_pos']
          final_table.at[[currIndex], POS_END] = contig_table_by_read.at[contig_row, 'r1_pos'] + juncbreak  # + mVal
          final_table.at[[currIndex], REV] = contig_table_by_read.at[contig_row, 'r1_is_reversed']
          final_table.at[[currIndex], CHR] = contig_table_by_read.at[contig_row, 'r1_chr']
          final_table.at[[currIndex], JUNC] = juncbreak
          final_table.at[[currIndex], STRAND] = strand
          # increment counter

          loop_var = loop_var + 1
          if (count > 1 and IG_found != "" and Gene_found != "" and gene_len > window_size and IG_len > window_size):
            print("count " + str(count) + " ig len " + str(IG_len) + " gene len " + str(gene_len))
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
  # @Bryce Should be changed to a tmp location
  tmp_file = out_path + "/r1counts.txt"
  status = call("zcat " + fastq + " | grep " + region + " > " + tmp_file, shell=True)
  if status < 0:
    print("### Cat Command Failed....now exiting!!")
    sys.exit(-1)
  else:
    with open(tmp_file, "r") as myfile:
      data = myfile.readlines()
      r1_count = len(data)
  return r1_count


###############################################
#  Simple approach to find reads at junction
#  Need to be replaced with samtools in next
#  update. NOT CURRENTLY USED
###############################################
def getAllMappedReadsatJunction(region, contig, readsam):
  r1_count = 0
  r2_count = 0

  region_rev = reverseComplement(region)

  # @Bryce Should be changed to a tmp location
  tmp_file = out_path + "/r1counts.txt"
  status = call(
    "egrep \"" + region + "|" + region_rev + "\" " + readsam + " | awk '{ print $1 }' | sort | uniq  > " + tmp_file,
    shell=True)
  if status < 0:
    print("### Cat Command Failed....now exiting!!")
    sys.exit(-1)
  else:
    with open(tmp_file, "r") as myfile:
      data = myfile.readlines()
      r1_count = len(data)
  return r1_count


##########################################
#
#
#
##########################################
def get_genomic_bp(cigar, matchlen, orient, pos_start):
  gen_bp = pos_start
  # find location of match for region in Q
  Mval = 0
  Sval = 0
  Hval = 0
  Junc = 0
  if (cigar != ''):
    tmp_cigar = cigar
    tmp_cigar = tmp_cigar.replace('H', 'M')
    tmp_cigar = tmp_cigar.replace('S', 'M')
    tmp_cigar = tmp_cigar.replace('I', 'M')
    tmp_cigar = tmp_cigar.replace('D', 'M')

    # if clipping is sbefore M  then Mval =Mval -length
    sindex = cigar.find('S')
    hindex = cigar.find('H')
    mindex = cigar.find('M')
    total_len = 0
    # split
    cigar_list = tmp_cigar.split('M')
    for values in cigar_list:
      if (values != '' and cigar.find(str(values) + 'M') != -1):
        if (int(values) > Mval):
          Mval = int(values)
      # Generall only one S and/or H, keep first S or H
      if (values != '' and cigar.find(str(values) + 'S') != -1 and Sval == 0):
        Sval = int(values)
      if (values != '' and cigar.find(str(values) + 'H') != -1 and Hval == 0):
        Hval = int(values)
    # if(values != '' and cigar.find(str(values)+'M') !=-1 and Hval ==0):
    #      Mval = int(values)

    # if clipping is sbefore M  then Mval =Mval -length
    print("cigar=" + cigar)
    sindex = cigar.find('S')
    hindex = cigar.find('H')
    mindex = cigar.find('M')
    print(" h in " + str(hindex) + " m index " + str(mindex) + " s index" + str(sindex))
    if (hindex < mindex and hindex != -1):  # Hard clipping first
      print("h first")
      gen_bp = pos_start
    elif (sindex < mindex and sindex != -1):  # add soft clip region to np
      gen_bp = pos_start  # + Sval
      print(" s fiesr")
    else:  # if(mindex < sindex and index < hindex):
      gen_bp = pos_start + Mval
      print(" m first=" + str(Mval))
  # split
  # count
  # if H dont add
  return gen_bp


##########################################
#
#  Check for individual gene calls
#
##########################################
def check_gene_call(table, gene, nreads, min_con_len):
  # initialize return vars
  found_flag = 0
  count_gene = 0
  call = 0
  contig_length = 0
  ig_breakpoint = 0
  breakpoint = 0
  source = ''
  gene_overlap = 0
  gene_cigar = ''
  ig_cigar = ''
  ig_overlap = 0
  der_gene = ''
  der_ig = ''
  pos_strand_der = 0
  neg_strand_der = 0
  pos_strand_list = (0, 0, 0, 0, 0, 0, 0, 0, 0)
  neg_strand_list = (0, 0, 0, 0, 0, 0, 0, 0, 0)
  res_gene = ()

  # filter Tx by window and contig overlap
  # Case 1 : Gene_1 is Target gene
  if table.empty:
    res_gene = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    emp_list = []
    for t_elem in res_gene:
      emp_list.append(t_elem)
    return emp_list

  #updated code
  #start with quey gene
  gg =  'Gene_'+str(1)

  table_by_gene = table #_c[(table_c.Gene_1 == gene)]
  table_by_gene = table_by_gene[(table_by_gene.IgTxCalled == 1)]
  #loop all rows
  for row in table_by_gene.index:
    print("In Table Gene1 is target")
    print(table_by_gene.at[row, 'name'])
    print(table_by_gene.at[row, gg])

    total_aligned_length = 0
    ig_breakpoint = ''
    ig_aligned_length = 0
    gene_strand =''
    all_gene_cigar = ''
    all_ig_cigar = ''
    gene_overlap = 0
    ig_overlap = 0

    #get max splits for contigs
    max_contig_splits = 1
    for col in table_by_gene.columns:
      #print(col)
    #assuming  max 10 splits per contigs which is highly unlikely
      for tmp_idx in range(1,10):
        tmp_gene = 'Gene_' + str(tmp_idx)
        if(tmp_gene == col):
          if(tmp_idx >= max_contig_splits):
            max_contig_splits = tmp_idx
    #for  loop range add 1
    max_contig_splits = max_contig_splits + 1

    print("We have total contig splits "+str(max_contig_splits))
    #loop through all columns of each row
    for gindex in range(1,max_contig_splits):
      gg =  'Gene_'+str(gindex)

      print("Gene is "+gg+" "+str(gindex))
      gindex2 =  gindex

      #set indices
      pos_col = 'pos_'+str(gindex)+'_start'
      cigar_col= 'cigar_'+str(gindex)
      align_col = 'aligned_length_'+str(gindex)
      contig_flag_col = 'Contig_reversed'+str(gindex)
      frag_col = 'fragments_at_junc_'+str(gindex)
      contig_col = 'Contig_length_'+str(gindex)
      strand_col = 'strand_'+str(gindex)

      if((table_by_gene.at[row,gg]==gene)):

        #get values for target
        print(pos_col+cigar_col)
        pos_start = table_by_gene.at[row, pos_col]
        gene_cigar = table_by_gene.at[row, cigar_col]
        if(all_gene_cigar ==''):
           all_gene_cigar = table_by_gene.at[row, cigar_col]
        else:
           all_gene_cigar = all_gene_cigar + table_by_gene.at[row, cigar_col]
        matchlen = table_by_gene.at[row, align_col]
        count_gene = table_by_gene.at[row, frag_col]
        contig_length = table_by_gene.at[row, contig_col]
        if(matchlen > gene_overlap):
       	  gene_overlap = matchlen
        total_aligned_length = total_aligned_length + matchlen
        breakpoint = get_genomic_bp(gene_cigar, matchlen, table_by_gene.at[row, contig_flag_col], pos_start)
        gene_strand = table_by_gene.at[row, strand_col]
        #print(" bb="+str(breakpoint)+" "+str(matchlen)+" "+gene_cigar+" "+str(pos_start))

      #if gg is IG
      elif('IG' in str(table_by_gene.at[row,gg])):
        #print("Found IG at "+gg)
        ig_breakpoint = get_genomic_bp(table_by_gene.at[row, cigar_col], table_by_gene.at[row, align_col],
                                        table_by_gene.at[row, contig_flag_col], table_by_gene.at[row, pos_col])
        gene_tmp = table_by_gene.at[row,gg]
        if(table_by_gene.at[row, align_col] > ig_overlap):
          ig_overlap = table_by_gene.at[row, align_col]
        source = ig_dict[gene_tmp]
        if(table_by_gene.at[row, align_col] > ig_aligned_length):
          ig_aligned_length = table_by_gene.at[row, align_col]

        ig_cigar = table_by_gene.at[row, cigar_col]
        if(all_ig_cigar ==''):
           all_ig_cigar = table_by_gene.at[row, cigar_col]
        else:
           all_ig_cigar = all_ig_cigar + table_by_gene.at[row, cigar_col]

     #check all alignments for the contig
      while(gindex2 < max_contig_splits-1):
        gindex2=  gindex2 + 1
        gg2 = 'Gene_'+str(gindex2)
        #set indices
        pos_col2 = 'pos_'+str(gindex2)+'_start'
        cigar_col2 = 'cigar_'+str(gindex2)
        align_col2 = 'aligned_length_'+str(gindex2)
        contig_flag_col2 = 'Contig_reversed'+str(gindex2)
        frag_col2 = 'fragments_at_junc_'+str(gindex2)

        #if gene is same as target gene
        if((table_by_gene.at[row,gg] == table_by_gene.at[row,gg2]) and (table_by_gene.at[row,gg]==gene)):
          print("Found double match"+gg2)
          align_col2 = 'aligned_length_'+str(gindex2)
          total_aligned_length = total_aligned_length + table_by_gene.at[row, align_col2]
          gene_cigar = gene_cigar+';'+table_by_gene.at[row, cigar_col2]
          if(table_by_gene.at[row, frag_col2] > count_gene):
            count_gene = table_by_gene.at[row, frag_col2]
          if( table_by_gene.at[row, align_col2] > gene_overlap):
            gene_overlap =  table_by_gene.at[row, align_col2]

        #multiple alignments of IG
        if(('IG' in str(table_by_gene.at[row,gg2])) and (table_by_gene.at[row,gg]==table_by_gene.at[row,gg2])):
          print("Found IG at "+gg2)
          ig_breakpoint2 = get_genomic_bp(table_by_gene.at[row, cigar_col2], table_by_gene.at[row, align_col2],
                                        table_by_gene.at[row, contig_flag_col2], table_by_gene.at[row, pos_col2])
          if(table_by_gene.at[row, align_col2] > ig_aligned_length):
            ig_aligned_length = table_by_gene.at[row, align_col2]
          ig_breakpoint = ig_breakpoint+';'+ig_breakpoint2
          ig_cigar = ig_cigar+';'+table_by_gene.at[row, cigar_col2]
    #test onky since we are missing fastq
    print("************************\n"+str(total_aligned_length) + " "+str(ig_aligned_length)+"*****************\n")
    #gene overlap is max alignment and total is total region aligned
    #switching total to gene to filer by max.
    #we keep both in case we would like to  get total length in later revisions
    total_aligned_length = gene_overlap
    #ig_overlap = ig_aligned_length
    if (total_aligned_length >= min_con_len and ig_aligned_length >= min_con_len and count_gene >= nreads):
      call = 1
      print("Tx pass")
      if (gene_strand == 'Pos'):
        print("Found pos")
        pos_strand_der = 1
        pos_strand_list = (
        pos_strand_der, contig_length, gene_overlap, all_gene_cigar, ig_overlap, all_ig_cigar, breakpoint, ig_breakpoint,
        count_gene)
      elif (gene_strand == 'Neg'):
        neg_strand_der = 1
        print("found neg")
        neg_strand_list = (
        neg_strand_der, contig_length, gene_overlap, all_gene_cigar, ig_overlap, all_ig_cigar, breakpoint, ig_breakpoint,
        count_gene)

##end of updated code

  print("call = " + str(call))
  # Collapse into one list
  if (call == 1):
    #print(pos_strand_list)
    #print(neg_strand_list)
    res_gene = [call, source]
    for x in pos_strand_list:
      res_gene.append(x)
    for x in neg_strand_list:
      res_gene.append(x)
    found_flag = 1
  # check for Gene 2
  if (found_flag == 0):
    res_gene = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
  print("res gene")
  print(res_gene)

  full_list = []
  for t_elem in res_gene:
    full_list.append(t_elem)
  return full_list


##########################################
# Generate summary table
##########################################
def gen_summ_table(filt_table, results_table, nreads, min_con_len, sample):
  # prepare final table
  list_of_genes = ['NSD2', 'CCND3', 'MYC', 'MAFA', 'CCND1', 'CCND2', 'MAF', 'MAFB']
  list_of_features = ['Call', 'Ig_Loci', 'PosStrand_Derivative', 'PosStrand_Contig_Length',
                      'PosStrand_Target_Match', 'PosStrand_Target_Match_Cigars', 'PosStrand_Ig_Match',
                      'PosStrand_Ig_Match_Cigars', 'PosStrand_Target_Breakpoint', 'PosStrand_Ig_Breakpoint',
                      'PosStrand_Junction_Fragments', 'NegStrand_Derivative',
                      'NegStrand_Contig_Length', 'NegStrand_Target_Match', 'NegStrand_Target_Match_Cigars',
                      'NegStrand_Ig_Match', 'NegStrand_Ig_Match_Cigars',
                      'NegStrand_Target_Breakpoint', 'NegStrand_Ig_Breakpoint', 'NegStrand_Junction_Fragments']

  column_names = ['Specimen']
  # Make list of column_names:
  for gene in list_of_genes:
    for feature in list_of_features:
      header = '_'.join([gene, feature])
      column_names.append(header)

  # check for NSD2
  nsd2_call = check_gene_call(filt_table, "NSD2", nreads, min_con_len)
  ccnd3_call = check_gene_call(filt_table, "CCND3", nreads, min_con_len)
  myc_call = check_gene_call(filt_table, "MYC", nreads, min_con_len)
  mafa_call = check_gene_call(filt_table, "MAFA", nreads, min_con_len)
  ccnd1_call = check_gene_call(filt_table, "CCND1", nreads, min_con_len)
  ccnd2_call = check_gene_call(filt_table, "CCND2", nreads, min_con_len)
  maf_call = check_gene_call(filt_table, "MAF", nreads, min_con_len)
  mafb_call = check_gene_call(filt_table, "MAFB", nreads, min_con_len)

  results = ((sample,), nsd2_call, ccnd3_call, myc_call, mafa_call, ccnd1_call, ccnd2_call, maf_call, mafb_call)

  #print(mafb_call)
  con_results = []
  con_results.append(sample)
  # con_results=con_results+nsd2_call+ccnd3_call+myc_call+mafa_call+ccnd1_call+ccnd2_call+maf_call+mafb_call
  for x in nsd2_call:
    con_results.append(x)  # ,ignore_index=True,verify_integrity=False)
  for x in ccnd3_call:
    con_results.append(x)
  for x in myc_call:
    con_results.append(x)
  for x in mafa_call:
    con_results.append(x)
  for x in ccnd1_call:
    con_results.append(x)
  for x in ccnd2_call:
    con_results.append(x)
  for x in maf_call:
    con_results.append(x)
  for x in mafb_call:
    con_results.append(x)

  print("full list")
  #print(con_results)
  translocationsTable = pd.DataFrame(columns=column_names)
  translocationsTable.loc[len(translocationsTable)] = con_results
  return translocationsTable


###############################################
#  Simple approach to find reads at junction
#  Need to be replaced with samtools in next
#  update
###############################################
def getFragsatJunction_samtools(location, contig, bam):
  r1_count = 0
  r2_count = 0
  # grep region1 from fastq
  # @Bryce Should be changed to a tmp location
  tmp_file = out_path + "/r1counts.txt"

  print("samtools view " + bam + " " + location + " | awk '{ print $1 }' | sort | uniq | wc -l ")
  status = call("samtools view " + bam + " " + location + " | awk '{ print $1 }' | sort | uniq  > " + tmp_file,
                shell=True)
  if status < 0:
    print("### Cat Command Failed....now exiting!!")
    sys.exit(-1)
  else:
    with open(tmp_file, "r") as myfile:
      data = myfile.readlines()
      r1_count = len(data)
  return r1_count


###################################################
#
# END FUNCTION DEFINITIONS
#
###################################################
print("Importing Data")

print("Reading SAM File")

# Create python object for the SAM file

# get Ig/parter bed file
ig_headers = ['chr', 'start', 'stop', 'name']
igregions = pd.read_csv(bed_file, sep="\t", header=None, names=ig_headers)

results_full_table = pd.DataFrame()
filt_table = pd.DataFrame()
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
print(results_table)
if not results_table.empty:
  filt_table = results_table[
    results_table.Gene_1.isin(['IGH', 'IGL', 'NSD2', 'CCND1', 'CCND2', 'CCND3', 'IGK', 'MYC', 'MAFA', 'MAFB', 'MAF'])]
  print("\nWriting outfile ****************")

  out_file = out_path + "/ContigResults.txt"
  out_file_filt = out_path + "/FilteredContigResults.txt"
  results_table.to_csv(out_file, sep="\t", index=False, na_rep=0, float_format='%.0f')

  if not filt_table.empty:
    filt_table.to_csv(out_file_filt, sep="\t", index=False, na_rep=0, float_format='%.0f')

# generate summary table
summ_table = gen_summ_table(filt_table, results_table, min_reads, window_size, sample_name)
#print(summ_table)
out_file_summ = out_path + "/DEX_IgTx_GA_Summary.txt"
summ_table.to_csv(out_file_summ, sep="\t", index=False, na_rep=0, float_format='%.0f')
print("Test Done")
