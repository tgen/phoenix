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
    name.append(file_name + ":" + read.query_name)  # +"_"+str(index))
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
    print(str(r1_cigar) + " and  " + str(name) + " mapq " + str(r1_mapq))
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


# Function to check contigs
def check_contigs(contig_table, fastq_list, igregions, window_size=400, min_reads=5):
  # for i in range(len(contig_table)) :
  #    print(contig_table.loc[i, 'r1_cigar'])
  list_names = []

  for row in contig_table.index:
    r1_cigar = contig_table.at[row, 'r1_cigar']
    name = contig_table.at[row, 'name']
    if (name not in list_names):
      list_names.append(name)
    # split cigar by M,S,H
    # print("LOOP.....")
    print(name + "=====r1_cigar=====")  # + r1_cigar)
    # mindex = r1_cigar.find('M')
    # sindex = r1_cigar.find('S')
    # hindex = r1_cigar.find('H')
    # get index of MSH
    Mval = 0
    if (r1_cigar != ''):
      tmp_cigar = r1_cigar
      tmp_cigar = tmp_cigar.replace('H', 'M')
      tmp_cigar = tmp_cigar.replace('S', 'M')
      tmp_cigar = tmp_cigar.replace('I', 'M')
      tmp_cigar = tmp_cigar.replace('D', 'M')
      print(tmp_cigar)
      cigar_list = tmp_cigar.split('M')
      for values in cigar_list:
        print("vals = " + values)
        if (values != '' and r1_cigar.find(str(values) + 'M') != -1):
          if (int(values) > Mval):
            Mval = int(values)
      if (Mval < window_size):
        print("Value lower than thr " + str(Mval))
  print("\n=====\n*************\n========\n")
  # next check if two names are on multiple lines
  for names in list_names:
    # extract by name
    print(names)
    contig_table_by_read = contig_table[(contig_table.name == names)]
    # if contigs aligns to multiple locations
    count = len(contig_table_by_read.index)
    if (count > 1):
      print(contig_table_by_read)
      # get breakpoint on contig?
      # check their locations
      for contig_row in contig_table_by_read.index:
        contig_chr = contig_table_by_read.at[contig_row, 'r1_chr']
        contig_pos = contig_table_by_read.at[contig_row, 'r1_pos']
        mygene = isKnownTx(contig_chr, int(contig_pos), igregions)
        print("contig_chr=" + contig_chr + "contig_pos " + str(contig_pos) + mygene)

        juncbreak = 334  # getBreakpoint()
        contig = contig_table_by_read.at[contig_row, 'seq']
        namesSplit = names.split(':')
        fastq = namesSplit[0]
        fastq = fastq.replace('Trinity_sorted.bam', '001.fastq.gz')
        # fastq = namesSplit #namesSplit[0]+ "_001_fastq.gz"
        ffastq = ''
        for fqs in fastq_list:
          if (fastq in fqs):
            ffastq = fqs
        print("\nfastq full is " + str(ffastq))
        # fastq_file = fqs[]
        print("\nfastq name is " + fastq)
        # print("\nfastq full is " + fastq)
        val = getReadsatjunction(juncbreak, contig, ffastq)


# Function to check if given  location (chr:pos) is in list of IgTx regions
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


def getReadsatjunction(juncbreak, contig, fastq):
  # extract 25mers left and right of contig
  jright = juncbreak + 25
  jleft = juncbreak - 25

  r1 = contig[juncbreak: jright: 1]
  r2 = contig[jleft: juncbreak: 1]

  # grep region1 from fastq
  status = call("zcat " + fastq + " | grep " + r1 + " | wc -l > r1counts.txt")
  if status < 0:
    print("### Cat Command Failed....now exiting!!")
    sys.exit(-1)
  else:
    with open("r1counts.txt", "r") as myfile:
      data = myfile.readlines()
      r1_count = int(data)
      print("read count is " + r1_count)

  # grep region2 from fastq
  # val2 = system("zcat "+ fastq + " | grep " + r2 + "| wc -l > r1counts.txt")

  # if count of region1 is  greater than > min_reads and count  of region2 > min_reads
  if (r1_count >= 1):  # min_reads and val2 >=  min_reads):
    return 1
  else:
    return 0


# Function to score and count dynamic window sizes from read1 counts
def call_r1_windows(result_row, input_table, orientation, windowname, windowchecksize=1000, readorder="r1"):
  position_column = readorder + "_pos"
  maxwindowcount = 0
  maxwindowwidth = 0
  maxwindowlocation = 0
  count = len(input_table.index)
  for row in input_table.index:
    # Get PositionA value, calculate +/- window, capture rows within the +/- window
    # Create variable for the read specific position field
    position_current = input_table.at[row, position_column]
    position_negative = (position_current - windowchecksize)
    position_positive = (position_current + windowchecksize)
    position_table = input_table[(input_table.r1_pos > position_negative) & (input_table.r1_pos < position_positive)]
    position_count = len(position_table.index)
    # Create Tracking Variables
    # if the current window cound is greater than existing Max, update Max and Next accordingly
    if position_count > maxwindowcount:
      nextlargestwindowcount = maxwindowcount
      maxwindowcount = position_count
      nextwindowwidth = maxwindowwidth
      nextwindowlocation = maxwindowlocation
      max_window_top = position_table.head(1)
      max_window_bottom = position_table.tail(1)
      if readorder == "r1":
        max_window_top_position = max_window_top.iat[0, 2]  # was 8
        max_window_bottom_position = max_window_bottom.iat[0, 2]  # was 8
      elif readorder == "r2":
        max_window_top_position = max_window_top.iat[0, 10]  # was 13
        max_window_bottom_position = max_window_bottom.iat[0, 10]  # was 13

      print("SN XXX window wdith " + str(max_window_bottom_position) + " top pos " + str(max_window_top_position))
      maxwindowwidth = max_window_bottom_position - max_window_top_position

      # Determine the proximal breakpoint location based on read orientation
      if orientation == "FR":
        maxwindowlocation = max_window_bottom_position
      elif orientation == "RF":
        maxwindowlocation = max_window_top_position
      elif orientation == "FF":
        maxwindowlocation = max_window_bottom_position
      elif orientation == "RR":
        maxwindowlocation = max_window_top_position

      # if the two bundle counts will be greater than total count reset the NextLargest to 0
      if (nextlargestwindowcount + maxwindowcount) > count:
        nextlargestwindowcount = 0
        nextwindowwidth = 0
        nextwindowlocation = 0
      windowend = position_positive
    # if the second bundle in order is greater than the current next but less than the max then capture current as a next
    elif position_count > nextlargestwindowcount and position_count <= maxwindowcount and position_negative >= windowend:
      nextlargestwindowcount = position_count
      max_window_top = position_table.head(1)
      max_window_bottom = position_table.tail(1)
      if readorder == "r1":
        max_window_top_position = max_window_top.iat[0, 2]  # was 8
        max_window_bottom_position = max_window_bottom.iat[0, 2]  # was 8
      elif readorder == "r2":
        max_window_top_position = max_window_top.iat[0, 10]  # was 13
        max_window_bottom_position = max_window_bottom.iat[0, 10]  # was 13

      nextwindowwidth = max_window_bottom_position - max_window_top_position

      # Determine the proximal breakpoint location based on read orientation
      if orientation == "FR":
        mnextwindowlocation = max_window_bottom_position
      elif orientation == "RF":
        nextwindowlocation = max_window_top_position
      elif orientation == "FF":
        nextwindowlocation = max_window_bottom_position
      elif orientation == "RR":
        nextwindowlocation = max_window_top_position

  # Define variables
  MAX_WINDOW_COUNT = orientation + "_" + windowname + "_" + readorder + '_maxwindowcount'
  MAX_WINDOW_WIDTH = orientation + "_" + windowname + "_" + readorder + '_maxwindowwidth'
  MAX_WINDOW_LOCATION = orientation + "_" + windowname + "_" + readorder + '_maxwindowlocation'
  NEXT_WINDOW_COUNT = orientation + "_" + windowname + "_" + readorder + '_nextlargestwindowcount'
  NEXT_WINDOW_WIDTH = orientation + "_" + windowname + "_" + readorder + '_nextwindowwidth'
  NEXT_WINDOW_LOCATION = orientation + "_" + windowname + "_" + readorder + '_nextwindowlocation'

  # Add values to data table
  tf_table.at[[result_row], MAX_WINDOW_COUNT] = maxwindowcount
  tf_table.at[[result_row], MAX_WINDOW_WIDTH] = maxwindowwidth
  tf_table.at[[result_row], MAX_WINDOW_LOCATION] = maxwindowlocation
  tf_table.at[[result_row], NEXT_WINDOW_COUNT] = nextlargestwindowcount
  tf_table.at[[result_row], NEXT_WINDOW_WIDTH] = nextwindowwidth
  tf_table.at[[result_row], NEXT_WINDOW_LOCATION] = nextwindowlocation

  # Message the results summary
  print(orientation + " " + windowname + "_" + readorder +
        " Max - " + str(maxwindowcount) + ", " + str(maxwindowwidth) + ", " + str(maxwindowlocation))
  print(orientation + " " + windowname + "_" + readorder +
        " Next - " + str(nextlargestwindowcount) + ", " + str(nextwindowwidth) + ", " + str(nextwindowlocation))


######
# Function to score and count dynamic window sizes from read1 counts
def call_r2_windows(result_row, input_table, orientation, windowname, windowchecksize=1000, readorder="r2"):
  position_column = readorder + "_pos"
  maxwindowcount = 0
  maxwindowwidth = 0
  maxwindowlocation = 0
  count = len(input_table.index)
  for row in input_table.index:
    # Get PositionA value, calculate +/- window, capture rows within the +/- window
    # Create variable for the read specific position field
    position_current = input_table.at[row, position_column]
    position_negative = (position_current - windowchecksize)
    position_positive = (position_current + windowchecksize)
    position_table = input_table[(input_table.r2_pos > position_negative) & (input_table.r2_pos < position_positive)]
    position_count = len(position_table.index)
    # Create Tracking Variables
    # if the current window cound is greater than existing Max, update Max and Next accordingly
    if position_count > maxwindowcount:
      nextlargestwindowcount = maxwindowcount
      maxwindowcount = position_count
      nextwindowwidth = maxwindowwidth
      nextwindowlocation = maxwindowlocation
      max_window_top = position_table.head(1)
      max_window_bottom = position_table.tail(1)
      if readorder == "r1":
        max_window_top_position = max_window_top.iat[0, 2]  # was 8
        max_window_bottom_position = max_window_bottom.iat[0, 2]  # was 8
      elif readorder == "r2":
        max_window_top_position = max_window_top.iat[0, 10]  # was 13
        max_window_bottom_position = max_window_bottom.iat[0, 10]  # was 13

      maxwindowwidth = max_window_bottom_position - max_window_top_position

      # Determine the proximal breakpoint location based on read orientation
      if orientation == "FR":
        maxwindowlocation = max_window_top_position
      elif orientation == "RF":
        maxwindowlocation = max_window_bottom_position
      elif orientation == "FF":
        maxwindowlocation = max_window_top_position
      elif orientation == "RR":
        maxwindowlocation = max_window_bottom_position

      # if the two bundle counts will be greater than total count reset the NextLargest to 0
      if (nextlargestwindowcount + maxwindowcount) > count:
        nextlargestwindowcount = 0
        nextwindowwidth = 0
        nextwindowlocation = 0
      windowend = position_positive
    # if the second bundle in order is greater than the current next but less than the max then capture current as a next
    elif position_count > nextlargestwindowcount and position_count <= maxwindowcount and position_negative >= windowend:
      nextlargestwindowcount = position_count
      max_window_top = position_table.head(1)
      max_window_bottom = position_table.tail(1)
      if readorder == "r1":
        max_window_top_position = max_window_top.iat[0, 2]  # was 8
        max_window_bottom_position = max_window_bottom.iat[0, 2]  # was 8
      elif readorder == "r2":
        max_window_top_position = max_window_top.iat[0, 10]  # was 13
        max_window_bottom_position = max_window_bottom.iat[0, 10]  # was 13

      nextwindowwidth = max_window_bottom_position - max_window_top_position

      # Determine the proximal breakpoint location based on read orientation
      if orientation == "FR":
        nextwindowlocation = max_window_top_position
      elif orientation == "RF":
        nextwindowlocation = max_window_bottom_position
      elif orientation == "FF":
        nextwindowlocation = max_window_top_position
      elif orientation == "RR":
        nextwindowlocation = max_window_bottom_position

  # Define variables
  MAX_WINDOW_COUNT = orientation + "_" + windowname + "_" + readorder + '_maxwindowcount'
  MAX_WINDOW_WIDTH = orientation + "_" + windowname + "_" + readorder + '_maxwindowwidth'
  MAX_WINDOW_LOCATION = orientation + "_" + windowname + "_" + readorder + '_maxwindowlocation'
  NEXT_WINDOW_COUNT = orientation + "_" + windowname + "_" + readorder + '_nextlargestwindowcount'
  NEXT_WINDOW_WIDTH = orientation + "_" + windowname + "_" + readorder + '_nextwindowwidth'
  NEXT_WINDOW_LOCATION = orientation + "_" + windowname + "_" + readorder + '_nextwindowlocation'

  # Add values to data table
  tf_table.at[[result_row], MAX_WINDOW_COUNT] = maxwindowcount
  tf_table.at[[result_row], MAX_WINDOW_WIDTH] = maxwindowwidth
  tf_table.at[[result_row], MAX_WINDOW_LOCATION] = maxwindowlocation
  tf_table.at[[result_row], NEXT_WINDOW_COUNT] = nextlargestwindowcount
  tf_table.at[[result_row], NEXT_WINDOW_WIDTH] = nextwindowwidth
  tf_table.at[[result_row], NEXT_WINDOW_LOCATION] = nextwindowlocation

  # Message the results summary
  print(orientation + " " + windowname + "_" + readorder +
        " Max - " + str(maxwindowcount) + ", " + str(maxwindowwidth) + ", " + str(maxwindowlocation))
  print(orientation + " " + windowname + "_" + readorder +
        " Next - " + str(nextlargestwindowcount) + ", " + str(nextwindowwidth) + ", " + str(nextwindowlocation))


######

# Function to call events
def call_structural_event(result_row, orientation,
                          window1_chr, window1_start, window1_end,
                          window2_chr, window2_start, window2_end):
  # Make tables for both possible derivatives
  # Read1 aligned to forward strand
  orientation_result_table = result_table[result_table.pair_orientation == orientation]
  orientation_count = len(orientation_result_table.index)
  print(orientation + " count - " + str(orientation_count))

  if orientation_count != 0:
    # Calculate the distribution of discordant reads on each window end
    # WINDOW-1 against read1
    orientation_result_table_window1 = orientation_result_table[
      (orientation_result_table.r1_chr_int == window1_chr) &
      (orientation_result_table.r1_pos >= window1_start) &
      (orientation_result_table.r1_pos <= window1_end)]
    orientation_result_table_window1_count = len(orientation_result_table_window1.index)
    print(orientation + " Window1 count - " + str(orientation_result_table_window1_count))

    if orientation_result_table_window1_count != 0:
      orientation_window1_table_r1 = orientation_result_table_window1.sort_values(['r1_chr_int', 'r1_pos'],
                                                                                  ascending=[1, 1])
      orientation_window1_top_r1 = orientation_window1_table_r1.head(1)
      orientation_window1_bottom_r1 = orientation_window1_table_r1.tail(1)
      orientation_window1_top_position_r1 = orientation_window1_top_r1.iat[0, 2]  # was 8
      orientation_window1_bottom_position_r1 = orientation_window1_bottom_r1.iat[0, 2]  # was 8
      orientation_window1_size_r1 = int(orientation_window1_bottom_position_r1) - int(
        orientation_window1_top_position_r1)
      print(orientation + " Window1_r1: " +
            str(orientation_window1_size_r1) +
            "bp " +
            str(orientation_window1_top_position_r1) +
            "-" +
            str(orientation_window1_bottom_position_r1))
      # Call the window calling function
      call_r1_windows(result_row=result_row, input_table=orientation_window1_table_r1,
                      orientation=orientation, windowname="win1")

      # WINDOW-1 against read2
      orientation_window1_table_r2 = orientation_result_table_window1.sort_values(['r2_chr_int', 'r2_pos'],
                                                                                  ascending=[1, 1])
      orientation_window1_top_r2 = orientation_window1_table_r2.head(1)
      orientation_window1_bottom_r2 = orientation_window1_table_r2.tail(1)
      orientation_window1_top_position_r2 = orientation_window1_top_r2.iat[0, 10]  # was 13
      orientation_window1_bottom_position_r2 = orientation_window1_bottom_r2.iat[0, 10]  # was 13
      orientation_window1_size_r2 = int(orientation_window1_bottom_position_r2) - int(
        orientation_window1_top_position_r2)
      print(orientation + " Window1_r2: " +
            str(orientation_window1_size_r2) +
            "bp " +
            str(orientation_window1_top_position_r2) +
            "-" +
            str(orientation_window1_bottom_position_r2))
      # Call the window calling function
      call_r2_windows(result_row=result_row, input_table=orientation_window1_table_r2,
                      orientation=orientation, windowname="win1")
    else:
      orientation_window1_size_r1 = 0
      orientation_window1_size_r2 = 0

    # WINDOW-2 against read1
    orientation_result_table_window2 = orientation_result_table[
      (orientation_result_table.r1_chr_int == window2_chr) &
      (orientation_result_table.r1_pos >= window2_start) &
      (orientation_result_table.r1_pos <= window2_end)]
    orientation_result_table_window2_count = len(orientation_result_table_window2.index)
    print(orientation + " Window2 count - " + str(orientation_result_table_window2_count))

    if orientation_result_table_window2_count != 0:
      orientation_window2_table_r1 = orientation_result_table_window2.sort_values(['r1_chr_int', 'r1_pos'],
                                                                                  ascending=[1, 1])
      orientation_window2_top_r1 = orientation_window2_table_r1.head(1)
      orientation_window2_bottom_r1 = orientation_window2_table_r1.tail(1)
      orientation_window2_top_position_r1 = orientation_window2_top_r1.iat[0, 2]  # was 8
      orientation_window2_bottom_position_r1 = orientation_window2_bottom_r1.iat[0, 2]  # was 8
      orientation_window2_size_r1 = int(orientation_window2_bottom_position_r1) - int(
        orientation_window2_top_position_r1)
      print(orientation + " Window2_r1: " +
            str(orientation_window2_size_r1) +
            "bp " +
            str(orientation_window2_top_position_r1) +
            "-" +
            str(orientation_window2_bottom_position_r1))
      # Call the window calling function
      print(orientation_window2_table_r1)
      orientation_window2_table_r1.to_csv("orientation_window2_table_r1.txt", sep="\t", index=False,
                                          float_format='%.0f')
      call_r1_windows(result_row=result_row, input_table=orientation_window2_table_r1,
                      orientation=orientation, windowname="win2")

      # WINDOW-2 against read2
      orientation_window2_table_r2 = orientation_result_table_window2.sort_values(['r2_chr_int', 'r2_pos'],
                                                                                  ascending=[1, 1])
      orientation_window2_top_r2 = orientation_window2_table_r2.head(1)
      orientation_window2_bottom_r2 = orientation_window2_table_r2.tail(1)
      orientation_window2_top_position_r2 = orientation_window2_top_r2.iat[0, 10]  # was 13
      orientation_window2_bottom_position_r2 = orientation_window2_bottom_r2.iat[0, 10]  # was 13
      orientation_window2_size_r2 = int(orientation_window2_bottom_position_r2) - int(
        orientation_window2_top_position_r2)
      print(orientation + " Window2_r2: " +
            str(orientation_window2_size_r2) +
            "bp " +
            str(orientation_window2_top_position_r2) +
            "-" +
            str(orientation_window2_bottom_position_r2))
      # Call the window calling function
      call_r2_windows(result_row=result_row, input_table=orientation_window2_table_r2,
                      orientation=orientation, windowname="win2")
    else:
      orientation_window2_size_r1 = 0
      orientation_window2_size_r2 = 0
  else:
    orientation_result_table_window1_count = 0
    orientation_window1_size_r1 = 0
    orientation_window1_size_r2 = 0
    orientation_result_table_window2_count = 0
    orientation_window2_size_r1 = 0
    orientation_window2_size_r2 = 0

  # Set variables so data is added to the data table
  ORIENTATION_COUNT = orientation + '_discordantFrag_Count'
  ORIENTATION_WIN1_COUNT = orientation + '_win1_count'
  ORIENTATION_WIN1_R1_SIZE = orientation + '_win1_r1_size'
  ORIENTATION_WIN1_R2_SIZE = orientation + '_win1_r2_size'
  ORIENTATION_WIN2_COUNT = orientation + '_win2_count'
  ORIENTATION_WIN2_R1_SIZE = orientation + '_win2_r1_size'
  ORIENTATION_WIN2_R2_SIZE = orientation + '_win2_r2_size'

  # Add values to data table
  tf_table.at[[row], ORIENTATION_COUNT] = orientation_count
  tf_table.at[[row], ORIENTATION_WIN1_COUNT] = orientation_result_table_window1_count
  tf_table.at[[row], ORIENTATION_WIN1_R1_SIZE] = orientation_window1_size_r1
  tf_table.at[[row], ORIENTATION_WIN1_R2_SIZE] = orientation_window1_size_r2
  tf_table.at[[row], ORIENTATION_WIN2_COUNT] = orientation_result_table_window2_count
  tf_table.at[[row], ORIENTATION_WIN2_R1_SIZE] = orientation_window2_size_r1
  tf_table.at[[row], ORIENTATION_WIN2_R2_SIZE] = orientation_window2_size_r2


# END FUNCTION DEFINITIONS

print("Importing TophatFusion Table")
# Import Collapsed Tophat-Fusion table
# tf_table_file = tf_file = sys.argv[1]
# tf_table = pd.read_csv(tf_table_file, sep="\t",keep_default_na=False, na_filter=False) #,na_values=[''])
# tf_table.replace('', '.',regex=True)

print("Reading SAM File")
# Create python object for the SAM file
bam_list_string = sys.argv[1]
bam_list = bam_list_string.split(',')
# get fastq
fastq_str = sys.argv[2]
fastq_list = fastq_str.split(',')
# get Ig/parter bed file
ig_headers = ['chr', 'start', 'stop', 'name']
igregions = pd.read_csv(sys.argv[4], sep="\t", header=None, names=ig_headers)

full_table = pd.DataFrame()
index_sam = 0
for assm_sam_file in bam_list:
  samfile = pysam.AlignmentFile(assm_sam_file, "r")

  print("Creating SAM File Table")
  # Call Function to convert SAM file into pandas dataframe

  table = bam_to_df(samfile, file_name=assm_sam_file)
  check_contigs(table, fastq_list, igregions, window_size, min_reads)
  if (index_sam == 0):
    full_table = table
  else:
    print("Merging Regional bam")
    frames = [full_table, table]
    full_table = pd.concat(frames)
    # full_table = full_table.merge(table, left_on='name', right_on='name')
  index_sam = index_sam + 1

full_table.reset_index(drop=True)
full_table.set_index('name')  # ["name"])
# generate common bam
# Sort the bam
# write bam to file
# pysam.sort()
print("****************")
full_table.to_csv("MergedTableAllRegions.txt", sep="\t", index=False, na_rep=0, float_format='%.0f')
print(full_table)
print("Test Done")
