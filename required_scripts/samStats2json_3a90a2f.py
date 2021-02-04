#!/usr/bin/env python3

import re
import json
import argparse
import os
import sys
import csv
from copy import deepcopy

currentDirectory = os.getcwd()

# defines valid filetypes users can enter as an argument. As more files are supported,
#       new filetypes must be added here.
fileTypes = {
    "picard_alignment_summary_metrics": {
        "picard": True,
        "picard_metrics": True,
        "picard_histogram": False,
        "metrics_rows_per_level": 3,
        "metrics_cols_to_concat": 1,
        "help": "\"bam_name.bam.alignment_summary_metrics\" file generated from Picard."
    },
    "picard_hs_metrics": {
        "picard": True,
        "picard_metrics": True,
        "picard_histogram": True,
        "hist_col_per_level": 2,
        "metrics_rows_per_level": 1,
        "metrics_cols_to_concat": 1,
        "help": "\"bam_name.bam.hs_metrics\" file generated from Picard."
    },
    "picard_insert_size_metrics": {
        "picard": True,
        "picard_metrics": True,
        "picard_histogram": True,
        "hist_col_per_level": 1,
        "metrics_rows_per_level": 1,
        "metrics_cols_to_concat": 1,
        "help": "\"bam_name.bam.insert_size_metrics\" file generated from Picard."
    },
    "picard_mark_duplicates_metrics": {
        "picard": True,
        "picard_metrics": True,
        "picard_histogram": True,
        "hist_col_per_level": 4,
        "metrics_rows_per_level": 1,
        "metrics_cols_to_concat": 1,
        "help": "\"bam_name.bam.insert_size_metrics\" file generated from Picard."
    },
    "picard_oxog_metrics": {
        "picard": True,
        "picard_metrics": True,
        "picard_histogram": False,
        "metrics_rows_per_level": 16,
        "metrics_cols_to_concat": 1,
        "help": "\"bam_name.bam.insert_size_metrics\" file generated from Picard."
    },
    "picard_quality_yield_metrics": {
        "picard": True,
        "picard_metrics": True,
        "picard_histogram": False,
        "metrics_rows_per_level": 1,
        "metrics_cols_to_concat": 1,
        "help": "\"bam_name.bam.quality_yield_metrics\" file generated from Picard."
    },
    "picard_rna_metrics": {
        "picard": True,
        "picard_metrics": True,
        "picard_histogram": True,
        "hist_col_per_level": 1,
        "metrics_rows_per_level": 1,
        "metrics_cols_to_concat": 1,
        "help": "\"bam_name.bam.rna_metrics\" file generated from GATK CollectRnaSeqMetrics."
    },
    "picard_wgs_wnzc_metrics": {
        "picard": True,
        "picard_metrics": True,
        "picard_histogram": True,
        "hist_col_per_level": 2,
        "metrics_rows_per_level": 2,
        "metrics_cols_to_concat": 1,
        "help": "\"bam_name.bam.wgs_wnzc_metrics\" file generated from GATK CollectWgsMetricsWithNonZeroCoverage."
    },
    "picard_bait_bias_detail_metrics": {
        "picard": True,
        "picard_metrics": True,
        "picard_histogram": False,
        "metrics_rows_per_level": 192,
        "metrics_cols_to_concat": 3,
        "help": "\"bam_name.bam.wgs_metrics\" file generated from GATK CollectMultipleMetrics."
    },
    "picard_bait_bias_summary_metrics": {
        "picard": True,
        "picard_metrics": True,
        "picard_histogram": False,
        "metrics_rows_per_level": 12,
        "metrics_cols_to_concat": 2,
        "help": "\"bam_name.bam.wgs_metrics\" file generated from GATK CollectMultipleMetrics."
    },
    "picard_base_distribution_by_cycle_metrics": {
        "picard": True,
        "picard_metrics": True,
        "picard_histogram": False,
        "metrics_rows_per_level": "All",
        "metrics_cols_to_concat": 2,
        "help": "\"bam_name.bam.wgs_metrics\" file generated from GATK CollectMultipleMetrics."
    },
    "picard_error_summary_metrics": {
        "picard": True,
        "picard_metrics": True,
        "picard_histogram": False,
        "metrics_rows_per_level": 6,
        "metrics_cols_to_concat": 2,
        "metrics_cols_to_skip": 2,
        "help": "\"bam_name.bam.wgs_metrics\" file generated from GATK CollectWgsMetrics."
    },
    "picard_gc_bias_summary_metrics": {
        "picard": True,
        "picard_metrics": True,
        "picard_histogram": False,
        "metrics_rows_per_level": 1,
        "metrics_cols_to_concat": 1,
        "metrics_cols_to_skip": 0,
        "help": "\"bam_name.bam.summary_metrics\" file generated from Picard."
        },
    "picard_pre_adapter_detail_metrics": {
        "picard": True,
        "picard_metrics": True,
        "picard_histogram": False,
        "metrics_rows_per_level": 192,
        "metrics_cols_to_concat": 3,
        "help": "\"bam_name.bam.wgs_metrics\" file generated from GATK CollectMultipleMetrics."
    },
    "picard_quality_by_cycle_metrics": {
        "picard": True,
        "picard_metrics": False,
        "picard_histogram": True,
        "metrics_rows_per_level": 1,
        "metrics_cols_to_concat": 1,
        "hist_col_per_level": 1,
        "help": "\"bam_name.bam.wgs_metrics\" file generated from GATK CollectMultipleMetrics."
    },
    "picard_quality_distribution_metrics": {
        "picard": True,
        "picard_metrics": False,
        "picard_histogram": True,
        "metrics_rows_per_level": 1,
        "metrics_cols_to_concat": 1,
        "hist_col_per_level": 1,
        "help": "\"bam_name.bam.wgs_metrics\" file generated from GATK CollectMultipleMetrics."
    },
    "tgen_mutation_burden": {
        "picard": True,
        "picard_metrics": True,
        "picard_histogram": False,
        "metrics_rows_per_level": 1,
        "metrics_cols_to_concat": 1,
        "help": "\"bam_name.annotate_flag.mutation_burden\" file generated from tgen_mutation_burden in a Picard style output."
    },
    "quality_yield_metrics": {
        "picard": True,
        "picard_metrics": True,
        "picard_histogram": False,
        "metrics_rows_per_level": 1,
        "metrics_cols_to_concat": 1,
        "help": "\"bam_name.bam.quality_yield_metrics\" file generated from Picard."
    },
    "samtools_idxstats": {
        "samtool": True,
        "help": "\"bam_name.bam.flagstats.txt\" file generated from Samtools."
    },
    "samtools_markdup": {
        "samtool": True,
        "help": "\"bam_name.bam.markdup.txt\" file generated from Samtools."
    },
    "samtools_flagstats": {
        "samtool": True,
        "help": "\"bam_name.bam.flagstats.txt\" file generated from Samtools."
    },
    "verifybamid": {
        "verifybamid": True,
        "help": "\"bam_name.bam.verifybamid2.selfSM\" file generated from verifyBamID"
    },
    "bt_cell_counts": {
        "btCellCounts": True,
        "help": "\"bam_name.bam.BTcell.loci.counts.txt\" file generated from rna_getBTcellLociCounts task in the "
                "phoenix pipeline "
    },
    "sex_check": {
        "sexCheck": True,
        "help": "\"bam_name.bam.sexCheck.txt\" file generated from freebayes_sex_check task in the phoenix pipeline."
    },
    "oxog_metrics_summary": {
        "oxogMetricsSummary": True,
        "help": "\"bam_name.bam.oxog_metrics_summary.tsv\" file generated from gatk_convertsequencingarrtifacttooxog "
                "task in the phoenix pipeline. "
    },
    "snpsniffer_summary": {
        "snpSnifferSummary": True,
        "help": "\"project_name.SnpSniffer_Mismatch_Summary.tsv\" file generated from snpsniffer_summary task in the "
                "phoenix pipeline. "
    },
    "cellranger_metrics": {
        "cellrangerMetrics": True,
        "help": "\"bam_name.metrics_summary.csv\" file generated from the cellranger task in the "
                "phoenix pipeline. "
    },
    "cellranger_vdj_metrics": {
        "cellrangerVDJMetrics": True,
        "help": "\"bam_name.metrics_summary.csv\" file generated from the cellranger task in the "
                "phoenix pipeline. "
    },
    "starsolo_metrics": {
        "starsoloMetrics": True,
        "help": "\"bam_name.Barcodes.stats\" file generated from the starsolo_count task in the "
                "phoenix pipeline. "
    },
    "samtools_baseQualityYield_summary": {
        "samtool": True,
        "help": "\"bam_name.samtools_baseQualityYield_summary.tsv\" file generated from the samtools_stats task in "
                "the phoenix pipeline. "
    },
    "samtools_coverage_summary": {
        "samtool": True,
        "help": "\"bam_name.samtools_coverage_summary.tsv\" file generated from the samtools_stats task in the "
                "phoenix pipeline. "
    },
    "samtools_insertSize_summary": {
        "samtool": True,
        "help": "\"bam_name.samtools_insertSize_summary.tsv\" file generated from the samtools_stats task in the "
                "phoenix pipeline. "
    },
    "samtools_summaryNumbers_summary": {
        "samtool": True,
        "help": "\"bam_name.samtools_summaryNumbers_summary.tsv\" file generated from the samtools_stats task in the "
                "phoenix pipeline. "
    },
    "samtools_markdup_summary": {
        "samtool": True,
        "help": "\"bam_name.samtools_markdup_summary.tsv\" file generated from the samtools_stats task in the "
                "phoenix pipeline. "
    }
}

# Additional picard files that may need support in future development
"""
    "bait_Bias_Detail_Metrics": {
        "picard_metrics": True,
        "help": "\"bam_name.bam.bait_Bias_Detail_Metrics\" file generated from Picard."
        },
    "bait_bias_summary_metrics": {
        "picard_metrics": True,
        "help": "\"bam_name.bam.bait_bias_summary_metrics\" file generated from Picard."
        },
    "base_distribution_by_cycle_metrics": {
        "picard_metrics": True,
        "help": "\"bam_name.bam.base_distribution_by_cycle_metrics\" file generated from Picard."
        },
    "gc_bias.detail_metrics": {
        "picard_metrics": True,
        "help": "\"bam_name.bam.detail_metrics\" file generated from Picard."
        },
    "pre_adapter_detail_metrics": {
        "picard_metrics": True,
        "help": "\"bam_name.bam.pre_adapter_detail_metrics\" file generated from Picard."
        },
    "pre_adapter_summary_metrics": {
        "picard_metrics": True,
        "help": "\"bam_name.bam.pre_adapter_summary_metrics\" file generated from Picard."
        },
    "quality_by_cycle_metrics": {
        "picard_metrics": True,
        "picard_histogram": True,
        "help": "\"bam_name.bam.quality_by_cycle_metrics\" file generated from Picard."
        },
    "quality_distribution_metrics": {
        "picard_metrics": True,
        "picard_histogram": True,
        "help": "\"bam_name.bam.quality_distribution_metrics\" file generated from Picard."
        }
"""


def description(ftype):
    """Create the epilog string for the command line help."""
    string = "NOTES: \n"
    string = string + "  -FileTypes"
    for key in ftype:
        string = string + "\n     " + key + "\t\t" + ftype[key]["help"]

    return string


def parse_arguments_and_validate():
    """Parse arguments, validate and return the args"""

    parser = argparse.ArgumentParser(
        description='Convert stats files generated by gatk or samtools into json format.',
        epilog=description(fileTypes),
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('statfile', metavar="statfile", help='the file that will be converted to json format')

    parser.add_argument('filetype', metavar='FileType', choices=fileTypes,
                        help='the file type that will be imported')

    parser.add_argument('-O', '--output', default=currentDirectory,
                        help='the output directory the json file will be saved to')

    parser.add_argument('-S', '--samplename', metavar="<arg>",
                        help='accumulated stats at the sample level')

    parser.add_argument('-L', '--libraryname', metavar="<arg>",
                        help='accumulated stats at the library level, requires sample name arguments')

    parser.add_argument('-R', '--readgroupname', metavar="<arg>",
                        help='accumulated stats at the readgroup level, requires library and sample arguments')

    parser.add_argument('-C', '--chrname', metavar="<arg>",
                        help='optional flag to rename chromosomes')

    # need to make dynamic still "2.0"
    parser.add_argument('--version', action='version', version='%(prog)s 2.0')

    # prints help message when 0 arguments are entered
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    parser_args = parser.parse_args()

    # Validate input arguments
    if parser_args.readgroupname:
        # If the readgoupname is provided then both the libraryname and samplename must be provided as well
        if not parser_args.libraryname:
            print("Error: Library name required.")
            sys.exit()
        if not parser_args.samplename:
            print("Error: Sample name required.")
            sys.exit()

    if parser_args.libraryname:
        # If the libraryname is provided then the samplename must be provided as well.
        if not parser_args.samplename:
            print("Error: Sample name required")
            sys.exit()

    return parser_args


# Used within output_file_location to determine if user has permissions to the keyed directory.
def permissions(y):
    """Exit if the user does not have permissions to write to the output directory."""
    if not os.access(y, os.W_OK):
        print("Error: User does not have permissions to write to output directory.")
        sys.exit()


# Returns the name and location the file will be saved to
# args(inputfile , outputfile)
def output_file_location(i, output):
    """Checks if the output is a file or a directory and returns the output name of the json. """

    # outputFile is an existing directory. Need to save the default filename to this directory
    if os.path.isdir(output):
        permissions(output)
        output = output + "/" + os.path.basename(i) + ".json"

    # outputFile is not an existing directory. Need to check if the base is an existing directory
    else:

        # outputFile is the name of the file with the correct directories
        if os.path.isdir(os.path.dirname(output)):
            permissions(os.path.dirname(output))

        # outputFile is the name of the file that is to be used.
        else:
            # outputFile directory does not exist
            if output.__contains__("/"):
                print("Error: output directory does not exist")
                sys.exit()

            # outputFile is only the name of the file with no directories
            else:
                permissions(currentDirectory)
                output = currentDirectory + "/" + output

    return output


# Get sample, library, and read_group name for picard files that have the columns
def get_level_name(data_list, row_index, data_head_all):
    # Set the sample, library, and read_group values if not none.
    if "SAMPLE" in data_head_all:
        samplename = data_list[row_index][data_head_all.index("SAMPLE")]

        if samplename == "" and row_index > 0:
            for i in range(row_index, -1, -1):
                samplename = data_list[i][data_head_all.index("SAMPLE")]
                if len(samplename) > 0:
                    break

    elif "SAMPLE_ALIAS" in data_head_all:
        samplename = data_list[row_index][data_head_all.index("SAMPLE_ALIAS")]
    else:
        samplename = None

    libraryname = data_list[row_index][data_head_all.index("LIBRARY")] if "LIBRARY" in data_head_all else None
    readgroupname = data_list[row_index][data_head_all.index("READ_GROUP")] if "READ_GROUP" in data_head_all else None

    return samplename, libraryname, readgroupname


# Convert file_string into list
def stats_string_to_list(file_string):
    """Convert file_string into list of rows and a list of col values"""

    # Split the contents of the file into a list using the newline character
    file_row_list = file_string.split("\n")

    # Create an empty list to append to
    file_col_list = []

    # If an item (row) in the stats_file_row_list is not empty then split on tab and append to stats_file_col_list.
    for file_row in file_row_list:
        if not file_row == '':
            row = file_row.split("\t")
            file_col_list.append(row)

    return file_row_list, file_col_list


# Create filtered string
def filter_string(string, start, end):
    """Returns a string that is contains all of the characters
    in between and including the start and end from string ."""
    m = re.compile(r'%s.*?%s' % (start, end), re.S)
    filt_string = m.search(string).group(0)

    return filt_string


# Provides location within the json structure for where the current row is to be inserted
def add_to_json(input_json, json_object_name, sample_name=None, library_name=None, readgroup_name=None):
    """Provides location within the json structure for where the current row is to be inserted."""
    # TODO comment code within add_to_json

    if not input_json:
        input_json = {
            "SAMPLES": [
                {
                    "LIBRARIES": [
                        {
                            "READGROUPS": [
                                {

                                }
                            ]
                        }
                    ]
                }
            ]
        }

    flag_sample = False
    flag_library = False
    flag_readgroup = False

    # Bam level location. return the string for the top level
    if not sample_name and not library_name and not readgroup_name:
        return json_object_name, input_json

    it1 = 0
    it2 = 0
    it3 = 0

    final1 = 0
    final2 = 0
    final3 = 0

    first_sample = False
    first_library = False
    first_readgroup = False

    for sample in input_json["SAMPLES"]:

        if 'SM' in sample:
            if sample["SM"] != sample_name:
                it1 += 1
                continue
            else:
                final1 = it1

        else:
            input_json["SAMPLES"][it1]["SM"] = sample_name
            first_sample = True

        flag_sample = True
        it2 = 0

        if not library_name:
            flag_library = True
            continue

        for library in sample["LIBRARIES"]:

            if 'LB' in library:
                if library["LB"] != library_name:
                    it2 += 1
                    continue
                else:
                    final2 = it2
            else:
                input_json["SAMPLES"][it1]["LIBRARIES"][it2]["LB"] = library_name
                first_library = True

            flag_library = True
            it3 = 0

            if not readgroup_name:
                flag_readgroup = True
                continue

            for readgroup in library["READGROUPS"]:

                if 'ID' in readgroup:
                    if readgroup["ID"] != readgroup_name:
                        it3 += 1
                        continue
                    else:
                        final3 = it3
                else:
                    input_json["SAMPLES"][it1]["LIBRARIES"][it2]["READGROUPS"][it3]['ID'] = readgroup_name
                    first_readgroup = True

                flag_readgroup = True

                if not first_readgroup:
                    it3 += 1
            if not first_library:
                it2 += 1
        if not first_sample:
            it1 += 1

    if flag_sample is False:
        temp = {"SM": sample_name, "LIBRARIES": []}
        input_json["SAMPLES"].append(temp)
        final1 = it1

    if flag_library is False and library_name:
        temp = {"LB": library_name, "READGROUPS": []}
        input_json["SAMPLES"][final1]["LIBRARIES"].append(temp)
        final2 = it2

    if flag_readgroup is False and readgroup_name:
        temp = {"ID": readgroup_name}
        input_json["SAMPLES"][final1]["LIBRARIES"][final2]["READGROUPS"].append(temp)
        final3 = it3

    if sample_name and library_name and readgroup_name:
        # readgroup level
        json_location = json_object_name + "[\"SAMPLES\"][" + str(final1) + "][\"LIBRARIES\"][" + \
                        str(final2) + "][\"READGROUPS\"][" + str(final3) + "]"
        return json_location, input_json
    elif sample_name and library_name:
        # library level
        json_location = json_object_name + "[\"SAMPLES\"][" + str(final1) + "][\"LIBRARIES\"][" + \
                        str(final2) + "]"
        return json_location, input_json
    elif sample_name:
        # sample level
        json_location = json_object_name + "[\"SAMPLES\"][" + str(final1) + "]"
        return json_location, input_json


def flat_file_add_to_output_dict(output_dict, location_in_json, data_dict):
    """Add key value pairs to the output_dict"""

    for key, value in data_dict.items():
        add_command = location_in_json + "[\"" + key + "\"] = " + value
        try:
            exec(add_command)
        except (NameError, SyntaxError):
            add_command = location_in_json + "[\"" + key + "\"] = \"" + value + "\""
            exec(add_command)

    return output_dict


# Extracts the command and date from the file header
def picard_header(file_col_list):
    """Extracts the command and date rows from a picard file that has been converted into a list of the columns."""

    # GATK picard stats files have a row in the files that includes the command used to make the file and a timestamp
    # the command was started. Both of these rows start with "# " with no other lines having this same pattern.
    # This step loops through each line looking for the pattern, removes it and returns the result.
    bam_headers = []
    for file_line in file_col_list:
        if file_line[0].startswith("# "):
            bam_headers.append(file_line[0].replace("# ", ""))

    return bam_headers


# Extracts the Metrics Class data and the BAM headers
def picard_metrics(file_string, skip_col=None):
    """
    Extract the picard specific stats and header row from the
    file_string and return the stats and header as lists.
    """

    filtered_string = filter_string(file_string, '## METRICS CLASS', '\n\n')

    picard_row_list, picard_col_list = stats_string_to_list(filtered_string)

    # Initialize the returned stat_data list
    stat_data = []

    # Initialize the returned header list
    data_head_list = []

    # Skipping over header rows that start with "##" the header and stats are extracted to separate lists and returned.
    for index, file_line in enumerate(picard_col_list):
        if file_line[0].__contains__("##"):
            pass
        elif index == 1:
            data_head_list = file_line
        else:
            stat_data.append(file_line)

    # Remove the SAMPLE, LIBRARY, and READ_GROUP columns from stat_data and data_head_list if they exists.
    stat_data_only = deepcopy(stat_data)
    data_head_list_stats_only = data_head_list[:]

    if "SAMPLE" in data_head_list:
        [stat_line.pop(data_head_list_stats_only.index("SAMPLE")) for stat_line in stat_data_only]
        data_head_list_stats_only.pop(data_head_list_stats_only.index("SAMPLE"))

    if "SAMPLE_ALIAS" in data_head_list:
        [stat_line.pop(data_head_list_stats_only.index("SAMPLE_ALIAS")) for stat_line in stat_data_only]
        data_head_list_stats_only.pop(data_head_list_stats_only.index("SAMPLE_ALIAS"))

    if "LIBRARY" in data_head_list:
        [stat_line.pop(data_head_list_stats_only.index("LIBRARY")) for stat_line in stat_data_only]
        data_head_list_stats_only.pop(data_head_list_stats_only.index("LIBRARY"))

    if "READ_GROUP" in data_head_list:
        [stat_line.pop(data_head_list_stats_only.index("READ_GROUP")) for stat_line in stat_data_only]
        data_head_list_stats_only.pop(data_head_list_stats_only.index("READ_GROUP"))

    if skip_col is not None:
        [stat_line.pop(skip_col) for stat_line in stat_data_only]
        data_head_list_stats_only.pop(skip_col)

    return stat_data, data_head_list, stat_data_only, data_head_list_stats_only


# Extracts the Histogram data from a picard file
def picard_histogram(file_string):
    """Extracts the Histogram data from a picard file."""

    filtered_string = filter_string(file_string, '## HIST', '\n\n')

    picard_hist_row_list, picard_hist_col_list = stats_string_to_list(filtered_string)

    # Initialize the returned hist_data list
    hist_data = []

    # Initialize the returned hist_header list
    hist_head = []

    # Skipping over header rows that start with "##" the header and histogram
    # stats extracted to separate lists and returned.
    for index, hist_split_line in enumerate(picard_hist_col_list):
        if hist_split_line[0].__contains__("##"):
            pass
        elif index == 1:
            hist_head_temp = hist_split_line

            # removes extra characters appended to column headers by Picard
            hist_head = []
            for col in hist_head_temp:
                col = col.replace(".fr", "")
                hist_head.append(col)
        else:
            hist_data.append(hist_split_line)

    return hist_data, hist_head


def picard_add_metrics_to_dict(data_list, data_list_stats, data_head_all,
                               data_head_stats, rows_per_level, cols_to_concat):
    """Function used when filetype argument is picards alignment summary metric (asm)."""

    # Initilize the output_dict, data_dict_temp, and row_accumulation_counter
    output_dict = {}
    data_dict_temp = {}
    row_accumulation_counter = 1

    if rows_per_level == "All":
        rows_per_level = len(data_list_stats)

    # Loops through each row within the input file excluding the SAMPLE, LIBRARY, and READ_GROUP columns
    for row_index, data_line in enumerate(data_list_stats):
        if rows_per_level != 1 or data_head_stats[0] == "CATEGORY":
            # Loops through each stat column modifying the key name by concatenating
            # the column name with the CATEGORY columns
            for key_name in data_head_stats[cols_to_concat:]:
                category = data_line[0]
                if cols_to_concat > 1:
                    for category_col in data_line[1:cols_to_concat]:
                        category = category + "_" + category_col

                data_dict_temp[key_name + "_" + category] = data_line[data_head_stats.index(key_name)]
        else:
            for key_name in data_head_stats[0:]:
                data_dict_temp[key_name] = data_line[data_head_stats.index(key_name)]

        # Set the sample, library, and read_group values if not none.
        samplename, libraryname, readgroupname = get_level_name(data_list, row_index, data_head_all)

        data_line_sample_name = args.samplename or samplename
        data_line_library_name = args.libraryname or libraryname
        data_line_read_group_name = args.readgroupname or readgroupname

        # The picard metrics files can have a different number of lines associated to one accumulation level
        # keep the number of rows_per_level and then add to the output_dict when reset.
        if row_accumulation_counter % rows_per_level == 0:
            dict_location, output_dict = add_to_json(
                output_dict,
                "output_dict",
                data_line_sample_name,
                data_line_library_name,
                data_line_read_group_name
            )

            for data_temp_line in data_dict_temp:
                # If the value of the of data_temp_line is not a float then wrap the value in double quotes so it will
                # get added to the dictionary correctly.
                try:
                    float(data_dict_temp[data_temp_line])
                except ValueError:
                    add_command = dict_location + "[\"" + data_temp_line + "\"] = " + \
                                  "\"" + data_dict_temp[data_temp_line] + "\""
                else:
                    add_command = dict_location + "[\"" + data_temp_line + "\"] = " + \
                                  data_dict_temp[data_temp_line]

                exec(add_command)

            # resets asm_temp for next 3 rows
            data_dict_temp = {}

        row_accumulation_counter += 1

    return output_dict


def picard_add_histogram_to_dict(input_dict, hist_data, hist_head, cols_per_level,
                                 sample_index, library_index, read_group_index, index):
    """Add the the output from picard_histogram to the appropriate location in the final json."""

    dict_location, input_dict = add_to_json(input_dict, "input_dict", sample_index, library_index, read_group_index)

    # Initilize the data_dict_temp and calculate the column to start and end at.
    data_dict_temp = {}
    col_start = index * cols_per_level + 1
    col_end = col_start + cols_per_level

    # Parse the keys from the hist_data from the first column
    keys = [x[0] for x in hist_data]

    # Parse the histogram data from hist_data for each of the columns in this level
    for col in range(col_start, col_end):
        try:
            values = [int(x[col]) for x in hist_data]
        except ValueError:
            values = [float(x[col]) for x in hist_data]

        data_dict_temp[hist_head[col]] = dict(zip(keys, values))

    if "HISTOGRAM" in eval(dict_location + ".keys()"):
        add_command = dict_location + "[\"HISTOGRAM\"][\"" + hist_head[0] + "\"] = data_dict_temp"
    else:
        add_command = dict_location + "[\"HISTOGRAM\"] = {\"" + hist_head[0] + "\": data_dict_temp}"

    exec(add_command)

    return input_dict


# Extracts the data from a samtool output file
def samtool_data_extract(file_row_list, samtools_file_type, sample_index, library_index, read_group_index):
    """Extract data from a samtools stats file."""

    output_dict = {}

    location_in_json, output_dict = add_to_json(
        output_dict,
        "output_dict",
        sample_index,
        library_index,
        read_group_index
    )

    if samtools_file_type == "samtools_markdup":
        # Initialize the data_list dictionary that will have key values added and returned.

        # Loop through each of the fileds less the last newline
        for file_row_list_row in file_row_list[:-1]:
            # The delimeter is ": " in the samtools markdup.txt file and is used to split each row.
            # except for the last row "ESTIMATED_LIBRARY_SIZE" which needs to be handled separately.
            if "ESTIMATED_LIBRARY_SIZE" in file_row_list_row:
                key, val = file_row_list_row.split(" ")
            else:
                key, val = file_row_list_row.split(": ")

            # The whitespace is removed from each key if it exists and replaced with an underscore.
            key = key.replace(" ", "_")

            # If the value of the of data_temp_line is not a float then wrap the value in double quotes so it will
            # get added to the dictionary correctly.
            try:
                float(val)
            except ValueError:
                # The whitespace is removed from each key if it exists and replaced with an underscore.
                add_command = location_in_json + "[\"" + key + "\"] = " + \
                              "\"" + val + "\""
            else:
                add_command = location_in_json + "[\"" + key + "\"] = " + \
                              val

            exec(add_command)

    elif samtools_file_type == "samtools_flagstats":
        # Samtools is not structured to make parsing easy. We check number of rows for 13 lines. Key names are hardcoded
        # since the file isn't easily parsable for key. Users may change key names in the data_list section.
        # support for providing custom keynames may be provided in the future.

        data_list = [
            "intotal",
            "secondary",
            "supplementary",
            "duplicates",
            "mapped",
            "pairedinsequencing",
            "read1",
            "read2",
            "properlypaired",
            "withitselfandmatemapped",
            "singletons",
            "withmatemappedtoadifferentchr",
            "withmatemappedtoadifferentchrmapQ_ge_5"
        ]

        # Store the QC-passed reads. Future support may be added to include the QC-failed reads.
        # However, Illumina bcl2fastq does not include these reads in the fastqs so for now support is not included.
        flag_stat = [x.split(" ")[0] for x in file_row_list[:-1]]

        # Validate the expected number of lines in the flagstat file
        if len(flag_stat) != 13:
            print("Error: Incorrect number of rows within flagstats file. Expected 13 lines but ", len(flag_stat),
                  "lines in the file")
            exit(1)

        data_dict = dict(zip(data_list, flag_stat))

        for key, value in data_dict.items():
            add_command = location_in_json + "[\"" + key + "\"] = " + value

            exec(add_command)

    elif samtools_file_type == "samtools_idxstats":

        for file_row_list_row in file_row_list[:-1]:
            split_list = file_row_list_row.split("\t")

            key = split_list[0]
            value = split_list[2]

            add_command = location_in_json + "[\"" + key + "\"] = " + value

            exec(add_command)

    elif samtools_file_type in ["samtools_baseQualityYield_summary",
                                "samtools_coverage_summary",
                                "samtools_insertSize_summary",
                                "samtools_summaryNumbers_summary",
                                "samtools_markdup_summary"]:

        keys = file_row_list[0].split("\t")
        values = file_row_list[1].split("\t")

        data_dict = dict(zip(keys, values))
        data_dict.pop("Sample")
        data_dict.pop("Library")
        data_dict.pop("Read_Group")
        data_dict.pop("BAM_File")

        for key, value in data_dict.items():
            add_command = location_in_json + "[\"" + key + "\"] = " + value

            exec(add_command)

    return output_dict


def verifybamid_data_extract(file_col_list, sample_index, library_index, read_group_index):
    """Extract data from a verifybamid stats file."""

    output_dict = {}

    location_in_json, output_dict = add_to_json(
        output_dict,
        "output_dict",
        sample_index,
        library_index,
        read_group_index
    )

    data_list = [
        "BAM",
        "PC1",
        "PC2",
        "alpha"
    ]

    data_dict = dict(zip(data_list, file_col_list[0]))
    del data_dict["BAM"]

    for key, value in data_dict.items():
        add_command = location_in_json + "[\"" + key + "\"] = " + value

        exec(add_command)

    return output_dict


def bt_cell_counts_data_extract(file_row_list, sample_index, library_index, read_group_index):
    """Extract data from a bt cell counts file."""

    output_dict = {}

    location_in_json, output_dict = add_to_json(
        output_dict,
        "output_dict",
        sample_index,
        library_index,
        read_group_index
    )

    loci = [x.split("\t")[3] for x in file_row_list[:-1]]
    counts = [x.split("\t")[4] for x in file_row_list[:-1]]

    data_dict = dict(zip(loci, counts))

    for key, value in data_dict.items():
        add_command = location_in_json + "[\"" + key + "\"] = " + value

        exec(add_command)

    return output_dict


def sex_check_data_extract(file_col_list, sample_index, library_index, read_group_index):
    """Extract data from a bt cell counts file."""

    output_dict = {}

    location_in_json, output_dict = add_to_json(
        output_dict,
        "output_dict",
        sample_index,
        library_index,
        read_group_index
    )

    data_dict = dict(zip(file_col_list[0], file_col_list[1]))
    del data_dict["BAM"]
    del data_dict["SAMPLE"]

    output_dict = flat_file_add_to_output_dict(output_dict, location_in_json, data_dict)

    return output_dict


def oxog_metrics_summary_data_extract(file_col_list, sample_index, library_index, read_group_index):
    """Extract data from a bt cell counts file."""

    output_dict = {}

    location_in_json, output_dict = add_to_json(
        output_dict,
        "output_dict",
        sample_index,
        library_index,
        read_group_index
    )

    data_dict = dict(zip(file_col_list[0], file_col_list[1]))

    output_dict = flat_file_add_to_output_dict(output_dict, location_in_json, data_dict)

    return output_dict


def snpsniffer_summary_data_extract(file_col_list):
    """Extract data from a bt cell counts file."""

    output_dict = {}

    header = file_col_list.pop(0)

    for row in file_col_list:
        sample_index = row[0]
        library_index = row[1]
        read_group_index = None

        location_in_json, output_dict = add_to_json(
            output_dict,
            "output_dict",
            sample_index,
            library_index,
            read_group_index
        )

        data_dict = dict(zip(header, row))
        del data_dict["SAMPLE"]
        del data_dict["LIBRARY"]

        output_dict = flat_file_add_to_output_dict(output_dict, location_in_json, data_dict)

    return output_dict


def cellranger_metrics_data_extract(dict_list_stats, sample_index, library_index, read_group_index):
    """Extract data from a cellranger metrics_summary.csv file."""

    output_dict = {}

    location_in_json, output_dict = add_to_json(
        output_dict,
        "output_dict",
        sample_index,
        library_index,
        read_group_index
    )

    data_dict = {}

    for key in dict_list_stats:
        if "%" in dict_list_stats[key]:
            value = float(dict_list_stats[key].replace("%", ""))
            value = value / 100
            value = float("%.3f" % round(value, 3))
        else:
            value = dict_list_stats[key].replace(",", "")

        field = key.replace(' ', '_').replace(',', '').replace('(', '').replace(')', '').replace('-', '_')
        data_dict[field] = value

    for key, value in data_dict.items():
        add_command = location_in_json + "[\"" + key + "\"] = " + str(value)

        exec(add_command)

    return output_dict


def starsolo_metrics_data_extract(file_row_list, sample_index, library_index, read_group_index):
    """Extract data from a bt cell counts file."""

    output_dict = {}

    location_in_json, output_dict = add_to_json(
        output_dict,
        "output_dict",
        sample_index,
        library_index,
        read_group_index
    )

    keys = [x.split()[0] for x in file_row_list[:-1]]
    values = [x.split()[1] for x in file_row_list[:-1]]

    data_dict = dict(zip(keys, values))

    for key, value in data_dict.items():
        add_command = location_in_json + "[\"" + key + "\"] = " + value

        exec(add_command)

    return output_dict


def rename_chromosomes(stats_file_string, chr_dictionary):
    for item in chr_dictionary.keys():
        stats_file_string = re.sub('^'+item, chr_dictionary[item], stats_file_string, flags=re.MULTILINE)
    return stats_file_string


if __name__ == '__main__':

    # Parse and validate arguments
    args = parse_arguments_and_validate()

    # Store the input stats file as a string
    dict_list = []

    csv_file = csv.DictReader(open(args.statfile, 'rt'))

    for line in csv_file:
        dict_list.append(line)

    with open(args.statfile, mode="rt") as file:
        stats_file_string = file.read()

    if args.chrname:
        chr_dict = {}
        with open(args.chrname, mode="rt") as f:
            for line in f:
                key, val = line.strip().split()
                chr_dict[key] = val.strip()
        stats_file_string = rename_chromosomes(stats_file_string, chr_dict)

    stats_file_row_list, stats_file_col_list = stats_string_to_list(stats_file_string)

    # Initialize the final_json
    final_json = {}

    if "samtool" in fileTypes[args.filetype]:
        final_json = samtool_data_extract(
            stats_file_row_list,
            args.filetype,
            args.samplename,
            args.libraryname,
            args.readgroupname
        )

    if "picard" in fileTypes[args.filetype]:
        # parse out the command and date started from the picard stats_file_col_list into a list
        bamHead = picard_header(stats_file_col_list)

        if args.filetype == 'picard_alignment_summary_metrics':
            if stats_file_col_list[6][0] == 'UNPAIRED':
                metrics_rows_per_level = 1
            else:
                metrics_rows_per_level = fileTypes[args.filetype]["metrics_rows_per_level"]
        else:
            metrics_rows_per_level = fileTypes[args.filetype]["metrics_rows_per_level"]

        if fileTypes[args.filetype]["picard_metrics"]:
            if "metrics_cols_to_skip" in fileTypes[args.filetype]:
                col_to_skip = fileTypes[args.filetype]["metrics_cols_to_skip"]
            else:
                col_to_skip = None
            # parse out the METRICS stat lines and columns headers from the stats_file_string
            data, data_head, data_stats_only, data_head_stats_only = picard_metrics(stats_file_string,
                                                                                    skip_col=col_to_skip)

            final_json = picard_add_metrics_to_dict(
                data,
                data_stats_only,
                data_head,
                data_head_stats_only,
                metrics_rows_per_level,
                fileTypes[args.filetype]["metrics_cols_to_concat"]
            )

        if fileTypes[args.filetype]["picard_histogram"] and fileTypes[args.filetype]["picard_metrics"]:
            # parse out the METRICS stat lines and columns headers from the stats_file_string
            data, data_head, data_stats_only, data_head_stats_only = picard_metrics(stats_file_string)

            histogram_data, histogram_data_head = picard_histogram(stats_file_string)

            index_end = int(len(data)/metrics_rows_per_level)

            for level_index in range(0, index_end):
                hist_sample, hist_library, hist_readgroup = get_level_name(data, level_index, data_head)

                final_json = picard_add_histogram_to_dict(
                    final_json,
                    histogram_data,
                    histogram_data_head,
                    fileTypes[args.filetype]["hist_col_per_level"],
                    args.samplename or hist_sample,
                    args.libraryname or hist_library,
                    args.readgroupname or hist_readgroup,
                    level_index
                )
        elif fileTypes[args.filetype]["picard_histogram"]:
            histogram_data, histogram_data_head = picard_histogram(stats_file_string)

            final_json = picard_add_histogram_to_dict(
                final_json,
                histogram_data,
                histogram_data_head,
                fileTypes[args.filetype]["hist_col_per_level"],
                args.samplename,
                args.libraryname,
                args.readgroupname,
                0
            )

        final_json["COMMAND"] = bamHead[0]
        final_json["DATE"] = bamHead[1]

    if "verifybamid" in fileTypes[args.filetype]:
        final_json = verifybamid_data_extract(
            stats_file_col_list,
            args.samplename,
            args.libraryname,
            args.readgroupname
        )

    if "btCellCounts" in fileTypes[args.filetype]:
        final_json = bt_cell_counts_data_extract(
            stats_file_row_list,
            args.samplename,
            args.libraryname,
            args.readgroupname
        )

    if "sexCheck" in fileTypes[args.filetype]:
        final_json = sex_check_data_extract(
            stats_file_col_list,
            args.samplename,
            args.libraryname,
            args.readgroupname
        )

    if "oxogMetricsSummary" in fileTypes[args.filetype]:
        final_json = oxog_metrics_summary_data_extract(
            stats_file_col_list,
            args.samplename,
            args.libraryname,
            args.readgroupname
        )

    if "snpSnifferSummary" in fileTypes[args.filetype]:
        final_json = snpsniffer_summary_data_extract(
            stats_file_col_list
        )

    if "cellrangerMetrics" in fileTypes[args.filetype]:
        final_json = cellranger_metrics_data_extract(
            dict_list[0],
            args.samplename,
            args.libraryname,
            args.readgroupname
        )

    if "cellrangerVDJMetrics" in fileTypes[args.filetype]:
        final_json = cellranger_metrics_data_extract(
            dict_list[0],
            args.samplename,
            args.libraryname,
            args.readgroupname
        )

    if "starsoloMetrics" in fileTypes[args.filetype]:
        final_json = starsolo_metrics_data_extract(
            stats_file_row_list,
            args.samplename,
            args.libraryname,
            args.readgroupname
        )

    outputFile = output_file_location(args.statfile, args.output)

    writeOutput = open(outputFile, "w+")
    writeOutput.write(json.dumps(final_json))
