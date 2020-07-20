#!/usr/bin/env python3

import os
import json
import sys
import requests as requests
import base64
import argparse
from random import randint
from time import sleep

# ADD THE PREFIX FOR EACH FILE TYPE BASED ON FIELD NAME WITHIN the TGen LIMS
fileTypes = {
    "picard_alignment_summary_metrics": {
        "prefix": "picardTools_alignmentStats_",
        "help": "\"bam_name.bam.alignment_Summary_Metrics\" file generated from Picard."
        },
    "picard_hs_metrics": {
        "prefix": "picardTools_hsStats_",
        "help": "\"bam_name.bam.hs_Metrics\" file generated from Picard."
        },
    "picard_insert_size_metrics": {
        "prefix": "picardTools_insertSizeStats_",
        "help": "\"bam_name.bam.insert_size_metrics\" file generated from Picard."
        },
    "picard_quality_yield_metrics": {
        "prefix": "picardTools_qualityYieldMetrics_",
        "help": "\"bam_name.bam.quality_yield_metrics\" file generated from Picard."
        },
    "picard_gc_bias_summary_metrics": {
        "prefix": "picardTools_gc_bias_",
        "help": "\"bam_name.bam.summary_metrics\" file generated from Picard."
        },
    "picard_rna_metrics": {
        "prefix": "picardTools_rnaMetrics_",
        "help": "\"bam_name.bam.rna_metrics\" file generated from Picard."
        },
    "picard_wgs_wnzc_metrics": {
        "prefix": "picardTools_wnczMetrics_",
        "help": "\"bam_name.bam.wgs_metrics\" file generated from Picard."
        },
    "picard_error_summary_metrics": {
        "prefix": "picardTools_errorMetrics_",
        "help": "\"bam_name.bam.error_summary_metrics\" file generated from Picard."
        },
    "samtools_idxstats": {
        "prefix": "idxStatsChrCount_",
        "help": "\"bam_name.bam.flagstats.txt\" file generated from Samtools."
    },
    "samtools_markdup": {
        "prefix": "samStats_markdup_",
        "help": "\"bam_name.bam.markdup.txt\" file generated from Samtools."
        },
    "samtools_flagstats": {
        "prefix": "samStats_reads_",
        "help": "\"bam_name.bam.flagstats.txt\" file generated from Samtools."
        },
    "bt_cell_counts": {
        "prefix": "samToolsIgCount_",
        "help": "\"bam_name.bam.BTcell.loci.counts.txt\" file generated from rna_getBTcellLociCounts task in the "
                "phoenix pipeline."
    },
    "verifybamid": {
        "prefix": "verifybamid_",
        "help": "\"bam_name.bam.verifybamid2.out\" file generated from verifybamid2."
    },
    "sex_check": {
        "prefix": "sexCheck_",
        "help": "\"bam_name.bam.sexCheck.out\" file generated from freebayes_sex_check task in the phoenix pipeline."
    },
    "oxog_metrics_summary": {
        "prefix": "picardTools_oxogMetrics_",
        "help": "\"bam_name.bam.oxog_metrics_summary.tsv\" file generated from gatk_convertsequencingarrtifacttooxog "
                "task in the phoenix pipeline. "
    },
    "snpsniffer_summary": {
        "prefix": "snpsniffer_summary_",
        "help": "\"project_name.SnpSniffer_Mismatch_Summary.tsv\" file generated from snpsniffer_summary "
                "task in the phoenix pipeline. "
    },
    "cellranger_metrics": {
        "prefix": "cellranger_metrics_",
        "help": "\".cellranger_count.bam.metrics_summary.csv\" file generated from the cellranger "
                "task in the phoenix pipeline. "
    },
    "cellranger_vdj_metrics": {
        "prefix": "cellranger_metrics_",
        "help": "\".cellranger_vdj.bam.metrics_summary.csv\" file generated from the cellranger_vdj_ "
                "task in the phoenix pipeline. "
    },
    "starsolo_metrics": {
        "prefix": "starsolo_metrics_",
        "help": "\"bam_name.Barcodes.stats\" file generated from the starsolo_count task in the phoenix pipeline. "
    },
    "samtools_baseQualityYield_summary": {
        "prefix": "samStats_baseQualityYield_",
        "help": "\"bam_name.Barcodes.stats\" file generated from the starsolo_count task in the phoenix pipeline. "
    },
    "samtools_coverage_summary": {
        "prefix": "samStats_coverage_",
        "help": "\"bam_name.Barcodes.stats\" file generated from the starsolo_count task in the phoenix pipeline. "
    },
    "samtools_insertSize_summary": {
        "prefix": "samStats_insertSize_",
        "help": "\"bam_name.Barcodes.stats\" file generated from the starsolo_count task in the phoenix pipeline. "
    },
    "samtools_summaryNumbers_summary": {
        "prefix": "samStats_summaryNumbers_",
        "help": "\"bam_name.Barcodes.stats\" file generated from the starsolo_count task in the phoenix pipeline. "
    }
}

# additional picard files that may need support in future development
"""
    "base_distribution_by_cycle_metrics": {
        "prefix": "",
        "help": "\"bam_name.bam.base_distribution_by_cycle_metrics\" file generated from Picard."
        },
    "gc_bias.detail_metrics": {
        "prefix": "",
        "help": "\"bam_name.bam.detail_metrics\" file generated from Picard."
        },
    "quality_by_cycle_metrics": {
        "prefix": "",
        "help": "\"bam_name.bam.quality_by_cycle_metrics\" file generated from Picard."
        },
    "quality_distribution_metrics": {
        "prefix": "",
        "help": "\"bam_name.bam.quality_distribution_metrics\" file generated from Picard."
        }
"""

# samtools idx contigs to push into the LIMS.
# These are the only ones that currently available.
'''
contig_list = [
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "chrX",
    "chrY",
    "chrM",
    "chrUn_GL000220v1",
    "chr22_KI270733v1_random",
    "chrUn_KN707925v1_decoy",
    "chr14_KI270726v1_random",
    "chrUn_JTFH01000119v1_decoy",
    "chrUn_JTFH01000204v1_decoy",
    "chrUn_JTFH01000678v1_decoy",
    "chrUn_JTFH01000497v1_decoy",
    "chrUn_KN707927v1_decoy",
    "chrUn_JTFH01000562v1_decoy",
    "chrUn_JTFH01000659v1_decoy",
    "chrUn_JTFH01001134v1_decoy",
    "chrUn_JTFH01000888v1_decoy",
    "chrUn_JTFH01000695v1_decoy",
    "chrUn_JTFH01000697v1_decoy"
]
'''


def description(ftype):

    string = "NOTES: \n"
    string = string + "  -FileTypes"

    for key in ftype:
        string = string + "\n     " + key + "\t\t" + ftype[key]["help"]

    return string


def server_request(request_type, url, data, header):
    """Request server response and try up to 20 time with a delay if unsuccessful"""

    for attempt in range(20):
        response = requests.request(request_type,
                                    url,
                                    data=data,
                                    headers=header,
                                    verify=CA)

        print("This is request attempt: " + str(attempt + 1))
        print(response)
        print(response.content)
        print(response.text)

        if response.status_code == 502:
            if attempt == 19:
                print("ERROR:\n")
                raise ValueError("The LIMS server has returned a https 502 error 20 consecutive times.")
            else:
                sleep(randint(4, 10))
                continue
        else:
            break

    parsed_response = json.loads(response.text)

    return response, parsed_response


# Command line arguments.
parser = argparse.ArgumentParser(description="Module used to push sample QC stats into the LIMS. Used in conjunction"
                                             " with samStats2json output files", epilog=description(fileTypes),
                                 formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('jsonfile', metavar="StatFile", help='path to the json that will be pushed to the LIMS')

parser.add_argument('filetype', metavar='FileType', choices=fileTypes,
                    help='the file type that will be imported')

parser.add_argument('isilonpath', metavar='IsilonPath',
                    help='the path where the results of the study from the pipeline were stored')

parser.add_argument('project', metavar='Project',  help='patient ID')

parser.add_argument('study', metavar='Study',  help='name of the study (project)')

parser.add_argument('-L', "--libraryID", metavar='LibraryID', help='name of the libraryID')

parser.add_argument('-c', '--contigList', metavar='STR', default='None', type=str,
                    help='comma separated list of contigs to push into the LIMS . (=\'None\')')

# prints help message when no arguments are entered
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

contig_list = args.contigList.split(',')

# Required cert for the TGen LIMS
CA = "/packages/python/3.6.0/share/TGenCA.cer"

# Set REST api https paths used on the private TGen network
urlStart = "https://tgen-fmp3.ad.tgen.org/fmi/data/v1/databases/KBase_V2/sessions"
urlFind = "https://tgen-fmp3.ad.tgen.org/fmi/data/v1/databases/KBase_V2/layouts/SamplesSequencing_REST/_find"
urlFind2 = "https://tgen-fmp3.ad.tgen.org/fmi/data/v1/databases/KBase_V2/layouts/DCLStudies_REST/_find"
urlPatch = "https://tgen-fmp3.ad.tgen.org/fmi/data/v1/databases/KBase_V2/layouts/SamplesSequencing_REST/records/"

userName = None

# Get token from the TGen LIMS to continue doing work
with open(os.path.expanduser("~/credentials"), mode='r') as p:
    for line in p:
        cred = line.split(":")
        if cred[0] == "jetstream":
            userName = line.rstrip()

b64Val = base64.b64encode(userName.encode('UTF-8')).decode("ascii")

payload = "{}"
headers = {
    'Content-Type': "application/json",
    'cache-control': "no-cache",
    'Authorization': 'Basic %s' % b64Val
    }

# Get the token from the response.
responseToken, parsed_responseToken = server_request("POST", urlStart, payload, headers)

# This step will produce a KeyError if a token was not provided in the response.
# Possible reasons:
#   - credentials are invalid
#   - server is down
#   - urlStart is incorrect
token = parsed_responseToken['response']['token']

# Set the header needed for requests
headers = {
            'Content-Type': "application/json",
            'cache-control': "no-cache",
            'Authorization': 'Bearer %s' % token
        }

# Open json file and store it as contents
try:
    with open(args.jsonfile, mode='rt') as file:
        contents = json.loads(file.read())

    # Loop through each sample
    sample_counter = 0

    for sample in contents["SAMPLES"]:

        # Loops through each library in the sample and uploads each to KBase.
        counter = 0
        libraryID = ""
        prefix = fileTypes[args.filetype]['prefix']

        for library in contents["SAMPLES"][sample_counter]["LIBRARIES"]:
            patchData = {"fieldData": {}}

            if args.libraryID:
                libraryData = {"query": [{"DCL Sample ID": "=" + args.libraryID}]}
            else:
                for x in library:
                    if x == "READGROUPS" or x == "HISTOGRAM":
                        pass
                    elif x == "LB":
                        libraryID = contents["SAMPLES"][sample_counter]["LIBRARIES"][counter][x]
                    else:
                        # Special filtering for idx stats
                        if prefix == "idxStatsChrCount_" and x in contig_list:
                            patchData["fieldData"][prefix + x] = \
                                contents["SAMPLES"][sample_counter]["LIBRARIES"][counter][x]
                        elif prefix != "idxStatsChrCount_":
                            patchData["fieldData"][prefix + x] = \
                                contents["SAMPLES"][sample_counter]["LIBRARIES"][counter][x]

                libraryData = {"query": [{"DCL Sample ID": "=" + libraryID}]}

            # POSTS
            # Initializing validation checks
            studyData = {"query": [{"Study Name Code": "=" + args.study}]}

            json_studyData = json.dumps(studyData)
            json_libraryData = json.dumps(libraryData)

            postResponse1, parsed_postResponse1 = server_request("POST", urlFind, json_libraryData, headers)

            if postResponse1.status_code != 200:
                print("ERROR:\n" + "Code: " + str(parsed_postResponse1["messages"][0]["code"]) +
                      "\n" + "Message: " + str(parsed_postResponse1["messages"][0]['message']))
                raise ValueError("Something went wrong finding the study name code and project name")

            db_study = parsed_postResponse1['response']['data'][0]['fieldData']['DCL_Studies::Study Name Code']
            db_project = parsed_postResponse1['response']['data'][0]['fieldData']['DCL Patient ID']

            if not args.study == db_study:
                print("UPLOAD FAILED: Study name code does not match.")
                print("Entered value: " + args.study)
                print("Value from KBase: " + db_study)
                raise Exception('The provided Study does not match the Study Name Code in the LIMS.')

            if not args.project == db_project:
                print("UPLOAD FAILED: Project name does not match.")
                print("Entered value: " + args.project)
                print("Value from KBase: " + db_project)
                raise Exception('The provided Project does not match the DCL Patient ID in the LIMS.')

            postResponse2, parsed_postResponse2 = server_request("POST", urlFind2, json_studyData, headers)

            if postResponse2.status_code != 200:
                print("ERROR:\n" + "Code: " + str(parsed_postResponse2["messages"][0]["code"]) +
                      "\n" + "Message: " + str(parsed_postResponse2["messages"][0]['message']))
                raise ValueError("Something went wrong finding the Isilon path")

            db_isilonpath = parsed_postResponse2['response']['data'][0]['fieldData']['Pipeline_FinalStoragePath']

            if not args.isilonpath == db_isilonpath:
                print("UPLOAD FAILED: Isilon path does not match. ")
                print("Entered value: " + args.isilonpath)
                print("Value from KBase: " + db_isilonpath)
                raise Exception('The provided IsilonPath does not match the Pipeline_FinalStoragePath in the LIMS.')

            # PATCH
            # Converts data to be patched into json format
            json_patchData = json.dumps(patchData)

            # The record ID used in the patch
            recordID = str([d['recordId'] for d in parsed_postResponse1['response']['data']][0])

            for i in range(20):
                patchResponse, parsed_patchResponse = server_request("PATCH",
                                                                     urlPatch + recordID,
                                                                     json_patchData,
                                                                     headers)

                if patchResponse.status_code != 200:
                    if i == 19:
                        print("ERROR:\n" + "Code: " + str(parsed_patchResponse["messages"][0]["code"]) +
                              "\n" + "Message: " + str(parsed_patchResponse["messages"][0]['message']))
                        raise ValueError("The library record has been in use by another user for 20 consecutive "
                                         "attempts.")
                    elif parsed_patchResponse["messages"][0]["code"] == '301':
                        sleep(randint(4, 10))
                        continue
                    else:
                        print("ERROR:\n" + "Code: " + str(parsed_patchResponse["messages"][0]["code"]) +
                              "\n" + "Message: " + str(parsed_patchResponse["messages"][0]['message']))
                        raise ValueError("Something went wrong with the PATCH.")
                else:
                    break

            counter += 1

        sample_counter += 1

finally:
    # Closes out opened session with FileMaker.
    close_headers = {
        'cache-control': "no-cache",
        'Authorization': 'Bearer %s' % token
    }

    close_payload = ""
    print("Closing connection...")

    close_response, parsed_close_response = server_request("DELETE",
                                                           urlStart + "/" + token,
                                                           close_payload,
                                                           close_headers)

    if close_response.status_code == 200:
        print("Connection Closed: " + str(close_response))
