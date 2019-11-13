#!/usr/bin/env python3

import os
import json
import sys
import requests
import base64
import argparse

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
    "picard_rna_metrics": {
        "prefix": "picardTools_rnaMetrics_",
        "help": "\"bam_name.bam.rna_metrics\" file generated from Picard."
        },
    "picard_wgs_metrics": {
        "prefix": "picardTools_wgsMetrics_",
        "help": "\"bam_name.bam.wgs_metrics\" file generated from Picard."
        },
    "samtools_markdup": {
        "prefix": "samStats_markdup_",
        "help": "\"bam_name.bam.markdup.txt\" file generated from Samtools."
        },
    "samtools_flagstats": {
        "prefix": "samStats_reads_",
        "help": "\"bam_name.bam.flagstats.txt\" file generated from Samtools."
        },

}

# additional picard files that may need support in future development
"""
    "bait_Bias_Detail_Metrics": {
        "prefix": "",
        "help": "\"bam_name.bam.bait_Bias_Detail_Metrics\" file generated from Picard."
        },
    "bait_bias_summary_metrics": {
        "prefix": "",
        "help": "\"bam_name.bam.bait_bias_summary_metrics\" file generated from Picard."
        },
    "base_distribution_by_cycle_metrics": {
        "prefix": "",
        "help": "\"bam_name.bam.base_distribution_by_cycle_metrics\" file generated from Picard."
        },
    "error_summary_metrics": {
        "prefix": "",
        "help": "\"bam_name.bam.error_summary_metrics\" file generated from Picard."
        },
    "gc_bias.detail_metrics": {
        "prefix": "",
        "help": "\"bam_name.bam.detail_metrics\" file generated from Picard."
        },
    "gc_bias.summary_metrics": {
        "prefix": "",
        "help": "\"bam_name.bam.summary_metrics\" file generated from Picard."
        },
    "pre_adapter_detail_metrics": {
        "prefix": "",
        "help": "\"bam_name.bam.pre_adapter_detail_metrics\" file generated from Picard."
        },
    "pre_adapter_summary_metrics": {
        "prefix": "",
        "help": "\"bam_name.bam.pre_adapter_summary_metrics\" file generated from Picard."
        },
    "quality_by_cycle_metrics": {
        "prefix": "",
        "help": "\"bam_name.bam.quality_by_cycle_metrics\" file generated from Picard."
        },
    "quality_distribution_metrics": {
        "prefix": "",
        "help": "\"bam_name.bam.quality_distribution_metrics\" file generated from Picard."
        },
    "quality_yield_metrics": {
        "prefix": "",
        "help": "\"bam_name.bam.quality_yield_metrics\" file generated from Picard."
    }
"""


def description(ftype):

    string = "NOTES: \n"
    string = string + "  -FileTypes"

    for key in ftype:
        string = string + "\n     " + key + "\t\t" + ftype[key]["help"]

    return string


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


# prints help message when no arguments are entered
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

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
responseToken = requests.request("POST", urlStart, data=payload, headers=headers, verify=CA)
parsed_responseToken = json.loads(responseToken.text)

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

    # Loops through each library in the sample and uploads each to KBase.
    counter = 0
    libraryID = ""
    prefix = fileTypes[args.filetype]['prefix']

    for library in contents["SAMPLES"][0]["LIBRARIES"]:
        patchData = {"fieldData": {}}
        if args.libraryID:
            libraryData = {"query": [{"DCL Sample ID": "=" + args.libraryID}]}
        else:
            for x in library:
                if x == "READGROUPS" or x == "HISTOGRAM":
                    pass
                elif x == "LB":
                    libraryID = contents["SAMPLES"][0]["LIBRARIES"][counter][x]
                else:
                    patchData["fieldData"][prefix + x] = contents["SAMPLES"][0]["LIBRARIES"][counter][x]

            libraryData = {"query": [{"DCL Sample ID": "=" + libraryID}]}

        # POSTS
        # Initializing validation checks
        studyData = {"query": [{"Study Name Code": "=" + args.study}]}

        json_studyData = json.dumps(studyData)
        json_libraryData = json.dumps(libraryData)

        postResponse1 = requests.request("POST", urlFind, data=json_libraryData, headers=headers, verify=CA)
        parsed_postResponse1 = json.loads(postResponse1.text)

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

        postResponse2 = requests.request("POST", urlFind2, data=json_studyData, headers=headers, verify=CA)
        parsed_postResponse2 = json.loads(postResponse2.text)

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
        json_patchData = json.dumps(patchData, indent=4)

        # The record ID used in the patch
        recordID = str([d['recordId'] for d in parsed_postResponse1['response']['data']][0])

        patchResponse = requests.request("PATCH", urlPatch + recordID, data=json_patchData, headers=headers, verify=CA)

        parsed_patchResponse = json.loads(patchResponse.text)

        if patchResponse.status_code != 200:
            print("ERROR:\n" + "Code: " + str(parsed_patchResponse["messages"][0]["code"]) +
                  "\n" + "Message: " + str(parsed_patchResponse["messages"][0]['message']))
            raise ValueError("Something went wrong with the PATCH.")

        counter += 1

finally:
    # Closes out opened session with FileMaker.
    close_headers = {
        'cache-control': "no-cache",
        'Authorization': 'Bearer %s' % token
    }

    close_payload = ""

    close_response = requests.request("DELETE",
                                      urlStart + "/" + token,
                                      data=close_payload,
                                      headers=close_headers,
                                      verify=CA)
    print(close_response)
