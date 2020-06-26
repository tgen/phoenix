#!/usr/bin/env/python3


# NOTES:
#  /// \brief Reads with MAPQ below this value are filtered out during SV locus generation and many subsequent
#   /// candidate creation and scoring steps
#   ///
#   unsigned minMapq = 15;
# os << "##FORMAT=<ID=PR,Number=.,Type=Integer,Description=\"Spanning paired-read support for the ref and alt alleles in the order listed\">\n";
# os << "##FORMAT=<ID=SR,Number=.,Type=Integer,Description=\"Split reads for the ref and alt alleles in the order listed, for reads where P(allele|read)>0.999\">\n";

import os
from sys import argv, exit  # Used to bring in the feature argv, variables or arguments
from multiprocessing import Process, Queue, cpu_count
import subprocess
import platform
import pysam
import getopt
import re
import logging

# LOGGING SETUP
# create a logger
__logfile__ = os.path.join(os.getcwd(), "log_prepare_sv_vcf.txt")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
# create a generic formatter
formatter = logging.Formatter('%(levelname)s %(asctime)-15s %(module)s %(lineno)d\t %(message)s')
# create file_handlers with different settings from generic handling
# file_handler = logging.FileHandler(__logfile__)
# file_handler.setLevel(logging.INFO)
# file_handler.setFormatter(formatter)

# create stream handler for printing to standard output
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)

# add the handler to the logger ; Allow to manage differently and customize some logs
# logger.addHandler(file_handler)
logger.addHandler(stream_handler)


class SearchDiscordantReads:
    """
    all information to be able to capture the discordant reads for the current variant record
    myvcf.header.samples.header.samples.header.samples[0]
    myvcf.header.samples.header.samples.header.samples[1]
    """

    def __init__(self, svtype, reference_genome, vcf_file, tumor_bam_fn, normal_bam_fn, variant,
                 slops, insert_size, sigma,
                 minmapq, tpl_reads_orientation=(-1, -1), coltumorsample=11, threads=1):
        self.svtype = svtype
        self.refgen = reference_genome
        self.vcf_file = vcf_file
        self.tumor_bam_fn = tumor_bam_fn
        self.normal_bam_fn = normal_bam_fn
        self.variant = variant  # pysam variant record
        self.mate_variant_rec = None  # pysam variant record
        self.chrom = variant.chrom
        self.pos = variant.pos
        self.slops = slops  # tuple of tuple of 2 slop values ((pos_left, pos_right),(end_left, end_right))
        self.insert_size = insert_size
        self.sigma = sigma
        self.chr2 = variant.info['CHR2']
        self.endpossv = variant.info['ENDPOSSV']
        self.svlen = variant.info['SVLEN']
        self.minmapq = minmapq
        self.tpl_reads_orientation = tpl_reads_orientation
        self.coltumorsample = coltumorsample
        self.threads = threads

    def print_coords_event_for_igv(self):
        logger.debug('slops = {}  ; read_orientation = {}'.format(str(self.slops), str(self.tpl_reads_orientation)))
        logger.debug("coordinates for igv: {}:{}  {}:{}".format(self.chrom, self.pos, self.chr2, self.endpossv))
        logger.debug(
            "coordinates for igv with slop added: {}:{}-{}  {}:{}-{}".format(self.chrom, self.pos - self.slops[0][0], self.pos + self.slops[0][1],
                                                                             self.chr2, self.endpossv - self.slops[1][0], self.endpossv + self.slops[1][1]))

    def print_for_debug(self):
        logger.debug("tpl_reads_orientation: {}\nminmapq: {}\nsvlen: {}\nthreads: {}".format(self.tpl_reads_orientation, self.minmapq, self.svlen, self.threads))
        self.print_coords_event_for_igv()

    def capture_mate_record_to_current_variant(self):
        if self.svtype == "BND":
            try:
                # This handle will be used to fetch the mate breakend
                vcfhandle_to_get_mate_record = pysam.VariantFile(self.vcf_file)
                logger.debug("{}:{} ; {}".format(self.chr2, self.endpossv, self.variant.info['MATEID']))
                # need here to get the handle again to get consumed variant in the list, as the mate could have been processed already and not being in myvcf handle
                variant_rec_mate = fetch_mate_variant_record(vcfhandle_to_get_mate_record, self.chr2, self.endpossv, self.variant.info['MATEID'][0])
                logger.debug("variant_rec_mate returned is: --> " + str(variant_rec_mate))

                read_orientation = get_read_orientation_from_alt(self.variant.alts)
                mate_read_orientation = get_read_orientation_from_alt(variant_rec_mate.alts)
                # we create a tuple for the read orientation of both the current and its mate breakend ; then we pass it to the function that capture the discordant reads
                self.tpl_reads_orientation = (read_orientation, mate_read_orientation)
                self.mate_variant_rec = variant_rec_mate
            except Exception as e:
                logger.error(e)


def get_sample_column_index_from_vcf_header(variant, sample_name):
    """
    return the index of the column in variant for the given sample name
    returned column index will be from 0 to N if sample name found otherwise raise error
    :param variant: variant record from pysam Variant
    :param sample_name: sample name which is expected to be in the VCF
    :return: integer index value
    """
    if len(variant.header.samples) > 2:
        raise Exception("ERROR: Number of sample found in VCF is greater than 2; expected somatic Constitutional versus Case")
    for i in range(0, len(variant.header.samples)):
        if str(sample_name) in variant.header.samples[i]:
            logger.info("TUMOR COLUMN Value is: {}   for {}".format(str(i),sample_name))
            return i
    raise Exception("ERROR: SAMPLE NOT FOUND in VCF")


def make_dico_contig_lengths_from_bam_handler(bamh):
    """
    create a dictionary of contig with their length as values
    :param bamh: bam handle from pysam bam/cram reader
    :return: dictionary of contigs:lengths
    """
    return dict(zip(bamh.references, bamh.lengths))


def get_contig_length(dict_contig_length, contig_name):
    """
    capture the contig and their length from a dictionary
    :param dict_contig_length: dictionary of contig with their length
    :param contig_name: contig name for which we want the length
    :return: length of the contig
    """
    if contig_name in dict_contig_length.keys():
        return dict_contig_length[contig_name]
    else:
        raise KeyError("ERROR: Contig name << {} >> NOT FOUND".format(contig_name))


def is_tool_in_path(name):
    """
    is_tool check if a tool is in the PATH
    :param name: tool name or executable
    :return: Boolean
    """

    try:
        devnull = open(os.devnull)
        subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            return False
    return True


def read_aln_file(aln_file, reference_genome_file, threads=2):
    """
    read the alignment file whether it is a SAM, BAM or CRAM file and returns the bam file handle
    :param threads: use this much threads when reading CRAM file; not use when BAM
    :param aln_file:
    :param reference_genome_file:
    :return: aln read handle (bamh)
    """

    extension = os.path.splitext(aln_file)[1]
    try:
        if extension == ".cram":
            return pysam.AlignmentFile(aln_file, mode='rc', reference_filename=reference_genome_file, threads=threads)
        elif extension == ".bam":
            return pysam.AlignmentFile(aln_file, mode='rb', threads=threads)
        elif extension == ".sam":
            return pysam.AlignmentFile(aln_file, mode='r')
        else:
            logger.debug("extension found: " + extension)
            raise Exception("EXPECTED EXTENSION for ALIGNMENT FILE NOT FOUND; must be either .cram, .bam or .sam")
            # exit(2)
    except FileNotFoundError as fnf:
        logger.exception(fnf)
        exit(2)
    except Exception as e:
        logger.exception(e)
        exit(2)


def write_new_vcf(vcf, newheader, LVAR, logfile=__logfile__, output_vcf=None, compress_vcf=False):
    """
    function to write the updated records in a new vcf_file
    :param vcf: name or path name of the input vcf_file to make the output_vcf name if output_vcf not provided
    :param newheader: New Header for the new VCF file
    :param LVAR: List of Variant Records in Pysam Format
    :param output_vcf: default value is None ; Otherwise value is provided by user
    :return: None
    """

    if output_vcf is None:
        output_vcf = vcf + ".flag.vcf"
    if os.path.splitext(output_vcf)[1] == '.gz':
        output_vcf = os.path.splitext(output_vcf)[0]
    try:
        with open(output_vcf, "w") as f:
            f.write(str(newheader))
            if LVAR is not None:
                for record in LVAR:
                    f.write(str(record))
    except IOError as ioe:
        logger.exception(ioe)
        logger.info("logfile is:" + logfile)
        exit(2)
    except Exception as e:
        logger.exception(e)
        logger.info("logfile is: " + logfile)
        exit(2)

    if compress_vcf:
        output_vcf_gz = output_vcf + ".gz"
        command_bcftools_view = "bcftools view -Oz -o {} {}".format(output_vcf_gz, output_vcf)
        logger.info(command_bcftools_view)
        command_bcftools_index = 'bcftools index --force --tbi {}'.format(output_vcf_gz)
        logger.info(command_bcftools_index)
        try:
            os.system(command_bcftools_view)
            os.system(command_bcftools_index)
        except IOError as ioe:
            logger.error(ioe)
        except Exception as e:
            logger.error(e)
        logger.info("vcf file written successfully")
        logger.info("output file is: " + str(output_vcf_gz))
    else:
        logger.info("vcf file written successfully")
        logger.info("output file is: " + str(output_vcf))


def update_vcf_header(vcf, slop, insert_size, arg_cmd):
    """
    update the VCF header with new lines, new fields and new flags when needed
    Add the new flags to the header in INFO and FORMAT fields
    :param vcf: pysam vcf_file object
    :param slop: value of the slop used to capture the discordant reads
    :param insert_size: value of the insert_size used by read aligner to estimate if reads are properly paired or not
    :param arg_cmd: arguments used in the command
    :return: updated pysam vcf_file object
    """

    # -----#
    # INFO #
    # -----#
    vcf.header.info.add("TRA", ".", "Flag",
                        "Variant is a TRANSLOCATION")
    vcf.header.info.add("INV", ".", "Flag",
                        "Variant is a BND Inversion")
    vcf.header.info.add("INVDUP", ".", "Flag",
                        "Variant is a BND Inverted Duplication")
    vcf.header.info.add("DUPTANDEM", ".", "Flag",
                        "Variant is a BND DUP-TANDEM")
    vcf.header.info.add("CHR2", "1", "String",
                        "Chromosome name of the associated breakend")
    vcf.header.info.add("ENDPOSSV", "1", "Integer",
                        "POSition of the associated breakend in CHR2")
    vcf.header.info.add("UNTESTED", "0", "Flag",
                        "indicate if we tried to recapture the discordant reads; apply to any INS, or to DEL with SVLEN within specific range")
    # --------#
    # FORMAT  #
    # --------#
    vcf.header.formats.add("RUT1", "1", "Float", "for POS region, Ratio of Unique Start positions over the Total of number of Start Positions (aka total number of reads)")
    vcf.header.formats.add("RUT2", "1", "Float",
                           "for ENDPOSSV region, Ratio of Unique Start positions over the Total of number of Start Positions (aka total number of reads)")
    vcf.header.formats.add("RDISTDISC1", "1", "Integer",
                           "Distribution size of Discordant pair at CHR:POS ; value 0 means no discordant reads in [POS+SLOP] interval")
    vcf.header.formats.add("RDISTDISC2", "1", "Integer",
                           "Distribution size of Discordant pair at CHR2:ENDPOSSV ; value 0 means no discordant reads in [POS+SLOP] interval")
    vcf.header.formats.add("RCDIS1", "1", "Integer",
                           "Number of Recaptured Discordant Pairs in Left Region (POS) with which we calculated the range distribution RDISTDISC1")
    vcf.header.formats.add("RCDIS2", "1", "Integer",
                           "Number of Recaptured Discordant Pairs in Right Region (ENDPOSSV) with which we calculated the range distribution RDISTDISC2")
    vcf.header.formats.add("DRNOISE1", "6", "Integer",
                           "Number of Recaptured Discordant Pairs not associated with current event; we may consider this as read noise within and around current left breakend; "
                           " idx1=DEL, idx2=DUP, idx3=INS, idx4=INV, idx5=TRA, idx6=TOTAL_READ_COUNTED_AS_NOISE")
    vcf.header.formats.add("DRNOISE2", "6", "Integer",
                           "Same as DRNOISE1 but for ENDPOSSV region")
    vcf.header.formats.add("EVALRCNOISE", "2", "Float",
                           "This variable represent a way of evaluating the noise related to the total number of reads counted for the current event. Is is calculated as follow: "
                           "(PR+SR)/(PR+SR+DRNOISE1) for first value and (RCDIS1/(RCDIS1+DRNOISE1)) for the second value")
    # ADDED LINES for COMMENTS
    vcf.header.add_line("##INSERTSIZE=" + str(insert_size) + ". Insert size value used to define whether we may have discordant reads associated to the variant")
    vcf.header.add_line("##SLOP=" + str(slop) + ". If CIPOS or/and CIEND exist in record, CIPOS is added to slop for POS, CIEND is added to slop for ENDPOSSV")
    if arg_cmd is not None:
        str_arg_cmd = ' '.join([str(x) for x in arg_cmd])
        vcf.header.add_line("##CommandPrepareSVvcf='{}'".format(str(str_arg_cmd)))
    return vcf


def capture_noise_around_pos(res, lrnoise, aln_rec, svlen, insert_size, sigma, factor_of_sigma=1):
    """
    check if the read pairs are considered as noise within the processed region of the current processed breakend aka variant record
    :param res: 4 code number for the current type of read orientation
    :param lrnoise: ordered list representing the count of reads considered as noise
    :param aln_rec:
    :param svlen: value only use in the logger; will be deprecated in later version
    :param insert_size: mean insert size
    :param sigma: sigma can be a sigma value or the standard deviation
    :param factor_of_sigma: if sigma value represents a standard deviation, this factor must be a value of 1
    :return:
    """
    try:
        a = aln_rec
        if (res == '1100' or res == '1001') and (a.reference_name != a.next_reference_name):
            logger.debug("TRANSLOCATION READ: " + str(str(aln_rec).split("\t")[0:8]))
            lrnoise[4] += 1  # TRA
            lrnoise[5] += 1  # TOT
        elif (res == '1100' or res == '1001') and (a.reference_name == a.next_reference_name):
            logger.debug("DUPLICATION READ: " + str(str(aln_rec).split("\t")[0:8]))
            lrnoise[1] += 1  # DUP
            lrnoise[5] += 1  # TOT
        elif res == '1110' or res == '1011' or res == '0100' or res == '0001':
            if a.reference_name != a.next_reference_name:
                logger.debug("TRANSLOCATION READ: " + str(str(aln_rec).split("\t")[0:8]))
                lrnoise[4] += 1  # TRA
            else:
                logger.debug("INVERSION READ " + str(str(aln_rec).split("\t")[0:8]))
                lrnoise[3] += 1  # INV
            lrnoise[5] += 1  # TOT
        elif res == '0110' or res == '0011':
            logger.debug("svlen = " + str(svlen) + " and INDEL isize = " + str(a.isize))
            if a.reference_name != a.next_reference_name:
                logger.debug("TRANSLOCATION READ: " + str(str(aln_rec).split("\t")[0:8]))
                lrnoise[4] += 1  # TRA
                lrnoise[5] += 1  # TOT
            elif a.isize < (insert_size - (factor_of_sigma * sigma)) or a.isize > insert_size + (factor_of_sigma * sigma):
                lrnoise[0] += 1  # DEL
                logger.debug("DELETION READ: " + str(str(aln_rec).split("\t")[0:8]))
                lrnoise[5] += 1  # TOT
            elif a.isize < (insert_size - (factor_of_sigma * sigma)) or a.isize > insert_size + (factor_of_sigma * sigma):
                lrnoise[2] += 1  # INS
                logger.debug("INSERTION READ: " + str(str(aln_rec).split("\t")[0:8]))
                lrnoise[5] += 1  # TOT
        else:
            logger.error("res_code: {}".format(str(res)))
            logger.error("No code assigned to the res object: Should Not have happened; Check with Authors")
    except Exception as e:
        logger.exception(e)
        exit(2)
    return lrnoise


def info_on_captured_discordant_reads(bamh, slop, chrom, pos):
    """
    FUNCTION NOT USE ; Was for testing only ; info_on_captured_discordant_reads
    :param bamh: bam handler
    :param slop: slop value to be added to pos
    :param chrom: contig name
    :param pos: current position of the event
    """

    items = bamh.fetch(contig=chrom, start=max(pos - slop, 1), end=pos + slop)
    # get the number of reads captured in that region
    tot_num_reads_in_region = sum(1 for x in items)
    logger.debug(tot_num_reads_in_region)
    logger.debug(sum(1 for x in items))


def is_read_in_expected_orientation(is_rev, expected_read_orientation):
    """
    check if the read has the expected orientation for the current event
    :param is_rev:  Boolean that tell us if the read is reverse (R)
    :param expected_read_orientation: Expected read orientation for the svtype and position we are dealing with
    :return: Boolean
    """
    logger.debug("function_name: {} ; is_rev: {}  and expected_read_orientation: {}".format("is_read_in_expected_orientation", str(is_rev), str(expected_read_orientation)))
    if is_rev and expected_read_orientation == "-":
        return True
    elif not is_rev and expected_read_orientation == "+":
        return True
    return False


def is_mate_read_in_expected_orientation(is_mate_rev, expected_read_orientation):
    """
        check if the mate read has the expected orientation for the current event
        :param is_mate_rev: Boolean that tell us if the mate read is reverse (R)
        :param expected_read_orientation: Expected read orientation for the svtype and position we are dealing with
        :return: Boolean
        """
    logger.debug("is_rev: {}  and expected_read_orientation: {}".format(str(is_mate_rev), str(expected_read_orientation)))
    if is_mate_rev and expected_read_orientation == "-":
        return True
    elif not is_mate_rev and expected_read_orientation == "+":
        return True
    return False


def extract_elements_from_bnd_alt_column(alt):
    """
    Extract read orientation and the mate chromosome for the current BND event
    we return a tuple of string where the first part is the chromosome and the second part is the read orientation
    let's have a code for read orientations
    :param alt: alt column as is from pysam method that captured it
    :return: components of the alt column as a list
    """

    logger.debug(str(alt))
    # return a list ; that list must always be of length 3 because in bnd_alt the number of square brackets must always be two
    lcomponents = list(re.split('\[|\]', alt[0]))
    logger.debug("lcomponents = " + str(type(lcomponents)))
    logger.debug("lcomponents = " + str(lcomponents))
    return lcomponents


def fetch_mate_variant_record(vcfhandle, chr_mate, pos_mate, mateid, count=0, slop=50):
    """
    We fetch the MateID variant Record for the breakend being process
    :param vcfhandle:
    :param chr_mate:
    :param pos_mate:
    :param mateid: must be a string and not a tuple
    :param count:  Normally the mate_record is found really quickly after first try because the mate_record is at expected pos_mate position; in some case not so
    we have to expand the search (hence the slop value defined below); we limit the expansion to three tries, after that we search the whole contig;
    It is defined here to be used recurrently with the method
    :param slop: this slop differs in purpose from the user given slop; the point here is to find the mate_record as fast as possible
    within the vcf record so this slop here is just a region size to lookup for the vcf record
    :return: a Unique VariantRecord (it must be one and one only) otherwise None or raise Error
    """

    res_fetch = vcfhandle.fetch(contig=str(chr_mate), start=(int(pos_mate) - slop) - 1, end=int(pos_mate) + slop)
    total_items_found = sum(1 for v in res_fetch)
    logger.debug("search region:  {}:{}-{}".format(str(chr_mate), str((int(pos_mate) - slop) - 1), str(int(pos_mate) + slop)))
    logger.debug("res_fetch ==> " + str(res_fetch))
    logger.debug("total_items_found ==> " + str(total_items_found))
    if count < 3:
        res_fetch = vcfhandle.fetch(contig=str(chr_mate), start=(int(pos_mate) - slop) - 1, end=int(pos_mate) + slop)
    else:
        res_fetch = vcfhandle.fetch(contig=str(chr_mate))
    try:
        if total_items_found >= 1:
            rec_found = None
            ## we check that the mate id is present in the search result
            for rec in res_fetch:
                logger.debug("mate rec captured by res_fetch  ==> " + str(rec))
                logger.debug(str(rec.chrom) + ":" + str(rec.pos))
                if 'MATEID' not in rec.info.keys():
                    # if we increase the slop, we might found records that are not BND breakends and therefore no MATEID is present in the INFO field
                    continue
                logger.debug("mateid we want ot find in captured/fetched records: " + str(mateid))
                # NOTE: rec.info['MATEID'] returns a tuple such as: ('MantaBND:2:254200:254201:0:0:0:1',)
                if str(rec.id) == str(mateid):
                    logger.debug("yeah mate id found ... returning rec --> " + str(rec))
                    rec_found = rec
                    break
            if rec_found is None:
                count += 1
                logger.debug("rec is still none")
                logger.info("broadening the search by increasing the slop around pos_mate b/c the value in pos_mate might not be equivalent to the value in pos_alt: loop_" + str(count))
                return fetch_mate_variant_record(vcfhandle, chr_mate, pos_mate, mateid, count=count, slop=1000 + slop)
            else:
                return rec_found
        else:
            count += 1
            logger.info("broadening the search by increasing the slop around pos_mate b/c the value in pos_mate might not be equivalent to the value in pos_alt: loop_" + str(count))
            return fetch_mate_variant_record(vcfhandle, chr_mate, pos_mate, mateid, count=count, slop=1000 + slop)
    except Exception as e:
        logger.error(e)
        logger.error("ERROR: MATE NOT FOUND; Check your VCF input to see if the id << " + str(mateid) + " >> exists.")
        exit(2)


def get_read_orientation_from_alt(alts):
    """
    capture the read orientation of the mate using the ALT column
    :param alts: pysam's variant alt field
    :return: read orientation [ possible values either + or - ]
    """
    alt_info = extract_elements_from_bnd_alt_column(alts)
    logger.debug("alt info = " + str(alt_info))
    # read orientation based on the location of the nucleotides string present in the ALT column of a BND
    # if string present on the left,  expected reads for current event are FORWARD (aka, F, aka, +) ; this means alt_info[0]!=""
    # if string present on the right,  expected reads for current event are REVERSE (aka, R, aka, -) ; this means alt_info[2]!=""
    if alt_info[0] != "":
        return "+"
    elif alt_info[2] != "":
        return "-"
    else:
        logger.error("ERROR: alt_info variable does not have the expected information to determine the read orientation of the current event")
        logger.error("ERROR: alt_info is: " + str(alt_info) + " and the current alt is: " + str(alts))
        raise Exception("ERROR: with alt_info variable in function get_read_orientation_from_alt")


def dec_to_bin(x):
    """
    convert decimal value into a binary value
    :param x: decimal value
    :return: binary value
    """
    return int(bin(x)[2:])


def build_res_code(is_rev, is_first, is_mate_rev, is_mate_first, isize):
    """
    build a 4 binary digit code value to represent the specific orientation of the reads in the pair
    :param is_rev: boolean value which tells if the read is on the reverse strand or not
    :param is_first: boolean value which tells if the read is on the first in pair or not
    :param is_mate_rev: boolean value which tells if the mate read is on the reverse strand or not
    :param is_mate_first: boolean value which tells if the mate read is on the first in pair or not
    :param isize: insert size captured from the bam alignment record
    :return: 4-digit-code value
    """

    res = ""
    if isize < 0:  # Damn I do not remember why exactly I had to check if isize<0 and reassigned the values as the opposite; I have an idea :-)
        temp_is_rev = is_mate_rev
        temp_is_first = is_mate_first
        is_mate_rev = is_rev
        is_mate_first = is_first
        is_rev = temp_is_rev
        is_first = temp_is_first
    if is_rev == 0 and is_first == 1:  # F1
        res = str(is_rev) + str(is_first)
        if is_mate_rev == 0 and is_mate_first == 0:  # F2
            res = res + str(is_mate_rev) + str(is_mate_first)  # F1F2
        elif is_mate_rev == 1 and is_mate_first == 0:  # R2
            res = res + str(is_mate_rev) + str(is_mate_first)  # F1R2
    elif is_rev == 1 and is_first == 1:  # R1
        res = str(is_rev) + str(is_first)
        if is_mate_rev == 0 and is_mate_first == 0:  # F2
            res = res + str(is_mate_rev) + str(is_mate_first)  # R1F2
        elif is_mate_rev and is_mate_first == 0:  # R2
            res = str(is_rev) + str(is_first) + str(is_mate_rev) + str(is_mate_first)  # R1R2
    elif is_rev == 0 and is_first == 0:  # F2
        res = str(is_rev) + str(is_first)
        if is_mate_rev == 0 and is_mate_first == 1:  # F1
            res = res + str(is_mate_rev) + str(is_mate_first)  # F2F1
        elif is_mate_rev == 1 and is_mate_first == 1:  # R1
            res = res + str(is_mate_rev) + str(is_mate_first)  # F2R1
    elif is_rev == 1 and is_first == 0:  # R2
        res = str(is_rev) + str(is_first)
        if is_mate_rev == 0 and is_mate_first == 1:  # F1
            res = res + str(is_mate_rev) + str(is_mate_first)  # R2F1
        elif is_mate_rev == 1 and is_mate_first == 1:  # R1
            res = res + str(is_mate_rev) + str(is_mate_first)  # R2R1
    else:
        logger.error("ERROR: No Code Build ; Missing Implementation; Ask The Author. Not good :-( ")
        exit(2)
    return res


def check_if_first_in_pair(flag_integer):
    """
    Check if the flag value embeds the information of the first in pair for the current process read
    :param flag_integer: value of the alignment flag found in column two of a bam alignment record
    :return: 1 or 0, aka, true or false
    """
    b = "{0:011}".format(dec_to_bin(flag_integer))
    logger.debug(str(flag_integer) + " --- " + str(b) + " -- " + str(b[-7]))
    if str(b[-7]) == "1":
        return 1
    return 0


def parse_cigartuple_to_sum_matches(cigartuple):
    """
    get the sum of matching based in current read from the cigar information provided as a tuple of data
    See pysam documentation to get the dictionary of matching values (https://readthedocs.org/projects/pysam/downloads/pdf/latest/)
    :param cigartuple: tuple of tuples from pysam methods (cigartuples)
    :return: sum of bases matching the reference genome template
    """

    summ = 0
    for t in cigartuple:
        if t[0] == 0:
            summ += t[1]
    return summ


def is_read_noise(svtype, lrnoise, aln_rec, vcf_svlen, insert_size, sigma, factor_sigma=3):
    """
    check if the read that does not fit the current event, i.e. that is still a discordant read
    but not involved in current event is considered as noise within the current region of either POS or END breakend
    Important Note: We assume here that the read is NOT properly paired
    :param factor_sigma: in case a factor has to be apply to sigma I added the variable; if value is 1, no impact on results
    :param svtype: not used yet but later
    :param lrnoise:
    :param aln_rec:
    :param vcf_svlen:
    :param insert_size:
    :param sigma:
    :return: updated lrnoise (list of reads considered as noise)
    """

    # LRNOISE = [0, 0, 0, 0, 0, 0]
    # each index represent a SVTYPE assuming we always have the same number of SVTYPE; otherwise will need update
    # let's consider the SVTYPE alphabetically here. idx1=DEL, idx2=DUP, idx3=INS, idx4=INV, idx5=TRA, idx6=TOTAL_READ_COUNTED_AS_NOISE)
    # default_indexes = {'1100':4, '1001':4, '1110':1, '1011':1, '0100':4, '0001':4, '011':1, '001':1 } for list LRNOISE
    # indexes 3 and 5 will depend respectively upon the svlen for INS and upon the fact that the reads are on to two different chromosomes for TRA
    a = aln_rec
    svlen = vcf_svlen[0]
    logger.debug("current event in process is {} ; ".format(str(svtype)))
    # init variables
    is_rev = 0  # is current read reverse read? 0 false, 1 true
    is_mrev = 0  # is mate read of current read reverse? 0 false, 1 true
    is_first_in_pair = check_if_first_in_pair(a.flag)  # 1 if first in pair, 0 if not
    logger.debug("FIRST in PAIR: " + str(is_first_in_pair))
    try:
        if a.is_reverse:
            is_rev = 1
        if a.mate_is_reverse:
            is_mrev = 1
        # let's create a code here to assign indexes
        res = build_res_code(is_rev, is_first_in_pair, is_mrev, int(not is_first_in_pair), a.isize)

        logger.debug("RES_CODE: " + res + "  ; flag = " + str(a.flag) + " --> readname : " + str(a.query_name))
        logger.debug("svlen = " + str(svlen) + " and i-size = " + str(a.isize))
        if (res == '1100' or res == '1001') and (a.reference_name != a.next_reference_name):
            logger.debug("TRANSLOCATION READ: " + str(str(aln_rec).split("\t")[0:8]))
            lrnoise[4] += 1  # TRA
            lrnoise[5] += 1  # TOT
        elif (res == '1100' or res == '1001') and (a.reference_name == a.next_reference_name):
            logger.debug("DUPLICATION READ: " + str(str(aln_rec).split("\t")[0:8]))
            lrnoise[1] += 1  # DUP
            lrnoise[5] += 1  # TOT
        elif res == '1110' or res == '1011' or res == '0100' or res == '0001':
            if a.reference_name != a.next_reference_name:
                logger.debug("TRANSLOCATION READ: " + str(str(aln_rec).split("\t")[0:8]))
                lrnoise[4] += 1  # TRA
            else:
                logger.debug("INVERSION READ " + str(str(aln_rec).split("\t")[0:8]))
                lrnoise[3] += 1  # INV
            lrnoise[5] += 1  # TOT
        elif res == '0110' or res == '0011':
            logger.debug("svlen = " + str(svlen) + " and in-del a.isize = " + str(a.isize))
            if a.reference_name != a.next_reference_name:
                logger.debug("TRANSLOCATION READ: " + str(str(aln_rec).split("\t")[0:8]))
                lrnoise[4] += 1  # TRA
                lrnoise[5] += 1  # TOT
            elif abs(a.isize) < (insert_size - (factor_sigma * sigma)) or abs(a.isize) > insert_size + (factor_sigma * sigma):
                logger.debug("ISIZE = {} ; insert_size={}|sigma={} ; res code is:  {} ; and aln_rec is {} ".format(str(a.isize), str(insert_size), str(sigma), str(res), str(a)))
                lrnoise[0] += 1  # DEL
                logger.debug("INDEL DELETION/INSERTION READ: " + str(str(aln_rec).split("\t")[0:8]))
                lrnoise[5] += 1  # TOT
        else:
            logger.error("No code assigned to the res object: Should Not have happened; Check with Authors")

    except Exception as e:
        logger.exception(e)
        exit(2)
    return lrnoise


def getDiscordantReadsDistribution_parallel_version(sdr, svtype, refgen, bam, chr, pos, slops, chr2, pos2, svlen, minmapq, tpl_reads_orientation, insert_size, sigma, my_queue):
    """
    see documentation of the methods "getDiscordantReadsDistribution2" that describe most of the arguments
    here we describe the additional argument only
    :param my_queue: queue object to be processed
    :return: tuple of data
    """
    res = getDiscordantReadsDistribution2(sdr, svtype, refgen, bam, chr, pos, slops, chr2, pos2, svlen, minmapq, tpl_reads_orientation, insert_size, sigma)
    my_queue.put(res)


def getDiscordantReadsDistribution(sdr):
    """
    same method as the method getDiscordantReadsDistribution2 but with only one argument embedding all the args of interest in the class object SearchDiscordantReads (sdr)
    :param sdr: object of class search_discordant_reads
    :return: tuple of values of interest (see getDiscordantReadsDistribution2)
    """
    return getDiscordantReadsDistribution2(sdr, sdr.svtype, sdr.refgen, sdr.bam,
                                          sdr.chromosome, sdr.pos, sdr.slops,
                                          sdr.chr2, sdr.pos2, sdr.svlen,
                                          sdr.minmapq, sdr.tpl_reads_orientation,
                                          sdr.insert_size, sdr.sigma)


def getDiscordantReadsDistribution2(sdr, svtype, ref_gen, bam, chrom, pos, slops, chr2, pos2, svlen, minmapq, tpl_reads_orientation, insert_size, sigma):
    """
    Capture the Discordant Pairs within given region of pos+/-slop and return some stats
    :param sdr: class object sdr [future use]
    :param insert_size: insert size
    :param sigma: standard deviation for insert size
    :param ref_gen: reference genome used to align reads
    :param bam: tumor bam file
    :param chrom: chromosome for current variant
    :param pos: position for current variant
    :param slops: tuple of slops (4 values) ; TODO make it an object from class slops
    :param chr2: mate chromosome of the current breakend
    :param pos2: mate position of the current breakend
    :param svlen: length of Mutation called by tool
    :param minmapq: value of the minimum mapping quality required for the read to be taken into account (if not the whole read pair is discarded)
    :param tpl_reads_orientation: orientation of the read given by the current variant (mostly fro BND as it is given in the ALT column)
    :return: stats about the spread to reads in the region; return currently 4 variables:
    1) Range(max(pos_start)-min(post_start),
    2) total number of pos_start, aka total number of discordant pairs for current region or break,
    3) RUT aka Ratio of Unique Start Position divided by the total number of start positions (or total number of reads)
    4) the read NOISE "seen" within the region of interest
    """

    logger.debug("##" * 50)
    logger.debug("SVTYPE is {} ; and tpl_reads_orientation_is: {} ".format(svtype, str(tpl_reads_orientation)))
    logger.debug("--" * 50)
    bamh = read_aln_file(bam, ref_gen)
    # getting the length of the chromosomes # TODO: move the creation of the dictionary to not be created each time even though it does not take lot of time
    dico_lengths_contigs = make_dico_contig_lengths_from_bam_handler(bamh)  # dict(zip(bamh.references, bamh.lengths))
    chr1len = dico_lengths_contigs[chrom]
    chr2len = dico_lengths_contigs[chr2]
    # min(pos + slops[0][1], chr1len)) allows to check if the slop does not give a value beyond the length of the chromosomes
    items = bamh.fetch(contig=chrom, start=max(pos - slops[0][0], 1), end=min(pos + slops[0][1], chr1len))
    logger.debug("Total number of reads FETCH within region: " + str(sum(1 for x in items)))
    items = bamh.fetch(contig=chrom, start=max(pos - slops[0][0], 1), end=min(pos + slops[0][1], chr1len))

    # why do we need the pos2 here? well because we want to keep the pairs that are within both expected regions Left and Right of the current variant type
    pos2lowend = min(pos2 - slops[1][0], chr2len)
    pos2highend = min(pos2 + slops[1][1], chr2len)

    logger.debug("BEGIN_FETCH_" * 30)
    logger.debug("bamh.fetch(contig=chrom, start=max(pos - slop, 1), end=pos + slop) -->   {}:{}-{}".format(chrom, max(pos - slops[0][0], 1), pos + slops[0][1]))
    logger.debug("mate_location: chr2:pos2lowend-pos2highend  -->  {}:{}-{}".format(chr2, pos2lowend, pos2highend))
    logger.debug(str(chrom) + ":" + str(pos) + ' <--> ' + str(chr2) + ":" + str(pos2))
    logger.debug(str(chrom) + " <--> " + str(pos - slops[0][0]) + " <--> " + str(pos) + " <--> " + str(pos + slops[0][1]))
    logger.debug(str(chr2) + " <--> " + str(pos2lowend) + " <--> " + str(pos2) + " <--> " + str(pos2highend))
    logger.debug("END_FETCH_" * 30)

    # init list of coordinates
    CROI = []  # Count Read Of Interest for Event of Interest, aka read not properly paired and in line with the SV type
    lrnoise = [0, 0, 0, 0, 0, 0]
    # each index represent a SVTYPE assuming we always have the same number of SVTYPE; otherwise will need update
    # let's consider the SVTYPE alphabetically here.
    # idx0=DEL, idx1=DUP, idx2=INS, idx3=INV, idx4=TRA, idx5=TOTAL_READ_COUNTED_AS_NOISE)
    count_properly_paired = 0
    count_all_discordant_reads_found_in_pos_region = 0
    # ------------------------------#
    # Loop over the captured reads
    # ------------------------------#
    for x in items:  # x is a read line or aka an alignment record
        logger.debug("Processing aln_rec: {} ::::: chrom:pos == {}:{}".format(str(str(x).split("\t")[0:8]), chrom, pos))
        mq_mate = get_mate_mq_value(x)
        mq = x.mapping_quality
        # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        # here the heart of the read selection for discordant pairs:
        # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if not x.is_proper_pair:
            if not x.is_unmapped \
                    and not x.mate_is_unmapped \
                    and not x.is_qcfail \
                    and not x.is_duplicate \
                    and not x.is_secondary \
                    and is_read_in_expected_orientation(x.is_reverse, tpl_reads_orientation[0]) \
                    and mq >= minmapq and mq_mate >= minmapq:
                # and not x.is_supplementary \
                logger.debug("%%%%%%%%%%" * 10)
                logger.debug("discordant pair found: " + str(str(x).split("\t")[0:8]))
                logger.debug("is it reverse read: {} \tExpected_reads_orientation: {}".format(str(x.is_reverse), str(tpl_reads_orientation)))
                logger.debug("bamh.fetch(contig=chrom, start=max(pos - slop, 1), end=pos + slop) [slops {} ]-->   {}:{}-{} ".format(slops, chrom, max(pos - slops[0][0], 1), pos + slops[0][1]))
                logger.debug("mate_location: chr2:pos2lowend-pos2highend  -->  {}:{}-{}".format(chr2, pos2lowend, pos2highend))
                logger.debug("%%%%%%%%%%" * 10)
                logger.debug(str(pos2lowend) + " <= " + str(x.next_reference_start) + " <= " + str(pos2highend))
                logger.debug(pos2lowend <= x.next_reference_start <= pos2highend)
                count_all_discordant_reads_found_in_pos_region += 1

                # another important section: check if the mate is also in the expected range for linked-break (CHR2:ENDPOSSV side)
                # and the read is in the expected orientation
                if pos2lowend <= x.next_reference_start <= pos2highend and is_mate_read_in_expected_orientation(x.mate_is_reverse, tpl_reads_orientation[1]):

                    logger.debug("MATE_DISCORDANT_READ PAIR matching OUR CURRENT CHR2:ENDPOSSV +/- slop region")
                    logger.debug("VALIDATED as discordant read from current EVENT : " + str(str(x).split("\t")[0:8]))
                    if "S" in x.cigarstring:
                        CROI.append(x.reference_start + parse_cigartuple_to_sum_matches(x.cigartuples))  # added x.tlen (template length) to take care of split read [ + x.tlen]
                        logger.debug("CROI content after appending: {} ; Added value:{} from read {} ; |||||  x.reference_start  {} <--> x.next_reference_start = {}".format(
                            str(len(CROI)), str(x.reference_start+x.tlen), str(x), str(x.reference_start), str(x.next_reference_start)))
                    else:
                        CROI.append(x.reference_start)
                        logger.debug("CROI content after appending: {} ; Added value:{} from read {} ; |||||  x.reference_start {} <--> x.next_reference_start = {}".format(
                            str(len(CROI)), str(x.reference_start), str(x), str(x.reference_start), str(x.next_reference_start)))
                else:
                    lrnoise = is_read_noise(svtype, lrnoise, x, svlen, insert_size, sigma)
                    logger.debug("NOISE: The Mate is NOT in any of expected regions [ POS +/- slop __ or__  ENDPOSSV +/- slop ] ...")
                    logger.debug("DISCORDANT BUT MATE NOT in EXPECTED REGION so NOISE" + str(str(x).split("\t")[0:8]))
                    # we know that the mate of the current READ is not in the expected region defined by ENDPOSSV+/-slop
                    # this mean the read belongs to another event; It may be a read belonging to a valid event ...
                    # or just a random read location which we refer as Noise ;
                    # with this implementation we cannot discriminate whether the read belong to another valid event or not # maybe a TODO
                    # How to implement:
                    # we need to count the number of Discordant Reads that are NOT related to the current event
                    # to do so we use the fact that some type of EVENT (SVTYPE such as BND (TRA or INV), DEL, DUP, INS) have different read orientation
                    # intraTRA (aka transposon) or DUP: RL (igv) or R1F2 or R2F1
                    # INV: LL/RR (igv) or F1F2/F2F1/R1R2/R2R1
                    # DEL or INS: LR(igv) or F1R2/F2R1
                    # we therefore can add to a list the reads that we consider as noise within current region of interest ; we know that the read already belongs to a discordant pair

        elif x.is_proper_pair:
        # It can be read from any type of event BUT TRA (TRA are always Discordant)
        # DEL, DUP, INS, and INV can be represented by properly paired reads.
        # DUP and INV are easy to capture properly paired reads (DUP: (-,+) or R1F2 or R2F1) ; INV: R1R2, R2R1, F1F2 or F2F1)
        # DEL and INS cannot be captured as noise when properly paired because they look like to normal reads (read orientation [ (+,-) or F1R2 or F2R1 ] ) --> So we skip them
        # this means we can have here any type of reads (from correctly paired to discordant reads but reads that are secondary alignment, supp aln, MQ of 0, or MateQ of 0, or duplicate)
            if not x.is_unmapped \
                    and not x.mate_is_unmapped \
                    and not x.is_qcfail \
                    and not x.is_duplicate \
                    and not x.is_secondary \
                    and mq >= minmapq and mq_mate >= minmapq:
                # and not x.is_supplementary \
                stranded_noise = True
                if stranded_noise \
                        and (tpl_reads_orientation == ('+', '+') or tpl_reads_orientation == ('-', '-'))\
                        and is_read_in_expected_orientation(x.is_reverse, tpl_reads_orientation[0])\
                        and is_mate_read_in_expected_orientation(x.is_reverse, tpl_reads_orientation[1]):
                    logger.debug("WWWWWWW"*20)
                    logger.debug("%%%%%%%%%%" * 10)
                    logger.debug("discordant pair found: " + str(str(x).split("\t")[0:8]))
                    logger.debug("is it reverse read: {} \tExpected_reads_orientation: {}".format(str(x.is_reverse), str(tpl_reads_orientation)))
                    logger.debug("bamh.fetch(contig=chrom, start=max(pos - slop, 1), end=pos + slop) [slops {} ]-->   {}:{}-{} ".format(slops, chrom, max(pos - slops[0][0], 1), pos + slops[0][1]))
                    logger.debug("mate_location: chr2:pos2lowend-pos2highend  -->  {}:{}-{}".format(chr2, pos2lowend, pos2highend))
                    logger.debug("%%%%%%%%%%" * 10)
                    lrnoise = is_read_noise(svtype, lrnoise, x, svlen, insert_size, sigma)
                # elif not stranded_noise \
                #         and not is_read_in_expected_orientation(x.is_reverse, tpl_reads_orientation[0])\
                #         and tpl_reads_orientation != ('+', '-')\
                #         and (x.reference_name != x.next_reference_name or x.tlen > insert_size+sigma):
                #     # Check we have the pair INSERT SIZE > than Expected
                #     lrnoise = is_read_noise(svtype, lrnoise, x, svlen)
                count_properly_paired += 1

    if not CROI:  # i.e no discordant pairs found
        logger.debug("No Discordant Pairs")
        # ------------------------------#
        # RDISTDISC:RCDIS:RUT:DRNOISE
        # ------------------------------#
        # This mean there was no discordant reads in that region ; we can return 0,0,0,lrnoise or 0,0,0,lrnoise
        return 0, 0, 1, lrnoise
    else:
        logger.debug("Final length of CROI: {}".format(str(len(CROI))))
        # print number of discordant reads within regions pos-slop and pos+slop
        # max(CR) - min(CR) # this is the difference between the left-most-hand and the right-most-hand starting positions
        # if less than 100 for Delly we use to exclude the variant; this information is capture for both ends of thebreakPoints
        # and added back to the VCF record by adding new flags to the FORMAT field as it is sample specific information ;
        # range is max(CR)-min(CR)
        # len(CR) should never be 0 here since we took car of it in the "if not CR" above
        # RUT=len(set(CR))/len(CR)  ; # RUT stand for Ratio_UniqueStartPos_over_TotalStartPost
        # we also can capture the UNIQUE number of START position and compare it to the total number of reads found
        # if the ratio of unique_start_positions/total_number_of_reads is not 1, this means some reads start at the same position
        # if the ratio is close to zero (0; well it can never reach 0 since at least one unique must be), or tend to
        # be closer to zero that 1, this means that most of the reads are piling on ; 0.5 mean fifty-fifty
        # so we also can return other statistics that show if the reads are piling on in order to filter the variants better
        # split read is a conceptual name for the situation, that your read is broken in 2 parts. And to say it is a split read you have a bit wise flag of 0x800 marking supplementary alignment.
        # It's a naming thing (another name for this is a chimeric read). It links with the fact, that 0x800 flag was added later than built in ability for aligners to find split-reads,
        # and some software which didn't have implemented this supplementary read flag give a flag of 0x100, marking that the read wasn't linear, but not saying was it true supplementary
        # (true split read) or not. So just a naming in this case
        try:
            # ------------------------------#
            # RDISTDISC:RCDIS:RUT:DRNOISE
            # ------------------------------#
            return max(CROI) - min(CROI), len(CROI), len(set(CROI)) / len(CROI), lrnoise
        except ZeroDivisionError as zde:
            logger.exception(zde)
        except Exception as e:
            logger.exception(e)


def get_mate_mq_value(aln_rec):
    """
    check if MQ exists in aln_rec or assign MQ flag and return value
    :param aln_rec:
    :return: mq_mate value
    """
    if aln_rec.has_tag('MQ'):  # some alignment records do not have the MQ value in the mate when pair aligned; BWA-known-issue
        mq_mate = aln_rec.get_tag('MQ')
    else:
        mq_mate = 60  # we assume that if the mate quality is not present, it is still a good one.
    return mq_mate


def process_insertion(sdr_obj):
    """
    Processing Insertions (INS); The reads involved in that type of event is not possible to discriminate from the normal reads
    Their orientation is the same F2R1 or F1R2 and is below the insert size
    :param sdr_obj: object of class searchDiscordantReads
    :return: updated variant record
    """
    sdr_obj.print_coords_event_for_igv()
    sdr_obj.variant.info['UNTESTED'] = True
    return get_discordant_info_for_each_sample_and_breaks(sdr_obj)
    # return sdr_obj.variant


def process_deletion(sdr_obj):
    """
    Processing deletions (DEL); The reads involved in that type of event is not possible to discriminate from the normal reads if the svlen of teh event is less than or within insert_size
    Their orientation is the same F2R1 or F1R2
    If the svlen is greater than insert_size, therefore, the event is described by NOT properly paired reads and they can be counted
    :param sdr_obj: object of class searchDiscordantReads
    :return: updated variant record
    """

    sdr_obj.print_coords_event_for_igv()
    if sdr_obj.variant.info['SVLEN'][0] <= sdr_obj.insert_size + sdr_obj.sigma:
        sdr_obj.variant.info['UNTESTED'] = True
    return get_discordant_info_for_each_sample_and_breaks(sdr_obj)
    # return sdr_obj.variant


def process_duplication(sdr_obj):
    """
    Processing duplications (DUP) where the expected type of reads is always the same whatever the svlen value, if svlen > insert_size+sigma
    The read orientation is expected to be: R1F2 or R2F1 ; R for POS breakend reads and F for END breakend reads
    :param sdr_obj: object of class searchDiscordantReads
    :return: updated variant record
    """
    sdr_obj.print_coords_event_for_igv()
    return get_discordant_info_for_each_sample_and_breaks(sdr_obj)
    # return sdr_obj.variant


def process_inversion(sdr_obj):
    """
    Processing inversions (INV). Inversion may be represented by two or four breakends, but 4 is now the norm
    The ALT column provided the information about the mate read orientation and from that information we can extrapolate the read orientation
    for the current pos as the inversion are always represented by the following read orientation combinations:
    R1R2, R2R1, F1F2 and F2F1;
    So if the read mate is reverse, the current read pos is reverse as well; If mate is Forward, the current read is Forward too
    :param sdr_obj: object of class searchDiscordantReads
    :return: updated variant record
    """
    sdr_obj.print_coords_event_for_igv()
    return get_discordant_info_for_each_sample_and_breaks(sdr_obj)
    # return sdr_obj.variant


def process_translocation(sdr_obj):
    """
    Processing inversions (INV). Inversion may be represented by two or four breakends, but 4 is now the norm
    The ALT column provided the information about the mate read orientation and from that information we can extrapolate the read orientation
    The read orientation is expected to be as the ones found with Duplication
    R1F2 or R2F1
    :param sdr_obj: object of class searchDiscordantReads
    :return: updated variant record
    """
    sdr_obj.print_coords_event_for_igv()
    return get_discordant_info_for_each_sample_and_breaks(sdr_obj)
    # return sdr_obj.variant


def pre_process_variant_to_get_chr2_endpossv_and_read_orientations(vcf_file, variant):
    """
    Capture the chr2 and endpossv from the information of the current variant
    1. check of its type (svtype)
    2. grab the read orientation from ALT if svtype is BND to make sure that the discordant read captured will be from the correct orientation
    3. return chr2, endpossv and read_orientation
    :param vcf_file: full pathname to the input vcf_file
    :param variant: current object variant to be processed VariantFile from pysam (for variant in myvcf_handle)
    :return: chr2, endpossv value (cannot be None) and the read orientations for the event (both read and its mate read)
    """

    CHR2, ENDPOSSV, tpl_reads_orientation, variant_rec_mate = None, None, None, None
    if "BND" in variant.info['SVTYPE']:
        # we deal here with either TRA or INV
        try:
            # NOTE: BND SVTYPEs are the only events containing square brackets within ALT column
            # We can therefore extract the necessary information all at once
            # 1) we need to extract the CHR and the POS of the MATE
            # 2) we need to know if the READ orientation
            list_alts_components = extract_elements_from_bnd_alt_column(variant.alts)
            chr_mate = list_alts_components[1].split(":")[0]
            pos_mate = int(list_alts_components[1].split(":")[1])
            mateid = variant.info['MATEID'][0]

            logger.debug("variants.alts: " + str(variant.alts))
            logger.debug("list_alts_components  --->  " + str(list_alts_components))
            logger.debug("chr_mate  --->  " + str(chr_mate))
            logger.debug("pos_mate  --->  " + str(pos_mate))

            # This handle will be used to fetch the mate breakend
            vcfhandle_to_get_mate_record = pysam.VariantFile(vcf_file)
            # need here to get the handle again to get consumed variant in the list, as the mate could have been processed already and not being in myvcf handle
            variant_rec_mate = fetch_mate_variant_record(vcfhandle_to_get_mate_record, chr_mate, pos_mate, mateid)
            logger.debug("variant_rec_mate returned is: --> " + str(variant_rec_mate))

            read_orientation = get_read_orientation_from_alt(variant.alts)
            mate_read_orientation = get_read_orientation_from_alt(variant_rec_mate.alts)
            # we create a tuple for the read orientation of both the current and its mate breakend ; then we pass it to the function that capture the discordant reads
            tpl_reads_orientation = (read_orientation, mate_read_orientation)

            logger.debug("variant_rec_mate " + str(variant_rec_mate))
            logger.debug("READ ORIENTATION " + str(read_orientation))
            logger.debug("MATE READ ORIENTATION " + str(mate_read_orientation))
            CHR2 = chr_mate
            ENDPOSSV = int(pos_mate)
            logger.debug("ASSOCIATED VARIANT in CHR2:ENDPOSSV --> " + str(CHR2) + ":" + str(ENDPOSSV))

            if CHR2 is None or ENDPOSSV is None:
                logger.debug(str(variant))
                raise Exception("ERROR: No Pattern Found in BND; printing the variant record to "
                                "check why the Regular expression did not catch the correct pattern")
        except TypeError as ve:
            logger.exception(ve)
        except Exception as e:
            logger.exception(str(variant))
            logger.exception("ERROR: NO Pattern Found in BND ; printing the variant record to " + str(variant) +
                             " ; check why the Regular expression did not catch the correct pattern")
            logger.exception(str(e))
            exit(2)
    else:
        # we deal with either DEL, DUP or INS here (at least with Manta version 1.6)
        # there is no information about the expected read orientation mentioned in the ALT column
        # but we know the type of expected orientation by technology Illumina
        try:
            CHR2 = variant.chrom
            ENDPOSSV = variant.stop
            if "DUP" in variant.info['SVTYPE']:
                tpl_reads_orientation = ('-', '+')  # actually, this order is only when we look at it from the left break point standpoint;
                # we can state that because we only have one record line in the vcf_file for that type of DUP event
            else:
                tpl_reads_orientation = ('+', '-')  # actually, this order is only when we look at it from the left break point standpoint
                # we can state that because we only have one record line in the vcf_file for that type of INDEL event
        except TypeError as ve:
            logger.exception("ERROR TYPE")
            logger.exception(ve)
        except Exception as e:
            logger.exception(str(e))
    return CHR2, ENDPOSSV, tpl_reads_orientation, variant_rec_mate


def set_slops_from_read_orientation(slop, tpl_reads_orientation, cipos, ciend, buffer=100):
    """
    set the slop boundaries according to the svtype and the reads orientation fro that svtype
    :param slop: minimum slop value
    :param tpl_reads_orientation: tuple of read orientation with string character + or - values
    :param cipos: CIPOS value from variant record
    :param ciend: CIEND value from variant record
    :param buffer: slop like value to which the author defined to default as 100 (value not based on any data); this value is added to the slop value no matter what;
    :return: tuple with 4 slop values
    """

    # :param svtype: svtype is important to know what side of the break the reads is (based on igv read orientation)
    # we first init the 4 slops variables
    slop_pos_left, slop_pos_right, slop_end_left, slop_end_right = 0, 0, 0, 0
    if tpl_reads_orientation[0] == "+" and tpl_reads_orientation[1] == "+":
        slop_pos_left = cipos + slop
        slop_pos_right = cipos
        slop_end_left = ciend + slop
        slop_end_right = ciend
    elif tpl_reads_orientation[0] == "-" and tpl_reads_orientation[1] == "-":
        slop_pos_left = cipos
        slop_pos_right = cipos + slop
        slop_end_left = ciend
        slop_end_right = ciend + slop
    elif tpl_reads_orientation[0] == "-" and tpl_reads_orientation[1] == "+":
        slop_pos_left = cipos
        slop_pos_right = cipos + slop
        slop_end_left = ciend + slop
        slop_end_right = ciend
    elif tpl_reads_orientation[0] == "+" and tpl_reads_orientation[1] == "-":
        slop_pos_left = cipos + slop
        slop_pos_right = cipos
        slop_end_left = ciend
        slop_end_right = ciend + slop
    logger.debug("slop_pos_left+buffer = {}, slop_pos_right+buffer = {}, slop_end_left+buffer = {}, slop_end_right+buffer = {}".format(
        slop_pos_left+buffer, slop_pos_right+buffer, slop_end_left+buffer, slop_end_right+buffer))
    # slops is a tuple of tuples. tuple[0] is for left breakend (aka POS), tuple[1] is for right breakend (aka ENDPOSSV)
    slops = ((slop_pos_left+buffer, slop_pos_right+buffer), (slop_end_left+buffer, slop_end_right+buffer))
    return slops


def reevaluate_slop_values(variant, mate_variant_rec, slop, tpl_reads_orientation):
    """
    Recalculate the left and right slop around the breakend based on variant type and read orientation
    :param variant: pysam variant record
    :param mate_variant_rec:  pysam variant record for the mate of current variant
    :param slop: user-given slop value or default value
    :param tpl_reads_orientation: read orientation expected for the current variant
    :return: tuple of 4 slops
    """

    ci_pos, ci_end = 0, 0  # assume as if the variant is PRECISE when CIPOS or CIEND do not exist
    logger.debug("tpl_reads_orientation = {}".format(tpl_reads_orientation))
    if 'CIPOS' in variant.info.keys() \
            and variant.info['CIPOS'] != "." \
            and isinstance(variant.info['CIPOS'], tuple) \
            and isinstance(variant.info['CIPOS'][1], int):
        ci_pos = int(variant.info['CIPOS'][1])
    logger.debug("ci_pos value: ".format(str(ci_pos)))

    if 'CIEND' in variant.info.keys() \
            and variant.info['CIEND'] != "." \
            and isinstance(variant.info['CIEND'], tuple) \
            and isinstance(variant.info['CIEND'][1], int):
        ci_end = int(variant.info['CIEND'][1])
    logger.debug("ci_end value: ".format(str(ci_end)))

    if "BND" in variant.info['SVTYPE']:
        if variant.chrom != variant.info['CHR2']:
            logger.debug("This is a ------ TRA ---------")
            variant.info['TRA'] = True
        else:
            # NOTE: BND is a catch-all for a generic breakpoint, so you can't assume all to be some trans-locations (according to LUMPY's authors)
            # NOTE: LUMPY never represents breakpoints with DEL or DUP orientations as BNDs. Only one sided inversions and inter-chromosomal events
            # NOTE: BND are mostly INV in Manta
            # The VCF standard describes two types of SV notations. One is by SV types, i.e. insertions, deletions, inversions, translocations, etc.
            # The other is by breakend notations, often labelled with SVTYPE=BND. To describe a SV with breakend notations, each SV has two positions, each captured by one breakend
            # (except for inversions, which have 4 separate records). Each breakend includes a genomic locus, as well as a half interval extending out to the partner breakend.
            # In VCF BND notations, the ALT field encodes directional information of the partner breakend.
            # 1 200 . N N[5:500[ partner breakend immediately after chr1:200, starting from chr5:500 and extending rightwards
            # 1 200 . N ]5:500]N partner breakend immediately before chr1:200, extending from the left and ending at chr5:500
            # 1 200 . N [5:500[N partner breakend immediately before chr1:200, starting from chr5:500 and extending rightwards
            # 1 200 . N N]5:500] partner breakend immediately after chr1:200, extending from the left and ending at chr5:500
            logger.debug("BND____tpl_reads_orientation = {}".format(tpl_reads_orientation))
            if tpl_reads_orientation == ("+", "-"):
                variant.info['INVDUP'] = True
            elif tpl_reads_orientation == ("+", "+") or tpl_reads_orientation == ("-", "-"):
                variant.info['INV'] = True
            elif tpl_reads_orientation == ("-", "+"):
                variant.info['DUPTANDEM'] = True
            else:
                logger.error("Should Have Never Reached Here; At least one of the four type of read orientation should be present")
        if 'CIPOS' in mate_variant_rec.info.keys() \
                and mate_variant_rec.info['CIPOS'] != "." \
                and isinstance(mate_variant_rec.info['CIPOS'], tuple) \
                and isinstance(mate_variant_rec.info['CIPOS'][1], int):
            ci_pos = int(mate_variant_rec.info['CIPOS'][1])

        # ci_end = ci_pos  # actually the ci_end MUST be captured from the twin-vcf_file-record
        if 'CIEND' in mate_variant_rec.info.keys() \
                and mate_variant_rec.info['CIEND'] != "." \
                and isinstance(mate_variant_rec.info['CIEND'], tuple) \
                and isinstance(mate_variant_rec.info['CIEND'][1], int):
            ci_end = int(mate_variant_rec.info['CIEND'][1])
        else:
            ci_end = ci_pos
    return set_slops_from_read_orientation(slop, tpl_reads_orientation, ci_pos, ci_end)


def Reverse(tuples):
    new_tup = ()
    for k in reversed(tuples):
        new_tup = new_tup + (k,)
    return new_tup


def get_discordant_info_for_each_sample_and_breaks(sdr):
    """
    capture the discordant reads for the current processed variant and update the variant record
    with appropriate flags and fields
    :param sdr: object of class SearchDiscordantRead
    :return: updated variant record
    """
    # we need at least 4 processors (even though some job profiling showed that 3 processors is enough)
    # Now we do have all the information to recapture the discordant pairs for the current variant call
    # -------------------------------------#
    # check if Parallel is possibly here
    # -------------------------------------#
    if cpu_count() > 3 >= sdr.threads and platform.mac_ver()[0] == "":
        my_queue1 = Queue()
        my_queue2 = Queue()
        my_queue3 = Queue()
        my_queue4 = Queue()

        p1 = Process(target=getDiscordantReadsDistribution_parallel_version,
                     args=(
                         sdr, sdr.svtype, sdr.refgen, sdr.tumor_bam_fn, sdr.variant.chrom, sdr.variant.pos, sdr.slops,
                         sdr.variant.info["CHR2"], sdr.variant.info["ENDPOSSV"], sdr.variant.info['SVLEN'], sdr.minmapq,
                         sdr.tpl_reads_orientation, sdr.insert_size, sdr.sigma, my_queue1))
        p2 = Process(target=getDiscordantReadsDistribution_parallel_version,
                     args=(
                         sdr, sdr.svtype, sdr.refgen, sdr.tumor_bam_fn, sdr.variant.info["CHR2"], sdr.variant.info["ENDPOSSV"], Reverse(sdr.slops),
                         sdr.variant.chrom, sdr.variant.pos, sdr.variant.info['SVLEN'],
                         sdr.minmapq, Reverse(sdr.tpl_reads_orientation), sdr.insert_size, sdr.sigma, my_queue2))
        p3 = Process(target=getDiscordantReadsDistribution_parallel_version,
                     args=(
                         sdr, sdr.svtype, sdr.refgen, sdr.normal_bam_fn, sdr.variant.chrom, sdr.variant.pos, sdr.slops,
                         sdr.variant.info["CHR2"], sdr.variant.info["ENDPOSSV"], sdr.variant.info['SVLEN'], sdr.minmapq,
                         sdr.tpl_reads_orientation, sdr.insert_size, sdr.sigma, my_queue3))
        p4 = Process(target=getDiscordantReadsDistribution_parallel_version,
                     args=(
                         sdr, sdr.svtype, sdr.refgen, sdr.normal_bam_fn, sdr.variant.info["CHR2"], sdr.variant.info["ENDPOSSV"], Reverse(sdr.slops),
                         sdr.variant.chrom, sdr.variant.pos, sdr.variant.info['SVLEN'],
                         sdr.minmapq, Reverse(sdr.tpl_reads_orientation), sdr.insert_size, sdr.sigma, my_queue4))
        list_queues = [my_queue1, my_queue2, my_queue3, my_queue4]
        list_p = [p1, p2, p3, p4]
        for p in list_p:
            p.start()
        for p in list_p:
            p.join()

        tumdistchr = my_queue1.get()
        tumdistchr2 = my_queue2.get()
        normdistchr = my_queue3.get()
        normdistchr2 = my_queue4.get()

        for q in list_queues:
            q.close()
    else:
        # serial version
        tumdistchr = getDiscordantReadsDistribution2(sdr, sdr.svtype, sdr.refgen, sdr.tumor_bam_fn, sdr.variant.chrom, sdr.variant.pos, sdr.slops,
                                                     sdr.variant.info['CHR2'], sdr.variant.info['ENDPOSSV'],
                                                     sdr.variant.info['SVLEN'], sdr.minmapq, sdr.tpl_reads_orientation, sdr.insert_size, sdr.sigma)
        tumdistchr2 = getDiscordantReadsDistribution2(sdr, sdr.svtype, sdr.refgen, sdr.tumor_bam_fn, sdr.variant.info['CHR2'], sdr.variant.info['ENDPOSSV'], Reverse(sdr.slops),
                                                      sdr.variant.chrom, sdr.variant.pos,
                                                      sdr.variant.info['SVLEN'], sdr.minmapq, Reverse(sdr.tpl_reads_orientation), sdr.insert_size, sdr.sigma)
        normdistchr = getDiscordantReadsDistribution2(sdr, sdr.svtype, sdr.refgen, sdr.normal_bam_fn, sdr.variant.chrom, sdr.variant.pos, sdr.slops,
                                                      sdr.variant.info['CHR2'], sdr.variant.info['ENDPOSSV'],
                                                      sdr.variant.info['SVLEN'], sdr.minmapq, sdr.tpl_reads_orientation, sdr.insert_size, sdr.sigma)
        normdistchr2 = getDiscordantReadsDistribution2(sdr, sdr.svtype, sdr.refgen, sdr.normal_bam_fn, sdr.variant.info['CHR2'], sdr.variant.info['ENDPOSSV'], Reverse(sdr.slops),
                                                       sdr.variant.chrom, sdr.variant.pos,
                                                       sdr.variant.info['SVLEN'], sdr.minmapq, Reverse(sdr.tpl_reads_orientation), sdr.insert_size, sdr.sigma)

    # assign the FLAG with correct values:
    # check column tumor sample
    idxN = 0 if sdr.coltumorsample == 1 else 1
    idxT = 1 if sdr.coltumorsample == 1 else 0
    logger.debug("captured index column Tumor sample:  " + str(sdr.coltumorsample) + " - sample indexes: Normal_" + str(idxN) + " --- Tumor_" + str(idxT))

    # adding 'DRNOISE' or the count of reads consider as noise within and around the current breakend
    sdr.variant.samples[idxT]['DRNOISE1'] = tumdistchr[3]
    sdr.variant.samples[idxT]['DRNOISE2'] = tumdistchr2[3]
    sdr.variant.samples[idxN]['DRNOISE1'] = normdistchr[3]
    sdr.variant.samples[idxN]['DRNOISE2'] = normdistchr2[3]

    # adding distance btw leftmost start and rightmost start
    sdr.variant.samples[idxT]['RDISTDISC1'] = tumdistchr[0]
    sdr.variant.samples[idxT]['RDISTDISC2'] = tumdistchr2[0]
    sdr.variant.samples[idxN]['RDISTDISC1'] = normdistchr[0]
    sdr.variant.samples[idxN]['RDISTDISC2'] = normdistchr2[0]
    # Adding Count of Discordant Reads
    sdr.variant.samples[idxT]['RCDIS1'] = tumdistchr[1]
    sdr. variant.samples[idxT]['RCDIS2'] = tumdistchr2[1]
    sdr.variant.samples[idxN]['RCDIS1'] = normdistchr[1]
    sdr.variant.samples[idxN]['RCDIS2'] = normdistchr2[1]
    # adding RUT (Ratio Unique start position over Total number of reads)
    sdr.variant.samples[idxT]['RUT1'] = tumdistchr[2]
    sdr. variant.samples[idxT]['RUT2'] = tumdistchr2[2]
    sdr.variant.samples[idxN]['RUT1'] = normdistchr[2]
    sdr.variant.samples[idxN]['RUT2'] = normdistchr2[2]
    try:
        if (sdr.variant.samples[idxT]['RCDIS1'] + sdr.variant.samples[idxT]['DRNOISE1'][5]) != 0:
            idxt_evalrcnoise_using_rcdis1 = (sdr.variant.samples[idxT]['RCDIS1'] / (sdr.variant.samples[idxT]['RCDIS1'] + sdr.variant.samples[idxT]['DRNOISE1'][5]))
        else:
            idxt_evalrcnoise_using_rcdis1 = 1
        if (sdr.variant.samples[idxN]['RCDIS1'] + sdr.variant.samples[idxN]['DRNOISE1'][5]) != 0:
            idxn_evalrcnoise_using_rcdis1 = (sdr.variant.samples[idxN]['RCDIS1'] / (sdr.variant.samples[idxN]['RCDIS1'] + sdr.variant.samples[idxN]['DRNOISE1'][5]))
        else:
            idxn_evalrcnoise_using_rcdis1 = 0

        if 'SR' in sdr.variant.info.keys():
            if (sdr.variant.samples[idxT]['PR'][1] + sdr.variant.samples[idxT]['SR'][1] + sdr.variant.samples[idxT]['DRNOISE1'][5]) !=0:
                sdr.variant.samples[idxT]['EVALRCNOISE'] = (
                (sdr.variant.samples[idxT]['PR'][1] + sdr.variant.samples[idxT]['SR'][1]) / (sdr.variant.samples[idxT]['PR'][1] + sdr.variant.samples[idxT]['SR'][1] + sdr.variant.samples[idxT]['DRNOISE1'][5]),
                idxt_evalrcnoise_using_rcdis1
                )
            else:
                sdr.variant.samples[idxT]['EVALRCNOISE'] = 1, idxt_evalrcnoise_using_rcdis1

            if (sdr.variant.samples[idxN]['PR'][1] + sdr.variant.samples[idxN]['SR'][1] + sdr.variant.samples[idxN]['DRNOISE1'][5]) != 0:
                sdr.variant.samples[idxN]['EVALRCNOISE'] = (
                (sdr.variant.samples[idxN]['PR'][1] + sdr.variant.samples[idxN]['SR'][1]) / (sdr.variant.samples[idxN]['PR'][1] + sdr.variant.samples[idxN]['SR'][1] + sdr.variant.samples[idxN]['DRNOISE1'][5]),
                idxn_evalrcnoise_using_rcdis1
                )
            else:
                sdr.variant.samples[idxN]['EVALRCNOISE'] = 0, idxn_evalrcnoise_using_rcdis1
        else:
            if sdr.variant.samples[idxT]['PR'][1] + sdr.variant.samples[idxT]['DRNOISE1'][5] != 0:
                sdr.variant.samples[idxT]['EVALRCNOISE'] = (
                    (sdr.variant.samples[idxT]['PR'][1]) / (
                                sdr.variant.samples[idxT]['PR'][1] + sdr.variant.samples[idxT]['DRNOISE1'][5]), idxt_evalrcnoise_using_rcdis1
                )
            else:
                sdr.variant.samples[idxT]['EVALRCNOISE'] = 1, idxt_evalrcnoise_using_rcdis1

            if sdr.variant.samples[idxN]['PR'][1] + sdr.variant.samples[idxN]['DRNOISE1'][5] != 0:
                sdr.variant.samples[idxN]['EVALRCNOISE'] = (
                    (sdr.variant.samples[idxN]['PR'][1]) / (
                                sdr.variant.samples[idxN]['PR'][1] + sdr.variant.samples[idxN]['DRNOISE1'][5]), idxn_evalrcnoise_using_rcdis1
                )
            else:
                sdr.variant.samples[idxN]['EVALRCNOISE'] = 0, idxn_evalrcnoise_using_rcdis1
    except ZeroDivisionError as zde:
        logger.error("DIV ZERO ERROR with variant: {}".format(str(sdr.variant)))
        logger.error(zde)
    return sdr.variant


def process_variant_record(variant, vcf_file, tumor_bam_file, normal_bam_file, ref_gen, slop, coltumorsample, minmapq, insert_size, sigma):
    """
    Process a variant record ; Capture necessary information and update the variant record with new flags and fields
    :param variant: pysam variant record
    :param vcf_file: vcf filename
    :param tumor_bam_file: tumor bam name
    :param normal_bam_file: normal bam name
    :param ref_gen: reference genome used to call variant
    :param slop: slop value integer to be added to positions of current variant in order to capture erads of interest
    :param coltumorsample: information about the column of the tumor sample in the somatic vcf
    :param minmapq: integer value of the minimum mapping quality the reas must have to consider the pair
    :param insert_size: mean insert size captured by samtools stats or any other tool that can estimate the sequenced fragment insert size; must be an integer
    :param sigma: stander deviation of the insert size (integer)
    :return: update variant record
    """

    # capture necessary information before further processing
    CHR2, ENDPOSSV, tpl_reads_orientation, mate_variant_record = pre_process_variant_to_get_chr2_endpossv_and_read_orientations(vcf_file, variant)
    # we update the INFO field with these two flags because we use them downstream
    variant.info['CHR2'] = CHR2
    variant.info['ENDPOSSV'] = ENDPOSSV

    # Dealing with SVLEN: This Field might not be present in every record so let's check here:
    if 'SVLEN' not in variant.info.keys():
        variant.info['SVLEN'] = 0

    # creation of an object to gather lots of variables and information across functions
    slops = None
    sdr = SearchDiscordantReads(variant.info['SVTYPE'], ref_gen, vcf_file, tumor_bam_file, normal_bam_file, variant,
                                slops, insert_size, sigma, minmapq, tpl_reads_orientation, coltumorsample)
    sdr.capture_mate_record_to_current_variant()
    sdr.slops = reevaluate_slop_values(variant, sdr.mate_variant_rec, slop, tpl_reads_orientation)

    # Check the type of EVENT we are dealing with and process accordingly
    if "INS" in variant.info['SVTYPE']:
        variant = process_insertion(sdr)
    elif "DEL" in variant.info['SVTYPE']:
        variant = process_deletion(sdr)
    elif "DUP" in variant.info['SVTYPE']:
        variant = process_duplication(sdr)
    elif "BND" in variant.info['SVTYPE']:
        # BND can be INV, INVDUP, DUPTANDEM or TRA by checking the CHR2
        if variant.chrom == variant.info['CHR2']:

            variant = process_inversion(sdr)
        elif variant.chrom != variant.info['CHR2']:
            variant = process_translocation(sdr)
        else:
            raise Exception("ERROR: issue with chromosome comparison; Reason: Unknown")
    else:
        raise Exception("ERROR: UNKNOWN SVTYPE; Check You Input VCF if format is compatible with current version of this tool")
    return variant


def start_processing(vcf_file, tumor_bam_file, normal_bam_file, ref_gen, slop, arg_cmd, minmapq, insert_size, sigma, tumor_sample_name_in_vcf):
    """
    looping over the variants and processing the captiure of discordant reads for each breakend listed in the vcf
    a breakend is represented by one record i.e. one line when dealing with well 4.3-spec formatted vcf
    :param vcf_file: vcf filename
    :param tumor_bam_file: tumor bam filename
    :param normal_bam_file: normal bam filename
    :param ref_gen: reference genome
    :param slop: slop value (default or set by user)
    :param arg_cmd: command used to run current script
    :param minmapq: minimum mapping quality required to take reads within read pair into account
    :param insert_size: insert_size used by the read aligner if possible in order to capture discordant reads correctly
    :param sigma: standard deviation of the mean insert size
    :param tumor_sample_name_in_vcf: tumor sample name found in VCF header in either column 10 or 11 as we must only deal with somatic call here
    :return: newHeader, LVAR which is the list of updated variant records with new Fields and Flags
    """

    # -----------------------------------------------------------------------------------------------------------------------#
    # when dealing with Illumina sequencing
    # tpl_reads_orientation will have:
    # DEL and INS should be (+,-) for POS and END respectively so F1R2 or F2R1
    # DUP should be (-,+) so R1F2 or R2F1
    # For INV and TRA we rely on the ALT column where the information about the expected read orientation of the mate is embedded in
    # for INV, if mate is on the Reverse strand, this mean the current POS is also on the Reverse strand, so (-,-), aka R1R1 or R2R1
    # for INV, if mate is on the Forward strand, this mean the current POS is also on the Forward strand, so (-,-), aka F1F2 or F2F1
    # for TRA, the orientation should always be (-, +), so R1F2 or R2F1
    # -----------------------------------------------------------------------------------------------------------------------#
    # tumor_bam_fn = tumor_bam_file  # example: MMRF_1816_1_BM_CD138pos_T1_KHWGL_L06685.bwa.final.bam
    # normal_bam_fn = normal_bam_file  # example: MMRF_1816_1_PB_Whole_C2_KHWGL_L06684.bwa.final.bam
    # slop_ori = slop  # deprecated

    # read the input VCF file using pysam
    myvcf = pysam.VariantFile(vcf_file)  # no mode implies we let pysam guessing the correct mode to use (see help)
    tot_variants_count = sum(1 for v in myvcf)
    logger.info("total number of variants in VCF: " + str(tot_variants_count))
    # regenerate the collection object myvcf after consumption
    myvcf.reset()
    # update vcf_file header
    myvcf = update_vcf_header(myvcf, slop, insert_size, arg_cmd)

    # Check if there are no variant, therefore we stop right here ; Need to do this here before any call to myvcf's records
    if tot_variants_count == 0:
        return myvcf.header, None

    # check if tumor_sample_name exist in vcf and if yes returns the index value for future use
    coltumorsample = get_sample_column_index_from_vcf_header(next(myvcf), tumor_sample_name_in_vcf)
    # we consumed one item from myvcf in method get_sample_column_index_from_vcf_header right above; need to reset it here
    myvcf.reset()

    # step_count is ~10% of tot_variants and round to the nearest nth value
    counter = 0
    step_count = max(int(round(tot_variants_count / 10, -(len(str(round(tot_variants_count / 10))) - 1))), 1)

    # According to cyvcf2's authors, pysam does not provide an iterable of variant records; b/c of this, using multiprocessing on myvcf is not possible
    # tuple_args = (vcf_file, tumor_bam_file, normal_bam_file, ref_gen, slop, coltumorsample, minmapq, insert_size, sigma)
    list_updated_variants = []  # init list for future updated variants
    for variant in myvcf:
        list_updated_variants.append(process_variant_record(variant, vcf_file, tumor_bam_file, normal_bam_file, ref_gen, slop, coltumorsample, minmapq, insert_size, sigma))
        counter += 1
        if counter % step_count == 0:
            logger.info("processed variants ... {}".format(counter))
        logger.debug("PROCESSING VARIANT --> " + str(variant.chrom) + ":" + str(variant.pos))
    return myvcf.header, list_updated_variants


def usage():
    """
    print usage for current script from a linux-based system
    :return:
    """
    return (
        print("USAGE: \n" + os.path.basename(__file__) + '\n-i|--vcf_file\t[M] <input file VCF file> '
                                                         '\n-t|--tumor_bam_file\t[M] <Tumor BAM/CRAM file> '
                                                         '\n-n|--normal_bam_file\t[M] <normal BAM/CRAM file> '
                                                         '\n--insert-size\t[M] <mean insert size  captured from samtools stats or picard stats>'
                                                         '\n--sigma\t[M] <standard deviation for mean insert size captured from samtools stats or picard stats>'
                                                         '\n--tumor-name\t[M] name of the tumor sample found in the vcf header line starting with ##CHROM'
                                                         '\n-s|--slop\t[O] <slop, default = 1000> '
                                                         '\n-m|--minmapq\t[O] <minimum mapping quality, default = 1> '
                                                         '\n-o|--output\t[O] output_filename(relative or full path), default will add extension "_addFlags.vcf_file" to input vcf_file or '
                                                         '"_addFlags.vcf_file.gz" if "-z" enabled'
                                                         '\n-z|--gz-out\t[O] compressed the vcf_file using bcftools view ; need bcftools in PATH, default is False ; '
                                                         'recommendation: DO NOT use if you later gene-annotate the vcf_file; '
                                                         '\n-l|--logfile\t[O] logfile, default is: "log_prepare_vcf_manta.txt"'
                                                         '\n-g|--refgen\t[O] reference_genome but [M] if CRAM file'
                                                         '\n-d|--debug\t[O] enable debug information'
                                                         '\n-c|--coltumorsample\t[O] (default = 11, column for the tumor sample ; DEPRECATED])')
    )


def main(list_argv):
    vcf = None
    tbam = None
    nbam = None
    tumor_sample_name_in_vcf = None
    slop = None
    insert_size_default = 500
    insert_size = None
    sigma_default = 100
    sigma = None
    minmapq = 1
    output_vcf = ""
    refgen = ""
    logfile = __logfile__
    debug = False
    coltumorsample = 11
    compress_vcf = False
    try:
        opts, args = getopt.getopt(list_argv, "c:di:g:hl:m:n:o:s:t:z", ["help", "debug", "vcf_file=", "tumor_bam_file=", "normal_bam_file=",
                                                                        "tumor-name=", "slop=",
                                                                        "minmapq=", "refgen=", "output=", "insert-size=",
                                                                        "sigma=", "logfile=", "coltumorsample", "gz-out"])
        logger.debug(opts)
        for opt, arg in opts:
            if opt in ('-h',):
                usage()
                exit()
            elif opt in ("-i", "--vcf_file"):
                vcf = str(arg)
            elif opt in ("-t", "--tumor_bam_file"):
                tbam = str(arg)
            elif opt in ("-n", "--normal_bam_file"):
                nbam = str(arg)
            elif opt in ("--tumor-name",):
                tumor_sample_name_in_vcf = str(arg)
            elif opt in ("-s", "--slop"):
                slop = int(float(str(arg)))
            elif opt in ("--insert-size",):
                insert_size = int(float(str(arg)))
            elif opt in ("--sigma",):
                sigma = int(float(str(arg)))
            elif opt in ("-m", "--minmapq"):
                minmapq = int(float(str(arg)))
            elif opt in ("-g", "--refgen"):
                refgen = str(arg)
            elif opt in ("-o", "--output"):
                output_vcf = str(arg)
            elif opt in ("-l", "--logfile"):
                logfile = str(arg)
            elif opt in ("-d", "--debug"):
                debug = True
            elif opt in ("-z", "--gz-out"):
                compress_vcf = True
                is_tool_in_path('bcftools')
            elif opt in ("-c", "--coltumorsample"):
                coltumorsample = int(arg)
                if coltumorsample != 10 and coltumorsample != 11:
                    raise Exception("coltumorsample MUST be either 10 or 11 in values as it represents the column number of the Tumor sample in the somatic VCF")
    except getopt.GetoptError:
        usage()
        exit(2)
    except TypeError as te:
        logger.exception(te)
        exit(2)
    try:
        if debug:
            logger.setLevel(logging.DEBUG)
        if logfile != "" and logfile is not None:
            file_handler = logging.FileHandler(logfile)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
        # Filename can be full path or relative assuming the right current folder
        logger.info('')
        logger.info('VCF is: \t' + vcf)
        logger.info('TBAM is: \t' + tbam)
        logger.info('NBAM is: \t' + nbam)
        logger.info('SLOP is: \t' + str(slop))
        logger.info('INSERT_SIZE is {} ... and SIGMA (aka std) is {}'.format(str(insert_size), str(sigma)))
        logger.info('MINMAPQ is: \t' + str(minmapq))
        logger.info('REFGEN is: \t' + refgen)
        logger.info('OUTFILE is: \t' + output_vcf)
        logger.info('LOGFILE is: \t' + str(logfile))
        logger.info('TUMOR_NAME is: \t' + str(tumor_sample_name_in_vcf))
        logger.info('Output a gz compressed vcf_file is: \t' + str(compress_vcf))
        logger.debug("WARNING: for Manta_1.6, we expect the NORMAL SAMPLE to be in the column 10 and the TUMOR sample in column 11; IF not wrong data will be added to the VCF; WARNING")

        # checking inputs
        for item in [vcf, tbam, nbam]:
            if item is None or item == '':
                raise Exception("ERROR: Mandatory input values are missing")
            if not os.path.exists(item):
                raise FileNotFoundError("ERROR: File NOT FOUND for: {}".format(item))
        if os.path.splitext(tbam)[1] == ".cram" or os.path.splitext(nbam)[1] == ".cram":
            if refgen is None or refgen == '':
                raise Exception("ERROR: Reference Genome is needed since at least one of aln file is in CRAM format ; please provide --refgen and its value ")
            elif not os.path.exists(refgen):
                raise FileNotFoundError("ERROR: File NOT FOUND for: {}".format(refgen))
        if insert_size is None:
            insert_size = insert_size_default
            sigma = sigma_default
        elif insert_size is not None and isinstance(insert_size, int):
            if sigma is None or not isinstance(sigma, int):
                raise TypeError("ERROR: sigma value for insert size MUST be provided and MUST be a numeric")
        if tumor_sample_name_in_vcf is None or not isinstance(tumor_sample_name_in_vcf, str):
            raise TypeError("ERROR: TUMOR sample name present in vcf must be provided")
        if slop is None:
            slop = insert_size + (4 * sigma)

        logger.info("INSERT_SIZE = {}, SIGMA = {}, SLOP = {}".format(insert_size, sigma, slop))
        logger.info("")
        logger.info("processing ...")
        new_header, list_var_add_fields = start_processing(vcf, tbam, nbam, refgen, slop, list_argv, minmapq, insert_size, sigma, tumor_sample_name_in_vcf)
        write_new_vcf(vcf, new_header, list_var_add_fields, output_vcf=output_vcf, logfile=logfile, compress_vcf=compress_vcf)
    except FileNotFoundError as fnf:
        logger.info("check if you assign all the correct inputs and PATH to your files")
        logger.exception(fnf)
    except Exception as e:
        logger.exception(e)


if __name__ == "__main__":
    # main(argv[0], argv[1:])
    main(argv[1:])
