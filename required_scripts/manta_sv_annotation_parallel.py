#!/usr/bin/env python3

import getopt
import os
import multiprocessing
import pysam
import time
import logging
from sys import argv, exit
from pybedtools import BedTool, cleanup, set_tempdir, get_tempdir


def getAnno(LR, RR, genes):
    """
    :param LR:  BedTool object with only one position defined for the Left Breakpoint of the current SV ; LR stands for Left Region
    :param RR: BedTool object with only one position defined for the Right Breakpoint of the current SV ; RR stands for Right Region
    :param genes: BedTool object with all the gene annotations in bedformat
    :return: tuple of strings (gene_name(s),dist,strand)
    """
    return ([getClosestUpstream(LR, genes), getIsec(LR, genes), getClosestDownstream(LR, genes),
             getClosestUpstream(RR, genes), getIsec(RR, genes), getClosestDownstream(RR, genes)])


def getIsec(locus, genes):
    """
    :param locus:  BedTool object with only one position defined
    :param genes: BedTool object with all the gene annotations in bedformat
    :return: tuple of strings (gene_name(s),dist,strand)
    """
    isec = locus.intersect(genes, wao=True)
    if isec is None or isec == "":
        return (".___.", "0", ".")  # print("NO Intersection with any gene")
    gene = set(isec.to_dataframe().iloc[0:1, 6]).pop()  # here we keep only the first gene in the list
    strand = set(isec.to_dataframe().iloc[0:1, 8]).pop()
    if gene == ".":
        return (".___.", "0", ".")
    return ((str(gene), str(0), str(strand)))


def getClosestUpstream(locus, genes):
    """
    :param locus:  BedTool object with only one position defined
    :param genes: BedTool object with all the gene annotations in bedformat
    :return: tuple of strings (gene_name(s),dist,strand)
    """
    theclosest = locus.closest(genes, io=True, fu=True, D="a", d=True, t="first")
    if theclosest is None or theclosest == "":
        return (".___.", "0", ".")  ## print("NO Closest Gene upstream with any gene")
    gene = set(theclosest.to_dataframe().iloc[0:1, 6]).pop()  ## here we keep only the first gene in the list
    if gene == ".":
        return (".___.", "0", ".")
    dist = set(theclosest.to_dataframe().iloc[0:1, 9]).pop()
    strand = set(theclosest.to_dataframe().iloc[0:1, 8]).pop()
    return (str(gene), str(dist), str(strand))


def getClosestDownstream(locus, genes):
    """
    :param locus:  BedTool object with only one position defined
    :param genes: BedTool object with all the gene annotations in bedformat
    :return: tuple of strings (gene_name(s),dist,strand)
    """
    theclosest = locus.closest(genes, io=True, fd=True, D="a", d=True, t="first")
    if theclosest is None or theclosest == "":
        # print("NO Closest Gene Downstream with any gene")
        return (".___.", str(0), ".")
    gene = set(theclosest.to_dataframe().iloc[0:1, 6]).pop()  ## here we keep only the first gene in the list
    if gene == ".":
        return (".___.", "0", ".")
    dist = set(theclosest.to_dataframe().iloc[0:1, 9]).pop()
    strand = set(theclosest.to_dataframe().iloc[0:1, 8]).pop()
    return (str(gene), str(dist), str(strand))


def getGenesOverlappingRegion(rec, genes):
    """
    :param rec:  pysam vcf record
    :param genes: BedTool object with all the gene annotations in bedformat
    :return: tuple of strings (gene_name(s))
    """

    NOGENE=str(set([".___."]))
    if rec.info['SVTYPE'] == "TRA" or ("TRA" in rec.info.keys() and rec.info['TRA']):
        return NOGENE
    chr1 = rec.chrom
    pos1 = rec.pos
    pos2 = rec.info['ENDPOSSV']

    if int(pos2) < int(pos1):
        pos1 = rec.info['ENDPOSSV']
        pos2 = rec.pos
    locus = BedTool(' '.join([chr1, str(pos1-1), str(pos2)]), from_string=True)
    isec = locus.intersect(genes, wao=True)
    if isec is None or isec == "":
        return NOGENE  # print("NO Intersection with any gene")
    gene = set(isec.to_dataframe().iloc[0::, 6])  # here we get ALL the Genes in column 6
    # strand = set(isec.to_dataframe().iloc[0:, 8]).pop()
    if gene == "." or gene == "{'.'}":
        return NOGENE
    return str(gene)


def getLengthRegion(rec):
    """
    :param rec:  pysam vcf record
    :return: length of the region of interest; pos2-pos1;
    """
    if rec.info['SVTYPE'] == "TRA" or ("TRA" in rec.info.keys() and rec.info['TRA']):
        return -1
    return abs(rec.info['ENDPOSSV']-rec.pos)


def newHeaders(myvcf, arg_cmd):
    """
    Update the VCF header with newly added fields and other information
    :param myvcf:  pysam vcf object VariantFile
    :return: updated myvcf object;
    """
    myvcf.header.info.add("RGENUPS", "3", "String",
                          "Gene,DistToBreakPoint,Strand for Upstream Gene to the Right BreakPoint")
    myvcf.header.info.add("RGENISEC", "3", "String",
                          "Gene,DistToBreakPoint,Strand for Gene Intersecting the Right BreakPoint")
    myvcf.header.info.add("RGENDNS", "3", "String",
                          "Gene,DistToBreakPoint,Strand for Downstream Gene to the Right BreakPoint")
    myvcf.header.info.add("LGENUPS", "3", "String",
                          "Gene,DistToBreakPoint,Strand for Upstream Gene to the Left BreakPoint")
    myvcf.header.info.add("LGENISEC", "3", "String",
                          "Gene,DistToBreakPoint,Strand for Gene Intersecting the Left BreakPoint")
    myvcf.header.info.add("LGENDNS", "3", "String",
                          "Gene,DistToBreakPoint,Strand for Downstream Gene to the Left BreakPoint")
    myvcf.header.info.add("GENELISTbtwBP", ".", "String",
                          "List of genes between the two defined breakpoint for DEL, DUP, INS, and INV")
    myvcf.header.info.add("LENSV", "1", "Integer",
                          "Length of the SV for DEL, DUP, INS, and INV; TRA will have lenght of 0")
    if arg_cmd is not None:
        str_arg_cmd = ' '.join([str(x) for x in arg_cmd[1:]])
        myvcf.header.add_line("##CommandAnnotationSVvcf='{}'".format(str(str_arg_cmd)))
    return myvcf


def write_new_vcf(vcf, newheader, lvar, out_filename=None):
    """
    :param vcf: full (or relative) path to input vcf file
    :param newheader: Pysam's Header object
    :param lvar: list of Pysam's Variant Records Objects
    :param out_filename: Name of the output file (not the path)
    :return: None
    :rtype: None
    """

    if out_filename is None or out_filename == "":
        out_filename = "".join([str(vcf), ".parallel.anno.vcf"])
    with open(out_filename, "w") as f:
        f.write(str(newheader))
        for record in lvar:
            f.write(str(record))


def annotateVariant_parallel(i):
    """
    Needed to run variant annotation in parallel
    :param i: index to loop over the list of variant
    :type i: integer
    :return: tuple of fields updated
    :rtype: tuple
    """
    return annotateVariant(i)


def annotateVariant(i):
    """
    capture and add annotation to a variant
    :param i: index
    :type i: integer
    :return: tuple of added or updated fields in the variant record
    :rtype: tuple
    """
    rec = my_list[i]
    chr1 = rec.chrom
    chr2 = rec.info['CHR2']
    pos1 = rec.pos
    pos2 = rec.info['ENDPOSSV']

    rec.info['LENSV'] = getLengthRegion(rec)
    if get_genes_btw_bps:
        rec.info['GENELISTbtwBP'] = getGenesOverlappingRegion(rec, genes).replace(" ", "")
    try:
        # creating the BedTool object with only one locus defined within (LR is Left Region of the SV, RR is Right Region)
        LR = BedTool(' '.join([chr1, str(pos1 - 1), str(pos1)]), from_string=True)
        RR = BedTool(' '.join([chr2, str(pos2 - 1), str(pos2)]), from_string=True)
        # getting the gene annotations for the isec, the upstream and the downstream genes and ADDING these to the VCF Record (rec)
        rec.info['RGENUPS'], rec.info['RGENISEC'], rec.info['RGENDNS'], rec.info['LGENUPS'], rec.info['LGENISEC'], rec.info['LGENDNS'] = getAnno(LR, RR, genes)
    except Exception as exc:
        print(exc)
    return str(rec)


class dataForAnno:
    """
    Class to capture all the data and input in one object in order to simplify the passage of arguments
    and to print the inputs and their respective values
    """
    def __init__(self, vcf_file, annof, outvcfname, get_genes_btw_bps, tempdir="None", threads=multiprocessing.cpu_count()):
        self.vcf_file = vcf_file
        self.annof = annof
        self.outvcfname = outvcfname
        self.get_genes_btw_bps = get_genes_btw_bps
        self.tempdir = tempdir
        self.threads = threads


    def print_inputs_summary(self):
        print('VCF\tis:\t', self.vcf_file)
        print('annof\tis:\t', self.annof)
        print('outvcfname\tis:\t', self.outvcfname)
        print('get_inbtw_genes\tis:\t', self.get_genes_btw_bps)
        print('threads\tis:\t', self.threads)
        print("pyBedtools Temporary Dir:\t{}".format(str(get_tempdir())))
        print("pyBedtools Temporary Dir:\t{}".format(str(self.tempdir)))


def usage(scriptname):
    """
    Show how to use and what options to the script are available
    :param scriptname: current name of the script
    :type scriptname: string
    :return: print usage
    :rtype: string
    """
    print(
        "USAGE: \n" + scriptname + ' -a <AnnotationFilename(BED-formatFile)[M] -i <inputfile[M]> -o <outputfilename[O]> '
                                   '-t <numberOfThreads[O]> --tempdir <TemporaryDirectoryForPybedtools[O/default=None]> --add-genes-in-between[O/Flag/default=False] \n')


def processArgs(scriptname, argv):
    """
    Process the given arguments and options and checks inputs
    :param scriptname: scriptname
    :type scriptname: string
    :param argv: list of arguments
    :type argv: list
    :return: object of Class dataForAnno
    :rtype: object Class dataForAnno
    """
    global vcf, vcf_file
    global annof
    global outvcfname
    global get_genes_btw_bps
    vcf = None
    annof = None
    outvcfname = None
    #tempdir = None
    get_genes_btw_bps = False
    threads = multiprocessing.cpu_count()

    try:
        opts, args = getopt.getopt(argv, "hi:a:o:t:",
                                   ["vcf=", "annof=", "outvcfname=", "threads=",
                                    "add-genes-in-between", "tempdir="])
        print(opts)
    except getopt.GetoptError:
        usage(scriptname)
        exit(2)
    for opt, arg in opts:
        if opt == '-h':
            usage(scriptname)
            exit()
        elif opt in ("-i", "--vcf"):
            vcf = str(arg)
        elif opt in ("-a", "--annof"):
            annof = str(arg)
        elif opt in ("-o", "--outvcfname"):
            outvcfname = str(arg)
        elif opt in ("--tempdir",):
            tempdir = str(arg)
        elif opt in ("--add-genes-in-between",):
            get_genes_btw_bps = True
        elif opt in ("-t", "--threads"):
            threads = int(arg)

    # check args values and input requirements
    if vcf is None or annof is None:
        usage(scriptname)
        raise FileNotFoundError("VCF and ANNOF files are MANDATORY")
        exit(2)
    if outvcfname is None or outvcfname == "":
        outvcfname = ''.join(str(vcf) + ".anno.vcf")

    odfa = dataForAnno(vcf, annof, outvcfname, get_genes_btw_bps, tempdir=tempdir, threads=threads)
    odfa.print_inputs_summary()
    return odfa


if __name__ == "__main__":
    # main(argv[0], argv[1:])
    # #@@@@@@@@
    # # MAIN ##
    # #@@@@@@@@
    global vcf, vcf_file
    global annof
    global outvcfname
    global get_genes_btw_bps
    global genes
    global my_list
    pool = None
    result_list = None
    dfa = processArgs(argv[0], argv[1:])

    try:
        # ---------------#
        # define logging
        logger = multiprocessing.log_to_stderr()
        formatter = logging.Formatter('%(levelname)s %(asctime)-15s %(module)s %(lineno)d\t %(message)s')
        # create file_handlers with different settings from generic handling
        file_handler = logging.FileHandler(os.path.join(os.path.dirname(dfa.vcf_file), "log_mp_annotate.out"))
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        # logger.setLevel(multiprocessing.SUBDEBUG)
        logger.setLevel(logging.INFO)
        # -----------------#

        if dfa.tempdir is not None and os.path.exists(dfa.tempdir):
            set_tempdir(dfa.tempdir)
        elif not os.path.exists(dfa.tempdir):
            raise FileNotFoundError("\nDIR NOT FOUND please create: " + str(dfa.tempdir) + "\n")
        print("sorting the annotation bed file ... ")
        genes = BedTool(dfa.annof).sort()
        # we read the entire vcf file at once or we stream the record and process the extraction of the required info to annotate the breakpoints
        vcf_recs = pysam.VariantFile(dfa.vcf_file, "r")
        # updating header with new info Tags
        vcf_recs = newHeaders(vcf_recs, argv)
        newheader = vcf_recs.header
        # init list for new Variants
        res = []
        my_list = [variant for variant in vcf_recs]
        my_range = [i for i in range(len(my_list))]
        pool = multiprocessing.Pool(processes=dfa.threads)  # start X worker processes
        result_list = pool.map(func=annotateVariant_parallel, iterable=my_range)
        time.sleep(1)
    except FileNotFoundError as fnf:
        print(str(fnf))
    except NotADirectoryError as nade:
        print(str(nade))
    except IOError as ioe:
        print("ERROR I/O: " + str(ioe))
    except Exception as e:
        print(e)
    finally:
        if pool is not None:
            pool.close()
            pool.join()
    try:
        from operator import itemgetter
        sorted(result_list, key=itemgetter(1, 2))
        # print([str(rec) for rec in result_list[0:2]])
        write_new_vcf(dfa.vcf_file, newheader, result_list, dfa.outvcfname)
        print("vcf outfile is: {}".format(str(dfa.outvcfname)))
    except KeyError as kerr:
        print(str(kerr))
    except FileNotFoundError as fnf:
        print(str(fnf))
    except Exception as e:
        print(e)
    finally:
        # remove all the temporary file created by bedtools
        cleanup(remove_all=True)

#  Example of vcf record
# X       96655464 DEL00008569 T <DEL> . PASS
# PRECISE;SVTYPE=DEL;SVMETHOD=EMBL.DELLYv0.7.5;CHR2=X;ENDPOSSV=96781481;INSLEN=0;HOMLEN=3;PE=21;
# MAPQ=60;CT=3to5;CIPOS=-4,4;CIEND=-4,4;SR=7;SRQ=1;
# CONSENSUS=ACAGTGTACTATGTGATGTTTTGACATATGTATACCAAATCCATTTAGCACTTGGTAACAAAAGGTAAGAATAGACATTGAATACTGTACTATTTTTA;
# CE=1.87385;RDRATIO=0.257492;SOMATIC
# GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV:RCALT:RDISTDISC1:RDISTDISC2:RCDIS1:RCDIS2
# 1/1:-17.8973,-1.80346,0:18:PASS:2439:1227:2350:1:2:21:0:6:27:690:704:21:21
# 0/0:0,-1.20172,-11.4976:12:LowQual:3039:6207:3199:2:16:0:4:0:0:-1:-1:0:0
