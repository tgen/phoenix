import getopt
from sys import argv, exit
import multiprocessing
import pysam
import pybedtools

## capture the arguments required for the annotations;
## 1) the annotation file ; 2) the VCF to be annotated ; 3) the name of the output file [O] otherwise, same name + anno.vcf
## if further arguments are necessary we will deal with that later

def usage(scriptname):
    print(
        "USAGE: \n" + scriptname + '-i <inputfile[M]> -o <outputfilename[O]> -a <AnnotationFilename(BED-formatFile)[M]>\n')


def processArgs(scriptname, argv):
    global vcf
    global annof
    global outvcfname;
    outvcfname = "";

    try:
        opts, args = getopt.getopt(argv, "hi:a:o:", ["vcf=", "annof=", "outvcfname="])
        print(opts)
    except getopt.GetoptError:
        usage(scriptname)
        exit(2)
    for opt, arg in opts:
        if opt == '-h':
            usage(scriptname);
            exit()
        elif opt in ("-i", "--vcf"):
            vcf = arg
        elif opt in ("-a", "--annof"):
            annof = arg
        elif opt in ("-o", "--outvcfname"):
            outvcfname = arg

    ## check args values and input requirements
    if vcf is None or annof is None:
        usage(scriptname);
        exit(2)
    if outvcfname is None or outvcfname == "":
        outvcfname = ''.join(str(vcf) + ".anno.vcf")

    print('VCF  is: ',
          vcf);  ## iput file is a list of the filename we need to process in the loop; Filename can be full path or relative assuming the right current folder
    print('annof is: ', annof)
    print('outvcfname is: ', outvcfname)


def getAnno(LR, RR, genes):
    """
    :param LR:  BedTool object with only one position defined for the Left Breakpoint of the current SV ; LR stands for Left Region
    :param RR: BedTool object with only one position defined for the Right Breakpoint of the current SV ; RR stands for Right Region
    :param genes: BedTool object with all the gene annotations in bedformat
    :return: tuple of strings (gene_name(s),dist,strand)
    """
    return ([getClosestUpstream(LR, genes), getIsec(LR, genes), getClosesetDownstream(LR, genes),
             getClosestUpstream(RR, genes), getIsec(RR, genes), getClosesetDownstream(RR, genes)])


def getIsec(locus, genes):
    """
    :param locus:  BedTool object with only one position defined
    :param genes: BedTool object with all the gene annotations in bedformat
    :return: tuple of strings (gene_name(s),dist,strand)
    """
    isec = locus.intersect(genes, wao=True)
    if isec is None or isec == "":
        return (".___.", "0", ".")  ## print("NO Intersection with any gene")
    gene = set(isec.to_dataframe().iloc[0:1, 6]).pop()  ## here we keep only the first gene in the list
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
        return ((".___.", "0", "."))  ## print("NO Closest Gene upstream with any gene")
    gene = set(theclosest.to_dataframe().iloc[0:1, 6]).pop()  ## here we keep only the first gene in the list
    if gene == ".":
        return (".___.", "0", ".")
    dist = set(theclosest.to_dataframe().iloc[0:1, 9]).pop()
    strand = set(theclosest.to_dataframe().iloc[0:1, 8]).pop()
    return ((str(gene), str(dist), str(strand)))


def getClosesetDownstream(locus, genes):
    """
    :param locus:  BedTool object with only one position defined
    :param genes: BedTool object with all the gene annotations in bedformat
    :return: tuple of strings (gene_name(s),dist,strand)
    """
    theclosest = locus.closest(genes, io=True, fd=True, D="a", d=True, t="first")
    if theclosest is None or theclosest == "":
        # print("NO Closest Gene Downstream with any gene")
        return ((".___.", str(0), "."))
    gene = set(theclosest.to_dataframe().iloc[0:1, 6]).pop()  ## here we keep only the first gene in the list
    if gene == ".":
        return (".___.", "0", ".")
    dist = set(theclosest.to_dataframe().iloc[0:1, 9]).pop()
    strand = set(theclosest.to_dataframe().iloc[0:1, 8]).pop()
    return ((str(gene), str(dist), str(strand)))


def getGenesOverlappinRegion(rec, genes):
    """
    :param rec:  vcf record pysam-formatted
    :param genes: BedTool object with all the gene annotations in bedformat
    :return: tuple of strings (gene_name(s))
    """
    NOGENE=str(set([".___."]))
    if rec.info['SVTYPE']=="TRA":
        return NOGENE
    chr1 = rec.chrom
    POS1 = rec.pos
    POS2 = rec.info['ENDPOSSV']

    if int(POS2) < int(POS1):
        POS1 = rec.info['ENDPOSSV'] ; POS2 = rec.pos ;
    locus = pybedtools.BedTool(' '.join([chr1, str(POS1-1), str(POS2)]), from_string=True)
    isec = locus.intersect(genes, wao=True)
    if isec is None or isec == "":
        return NOGENE  ## print("NO Intersection with any gene")
    gene = set(isec.to_dataframe().iloc[0::, 6])  ## here we get ALL the Genes in column 6
    # strand = set(isec.to_dataframe().iloc[0:, 8]).pop()
    if gene == "." or gene == "{'.'}":
        return NOGENE
    return (str(gene))


def getLenghtRegion(rec):
    """
    :param rec:  vcf record pysam-formatted
    :return: lenght of the region of interest; pos2-pos1;
    """
    if rec.info['SVTYPE']=="TRA":
        return (-1);
    return abs(rec.info['ENDPOSSV']-rec.pos)


def newHeaders(myvcf):
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
    return myvcf

def write_new_vcf(vcf, newheader, LVAR, outfilename=""):
    if outfilename is None or outfilename == "":
        outfilename = "".join([str(vcf), ".parallel.anno.vcf"])
    with open(outfilename, "w") as f:
        f.write(str(newheader))
        print("in write_new_vcf")
        for record in LVAR:
            f.write(str(record))


def main(scriptname, argv):
    processArgs(scriptname, argv)
    import time
    from pybedtools import BedTool
    global genes
    genes = BedTool(annof).sort()
    ## we read the entire vcf file at once or we stream the record and process the extraction of the required info to annotate the breakpoints
    vcfr = pysam.VariantFile(vcf, "r")
    ## updating header with new info Tags
    vcfr = newHeaders(vcfr)

    my_list = [variant for variant in vcfr]
    from multiprocessing import Pool
    import logging
    logger = multiprocessing.log_to_stderr()
    logger.setLevel(multiprocessing.SUBDEBUG)
    pool = Pool(processes=12)              # start 4 worker processes
    result_list = pool.map(func=annotateVariant_parallel, iterable=iter(my_list))
    print(logger)
    with open('outFile.vcf', 'w') as out_file:
        out_file.writelines(result_list)

def annotateVariant_parallel(i):
    res = annotateVariant(i)
    return res


def annotateVariant(i):
    rec = my_list[i]
    chr1 = rec.chrom
    chr2 = rec.info['CHR2']
    POS1 = rec.pos
    POS2 = rec.info['ENDPOSSV']

    rec.info['LENSV'] = getLenghtRegion(rec)
    rec.info['GENELISTbtwBP'] = getGenesOverlappinRegion(rec, genes).replace(" ","")

    # creating the BedTool object with only one locus defined within (LR is Left Region of the SV, RR is Right Region)
    LR = pybedtools.BedTool(' '.join([chr1, str(POS1 - 1), str(POS1)]), from_string=True)
    RR = pybedtools.BedTool(' '.join([chr2, str(POS2 - 1), str(POS2)]), from_string=True)
    ## getting the gene annotations for the isec, the upstream and the downstream genes and ADDING these to the VCF Record (rec)
    rec.info['RGENUPS'], rec.info['RGENISEC'], rec.info['RGENDNS'], \
    rec.info['LGENUPS'], rec.info['LGENISEC'], rec.info['LGENDNS'] = getAnno(LR, RR, genes)
    return str(rec)

##@@@@@@@@
## MAIN ##
##@@@@@@@@
if __name__ == "__main__":
    # main(argv[0], argv[1:])
    # exit()

    processArgs(argv[0], argv[1:])
    import time
    from pybedtools import BedTool
    global genes
    print("sorting the annotation bedfile ... ")
    genes = BedTool(annof).sort()
    ## we read the entire vcf file at once or we stream the record and process the extraction of the required info to annotate the breakpoints
    vcfr = pysam.VariantFile(vcf, "r")
    ## updating header with new info Tags
    vcfr = newHeaders(vcfr)
    newheader = vcfr.header
    ## init list for new Variants
    res = []
    global my_list
    my_list = [variant for variant in vcfr]

    my_range = [i for i in range(len(my_list))]
    from multiprocessing import Pool
    import logging
    try:
        logger = multiprocessing.log_to_stderr()
        logger.setLevel(multiprocessing.SUBDEBUG)
        ncpus = multiprocessing.cpu_count()
        pool = Pool(processes=ncpus)  # start X worker processes
        result_list = pool.map(func=annotateVariant_parallel, iterable=my_range)
        time.sleep(2)
        # print(logger)
    finally:
        pool.close()
        pool.join()
    # print("results list")
    # print(type(result_list))
    # print(len(result_list))
    from operator import itemgetter
    sorted(result_list,key=itemgetter(1,2))
    # print([str(rec) for rec in result_list[0:2]])
    write_new_vcf(vcf, newheader, result_list, outvcfname)


## Example of vcf record
# X       96655464 DEL00008569 T <DEL> . PASS
# PRECISE;SVTYPE=DEL;SVMETHOD=EMBL.DELLYv0.7.5;CHR2=X;ENDPOSSV=96781481;INSLEN=0;HOMLEN=3;PE=21;
# MAPQ=60;CT=3to5;CIPOS=-4,4;CIEND=-4,4;SR=7;SRQ=1;
# CONSENSUS=ACAGTGTACTATGTGATGTTTTGACATATGTATACCAAATCCATTTAGCACTTGGTAACAAAAGGTAAGAATAGACATTGAATACTGTACTATTTTTA;
# CE=1.87385;RDRATIO=0.257492;SOMATIC
# GT:GL:GQ:FT:RCL:RC:RCR:CN:DR:DV:RR:RV:RCALT:RDISTDISC1:RDISTDISC2:RCDIS1:RCDIS2
# 1/1:-17.8973,-1.80346,0:18:PASS:2439:1227:2350:1:2:21:0:6:27:690:704:21:21
# 0/0:0,-1.20172,-11.4976:12:LowQual:3039:6207:3199:2:16:0:4:0:0:-1:-1:0:0
