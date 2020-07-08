
from multiprocessing import Process, Queue, cpu_count
import pysam, sys
import getopt
import fileinput
from sys import argv       # Used to bring in the feature argv, variables or arguments


def main(scriptname, argv):
    vcf = '' ; tbam = '' ; nbam = '' ; slop = 1000;  ## Default values
    try:
        opts, args = getopt.getopt(argv,"hi:t:n:s:",["vcf=","tbam=", "nbam=", "slop="])
        print(opts)
    except getopt.GetoptError:
        print("USAGE: \n" + scriptname + '-i <inputfile VCF file[M]> -t <Tumor BAM file [M]> -n <normal BAM file [M]> -s <slop[O]>\n')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print("USAGE: \n" +  scriptname + ' -i <inputfile VCF file[M]> -t <Tumor BAM file [M]> -n <normal BAM file [M]> -s <slop[O]>\n')
            sys.exit()
        elif opt in ("-i", "--vcf"):
            vcf = arg
        elif opt in ("-t", "--tbam"):
            tbam = arg
        elif opt in ("-n", "--nbam"):
            nbam = arg
        elif opt in ("-s", "--slop"):
            slop = int(arg)
    ## Filename can be full path or relative assuming the right current folder
    print('VCF  is: ', vcf ) ; 
    print('TBAM is: ', tbam)
    print('NBAM is: ', nbam)
    print('SLOP is: ', slop)

    newheader,LVAR = process_files_and_get_dist_of_discordant_pairs_and_return_new_vcf_records(vcf, tbam, nbam, slop)
    write_new_vcf(vcf, newheader, LVAR)

def getDiscordantReadsDistribution_parallel_version(bam, chr, pos, slop, chr2, pos2, myqueue):
    res = getDiscordantReadsDistribution(bam, chr, pos, slop, chr2, pos2)
    myqueue.put((res))

def getDiscordantReadsDistribution(bam, chr, pos, slop, chr2, pos2):
    bamf = pysam.AlignmentFile(bam, "rb")
    iter = bamf.fetch(chr, max(pos-slop,1), pos+slop) ## TODO check the max length of the chr and pos+slop should not go out of range
    pos2lowend=pos2-slop
    pos2highend=pos2+slop
#    print(str(chr) + " <--> " + str(pos-slop) + " <--> " + str(pos) + " <--> " + str(pos+slop))
#    print(str(chr) + ' <--> ' + str(chr2))
    ## init list of coordinates
    CR = []
    for x in iter:  ## x is a read line or aka an alignment record
        if x.has_tag('MQ'):	## some alignemnt records do not have the MQ value in the mate when pair aligned; BWA-known-issue
            mq_mate=x.get_tag('MQ')
        else:
            mq_mate=60	## we assume that if the mate quality is not present it is still a good one. TODO: Re-Calculate MQ if possible
        mq=x.mapping_quality
        ## here the heart of the read selection:
        if not x.is_unmapped and not x.mate_is_unmapped and not x.is_proper_pair and not x.is_duplicate and mq>=1 and mq_mate>=1:
#            print(str(x))
#            print(str(pos2lowend), "<=", str(x.next_reference_start),"<=",str(pos2highend))
#            print(pos2lowend <= x.next_reference_start <= pos2highend)
            if pos2lowend <= x.next_reference_start <= pos2highend:
                CR.append(x.reference_start)
    if not CR:
#        print("No Discordant Pairs")
        return(-1,0)	## This mean there was no discordant reads in that region ; we can return 0,0 or -1,0
    else:
        return(max(CR)-min(CR),len(CR))


def process_files_and_get_dist_of_discordant_pairs_and_return_new_vcf_records(vcf, tbam, nbam, slop):
    tumor_bam_fn=tbam ; #"MMRF_1816_1_BM_CD138pos_T1_KHWGL_L06685.bwa.final.bam"
    normal_bam_fn=nbam ; # "MMRF_1816_1_PB_Whole_C2_KHWGL_L06684.bwa.final.bam"
    #read the input VCF file using pysam (note a bug was reported on github (issue #350) but does not impact us here as we do not use the full record information)
    myvcf=pysam.VariantFile(vcf,"r")
    # Add the new fields to the header.
    myvcf.header.formats.add("RCALT","1","Integer","Read Count supporting ALT SV RC=(DV+RV) captured by Delly")
    myvcf.header.formats.add("RDISTDISC1","1","Integer","Distribution size of Discordant pair at CHR:POS ; value -1 means no discordant reads found in captured interval")
    myvcf.header.formats.add("RDISTDISC2","1","Integer","Distribution size of Discordant pair at CHR2:THEEND ; value -1 means no discordant reads found in captured interval")
    myvcf.header.formats.add("RCDIS1","1","Integer","Number of Recaptured Discordant Pairs in Left Region with which we calculated the range distribution RDISTDISC1")
    myvcf.header.formats.add("RCDIS2","1","Integer","Number of Recaptured Discordant Pairs in Right Region with which we calculated the range distribution RDISTDISC2")

    LVAR = []  ## init list for updated variants
    slop_ori = slop  ## we need to keep the otiginal value of the slop variable
    for variant in myvcf:	# loop over the list of variants records
        variant.samples[0]['RCALT']=variant.samples[0]['DV']+variant.samples[0]['RV']
        variant.samples[1]['RCALT']=variant.samples[1]['DV']+variant.samples[1]['RV']
#        print(variant)
        ## we check if there is an overlap between the left and right region, if so we reduce the slop
        ## but we can not have the slop less than 1
        if variant.chrom ==  variant.info['CHR2']:
            if variant.alts[0]!="<INV>":
                while variant.info['ENDPOSSV']-slop <= variant.pos+slop:
                    ## TODO: avoid the out of range error here by putting: max(variant.info['ENDPOSSV']-slop,1)
                    slop=max(slop-1,1)
            else:
                while variant.info['ENDPOSSV']+slop <= variant.pos-slop:
                    ## TODO: avoid the out of range error here by putting: max(variant.pos-slop,1)
                    slop=max(slop-1,1)

        if cpu_count() > 4:
            myqueue1 = Queue() ; myqueue2 = Queue() ;myqueue3 = Queue() ;myqueue4 = Queue()
            p1 = Process(target=getDiscordantReadsDistribution_parallel_version, args=(tumor_bam_fn,variant.chrom,variant.pos,slop,variant.info["CHR2"],variant.info["ENDPOSSV"],myqueue1))
            p2 = Process(target=getDiscordantReadsDistribution_parallel_version, args=(tumor_bam_fn,variant.info["CHR2"],variant.info["ENDPOSSV"],slop,variant.chrom,variant.pos,myqueue2))
            p3 = Process(target=getDiscordantReadsDistribution_parallel_version, args=(normal_bam_fn,variant.chrom,variant.pos,slop,variant.info["CHR2"],variant.info["ENDPOSSV"],myqueue3))
            p4 = Process(target=getDiscordantReadsDistribution_parallel_version, args=(normal_bam_fn,variant.info["CHR2"],variant.info["ENDPOSSV"],slop,variant.chrom,variant.pos,myqueue4))
            p1.start() ; p2.start() ; p3.start() ; p4.start();
            p1.join() ; p2.join() ; p3.join() ; p4.join() ;
            TUMDISTCHR = myqueue1.get()
            TUMDISTCHR2 = myqueue2.get()
            NORMDISTCHR = myqueue3.get()
            NORMDISTCHR2 = myqueue4.get()
            myqueue1.close() ; myqueue2.close() ; myqueue3.close() ; myqueue4.close()
        else:
            TUMDISTCHR = getDiscordantReadsDistribution(tumor_bam_fn, variant.chrom, variant.pos, slop, variant.info['CHR2'], variant.info['ENDPOSSV'])
            TUMDISTCHR2 = getDiscordantReadsDistribution(tumor_bam_fn, variant.info['CHR2'], variant.info['ENDPOSSV'], slop, variant.chrom, variant.pos)
            NORMDISTCHR = getDiscordantReadsDistribution(normal_bam_fn, variant.chrom, variant.pos, slop, variant.info['CHR2'], variant.info['ENDPOSSV'])
            NORMDISTCHR2 = getDiscordantReadsDistribution(normal_bam_fn, variant.info['CHR2'], variant.info['ENDPOSSV'], slop, variant.chrom, variant.pos)

        variant.samples[0]['RDISTDISC1']=TUMDISTCHR[0]
        variant.samples[0]['RDISTDISC2']=TUMDISTCHR2[0]
        variant.samples[1]['RDISTDISC1']=NORMDISTCHR[0]
        variant.samples[1]['RDISTDISC2']=NORMDISTCHR2[0]
        variant.samples[0]['RCDIS1']=TUMDISTCHR[1]
        variant.samples[0]['RCDIS2']=TUMDISTCHR2[1]
        variant.samples[1]['RCDIS1']=NORMDISTCHR[1]
        variant.samples[1]['RCDIS2']=NORMDISTCHR2[1]
        LVAR.append(variant)
        slop = slop_ori ## we reinit here the slop here in case we entered one of the while loop under the << if variant.chrom ==  variant.info['CHR2']: >>
    return(myvcf.header, LVAR)


def write_new_vcf(vcf, newheader,LVAR):
    """ function to write the updated records in a new vcf"""
    with open(vcf+"_addDist.vcf", "w") as f:
        f.write(str(newheader))
        for record in LVAR:
            f.write(str(record))


if __name__ == "__main__":
    main(argv[0], argv[1:])

