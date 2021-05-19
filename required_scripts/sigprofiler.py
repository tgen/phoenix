#!/usr/bin/env python3

import os
import sys
import argparse

currentDirectory = os.getcwd()

def parse_arguments():
    """Parse arguments, validate and return the args"""

    parser = argparse.ArgumentParser(
        description='Compute mutational signatures.',
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-o', '--output', default=currentDirectory,
                        help='Path to and name of output directory - Can be relative or full path')

    parser.add_argument('-i', '--vcfpath',
                        help='Path to directory containing vcf file(s) - Can be relative or full path')

    parser.add_argument('-g', '--genome', default='GRCh38',
                        help='Optional definition of genome, defaults to GRCh38')

    parser.add_argument('-p', '--project', metavar="<arg>",
                        help='Name of the project')

    parser.add_argument('-e', '--exome', default=False, action="store_true",
                        help='Set if input is from exome')

    parser.add_argument('-t', '--threads', default=-1, type=int,
                        help='Set number of threads, default will use all threads')

    # prints help message when 0 arguments are entered
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    parser_args = parser.parse_args()

    return parser_args

if __name__ == '__main__':

    # Parse and validate arguments
    args = parse_arguments()

    from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

    matrices = matGen.SigProfilerMatrixGeneratorFunc(args.project, args.genome, args.vcfpath, plot=True, exome=args.exome, bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)

    from SigProfilerExtractor import sigpro as sig

    #
    # SBS Extraction and decomposition
    #
    if args.exome:
        input_data = args.vcfpath+"/output/SBS/"+args.project+".SBS96.exome"
    else:
        input_data = args.vcfpath+"/output/SBS/"+args.project+".SBS96.all"

    sig.sigProfilerExtractor("matrix", args.output, input_data, reference_genome=args.genome, opportunity_genome = args.genome, minimum_signatures=1, maximum_signatures=10, cpu=args.threads)


    from SigProfilerExtractor import decomposition as decomp
    signatures = args.output+"/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt"
    activities = args.output+"/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Activities/SBS96_De-Novo_Activities_refit.txt"
    samples= args.output+"/SBS96/Samples.txt"

    #to get all cosmic signatures without filtering
    decomp.decompose(signatures, activities, samples, args.output, genome_build=args.genome, verbose=False, nnls_add_penalty=0.0, nnls_remove_penalty=0.0, initial_remove_penalty=0.0, de_novo_fit_penalty=0.02)

    #
    # DBS Extraction and decomposition
    #
    if args.exome:
        input_data = args.vcfpath+"/output/DBS/"+args.project+".DBS78.exome"
    else:
        input_data = args.vcfpath+"/output/DBS/"+args.project+".DBS78.all"

    sig.sigProfilerExtractor("matrix", args.output, input_data, reference_genome=args.genome, opportunity_genome = args.genome, context_type = "78", minimum_signatures=1, maximum_signatures=10, cpu=args.threads)


    from SigProfilerExtractor import decomposition as decomp
    signatures = args.output+"/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Signatures/DBS78_De-Novo_Signatures.txt"
    activities = args.output+"/DBS78/Suggested_Solution/DBS78_De-Novo_Solution/Activities/DBS78_De-Novo_Activities_refit.txt"
    samples= args.output+"/DBS78/Samples.txt"

    #to get all cosmic signatures without filtering
    decomp.decompose(signatures, activities, samples, args.output, genome_build=args.genome, verbose=False, nnls_add_penalty=0.0, nnls_remove_penalty=0.0, initial_remove_penalty=0.0, de_novo_fit_penalty=0.02)

    #
    # ID Extraction and decomposition
    #
    if args.exome:
        input_data = args.vcfpath+"/output/ID/"+args.project+".ID83.exome"
    else:
        input_data = args.vcfpath+"/output/ID/"+args.project+".ID83.all"

    sig.sigProfilerExtractor("matrix", args.output, input_data, reference_genome=args.genome, opportunity_genome = args.genome, context_type = "83", minimum_signatures=1, maximum_signatures=10, cpu=args.threads)


    from SigProfilerExtractor import decomposition as decomp
    signatures = args.output+"/ID83/Suggested_Solution/ID83_De-Novo_Solution/Signatures/ID83_De-Novo_Signatures.txt"
    activities = args.output+"/ID83/Suggested_Solution/ID83_De-Novo_Solution/Activities/ID83_De-Novo_Activities_refit.txt"
    samples= args.output+"/ID83/Samples.txt"

    #to get all cosmic signatures without filtering
    decomp.decompose(signatures, activities, samples, args.output, genome_build=args.genome, verbose=False, nnls_add_penalty=0.0, nnls_remove_penalty=0.0, initial_remove_penalty=0.0, de_novo_fit_penalty=0.02)
