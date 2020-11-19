library(dplyr)
library(tidyverse)
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))

#Read in the individual flatFiles

option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-o", "--outfile"), action="store", default="combinedIgTxCall.txt",
              help="output results file"),
  make_option(c("-s", "--specimen"),action="store", default="SAMPLE"),
  make_option(c("-p", "--pairoscope_file"),action="store",
              help="Flat file containing Ig Calls using Pairoscope"),
  make_option(c("-m", "--manta_file"), action="store",
              help="Flat file containing Ig Calls using manta"),
  make_option(c("-g", "--gammit_file"), action="store",
              help="Flat file containing Ig Calls using Gammit"),
  make_option(c("-c", "--count"), action="store",type="integer", default=2,
              help="Minimum caller count [default %default]",
              metavar="number")
   )

opt <- parse_args(OptionParser(option_list=option_list))

#Test case/example output from 3 Ig callers
#opt$gammit_file="/Volumes/MMRF/cell_lines/phoenix_results/KMS34_GM12878/exome/tumor_only_structural_calls/gammit/KMS34_GM12878_p2_CL_1_T1_TWTK1/KMS34_GM12878_p2_CL_1_T1_TWTK1_DEX_IgTx_GA_Summary.txt"
#opt$pairoscope_file="/Volumes/MMRF/cell_lines/phoenix_results/KMS34_GM12878/exome/somatic_structural_calls/pairoscope/KMS34_GM12878_p2_CL_1_T1_TWTK1/KMS34_GM12878_p2_CL_1_T1_TWTK1.bwa_pairoscope_igtx_calls.txt"
#opt$manta_file="/Volumes/MMRF/cell_lines/phoenix_results/KMS34_GM12878/exome/tumor_only_structural_calls/manta/KMS34_GM12878_p2_CL_1_T1_TWTK1/KMS34_GM12878_p2_CL_1_T1_TWTK1.bwa_mantaTxCalls.txt"
#opt$outfile = "/Users/snasser/Desktop/ttsjk.txt"
#opt$count = 2
#opt$specimen = "KMS34"

write("Check input files...\n", stderr())
if(is.null(opt$pairoscope_file) || !file.exists(as.character(opt$pairoscope_file)))
{
     write("Please provide Pairoscope Ig Tx files to procees..\n", stderr())
} else if (is.null(opt$manta_file) || !file.exists(as.character(opt$manta_file)))
{
    write("Please provide manta Ig Tx files to procees..\n", stderr())
} else if (is.null(opt$gammit_file) || !file.exists(as.character(opt$gammit_file)))
{
  write("Please provide Gammit Ig Tx files to procees..\n", stderr())
} else
{
  write("Process Data...\n", stderr())
  combined_calls<-NULL
  #read from from  file
  pairoscope= read.table(file=opt$pairoscope_file, header = TRUE,sep = '\t')
  pair_calls=pairoscope %>% select(ends_with("Call"))
  pair_calls <- pair_calls %>% rename_all(list(~ str_replace(., "CALL", "CALL_Pairoscope")))

  specimen = tibble(Specimen=opt$specimen)
  #manta
  manta=read.table(file=opt$manta_file,header = TRUE,sep = '\t')
  manta_calls = manta %>% select(ends_with("Called"))
  manta_calls <- manta_calls %>% rename_all(list(~ str_replace(.,"Target_Called", "CALL_Manta")))

  #gammit
  gammit = read.table(file=opt$gammit_file,header = TRUE,sep = '\t')
  gammit_calls = gammit %>% select(ends_with("Call"))
  gammit_calls <- gammit_calls %>% rename_all(list(~ str_replace(.,"Call", "CALL_Gammit")))

  #merge

  combined_calls=vctrs::vec_cbind(!!!list(specimen,gammit_calls, manta_calls,pair_calls))

  combined_calls= combined_calls %>%
    mutate (NSD2_CALLER_COUNT = combined_calls %>% select(starts_with("NSD2"))  %>% sum(),
      NSD2_Summary_CALL = if_else(NSD2_CALLER_COUNT >= opt$count, 1, 0),
      MAF_CALLER_COUNT = combined_calls %>% select(starts_with("MAF"))  %>% sum(),
      MAF_Summary_CALL = if_else(MAF_CALLER_COUNT >= opt$count, 1, 0),
      MAFA_CALLER_COUNT = combined_calls %>% select(starts_with("MAFA"))  %>% sum(),
      MAFA_Summary_CALL = if_else(MAFA_CALLER_COUNT >= opt$count, 1, 0),
      MAFB_CALLER_COUNT = combined_calls %>% select(starts_with("MAFB"))  %>% sum(),
      MAFB_Summary_CALL = if_else(MAFB_CALLER_COUNT >= opt$count, 1, 0),
      MYC_CALLER_COUNT = combined_calls %>% select(starts_with("MYC"))  %>% sum(),
      MYC_Summary_CALL = if_else(MYC_CALLER_COUNT >= opt$count, 1, 0),
      CCND1_CALLER_COUNT = combined_calls %>% select(starts_with("CCND1"))  %>% sum(),
      CCND1_Summary_CALL = if_else(CCND1_CALLER_COUNT >= opt$count, 1, 0),
      CCND2_CALLER_COUNT = combined_calls %>% select(starts_with("CCND2"))  %>% sum(),
      CCND2_Summary_CALL = if_else(CCND2_CALLER_COUNT >= opt$count, 1, 0),
      CCND3_CALLER_COUNT = combined_calls %>% select(starts_with("CCND3"))  %>% sum(),
      CCND3_Summary_CALL = if_else(CCND3_CALLER_COUNT >= opt$count, 1, 0)
        )
  combined_calls=combined_calls[,order(colnames(combined_calls), decreasing = TRUE)]

  # combined_calls=combined_calls %>% relocate(Specimen)
  write("Save results...\n", stderr())
  write_tsv(combined_calls, opt$outfile, append = FALSE,na="NA")
  write("Done.\n", stderr())
  }
