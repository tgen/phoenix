library(dplyr)
library(tidyverse)
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))

# Initialize list of matrix object
genes <- c("NSD2", "CCND3", "MYC", "MAFA", "CCND1", "CCND2", "MAF", "MAFB")
callers <- c()
header <- c("Specimen")

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
              metavar="number"),
  make_option(c("--genes"), action="store",
              help="Flat file containing genes of interest (not required)")
)

opt <- parse_args(OptionParser(option_list=option_list))

write("Check input files...\n", stderr())
if(!is.null(opt$pairoscope_file))
{
  callers <- append(callers, "Pairoscope")
  if(!file.exists(as.character(opt$pairoscope_file)))
  {
    write("Pairoscope Ig Tx file not found...\n", stderr())
  }
}

if (!is.null(opt$manta_file))
{
  callers <- append(callers, "Manta")
  if (!file.exists(as.character(opt$manta_file)))
  {
    write("Manta Ig Tx file not found...\n", stderr())
  }
}

if (!is.null(opt$gammit_file))
{
  callers <- append(callers, "Gammit")
  if (!file.exists(as.character(opt$gammit_file)))
  {
    write("Gammit Ig Tx file not found...\n", stderr())
  }
}

if (!is.null(opt$genes))
{
  genes = read.table(file=opt$genes)
}

for (gene in genes) {
  header <- append(header, paste(gene, "Summary_CALL", sep = "_"))
  header <- append(header, paste(gene, "IgSource", sep = "_"))
  header <- append(header, paste(gene, "CALLER_COUNT", sep = "_"))
  for (caller in callers) {
    header <- append(header, paste(gene, "CALL", caller, sep = "_"))
  }
}

init_matrix <- matrix(0, ncol = length(header), nrow = 1)
colnames(init_matrix) <- header
init_table <- as_tibble(init_matrix)
init_table <- init_table %>% mutate_at(1, as.character)

write("Processing Data...\n", stderr())
specimen = tibble(Specimen=opt$specimen)
combined_calls<-NULL
call_list <- list(specimen)

#pairoscope
if(!is.null(opt$pairoscope_file) && file.exists(as.character(opt$pairoscope_file)))
{
  pairoscope= read.table(file=opt$pairoscope_file, header = TRUE,sep = '\t')
  pair_calls=pairoscope %>% select(ends_with("Call"))
  pair_calls <- pair_calls %>% rename_all(list(~ str_replace(., "CALL", "CALL_Pairoscope")))
  pair_source = pairoscope %>% select(ends_with("IGSOURCE"))
  pair_source <- pair_source %>% rename_all(list(~ str_replace(.,"IGSOURCE", "IgSource")))
  pair_source$NSD2_IgSource = ifelse(pair_calls$NSD2_CALL_Pairoscope==1, pair_source$NSD2_IgSource,0)
  pair_source$CCND1_IgSource = ifelse(pair_calls$CCND1_CALL_Pairoscope==1, pair_source$CCND1_IgSource,0)
  pair_source$CCND2_IgSource = ifelse(pair_calls$CCND2_CALL_Pairoscope==1, pair_source$CCND2_IgSource,0)
  pair_source$CCND3_IgSource = ifelse(pair_calls$CCND3_CALL_Pairoscope==1, pair_source$CCND3_IgSource,0)
  pair_source$MYC_IgSource = ifelse(pair_calls$MYC_CALL_Pairoscope==1, pair_source$MYC_IgSource,0)
  pair_source$MAF_IgSource = ifelse(pair_calls$MAF_CALL_Pairoscope==1, pair_source$MAF_IgSource,0)
  pair_source$MAFA_IgSource = ifelse(pair_calls$MAFA_CALL_Pairoscope==1, pair_source$MAFA_IgSource,0)
  pair_source$MAFB_IgSource = ifelse(pair_calls$MAFB_CALL_Pairoscope==1, pair_source$MAFB_IgSource,0)

  call_list <- append(call_list, pair_calls)
}

#manta
if(!is.null(opt$manta_file) && file.exists(as.character(opt$manta_file)))
{
  manta=read.table(file=opt$manta_file,header = TRUE,sep = '\t')
  manta_calls = manta %>% select(ends_with("Called"))
  manta_calls <- manta_calls %>% rename_all(list(~ str_replace(.,"Target_Called", "CALL_Manta")))
  manta_source = manta %>% select(ends_with("Ig_Loci"))
  manta_source <- manta_source %>% rename_all(list(~ str_replace(.,"Ig_Loci", "IgSource")))
  call_list <- append(call_list, manta_calls)
}

#gammit
if(!is.null(opt$gammit_file) && file.exists(as.character(opt$gammit_file)))
{
  gammit = read.table(file=opt$gammit_file,header = TRUE,sep = '\t')
  gammit_calls = gammit %>% select(ends_with("Call"))
  gammit_calls <- gammit_calls %>% rename_all(list(~ str_replace(.,"Call", "CALL_Gammit")))

  gammit_source = gammit %>% select(ends_with("Ig_Loci"))
  gammit_source <- gammit_source %>% rename_all(list(~ str_replace(.,"Ig_Loci", "IgSource")))
  call_list <- append(call_list, gammit_calls)
}

#merge
combined_calls=vctrs::vec_cbind(!!!call_list)

combined_calls= combined_calls %>%
  mutate (NSD2_CALLER_COUNT = combined_calls %>% select(starts_with("NSD2_"))  %>% sum(),
          NSD2_Summary_CALL = if_else(NSD2_CALLER_COUNT >= opt$count, 1, 0),
          MAF_CALLER_COUNT = combined_calls %>% select(starts_with("MAF_"))  %>% sum(),
          MAF_Summary_CALL = if_else(MAF_CALLER_COUNT >= opt$count, 1, 0),
          MAFA_CALLER_COUNT = combined_calls %>% select(starts_with("MAFA_"))  %>% sum(),
          MAFA_Summary_CALL = if_else(MAFA_CALLER_COUNT >= opt$count, 1, 0),
          MAFB_CALLER_COUNT = combined_calls %>% select(starts_with("MAFB_"))  %>% sum(),
          MAFB_Summary_CALL = if_else(MAFB_CALLER_COUNT >= opt$count, 1, 0),
          MYC_CALLER_COUNT = combined_calls %>% select(starts_with("MYC_"))  %>% sum(),
          MYC_Summary_CALL = if_else(MYC_CALLER_COUNT >= opt$count, 1, 0),
          CCND1_CALLER_COUNT = combined_calls %>% select(starts_with("CCND1_"))  %>% sum(),
          CCND1_Summary_CALL = if_else(CCND1_CALLER_COUNT >= opt$count, 1, 0),
          CCND2_CALLER_COUNT = combined_calls %>% select(starts_with("CCND2_"))  %>% sum(),
          CCND2_Summary_CALL = if_else(CCND2_CALLER_COUNT >= opt$count, 1, 0),
          CCND3_CALLER_COUNT = combined_calls %>% select(starts_with("CCND3_"))  %>% sum(),
          CCND3_Summary_CALL = if_else(CCND3_CALLER_COUNT >= opt$count, 1, 0)
  )
#Add IG source if matches across all caller
NSD2_IgSource <- c(if(exists("gammit_source")){ gammit_source$NSD2_IgSource },
                   if(exists("manta_source")){ manta_source$NSD2_IgSource },
                   if(exists("pair_source")){ pair_source$NSD2_IgSource })
CCND1_IgSource <- c(if(exists("gammit_source")){ gammit_source$CCND1_IgSource },
                    if(exists("manta_source")){ manta_source$CCND1_IgSource },
                    if(exists("pair_source")){ pair_source$CCND1_IgSource })
CCND2_IgSource <- c(if(exists("gammit_source")){ gammit_source$CCND2_IgSource },
                    if(exists("manta_source")){ manta_source$CCND2_IgSource },
                    if(exists("pair_source")){ pair_source$CCND2_IgSource })
CCND3_IgSource <- c(if(exists("gammit_source")){ gammit_source$CCND3_IgSource },
                    if(exists("manta_source")){ manta_source$CCND3_IgSource },
                    if(exists("pair_source")){ pair_source$CCND3_IgSource })
MYC_IgSource <- c(if(exists("gammit_source")){ gammit_source$MYC_IgSource },
                  if(exists("manta_source")){ manta_source$MYC_IgSource },
                  if(exists("pair_source")){ pair_source$MYC_IgSource })
MAF_IgSource <- c(if(exists("gammit_source")){ gammit_source$MAF_IgSource },
                  if(exists("manta_source")){ manta_source$MAF_IgSource },
                  if(exists("pair_source")){ pair_source$MAF_IgSource })
MAFA_IgSource <- c(if(exists("gammit_source")){ gammit_source$MAFA_IgSource },
                   if(exists("manta_source")){ manta_source$MAFA_IgSource },
                   if(exists("pair_source")){ pair_source$MAFA_IgSource })
MAFB_IgSource <- c(if(exists("gammit_source")){ gammit_source$MAFB_IgSource },
                   if(exists("manta_source")){ manta_source$MAFB_IgSource },
                   if(exists("pair_source")){ pair_source$MAFB_IgSource })

NSD2_IgSource <- NSD2_IgSource[NSD2_IgSource != 0]
CCND1_IgSource <- CCND1_IgSource[CCND1_IgSource != 0]
CCND2_IgSource <- CCND2_IgSource[CCND2_IgSource != 0]
CCND3_IgSource <- CCND3_IgSource[CCND3_IgSource != 0]
MYC_IgSource <- MYC_IgSource[MYC_IgSource != 0]
MAF_IgSource <- MAF_IgSource[MAF_IgSource != 0]
MAFA_IgSource <- MAFA_IgSource[MAFA_IgSource != 0]
MAFB_IgSource <- MAFB_IgSource[MAFB_IgSource != 0]

combined_calls$NSD2_IgSource = ifelse(length(unique(NSD2_IgSource))==1, unique(NSD2_IgSource),0)
combined_calls$CCND1_IgSource = ifelse(length(unique(CCND1_IgSource))==1, unique(CCND1_IgSource),0)
combined_calls$CCND2_IgSource = ifelse(length(unique(CCND2_IgSource))==1, unique(CCND2_IgSource),0)
combined_calls$CCND3_IgSource = ifelse(length(unique(CCND3_IgSource))==1, unique(CCND3_IgSource),0)
combined_calls$MYC_IgSource = ifelse(length(unique(MYC_IgSource))==1, unique(MYC_IgSource),0)
combined_calls$MAF_IgSource = ifelse(length(unique(MAF_IgSource))==1, unique(MAF_IgSource),0)
combined_calls$MAFA_IgSource = ifelse(length(unique(MAFA_IgSource))==1, unique(MAFA_IgSource),0)
combined_calls$MAFB_IgSource = ifelse(length(unique(MAFB_IgSource))==1, unique(MAFB_IgSource),0)

#combined_calls=combined_calls[,order(colnames(combined_calls), decreasing = TRUE)]

summary_table <- init_table %>%
    right_join(combined_calls) %>%
    replace(is.na(.), 0)

summary_table=summary_table[,order(colnames(summary_table), decreasing = TRUE)]

# combined_calls=combined_calls %>% relocate(Specimen)
write("Save results...\n", stderr())
write_tsv(summary_table, opt$outfile, append = FALSE, na="NA")
write("Done.\n", stderr())
