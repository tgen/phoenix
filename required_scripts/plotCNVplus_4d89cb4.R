#!/usr/bin/env Rscript --vanilla

# This code was copied and modified from the github 
# broadinstitute/gatk and tgen/tCoNut repositories.

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))

#### Preset Options For Readability ####

# Homo_sapiens_assembly38 primary contigs and lengths
contig_names <- c("chr1",
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
                 "chrY")
contig_lengths <- c(248956422,
                   242193529, 
                   198295559, 
                   190214555, 
                   181538259, 
                   170805979, 
                   159345973, 
                   145138636, 
                   138394717, 
                   133797422, 
                   135086622, 
                   133275309, 
                   114364328, 
                   107043718, 
                   101991189, 
                   90338345, 
                   83257441, 
                   80373285, 
                   58617616, 
                   64444167, 
                   46709983, 
                   50818468, 
                   156040895, 
                   57227415)

#### Argument parsing ####
option_list <- list(
  # Required
  make_option(c("--sample_name", "-s"), dest="sample_name", action="store", default = NA, type = 'character',
              help = "[Required] The title used for each plot"),
  make_option(c("--output_prefix", "-o"), dest="output_prefix", action="store", default = NA, type = 'character',
              help = "[Required] Prefix for each of the plots"),
  make_option(c("--plots_directory", "-O"), dest="plots_directory", action="store", default = NA, type = 'character',
              help = "[Required] Directory for each of the plots to be writen too"),
  make_option(c("--denoised_copy_ratios_file", "-r"), dest="denoised_copy_ratios_file", action="store",
              help = "[Required] .denoisedCR.tsv output from gatk DenoiseReadCounts"),
  make_option(c("--allelic_counts_file", "-a"), dest="allelic_counts_file", action="store",
              help = "[Required] .hets.tsv output from gatk ModelSegments"),
  make_option(c("--modeled_segments_file", "-m"), dest="modeled_segments_file", action="store",
              help = "[Required] .modelFinal.seg output from gatk ModelSegments"),
  
  # Optional With Defaults
  make_option(c("--contig_names_string", "-t"), dest="contig_names_string", action="store", type = 'character', default = NA,          #string with elements separated by "CONTIG_DELIMITER"
              help = "[Optional] \"CONTIG_DELIMITER\" seperated list of contig names to plot [default chr1CONTIG_DELIMITERchr2CONTIG_DELIMITER...]"),
  make_option(c("--contig_lengths_string", "-T"), dest="contig_lengths_string", action="store", type = 'character', default = NA,      #string with elements separated by "CONTIG_DELIMITER"
              help = "[Optional] \"CONTIG_DELIMITER\" seperated list of contig lengths in the same order as contig_names [default hg38 lengths 248956422CONTIG_DELIMITER242193529CONTIG_DELIMITER...]"),
  make_option(c("--re_center_CNA", "-e"), dest="re_center_CNA", action="store", default=FALSE,
              help = "[Optional] Re-Center all copy number values based on the distribution of heterozygous Log ratio [default %default]"),
  make_option(c("--re_centered_seg_directory", "-E"), dest="re_centered_seg_directory", action="store", default=NA, type = 'character',
              help = "[Optional] Directory to write the recentered seg file to. Required if --re_center_CNA is TRUE [default %default]"),
  make_option(c("--CNlossColor", "-c"), dest="CNlossColor", action="store", default="#018571",
              help = "[Optional] Color of points to use for CN losses [default %default]"),
  make_option(c("--CNgainColor", "-C"), dest="CNgainColor", action="store", default="#bf812d",
              help = "[Optional] Color of points to use for CN gains [default %default]"),
  make_option(c("--CNhetColor", "-n"), dest="CNhetColor", action="store", default="#00A8FF",
              help = "[Optional] Color of points to use for hets [default %default]"),
  make_option(c("--BAFcolor", "-b"), dest="BAFcolor", action="store", default="#018571",
              help = "[Optional] Color of points to use for ALT-Alleles-Frequency [default %default]"),
  make_option(c("--SEGcolor", "-g"), dest="SEGcolor", action="store", default="#00A8FF",
              help = "[Optional] Color of lines to denote segements [default %default]"),
  make_option(c("--CNlossLim", "-l"), dest="CNlossLim", action="store", default=-0.4150375,
              help = "[Optional] Log2 value used to denote CN losses [default %default]"),
  make_option(c("--CNgainLim", "-L"), dest="CNgainLim", action="store", default=0.3219281,
              help = "[Optional] Log2 value used to denote CN gains [default %default]"),
  make_option(c("--hetDPfilter", "-d"), dest="hetDPfilter", action="store", default=20,
              help = "[Optional] Depth requirement to filter hets [default %default]"),
  make_option(c("--hetAFlow", "-f"), dest="hetAFlow", action="store", default=0.4,
              help = "[Optional] Min allele frequency to keep hets [default %default]"),
  make_option(c("--hetAFhigh", "-F"), dest="hetAFhigh", action="store", default=0.6,
              help = "[Optional] Max allele frequency to keep hets [default %default]"),
  make_option(c("--hetMAFposteriorOffset", "-z"), dest="hetMAFposteriorOffset", action="store", default=0.01,
              help = "[Optional] Value to be added/subtracted to MINOR_ALLELE_FRACTION_POSTERIOR_10 and MINOR_ALLELE_FRACTION_POSTERIOR_90 from modeled_segments_file for filter hets. [default %default]"),
  make_option(c("--point_size", "-p"), dest="point_size", action="store", default=2,
              help = "[Optional] Size of points to be ploted [default %default]"),
  make_option(c("--lowerCNvalidatePeakOffset", "-j"), dest="lowerCNvalidatePeakOffset", action="store", default=0.125,
              help = "[Optional] The lower real copy number offset used for the selection of hets around a mode in the copy number density plot of hets used for validating the mode to be used for centering. [default %default]"),
  make_option(c("--UpperCNvalidatePeakOffset", "-J"), dest="UpperCNvalidatePeakOffset", action="store", default=0.125,
              help = "[Optional] The upper real copy number offset used for the selection of hets around a mode in the copy number density plot of hets used for validating the mode to be used for centering. [default %default]"),
  make_option(c("--lowerCNcenteringPeakOffset", "-i"), dest="lowerCNcenteringPeakOffset", action="store", default=0.25,
              help = "[Optional] The lower real copy number offset used for the selection of hets around a mode in the copy number density plot of hets used for centering of all segments. [default %default]"),
  make_option(c("--UpperCNcenteringPeakOffset", "-I"), dest="UpperCNcenteringPeakOffset", action="store", default=0.25,
              help = "[Optional] The upper real copy number offset used for the selection of hets around a mode in the copy number density plot of hets used for centering of all segments. [default %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Required
sample_name <- opt[["sample_name"]]
output_prefix <- opt[["output_prefix"]]
plots_directory <- opt[["plots_directory"]]

denoised_copy_ratios_file <- opt[["denoised_copy_ratios_file"]]
allelic_counts_file <- opt[["allelic_counts_file"]]
modeled_segments_file <- opt[["modeled_segments_file"]]

# Optional With Defaults

if (!is.na(opt[["contig_names_string"]])) {
  contig_names_string <- opt[["contig_names_string"]]
  contig_names <- as.list(strsplit(contig_names_string, "CONTIG_DELIMITER")[[1]])
}
if (!is.na(opt[["contig_lengths_string"]])) {
  contig_lengths_string <- opt[["contig_lengths_string"]]
  contig_lengths <- as.numeric(as.list(strsplit(contig_lengths_string, "CONTIG_DELIMITER")[[1]]))
}
re_center_CNA <- opt[["re_center_CNA"]]
re_centered_seg_directory <- opt[["re_centered_seg_directory"]]
CNgainColor <- opt[["CNgainColor"]]
CNlossColor <- opt[["CNlossColor"]]
CNhetColor <- opt[["CNhetColor"]]
BAFcolor <- opt[["BAFcolor"]]
SEGcolor <- opt[["SEGcolor"]]
CNgainLim <- opt[["CNgainLim"]]
CNlossLim <- opt[["CNlossLim"]]
hetDPfilter <- opt[["hetDPfilter"]]
hetAFlow <- opt[["hetAFlow"]]
hetAFhigh <- opt[["hetAFhigh"]]
hetMAFposteriorOffset <- opt[["hetMAFposteriorOffset"]]
point_size <- opt[["point_size"]]
lowerCNvalidatePeakOffset <- opt[["lowerCNvalidatePeakOffset"]]
UpperCNvalidatePeakOffset <- opt[["UpperCNvalidatePeakOffset"]]
lowerCNcenteringPeakOffset <- opt[["lowerCNcenteringPeakOffset"]]
UpperCNcenteringPeakOffset <- opt[["UpperCNcenteringPeakOffset"]]

#### Validation ####

if (is.na(opt[["sample_name"]])) {
  print("Error: --sample_name is required to generate plots.")
  quit(save="no", status=1, runLast=FALSE)
}
if (is.na(opt[["output_prefix"]])) {
  print("Error: --output_prefix is required to generate plots.")
  quit(save="no", status=1, runLast=FALSE)
}
if (!dir.exists(plots_directory)) {
  print("Error: --plots_directory is required to generate plots.")
  quit(save="no", status=1, runLast=FALSE)
}
if (re_center_CNA == TRUE) {
  if (!dir.exists(re_centered_seg_directory)) {
    print("Error: --re_centered_seg_directory is required when re_center_CNA is set to TRUE.")
    quit(save="no", status=1, runLast=FALSE)
  }
}
if (!file.exists(denoised_copy_ratios_file)) {
  print("Error: --denoised_copy_ratios_file is required to generate plots.")
  quit(save="no", status=1, runLast=FALSE)
}
if (!file.exists(allelic_counts_file)) {
  print("Error: --allelic_counts_file is required to generate plots.")
  quit(save="no", status=1, runLast=FALSE)
}
if (!file.exists(modeled_segments_file)) {
  print("Error: --modeled_segments_file is required to generate plots.")
  quit(save="no", status=1, runLast=FALSE)
}

#### Read TSV Function ####

ReadTSV <- function(tsv_file) {
  temp_file <- tempfile()
  system(sprintf('grep -v ^@ "%s" > %s', tsv_file, temp_file))
  return(suppressWarnings(fread(temp_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE, data.table=FALSE, showProgress=FALSE, verbose=FALSE)))
}

#### dLRs Function ####
dLRs <- function(x) {
  return(IQR(diff(na.omit(x))) / (4 * qnorm((1 + 0.5) / 2) / sqrt(2)))
}

#### Read In Files ####
denoised_copy_ratios <- ReadTSV(denoised_copy_ratios_file)
allelic_counts <- ReadTSV(allelic_counts_file)
modeled_segments <- ReadTSV(modeled_segments_file)

setDT(denoised_copy_ratios)
setDT(allelic_counts)
setDT(modeled_segments)

#### Process Files, lists, and perform recentering if specified  ####

contig_ends <- cumsum(contig_lengths)
contig_starts <- c(0, head(contig_ends, -1))

denoised_copy_ratios$MB <- denoised_copy_ratios$START/1000000

allelic_counts$DP <- allelic_counts$REF_COUNT + allelic_counts$ALT_COUNT
allelic_counts$AF <- allelic_counts$ALT_COUNT/allelic_counts$DP
allelic_counts$MB <- allelic_counts$POSITION/1000000

# Filter the hets in the allelic counts file by depth, and min and max allele frequency
hets <- allelic_counts[allelic_counts$DP >= hetDPfilter & allelic_counts$AF >= hetAFlow & allelic_counts$AF <= hetAFhigh,]

# Filter out hets that are not MINOR_ALLELE_FRACTION_POSTERIOR_10 < AF < MINOR_ALLELE_FRACTION_POSTERIOR_90
hets$MAF <- ifelse(hets$AF >= 0.5, 1-hets$AF, hets$AF)

hets <- hets[modeled_segments, on = .(CONTIG, POSITION >= START, POSITION <= END), MINOR_ALLELE_FRACTION_POSTERIOR_10 := i.MINOR_ALLELE_FRACTION_POSTERIOR_10 ][]
hets <- hets[modeled_segments, on = .(CONTIG, POSITION >= START, POSITION <= END), MINOR_ALLELE_FRACTION_POSTERIOR_90 := i.MINOR_ALLELE_FRACTION_POSTERIOR_90 ][]

hets <- hets[hets$MAF >= hets$MINOR_ALLELE_FRACTION_POSTERIOR_10 - hetMAFposteriorOffset & hets$MAF <= hets$MINOR_ALLELE_FRACTION_POSTERIOR_90 + hetMAFposteriorOffset,]

# Annotate each het with it's coresponding copy number ratio
hets <- hets[denoised_copy_ratios, on = .(CONTIG, POSITION >= START, POSITION <= END), LOG2_COPY_RATIO := i.LOG2_COPY_RATIO ][]
hets <- na.omit(hets)

setDF(hets)
setDF(denoised_copy_ratios)
setDF(allelic_counts)
setDF(modeled_segments)

# Re-Center copy number values
if (re_center_CNA == TRUE) {
  
  # Function for finding the modes of the distribution of copy ratios of the hets
  all_modes <- function(x) {
    mode <- NULL
    for ( CN in 2:(length(x)-1)) {
      if (( x[CN] > x[CN-1]) & (x[CN] > x[CN+1])) {
        mode <- c(mode,CN)
      }
    }
    if (length(mode) == 0) {
      mode <- 'Fail'
    }
    return(mode)
  }
  
  modes_values <- all_modes(density(hets$LOG2_COPY_RATIO, adjust = 1.4)$y)
  
  if (length(density(hets$LOG2_COPY_RATIO, adjust = 1.4)$x[modes_values]) == 1) {
    
    MAINPEAK <- density(hets$LOG2_COPY_RATIO, adjust = 1.4)$x[modes_values[1]]
    MAINPEAKlinear <- (2^MAINPEAK)*2
    
    CUTOFF <- log2((MAINPEAKlinear - lowerCNcenteringPeakOffset)/2)
    CUTOFF2 <- log2((MAINPEAKlinear + UpperCNcenteringPeakOffset)/2)
    
    POSITIONS <- hets[ hets$LOG2_COPY_RATIO > CUTOFF & hets$LOG2_COPY_RATIO < CUTOFF2 ,]
    
    AVE_CR <- mean(POSITIONS$LOG2_COPY_RATIO)
    
  } else {
    
    densum <- sum(density(hets$LOG2_COPY_RATIO, adjust = 1.4)$y)
    
    MAINPEAK <- density(hets$LOG2_COPY_RATIO, adjust = 1.4)$x[modes_values[1]]
    MAINPEAKlinear <- (2^MAINPEAK)*2
    
    CUTOFF <- log2((MAINPEAKlinear - lowerCNcenteringPeakOffset)/2)
    CUTOFF2 <- log2((MAINPEAKlinear + UpperCNcenteringPeakOffset)/2)
    CUTOFF3 <- log2((MAINPEAKlinear - lowerCNvalidatePeakOffset)/2)
    CUTOFF4 <- log2((MAINPEAKlinear + UpperCNvalidatePeakOffset)/2)
    
    densityList <- density(hets$LOG2_COPY_RATIO, adjust = 1.4)$x
    densityListFilt <- which(densityList >= CUTOFF3 & densityList <= CUTOFF4)
    densityListFiltSum <- sum(density(hets$LOG2_COPY_RATIO, adjust = 1.4)$y[densityListFilt])
    
    POSITIONS <- hets[ hets$LOG2_COPY_RATIO > CUTOFF & hets$LOG2_COPY_RATIO < CUTOFF2 ,]
    AVE_CR <- mean(POSITIONS$LOG2_COPY_RATIO)
    
    percent <- densityListFiltSum/densum
    
    if (percent > 0.025) {
      
    } else {
      
      INCVAR <- 1
      
      for (peak in density(hets$LOG2_COPY_RATIO, adjust = 1.4)$y[modes_values]) {
        
        MAINPEAK <- density(hets$LOG2_COPY_RATIO, adjust = 1.4)$x[modes_values[INCVAR]]
        MAINPEAKlinear <- (2^MAINPEAK)*2
        
        CUTOFF <- log2((MAINPEAKlinear - lowerCNcenteringPeakOffset)/2)
        CUTOFF2 <- log2((MAINPEAKlinear + UpperCNcenteringPeakOffset)/2)
        CUTOFF3 <- log2((MAINPEAKlinear - lowerCNvalidatePeakOffset)/2)
        CUTOFF4 <- log2((MAINPEAKlinear + UpperCNvalidatePeakOffset)/2)
        
        densityList <- density(hets$LOG2_COPY_RATIO, adjust = 1.4)$x
        densityListFilt <- which(densityList >= CUTOFF3 & densityList <= CUTOFF4)
        densityListFiltSum <- sum(density(hets$LOG2_COPY_RATIO, adjust = 1.4)$y[densityListFilt])
        
        POSITIONS <- hets[ hets$LOG2_COPY_RATIO > CUTOFF & hets$LOG2_COPY_RATIO < CUTOFF2 ,]
        
        AVE_CR <- mean(POSITIONS$LOG2_COPY_RATIO)
        
        percent <- densityListFiltSum/densum
        
        if (percent > 0.025 & INCVAR < length(density(hets$LOG2_COPY_RATIO, adjust = 1.4)$x[modes_values])) {
          break
        } else if (percent > 0.025 & INCVAR == length(density(hets$LOG2_COPY_RATIO, adjust = 1.4)$x[modes_values])) {
          break
        } else {
          INCVAR <- INCVAR + 1
        }
      }
    }
  }
  
  if (length(density(hets$LOG2_COPY_RATIO, adjust = 1.4)$x[modes_values]) == 1) {
    png(filename = paste0(plots_directory, "/", output_prefix, ".hets.density.png"), width = 648, height = 368, units = "mm", res = 300)
    plot(density(hets$LOG2_COPY_RATIO, adjust = 1.4)) + abline(v=CUTOFF) + abline(v=CUTOFF2) + abline(v=MAINPEAK)
    dev.off()
  } else {
    png(filename = paste0(plots_directory, "/", output_prefix, ".hets.density.png"), width = 648, height = 368, units = "mm", res = 300)
    plot(density(hets$LOG2_COPY_RATIO, adjust = 1.4)) + abline(v=CUTOFF) + abline(v=CUTOFF2) + abline(v=CUTOFF3, lty=2) + abline(v=CUTOFF4, lty=2) + abline(v=MAINPEAK)
    dev.off()
  }
  
  # Adjust all of the input tables with COPY RATIO based on the AVE_CR
  hets$LOG2_COPY_RATIO <- hets$LOG2_COPY_RATIO - AVE_CR
  denoised_copy_ratios$LOG2_COPY_RATIO <- denoised_copy_ratios$LOG2_COPY_RATIO - AVE_CR
  modeled_segments$LOG2_COPY_RATIO_POSTERIOR_10 <- modeled_segments$LOG2_COPY_RATIO_POSTERIOR_10 - AVE_CR
  modeled_segments$LOG2_COPY_RATIO_POSTERIOR_50 <- modeled_segments$LOG2_COPY_RATIO_POSTERIOR_50 - AVE_CR
  modeled_segments$LOG2_COPY_RATIO_POSTERIOR_90 <- modeled_segments$LOG2_COPY_RATIO_POSTERIOR_90 - AVE_CR
  
  # Write out the seg file
  SEG <- data.frame("Sample" = rep(sample_name, length(modeled_segments$CONTIG)),
                    "Chromosome" = modeled_segments$CONTIG, 
                    "Start" = modeled_segments$START,
                    "End" = modeled_segments$END, 
                    "Num_Probes" = modeled_segments$NUM_POINTS_COPY_RATIO, 
                    "Segment_Mean" = modeled_segments$LOG2_COPY_RATIO_POSTERIOR_50, 
                    stringsAsFactors = TRUE)
  
  write.table(SEG, file = paste0(re_centered_seg_directory, "/", output_prefix, ".re_centered.cr.igv.seg"), quote = FALSE, row.names = FALSE, sep = "\t")
  
}

# Create amplification and deletion dataframes
amp <- denoised_copy_ratios[denoised_copy_ratios$LOG2_COPY_RATIO > CNgainLim,]
del <- denoised_copy_ratios[denoised_copy_ratios$LOG2_COPY_RATIO < CNlossLim,]

# transform to linear copy ratio
denoised_copy_ratios[["COPY_RATIO"]] <- (2^denoised_copy_ratios[["LOG2_COPY_RATIO"]])*2

# determine copy-ratio midpoints
denoised_copy_ratios[["MIDDLE"]] <- round((denoised_copy_ratios[["START"]] + denoised_copy_ratios[["END"]]) / 2)

# Calculate dLRs and write to stats file
dLRs_value <- dLRs(denoised_copy_ratios[["LOG2_COPY_RATIO"]])
dLRs_value <- data.frame("Sample" = sample_name,
                         "dLRs" = dLRs_value,
                         stringsAsFactors = TRUE)
write.table(dLRs_value, file = paste0(re_centered_seg_directory, "/", output_prefix, ".dLRs.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")


#### CNV Multi CHR Plot Function ####

plotCNVgraph <- function(plots_directory, output_prefix, denoised_copy_ratios, amp, del, hets, CNgainColor=CNgainColor, CNlossColor=CNlossColor, CNhetColor=CNhetColor, plotHets=TRUE ) {
  
  #create PNG file
  if (plotHets) {
    fname <- paste0(plots_directory, "/", output_prefix, '_cna_withhets.png')
  } else {
    fname <- paste0(plots_directory, "/", output_prefix, '_cna.png')
  }
  
  png(file=fname,width=500*6*3,height=500*4*3,res=300)
  par(mfrow=c(6,4))
  
  for (i in contig_names){  
    tmpCNR <- subset(denoised_copy_ratios, CONTIG == i)
    tmpAMP <- subset(amp, CONTIG == i)
    tmpDEL <- subset(del, CONTIG == i)
    tmpHETS <- subset(hets, CONTIG == i)
    
    if ( length(rownames(tmpCNR)) == 0 ) {
      next()
    }
    plot(tmpCNR$MB, 
         tmpCNR$LOG2_COPY_RATIO,
         type="l",
         ylim=c(-2,2),
         col="black",
         main=i,
         ylab="Log2(T/N)",
         xlab="Physical Position (Mb)",
         xlim=c(floor(min(tmpCNR$MB)),ceiling(max(tmpCNR$MB))),
         xaxt="n",
         cex.lab=1.35,
         cex.main=1.5,
         cex.axis=1.15)
    points(tmpAMP$MB,
           tmpAMP$LOG2_COPY_RATIO,
           col=CNgainColor)
    points(tmpDEL$MB,
           tmpDEL$LOG2_COPY_RATIO,
           col=CNlossColor)
    if (plotHets) {
      points(tmpHETS$MB,
             tmpHETS$LOG2_COPY_RATIO,
             col=CNhetColor,
             pch=19,
             cex=1.5)
    }

    x <- seq(0,250,25)
    axis(1,at=x)
  }
  
  invisible(dev.off())
}

#### BAF Multi CHR Plot Function ####

plotBAFgraph <- function(plots_directory, output_prefix, allelic_counts, BAFcolor=BAFcolor) {
  
  #create PNG file
  png(file=paste0(plots_directory, "/", output_prefix,'_baf.png'),width=500*6*3,height=500*4*3, res=300)
  par(mfrow=c(6,4))
  
  for (i in contig_names){
    tmpBAF <- subset(allelic_counts, CONTIG == i)
    
    if (length(tmpBAF$CONTIG) == 0) {
      next()
    }
    
    plot(tmpBAF$MB,
         tmpBAF$AF,
         type="p",
         col=BAFcolor,
         bg=BAFcolor,
         pch=16,
         ylim=c(0,1),
         yaxt="n",
         main=i,
         ylab="Frequencies",
         xlab="Physical Position (Mb)",
         xlim=c(floor(min(tmpBAF$MB)),ceiling(max(tmpBAF$MB))),
         xaxt="n",
         lwd=3.5,
         cex.lab=1.5,
         cex.main=1.5
    )
    abline(h=0.5,lty=2)
    
    y <- seq(0,1,0.5)
    axis(2,at=y)
    x <- seq(0,250,25)
    axis(1,at=x)
  }
  invisible(dev.off())

}

#### CNV/BAF Combined Plot Functions ####

SetUpPlot <- function(sample_name, y.lab, y.min, y.max, contig_names, contig_starts, contig_ends, do_label_contigs) {
  num_contigs <- length(contig_names)
  contig_centers <- (contig_starts + contig_ends) / 2
  genome_length <- contig_ends[num_contigs]
  
  if (do_label_contigs) {
    # Bottom plot. Do not plot Main title.
    suppressWarnings(par(mar=c(2.1, 3.6, 0, 0), mgp=c(2, -0.2, -1.1)))
    
    plot(0, type="n", bty="n", xlim=c(0, genome_length), ylim=c(0, y.max), xlab="", ylab="", main=NULL, xaxt="n")
    
    mtext(side=2.2, line=1.5, y.lab, cex=0.75, las=FALSE, outer=FALSE)
    mtext(text=contig_names[1:num_contigs], side=1, line=ifelse(1:num_contigs %% 2 == 1, -0.45, 0.0),
          at = contig_centers[1:num_contigs], las=1, cex=par("cex.axis") * par("cex") * 0.7)
  } else {
    # Top plot. Do not plot x.lab or contig lables
    suppressWarnings(par(mar=c(0, 3.6, 3.6, 0), mgp=c(2, -0.2, -1.1)))
    
    plot(0, type="n", bty="n", xlim=c(0, genome_length), ylim=c(0, y.max), xlab="", ylab="", main=sample_name, xaxt="n")
    
    mtext(side=2.2, line=1.5, y.lab, cex=0.75, las=FALSE, outer=FALSE)
    
  }
  
  for (i in 1:num_contigs) {
    if (num_contigs > 1) {
      use.col <- ifelse(i %% 2 == 1, "grey90", "white")
    } else {
      #use.col = "white"
      use.col <- "grey90"
    }
    rect(xleft=contig_starts[i], ybottom=y.min, xright=contig_ends[i], ytop=y.max, col=use.col, border=NA)
  }
}

PlotCopyRatiosWithModeledSegments <- function(denoised_copy_ratios, modeled_segments, contig_names, contig_starts, point_size=0.2, CNgainColor=CNgainColor, CNlossColor=CNlossColor, CNgainLim=CNgainLim, CNlossLim=CNlossLim, chr_plots=TRUE, SEGlwd=2, totLength=sum(contig_lengths) ) {
  points_start_index <- 1
  for (s in seq_len(nrow(modeled_segments))) {
    #skip segments with no points
    num_points <- modeled_segments[s, "NUM_POINTS_COPY_RATIO"]
    if (num_points == 0) {
      next
    }
    
    contig <- modeled_segments[s, "CONTIG"]
    seg_start <- modeled_segments[s, "START"]
    seg_stop <- modeled_segments[s, "END"]
    offset <- contig_starts[match(contig, contig_names)]

    denoised_copy_ratios_range <- denoised_copy_ratios[denoised_copy_ratios$CONTIG == contig & denoised_copy_ratios$START >= seg_start & denoised_copy_ratios$END <= seg_stop, ]
    
    genomic_coordinates_max <- offset + denoised_copy_ratios_range[denoised_copy_ratios_range$LOG2_COPY_RATIO >= 1.584963, "MIDDLE"]
    genomic_coordinates_gain <- offset + denoised_copy_ratios_range[denoised_copy_ratios_range$LOG2_COPY_RATIO >= CNgainLim & denoised_copy_ratios_range$LOG2_COPY_RATIO < 1.584963, "MIDDLE"]
    genomic_coordinates_neutral <- offset + denoised_copy_ratios_range[denoised_copy_ratios_range$LOG2_COPY_RATIO < CNgainLim & denoised_copy_ratios_range$LOG2_COPY_RATIO >= CNlossLim, "MIDDLE"]
    genomic_coordinates_loss <- offset + denoised_copy_ratios_range[denoised_copy_ratios_range$LOG2_COPY_RATIO < CNlossLim , "MIDDLE"]
    
    denoised_copy_ratios_max <- denoised_copy_ratios_range[denoised_copy_ratios_range$LOG2_COPY_RATIO >= 1.584963, "COPY_RATIO"]
    denoised_copy_ratios_gain <- denoised_copy_ratios_range[denoised_copy_ratios_range$LOG2_COPY_RATIO >= CNgainLim & denoised_copy_ratios_range$LOG2_COPY_RATIO < 1.584963, "COPY_RATIO"]
    denoised_copy_ratios_neutral <- denoised_copy_ratios_range[denoised_copy_ratios_range$LOG2_COPY_RATIO < CNgainLim & denoised_copy_ratios_range$LOG2_COPY_RATIO >= CNlossLim, "COPY_RATIO"]
    denoised_copy_ratios_loss <- denoised_copy_ratios_range[denoised_copy_ratios_range$LOG2_COPY_RATIO < CNlossLim , "COPY_RATIO"]

    if (length(genomic_coordinates_max) > 0 ) {
      if( !chr_plots ) {
        denoised_copy_ratios_max[] <- 6
        points(x=genomic_coordinates_max, y=denoised_copy_ratios_max, col=CNgainColor, pch=".", cex=point_size + 3 )
      } else {
        points(x=genomic_coordinates_max, y=denoised_copy_ratios_max, col=CNgainColor, pch=".", cex=point_size )
      }
    }
    if (length(genomic_coordinates_gain) > 0 ) {
      points(x=genomic_coordinates_gain, y=denoised_copy_ratios_gain, col=CNgainColor, pch=".", cex=point_size)
    }
    if (length(genomic_coordinates_neutral) > 0 ) {
      points(x=genomic_coordinates_neutral, y=denoised_copy_ratios_neutral, col="black", pch=".", cex=point_size)
    }
    if (length(genomic_coordinates_loss) > 0 ) {
      points(x=genomic_coordinates_loss, y=denoised_copy_ratios_loss, col=CNlossColor, pch=".", cex=point_size)
    }
    
    points_start_index <- points_start_index + num_points
  }

  points_start_index <- 1
  for (s in seq_len(nrow(modeled_segments))) {
    #skip segments with no points
    num_points <- modeled_segments[s, "NUM_POINTS_COPY_RATIO"]
    if (num_points == 0) {
      next
    }

    contig <- modeled_segments[s, "CONTIG"]
    offset <- contig_starts[match(contig, contig_names)]
    segment_start <- offset + modeled_segments[s, "START"]
    segment_end <- offset + modeled_segments[s, "END"]

    copy_ratio_posterior_10 <- (2^modeled_segments[s, "LOG2_COPY_RATIO_POSTERIOR_10"])*2
    copy_ratio_posterior_50 <- (2^modeled_segments[s, "LOG2_COPY_RATIO_POSTERIOR_50"])*2
    copy_ratio_posterior_90 <- (2^modeled_segments[s, "LOG2_COPY_RATIO_POSTERIOR_90"])*2

    if( !chr_plots ) {
      if (copy_ratio_posterior_10 > 6 ) {
        copy_ratio_posterior_10 <- 6
      }
      if (copy_ratio_posterior_50 > 6 ) {
        copy_ratio_posterior_50 <- 6
      }
      if (copy_ratio_posterior_90 > 6 ) {
        copy_ratio_posterior_90 <- 6
      }
    }
    
    segments(x0=segment_start, y0=copy_ratio_posterior_50, x1=segment_end, y1=copy_ratio_posterior_50, col=SEGcolor, lwd=SEGlwd, lty=1)
    rect(xleft=segment_start, ybottom=copy_ratio_posterior_10, xright=segment_end, ytop=copy_ratio_posterior_90, col=SEGcolor, border=SEGcolor, lwd=1, lty=1)

    points_start_index <- points_start_index + num_points
  }

  segments(x0 = 0, y0=2, x1=totLength, y1=2, col="gray70")
}

PlotAlternateAlleleFractionsWithModeledSegments <- function(allelic_counts, modeled_segments, contig_names, contig_starts, point_size=0.4, BAFcolor=BAFcolor) {
  points_start_index <- 1
  for (s in seq_len(nrow(modeled_segments))) {
    #skip segments with no points
    num_points <- modeled_segments[s, "NUM_POINTS_ALLELE_FRACTION"]
    if (num_points == 0) {
      next
    }
    
    contig <- modeled_segments[s, "CONTIG"]
    offset <- contig_starts[match(contig, contig_names)]
    seg_start <- modeled_segments[s, "START"]
    seg_stop <- modeled_segments[s, "END"]
    segment_start <- offset + modeled_segments[s, "START"]
    segment_end <- offset + modeled_segments[s, "END"]

    genomic_coordinates <- offset + allelic_counts[allelic_counts$CONTIG == contig & allelic_counts$POSITION >= seg_start & allelic_counts$POSITION <= seg_stop, "POSITION"]

    ref_counts <- allelic_counts[allelic_counts$CONTIG == contig & allelic_counts$POSITION >= seg_start & allelic_counts$POSITION <= seg_stop, "REF_COUNT"]
    alt_counts <- allelic_counts[allelic_counts$CONTIG == contig & allelic_counts$POSITION >= seg_start & allelic_counts$POSITION <= seg_stop, "ALT_COUNT"]
    alternate_allele_fractions <- alt_counts / (alt_counts + ref_counts)
    
    points(x=genomic_coordinates, y=alternate_allele_fractions, col=BAFcolor, pch=".", cex=point_size)
    
    minor_allele_fraction_posterior_10 <- modeled_segments[s, "MINOR_ALLELE_FRACTION_POSTERIOR_10"]
    minor_allele_fraction_posterior_50 <- modeled_segments[s, "MINOR_ALLELE_FRACTION_POSTERIOR_50"]
    minor_allele_fraction_posterior_90 <- modeled_segments[s, "MINOR_ALLELE_FRACTION_POSTERIOR_90"]
    
    segments(x0=segment_start, y0=minor_allele_fraction_posterior_50, x1=segment_end, y1=minor_allele_fraction_posterior_50, col="black", lwd=2, lty=1)
    rect(xleft=segment_start, ybottom=minor_allele_fraction_posterior_10, xright=segment_end, ytop=minor_allele_fraction_posterior_90, lwd=1, lty=1)
    
    major_allele_fraction_posterior_10 <- 1 - minor_allele_fraction_posterior_10
    major_allele_fraction_posterior_50 <- 1 - minor_allele_fraction_posterior_50
    major_allele_fraction_posterior_90 <- 1 - minor_allele_fraction_posterior_90
    
    segments(x0=segment_start, y0=major_allele_fraction_posterior_50, x1=segment_end, y1=major_allele_fraction_posterior_50, col="black", lwd=2, lty=1)
    rect(xleft=segment_start, ybottom=major_allele_fraction_posterior_90, xright=segment_end, ytop=major_allele_fraction_posterior_10, lwd=1, lty=1)
    
    points_start_index <- points_start_index + num_points
  }
}

#### Build Segmented Plot Function ####

WriteModeledSegmentsPlot <- function(sample_name, allelic_counts, denoised_copy_ratios, modeled_segments, contig_names, contig_lengths, output_prefix, plots_directory) {
  
  output_file <- paste0(plots_directory, "/", output_prefix, '.png')
  num_plots <- 2
  
  png(output_file, 12, 3.5 * num_plots, units="in", type="cairo", res=300, bg="white")
  par(mfrow=c(num_plots, 1), cex=0.75, las=1)
  
  SetUpPlot(sample_name = sample_name, 
            y.lab = "Estimated Copy Number", 
            y.min = 0, 
            y.max = 6,
            contig_names = contig_names, 
            contig_starts = contig_starts, 
            contig_ends = contig_ends, 
            do_label_contigs = FALSE)
  
  PlotCopyRatiosWithModeledSegments(denoised_copy_ratios = denoised_copy_ratios, 
                                    modeled_segments = modeled_segments, 
                                    contig_names = contig_names, 
                                    contig_starts = contig_starts, 
                                    point_size = point_size, 
                                    CNgainColor = CNgainColor, 
                                    CNlossColor = CNlossColor, 
                                    CNgainLim = CNgainLim, 
                                    CNlossLim = CNlossLim,
                                    chr_plots = FALSE,
                                    SEGlwd=2)
  
  SetUpPlot(sample_name = sample_name, 
            y.lab = "ALT-Allele Frequency", 
            y.min = 0, 
            y.max = 1.0,
            contig_names = contig_names, 
            contig_starts = contig_starts, 
            contig_ends = contig_ends, 
            do_label_contigs = TRUE)
  
  PlotAlternateAlleleFractionsWithModeledSegments(allelic_counts = allelic_counts, 
                                                  modeled_segments = modeled_segments, 
                                                  contig_names = contig_names, 
                                                  contig_starts = contig_starts,
                                                  point_size = point_size,
                                                  BAFcolor = BAFcolor)
  
  dev.off()
  
  # Make individual contig plots
  for(i in seq_along(contig_names)) {

    contig_name <- contig_names[i]
    contig_length <- contig_lengths[i]
    
    output_file <- paste0(plots_directory, "/", output_prefix, '_', contig_name, '.png')
    
    modeled_segments_filt <- modeled_segments[modeled_segments$CONTIG == contig_name, ]
    denoised_copy_ratios_filt <- denoised_copy_ratios[denoised_copy_ratios$CONTIG == contig_name, ]
    allelic_counts_filt <- allelic_counts[allelic_counts$CONTIG == contig_name, ]
    
    if (length(denoised_copy_ratios_filt$COPY_RATIO) == 0) {
      next()
    }
    
    ymax <- ceiling(max(denoised_copy_ratios_filt$COPY_RATIO))
    
    if (ymax < 4) {
      ymax <- 4
    }
    
    png(output_file, 12, 3.5 * num_plots, units="in", type="cairo", res=300, bg="white")
    par(mfrow=c(num_plots, 1), cex=0.75, las=1)
    
    # Delete the x axis title contig
    
    SetUpPlot(sample_name = sample_name, 
              y.lab = "Estimated Copy Number", 
              y.min = 0, 
              y.max = ymax, 
              contig_names = contig_name, 
              contig_starts = 0, 
              contig_ends = contig_length, 
              do_label_contigs = FALSE)
    
    PlotCopyRatiosWithModeledSegments(denoised_copy_ratios = denoised_copy_ratios_filt, 
                                      modeled_segments = modeled_segments_filt, 
                                      contig_names = contig_name, 
                                      contig_starts = 0, 
                                      point_size = point_size + 2,
                                      CNgainColor = CNgainColor, 
                                      CNlossColor = CNlossColor, 
                                      CNgainLim = CNgainLim, 
                                      CNlossLim = CNlossLim,
                                      chr_plots = TRUE,
                                      SEGlwd=4,
                                      totLength=contig_length)
    
    SetUpPlot(sample_name = sample_name, 
              y.lab = "ALT-Allele Frequency",
              y.min = 0, 
              y.max = 1.0,
              contig_names = contig_name, 
              contig_starts = 0, 
              contig_ends = contig_length, 
              do_label_contigs = TRUE)
    
    PlotAlternateAlleleFractionsWithModeledSegments(allelic_counts = allelic_counts_filt, 
                                                    modeled_segments = modeled_segments_filt, 
                                                    contig_names = contig_name, 
                                                    contig_starts = 0,
                                                    point_size = point_size + 2,
                                                    BAFcolor = BAFcolor)
    
    dev.off()
  }
}

#### Call plot functions ####

plotCNVgraph(plots_directory = plots_directory, 
             output_prefix = output_prefix,
             denoised_copy_ratios = denoised_copy_ratios,
             amp = amp, 
             del = del, 
             hets = hets, 
             CNgainColor = CNgainColor, 
             CNlossColor = CNlossColor, 
             CNhetColor = CNhetColor, 
             plotHets = TRUE)

plotCNVgraph(plots_directory = plots_directory, 
             output_prefix = output_prefix,
             denoised_copy_ratios = denoised_copy_ratios,
             amp = amp, 
             del = del, 
             hets = hets, 
             CNgainColor = CNgainColor, 
             CNlossColor = CNlossColor, 
             CNhetColor = CNhetColor, 
             plotHets = FALSE)

plotBAFgraph(plots_directory = plots_directory, 
             output_prefix = output_prefix,
             allelic_counts = allelic_counts,
             BAFcolor = BAFcolor)

WriteModeledSegmentsPlot(sample_name = sample_name,
                         allelic_counts = allelic_counts,
                         denoised_copy_ratios = denoised_copy_ratios,
                         modeled_segments = modeled_segments,
                         contig_names = contig_names,
                         contig_lengths = contig_lengths,
                         output_prefix = output_prefix,
                         plots_directory = plots_directory)

