#!/usr/bin/env Rscript --vanilla

# Load required modules
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(optparse))

# Define Options
option_list = list(
  make_option(c("-i", "--samtoolsStatsFile"),
              type="character",
              default=NULL,
              help="samtools stats summary numbers table",
              metavar="filename"),
  make_option(c("-d", "--samtoolsDuplicatesFile"),
              type="character",
              default=NULL,
              help="samtools markdup stats table",
              metavar="filename"),
  make_option(c("-b", "--bam"),
              type="character",
              default="NA",
              help="BAM filename (For Graph Title) [NA]",
              metavar="file.bam"),
  make_option(c("-s", "--sample"),
              type="character",
              default="NA",
              help="Unique Sample Identifier (ie. Bam spec: RG_SM) [NA]",
              metavar="RG_SM"),
  make_option(c("-l", "--library"),
              type="character",
              default="NA",
              help="Unique Library Identifier(ie. Bam spec: RG_LB) [NA]",
              metavar="RG_LB"),
  make_option(c("-r", "--readgroup"),
              type="character",
              default=NA,
              help="Unique Read Group Identifier(ie. Bam spec: RG_ID) [NA]",
              metavar="RG_ID"),
  make_option(c("-f", "--readformat"),
              type="character",
              default="PairedEnd",
              help="Sequencing Read Format (PairedEnd or SingleEnd) [PairedEnd]",
              metavar="Format PE or SE")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

####################################
## Validate Required Options are Provided
####################################

if (is.null(opt$bam)){
  print_help(opt_parser)
  stop("You must provide the source BAM filename to -b/--bam", call.=FALSE)
}

####################################
## Define Functions
####################################

# Summarize base distiribution by cycle and read
baseDistribution_summary <- function(input, bam, rgsm, rglb, rgid) {
  # Extract first read data table
  fbc <- grep("^FBC",input, value=TRUE)
  fbc <- separate(tibble(fbc),
                         col=1,
                         into=c("ID", "Cycle", "A_Bases", "C_Bases", "G_Bases", "T_Bases", "N_Bases", "NoBases"),
                         sep="\t") %>%
    type_convert() %>%
    select(-ID, -NoBases)

  # Pivot for graphing
  fbc_long <- fbc %>% pivot_longer(-Cycle, names_to = "Bases", values_to = "Percentage") %>% add_column(Source = "First Read - R1")

  # Extract second read data table
  lbc <- grep("^LBC",input, value=TRUE)
  lbc <- separate(tibble(lbc),
                  col=1,
                  into=c("ID", "Cycle", "A_Bases", "C_Bases", "G_Bases", "T_Bases", "N_Bases", "NoBases"),
                  sep="\t") %>%
    type_convert() %>%
    select(-ID, -NoBases)

  # Pivot for graphing
  lbc_long <- lbc %>% pivot_longer(-Cycle, names_to = "Bases", values_to = "Percentage") %>% add_column(Source = "Last Read - R2")

  # Join tables
  bc_long <- bind_rows(fbc_long, lbc_long)

  # Graph
  ggplot(bc_long, aes(Cycle, Percentage, color = Bases)) +
    geom_line(stat = "identity", position = position_jitterdodge(jitter.height = 0.1) ) +
    scale_y_continuous() +
    scale_color_manual(values = c("Green", "Blue", "Gold", "Black", "Red"),
                       name = "Base",
                       breaks = c("A_Bases", "C_Bases", "G_Bases", "T_Bases", "N_Bases"),
                       labels = c("A", "C", "G", "T", "N"), ) +
    ggtitle(label = bam) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="top") +
    facet_grid(. ~ Source)
  ggsave(paste(bam, "samtools_baseDistribution_linePlot.png", sep = "_"))

}

baseQuality_summary <- function(input, bam, rgsm, rglb, rgid, format) {
  ## Import first fragment base qualities and distirbution
  ffbq <- grep("^FFQ",input, value=TRUE)
  ffbq <- separate(tibble(ffbq),
                   col=1,
                   into = c("ID","Cycle", seq(from = 0, to = 40, by = 1)),
                   sep="\t") %>%
    type_convert() %>%
    select(-ID)
  # Convert to longformat to make graphing easier and add column to indicate read of origin
  ffbq_long <- ffbq %>%
    pivot_longer(-Cycle, names_to = "BaseQuality", values_to = "Bases") %>%
    add_column(Read = "First") %>% type_convert()
  # Create new column with the total quality value
  ffbq_long <- ffbq_long %>%
    select(Read, Cycle, BaseQuality, Bases) %>%
    mutate(Total_BaseQuality = BaseQuality * Bases)

  if (format == "PairedEnd") {
    ## Import second fragment base qualities and distirbution
    sfbq <- grep("^LFQ",input, value=TRUE)
    sfbq <- separate(tibble(sfbq),
                     col=1,
                     into = c("ID","Cycle", seq(from = 0, to = 40, by = 1)),
                     sep="\t") %>%
      type_convert() %>%
      select(-ID)
    # Convert to longformat to make graphing easier and add column to indicate read of origin
    sfbq_long <- sfbq %>%
      pivot_longer(-Cycle, names_to = "BaseQuality", values_to = "Bases") %>%
      add_column(Read = "Second") %>% type_convert()
    # Create new column with the total quality value
    sfbq_long <- sfbq_long %>%
      select(Read, Cycle, BaseQuality, Bases) %>%
      mutate(Total_BaseQuality = BaseQuality * Bases)

    # Join the two tables together, remove lines with Bases = NA
    bq_table <- bind_rows(ffbq_long, sfbq_long) %>% filter(!is.na(Bases))
  } else {
    bq_table <- ffbq_long %>% filter(!is.na(Bases))
  }

  # Summarize the overall yield
  bq_summary <- bq_table %>%
    mutate(Q20_Base_Count = if_else(BaseQuality >= 20, Bases, 0)) %>%
    mutate(Q30_Base_Count = if_else(BaseQuality >= 30, Bases, 0)) %>%
    mutate(Total_Quality_Yield = BaseQuality * Bases) %>%
    summarise(Total_Bases = sum(Bases),
              Q20_Bases = sum(Q20_Base_Count),
              Q30_Bases = sum(Q30_Base_Count),
              Q20_Equivalent_Yield = sum(Total_Quality_Yield) / 20) %>%
    add_column(Sample = rgsm) %>%
    add_column(Library = rglb) %>%
    add_column(Read_Group = rgid) %>%
    add_column(BAM_File = bam)

  # Save the data
  write_tsv(bq_summary,paste(bam, "samtools_baseQualityYield_summary.tsv", sep = "_"))

  # Plot #1 - Pecentage of Bases at each Quality Value (x=Q-value y=percent) facet by read
  plot1_data <- bq_table %>%
    group_by(Read, BaseQuality) %>%
    summarise(TotalBases = sum(Bases))

  # Get sum of total bases for each read
  if (format == "PairedEnd") {
    r1_base_sum <- plot1_data %>% filter(Read == "First") %>% summarise(Bases_Sum = sum(TotalBases)) %>% pull(var = Bases_Sum)
    r2_base_sum <- plot1_data %>% filter(Read == "Second") %>% summarise(Bases_Sum = sum(TotalBases)) %>% pull(var = Bases_Sum)
    # Calculate the percentage of each base quality by cycle
    plot1_data <- plot1_data %>%
      mutate(Percent_Bases = case_when(Read == "First" ~ TotalBases / r1_base_sum,
                                       Read == "Second" ~ TotalBases / r2_base_sum)) %>%
      add_column(Sample = rgsm) %>%
      add_column(Library = rglb) %>%
      add_column(Read_Group = rgid) %>%
      add_column(BAM_File = bam)
  } else {
    r1_base_sum <- plot1_data %>% filter(Read == "First") %>% summarise(Bases_Sum = sum(TotalBases)) %>% pull(var = Bases_Sum)
    # Calculate the percentage of each base quality by cycle
    plot1_data <- plot1_data %>%
      mutate(Percent_Bases = case_when(Read == "First" ~ TotalBases / r1_base_sum)) %>%
      add_column(Sample = rgsm) %>%
      add_column(Library = rglb) %>%
      add_column(Read_Group = rgid) %>%
      add_column(BAM_File = bam)
  }

  # Save the data
  write_tsv(plot1_data,paste(bam, "samtools_baseQualityDistribution_histogram.tsv", sep = "_"))

  # Graph
  ggplot(plot1_data, aes(BaseQuality, Percent_Bases, fill = Read)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_y_continuous(breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
                       name = "Percentage of Bases") +
    scale_fill_manual(values = c("#d8b365", "#5ab4ac"),
                       name = "Sequencing Read",
                       labels = c("First - R1","Last - R2")) +
    ggtitle(bam) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="top")
  ggsave(paste(bam, "samtools_baseQualityDistribution_histogram.png", sep = "_"))



  # Plot #2 - Mean Quality Value at each Cycle (x=cycle y=meanQualityByCycle) facet by read
  plot2_data <- bq_table %>%
    mutate(Total_Quality = BaseQuality * Bases) %>%
    group_by(Read, Cycle) %>%
    summarise(MeanBaseQuality = sum(Total_Quality) / sum(Bases)) %>%
    add_column(Sample = rgsm) %>%
    add_column(Library = rglb) %>%
    add_column(Read_Group = rgid) %>%
    add_column(BAM_File = bam)

  # Save the data
  write_tsv(plot2_data,paste(bam, "samtools_meanBaseQualityByCycle_histogram.tsv", sep = "_"))

  # Graph
  ggplot(plot2_data, aes(Cycle, MeanBaseQuality, color = Read)) +
    geom_line(stat = "identity") +
    scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40),
                       limits = c(20,40)) +
    scale_color_manual(values = c("#d8b365", "#5ab4ac"),
                       name = "Sequencing Read",
                       labels = c("First - R1","Last - R2")) +
    ggtitle(bam) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="top")
  ggsave(paste(bam, "samtools_meanBaseQualityByCycle_lineplot.png", sep = "_"))
}

#### Coverage Summary

coverage_summary <- function(input, bam, rgsm, rglb, rgid) {
  # Import data table
  coverage <- grep("^COV",input, value=TRUE)
  coverage <- separate(tibble(coverage),
                       col=1,
                       into=c("ID","Range", "Depth", "Bases"),
                       sep="\t") %>%
    type_convert() %>%
    select(-ID)

  # Add meta-data
  coverage <- coverage %>%
    add_column(Sample = rgsm) %>%
    add_column(Library = rglb) %>%
    add_column(Read_Group = rgid) %>%
    add_column(BAM_File = bam)

  # Save histogram table
  write_tsv(coverage, paste(bam, "samtools_coverage_histogram.tsv", sep = "_"))
  #coverage <- read_tsv(paste(bam, "samtools_coverage_histogram.tsv", sep = "_"))

  # Add column for total depth
  coverage <- coverage %>%
    mutate(Last_Ordered_Base = cumsum(Bases)) %>%
    mutate(First_Ordered_Base = Last_Ordered_Base - Bases + 1) %>%
    mutate(Pct_Total_Bases = Bases / sum(Bases)) %>%
    mutate(Total_Depth = Depth * Bases) %>%
    mutate(Bases_2x = case_when(Depth >= 2 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_5x = case_when(Depth >= 5 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_10x = case_when(Depth >= 10 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_15x = case_when(Depth >= 15 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_20x = case_when(Depth >= 20 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_25x = case_when(Depth >= 25 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_30x = case_when(Depth >= 30 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_35x = case_when(Depth >= 35 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_40x = case_when(Depth >= 40 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_50x = case_when(Depth >= 50 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_60x = case_when(Depth >= 60 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_70x = case_when(Depth >= 70 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_80x = case_when(Depth >= 80 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_90x = case_when(Depth >= 90 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_100x = case_when(Depth >= 100 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_110x = case_when(Depth >= 110 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_120x = case_when(Depth >= 120 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_130x = case_when(Depth >= 130 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_140x = case_when(Depth >= 140 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_150x = case_when(Depth >= 150 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_200x = case_when(Depth >= 200 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_300x = case_when(Depth >= 300 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_400x = case_when(Depth >= 400 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_500x = case_when(Depth >= 500 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_600x = case_when(Depth >= 600 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_700x = case_when(Depth >= 700 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_800x = case_when(Depth >= 800 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_900x = case_when(Depth >= 900 ~ Bases, TRUE ~ 0)) %>%
    mutate(Bases_1000x = case_when(Depth >= 1000 ~ Bases, TRUE ~ 0))

  # Calculate MEAN, for SD calculation
  fast_mean <- coverage %>%
    summarise(Bases_Tested = sum(Bases),
              Tested_Total_Depth = sum(Total_Depth)) %>%
    mutate(MEAN = Tested_Total_Depth / Bases_Tested) %>%
    pull(var = MEAN)

  # Add columne to coverage table SD_lineNumerator = ((Depth - MEAN)^2) * Bases
  coverage <- coverage %>% mutate(SD_lineNumerator = ((Depth - fast_mean)^2) * Bases)
  ## Then to calculate Variance = (1 / Bases_Tested) * SD_lineNumerator
  ## Then to calcualte SD = sqrt(Variance)

  # Summarize
  summary <- coverage %>%
    summarise(Bases_Tested = sum(Bases),
              Tested_Total_Depth = sum(Total_Depth),
              Bases_2x = sum(Bases_2x),
              Bases_5x = sum(Bases_5x),
              Bases_10x = sum(Bases_10x),
              Bases_15x = sum(Bases_15x),
              Bases_20x = sum(Bases_20x),
              Bases_25x = sum(Bases_25x),
              Bases_30x = sum(Bases_30x),
              Bases_35x = sum(Bases_35x),
              Bases_40x = sum(Bases_40x),
              Bases_50x = sum(Bases_50x),
              Bases_60x = sum(Bases_60x),
              Bases_70x = sum(Bases_70x),
              Bases_80x = sum(Bases_80x),
              Bases_90x = sum(Bases_90x),
              Bases_100x = sum(Bases_100x),
              Bases_110x = sum(Bases_110x),
              Bases_120x = sum(Bases_120x),
              Bases_130x = sum(Bases_130x),
              Bases_140x = sum(Bases_140x),
              Bases_150x = sum(Bases_150x),
              Bases_200x = sum(Bases_200x),
              Bases_300x = sum(Bases_300x),
              Bases_400x = sum(Bases_400x),
              Bases_500x = sum(Bases_500x),
              Bases_600x = sum(Bases_600x),
              Bases_700x = sum(Bases_700x),
              Bases_800x = sum(Bases_800x),
              Bases_900x = sum(Bases_900x),
              Bases_1000x = sum(Bases_1000x),
              SD_lineNumerator = sum(SD_lineNumerator))

  # Get position information for coverge percentiles
  median_position <- summary %>% pull(var = Bases_Tested) * 0.5
  median_value <- coverage %>% filter(First_Ordered_Base <= median_position) %>% filter(Last_Ordered_Base > median_position) %>% pull(var = Depth)
  percentile1_position <- summary %>% pull(var = Bases_Tested) * 0.01
  percentile1_value <- coverage %>% filter(First_Ordered_Base <= percentile1_position) %>% filter(Last_Ordered_Base > percentile1_position) %>% pull(var = Depth)
  percentile5_position <- summary %>% pull(var = Bases_Tested) * 0.05
  percentile5_value <- coverage %>% filter(First_Ordered_Base <= percentile5_position) %>% filter(Last_Ordered_Base > percentile5_position) %>% pull(var = Depth)
  percentile10_position <- summary %>% pull(var = Bases_Tested) * 0.10
  percentile10_value <- coverage %>% filter(First_Ordered_Base <= percentile10_position) %>% filter(Last_Ordered_Base > percentile10_position) %>% pull(var = Depth)
  percentile25_position <- summary %>% pull(var = Bases_Tested) * 0.25
  percentile25_value <- coverage %>% filter(First_Ordered_Base <= percentile25_position) %>% filter(Last_Ordered_Base > percentile25_position) %>% pull(var = Depth)
  percentile75_position <- summary %>% pull(var = Bases_Tested) * 0.75
  percentile75_value <- coverage %>% filter(First_Ordered_Base <= percentile75_position) %>% filter(Last_Ordered_Base > percentile75_position) %>% pull(var = Depth)
  percentile90_position <- summary %>% pull(var = Bases_Tested) * 0.90
  percentile90_value <- coverage %>% filter(First_Ordered_Base <= percentile90_position) %>% filter(Last_Ordered_Base > percentile90_position) %>% pull(var = Depth)
  percentile95_position <- summary %>% pull(var = Bases_Tested) * 0.95
  percentile95_value <- coverage %>% filter(First_Ordered_Base <= percentile95_position) %>% filter(Last_Ordered_Base > percentile95_position) %>% pull(var = Depth)
  percentile99_position <- summary %>% pull(var = Bases_Tested) * 0.99
  percentile99_value <- coverage %>% filter(First_Ordered_Base <= percentile99_position) %>% filter(Last_Ordered_Base > percentile99_position) %>% pull(var = Depth)

  # Determine the MAD
  mad_cov <- coverage %>%
    select(Depth, Bases) %>%
    mutate(Abs_Deviation = abs(Depth - median_value)) %>%
    arrange(desc(Abs_Deviation)) %>%
    mutate(Last_Ordered_Base = cumsum(Bases)) %>%
    mutate(First_Ordered_Base = Last_Ordered_Base - Bases + 1) %>%
    filter(First_Ordered_Base <= median_position) %>%
    filter(Last_Ordered_Base > median_position) %>%
    pull(var = Abs_Deviation)

  # Add Summary Calculations
  summary <- summary %>%
    mutate(Mean_Coverage = Tested_Total_Depth / Bases_Tested) %>%
    mutate(Variance = (1 / Bases_Tested) * SD_lineNumerator) %>%
    mutate(SD = sqrt(Variance)) %>%
    mutate(Median = median_value) %>%
    mutate(MAD = mad_cov) %>%
    mutate(Percentile_1 = percentile1_value) %>%
    mutate(Percentile_5 = percentile5_value) %>%
    mutate(Percentile_10 = percentile10_value) %>%
    mutate(Percentile_25 = percentile25_value) %>%
    mutate(Percentile_75 = percentile75_value) %>%
    mutate(Percentile_90 = percentile90_value) %>%
    mutate(Percentile_95 = percentile95_value) %>%
    mutate(Percentile_99 = percentile99_value) %>%
    mutate(IQR = Percentile_75 - Percentile_25) %>%
    mutate(Pct_2x = Bases_2x / Bases_Tested) %>%
    mutate(Pct_5x = Bases_5x / Bases_Tested) %>%
    mutate(Pct_10x = Bases_10x / Bases_Tested) %>%
    mutate(Pct_15x = Bases_15x / Bases_Tested) %>%
    mutate(Pct_20x = Bases_20x / Bases_Tested) %>%
    mutate(Pct_25x = Bases_25x / Bases_Tested) %>%
    mutate(Pct_30x = Bases_30x / Bases_Tested) %>%
    mutate(Pct_35x = Bases_35x / Bases_Tested) %>%
    mutate(Pct_40x = Bases_40x / Bases_Tested) %>%
    mutate(Pct_50x = Bases_50x / Bases_Tested) %>%
    mutate(Pct_60x = Bases_60x / Bases_Tested) %>%
    mutate(Pct_70x = Bases_70x / Bases_Tested) %>%
    mutate(Pct_80x = Bases_80x / Bases_Tested) %>%
    mutate(Pct_90x = Bases_90x / Bases_Tested) %>%
    mutate(Pct_100x = Bases_100x / Bases_Tested) %>%
    mutate(Pct_110x = Bases_110x / Bases_Tested) %>%
    mutate(Pct_120x = Bases_120x / Bases_Tested) %>%
    mutate(Pct_130x = Bases_130x / Bases_Tested) %>%
    mutate(Pct_140x = Bases_140x / Bases_Tested) %>%
    mutate(Pct_150x = Bases_150x / Bases_Tested) %>%
    mutate(Pct_200x = Bases_200x / Bases_Tested) %>%
    mutate(Pct_300x = Bases_300x / Bases_Tested) %>%
    mutate(Pct_400x = Bases_400x / Bases_Tested) %>%
    mutate(Pct_500x = Bases_500x / Bases_Tested) %>%
    mutate(Pct_600x = Bases_600x / Bases_Tested) %>%
    mutate(Pct_700x = Bases_700x / Bases_Tested) %>%
    mutate(Pct_800x = Bases_800x / Bases_Tested) %>%
    mutate(Pct_900x = Bases_900x / Bases_Tested) %>%
    mutate(Pct_1000x = Bases_1000x / Bases_Tested) %>%
    add_column(Sample = rgsm) %>%
    add_column(Library = rglb) %>%
    add_column(Read_Group = rgid) %>%
    add_column(BAM_File = bam)

  # TO ADD: HET_SNP_SENSITIVITY, HET_SNP_Q
  # Can't do equivalent as done on the fly - PCT_EXC_ADAPTER PCT_EXC_MAPQ, PCT_EXC_DUPE, PCT_EXC_UNPAIRED, PCT_EXC_BASEQ, PCT_EXC_OVERLAP, PCT_EXC_CAPPED, PCT_EXC_TOTAL

  # Save summary table
  write_tsv(summary, paste(bam, "samtools_coverage_summary.tsv", sep = "_"))

  # Define ploting positions for text blob
  max_depth <- coverage %>% filter(Depth == max(Depth)) %>% filter(Bases == min(Bases)) %>% pull(var = Depth)
  max_pct_bases <- coverage %>% filter(Pct_Total_Bases == max(Pct_Total_Bases)) %>% pull(var = Pct_Total_Bases)

  mean_cov <- pull(summary, var = Mean_Coverage)
  sd_cov <- pull(summary, var = SD) %>% round(digits = 3)
  pct10x_cov <- pull(summary, var = Pct_10x) %>% round(digits = 3) * 100
  pct20x_cov <- pull(summary, var = Pct_20x) %>% round(digits = 3) * 100
  pct30x_cov <- pull(summary, var = Pct_30x) %>% round(digits = 3) * 100
  pct100x_cov <- pull(summary, var = Pct_100x) %>% round(digits = 3) * 100
  pct150x_cov <- pull(summary, var = Pct_150x) %>% round(digits = 3) * 100
  pct200x_cov <- pull(summary, var = Pct_200x) %>% round(digits = 3) * 100
  coverage_label <- paste(paste("Mean(SD): ", round(mean_cov, digits = 2), " +/- ", sd_cov, sep = ""),
                          paste("Median(MAD): ", median_value, " +/- ", mad_cov, sep = ""),
                          paste("10x: ", pct10x_cov, "%", sep = ""),
                          paste("20x: ", pct20x_cov, "%", sep = ""),
                          paste("30x: ", pct30x_cov, "%", sep = ""),
                          paste("100x: ", pct100x_cov, "%", sep = ""),
                          paste("150x: ", pct150x_cov, "%", sep = ""),
                          paste("200x: ", pct200x_cov, "%", sep = ""),
                          sep = "\n")

  # Plot
  ggplot(coverage, aes(Depth, Pct_Total_Bases)) +
    geom_density(stat="identity", fill = "#d8b365") +
    geom_vline(xintercept = mean_cov, color = "#5ab4ac", linetype = "dashed") +
    annotate(geom = "text", x=mean_cov * 3.25, y=max_pct_bases * 0.85, label = coverage_label) +
    scale_x_continuous(name = "Coverage Depth", limits = c(1,mean_cov * 4)) +
    scale_y_continuous(name = "Percent Bases", labels = scales::percent) +
    ggtitle(bam) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste(bam, "samtools_coverage_histogram.png", sep = "_"))

  # Return the mean insert size to pass to the GC effect on coverage function
  return(median_value)
}

#### GC Depth Summary
gcdepth_summary <- function(input, median_coverage, bam, rgsm, rglb, rgid) {

  ## Extract data table for relationship between GC percentage and coverage
  gcDepth <- grep("^GCD",input, value=TRUE)
  gcDepth <- separate(tibble(gcDepth),
                         col=1,
                         into=c("ID", "PercentGC", "UniquePercentiles", "Percentile10", "Percentile25", "Percentile50", "Percentile75", "Percentile90"),
                         sep="\t") %>%
    type_convert() %>%
    select(-ID)

  # Add meta-data to insert size table (histogram)
  gcDepth <- gcDepth %>%
    add_column(Sample = rgsm) %>%
    add_column(Library = rglb) %>%
    add_column(Read_Group = rgid) %>%
    add_column(BAM_File = bam)

  # Add percent (NEED TO DO FIRST ON CORRECT)
  gcDepth <- gcDepth %>%
    mutate(Pct_Genome = if_else(row_number()==1, UniquePercentiles, UniquePercentiles - lag(UniquePercentiles)))

  # Set colors for forced legend
  colors <- c("GC Bin Frequency" = "#d8b365", "Median Coverage (+/-IQR)" = "#5ab4ac")
  # Calculate a scale factor for the GC bins
  max_bin_freq <- max(gcDepth$Pct_Genome)
  scale_factor <- (median_coverage / 2 ) / max_bin_freq
  # Plot the frequency of the genome territory %GC bins and the median +/- ICR
  ggplot(gcDepth, aes(PercentGC, Percentile50)) +
    geom_hline(yintercept = median_coverage, linetype = "longdash", color = "#a6611a") +
    geom_area(aes(PercentGC, Pct_Genome * scale_factor), stat = "identity", fill = "#d8b365") +
    geom_point(aes(PercentGC, Pct_Genome * scale_factor, color = "GC Bin Frequency"), fill = "#d8b365") +
    geom_pointrange(aes(ymin = Percentile25, ymax = Percentile75, color = "Median Coverage (+/-IQR)"), fill = "#5ab4ac") +
    annotate(geom = "text", x = 5, y= median_coverage + (median_coverage*0.05), label = "Median") +
    scale_y_continuous(name = "Coverage", limits = c(0,median_coverage * 2), sec.axis = sec_axis(~ . / scale_factor, name = "GC Bin Frequency")) +
    scale_x_continuous(name = "Percent GC") +
    labs(color = "Legend") +
    scale_color_manual(values = colors) +
    ggtitle(bam) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.y.right = element_text(color = "#d8b365"),
          axis.text.y.right = element_text(color = "#d8b365"),
          axis.ticks.y.right = element_line(color = "#d8b365"),
          legend.position = "top")
  ggsave(paste(bam, "samtools_gcDepth_plot.png", sep = "_"))

  # Save insert size histogram table
  write_tsv(gcDepth, paste(bam, "samtools_gcDepth_histogram.tsv", sep = "_"))
}

# Summarize indels by cycle and size
indel_summary <- function(input, bam, rgsm, rglb, rgid) {
  # Extract indel size distribution
  indel_size <- grep("^ID",input, value=TRUE)
  indel_size <- separate(tibble(indel_size),
                         col=1,
                         into=c("ID", "Length", "Insertion_Count", "Deletion_Count"),
                         sep="\t") %>%
    type_convert() %>%
    select(-ID)

  # Pivot for graphing
  indel_size_long <- indel_size %>% pivot_longer(-Length, names_to = "Source", values_to = "Count")

  # Graph
  ggplot(indel_size_long, aes(Length, Count, color = Source)) +
    geom_line() +
    scale_y_log10(name = "INDEL Count",
                  breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000, 10000000),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)), ) +
    annotation_logticks(sides = "l") +
    scale_color_manual(values = c("#d8b365", "#5ab4ac"), labels = c("Deletion", "Insertion")) +
    ggtitle(label = bam) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="top")
  ggsave(paste(bam, "samtools_indelSize_linePlot.png", sep = "_"))

  # Extreact indel cycle distribution
  indel_cycle <- grep("^IC",input, value=TRUE)
  indel_cycle <- separate(tibble(indel_cycle),
                         col=1,
                         into=c("ID", "Cycle", "Insertion_Forward_Count", "Insertion_Reverse_Count", "Deletion_Forward_Count", "Deletion_Reverse_Count"),
                         sep="\t") %>%
    type_convert() %>%
    select(-ID)

  # Pivot for graphing
  indel_cycle_long <- indel_cycle %>% pivot_longer(-Cycle, names_to = "Source", values_to = "Count")

  # Graph
  ggplot(indel_cycle_long, aes(Cycle, Count, color = Source)) +
    geom_line() +
    scale_y_continuous(name = "INDEL Count",
                       labels = scales::comma) +
    scale_color_manual(values = c("#a6611a", "#dfc27d", "#80cdc1", "#018571"), labels = c("Deletion-Forward", "Deletion-Reverse", "Insertion-Forward", "Insertion-Reverse")) +
    ggtitle(label = bam) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="top")
  ggsave(paste(bam, "samtools_indelDistByCycle_linePlot.png", sep = "_"))

}

#### Insert Summary
insertSize_summary <- function(input, bam, rgsm, rglb, rgid) {
  # Extract the firstRead length
  read1_length <- grep("^FRL",input, value=TRUE)
  read1_length <- separate(tibble(read1_length),
                           col=1,
                           into=c("ID", "Read_Length", "Count"),
                           sep="\t") %>%
    type_convert()
  # Because there can be multiple read lengths used we capture the max read length value
  read1_length <- max(read1_length$Read_Length)

  # Extract the lastRead length
  read2_length <- grep("^LRL",input, value=TRUE)
  read2_length <- separate(tibble(read2_length),
                           col=1,
                           into=c("ID", "Read_Length", "Count"),
                           sep="\t") %>%
    type_convert()
  # Because there can be multiple read lengths used we capture the max read length value
  read2_length <- max(read2_length$Read_Length)

  # Import inserts size data table
  insertSize <- grep("^IS",input, value=TRUE)
  insertSize <- separate(tibble(insertSize),
                         col=1,
                         into=c("ID", "Insert_Size", "Total_Pairs", "Inward_Pairs", "Outward_Pairs", "Other_Pairs"),
                         sep="\t") %>%
    type_convert() %>%
    select(-ID)

  # Add meta-data to insert size table (histogram)
  insertSize <- insertSize %>%
    add_column(Sample = rgsm) %>%
    add_column(Library = rglb) %>%
    add_column(Read_Group = rgid) %>%
    add_column(BAM_File = bam)

  # Save insert size histogram table
  write_tsv(insertSize, paste(bam, "samtools_insertSize_histogram.tsv", sep = "_"))
  #insertSize <- read_tsv(paste(bam, "samtools_insertSize_histogram.tsv", sep = "_"))

  # Calculate the combined read length
  combined_read_length <- read1_length + read2_length

  # Add column to table to indicate if below combined read length
  insertSize <- insertSize %>%
    mutate(Last_Ordered_Pair = cumsum(Total_Pairs)) %>%
    mutate(First_Ordered_Pair = Last_Ordered_Pair - Total_Pairs + 1) %>%
    mutate(Below_Paired_Lenth = if_else(Insert_Size < combined_read_length, Total_Pairs, 0, missing = NULL)) %>%
    mutate(Bases_Lost = if_else(Below_Paired_Lenth > 0, (combined_read_length - Insert_Size) * Total_Pairs, 0, missing = NULL)) %>%
    mutate(Total_Pairs_Thousands = Total_Pairs / 1000) %>%
    mutate(Inward_Pairs_Thousands = Inward_Pairs / 1000) %>%
    mutate(Outward_Pairs_Thousands = Outward_Pairs / 1000) %>%
    mutate(Total_Insert_Space = Insert_Size * Total_Pairs)

  # Calculate MEAN, for SD calculation
  fast_mean <- insertSize %>%
    summarise(Pairs = sum(Total_Pairs),
              InsertSpace = sum(Total_Insert_Space)) %>%
    mutate(MEAN = InsertSpace / Pairs) %>%
    pull(var = MEAN)

  # Add column to insertSize table SD_lineNumerator = ((InsertSize - MEAN)^2) * Total_Pairs
  insertSize <- insertSize %>%
    mutate(SD_lineNumerator = ((Insert_Size - fast_mean)^2) * Total_Pairs)

  # Summarize InsertSize table
  summary <- insertSize %>%
    summarise(Total_Pairs = sum(Total_Pairs),
              Total_Insert_Length = sum(Total_Insert_Space),
              R1_Length = read1_length,
              R2_Length = read2_length,
              Combined_Read_Length = combined_read_length,
              Pairs_Below_Read_Length = sum(Below_Paired_Lenth),
              Bases_Lost_By_Overlap = sum(Bases_Lost),
              SD_lineNumerator = sum(SD_lineNumerator))

  # Get position information for insertSize percentiles
  median_position <- round(summary %>% pull(var = Total_Pairs) * 0.5)
  median_value <- insertSize %>% filter(First_Ordered_Pair <= median_position) %>% filter(Last_Ordered_Pair >= median_position) %>% pull(var = Insert_Size)
  percentile1_position <- round(summary %>% pull(var = Total_Pairs) * 0.01)
  percentile1_value <- insertSize %>% filter(First_Ordered_Pair <= percentile1_position) %>% filter(Last_Ordered_Pair >= percentile1_position) %>% pull(var = Insert_Size)
  percentile5_position <- round(summary %>% pull(var = Total_Pairs) * 0.05)
  percentile5_value <- insertSize %>% filter(First_Ordered_Pair <= percentile5_position) %>% filter(Last_Ordered_Pair >= percentile5_position) %>% pull(var = Insert_Size)
  percentile10_position <- round(summary %>% pull(var = Total_Pairs) * 0.10)
  percentile10_value <- insertSize %>% filter(First_Ordered_Pair <= percentile10_position) %>% filter(Last_Ordered_Pair >= percentile10_position) %>% pull(var = Insert_Size)
  percentile25_position <- round(summary %>% pull(var = Total_Pairs) * 0.25)
  percentile25_value <- insertSize %>% filter(First_Ordered_Pair <= percentile25_position) %>% filter(Last_Ordered_Pair >= percentile25_position) %>% pull(var = Insert_Size)
  percentile75_position <- round(summary %>% pull(var = Total_Pairs) * 0.75)
  percentile75_value <- insertSize %>% filter(First_Ordered_Pair <= percentile75_position) %>% filter(Last_Ordered_Pair >= percentile75_position) %>% pull(var = Insert_Size)
  percentile90_position <- round(summary %>% pull(var = Total_Pairs) * 0.90)
  percentile90_value <- insertSize %>% filter(First_Ordered_Pair <= percentile90_position) %>% filter(Last_Ordered_Pair >= percentile90_position) %>% pull(var = Insert_Size)
  percentile95_position <- round(summary %>% pull(var = Total_Pairs) * 0.95)
  percentile95_value <- insertSize %>% filter(First_Ordered_Pair <= percentile95_position) %>% filter(Last_Ordered_Pair >= percentile95_position) %>% pull(var = Insert_Size)
  percentile99_position <- round(summary %>% pull(var = Total_Pairs) * 0.99)
  percentile99_value <- insertSize %>% filter(First_Ordered_Pair <= percentile99_position) %>% filter(Last_Ordered_Pair >= percentile99_position) %>% pull(var = Insert_Size)

  # Determin the MAD
  mad_insertSize <- insertSize %>%
    select(Insert_Size, Total_Pairs) %>%
    mutate(Abs_Deviation = abs(Insert_Size - median_value)) %>%
    arrange(desc(Abs_Deviation)) %>%
    mutate(Last_Ordered_Pair = cumsum(Total_Pairs)) %>%
    mutate(First_Ordered_Pair = Last_Ordered_Pair - Total_Pairs + 1) %>%
    filter(First_Ordered_Pair <= median_position) %>%
    filter(Last_Ordered_Pair > median_position) %>%
    pull(var = Abs_Deviation)

  # Add Summary Calculations
  summary <- summary %>%
    mutate(Mean = Total_Insert_Length / Total_Pairs) %>%
    mutate(Variance = (1 / Total_Pairs) * SD_lineNumerator) %>%
    mutate(SD = sqrt(Variance)) %>%
    mutate(Median = median_value) %>%
    mutate(MAD = mad_insertSize) %>%
    mutate(Percentile_1 = percentile1_value) %>%
    mutate(Percentile_5 = percentile5_value) %>%
    mutate(Percentile_10 = percentile10_value) %>%
    mutate(Percentile_25 = percentile25_value) %>%
    mutate(Percentile_75 = percentile75_value) %>%
    mutate(Percentile_90 = percentile90_value) %>%
    mutate(Percentile_95 = percentile95_value) %>%
    mutate(Percentile_99 = percentile99_value) %>%
    mutate(IQR = Percentile_75 - Percentile_25) %>%
    mutate(Total_Bases = Total_Pairs * Combined_Read_Length) %>%
    mutate(Pct_Below_PairedReadLength = Pairs_Below_Read_Length / Total_Pairs) %>%
    mutate(Pct_Bases_Lost = Bases_Lost_By_Overlap / Total_Bases ) %>%
    add_column(Sample = rgsm) %>%
    add_column(Library = rglb) %>%
    add_column(Read_Group = rgid) %>%
    add_column(BAM_File = bam)

  # Write out the summary table
  write_tsv(summary, paste(bam, "samtools_insertSize_summary.tsv", sep = "_"))

  # Define ploting positions for text blob
  max_total_pairs <- insertSize %>% filter(Total_Pairs == max(Total_Pairs)) %>% pull(var = Total_Pairs_Thousands)
  max_insert_length <- insertSize %>% filter(Insert_Size == max(Insert_Size)) %>% pull(var = Insert_Size)
  if(max_insert_length < 1000) {
    x_limit <- 1000
    max_insert_length <- 950
  } else {
    x_limit <- max_insert_length + 200
  }

  # Capture variables needed for plot
  mean_insert_length <- pull(summary, var = Mean) %>% round(digits = 3)
  standard_deviation <- pull(summary, var = SD) %>% round(digits = 3)
  median_insert_length <- pull(summary, var = Median) %>% round(digits = 3)
  mad_insert_length <- pull(summary, var = MAD) %>% round(digits = 3)
  pct_below_pairedlength <- pull(summary, var = Pct_Below_PairedReadLength) %>% round(digits = 3)
  pct_bases_lost <- pull(summary, var = Pct_Bases_Lost) %>% round(digits = 3)

  # Generate a summary text blob to add to plot
  summary_label <- paste(paste("Mean(SD): ", mean_insert_length, " +/- ", standard_deviation, sep = ""),
                         paste("Median(MAD): ", median_insert_length, " +/- ", mad_insert_length, sep = ""),
                         paste("Overlapping Pairs: ", pct_below_pairedlength * 100, "%", sep = ""),
                         paste("Bases Lost: ", pct_bases_lost * 100, "%", sep = ""),
                         sep="\n")

  # Make long format file for graphing inward and outward pairs
  long_insertSize <- insertSize %>%
    select(Insert_Size, Inward_Pairs_Thousands, Outward_Pairs_Thousands) %>%
    pivot_longer(-Insert_Size, names_to = "Pair_Type", values_to = "Count")

  # Plot
  ggplot(long_insertSize, aes(x = Insert_Size, y = Count, fill = Pair_Type)) +
    geom_density(stat="identity", alpha = 0.8) +
    geom_vline(xintercept = mean_insert_length, color = "#5ab4ac", linetype = "dashed") +
    geom_vline(xintercept = combined_read_length, color = "#5ab4ac", linetype = "longdash") +
    annotate(geom = "text", x=max_insert_length * 0.85, y=max_total_pairs * 0.85, label = summary_label) +
    scale_x_continuous(name="Insert Size (bp)", limits = c(0, x_limit)) +
    scale_y_continuous(name="Total Pairs (Thousands)") +
    scale_fill_manual(values = c("#d8b365", "#5ab4ac"),
                      name = "Pair Type",
                      labels = c("Inward Pairs","Outward Pairs")) +
    ggtitle(bam) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="top")
  ggsave(paste(bam, "samtools_insertSize_histogram.png", sep = "_"))
}

#### Samtools markdups stats table processed as follows is the the expected input:
##
markdup_summary <- function(file, bam, rgsm, rglb, rgid) {
  # Import data table
  mdups <- read_lines(file)
  mdups <- separate(tibble(mdups),
                    col=1,
                    into=c("Name","Value"),
                    sep=": ")  %>%
    type_convert() %>%
    filter(Name != "COMMAND") %>%
    mutate(Name = str_replace_all(Name, " ", "_"))

  flipped_mdups <- mdups %>%
    pivot_wider(names_from = Name, values_from = Value) %>%
    type_convert()

  ## NOTES:
  ### DUPLICATE TOTAL = DUPLICATE PAIR + DUPLICATE SINGLE + DUPLICATE NON PRIMARY
  ### DUPLICATE PRIMARY TOTAL = DUPLICATE PAIR + DUPLICATE SINGLE
  flipped_mdups <- flipped_mdups %>%
    mutate(PERCENT_TOTAL_DUPLICATES = DUPLICATE_TOTAL / EXAMINED) %>%
    mutate(PERCENT_PRIMARY_PLATFORM_DUPLICATES = (DUPLICATE_PAIR_OPTICAL + DUPLICATE_SINGLE_OPTICAL) / DUPLICATE_PRIMARY_TOTAL) %>%
    mutate(PERCENT_TOTAL_PLATFORM_DUPLICATES = (DUPLICATE_PAIR_OPTICAL + DUPLICATE_SINGLE_OPTICAL + DUPLICATE_NON_PRIMARY_OPTICAL) / DUPLICATE_TOTAL)

  # Add columns with sample meta data
  flipped_mdups <- flipped_mdups %>%
    add_column(Sample = rgsm) %>%
    add_column(Library = rglb) %>%
    add_column(Read_Group = rgid) %>%
    add_column(BAM_File = bam)

  # Save samtools markdups summary table
  write_tsv(flipped_mdups, paste(bam, "samtools_markdup_summary.tsv", sep = "_"))
}

# Function to summarize samtools stats summary numbers
summaryNumbers_summary <- function(input, bam, rgsm, rglb, rgid, format) {
  # Import data table
  sn <- grep("^SN",input, value=TRUE)
  sn <- separate(tibble(sn),
                 col=1,
                 into=c("ID", "Name","Value"), sep="\t") %>%
    type_convert() %>%
    select(-ID)

  # Clean up names
  sn <- sn %>%
    mutate(Name = str_replace_all(Name, " ", "_")) %>%
    mutate(Name = str_replace_all(Name, ":", "")) %>%
    mutate(Name = str_replace_all(Name, "-", "_")) %>%
    mutate(Name = str_replace_all(Name, "_\\(%\\)", "")) %>%
    mutate(Name = str_replace_all(Name, "\\(", "")) %>%
    mutate(Name = str_replace_all(Name, "\\)", "")) %>%
    mutate(Name = str_replace_all(Name, "_>_0", "_above0x"))

  # Transpose and ensure order matches original vertical table order and update types from character
  flipped_sn <- sn %>% spread(Name, Value) %>% select(sn$Name) %>% type_convert()

  # Add column for chimeric reads (if single-end, chimeras cannot be determined so force to 0)
  if (format == "PairedEnd") {
    flipped_sn <- flipped_sn %>%
      mutate(PERCENT_CHIMERA = pairs_on_different_chromosomes / reads_mapped_and_paired)
  } else {
    flipped_sn <- flipped_sn %>%
      mutate(PERCENT_CHIMERA = 0)
  }


  # Add columns with sample meta data
  flipped_sn <- flipped_sn %>%
    add_column(Sample = rgsm) %>%
    add_column(Library = rglb) %>%
    add_column(Read_Group = rgid) %>%
    add_column(BAM_File = bam)

  # Save samtools markdups summary table
  write_tsv(flipped_sn, paste(bam, "samtools_summaryNumbers_summary.tsv", sep = "_"))
}


####################################
## Determine Base Directory
####################################

# Capture the directory of this script as functions are in files with relative locations to the primary script
## This is harder than I expected
## (https://stackoverflow.com/questions/1815606/determine-path-of-the-executing-script)
argv <- commandArgs(trailingOnly = FALSE)
# base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))

####################################
## Read in the samtools stats output file
####################################

stats_file <- read_lines(opt$samtoolsStatsFile)

####################################
## Call Functions
####################################

## Depending on input execute the appropriate summary functions

# Execute summaries of samtools stats output
if (!is.null(opt$samtoolsStatsFile)) {

  # Extract summary stats
  print("Summarizing Samtools stats Summary Numbers Statistics:")
  # Call Coverage Summary
  summaryNumbers_summary(stats_file, opt$bam, opt$sample, opt$library, opt$readgroup, opt$readformat)

  # Calculate coverage
  print("Summarizing Samtools stats Coverage Statistics:")
  # Call Coverage Summary
  median_cov <- coverage_summary(stats_file, opt$bam, opt$sample, opt$library, opt$readgroup)

  if ( opt$readformat == "PairedEnd") {
    # Calculate insert size and fraction of reads crossing over
    print("Summarizing Samtools stats Insert Size Statistics:")
    # Call the insertSize function
    insertSize_summary(stats_file, opt$bam, opt$sample, opt$library, opt$readgroup)
  }

  # Calcualte per cycle and overal base quality distributions
  print("Summarizing Samtools stats Base Quality Statistics:")
  # Call the base quality function
  baseQuality_summary(stats_file, opt$bam, opt$sample, opt$library, opt$readgroup, opt$readformat)
  
  # Plot Base distribution by Cycle
  print("Summarizing Samtools stats Base Distribution per Cycle:")
  # Call the base distribution function
  baseDistribution_summary(stats_file, opt$bam, opt$sample, opt$library, opt$readgroup)
  
  # Summarize GC data
  print("Summarizing Samtools stats GC Effect on Coverage:")
  # Call the base distribution function (mean_cov is returned from the coverage function)
  gcdepth_summary(stats_file, median_cov, opt$bam, opt$sample, opt$library, opt$readgroup)

  # Summarize INDEL data
  print("Summarizing Samtools stats INDEL size and distribution:")
  # Call the base distribution function (mean_cov is returned from the coverage function)
  indel_summary(stats_file, opt$bam, opt$sample, opt$library, opt$readgroup)
}

# Execute summaries of samtools markduplicates output
if (!is.null(opt$samtoolsDuplicatesFile)) {
  print("Summarizing Samtools markdups Statistics:")
  # Call Coverage Summary
  markdup_summary(opt$samtoolsDuplicatesFile, opt$bam, opt$sample, opt$library, opt$readgroup)
}

