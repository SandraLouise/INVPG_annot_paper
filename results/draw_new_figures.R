################################################################################
## Import packages
library(ggplot2)
library(magrittr)
library(dplyr)

################################################################################
## Set input paths

# Simulated data
found_missing_path = "../simulated_CHM13_chr21/merged_2haps_div1.found_missing.txt"

# Human data
hgsvc2_chr7_path = "../INVPG-annot/0_Data/HUMAN/HGSVC2_hg38_chr7_inv.bed"
hgsvc2_chrX_path = "../INVPG-annot/0_Data/HUMAN/HGSVC2_hg38_chrX_inv.bed"

inv_hg_chr7_path = "../GRCh38/HG_chr7_merged.bed"
inv_hg_chrX_path = "../GRCh38/HG_chrX_merged.bed"

################################################################################
## Data functions
add_size_category_col <- function(df){
  df$INV_SIZE = "50 bp - 1 Kb"
  df[(df$END - df$START + 1) >= 1000,]$INV_SIZE = "1 - 5 Kb"
  df[(df$END - df$START + 1) >= 5000,]$INV_SIZE = "5 Kb - 100 Kb"
  df[(df$END - df$START + 1) >= 100000,]$INV_SIZE = "> 100 Kb"
  df$INV_SIZE = factor(df$INV_SIZE, levels = c("50 bp - 1 Kb", "1 - 5 Kb", "5 Kb - 100 Kb", "> 100 Kb"))
  return(df)
}

add_redundancy_col <- function(df){
  df$REDUNDANT = F
  for (i in 2:nrow(df)){
    if (df[i,]$PG == df[i-1,]$PG & df[i,]$START < df[i-1,]$END + 1){
      df[i,]$REDUNDANT = T
    }
  }
  return(df)
}

################################################################################
## Plotting functions

################################################################################
## Import data

# Simulated data

found_missing = read.table(found_missing_path, h=T, sep="\t")
recall <- found_missing[found_missing$quality == "Precise",]
recall$type = factor(recall$type, levels = c("Path-explicit", "Alignment-rescued", "Mixed"))

# Human data
# truth set
hgsvc2_chr7 = read.table(hgsvc2_chr7_path, h=F, sep="\t")
hgsvc2_chrX = read.table(hgsvc2_chrX_path, h=F, sep="\t")
colnames(hgsvc2_chr7) = c("Chr", "START", "END")
colnames(hgsvc2_chrX) = c("Chr", "START", "END")

hgsvc2_chr7 = add_size_category_col(hgsvc2_chr7)
hgsvc2_chrX = add_size_category_col(hgsvc2_chrX)

hgsvc2_merged = rbind(hgsvc2_chr7, hgsvc2_chrX)

input_nINV_merged <- hgsvc2_merged %>%
  group_by(INV_SIZE) %>%
  summarize(Count = n())

# annotated set
inv_hg_chr7 = read.table(inv_hg_chr7_path, h=F, sep="\t")
inv_hg_chrX = read.table(inv_hg_chrX_path, h=F, sep="\t")
colnames(inv_hg_chr7) = c("PG", "Chr", "START", "END", "Annot")
colnames(inv_hg_chrX) = c("PG", "Chr", "START", "END", "Annot")

inv_hg_chr7 = add_size_category_col(inv_hg_chr7)
inv_hg_chrX = add_size_category_col(inv_hg_chrX)
inv_hg_chr7 = add_redundancy_col(inv_hg_chr7)
inv_hg_chrX = add_redundancy_col(inv_hg_chrX)

inv_hg_merged = rbind(inv_hg_chr7, inv_hg_chrX)

################################################################################
## Draw plots

# Barplot Recall

barplot_recall <- recall %>%
  group_by(pg, type) %>%
  summarize(Count = n()) %>%
  ggplot(aes(x=pg, y=Count, fill=type)) +
  geom_bar(stat='identity', position= "stack", width = 0.5) +
  ylim(ymin=0, ymax=100) +
  labs(title="A", x="Pangenome graph pipeline", y = "Recall (%)", fill = "Annotation signal type") +
  scale_color_manual(values=c("#D1495B","#EDAE49"), aesthetics = "fill") +
  theme_minimal() +
  theme(legend.position = 'bottom')

# Barplot Number of TP

barplot_TP <- inv_hg_merged[inv_hg_merged$REDUNDANT == F,] %>%
  group_by(PG, INV_SIZE) %>%
  summarize(Count = n()) %>%
  ggplot(aes(x=INV_SIZE, y=Count, fill=PG)) +
  geom_rect(xmin=1.6, xmax=2.4, ymin=0, ymax=input_nINV_merged[input_nINV_merged$INV_SIZE == levels(input_nINV_merged$INV_SIZE)[2],]$Count,
            fill='grey', alpha=0.05) +
  geom_rect(xmin=2.6, xmax=3.4, ymin=0, ymax=input_nINV_merged[input_nINV_merged$INV_SIZE == levels(input_nINV_merged$INV_SIZE)[3],]$Count,
            fill='grey', alpha=0.05) +
  geom_rect(xmin=3.6, xmax=4.4, ymin=0, ymax=input_nINV_merged[input_nINV_merged$INV_SIZE == levels(input_nINV_merged$INV_SIZE)[4],]$Count,
            fill='grey', alpha=0.05) +
  geom_bar(stat='identity', position= position_dodge(), width = 0.5) +
  # ylim(ymin=0, ymax=100) +
  labs(title="A", x="Inversion size category", y = "# Annotated inversions", fill = "Pangenome graph pipeline") +
  theme_minimal() +
  theme(legend.position = 'right', plot.title = element_text(face = "bold"),
        axis.title.x = element_blank())
