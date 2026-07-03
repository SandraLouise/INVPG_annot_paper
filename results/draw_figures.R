#===============================================================================
## Import packages
#===============================================================================
library(ggplot2)
library(magrittr)
library(dplyr)
library(tidyr)
library(knitr)
library(ggvenn)
library(ggrepel)
library(scales)

#===============================================================================
# Functions
#===============================================================================
add_size_category_col <- function(df){
  df$INV_SIZE = "50 bp - 1 Kb"
  df[(df$END_truth - df$START_truth + 1) >= 1000,]$INV_SIZE = "1 - 5 Kb"
  df[(df$END_truth - df$START_truth + 1) >= 5000,]$INV_SIZE = "5 Kb - 100 Kb"
  df[(df$END_truth - df$START_truth + 1) >= 100000,]$INV_SIZE = "> 100 Kb"
  df$INV_SIZE = factor(df$INV_SIZE, levels = c("50 bp - 1 Kb", "1 - 5 Kb", "5 Kb - 100 Kb", "> 100 Kb"))
  return(df)
}

add_redundancy_col <- function(df){
  df$REDUNDANT = F
  for (i in 2:nrow(df)){
    if (df[i,]$PG == df[i-1,]$PG & df[i,]$START_truth < df[i-1,]$END_truth + 1){
      df[i,]$REDUNDANT = T
    }
  }
  return(df)
}

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

df_recall_accuracy <- function(recall_bkpt_accur, array_threshold, n_haplotypes){
  
  # array_threshold = offset_thresholds
  
  recall_bkpt_accur$START = recall_bkpt_accur$START_annot
  recall_bkpt_accur$END = recall_bkpt_accur$END_annot
  recall_bkpt_accur <- subset(recall_bkpt_accur, REDUNDANT == F)
  
  tmp_recall_bkpt_accur <- recall_bkpt_accur
  
  new_df <- tmp_recall_bkpt_accur %>% 
    group_by(HAPS, PG) %>% 
    summarize(Count = n())
  
  new_df$THRESHOLD = 100 # placeholder name for no threshold
  
  for (i_threshold in array_threshold) {
    
    tmp_recall_bkpt_accur <- recall_bkpt_accur
    tmp_recall_bkpt_accur$PASS_BKPT_ACCUR = F
    tmp_recall_bkpt_accur[tmp_recall_bkpt_accur$START_OFFSET < i_threshold & tmp_recall_bkpt_accur$END_OFFSET < i_threshold,]$PASS_BKPT_ACCUR = T
    
    tmp_recall_bkpt_accur <- subset(tmp_recall_bkpt_accur, PASS_BKPT_ACCUR == T)
    
    tmp_df <- tmp_recall_bkpt_accur %>% 
      group_by(HAPS, PG) %>% 
      summarize(Count = n())
    
    tmp_df$THRESHOLD = i_threshold
    
    new_df <- rbind(new_df, tmp_df)
    
  }
  
  return(new_df)
}

#===============================================================================
# Set input paths
#===============================================================================
#-------------------------------------------------------------------------------
# MASTER PATH TO INPUT FILES FOLDER
#-------------------------------------------------------------------------------
MASTER_PATH = "~/Documents/Paper/revisions/Final_figures/"

#-------------------------------------------------------------------------------
# Simulated data (CHM13)
#-------------------------------------------------------------------------------
#---------------------------------------
# Truth set
#---------------------------------------
sim_truth_path = paste0(MASTER_PATH, "Main/simulated_data_CHM13/truth_set/chr21_sim_100inv.bed")

#---------------------------------------
# Evaluation statistics
#---------------------------------------
# Path to the folder containing the eval files
eval_2haps_path = paste0(MASTER_PATH, "Main/simulated_data_CHM13/eval_stats/")

#---------------------------------------
# Recall intersect
#---------------------------------------
# Path to the folder containing the intersect files
intersect_2haps_path = paste0(MASTER_PATH, "Main/simulated_data_CHM13/recall_intersect/")
intersect_10haps_path = paste0(MASTER_PATH, "Main/simulated_data_CHM13/recall_intersect/")

#---------------------------------------
# Inv density
#---------------------------------------
inv_density_path = paste0(MASTER_PATH, "Main/simulated_data_CHM13/inv_density/simHG_Ninv_recall.txt")
inv_density = read.table(inv_density_path, h=T, sep="\t")
inv_density$Recall = inv_density$PassingRecall / inv_density$N_inv * 100
inv_density$PG = factor(inv_density$PG, levels = c("Cactus", "Minigraph", "Minigraph-Cactus", "PGGB"))

#---------------------------------------
# Div and cov parameters
#---------------------------------------
div_parameter_path = paste0(MASTER_PATH, "Supplementary/cov_div_parameters/2hap_div1_max_div.txt")
cov_parameter_path = paste0(MASTER_PATH, "Supplementary/cov_div_parameters/2hap_div1_min_cov.txt")

#-------------------------------------------------------------------------------
# Real data (GRCh38)
#-------------------------------------------------------------------------------
#---------------------------------------
# Truth set
#---------------------------------------
# genouest: /scratch/sromain/paper_PG_INV/data/human_data/HGSVC2_hg38_chr7_inv.bed
hgsvc2_chr7_path = paste0(MASTER_PATH, "Main/real_data_GRCh38/truth_set/HGSVC2_hg38_chr7_inv.bed")
# genouest: /scratch/sromain/paper_PG_INV/data/human_data/HGSVC2_hg38_chrX_inv.bed
hgsvc2_chrX_path = paste0(MASTER_PATH, "Main/real_data_GRCh38/truth_set/HGSVC2_hg38_chrX_inv.bed")

# INVPG-annot results
# genouest: /scratch/sromain/paper_PG_INV/results/GRCh38_chr7/inversions/HG_chr7_[pg].bed
inv_hg_chr7_path = list(
  "Cactus" = paste0(MASTER_PATH, "Main/real_data_GRCh38/invpg_beds/HG_chr7_cactus.bed"),
  "Minigraph" = paste0(MASTER_PATH, "Main/real_data_GRCh38/invpg_beds/HG_chr7_minigraph.bed"),
  "Minigraph-Cactus" = paste0(MASTER_PATH, "Main/real_data_GRCh38/invpg_beds/HG_chr7_mgc.bed"),
  "PGGB" = paste0(MASTER_PATH, "Main/real_data_GRCh38/invpg_beds/HG_chr7_pggb_v0.7.4_HGparam.bed")
)
# genouest: /scratch/sromain/paper_PG_INV/results/GRCh38_chr7/inversions/HG_chrX_[pg].bed
inv_hg_chrX_path = list(
  "Cactus" = paste0(MASTER_PATH, "Main/real_data_GRCh38/invpg_beds/HG_chrX_cactus.bed"),
  "Minigraph" = paste0(MASTER_PATH, "Main/real_data_GRCh38/invpg_beds/HG_chrX_minigraph.bed"),
  "Minigraph-Cactus" = paste0(MASTER_PATH, "Main/real_data_GRCh38/invpg_beds/HG_chrX_mgc.bed"),
  "PGGB" = paste0(MASTER_PATH, "Main/real_data_GRCh38/invpg_beds/HG_chrX_pggb_v0.7.4_HGparam.bed")
)

#---------------------------------------
# Intersect
#---------------------------------------
# genouest: /scratch/sromain/paper_PG_INV/results/GRCh38_chr7/inversions/HG_chr7_[pg].bed.recall_intersect
intersect_hg_chr7_path = list(
  "Cactus" = paste0(MASTER_PATH, "Main/real_data_GRCh38/intersect_beds/HG_chr7_cactus.bed.recall_intersect"),
  "Minigraph" = paste0(MASTER_PATH, "Main/real_data_GRCh38/intersect_beds/HG_chr7_minigraph.bed.recall_intersect"),
  "Minigraph-Cactus" = paste0(MASTER_PATH, "Main/real_data_GRCh38/intersect_beds/HG_chr7_mgc.bed.recall_intersect"),
  "PGGB" = paste0(MASTER_PATH, "Main/real_data_GRCh38/intersect_beds/HG_chr7_pggb_v0.7.4_HGparam.bed.recall_intersect")
)
# genouest: /scratch/sromain/paper_PG_INV/results/GRCh38_chr7/inversions/HG_chrX_[pg].bed.recall_intersect
intersect_hg_chrX_path = list(
  "Cactus" = paste0(MASTER_PATH, "Main/real_data_GRCh38/intersect_beds/HG_chrX_cactus.bed.recall_intersect"),
  "Minigraph" = paste0(MASTER_PATH, "Main/real_data_GRCh38/intersect_beds/HG_chrX_minigraph.bed.recall_intersect"),
  "Minigraph-Cactus" = paste0(MASTER_PATH, "Main/real_data_GRCh38/intersect_beds/HG_chrX_mgc.bed.recall_intersect"),
  "PGGB" = paste0(MASTER_PATH, "Main/real_data_GRCh38/intersect_beds/HG_chrX_pggb_v0.7.4_HGparam.bed.recall_intersect")
)

#---------------------------------------
# Reverse intersect
#---------------------------------------
rev_intersect_hg_chr7_path = list(
  "Cactus" = paste0(MASTER_PATH, "Main/real_data_GRCh38/intersect_beds/HG_chr7_cactus.bed.reverse.recall_intersect"),
  "Minigraph" = paste0(MASTER_PATH, "Main/real_data_GRCh38/intersect_beds/HG_chr7_minigraph.bed.reverse.recall_intersect"),
  "Minigraph-Cactus" = paste0(MASTER_PATH, "Main/real_data_GRCh38/intersect_beds/HG_chr7_mgc.bed.reverse.recall_intersect"),
  "PGGB" = paste0(MASTER_PATH, "Main/real_data_GRCh38/intersect_beds/HG_chr7_pggb.bed.reverse.recall_intersect")
)
rev_intersect_hg_chrX_path = list(
  "Cactus" = paste0(MASTER_PATH, "Main/real_data_GRCh38/intersect_beds/HG_chrX_cactus.bed.reverse.recall_intersect"),
  "Minigraph" = paste0(MASTER_PATH, "Main/real_data_GRCh38/intersect_beds/HG_chrX_minigraph.bed.reverse.recall_intersect"),
  "Minigraph-Cactus" = paste0(MASTER_PATH, "Main/real_data_GRCh38/intersect_beds/HG_chrX_mgc.bed.reverse.recall_intersect"),
  "PGGB" = paste0(MASTER_PATH, "Main/real_data_GRCh38/intersect_beds/HG_chrX_pggb.bed.reverse.recall_intersect")
)

#===============================================================================
# Import data
#===============================================================================
color_palette = list(
  "topo" = c("Path-explicit" = "#D1495B",
             "Alignment-rescued" = "#EDAE49")
)

#-------------------------------------------------------------------------------
# Simulated data (CHM13)
#-------------------------------------------------------------------------------
#---------------------------------------
# Truth set
#---------------------------------------
sim_truth = read_tsv(sim_truth_path, col_names=F) %>%
  select(CHR = X1, START_truth = X2, END_truth = X3)

sim_truth = add_size_category_col(sim_truth)

#---------------------------------------
# Evaluation statistics
#---------------------------------------
# 2 haps
eval_2haps = rbind(
  read.table(paste0(eval_2haps_path, "2hap_div0.1_cactus.eval"), h=T, sep="\t"),
  read.table(paste0(eval_2haps_path, "2hap_div1_cactus.eval"), h=T, sep="\t"),
  read.table(paste0(eval_2haps_path, "2hap_div5_cactus.eval"), h=T, sep="\t")
)
for (pg in c("minigraph", "mgc", "pggb")){
  for (div in c("0.1", "1", "5")){
    eval_2haps = rbind(
      eval_2haps,
      read.table(paste0(eval_2haps_path, "2hap_div", div, "_", pg, ".eval", collapse=""), 
                 h=T, sep="\t")
    )
  }
}

eval_2haps = eval_2haps %>%
  mutate(
    PG = case_when(
      PG == "cactus" ~ "Cactus",
      PG == "minigraph" ~ "Minigraph",
      PG == "mgc" ~ "Minigraph-Cactus",
      PG == "pggb" ~ "PGGB"
    ),
    DIV = as.numeric(gsub(pattern="div", replacement="", DIV))
  )

eval_2haps$PG = factor(
  eval_2haps$PG, 
  levels = c("Cactus", "Minigraph", "Minigraph-Cactus", "PGGB")
)

# 10 haps
eval_10haps = rbind(
  read.table(paste0(eval_2haps_path, "10hap_div1_cactus.eval"), h=T, sep="\t"),
  read.table(paste0(eval_2haps_path, "10hap_div1_minigraph.eval"), h=T, sep="\t"),
  read.table(paste0(eval_2haps_path, "10hap_div1_mgc.eval"), h=T, sep="\t"),
  read.table(paste0(eval_2haps_path, "10hap_div1_pggb.eval"), h=T, sep="\t")
)

eval_10haps = eval_10haps %>%
  mutate(
    PG = case_when(
      PG == "cactus" ~ "Cactus",
      PG == "minigraph" ~ "Minigraph",
      PG == "mgc" ~ "Minigraph-Cactus",
      PG == "pggb" ~ "PGGB"
    ),
    DIV = as.numeric(gsub(pattern="div", replacement="", DIV))
  )

eval_10haps$PG = factor(
  eval_10haps$PG, 
  levels = c("Cactus", "Minigraph", "Minigraph-Cactus", "PGGB")
)

#---------------------------------------
# Recall intersect - 2 haps
#---------------------------------------
intersect_2haps = rbind(
  read.table(paste0(intersect_2haps_path, "2hap_div0.1_cactus.bed.recall_intersect"), h=F, sep="\t"),
  read.table(paste0(intersect_2haps_path, "2hap_div1_cactus.bed.recall_intersect"), h=F, sep="\t"),
  read.table(paste0(intersect_2haps_path, "2hap_div5_cactus.bed.recall_intersect"), h=F, sep="\t")
)
for (pg in c("minigraph", "mgc", "pggb")){
  for (div in c("0.1", "1", "5")){
    intersect_2haps = rbind(
      intersect_2haps,
      read.table(paste0(intersect_2haps_path, "2hap_div", div, "_", pg, ".bed.recall_intersect", collapse=""), 
                 h=F, sep="\t")
    )
  }
}

colnames(intersect_2haps) = c("PG", "DIV", "START_truth", "END_truth",
                              "START_annot", "END_annot", "Annot")

# Rename PG, DIV for display + extract topology from Annot
intersect_2haps = intersect_2haps %>%
  mutate(
    PG = case_when(
      PG == "cactus" ~ "Cactus",
      PG == "minigraph" ~ "Minigraph",
      PG == "mgc" ~ "Minigraph-Cactus",
      PG == "pggb" ~ "PGGB"
    ),
    DIV = as.numeric(gsub(pattern="div", replacement="", DIV)),
    ANNOT_TOPO = case_when(
      grepl("INV:path", Annot) ~ "Path-explicit",
      grepl("INV:aln", Annot) ~ "Alignment-rescued",
      T ~ NA
    )
  )

# Order PG pipelines
intersect_2haps$PG = factor(
  intersect_2haps$PG, 
  levels = c("Cactus", "Minigraph", "Minigraph-Cactus", "PGGB")
)
# Order topologies
intersect_2haps$ANNOT_TOPO = factor(
  intersect_2haps$ANNOT_TOPO, 
  levels = c("Path-explicit", "Alignment-rescued")
)

intersect_2haps = add_redundancy_col(intersect_2haps)
intersect_2haps = add_size_category_col(intersect_2haps)

#---------------------------------------
# Recall intersect - 10 haps
#---------------------------------------
intersect_10haps = rbind(
  read.table(paste0(intersect_10haps_path, "10hap_div1_cactus.bed.recall_intersect"), h=F, sep="\t"),
  read.table(paste0(intersect_10haps_path, "10hap_div1_minigraph.bed.recall_intersect"), h=F, sep="\t"),
  read.table(paste0(intersect_10haps_path, "10hap_div1_mgc.bed.recall_intersect"), h=F, sep="\t"),
  read.table(paste0(intersect_10haps_path, "10hap_div1_pggb.bed.recall_intersect"), h=F, sep="\t")
)
colnames(intersect_10haps) = c("PG", "DIV", "START_truth", "END_truth",
                              "START_annot", "END_annot", "Annot")

# Rename PG, DIV for display + extract topology from Annot
intersect_10haps = intersect_10haps %>%
  mutate(
    PG = case_when(
      PG == "cactus" ~ "Cactus",
      PG == "minigraph" ~ "Minigraph",
      PG == "mgc" ~ "Minigraph-Cactus",
      PG == "pggb" ~ "PGGB"
    ),
    DIV = as.numeric(gsub(pattern="div", replacement="", DIV)),
    ANNOT_TOPO = case_when(
      grepl("INV:path", Annot) ~ "Path-explicit",
      grepl("INV:aln", Annot) ~ "Alignment-rescued",
      T ~ NA
    )
  )

# Order PG pipelines
intersect_10haps$PG = factor(
  intersect_10haps$PG, 
  levels = c("Cactus", "Minigraph", "Minigraph-Cactus", "PGGB")
)
# Order topologies
intersect_10haps$ANNOT_TOPO = factor(
  intersect_10haps$ANNOT_TOPO, 
  levels = c("Path-explicit", "Alignment-rescued")
)

intersect_10haps = add_redundancy_col(intersect_10haps)
intersect_10haps = add_size_category_col(intersect_10haps)

#---------------------------------------
# Breakpoint precision
#---------------------------------------
n_true_inv = 100
offset_thresholds = c(10, 25, 50)
n_haplotypes = c(2, 10)

bkpt_precision = rbind(
  intersect_2haps %>% 
    subset(DIV == 1) %>%
    mutate(HAPS = 2),
  intersect_10haps %>% 
    subset(DIV == 1) %>%
    mutate(HAPS = 10)
) %>%
  mutate(
    START_OFFSET = abs(START_annot - START_truth),
    END_OFFSET = abs(END_annot - END_truth),
    SUM_OFFSET = START_OFFSET + END_OFFSET,
    SUM_OFFSET = case_when(
      SUM_OFFSET == 0 ~ 1,
      T ~ SUM_OFFSET
    )
  ) %>%
  subset(REDUNDANT == F & Annot != ".")

bkpt_precision = df_recall_accuracy(bkpt_precision, offset_thresholds, n_haplotypes)
bkpt_precision = bkpt_precision %>%
  mutate(
    THRESHOLD = as.character(THRESHOLD),
    THRESHOLD = case_when(
      THRESHOLD == "100" ~ "no threshold",
      T ~ THRESHOLD
    ),
    HAPS = as.character(HAPS)
  )
bkpt_precision$THRESHOLD = factor(
  bkpt_precision$THRESHOLD, 
  levels=c("no threshold", "50", "25", "10")
)
bkpt_precision$HAPS = factor(
  bkpt_precision$HAPS, 
  levels=c("2", "10")
)

#---------------------------------------
# Div and cov parameters
#---------------------------------------
div_parameter <- read.table(div_parameter_path, h = FALSE, sep = "\t", 
                            stringsAsFactors = FALSE)
cov_parameter <- read.table(cov_parameter_path, h = FALSE, sep = "\t", 
                            stringsAsFactors = FALSE)

# Noms des colonnes
colnames(div_parameter) <- c("tool", "dataset", "max_div", "metric", "value")
colnames(cov_parameter) <- c("tool", "dataset", "min_cov", "metric", "value")

# Garder uniquement Recall
div_parameter_recall <- div_parameter %>%
  filter(metric == "Recall") %>%
  mutate(tool = recode(tool,
                       "cactus" = "Cactus",
                       "mgc" = "Minigraph-Cactus",
                       "minigraph" = "Minigraph",
                       "pggb" = "PGGB")) %>%
  # S'assurer que max_div est numérique
  mutate(
    max_div = as.numeric(max_div),
    value = as.numeric(value)
  )

# Garder uniquement Recall
cov_parameter_recall <- cov_parameter %>%
  filter(metric == "Recall") %>%
  mutate(tool = recode(tool,
                       "cactus" = "Cactus",
                       "mgc" = "Minigraph-Cactus",
                       "minigraph" = "Minigraph",
                       "pggb" = "PGGB")) %>%
  # S'assurer que max_div est numérique
  mutate(
    min_cov = as.numeric(min_cov),
    value = as.numeric(value)
  )

#-------------------------------------------------------------------------------
# Real data (GRCh38)
#-------------------------------------------------------------------------------
# Truth set
hgsvc2_chr7 = read.table(hgsvc2_chr7_path, h=F, sep="\t")
hgsvc2_chrX = read.table(hgsvc2_chrX_path, h=F, sep="\t")
colnames(hgsvc2_chr7) = c("Chr", "START_truth", "END_truth")
colnames(hgsvc2_chrX) = c("Chr", "START_truth", "END_truth")

hgsvc2_chrX_merged = rbind(
  hgsvc2_chr7,
  hgsvc2_chrX
)

hgsvc2_chrX_merged = add_size_category_col(hgsvc2_chrX_merged) %>%
  group_by(INV_SIZE) %>%
  summarise(Count = n())


# INVPG-annot results
inv_hg_merged_all = rbind(
  # Chr7
  read.table(intersect_hg_chr7_path[[names(intersect_hg_chr7_path)[1]]], h=F, sep="\t") %>%
    mutate(PG = names(intersect_hg_chr7_path)[1]),
  read.table(intersect_hg_chr7_path[[names(intersect_hg_chr7_path)[2]]], h=F, sep="\t") %>%
    mutate(PG = names(intersect_hg_chr7_path)[2]),
  read.table(intersect_hg_chr7_path[[names(intersect_hg_chr7_path)[3]]], h=F, sep="\t") %>%
    mutate(PG = names(intersect_hg_chr7_path)[3]),
  read.table(intersect_hg_chr7_path[[names(intersect_hg_chr7_path)[4]]], h=F, sep="\t") %>%
    mutate(PG = names(intersect_hg_chr7_path)[4]),
  # ChrX
  read.table(intersect_hg_chrX_path[[names(intersect_hg_chrX_path)[1]]], h=F, sep="\t") %>%
    mutate(PG = names(intersect_hg_chrX_path)[1]),
  read.table(intersect_hg_chrX_path[[names(intersect_hg_chrX_path)[2]]], h=F, sep="\t") %>%
    mutate(PG = names(intersect_hg_chrX_path)[2]),
  read.table(intersect_hg_chrX_path[[names(intersect_hg_chrX_path)[3]]], h=F, sep="\t") %>%
    mutate(PG = names(intersect_hg_chrX_path)[3]),
  read.table(intersect_hg_chrX_path[[names(intersect_hg_chrX_path)[4]]], h=F, sep="\t") %>%
    mutate(PG = names(intersect_hg_chrX_path)[4])
) %>%
  select(-V1)
colnames(inv_hg_merged_all) = c("CHR", "START_truth", "END_truth", "START_annot", "END_annot", "Annot", "PG")

inv_hg_merged_all = subset(inv_hg_merged_all, Annot != ".")
inv_hg_merged_all = add_size_category_col(inv_hg_merged_all)
inv_hg_merged_all = add_redundancy_col(inv_hg_merged_all)


#===============================================================================
# MAIN
#===============================================================================
#-------------------------------------------------------------------------------
# Fig 3: Simulated data - Recall
#-------------------------------------------------------------------------------
#---------------------------------------
# Fig 3A
#---------------------------------------
plot1 <- intersect_2haps %>%
  subset(Annot != "." & REDUNDANT == F & DIV == 1) %>%
  group_by(PG, ANNOT_TOPO) %>%
  summarize(Count = n()) %>%
  ggplot(aes(x=PG, y=Count, fill=ANNOT_TOPO)) +
  geom_bar(stat='identity', position= "stack", width = 0.5) +
  ylim(ymin=0, ymax=100) +
  labs(title="A", x="Pangenome graph pipeline", y = "Recall (%)", fill = "Annotation signal type") +
  scale_color_manual(values=color_palette$topo, aesthetics = "fill") +
  theme_minimal() +
  theme(legend.position = 'bottom',
        plot.title = element_text(face = "bold"))

#---------------------------------------
# Fig 3B
#---------------------------------------
tot_truth_size_cat = sim_truth %>%
  group_by(INV_SIZE) %>%
  summarise(tot_truth = n())

plot2 <- intersect_2haps %>%
  # Organize and complete table
  subset(Annot != "." & REDUNDANT == F & DIV == 1) %>%
  group_by(PG, INV_SIZE) %>%
  summarise(count = n()) %>%
  left_join(tot_truth_size_cat, by="INV_SIZE") %>%
  mutate(recall = count / tot_truth * 100) %>%
  # Plot
  ggplot(aes(x=INV_SIZE, y=recall, group=PG)) +
  geom_line(aes(color=PG), linewidth=0.9, alpha=0.6)+
  geom_point(aes(color=PG, shape=PG), size=3, alpha=1) +
  ylim(ymin=0, ymax=100) +
  labs(title="B", x="Inversion size category", y = "Recall (%)", color = "Pangenome graph pipeline", shape = "Pangenome graph pipeline") +
  theme_minimal() +
  theme(legend.position = 'top',
        plot.title = element_text(face = "bold"))

#---------------------------------------
# Fig 3C --> TODO
#---------------------------------------

plot3 <- ggplot(inv_density, aes(x=Density_inv, y=Recall, group=PG)) +
  geom_line(aes(color=PG), linewidth=0.9, alpha=0.6) +
  geom_point(aes(color=PG, shape=PG), size=3, alpha=1) +
  ylim(60, 100) +
  labs(title="C", x="Inversion density (% of bp)", y = "Recall (%)") +
  theme_minimal() +
  theme(legend.position = 'bottom',
        plot.title = element_text(face = "bold")) +
  geom_label_repel(
    data = subset(inv_density, PG == "Minigraph"),
    aes(label = N_inv),
    size = 3.5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.01, "lines",),
    alpha = 0.6
  )

#---------------------------------------
# Fig 3D
#---------------------------------------
plot4 <- bkpt_precision %>%
  ggplot(aes(x=THRESHOLD, y=Count, group=interaction(PG, HAPS))) +
  geom_line(aes(color=PG, linetype = HAPS), linewidth=0.9, alpha=0.6)+
  geom_point(aes(color=PG, shape=PG), size=3, alpha=1) +
  ylim(ymin=40, ymax=100) +
  labs(title="D", x="Precision threshold (bp)", y = "Recall (%)", linetype = "Number of haplotypes") +
  theme_minimal() +
  guides(color = "none", shape = "none") +
  theme(legend.position = 'bottom',
        plot.title = element_text(face = "bold"))

#---------------------------------------
# Fig 3E
#---------------------------------------
redundancy_div1 = rbind(
  eval_2haps %>%
    subset(DIV == 1 & (Stat_name == "Recall" | Stat_name == "NonRedundant_bubbles")) %>%
    mutate(HAPS = "2"),
  eval_10haps %>%
    subset(DIV == 1 & (Stat_name == "Recall" | Stat_name == "NonRedundant_bubbles"))%>%
    mutate(HAPS = "10")
) %>%
  mutate(
    Stat_name = case_when(
      Stat_name == "NonRedundant_bubbles" ~ "No",
      Stat_name == "Recall" ~ "Yes"
    )
  )
redundancy_div1$HAPS = factor(redundancy_div1$HAPS,
                              levels=c("2","10"))
redundancy_div1$Stat_name = factor(redundancy_div1$Stat_name,
                              levels=c("Yes","No"))

plot5 <- redundancy_div1 %>% 
  ggplot(aes(x=HAPS, y=Stat_value, group=interaction(Stat_name, PG), color=PG)) +
  geom_line(aes(linetype=Stat_name), linewidth=0.9, alpha=0.6)+
  geom_point(aes(shape=PG), size=3, alpha=1) +
  ylim(ymin=25, ymax=100) +
  labs(title="E", x="Number of haplotypes", y = "Recall (%)", linetype = "Count redundant INV") +
  facet_grid(~ PG) +
  guides(color = "none", shape = "none") +
  theme_minimal() +
  theme(legend.position = 'bottom',
        plot.title = element_text(face = "bold"))

#---------------------------------------
# Arrange into Figure 3
#---------------------------------------
legend_annotType <- get_legend(plot1)
legend_PG <- get_legend(plot2)
# legend3_4 <- get_legend(plot3)
plot1 <- plot1 + theme(legend.position="none")
plot2 <- plot2 + theme(legend.position="none")
plot3 <- plot3 + theme(legend.position="none")
# plot4 <- plot4 + theme(legend.position="none")
# plot5 <- plot5 + theme(legend.position="none")
pdf(file=paste0(MASTER_PATH, "Main/Figure3.pdf"), 
    width=8, height=8)
grid.arrange(plot1, plot2, plot3, plot4, plot5, 
             legend_annotType, legend_PG, 
             ncol=2, 
             layout_matrix = cbind(c(6,1,7,2,4), c(6,1,7,3,5)),
             heights = c(0.5,3,0.8,3,3.5))
dev.off()

#-------------------------------------------------------------------------------
# Fig 4: Simulated data - Topology
#-------------------------------------------------------------------------------
#---------------------------------------
# Fig 4A
#---------------------------------------
aln_rescued_2haps = eval_2haps %>%
  subset(Stat_name == "Alignment-rescued")

plot_AlnResc_div <- ggplot(aln_rescued_2haps, 
                           aes(x=DIV, y=Stat_value, group=PG, color=PG)) +
  geom_line(linewidth=0.9, alpha=0.8) +
  geom_point(aes(color=PG, shape=PG), size=3, alpha=0.8) +
  labs(title = "A", 
       x = "Nucleotide polymorphism (%)", y = "Alignment-rescued (n)",
       color = "Pangenome graph pipeline", shape = "Pangenome graph pipeline") +
  # force the x ticks to be at 0.1, 1, 5
  scale_x_continuous(breaks = c(0.1, 1, 5)) +
  ylim(0, 100) +
  theme_minimal() +
  theme(legend.position = 'top', 
        plot.title = element_text(face = "bold"),
        panel.grid.minor.x = element_blank())

#---------------------------------------
# Fig 4B
#---------------------------------------
path_explicit_2haps = eval_2haps %>%
  subset(Stat_name == "Path-explicit")

plot_PathExp_div <- ggplot(path_explicit_2haps, 
                           aes(x=DIV, y=Stat_value, group=PG, color=PG)) +
  geom_line(linewidth=0.9, alpha=0.8) +
  geom_point(aes(color=PG, shape=PG), size=3, alpha=0.8) +
  labs(title = "B", 
       x = "Nucleotide polymorphism (%)", y = "Path-explicit (n)",
       color = "Pangenome graph pipeline", shape = "Pangenome graph pipeline") +
  # force the x ticks to be at 0.1, 1, 5
  scale_x_continuous(breaks = c(0.1, 1, 5)) +
  ylim(0, 100) +
  theme_minimal() +
  theme(legend.position = 'top', 
        plot.title = element_text(face = "bold"),
        panel.grid.minor.x = element_blank())

#---------------------------------------
# Fig 4C
#---------------------------------------
boxplot_topo_size <- intersect_2haps %>%
  subset(Annot != ".") %>%
  mutate(ANNOT_SIZE = END_annot - START_annot + 1) %>%
  mutate(DIV_cat = factor(DIV, levels=c("0.1", "1", "5"))) %>%
  ggplot(aes(x = DIV_cat, y = ANNOT_SIZE, color = ANNOT_TOPO)) +
  geom_boxplot(position = position_dodge(preserve = "single"), alpha = 0.8) +
  scale_y_continuous(transform = "log10") +
  scale_color_manual(values=color_palette$topo) +
  facet_grid(. ~ PG) +
  theme_minimal() +
  labs(title = "C", 
       x = "Nucleotide polymorphism (%)", y = "Bubble size (bp)",
       color = "Bubble topology") +
  theme(legend.position = 'bottom', plot.title = element_text(face = "bold"))

#---------------------------------------
# Arrange subfigs into Fig 4
#---------------------------------------
legend_pg_color = get_legend(plot_AlnResc_div)
legend_topo_color = get_legend(boxplot_topo_size)

plot_AlnResc_div = plot_AlnResc_div + theme(legend.position = 'none')
plot_PathExp_div = plot_PathExp_div + theme(legend.position = 'none')
boxplot_topo_size = boxplot_topo_size + theme(legend.position = 'none')

pdf(
  file=paste0(MASTER_PATH, "Main/Figure4.pdf"),
  width=7, height=6
)
grid.arrange(
  #1
  legend_pg_color,
  #2,3
  plot_AlnResc_div, plot_PathExp_div,
  #4
  boxplot_topo_size,
  #5
  legend_topo_color,
  ncol=2,
  layout_matrix = rbind(
    c(1,1),
    c(2,3),
    c(4,4),
    c(5,5)
  ),
  heights=c(
    0.2,
    1,
    1,
    0.1
  )
)
dev.off()

#-------------------------------------------------------------------------------
# Fig 5: Real data
#-------------------------------------------------------------------------------
#---------------------------------------
# Barplot Number of TP
#---------------------------------------
inv_hg_merged_all$INV_SIZE = factor(
  inv_hg_merged_all$INV_SIZE,
  levels=c("1 - 5 Kb", "5 Kb - 100 Kb", "> 100 Kb")
)

barplot_recall <- inv_hg_merged_all %>%
  subset(REDUNDANT == F) %>%
  group_by(PG, INV_SIZE) %>%
  summarize(Count = n()) %>%
  ungroup() %>%
  complete(INV_SIZE, PG, fill = list(Count = 0)) %>%
  ggplot(aes(x=INV_SIZE, y=Count, fill=PG)) +
  geom_rect(xmin=0.6, xmax=1.4, ymin=0, ymax=hgsvc2_chrX_merged[hgsvc2_chrX_merged$INV_SIZE == levels(hgsvc2_chrX_merged$INV_SIZE)[2],]$Count,
            fill='grey', alpha=0.05) +
  geom_rect(xmin=1.6, xmax=2.4, ymin=0, ymax=hgsvc2_chrX_merged[hgsvc2_chrX_merged$INV_SIZE == levels(hgsvc2_chrX_merged$INV_SIZE)[3],]$Count,
            fill='grey', alpha=0.05) +
  geom_rect(xmin=2.6, xmax=3.4, ymin=0, ymax=hgsvc2_chrX_merged[hgsvc2_chrX_merged$INV_SIZE == levels(hgsvc2_chrX_merged$INV_SIZE)[4],]$Count,
            fill='grey', alpha=0.05) +
  geom_bar(stat='identity', position="dodge", width=0.6) +
  # ylim(ymin=0, ymax=100) +
  labs(title="A", x="Inversion size category", y = "Annotated inversions (n)", fill = "Pangenome graph pipeline") +
  theme_minimal() +
  theme(legend.position = 'top', plot.title = element_text(face = "bold"),
        axis.title.x = element_blank())

# inv_hg_merged_count = inv_hg_merged[inv_hg_merged$REDUNDANT == F,] %>%
#   group_by(PG, INV_SIZE) %>%
#   summarize(Count = n())
# inv_hg_merged_count = left_join(inv_hg_merged_count, input_nINV_merged, by = "INV_SIZE")
# inv_hg_merged_count$Recall = inv_hg_merged_count$Count.x / inv_hg_merged_count$Count.y * 100

#---------------------------------------
# Venn diagram of TP inversions
#---------------------------------------
hg_merged_rev_intersect = rbind(
  # Chr7
  read.table(rev_intersect_hg_chr7_path[[names(rev_intersect_hg_chr7_path)[1]]], h=F, sep="\t") %>%
    mutate(PG = names(rev_intersect_hg_chr7_path)[1]),
  read.table(rev_intersect_hg_chr7_path[[names(rev_intersect_hg_chr7_path)[2]]], h=F, sep="\t") %>%
    mutate(PG = names(rev_intersect_hg_chr7_path)[2]),
  read.table(rev_intersect_hg_chr7_path[[names(rev_intersect_hg_chr7_path)[3]]], h=F, sep="\t") %>%
    mutate(PG = names(rev_intersect_hg_chr7_path)[3]),
  read.table(rev_intersect_hg_chr7_path[[names(rev_intersect_hg_chr7_path)[4]]], h=F, sep="\t") %>%
    mutate(PG = names(rev_intersect_hg_chr7_path)[4]),
  # ChrX
  read.table(rev_intersect_hg_chrX_path[[names(rev_intersect_hg_chrX_path)[1]]], h=F, sep="\t") %>%
    mutate(PG = names(rev_intersect_hg_chrX_path)[1]),
  read.table(rev_intersect_hg_chrX_path[[names(rev_intersect_hg_chrX_path)[2]]], h=F, sep="\t") %>%
    mutate(PG = names(rev_intersect_hg_chrX_path)[2]),
  read.table(rev_intersect_hg_chrX_path[[names(rev_intersect_hg_chrX_path)[3]]], h=F, sep="\t") %>%
    mutate(PG = names(rev_intersect_hg_chrX_path)[3]),
  read.table(rev_intersect_hg_chrX_path[[names(rev_intersect_hg_chrX_path)[4]]], h=F, sep="\t") %>%
    mutate(PG = names(rev_intersect_hg_chrX_path)[4])
) %>%
  select(-V5)

colnames(hg_merged_rev_intersect) = c("CHR", "START_annot", "END_annot", "Annot", "START_truth", "END_truth", "OVERLAP", "PG")
hg_merged_rev_intersect = add_redundancy_col(hg_merged_rev_intersect)

venn_PG_comp <- function(PG_annot_inv, subfig_title){
  filtered_PG_annot_inv = PG_annot_inv %>%
    subset(OVERLAP > 0 & REDUNDANT == F) %>%
    mutate(PG = case_when(
      PG == "Cactus" ~ "A",
      PG == "Minigraph" ~ "B",
      PG == "Minigraph-Cactus" ~ "C",
      PG == "PGGB" ~ "D"
    ))
  
  PG_comp = list()
  for (pg in unique(filtered_PG_annot_inv$PG)){
    PG_comp[[pg]] = subset(filtered_PG_annot_inv, PG == pg)$START_truth
  }
  
  return(PG_comp)
}

merged_PG_comp <- venn_PG_comp(hg_merged_rev_intersect)

venn_recall <- ggvenn(
  merged_PG_comp, 
  fill_color = scales::hue_pal()(4),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = FALSE
) + 
  labs(title="B") + 
  theme(plot.title = element_text(face = "bold"))

#---------------------------------------
# Barplot topology
#---------------------------------------
# Add topology column
inv_hg_merged_topo = inv_hg_merged_all %>%
  mutate(
    TOPO = case_when(
      grepl("INV:path", Annot) ~ "Path-explicit",
      grepl("INV:aln", Annot) ~ "Alignment-rescued"
    )
  ) %>%
  group_by(PG, TOPO) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  complete(PG, TOPO, fill = list(Count = 0))

inv_hg_merged_topo$TOPO = factor(
  inv_hg_merged_topo$TOPO,
  levels=c("Path-explicit", "Alignment-rescued")
)

barplot_topo <- inv_hg_merged_topo %>%
  ggplot(
    aes(
      x=PG, 
      y=count, 
      fill=TOPO
    )
  ) +
  geom_bar(
    stat="identity", 
    position="dodge",
    width=0.7
  ) +
  scale_fill_manual(
    values = color_palette$topo
  ) +
  labs(
    x="",
    y="Annotated inversions (n)",
    title="C",
    fill="Bubble topology"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

#---------------------------------------
# Arrange plots into Fig 5
#---------------------------------------
legend_pg <- get_legend(barplot_recall)
barplot_recall <- barplot_recall + theme(legend.position="none")
legend_topo <- get_legend(barplot_topo)
barplot_topo <- barplot_topo + theme(legend.position="none")

pdf(file=paste0(MASTER_PATH, "Main/Figure5.pdf"),
    width=7, height=6)
grid.arrange(
  #1
  legend_pg,
  #2,3
  barplot_recall, venn_recall,
  #4
  barplot_topo,
  #5
  legend_topo,
  ncol=2,
  layout_matrix = rbind(
    c(1,1),
    c(2,3),
    c(4,4),
    c(5,5)
  ),
  heights=c(
    0.2,
    1,
    1,
    0.1
  )
)
dev.off()

#===============================================================================
# SUPPLEMENTARY
#===============================================================================
#-------------------------------------------------------------------------------
# Sup Fig 1: Div parameter setting on recall
#-------------------------------------------------------------------------------
Sup_div <- ggplot(div_parameter_recall,
       aes(x = max_div,
           y = value,
           color = tool,
           shape = tool,
           group = tool)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(50, 100)) +
  # force the x ticks
  scale_x_continuous(
    breaks = c(0, 1, 2, 5, 10, 100),
    trans = pseudo_log_trans(base = 2)
  ) +
  labs(
    x = "d parameter value",
    y = "Recall (%)",
    color = "Pangenome graph pipeline",
    shape = "Pangenome graph pipeline"
  ) +
  # to split the legend on 2 rows
  guides(
    colour = guide_legend(nrow = 2),
    fill = guide_legend(nrow = 2)
  ) +
  theme_bw(base_size = 15) + 
  theme(legend.position = "bottom",
        panel.grid.minor.x = element_blank())

pdf(file="~/Documents/Paper/revisions/Final_figures/Supplementary/max_div_div1.pdf",
    width=7, height=5)
plot(Sup_div)
dev.off()

#-------------------------------------------------------------------------------
# Sup Fig 2: Cov parameter setting on recall
#-------------------------------------------------------------------------------
Sup_cov <- ggplot(cov_parameter_recall,
       aes(x = min_cov,
           y = value,
           color = tool,
           shape = tool,
           group = tool)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0, 100)) +
  # force the x ticks
  scale_x_continuous(
    breaks = c(0.1, 0.25, 0.50, 0.75, 0.9, 0.95, 0.99, 1)
  ) +
  labs(
    x = "cov parameter value",
    y = "Recall (%)",
    color = "Pangenome graph pipeline",
    shape = "Pangenome graph pipeline"
  ) +
  # to split the legend on 2 rows
  guides(
    colour = guide_legend(nrow = 2),
    fill = guide_legend(nrow = 2)
  ) +
  theme_bw(base_size = 15) + 
  theme(legend.position = "bottom",
        panel.grid.minor.x = element_blank())

pdf(file="~/Documents/Paper/revisions/Final_figures/Supplementary/min_cov_div1.pdf",
    width=7, height=5)
plot(Sup_cov)
dev.off()

#-------------------------------------------------------------------------------
# Sup Fig 3: Overlap parameter setting on recall
#-------------------------------------------------------------------------------
# Import results
results = read.table(
  file="~/Documents/Paper/revisions/Final_figures/Supplementary/overlap_parameter/CHM13_2hap_div1_recall_stats.txt",
  header = F,
  sep = " "
)

# Format table
pg = list(
  "cactus" = "Cactus",
  "minigraph" = "Minigraph",
  "mgc" = "Minigraph-Cactus",
  "pggb" = "PGGB"
)

results_toPlot = results %>%
  select(filename = V1, recall = V3, -V2) %>%
  mutate(
    pg = case_when(
      grepl(names(pg)[1], filename) ~ pg[[names(pg)[1]]],
      grepl(names(pg)[2], filename) ~ pg[[names(pg)[2]]],
      grepl(names(pg)[3], filename) ~ pg[[names(pg)[3]]],
      grepl(names(pg)[4], filename) ~ pg[[names(pg)[4]]]
    ),
    over = case_when(
      grepl("25.bed", filename) ~ 25,
      grepl("50.bed", filename) ~ 50,
      grepl("75.bed", filename) ~ 75,
      grepl("90.bed", filename) ~ 90
    )
  )

# Set pg order
results_toPlot$pg = factor(
  results_toPlot$pg,
  levels = c("Cactus", "Minigraph", "Minigraph-Cactus", "PGGB")
)

# Plot
Sup_overlap = ggplot(results_toPlot, aes(x=over, y=recall, group=pg)) +
  geom_line(aes(color=pg), linewidth=0.9, alpha=0.6)+
  geom_point(aes(color=pg, shape=pg), size=3, alpha=1) +
  ylim(ymin=50, ymax=100) +
  # force the x ticks to be at 25, 50, 75, 90
  scale_x_continuous(breaks = c(25, 50, 75, 90)) +
  labs(x="Minimum reciprocal overlap (%)", y = "Recall (%)", color = "Pangenome graph pipeline", shape = "Pangenome graph pipeline") +
  theme_bw(base_size = 15) +
  # to split the legend on 2 rows
  guides(
    colour = guide_legend(nrow = 2),
    fill = guide_legend(nrow = 2)
  ) +
  theme(legend.position = 'bottom',
        panel.grid.minor.x = element_blank())

pdf(file="~/Documents/Paper/revisions/Final_figures/Supplementary/Sup_overlap_parameters.pdf",
    width=7, height=5)
plot(Sup_overlap)
dev.off()


