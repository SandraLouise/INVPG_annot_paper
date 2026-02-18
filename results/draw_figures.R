#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-e", "--evalFile"), type="character", default=NULL, 
              help="eval file name", metavar="character"),
  make_option(c("-i", "--intersectFile"), type="character", default=NULL, 
              help="intersect file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="./", 
              help="output directory name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

library(gridExtra)
library(ggplot2)

results = read.table(opt$evalFile, h=T, sep="\t")
intersect = read.table(opt$intersectFile, h=F, sep="\t")

# FORMATING

results[results$PG == "cactus",]$PG = "Cactus" 
results[results$PG == "minigraph",]$PG = "Minigraph" 
results[results$PG == "mgc",]$PG = "Minigraph-Cactus" 
results[results$PG == "pggb",]$PG = "PGGB" 
results$PG = factor(results$PG, levels = c("Cactus", "Minigraph", "Minigraph-Cactus", "PGGB"))
results$DIV = as.factor(results$DIV)
# results$Stat_value[results$Stat_name == ""] = as.factor(results$SNP_pct)

colnames(intersect) = c("PG","DIV","TRUTH_START","TRUTH_END","ANNOT_START","ANNOT_END","ANNOT")
intersect[intersect$PG == "cactus",]$PG = "Cactus" 
intersect[intersect$PG == "minigraph",]$PG = "Minigraph" 
intersect[intersect$PG == "mgc",]$PG = "Minigraph-Cactus" 
intersect[intersect$PG == "pggb",]$PG = "PGGB" 
intersect$ANNOT_TOPO = "Unannotated"
for (i in 1:nrow(intersect)){
  if (grepl("INV:path", intersect[i,]$ANNOT, fixed = TRUE)){
    intersect[i,]$ANNOT_TOPO = 'Path-explicit'
  }
  else if (grepl("INV:aln", intersect[i,]$ANNOT, fixed = TRUE)){
    intersect[i,]$ANNOT_TOPO = 'Alignment-rescued'
  }
}
intersect$TRUTH_SIZE = intersect$TRUTH_END - intersect$TRUTH_START + 1

# PREPARE
sub_fig3A = results[results$Stat_name == "Recall",]
sub_fig3B = results[results$Stat_name == "Path-explicit" | results$Stat_name == "Alignment-rescued",]
sub_fig3B$Stat_name = factor(sub_fig3B$Stat_name, levels = c("Path-explicit", "Alignment-rescued"))
sub_fig3C = intersect
sub_fig3C$PG = factor(sub_fig3C$PG, levels = c("Cactus", "Minigraph", "Minigraph-Cactus", "PGGB"))
sub_fig3C$DIV = as.factor(sub_fig3C$DIV)
sub_fig3C$ANNOT_TOPO = factor(sub_fig3C$ANNOT_TOPO, levels = c("Path-explicit", "Alignment-rescued", "Unannotated"))

sub_fig4 = intersect[intersect$ANNOT_TOPO!="Unannotated",]
sub_fig4$START_OFFSET = abs(sub_fig4$ANNOT_START - sub_fig4$TRUTH_START)
sub_fig4$END_OFFSET = abs(sub_fig4$ANNOT_END - sub_fig4$TRUTH_END)
sub_fig4[sub_fig4$START_OFFSET == 0,]$START_OFFSET = 1
sub_fig4[sub_fig4$END_OFFSET == 0,]$END_OFFSET = 1
combined_offset <- data.frame(matrix(0, ncol = 4, nrow = nrow(sub_fig4)*2))
colnames(combined_offset) = c("PG", "DIV", "ANNOT_TOPO", "OFFSET")
combined_offset[1:nrow(sub_fig4),]$PG = sub_fig4$PG
combined_offset[(nrow(sub_fig4)+1):nrow(combined_offset),]$PG = sub_fig4$PG
combined_offset[1:nrow(sub_fig4),]$DIV = sub_fig4$DIV
combined_offset[(nrow(sub_fig4)+1):nrow(combined_offset),]$DIV = sub_fig4$DIV
combined_offset[1:nrow(sub_fig4),]$ANNOT_TOPO = sub_fig4$ANNOT_TOPO
combined_offset[(nrow(sub_fig4)+1):nrow(combined_offset),]$ANNOT_TOPO = sub_fig4$ANNOT_TOPO
combined_offset[1:nrow(sub_fig4),]$OFFSET = sub_fig4$START_OFFSET
combined_offset[(nrow(sub_fig4)+1):nrow(combined_offset),]$OFFSET = sub_fig4$END_OFFSET
combined_offset$PG = factor(combined_offset$PG, levels = c("Cactus", "Minigraph", "Minigraph-Cactus", "PGGB"))
combined_offset$DIV = as.factor(combined_offset$DIV)
combined_offset$ANNOT_TOPO = factor(combined_offset$ANNOT_TOPO, levels = c("Path-explicit", "Alignment-rescued"))

# FIGURE 3

# Recall
plot1 = ggplot(data=sub_fig3A, aes(x=DIV, y=Stat_value)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.6, fill="#1F78B4") +
  ylim(ymin=0, ymax=100) +
  facet_grid(~ PG) +
  labs(title="A", x="SNP %", y = "Recall (%)") +
  theme(legend.position="right")
  
# Topology count
plot2 = ggplot(data=sub_fig3B, aes(x=DIV, y=Stat_value, fill=Stat_name)) +
  geom_bar(stat="identity", position="stack", width=0.6, ) +
  ylim(ymin=0, ymax=100) +
  facet_grid(~ PG) +
  labs(title="B", x="SNP %", y = "# inversion annotations", fill = "Annotation signal type") +
  scale_color_manual(values=c("#D1495B","#EDAE49"), aesthetics = "fill") +
  theme(legend.position="top")
  
# Size vs. topology
plot3 = ggplot(data=sub_fig3C, aes(x=DIV, y=TRUTH_SIZE, color=ANNOT_TOPO)) +
  geom_boxplot(width=0.6) +
  scale_y_continuous(trans='log10') +
  # scale_color_manual(values=c("#d1495b","#edae49","#999999"))+
  scale_color_manual(values=c("#D1495B","#EDAE49","darkgrey"))+ # for coeno data
  facet_grid(~ PG) +
  labs(title="C", x="SNP %", y = "Inversion size (nt)", color = "Inversion category") +
  theme(legend.position="top")

# grid.arrange(plot1, plot2, plot3, ncol=1)

grob <- arrangeGrob(plot1, plot2, plot3, ncol = 1)
ggsave(paste(opt$out, "figure_3.pdf", sep="/"), grob, width = 7, height = 8)
ggsave(paste(opt$out, "figure_3.png", sep="/"), grob, width = 7, height = 8)

# FIGURE 4
ggplot(data=combined_offset, aes(x=DIV, y=OFFSET, color=ANNOT_TOPO)) +
  geom_boxplot(width=0.6) +
  scale_y_continuous(trans='log10') +
  scale_color_manual(values=c("#D1495B","#EDAE49")) +
  facet_grid(~ PG) +
  labs(x="SNP %", y = "Breakpoint position offset (nt)", color = "Annotation signal type") +
  theme(legend.position="bottom")

ggsave(paste(opt$out, "figure_4.pdf", sep="/"), width = 7, height = 4)
ggsave(paste(opt$out, "figure_4.png", sep="/"), width = 7, height = 4)
