library(data.table)
library(R.utils)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
genescore_dir <- toString(args[1])
bed_dir <- toString(args[2])
annot_name <- toString(args[3])


score_file = paste0(genescore_dir, "/", annot_name, ".txt")
gene_scores = read.delim(score_file, header=F)

if(!dir.exists(paste0(bed_dir, "/", annot_name))){
  dir.create(paste0(bed_dir, "/", annot_name))
}

source("all_bedgraph_methods.R")

gene_scores = read.delim(score_file, header=F)

scores = gene_scores[,2]
names(scores) = gene_scores[,1]

out1 = ABC_intergenic_bedgraph_calc(scores,
                                    output_cell = paste0(bed_dir, "/", annot_name),
                                    output_bed = paste0("ABC.bed"))

out4 = Roadmap_Enhancer_bedgraph_calc(scores,
				    output_cell = paste0(bed_dir, "/", annot_name),
                                    output_bed = paste0("Roadmap_Enhancer.bed"))

out5 = PCHiC_bedgraph_calc(scores,
                           output_cell = paste0(bed_dir, "/", annot_name),
                           output_bed = "PCHiC.bed")

out6 = Yoshida_bedgraph_calc(scores,
                             output_cell = paste0(bed_dir, "/", annot_name),
                             output_bed = "Yoshida.bed")

out7 = eQTLCPP_bedgraph_calc(scores,
                             output_cell = paste0(bed_dir, "/", annot_name),
                             output_bed = "FinemapBloodeQTL.bed")


df = read.table("../../processed_data/Gene_100kb.txt", header=T)
df[which(df[,2] < 0), 2] = 0
score = gene_scores[match(df$gene, gene_scores[,1]), 2]
score[is.na(score)] = 0

temp = cbind.data.frame(df[,1:3], score)
temp2 = temp[which(temp$score != 0),]
write.table(temp2, file = paste0(bed_dir, "/", annot_name, "/", "100kb.bed"),
            quote=F, sep = "\t", row.names=F, col.names=F)


df = read.table("../../processed_data/Gene_5kb.txt", header=T)
df[which(df[,2] < 0), 2] = 0
score = gene_scores[match(df$gene, gene_scores[,1]), 2]
score[is.na(score)] = 0

temp = cbind.data.frame(df[,1:3], score)
temp2 = temp[which(temp$score != 0),]
write.table(temp2, file = paste0(bed_dir, "/", annot_name, "/", "5kb.bed"),
            quote=F, sep = "\t", row.names=F, col.names=F)



