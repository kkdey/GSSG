
############################  ABC  gene score calculation (Genic + Intergenic)  ##################################

#### Input: processed ABCv9 scores taken from (https://www.engreitzlab.org/resources) for blood/immune cell-types only
###  We add a 150bp window surrounding the ABC region to correct for any biases in calling the Hi-C or enhancer peaks.
###

#REF:Fulco, C.P. et al, 2019. Activity-by-contact model of enhancer–promoter regulation from thousands of
#CRISPR perturbations. Nature genetics, 51(12), pp.1664-1669.

library(data.table)
df = data.frame(fread(paste0("~/GSSG/processed_data/ABCpaper_NasserFulcoEngreitz2020_Blood_AvgHiC.txt.gz")))
df2 = df[which(df$ABC.Score > 0.015), ]

tabb = t(xtabs(~class + TargetGene, data = df2))

all_genes = read.table("~/GSSG/genesets/All_genes.txt",
                       header=F)

# ABC9_I_top10_0_015

score = (tabb[,2])
score2 = score[names(score) %in% all_genes[,1]]
outdf = cbind.data.frame(names(score2)[order(score2, decreasing = T)[1:2200]], 1)
write.table(outdf, file = "~/GSSG/genesets//ABC9_I_top10_0_015.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")

#ABC9_G_top10_0_015

score = (tabb[,2]+tabb[,1])
score2 = score[names(score) %in% all_genes[,1]]
outdf = cbind.data.frame(names(score2)[order(score2, decreasing = T)[1:2200]], 1)
write.table(outdf, file = "~/GSSG/genesets//ABC9_I_top10_0_015.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")


############################  Enhancer Domain Score (EDS)   ##################################
####  We compute top genes with highest blood specific Enhancer Domain Score (EDS): reference


library(data.table)

dat1 = data.frame(fread("~/GSSG/processed_data/EDSpaper_WangGoldstein2020_bytissue.txt"))
annot = cbind.data.frame(dat1$ENSG, dat1$BLOOD_activitylinking_cons)
geneanno = data.frame(data.table::fread("~/GSSG/processed_data/geneanno.csv"))
common_genes = intersect(geneanno$id, annot[,1])


annot2 = cbind.data.frame(geneanno[match(common_genes, geneanno$id), 2], annot[match(common_genes, annot[,1]), 2])
colnames(annot2) = c("Gene", "Score")

all_genes = read.table("~/GSSG/genesets/All_genes.txt",
                       header=F)
annot3 = annot2[annot2[,1] %in% all_genes[,1], ]

top_genes = as.vector(annot3[order(annot3[,2], decreasing = T)[1:2200], 1])
df = cbind(top_genes, 1)
write.table(df, file = "~/GSSG/genesets/EDS_Binary_top10.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")


############################  Promoter Capture HiC scores   ##################################
########## Immune cell-types Hi-C contacts  ########################
## REF:Javierre, B.M., et al., 2016. Lineage-specific genome architecture links enhancers
## and non-coding disease variants to target gene promoters. Cell, 167(5), pp.1369-1384.


pchic_data = data.frame(fread("~/GSSG/processed_data/PCHiC_peak_matrix_cutoff5.tsv", header=T))
names1 = lapply(pchic_data$baitName, function(x) return(strsplit(as.character(x), ";")[[1]]))
genes_tabb = table(unlist(names1))
top_genes = names(genes_tabb)[order(genes_tabb, decreasing = T)[1:2303]]
df = cbind.data.frame(top_genes, 1)
df = df[-(1:3), ]

df = read.table("~/GSSG/genesets/PCHiC_EDS.txt", header = T)

all_genes = read.table("~/GSSG/genesets/All_genes.txt",
                       header=F)
df2 = df[df[,1] %in% all_genes[,1], ]

top_genes = df2[order(df2[,2], decreasing = T)[1:2200], 1]
tab = cbind.data.frame(top_genes, 1)
write.table(tab, file = "~/GSSG/genesets/PCHiC_binary.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")


############################  Specifically Expressed Genes (Finucane et al 2018)   ##################################

## REF: Finucane, H.K. et al, 2018. Heritability enrichment of
## specifically expressed genes identifies disease-relevant tissues and cell types. Nature genetics, 50(4), pp.621-629.

seg_scores = read.delim("~/GSSG/processed_data/GTEx.tstat.tsv")
geneanno = data.frame(data.table::fread("~/GSSG/data/geneanno.csv"))
common_genes = intersect(seg_scores$ENSGID, geneanno$id)
annot2 = cbind.data.frame(geneanno[match(common_genes, geneanno$id), 2],
                          seg_scores[match(common_genes, seg_scores$ENSGID), "Whole_Blood"])
colnames(annot2) = c("Gene", "Score")

all_genes = read.table("~/GSSG/genesets/All_genes.txt",
                       header=F)
annot3 = annot2[annot2[,1] %in% all_genes[,1], ]

top_genes = as.vector(annot3[order(annot3[,2], decreasing = T)[1:2200], 1])
df = cbind(top_genes, 1)
write.table(df, file = "~/GSSG/genesets/SEG_GTEx_top10.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")

############################  ATAC-distal gene score   ##################################
########  From Yoshida et al 2019 Cell
## REF:Yoshida, H., et l. 2019.
## The cis-regulatory atlas of the mouse immune system. Cell, 176(4), pp.897-912.


temp = read.delim("~/GSSG/processed_data/HOMOD_HOMOP_grades2.txt", header = F)
all_genes = read.table("~/GSSG/genesets/All_genes.txt",
                       header=F)
temp2 = temp[temp[,1] %in% all_genes[,1], ]

genes = as.character(temp[order(temp[,2], decreasing = T)[1:2200], 1])
df = cbind(genes, 1)


write.table(df, file = "~/GSSG/genesets/HOMOD_top10.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")

write.table(temp2, file = "~/GSSG/genesets/HOMOD_prob.txt",
            row.names = F, col.names = F, quote = F, sep = "\t")


############################  Expecto-MVP gene score   ##################################
########  From Zhou et al 2018 Nature Genetics
## REF: Zhou, J. et al. 2018. Deep learning sequence-based ab initio
## prediction of variant effects on expression and disease risk. Nature genetics, 50(8), pp.1171-1179.


dat = data.frame(fread("~/GSSG/processed_data/variation_potential.positive_constraint.probability.txt"))

immune = c("V1", "Whole.Blood", "X4Star.hESC", "CD4.Memory.Primary.Cells", "CD4.Naive.Primary.Cells", "CD8.Naive.Primary.Cells",
           "GM12878...ENCODE", "Mobilized.CD34.Primary.Cells.Female")

dat2 = dat[, match(immune, colnames(dat))]
all_genes = read.table("~/GSSG/genesets/All_genes.txt",
                       header=F)
dat3 = dat2[dat2[,1] %in% all_genes[,1], ]

scores = apply(dat3[,-1], 1, max)
gene_names = dat3[order(scores, decreasing = T)[1:2200], 1]
df = cbind.data.frame(gene_names, 1)
write.table(df, file = "~/GSSG/genesets/Expecto_MVP.txt",
            row.names=F, col.names=F, quote=F, sep = "\t")


############################  eQTL-CTS gene score   ##################################
#################  DICE project ###############################

df = read.delim("~/GSSG/processed_data/DICE_PROP_TISSUE_SPECIFIC_QTL.txt", header=F)
genes = df[order(df[,2], decreasing = T)[1:2200], 1]
df2 = cbind.data.frame(genes, 1)

write.table(df2, file = "~/GSSG/genesets/eQTL_CTS_top10.txt",
            row.names=F, col.names=F, quote=F, sep = "\t")

all_genes = read.table("~/GSSG/genesets/All_genes.txt",
                       header=F)

df3 = df[df[,1] %in% all_genes[,1], ]
write.table(df3, file = "~/GSSG/genesets/eQTL_CTS_prob.txt",
            row.names=F, col.names=F, quote=F, sep = "\t")

