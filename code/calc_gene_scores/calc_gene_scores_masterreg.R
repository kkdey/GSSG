
############################  Master-regulator gene score   ##################################
## REF: Vosa et al bioRxiv
## Unraveling the polygenic architecture of complex traits using blood eQTL metaanalysis

library(data.table)
cis_sig = data.frame(fread("~/GSSG/processed_data/cis-eQTL_significant_20181017.txt.gz"))
trans_sig = data.frame(fread("~/GSSG/processed_data/trans-eQTL_significant_20181017.txt.gz"))
common_snps = intersect(cis_sig$SNP, trans_sig$SNP)
gene_names = unique(cis_sig$GeneSymbol)

score = rep(0, length(gene_names))
for(num in 1:length(gene_names)){
  idx = which(cis_sig$GeneSymbol == gene_names[num])
  snp_ids = cis_sig$SNP[idx]
  idx2 = match(trans_sig$SNP, snp_ids)
  if(length(idx2[!is.na(idx2)]) == 0){
    score[num] = 0
  }else{
    trans_id = which(!is.na(idx2))
    score[num] = length(unique(trans_sig$GeneSymbol[trans_id]))
  }
  cat("We are at gene", num, "\n")
}

names(score) = gene_names

all_genes = read.table("~/GSSG/genesets/All_genes.txt",
                       header=F)
score2 = score[names(score) %in% all_genes[,1]]

final_genes = names(score2[which(score2 >= 3)])
df = cbind.data.frame(final_genes, 1)
write.table(df, file = "~/GSSG/genesets/master_regulator.txt",
            row.names=F, col.names=F, sep = "\t", quote=F)


#######################  Transcription Factor (TF) genes   #########################################################


library(readxl)
tf_data = read_excel("~/GSSG/processed_data/TF_genes.xlsx", sheet = 2)
tf_data2 = data.frame(tf_data[-1, ])
genes = tf_data2[,2]
genes2 = tf_data2[which(tf_data2[,4] == "Yes"),2]

df = cbind.data.frame(genes2, 1)
write.table(df, file = "~/GSSG/genesets/TF_genes_curated.txt",
            row.names = F, col.names = F, sep = "\t", quote=F)
