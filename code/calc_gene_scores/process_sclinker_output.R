suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--inputcell", type="character", default="sclinker/", help="Directory where you have sclinker output files"),
  make_option("--prefix", type="character", default="Alzheimers", help="Name of the file prefix"),
  make_option("--outcell", type="character", default="sclinker_genescores/", help="Output directory with gene program files")
)

opt <- parse_args(OptionParser(option_list=option_list))
dput(opt)

suppressPackageStartupMessages(library(data.table))


filename=paste0(opt$inputcell, "/", opt$prefix, "_score.csv")
pfilename = paste0(opt$inputcell, "/", opt$prefix, "_pval.csv")
df2 = data.frame(fread(filename))
pdf2 = data.frame(fread(pfilename))
temp = df2[,-1]
temp[temp < 0] = 0

pval_logfold2 = apply(temp, 2, function(x) return(2*pnorm(abs(x), 0, 1, lower.tail = F)))

qq_logfold2 = apply(pval_logfold2, 2, function(y){
  z = -2*log(y+1e-08)
  PR = (z - min(z))/(max(z) - min(z))
  return(PR)
})

if(!dir.exists(paste0(opt$outcell, "/", opt$prefix))){
  dir.create(paste0(opt$outcell, "/", opt$prefix))
}

for(mm in 1:ncol(qq_logfold2)){
  df = cbind.data.frame(df2[,1], qq_logfold2[, mm])
  write.table(df, file = paste0(opt$outcell, "/", opt$prefix,  "/",
                                colnames(qq_logfold2)[mm], ".txt"), row.names = F, col.names = F, sep = "\t", quote=F)

}

genes = df2[,1]
tab = cbind.data.frame(genes, 1)
write.table(tab, file = paste0(opt$outcell, "/", opt$prefix,  "/", "ALL.txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)

