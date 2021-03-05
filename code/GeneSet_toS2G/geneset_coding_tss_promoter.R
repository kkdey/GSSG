options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
annotation_cell <- toString(args[1])
annotation_folder <- toString(args[2])
baseline_cell <- toString(args[3])


library(data.table)

for(numchr in 1:22){
  df = data.frame(fread(paste0("zcat ", annotation_cell, "/", annotation_folder, "/",
                               "5kb/", "5kb", ".", numchr, ".annot.gz")))
  base = data.frame(fread(paste0("zcat ", baseline_cell, "/",
                                 "baselineLD_v2.1", "/", "baselineLD.", numchr, ".annot.gz")))

  newdf1 = cbind.data.frame(df[,1:4], df[,5]*base$TSS_Hoffman)
  colnames(newdf1) = c(colnames(base)[1:4], "AN")
  newdf2 = cbind.data.frame(df[,1:4], df[,5]*base$Promoter_UCSC)
  colnames(newdf2) = c(colnames(base)[1:4], "AN")
  newdf3 = cbind.data.frame(df[,1:4], df[,5]*base$Coding_UCSC)
  colnames(newdf3) = c(colnames(base)[1:4], "AN")

  if(!dir.exists(paste0(annotation_cell, "/", annotation_folder, "/", "TSS_Hoffman"))){
    dir.create(paste0(annotation_cell, "/", annotation_folder, "/", "TSS_Hoffman"))
  }

  write.table(newdf1,
              file = gzfile(paste0(annotation_cell, "/", annotation_folder, "/",
                                   "TSS_Hoffman/",
                                   "TSS_Hoffman.", numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)


  if(!dir.exists(paste0(annotation_cell, "/", annotation_folder, "/", "PROMOTER"))){
    dir.create(paste0(annotation_cell, "/", annotation_folder, "/", "PROMOTER"))
  }

  write.table(newdf2,
              file = gzfile(paste0(annotation_cell, "/", annotation_folder, "/",
                                   "PROMOTER/",
                                   "PROMOTER.", numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)

  if(!dir.exists(paste0(annotation_cell, "/", annotation_folder, "/", "CODING"))){
    dir.create(paste0(annotation_cell, "/", annotation_folder, "/", "CODING"))
  }

  write.table(newdf3,
              file = gzfile(paste0(annotation_cell, "/", annotation_folder, "/",
                                   "CODING/",
                                   "CODING.", numchr, ".annot.gz")),
              quote=FALSE, row.names=FALSE)

  cat("We are at chr:", numchr, "\n")

}

