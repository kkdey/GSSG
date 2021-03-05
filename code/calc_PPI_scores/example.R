
source("ppi_RWR.R")
gg_edgelist = get(load("GTEx_PANDA_tissues.RData"))
gene_scores = get(load("gene_scores_mat_trial.rda"))
res = ppi_RWR(gg_edgelist, gene_scores)
score = res$score ## gene prioritization