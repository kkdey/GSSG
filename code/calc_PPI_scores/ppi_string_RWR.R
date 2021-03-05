#' @param gene_scores  matrix of gene scores with genes along the rows and gene sets along the columns
#' @param restart_prob The probability of restart for RWR method (default is 0.5).
#' @param thresh_score An upper threshold for the scores in the gene score matrix.
#' @param edgelist_cutoff The cut-off on the edge connection strength. The default is 0.10 but the user is recommended to try other
#'                        options
#' @param NUM_RUNS The number of RWR iterations with random initialization based on gene scores.
#' @param NUM_SUBGENES The number of genes to keep in each sub-sample per run for a gene score. If NULL, we take 50% of genes per gene
#'                     score, thresholded below at 100.
#' @results:  A list with two entries : 'score' : The prioritization score between 0 and 1 based on the PPI RWR. We recommend thresholding it
#' to top X% genes. 'score_sd': Standard deviation of the score.
#' @import dnet, igraph, STRINGdb
#'
#'


ppi_string_RWR = function(gene_scores,
                   restart_prob=0.5,
                   thresh_score=5,
                   edgelist_cutoff = 0.10,
                   NUM_RUNS = 5,
                   NUM_SUBGENES=NULL){

  library(dnet)
  library(igraph)
  library(STRINGdb)

  if(is.null(rownames(gene_scores))){
    warning("The gene_scores matrix row names were not provided. We assign them indices that may not match the entries in gg_edgelist
            genes.")
    rownames(gene_scores) = 1:nrow(gene_scores)
  }
  if(is.null(colnames(gene_scores))){
    warning("The names of gene sets in columns of gene_scores not provided. We assign them names V1, V2....")
    colnames(gene_scores) = paste0("V", 1:ncol(gene_scores))
  }


  #############################  Create a String database graph ##############################################
  string_db <- STRINGdb$new( version="11", species=9606,
                             score_threshold=400, input_directory="" )
  string_graph = string_db$get_graph()

  #############################  Get all gene aliases for all proteins (may not be in graph)  #####################

  gene_aliases = string_db$get_aliases(V(string_graph))

  #####################  Subset the gene aliases to proteins that occur in the graph  #######################

  idx = which(!is.na(match(gene_aliases$STRING_id, as_ids(V(string_graph)))))
  gene_aliases_2 = gene_aliases[idx, ]

  common_genes = intersect(rownames(gene_scores), gene_aliases_2$alias)
  gene_scores2 = gene_scores[match(common_genes, rownames(gene_scores)), ]

  if(nrow(gene_scores2) < 100){
    stop("The number of genes matched in the gene scores matrix is too few (<100). Either check for match discepancy between
         gg_edgelist and gene_scores input genes or include more genes in the gene scores file.")
  }

  if(is.null(NUM_SUBGENES)){
    numgenes1 = floor(nrow(gene_scores2)*0.5)
    if(numgenes1 < 100){
      NUM_SUBGENES = 100
    }else{
      NUM_SUBGENES = numgenes1
    }
  }

  priority_list = list()

  for(num_iter in 1:NUM_RUNS){

    #############  Create seed genes for RWR  method    ###################################

    ll = list()
    for(cc in 1:ncol(gene_scores2)){
      tmp = gene_scores2[,cc]
      tmp[tmp > thresh_score] = thresh_score
      xx = (tmp+0.05)
      probx = xx/sum(xx)
      ll[[cc]] = sample(rownames(gene_scores2), NUM_SUBGENES, prob = probx, replace = T)
    }

    seed_genes = Reduce(union, ll)
    all_seeds = matrix(0, length(seed_genes), length(ll))
    for(cc in 1:ncol(all_seeds)){
      all_seeds[match(ll[[cc]], seed_genes),cc] = 1
    }
    rownames(all_seeds) = seed_genes
    colnames(all_seeds) = colnames(gene_scores2)

    common_genes = intersect(rownames(all_seeds), gene_aliases_2$alias)
    gene_aliases_3 = gene_aliases_2[which(!is.na(match(gene_aliases_2$alias, common_genes))),]
    all_seeds2 = all_seeds[match(common_genes, rownames(all_seeds)),]

    tmp_names = gene_aliases_3[match(rownames(all_seeds2), gene_aliases_3[,2]), 1]

    mods_df = c()
    for(mm in 1:ncol(all_seeds2)){
      tmp = tapply(all_seeds2[,mm], tmp_names, mean)
      mods_df = cbind(mods_df, as.numeric(tmp))
      cats = names(tmp)
    }
    rownames(mods_df) = cats

    all_seeds3 = mods_df


    PTmatrix <- dRWR(g=string_graph, normalise="laplacian", setSeeds=all_seeds3,
                     restart=0.5, parallel=TRUE)
    rownames(PTmatrix) = as_ids(V(string_graph))

    Step1_matrix = apply(PTmatrix, 2, function(x) {
      y = 1-ecdf(x)(x)
      return(y)
    })

    Step2_vec = apply(Step1_matrix, 1, function(x) {
      y = -2*sum(log(x+1e-08))
      return(y)
    })

    combined_p = pchisq(Step2_vec, 2*ncol(Step1_matrix), lower.tail = F)
    unscaled_score = -log(combined_p)
    PR = (unscaled_score - min(unscaled_score))/(max(unscaled_score) - min(unscaled_score))
    names(PR) = as_ids(V(string_graph))

    ids = match(rownames(all_seeds3), names(PR))
    ids = ids[!is.na(ids)]
    PR2 = PR[ids]
    PR_gene_names = gene_aliases_3[match(names(PR2), gene_aliases_3[,1]), 2]

    names(PR2) = PR_gene_names

    priority_list[[num_iter]] = PR2
    cat("We have performed RWR for Iteration:", num_iter, "\n")
  }

  union_genes = c()
  for(mm in 1:length(priority_list)){
    union_genes = c(union_genes, names(priority_list[[mm]]))
  }
  union_genes = unique(union_genes)

  priority_matt = matrix(0, length(union_genes), NUM_RUNS)
  for(mm in NUM_RUNS){
    priority_matt[match(names(priority_list[[mm]]), union_genes), mm] = as.numeric(priority_list[[mm]])
  }
  score = rowMeans(priority_matt)
  names(score) = union_genes
  score_sd = apply(priority_matt, 1, sd)
  outlist = list("score" = score, "score_sd" = score_sd)
  return(outlist)
}
