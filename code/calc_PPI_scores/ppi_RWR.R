#' @param gg_edgelist A data frame or matrix with 3 columns of the form (Gene1, Gene2, Weight)
#' @param gene_scores  matrix of gene scores with genes along the rows and gene sets along the columns
#' @param restart_prob The probability of restart for RWR method (default is 0.5).
#' @param edgelist_cutoff The cut-off on the edge connection strength. The default is 0.10 but the user is recommended to try other
#'                        options
#' @param NUM_RUNS The number of RWR iterations with random initialization based on gene scores.
#' @param NUM_SUBGENES The number of genes to keep in each sub-sample per run for a gene score. If NULL, we take 50% of genes per gene
#'                     score, thresholded below at 100.
#' @results:  A list with two entries : 'score' : The prioritization score between 0 and 1 based on the PPI RWR. We recommend thresholding it 
              to top X% genes. 'score_sd': Standard deviation of the score. 
#' @import dnet, igraph


ppi_RWR = function(gg_edgelist,
                   gene_scores,
                   restart_prob=0.5,
                   edgelist_cutoff = 0.10,
                   NUM_RUNS = 5,
                   NUM_SUBGENES=NULL){

  library(dnet)
  library(igraph)

  if(ncol(gg_edgelist) != 3){
    stop("The gg_edgelist input must be a data frame or matrix of the form (Gene1, Gene2, Weight). If weight not provided,
         add a third column with all 1.")
  }
  if(is.null(rownames(gene_scores))){
    warning("The gene_scores matrix row names were not provided. We assign them indices that may not match the entries in gg_edgelist
            genes.")
    rownames(gene_scores) = 1:nrow(gene_scores)
  }
  if(is.null(colnames(gene_scores))){
    warning("The names of gene sets in columns of gene_scores not provided. We assign them names V1, V2....")
    colnames(gene_scores) = paste0("V", 1:ncol(gene_scores))
  }


  ############### use cut-off to filter the edgelist matrix and convert it to igraph format  ###########################

  net = gg_edgelist
  net = net[which(net[,3] > edgelist_cutoff ),]
  gg = graph_from_edgelist(as.matrix(net[,1:2]), directed=F)
  graph_nodes = as_ids(V(gg))


  common_genes = intersect(rownames(gene_scores), graph_nodes)
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

    #####################  Perform Random Walk with Restart   ########################################


    PTmatrix <- dRWR(g=gg, normalise="laplacian", setSeeds=all_seeds,
                     restart=0.5, parallel=TRUE)
    rownames(PTmatrix) = graph_nodes

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
    names(PR) = as_ids(V(gg))
    priority_list[[num_iter]] = PR
    cat("We have performed RWR for Iteration:", num_iter, "\n")
  }

  priority_matt = do.call(cbind, priority_list)
  score = rowMeans(priority_matt)
  score_sd = apply(priority_matt, 1, sd)

  outlist = list("score" = score, "score_sd" = score_sd)
  return(outlist)
}




####################################  Examples  #############################################

#gg_edgelist = get(load("~/Documents/Enhancer_MasterReg/data/GTEx_PANDA_tissues.RData"))
#gene_scores = get(load("/Users/kushaldey/Documents/Mouse_Humans/output/gene_scores_mat_trial.rda"))
#res = ppi_RWR(gg_edgelist, gene_scores)
