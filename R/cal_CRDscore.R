#'
#' Calculate Circadian rhythm disruption score using scRNA or bulk transcriptomic data
#'
#'
#'
#' @param expr log2 transformed TPM matrix derived from scRNAseq, or normalized gene expression matrix from bulk samples; rows are genes, cols are cells or samples.
#' @param circadians The character vectors containing selected circadian features.
#' @param study.type A string value to indicate the study type for calculation. Allowed values contain c('scRNAseq', 'bulk_RNAseq').
#' @param num.rounds A integer value to indicate the permutation times of iteration; 1000 by default and 10000 will be better for reproducibility.
#' @param n.bins A integer value giving expression bins within genes, which were binned according to their average expression across cells or samples. Set to 50 by default.
#' @param seed A logical value to indicate if set random seed for reproducible results. TRUE by default.
#'
#' @author Leihe, Yixian Fan.
#'
#' @return Continuous numerical variables with CRDscore, named by cells or samples.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("exprset")
#' data("circadians_selected")
#' CRDscore <- cal_CRDscore(expr = exprset, circadians = circadians_selected, study.type = "scRNAseq")
#' }



cal_CRDscore <- function(expr, circadians = NULL, study.type = NULL, num.rounds = 1000,  n.bins = 50, seed = TRUE){

  options(warn = -1)
  if(is.null(study.type)){stop("please set your study types as 'scRNAseq' or 'bulk_RNAseq'")}

  # create input data
  input.list <- list(expr = as.matrix(expr), genes = rownames(expr))

  # set study types
  # scRNAseq
  if(study.type == "scRNAseq"){
    input.list$genes.mean <- rowMeans(input.list$expr)
    input.list$expr.scaled <- sweep(input.list$expr, 1, input.list$genes.mean) #scaled expr by mean(genes)

    dist <- 10*((2^input.list$expr)-1)
    input.list$genes.dist <- log2(rowMeans(dist, na.rm = T) + 1)
  }

  # bulk_RNAseq
  if(study.type == "bulk_RNAseq"){
    input.list$genes.dist <- rowMeans(input.list$expr)
    input.list$expr.scaled <- sweep(input.list$expr, 1, input.list$genes.dist) #scaled expr by mean(genes)
  }

  #set bins, n.bins = 50
  input.list$genes.dist.bins <- arules::discretize(input.list$genes.dist, method = "frequency", breaks = n.bins) %>%
    plyr::mapvalues(levels(.), 1:length(levels(.))) %>%
    as.numeric() %>%
    as.matrix()

  # calculate CRDscore
  b.sign <- is.element(input.list$genes, circadians)

  if(sum(b.sign) > 5){
    if(seed){
      set.seed(123)
      r.scores <- get_random_CRD(input.list, input.list$genes.dist.bins, b.sign, num.rounds = num.rounds)
    }else{
      r.scores <- get_random_CRD(input.list, input.list$genes.dist.bins, b.sign, num.rounds = num.rounds)
    }

    raw.scores <- colMeans(input.list$expr.scaled[b.sign,])
    input.list$CRD.scores <- r.scores-raw.scores

    #input.list$CRD.scores.raw <- raw.scores
  }else{
    stop("Non-enough-overlapping genes to calculate CRDscore")
  }
  return(input.list$CRD.score)
}
