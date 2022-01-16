#' Get random values from permutation test,  the mean of which as baseline to calculate CRD score.
#'
#' @param input.list A list of variables containing scaled gene expression matrix, giving number of bins.
#' @param genes.dist.bins A matrix of giving expression bins within genes, which were binned according to their average expression across cells or samples.
#' @param b.sign A logical value to indicate the overlapped features between target features and the features derived from gene expressio matrix.
#' @param num.rounds  A integer value to indicate the permutation times of iteration; 1000 by default and 10000 will be better for reproducibility.
#'
#' @return Continuous numerical variables with values, named by cells or samples.
#'
#'
#' @examples
#'
#' \dontrun{
#' r.scores <- get_random_CRD(input.list, input.list$genes.dist.bins, b.sign, num.rounds = num.rounds)
#' }
#'


get_random_CRD <- function(input.list, genes.dist.bins, b.sign, num.rounds = 1000){
  sign.bins <- as.matrix(table(genes.dist.bins[b.sign]))
  q <-rownames(sign.bins)
  bg <- matrix(data = F, nrow = length(genes.dist.bins), ncol = num.rounds)
  for (i in 1:nrow(sign.bins)){
    num.genes <- sign.bins[i]
    if(num.genes > 0){
      idx <- which(is.element(genes.dist.bins, q[i]))
      for (j in 1:num.rounds){
        idxj <- sample(idx, num.genes)
        bg[idxj, j] <- T
      }
    }
  }
  r.scores <- apply(bg, 2, function(x)colMeans(input.list$expr.scaled[x,]))
  r.scores <- rowMeans(r.scores)
  return(r.scores)
}
