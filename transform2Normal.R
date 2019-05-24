#This function was taken from the eQTL-BMA package developed by Tim Flutre

#Force normal distribution
transform2Normal <- function(mat, break.ties.rand=TRUE,
                             seed=8988){
  if(break.ties.rand && ! is.null(seed))
    set.seed(seed)
  
  mat.qn <- t(apply(mat, 2, function(exp.per.gene){
    if(break.ties.rand){
      idx <- sample(length(exp.per.gene))
      tmp <- qqnorm(exp.per.gene[idx], plot.it=FALSE)$x
      tmp[sort(idx, index.return=TRUE)$ix]
    } else
      qqnorm(exp.per.gene, plot.it=FALSE)$x
  }))
  colnames(mat.qn) <- colnames(mat)
  
  return(mat.qn)
}
