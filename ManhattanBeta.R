manhattan.Beta <- function(dataframe, colors = c("gray10", "gray50"), ymax = "max", xaxis.cex = 1, yaxis.cex = 1, limitchromosomes = 1:23, suggestiveline = NULL, genomewideline = NULL, annotate=NULL, Title, ...) {
  
  d=dataframe
  ymax=max(d$Beta)*1.1
  ymin=min(d$Beta)
  
  #throws error if you don't have columns named CHR, BP, and P in your data frame.
  if (!("CHR" %in% names(d) & "BP" %in% names(d) & "Beta" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and Beta")
  
  # limits chromosomes to plot. (23=x, 24=y, 25=par?, 26=mito?)
  if (any(limitchromosomes)) d=d[d$CHR %in% limitchromosomes, ]
  
  # remove na's, sort by CHR and BP, and keep snps where 0<P<=1
  d = d[order(d$CHR, d$BP), ]
  
  # sets colors based on colors argument.
  colors <- rep(colors,max(d$CHR))[1:max(d$CHR)]
  
  # sets the maximum value on the y axis
  #if (ymax == "max") ymax<-ceiling(max(d$Beta))
  
  # creates continuous position markers for x axis for entire chromosome. also creates tick points.
  d$pos = NA
  ticks = NULL
  lastbase = 0
  numchroms = length(unique(d$CHR))
  for (i in unique(d$CHR)) {
    if (i == 1) {
      d[d$CHR == i, ]$pos = d[d$CHR == i, ]$BP
    } else {
      lastbase=lastbase+tail(subset(d,CHR==i-1)$BP, 1)
      d[d$CHR == i, ]$pos = d[d$CHR == i, ]$BP+lastbase
    }
    ticks=c(ticks, d[d$CHR == i, ]$pos[floor(length(d[d$CHR == i, ]$pos)/2)+1])
  }

  
  # create the plot
  # creates a blank plot
  with(d, plot(pos, Beta, ylim = c(0,ymax), ylab = expression("|" ~ beta ~ "|"), xlab = "Chromosome", xaxt = "n", type = "n", cex = 0.3, yaxt = "n", main = Title, ...))
  # then make an axis that has chromosome number instead of position
  axis(1, at = ticks, lab = unique(d$CHR), cex.axis = xaxis.cex)
  axis(2, cex.axis = yaxis.cex)
  icol=1
  for (i in unique(d$CHR)) {
    with(d[d$CHR==i, ],points(pos, Beta, col=colors[icol], cex=0.3, ...))
    icol = icol+1
  }
  
  # create a new data frame with rows from the original data frame where SNP is in annotate character vector.
  # then plot those points over the original graph, but with a larger point size and a different color.
  if (!is.null(annotate)) {
    d.annotate=d[which(d$SNP %in% annotate), ]
    icol=1
    for (i in unique(d.annotate$CHR)) {
      with(d.annotate[d.annotate$CHR==i, ], points(pos, Beta, col = "red", cex=0.5, pch = 20, ...))
      icol = icol+1
    }
  }
  
  # add threshold lines
  if (!is.null(suggestiveline)) abline(h=suggestiveline, col="blue")
  if (!is.null(genomewideline)) abline(h=genomewideline, col="red")
}
