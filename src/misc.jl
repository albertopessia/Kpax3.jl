function computecounts(d::Array{UInt8, 2},
                       partition::Array{Int, 1})
  f <- function(k) {
      return(as.integer(colSums(d[partition == k, , drop=FALSE])))
  }

  clust <- sort(unique(partition))
  m <- ncol(d)

  vapply(clust, FUN=f, FUN.VALUE=integer(length=m), USE.NAMES=FALSE)
end
