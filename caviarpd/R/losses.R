#' Accuracy for a List of Clustering Estimates through Loss Functions
#'
#' Computes VI and Binder loss for any number of given clustering estimates. Estimates must be passed to function in list form.
#'
#' @param estimate.list An object of type 'list' containing any number of clustering estimates. Estimates can be from hierarchical clustering, caviarPD, or another method, but must be in list form.
#' @param truth The true partition of the data against which each estimate is compared, in numeric form.
#' @param order.by Since loss is computed for each estimate under both binder and VI, this argument specifies which loss function values the estimates are sorted by in the returned data frame.
#'
#' @return A data.frame object which contains each estimate's name (from the original list), number of clusters, VI loss, and Binder loss.
#'
#' @examples
#' iris.dis <- dist(iris[,-5])
#' iris.truth <- as.numeric(iris[,5])
#' est1 <- caviarPD(distance=iris.dis, mass=0.9, nSamples=1000, loss='binder')
#' est2 <- caviarPD(distance=iris.dis, mass=2.0, loss="VI")
#' loss.indexes(list(est1, est2), iris.truth)
#'
#' @export
#'
loss.indexes <- function(estimate.list, truth, order.by='binder') {
  if (class(estimate.list) != 'list') stop(" estimates must be in list form ")
  else if ( length(unique(sapply(estimate.list, length))) != 1 ) stop (" all estimates in estimate.list must be of equal length ")
  n <- length(estimate.list)
  mat <- matrix(0, nrow=n, ncol=3)
  for (i in 1:n) {
    mat[i,1] <- length(unique(estimate.list[[i]]))
    binder.loss <- binder(estimate.list[[i]], truth)
    vi.loss <- VI(estimate.list[[i]], truth)
    mat[i,2:3] <- c(binder.loss, vi.loss)
  }
  colnames(mat) <- c("N-Clusters", "Binder", "VI")
  rownames(mat) <- names(estimate.list)
  if (order.by == 'binder') {
    mat <- mat[order(mat[,2],decreasing=FALSE),]
  } else if (order.by == 'VI') {
    mat <- mat[order(mat[,3],decreasing=FALSE),]
  } else {
    stop(" must order results by either binder or VI loss ")
  }
  mat
}

