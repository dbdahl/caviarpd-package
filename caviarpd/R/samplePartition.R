#' @useDynLib caviarpd .samplePartition
#'
samplePartition <- function(distr, nSamples, randomizePermutation=FALSE) {
  UseMethod("samplePartition")
}

samplePartition.default <- function(distr, nSamples, randomizePermutation=FALSE) {
  engine <- function(priorID, p) {
    .Call(.samplePartition, nSamples, distr$nItems, seed4rust(), priorID, p, randomizePermutation)
  }
  partitionDispatch(engine, distr, c("CenteredPartition","PowerPartition","PeggedTimeDependentPartition"))
}
