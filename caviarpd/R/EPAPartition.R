EPAPartition <- function(similarity, permutation, mass, discount=0) {
  checkSimilarity(similarity)
  nItems <- nrow(similarity)
  if ( nItems < 1 ) stop("The number of rows in 'similarity' must be at least one.")
  checkPermutation(permutation)
  if ( length(permutation) != nItems ) stop("The length of 'permutation' must equal the number of rows in 'similarity'.")
  checkMassDiscount(mass, discount)
  result <- list(nItems=nItems, similarity=similarity, permutation=permutation-1L, mass=mass, discount=discount)
  class(result) <- c("EPAPartition", "partitionDistribution")
  result
}

