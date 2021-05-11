seed4rust <- function() sapply(1:32, function(i) sample.int(256L,1L)-1L)

checkSimilarity <- function(similarity) {
  if ( ! is.matrix(similarity) ) stop("'similarity' must be a symmetric matrix of strictly positive enteries.")
  if ( ! isSymmetric(similarity) ) stop("'similarity' must be a symmetric matrix of strictly positive enteries.")
  if ( any( similarity <= 0 ) ) stop("'similarity' must be a symmetric matrix of strictly positive enteries.")
}

checkPermutation <- function(permutation) {
  if ( is.null(permutation) ) stop("'permutation' must be non-null.")
  if ( ! is.vector(permutation) ) stop("'permutation' must a vector.")
  if ( ! all(permutation %% 1 == 0) ) stop("'permutation' must only contain integers.")
  n <- length(permutation)
  if ( ( min(permutation) < 1 ) || ( max(permutation) > n ) || ( length(unique(permutation)) != n ) ) stop("'permutation' is not a valid permutation.")
}

checkMassDiscount <- function(mass, discount) {
  if ( ( discount < 0.0 ) || ( discount >= 1 ) ) stop("'discount' must be in [0,1).")
  if ( mass <= -discount ) stop("'mass' must be greater than -'discount'.")
}

canonicalForm <- function(labels) {
  temp <- integer(length(labels))
  i <- 1
  for (s in unique(labels)) {
    temp[which(labels == s)] <- i
    i <- i + 1
  }
  temp
}

