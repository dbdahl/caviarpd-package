#' Cluster Analysis via Random Partition Distributions
#'
#' Returns a clustering estimate given pairwise distances. Users can also specify parameters for the mass, discount, temperature, loss function, or number of samples.
#' @param distance A pairwise distance matrix of class 'dist'.
#' @param mass The main tuning parameter; higher mass leads to more clusters.
#' @param temperature A positive number that accentuates or dampens distance between observations.
#' @param discount Typically 0, controls the distribution of subset sizes.
#' @param loss The salso method aims to estimate this loss function when searching the partition space for an optimal estimate, must be specified as either "binder" or "VI".
#' @param nSamples The number of samples used to estimate the loss function in the salso method.
#' @param maxSize Restriction parameter for the maximum number of clusters in the clustering estimate.
#' @param samplesOnly If TRUE, returns only the samples generated for a given mass, temperature, and discount rather than an actual clustering estimate.
#' @param distr The random partition distribution used to generate samples, must be specified as either "EPA" or "ddCRP".
#'
#' @return A list containing two elements: the clustering estimate and the summary. The summary contains the pairwise probabilities for all samples and is used to create the confidence plots.
#'
#' @examples
#' iris.dis <- dist(iris[,-5])
#' caviarPD(distance=iris.dis)
#' caviarPD(distance=iris.dis, mass=0.75, loss="binder", nSamples=10, maxSize=3)
#'
#' @export
#' @importFrom salso salso binder VI
#' @importFrom pumpkin EPAPartition DDCRPPartition samplePartition
#'
caviarPD <- function(distance, temperature=10.0, mass=1.0, discount=0.0, loss="VI",
                     nSamples=100, maxSize=0, samplesOnly=FALSE, distr="EPA") {

  ### ERROR CHECKING ###
  if (class(distance) != 'dist') stop(" 'distance' argument must be an object of class 'dist' ")
  else if(class(mass) != 'numeric' ) stop(" 'mass' argument must be an object of class 'numeric' ")
  else if(class(nSamples) != 'numeric' ) stop(" 'nSamples' argument must be an object of class 'numeric' ")
  else if (loss != "VI" & loss != 'binder') stop(" 'loss' argument must be specified as either
                                                'binder' or 'VI' ")
  else if (samplesOnly != TRUE & samplesOnly != FALSE) stop(" 'samplesOnly' argument is not interpretable as logical ")
  else if(nSamples <= 0 | nSamples %% 1 !=0) stop (" must specify a positive integer for the 'nsamples' argument ")


  # Is temp restricted to be positive? What is discount restricted to?

  #else if(class(temperature) != 'numeric' ) stop(" 'temperature' argument must be an object of class 'numeric' ")
  #else if(class(discount) != 'numeric' ) stop(" 'discount' argument must be an object of class 'numeric' ")
  #else if(class(maxSize) != 'numeric' ) stop(" 'maxSize' argument must be an object of class 'numeric' ")


  similarity <- exp( -temperature * as.matrix(distance) )
  if (distr=="EPA") {
    distr <- EPAPartition(similarity=similarity, mass=mass, discount=discount, permutation=seq_len(nrow(similarity)))
  } else if (distr=="ddCRP") {
    distr <- DDCRPPartition(similarity=similarity, mass=mass)
  } else {
    stop("partition distribution must be specified as either 'EPA' or 'ddCRP' in the 'distr' argument ")
  }
  samples <- samplePartition(distr, nSamples, randomizePermutation=TRUE)
  if(samplesOnly == TRUE) return(samples)
  estimate <- salso(samples, loss=loss, maxNClusters = maxSize)
  list(estimate=estimate, summary=summary(estimate))
}
