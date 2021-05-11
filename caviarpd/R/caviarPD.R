#' Cluster Analysis via Random Partition Distributions
#'
#' Returns a clustering estimate given pairwise distances. Users can also specify parameters for the mass, discount, temperature, loss function, or number of samples.
#' @param distance A pairwise distance matrix of class 'dist'.
#' @param mass The main tuning parameter; higher mass leads to more clusters.
#' @param temperature A positive number that accentuates or dampens distance between observations.
#' @param discount Typically 0, controls the distribution of subset sizes.
#' @param loss The salso method aims to estimate this loss function when searching the partition space for an optimal estimate, must be specified as either "binder" or "VI".
#' @param nSamples The number of samples used to estimate the loss function in the salso method.
#' @param maxNClusters Restriction parameter for the maximum number of clusters in the clustering estimate.
#' @param samplesOnly If TRUE, returns only the samples generated for a given mass, temperature, and discount rather than an actual clustering estimate.
#' @param distr The random partition distribution used to generate samples, must be specified as either "EPA" or "ddCRP".
#' @param nCores The number of CPU cores to use. A value of zero indicates to use all cores on the system.
#'
#' @return A list containing two elements: the clustering estimate and the summary. The summary contains the pairwise probabilities for all samples and is used to create the confidence plots.
#'
#' @examples
#' iris.dis <- dist(iris[,-5])
#' caviarPD(distance=iris.dis, nSamples=10)
#' caviarPD(distance=iris.dis, mass=0.75, loss="binder", nSamples=10, maxNClusters=3)
#' # In practice the user should use at least 100 samples, but for ease of testing we use less here.
#'
#' @export
#' @importFrom salso salso binder VI
#'
caviarPD <- function(distance, temperature=10.0, mass=1.0, discount=0.0, loss="VI",
                     nSamples=100, maxNClusters=0, samplesOnly=FALSE, distr="EPA", nCores=0) {

  ### ERROR CHECKING ###
  if (class(distance) != 'dist') stop(" 'distance' argument must be an object of class 'dist' ")
  else if( !is.numeric(mass) ) stop(" 'mass' must be numeric and greater than -'discount' ")
  else if (loss != "VI" & loss != 'binder') stop(" 'loss' argument must be specified as either 'binder' or 'VI' ")
  else if (samplesOnly != TRUE & samplesOnly != FALSE) stop(" 'samplesOnly' argument is not interpretable as logical ")
  else if(nSamples <= 0 | nSamples %% 1 !=0 | !is.numeric(nSamples)) stop (" must specify a positive integer for the 'nsamples' argument ")
  else if (temperature < 0) stop(" 'temperature' must be nonnegarive ")
  else if(maxNClusters < 0 | maxNClusters %% 1 !=0 | !is.numeric(maxNClusters)) warning (" 'maxNClusters' argument should be a positive integer, ignoring constraint ")


  similarity <- exp( -temperature * as.matrix(distance) )
  samples <- if (distr=="EPA") {
    sample_epa(nSamples, similarity, mass, discount, nCores)
  } else if (distr=="ddCRP") {
    samplePartition(DDCRPPartition(similarity=similarity, mass=mass), nSamples, randomizePermutation=TRUE)
  } else {
    stop("partition distribution must be specified as either 'EPA' or 'ddCRP' in the 'distr' argument ")
  }
  if(samplesOnly == TRUE) return(samples)
  salso(samples, loss=loss, maxNClusters = maxNClusters)
}
