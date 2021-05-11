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
#' @param nCores: The number of CPU cores to use. A value of zero indicates to use all cores on the system.
#'
#' @return A list containing two elements: the clustering estimate and the summary. The summary contains the pairwise probabilities for all samples and is used to create the confidence plots.
#'
#' @examples
#' iris.dis <- dist(iris[,-5])
#' caviarPD(distance=iris.dis, nSamples=10)
#' caviarPD(distance=iris.dis, mass=0.75, loss="binder", nSamples=10, maxSize=3)
#' # In practice the user should use at least 100 samples, but for ease of testing we use less here.
#'
#' @export
#' @importFrom salso salso binder VI
#'
caviarPD <- function(distance, temperature=10.0, mass=1.0, discount=0.0, loss="VI",
                     nSamples=100, maxSize=0, samplesOnly=FALSE, distr="EPA", nCores=0) {

  ### ERROR CHECKING ###
  if (class(distance) != 'dist') stop(" 'distance' argument must be an object of class 'dist' ")
  else if( !is.numeric(mass) ) stop(" 'mass' must be numeric and greater than -'discount' ")
  else if (loss != "VI" & loss != 'binder') stop(" 'loss' argument must be specified as either 'binder' or 'VI' ")
  else if (samplesOnly != TRUE & samplesOnly != FALSE) stop(" 'samplesOnly' argument is not interpretable as logical ")
  else if(nSamples <= 0 | nSamples %% 1 !=0 | !is.numeric(nSamples)) stop (" must specify a positive integer for the 'nsamples' argument ")
  else if (temperature < 0) stop(" 'temperature' must be nonnegarive ")
  else if(maxSize < 0 | maxSize %% 1 !=0 | !is.numeric(maxSize)) warning (" 'maxSize' argument should be a positive integer, ignoring constraint ")


  similarity <- exp( -temperature * as.matrix(distance) )
  if (distr=="EPA") {
    distr <- EPAPartition(similarity=similarity, mass=mass, discount=discount,
                          permutation=seq_len(nrow(similarity)))
    if(samplesOnly == TRUE) return(sample_epa(nSamples, similarity, mass, discount, 0))
    caviarpd(nSamples, similarity, mass, discount, loss=="VI", 16, maxSize, 0)
  } else if (distr=="ddCRP") {
    distr <- DDCRPPartition(similarity=similarity, mass=mass)
    samples <- samplePartition(distr, nSamples, randomizePermutation=TRUE)
    if(samplesOnly == TRUE) return(samples)
    salso(samples, loss=loss, maxNClusters = maxSize)
  } else {
    stop("partition distribution must be specified as either 'EPA' or 'ddCRP' in the 'distr' argument ")
  }
}
