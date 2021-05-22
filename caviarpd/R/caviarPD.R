#' Cluster Analysis via Random Partition Distributions
#'
#' Returns a clustering estimate given pairwise distances. Users can also specify parameters for the mass, discount, temperature, loss function, or number of samples.
#'
#' @param distance A pairwise distance matrix of class 'dist'.
#' @param mass The main tuning parameter for the number of clusters. A higher mass value tends to lead to more clusters.
#' @param loss The SALSO method (Dahl, Johnson, Müller, 2021) tries to minimize this expected loss when searching the partition space for an optimal estimate. This must be either "binder" or "VI".
#' @param nSamples Number of samples used to estimate the partition confidence and variance ratio for each mass value obtained from the search algorithm; only applicable if single=TRUE.
#' @param samplesOnly If TRUE, the function only returns the samples generated for a given mass, temperature, and discount rather than an actual clustering estimate.
#' @param distr The random partition distribution used to generate samples.  This must be specified as either "EPA" or "ddCRP".
#' @param temperature A positive number that accentuates or dampens distance between observations.
#' @param discount A number in \eqn{[0,1)} giving the discount parameter to control the distribution of subset sizes.
#' @param maxNClusters The maximum number of clusters that can be considered by the SALSO method.
#' @param nCores The number of CPU cores to use. A value of zero indicates to use all cores on the system.
#'
#' @return A list containing two elements: the clustering estimate and the summary. The summary contains the pairwise probabilities for all samples and is used to create the confidence plots.
#'
#' @references
#'
#' D. B. Dahl, D. J. Johnson, and P. Müller (2021), Search Algorithms and Loss
#' Functions for Bayesian Clustering, <arXiv:2105.04451>.
#'
#' @examples
#' iris.dis <- dist(iris[,-5])
#' # In practice the user should use at least 100 samples, but for ease of testing we use less here.
#' caviarPD(distance=iris.dis, nSamples=10, nCores=1)
#' caviarPD(distance=iris.dis, mass=0.75, loss="binder", nSamples=10, maxNClusters=3, nCores=1)
#'
#' @export
#' @importFrom salso salso binder VI
#'
caviarPD <- function(distance, mass=1.0, loss="binder", nSamples=1000, samplesOnly=FALSE,
                     distr="EPA", temperature=10.0, discount=0.0, maxNClusters=0, nCores=0) {

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
    .Call(.sample_epa,nSamples, similarity, mass, discount, nCores, mkSeed())
  } else if (distr=="ddCRP") {
    samplePartition(DDCRPPartition(similarity=similarity, mass=mass), nSamples, randomizePermutation=TRUE)
  } else {
    stop("partition distribution must be specified as either 'EPA' or 'ddCRP' in the 'distr' argument ")
  }
  if(samplesOnly == TRUE) return(samples)
  salso(samples, loss=loss, maxNClusters = maxNClusters)
}
