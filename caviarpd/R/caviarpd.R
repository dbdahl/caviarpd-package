#' Cluster Analysis via Random Partition Distributions
#'
#' Returns a clustering estimate given pairwise distances using the CaviarPD method.
#'
#' @param distance A pairwise distance matrix of class 'dist'.
#' @param nClusters A numeric vector giving the range for the number of clusters to consider when search for the mass parameter. Should be unset if the \code{mass} argument is used.
#' @param mass The main tuning parameter governing the number of clusters, with higher values tends toward more clusters. Should be unset if the \code{nClusters} argument is provided.
#' @param nSamples The number of samples used to generate the clustering estimate.
#' @param gridLength The length of the grid search for an optimal mass parameter. Only applicable if multiple values are provided for \code{nClusters}.
#' @param samplesOnly If TRUE, the function only returns the samples generated for a given mass, temperature, and discount rather than an actual clustering estimate.
#' @param loss The SALSO method (Dahl, Johnson, Müller, 2021) tries to minimize this expected loss when searching the partition space for an optimal estimate. This must be either "binder" or "VI".
#' @param distr The random partition distribution used to generate samples.  This must be specified as either "EPA" or "ddCRP".
#' @param temperature A positive number that accentuates or dampens distance between observations.
#' @param discount A number in \eqn{[0,1)} giving the discount parameter to control the distribution of subset sizes.
#' @param sd Number of standard deviations away from the expectation to be considered in finding boundary mass values.
#' @param maxNClusters The maximum number of clusters that can be considered by the SALSO method.
#' @param nCores The number of CPU cores to use. A value of zero indicates to use all cores on the system.
#'
#' @return A object of class \code{salso.estimate}, which provides a clustering estimate that can be displayed and plotted.
#'
#' @references
#'
#' D. B. Dahl, D. J. Johnson, and P. Müller (2021), Search Algorithms and Loss
#' Functions for Bayesian Clustering, [arXiv:2105.04451](https://arxiv.org/abs/2105.04451).
#'
#' @examples
#' iris.dis <- dist(iris[,-5])
#' # To reduce load on CRAN servers, limit the number of samples and CPU cores.
#' est <- caviarpd(distance=iris.dis, mass=1, nSamples=200, nCores=2)
#' est <- caviarpd(distance=iris.dis, nClusters=3:5, nSamples=200, nCores=2)
#' summ <- summary(est, orderingMethod=2)
#' plot(summ, type="heatmap")
#' plot(summ, type="mds")
#'
#' @export
#' @importFrom salso salso psm
#' @importFrom cluster silhouette
#' @importFrom stats uniroot
#'
caviarpd <- function(distance, nClusters, mass, nSamples=1000, gridLength=10, samplesOnly=FALSE,
                     loss="binder", distr="EPA", temperature=10.0, discount=0.0, sd=3, maxNClusters=0, nCores=0) {

  ### ERROR CHECKING ###
  if (class(distance) != 'dist') stop("'distance' argument must be an object of class 'dist'")
  if (loss != "VI" && loss != 'binder') stop("'loss' argument must be specified as either 'binder' or 'VI'")
  if (samplesOnly != TRUE && samplesOnly != FALSE) stop("'samplesOnly' argument is not interpretable as logical")
  if (nSamples <= 0 || nSamples %% 1 !=0 || !is.numeric(nSamples)) stop("must specify a positive integer for the 'nsamples' argument")
  if (temperature < 0) stop(" 'temperature' must be nonnegative ")
  if (maxNClusters < 0 || maxNClusters %% 1 !=0 || !is.numeric(maxNClusters)) warning("'maxNClusters' argument should be a positive integer, ignoring constraint")
  if (distr != "EPA" && distr != "ddCRP") stop("partition distribution must be specified as either 'EPA' or 'ddCRP' in the 'distr' argument")
  if ( !is.numeric(discount) || length(discount) != 1 || discount < 0 || discount >= 1.0 ) stop("'discount' must be in [0,1)")
  if (samplesOnly == TRUE && (missing(mass) || length(mass) > 1)) stop( "must specify single mass parameter in order to obtain samples" )
  if ( length(nSamples) == 1 ) {
    nSamplesSearch <- 0.5 * nSamples
  } else if ( length(nSamples) == 2 ) {
    nSamplesSearch <- min(nSamples)
    nSamples <- max(nSamples)
  } else {
    stop("no more than 2 sample sizes may be specified in the 'nsamples' argument")
  }

  similarity <- exp( -temperature * as.matrix(distance) )

  single <- function(grid) {
    n <- length(grid)
    sils <-  numeric(n)
    for ( i in 1:n ) {
      samples <- if ( distr=="EPA" ) {
        .Call(.sample_epa, nSamples, similarity, grid[i], discount, nCores)
      } else {
        samplePartition(DDCRPPartition(similarity=similarity, mass=grid[i]), nSamples, randomizePermutation=TRUE)
      }
      est <- suppressWarnings(salso(samples, loss=loss, maxNClusters=maxNClusters))
      if ( length(unique(est)) == 1 ) {
        sils[i] <- NA
      } else {
        s <- summary(silhouette(est, 1-psm(samples)))
        sils[i] <- s$avg.width
      }
    }
    if ( all(is.na(sils)) ) {
      return(grid[n])
    } else {
      return(grid[which.max(sils)])
    }
  }

  rootfinder <- function(ncl) {
    if ( ncl == 1 ) return(-discount + .Machine$double.eps)
    n <- nrow(similarity)
    nsubsets.average <- function(mass, n) sum(mass / (mass + 1:n - 1))
    nsubsets.variance <- function(mass, n) sum((mass * (1:n - 1)) / (mass + 1:n - 1)^2)
    nclust <- function(mass) {
      .Call(.caviarpd_n_clusters, nSamplesSearch, similarity, mass, discount, loss=="VI", 16, maxNClusters, nCores)
    }

    sd.lwr <- if ( loss=='VI' ) { .5 * sd } else sd

    # Function to find the lower mass bound for a given cluster count
    lwr <- function(nClusters) { function(m) { nsubsets.average(m, n) +
        sd.lwr*sqrt(nsubsets.variance(m, n)) - nClusters }}

    # Function to find the upper mass bound for a given cluster count
    upr <- function(nClusters) { function(m) { nsubsets.average(m, n) -
        sd*sqrt(nsubsets.variance(m, n)) - nClusters }}

    # Get mass bounds and variance
    LB <- uniroot(lwr(ncl), c(.0001, 1000))$root
    UB <- uniroot(upr(ncl), c(.0001, 1000))$root
    bounds <- c(LB, UB)
    func <- function(mass) { nclust(mass) - ncl }

    # Find best mass
    tryCatch(uniroot(func, bounds)$root, error=function(e) NA)
  }

  # Mass selection
  if ( missing(mass) && !missing(nClusters) ) {

    if ( !is.numeric(nClusters) ) stop("'nClusters' argument must be numeric")
    nClusters <- as.integer(nClusters)

    if ( length(nClusters) == 1 ) {

      best <- rootfinder(nClusters)

    } else {

      mass.lwr <- rootfinder(min(nClusters))
      mass.upr <- rootfinder(max(nClusters))
      if ( is.na(mass.lwr) ) mass.lwr <- 0.25
      if ( is.na(mass.upr) ) mass.upr <- 5
      massGrid <- seq(mass.lwr, mass.upr, length=gridLength)
      best <- single(massGrid)

    }

  } else if ( missing(nClusters) && !missing(mass) ) {

    if ( !is.numeric(mass) ||  ( mass <= -discount ) ) stop("'mass' must be numeric and greater than -'discount'")
    if ( length(mass) == 1 ) {
      best <- mass
    } else {
      best <- single(mass)
    }

  } else if (missing(mass) && missing(nClusters)) {

    stop("must specify exactly one of 'mass' or 'nClusters' argument")

  } else {

    stop("cannot specify both 'mass' and 'nClusters' arguments")

  }

  mass <- best
  samples <- if (distr=="EPA") {
    .Call(.sample_epa, nSamples, similarity, mass, discount, nCores)
  } else if (distr=="ddCRP") {
    samplePartition(DDCRPPartition(similarity=similarity, mass=mass), nSamples, randomizePermutation=TRUE)
  } else {
    stop("partition distribution must be specified as either 'EPA' or 'ddCRP' in the 'distr' argument ")
  }
  if ( isTRUE(samplesOnly) ) return(samples)
  suppressWarnings(salso(samples, loss=loss, maxNClusters=maxNClusters))
}
