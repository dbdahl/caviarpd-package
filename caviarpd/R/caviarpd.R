#' Cluster Analysis via Random Partition Distributions
#'
#' Returns a clustering estimate given pairwise distances using the CaviarPD method.
#'
#' @param distance An object of class 'dist' or a pairwise distance matrix.
#' @param nClusters A numeric vector that specifies the range for the number of clusters to consider in the search for a clustering estimate. Should be missing if the \code{mass} argument is used. See `Details`.
#' @param mass A numeric vector of mass values to consider in the search for a clustering estimate. Should be missing if the \code{nClusters} argument is used. See `Details`.
#' @param nSamples The number of samples used to generate the clustering estimate.
#' @param gridLength The length of the grid search for an optimal mass parameter. Only applicable if a range of values are provided for \code{nClusters}.
#' @param samplesOnly If TRUE, the function only returns the samples generated for a given mass, temperature, and discount rather than an actual clustering estimate.
#' @param loss The SALSO method (Dahl, Johnson, Müller, 2021) tries to minimize this expected loss when searching the partition space for an optimal estimate. This must be either "binder" or "VI".
#' @param distr The random partition distribution used to generate samples.  This must be specified as either "EPA" or "ddCRP".
#' @param temperature A positive number that accentuates or dampens distance between observations.
#' @param similarity Either \code{"exponential"} or \code{"reciprocal"} to indicate the desired similarity function.
#' @param discount A number in \eqn{[0,1)} giving the discount parameter to control the distribution of subset sizes.
#' @param sd Number of standard deviations away from the expectation to be considered in finding boundary mass values.
#' @param maxNClusters The maximum number of clusters that can be considered by the SALSO method.
#' @param nCores The number of CPU cores to use. A value of zero indicates to use all cores on the system.
#'
#' @details
#' The \code{mass} argument is the main tuning parameter governing the number of clusters,
#'  with higher values tending toward more clusters. The \code{mass} is a real number bounded
#'  below by \eqn{-}\code{discount}. When a vector of mass values is supplied, a clustering
#'  estimate for each mass value is generated and the best clustering estimate is returned.
#'
#' Alternatively, a range for the number of clusters to be considered can be supplied with the
#'  \code{nClusters} argument. Mass values that return a clustering estimate with the minimum and
#'  maximum value of the range will estimated. A grid of mass values (of length \code{gridLength}) between
#'  the estimated min and max cluster mass values will be considered in the search for a clustering
#'  estimate. If \code{nClusters} is a single integer, then a clustering estimate with \code{nClusters}
#'  clusters will be returned.
#'
#' @return A object of class \code{salso.estimate}, which provides a clustering estimate (a vector of cluster labels) that can be displayed and plotted.
#'
#' @references
#'
#' D. B. Dahl, D. J. Johnson, and P. Müller (2021), Search Algorithms and Loss
#' Functions for Bayesian Clustering, [arXiv:2105.04451](https://arxiv.org/abs/2105.04451).
#'
#' @examples
#' # To reduce load on CRAN servers, limit the number of samples, grid length, and CPU cores.
#' set.seed(34)
#' iris.dis <- dist(iris[,-5])
#' est <- caviarpd(distance=iris.dis, mass=c(1, 2), nSamples=50, nCores=1)
#' samples <- caviarpd(distance=iris.dis, mass=1, nSamples=50, samplesOnly=TRUE, nCores=1)
#' est <- caviarpd(distance=iris.dis, nClusters=3, nSamples=50, nCores=1)
#' est <- caviarpd(distance=iris.dis, nClusters=3:5, nSamples=50, gridLength=5, nCores=1)
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
                     loss="binder", distr="EPA", temperature=10.0, similarity=c("exponential","reciprocal")[1], discount=0.0, sd=3, maxNClusters=0, nCores=0) {
  if ( is.matrix(distance) ) {
    if ( !isSymmetric(distance) || !is.numeric(distance) ) stop("'distance' is not a symmetric numerical matrix.")
  } else if ( class(distance) == 'dist' ) {
    distance <- as.matrix(distance)
  } else stop("'distance' argument must be an object of class 'dist' or a symmetric numerical matrix.")
  if ( !missing(nClusters) && (!is.numeric(nClusters) || !all(is.finite(nClusters)) || any(nClusters<1)) ) stop("'nClusters' must a numeric vector of finite values not less than 1")
  if ( !is.numeric(discount) || length(discount) != 1 || discount < 0 || discount >= 1.0 ) stop("'discount' must be in [0,1) and length 1")
  if ( !missing(mass) && (!is.numeric(mass) || !all(is.finite(mass)) || any( mass <= -discount )) ) stop("if supplied, 'mass' must be a numeric vector of finite values greater than -'discount'")
  if ( !is.numeric(nSamples) || ! length(nSamples) %in% c(1,2) || any(nSamples <= 0) || any(nSamples %% 1 != 0) ) stop("'nSamples' must be a strictly positive and length 1 or 2")
  if ( !is.numeric(gridLength) || length(gridLength) != 1 || gridLength < 2 || gridLength %% 1 != 0 ) stop("'gridLength' must be a strictly positive integer not less than 2")
  if ( !is.logical(samplesOnly) || !is.vector(samplesOnly) || length(samplesOnly) != 1 || ! samplesOnly %in% c(TRUE,FALSE) ) stop("'samplesOnly' must be a TRUE or FALSE")
  if ( length(loss) != 1 || ! loss %in% c("VI","binder") ) stop("'loss' must be either 'binder' or 'VI'")
  if ( distr != "EPA" && distr != "ddCRP" ) stop("'distr' must be either 'EPA' or 'ddCRP'")
  if ( !is.numeric(temperature) || !is.vector(temperature) || length(temperature) != 1 || temperature < 0 ) stop("'temperature' must be nonnegative and length 1")
  if ( !is.character(similarity) || length(similarity) != 1 || ! similarity %in% c("exponential","reciprocal") ) stop("'similarity' must be either 'exponential' or 'reciprocal'")
  if ( !is.numeric(sd) || length(sd) != 1 || sd < 1 || sd > 10 ) stop("'sd' must be nonnegative and between, say, 1 and 10")
  if ( !is.numeric(maxNClusters) || length(maxNClusters) != 1 || maxNClusters < 0 || maxNClusters %% 1 != 0 ) stop("'maxNClusters' must be 0 or a positive integer")
  if ( !is.numeric(nCores) || length(nCores) != 1 || nCores < 0 || nCores %% 1 != 0 ) stop("'nCores' must be 0 or a positive integer")
  if ( isTRUE(samplesOnly) && (missing(mass) || length(mass) > 1) ) stop( "must specify single mass parameter in order to obtain samples" )
  if ( length(nSamples) == 1 ) {
    nSamplesSearch <- max(1, 0.5 * nSamples)
  } else if ( length(nSamples) == 2 ) {
    nSamplesSearch <- min(nSamples)
    nSamples <- max(nSamples)
  } else {
    stop("no more than 2 sample sizes may be specified in 'nSamples'")
  }
  similarity <- if ( similarity == "exponential" ) {
    exp( -temperature * distance )
  } else if ( similarity == "reciprocal" ) {
    1/distance^temperature
  } else stop("Unsupported similarity")
  if ( ! all(is.finite(similarity)) ) stop("'distance', 'temperature', and/or 'similarity' yield similarity with nonfinite values")

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
        s <- summary(silhouette(est, 1-psm(samples, nCores)))
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
    if ( length(mass) == 1 ) {
      best <- mass
    } else {
      best <- single(mass)
    }
  } else if (missing(mass) && missing(nClusters)) {
    stop("must specify exactly one of 'mass' or 'nClusters'")
  } else {
    stop("cannot specify both 'mass' and 'nClusters'")
  }

  mass <- best
  samples <- if (distr=="EPA") {
    .Call(.sample_epa, nSamples, similarity, mass, discount, nCores)
  } else if (distr=="ddCRP") {
    samplePartition(DDCRPPartition(similarity=similarity, mass=mass), nSamples, randomizePermutation=TRUE)
  } else {
    stop("<impossible to get here because of previous check>")
  }
  if ( isTRUE(samplesOnly) ) return(samples)
  suppressWarnings(salso(samples, loss=loss, maxNClusters=maxNClusters, nCores=nCores))
}
