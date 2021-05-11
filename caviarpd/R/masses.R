#' Mass Parameter Selection for the CaviarPD Procedure
#'
#' User inputs a range of clustering sizes to obtain a mass value that can viably correspond to each cluster count.
#'
#' @param distance Pairwise distance matrix of class 'dist'.
#' @param ncl.range A vector of two values representing the minimum and maximum number of clusters that the user wishes to consider in selecting mass values.
#' @param single If TRUE, the algorithm returns both a list of masses for each targeted cluster count as well as a best overal mass selected from that list.
#' @param nSD Number of standard deviations in the EPA distribution to consider for a range of mass values; lower values allow the algorithm to run more quickly but risk an error if deemed insufficient.
#' @param temperature A positive number that accentuates or dampens distance between observations.
#' @param discount Typically 0, controls the distribution of subset sizes.
#' @param loss The salso method aims to estimate this loss function when searching the partition space for an optimal estimate, must be specified as either "binder" or "VI".
#' @param nSamplesA Number of samples used by the salso method to repeatedly calculate the cluster count for each attempted mass value.
#' @param nSamplesB Number of samples used by the salso method to obtain the optimal mass; only applicable if single=TRUE.
#' @param w Weights for selecting a single mass. The first weight is attached to the partition confidence, the second weight is attached to the variance trio, and the last is attached to the nunber of clusters.
#' @param nCores The number of CPU cores to use. A value of zero indicates to use all cores on the system.
#'
#' @return If single==FALSE, returns a list with two elements: a chain of mass values corresponding to each cluster count,
#' and the mass value that had the best overall confidence plot. If single==TRUE, returns only the best overall mass value.
#'
#' @examples
#' tooth.dis <- dist(scale(ToothGrowth[,-2]))
#' # In practice, use at least 100 samples and multiple cores. Less here for fast-running examples.
#' select.masses(tooth.dis, ncl.range=c(2,4), nSamplesA=10, nSamplesB=10, nCores=1)
#' iris.dis <- dist(iris[,-5])
#' select.masses(iris.dis, ncl.range=c(3,6), single=TRUE, nSamplesA=10, nSamplesB=10, nCores=1)
#'
#' @export
#' @importFrom stats dist uniroot var median
#'
select.masses <- function(distance, ncl.range, single=FALSE, nSD=3, discount=0.0, temperature=10.0,
                          loss='binder', nSamplesA=500, nSamplesB=1000, w=c(1,1,0), nCores=0) {

  ### ERROR CHECKING ###
  if (class(distance) != 'dist') stop(" 'distance' argument must be an object of class 'dist' ")
  else if (length(ncl.range) != 2 | !is.numeric(ncl.range)) stop(" 'ncl.range' argument must be numeric vector of 2 elements ")
  else if (!is.numeric(nSD) | nSD %% 1 !=0 | nSD < 1) stop(" 'nSD' argument must be a positive integer >= 1")
  else if (length(w) != 3) stop(" weights must have 3 elements ")
  else if (!is.numeric(w) | !(all(w >= 0))) stop(" weights must be numeric and nonnegative ")
  if (ncl.range[1] > ncl.range[2]) ncl.range <- ncl.range[c(2,1)]

  nsubsets.average <- function(mass, n) sum(mass / (mass + 1:n - 1))
  nsubsets.variance <- function(mass, n) sum((mass * (1:n - 1)) / (mass + 1:n - 1)^2)
  nclust <- function(mass, distance, loss='binder') { length(unique(caviarPD(distance=distance, mass, temperature=temperature, discount=discount, loss=loss, nSamples=nSamplesA, nCores=nCores))) }

  bounds.by.ncl <- function(ncl.range, nSD=3) {

    if (loss=='VI') nSD.lwr <- .5 * nSD else nSD.lwr <- nSD
    # Function to find the lower mass bound for a given cluster count
    lwr <- function(nClusters) { function(m) { nsubsets.average(m, 100) +
        nSD.lwr*sqrt(nsubsets.variance(m, 100)) - nClusters }}

    # Function to find the upper mass bound for a given cluster count
    upr <- function(nClusters) { function(m) { nsubsets.average(m, 100) -
        nSD*sqrt(nsubsets.variance(m, 100)) - nClusters }}

    # Initialize vectors for loop
    ncl <- seq(ncl.range[1], ncl.range[2], by=1)
    UB <- numeric(length(ncl))
    LB <- numeric(length(ncl))
    mid <- numeric(length(ncl))
    v <- numeric(length(ncl))

    # Get mass bounds and variance for each cluster count
    for (i in 1:length(ncl)) {
      LB[i] <- uniroot(lwr(ncl[i]), c(.001, 100))$root
      UB[i] <- uniroot(upr(ncl[i]), c(.001, 100))$root
      mid[i]  <- mean(c(LB[i], UB[i]))
    }

    # Return data
    data.frame(NClusters = ncl, Lower = round(LB,2), Upper = round(UB,2), Estimate=round(mid,2))
  }

  # Get an upper and lower bound of mass values to test for each number of clusters
  df <- bounds.by.ncl(ncl.range, nSD)
  ncls <- df$NClusters
  masses <- numeric(length(ncls))

  # Changes below
  ################################################

  n <- nrow(df)
  func2 <- function(ncl, distance, loss) { function(mass) { nclust(mass, distance, loss) - ncl } }
  boundsA <- c(df$Lower[1], df$Upper[1])
  boundsB <- c(df$Lower[n], df$Upper[n])
  masses[1] <- tryCatch( uniroot(func2(ncls[1], distance, loss), boundsA)$root, error=function(e) NA)
  masses[n] <- tryCatch( uniroot(func2(ncls[n], distance, loss), boundsB)$root, error=function(e) NA)
  for (i in 2:(n/2)) {
    bounds <- c(max(df$Lower[i], masses[i-1], na.rm=TRUE), min(df$Upper[i], masses[n-i+2], na.rm=TRUE))
    masses[i] <- tryCatch( uniroot(func2(ncls[i], distance, loss), bounds)$root, error=function(e) NA)
    masses[n-i+1] <- tryCatch( uniroot(func2(ncls[n-i+1], distance, loss), bounds)$root, error=function(e) NA)
  }
  if (n %% 2 == 1) {
    k <- median(1:length(ncls))
    masses[k] <- tryCatch( uniroot(func2(ncls[k], distance, loss), boundsA)$root, error=function(e) NA)
  }


  ################################################


  # Run the single,mass function if an optimal mass value is needed, otherwise, just return the sequence of masses
  mat <- matrix(sort(masses), nrow=1)
  colnames(mat) <- ncls[!is.na(masses)]
  if (single==FALSE) {
    return(mat)
  } else {
    final_mass <- single.mass(masses[!is.na(masses)], distance, temperature=temperature, discount=discount, nSamples=nSamplesB, w=w, loss=loss)
    return(list(masses=mat, best=final_mass))
  }
}


#' Single Mass Parameter Selection for the CaviarPD Procedure
#'
#' Calculates the partition confidence and variance ratios for each mass value to find the best mass. User can input masses from the select.masses function or supply their own.
#'
#' @param masses A vector of mass values from which to select the best mass. This can be a simple sequence for some range or a list of masses generated from the select.masses function.
#' @param distance Pairwise distance matrix of class 'dist'.
#' @param temperature A positive number that accentuates or dampens distance between observations.
#' @param discount Typically 0, controls the distribution of subset sizes.
#' @param nSamples The number of samples used to estimate the loss function in the salso method.
#' @param w Weights for selecting a single mass. The first weight is attached to the partition confidence, the second weight is attached to the variance trio, and the last is attached to the nunber of clusters.
#' @param loss The salso method aims to estimate this loss function when searching the partition space for an optimal estimate, must be specified as either "binder" or "VI".
#' @param nCores The number of CPU cores to use. A value of zero indicates to use all cores on the system.
#'
#' @return The value of the mass parameter that should be used in the caviarPD method for the given pairwise distance matrix.
#'
#' @examples
#' iris.dis <- dist(iris[,-5])
#' # In practice, use at least 100 samples and multiple cores. Less here for fast-running examples.
#' iris.masses <- select.masses(iris.dis, ncl.range=c(3,6), nSamplesA=10, nSamplesB=10, nCores=1)
#' single.mass(masses=iris.masses, distance=iris.dis, nSamples=10, nCores=1)
#' single.mass(masses=seq(.5, 2, by=.25), distance=iris.dis, nSamples=10, nCores=1)
#'
#' @importFrom salso salso psm
#'
#' @export
#'
single.mass <- function(masses, distance, temperature=10.0, discount=0.0, nSamples=1000, w=c(1,1,0), loss='binder', nCores=0) {

  ### ERROR CHECKING ###
  if (class(distance) != 'dist') stop(" 'distance' argument must be an object of class 'dist' ")
  else if( !is.numeric(masses) ) stop(" 'masses' must all be numeric and greater than -'discount' ")
  else if (length(w) != 3) stop(" weights must have 3 elements ")
  else if (!is.numeric(w) | !(all(w >= 0))) stop(" weights must be numeric and nonnegative ")

  pc <- numeric(length(masses))
  wc_var <- numeric(length(masses))
  total_var <- numeric(length(masses))
  ncl <- numeric(length(masses))
  for (i in 1:length(masses)) {
    b <- caviarPD(distance, masses[i], temperature=temperature, discount=discount, nSamples=nSamples, samplesOnly=TRUE, nCores=nCores)
    x <- salso(b, loss=loss)
    psmat <- psm(b)
    ncl[i] <- length(unique(x))
    pc_num <- sum(sapply(unique(x), function(label) {
      w <- which(x==label)
      m <- psmat[w,w]
      if (!is.null(nrow(m))) { diag(m) <- 0; sum(m)/2 } else { 0 } # Sometimes m just comes out to be a scalar so it has no "diagonals"
    }))
    pc_den <- sum(sapply(unique(x), function(label) choose(sum(x==label),2)))
    pc[i] <- pc_num/pc_den

    wcv_num <- sum( sapply(unique(x),  # Set NA values to zero
                           function(label) {
                             w <- which(x==label)
                             m <- as.matrix(psmat[w,w])
                             mvar <- if(nrow(m) <= 2) 0 else var(m[lower.tri(m)])
                             choose(sum(x==label),2) * mvar
                           }))
    wcv_den <- sum(sapply(unique(x), function(label) choose(sum(x==label),2)))
    wc_var[i] <- wcv_num / wcv_den
    total_var[i] <- var(psmat[lower.tri(psmat)])
  }
  vr <- wc_var / total_var

  final_weighted_avg <- numeric(length(masses))
  rank_pc <- rank(-pc)
  rank_vr <- rank(vr)
  rank_ncl <- rank(ncl)
  for (i in 1:length(masses)) { final_weighted_avg[i] <- w[1]*rank(-pc)[i] + w[2]*rank(vr)[i] + w[3]*rank(ncl)[i] }
  return(masses[which.min(final_weighted_avg)])
}


