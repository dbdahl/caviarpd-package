#' Mass Parameter Selection for the CaviarPD Procedure
#'
#' User inputs a range of clustering sizes to obtain a mass value that can viably correspond to each cluster count.
#'
#' @inheritParams caviarPD
#' @param ncl.range A vector of two values representing the minimum and maximum number of clusters that the user wishes to consider in selecting mass values.
#' @param single If TRUE, the algorithm returns both a list of masses for each targeted cluster count as well as a best overall mass selected from that list.
#' @param nSD Number of standard deviations in the EPA distribution to consider for a range of mass values; lower values allow the algorithm to run more quickly, but also risk not finding the mass for some cluster counts.
#' @param nSamplesSearch Number of samples used for each iteration of the search algorithm to find the mass value corresponding to a desired number of clusters. Setting this argument too high can result in long computations.
#' @param w Weights for selecting a single mass. The first weight is attached to the partition confidence, the second weight is attached to the variance trio, and the last is attached to the number of clusters.
#'
#' @return If single==FALSE, returns a list with two elements: a chain of mass values corresponding to each cluster count
#' and the mass value that had the best overall confidence plot. If single==TRUE, returns only the best overall mass value.
#'
#' @references
#'
#' D. B. Dahl, D. J. Johnson, and P. Müller (2021), Search Algorithms and Loss
#' Functions for Bayesian Clustering, <arXiv:2105.04451>.
#'
#' @example man/examples/select.masses.R
#' @export
#' @importFrom stats dist uniroot var median
#'
select.masses <- function(distance, ncl.range, single=FALSE, nSD=3, discount=0.0, temperature=10.0,
                          loss='binder', maxNClusters=0, nSamplesSearch=500, nSamples=1000, w=c(1,1,0), nCores=0) {

  ### ERROR CHECKING ###
  if (class(distance) != 'dist') stop(" 'distance' argument must be an object of class 'dist' ")
  else if (length(ncl.range) != 2 | !is.numeric(ncl.range)) stop(" 'ncl.range' argument must be numeric vector of 2 elements ")
  else if (!is.numeric(nSD) | nSD %% 1 !=0 | nSD < 1) stop(" 'nSD' argument must be a positive integer >= 1")
  else if (length(w) != 3) stop(" weights must have 3 elements ")
  else if (!is.numeric(w) | !(all(w >= 0))) stop(" weights must be numeric and nonnegative ")
  if (ncl.range[1] > ncl.range[2]) ncl.range <- ncl.range[c(2,1)]

  nsubsets.average <- function(mass, n) sum(mass / (mass + 1:n - 1))
  nsubsets.variance <- function(mass, n) sum((mass * (1:n - 1)) / (mass + 1:n - 1)^2)
  similarity <- exp( -temperature * as.matrix(distance) )
  nclust <- function(mass) {
    caviarpd_n_clusters(nSamplesSearch, similarity, mass, discount, loss=="VI", 16, maxNClusters, nCores, mkSeed())
  }

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
  func <- function(ncl) { function(mass) { nclust(mass) - ncl } }
  boundsA <- c(df$Lower[1], df$Upper[1])
  boundsB <- c(df$Lower[n], df$Upper[n])
  masses[1] <- tryCatch( uniroot(func(ncls[1]), boundsA)$root, error=function(e) NA)
  masses[n] <- tryCatch( uniroot(func(ncls[n]), boundsB)$root, error=function(e) NA)
  for (i in 2:(n/2)) {
    bounds <- c(max(df$Lower[i], masses[i-1], na.rm=TRUE), min(df$Upper[i], masses[n-i+2], na.rm=TRUE))
    masses[i] <- tryCatch( uniroot(func(ncls[i]), bounds)$root, error=function(e) NA)
    masses[n-i+1] <- tryCatch( uniroot(func(ncls[n-i+1]), bounds)$root, error=function(e) NA)
  }
  if (n %% 2 == 1) {
    k <- median(1:length(ncls))
    masses[k] <- tryCatch( uniroot(func(ncls[k]), boundsA)$root, error=function(e) NA)
  }

  ################################################

  # Run the single,mass function if an optimal mass value is needed, otherwise, just return the sequence of masses
  mat <- matrix(sort(masses), nrow=1)
  colnames(mat) <- ncls[!is.na(masses)]
  if (single==FALSE) {
    return(mat)
  } else {
    final_mass <- single.mass(masses[!is.na(masses)], distance, temperature=temperature, discount=discount, nSamples=nSamples, w=w, loss=loss)
    return(list(masses=mat, best=final_mass))
  }
}


#' Single Mass Parameter Selection for the CaviarPD Procedure
#'
#' Calculates the partition confidence and variance ratios for each mass value to find the best mass. Users can input masses from the \code{\link{select.masses}()} function or supply their own.
#'
#' @inheritParams select.masses
#' @param masses A vector of mass values from which to select the best mass. This can be a simple sequence for some range or a list of masses generated from the select.masses function.
#'
#' @return The value of the mass parameter that should be used in the CaviarPD method for the given pairwise distance matrix.
#'
#' @example man/examples/single.mass.R
#' @importFrom salso salso psm
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
  similarity <- exp( -temperature * as.matrix(distance) )
  for (i in 1:length(masses)) {
    b <- sample_epa(nSamples, similarity, masses[i], discount, nCores, mkSeed())
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


