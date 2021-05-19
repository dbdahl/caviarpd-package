#' Jaccard Distance for Categorical Attributes
#'
#' Computes the Jaccard pairwise distance matrix for categorical variables.
#'
#' @param data The data set of categorical variables for which the Jaccard distance is computed.
#'
#' @return An object of class "dist", which can be converted to an symmetric \eqn{n*n}
#' matrix. See the \code{\link[stats]{dist}} function documentation for details.
#'
#' @examples
#' jaccard(npk[,-5])
#' jaccard(warpbreaks[,-1])
#'
#' @export
#'
jaccard <- function(data) {
  levels <- lapply(1:ncol(data), function(i) sort(unique(data[,i])))
  expand <- t(apply(data, 1, function(x) Reduce(c,lapply(seq_along(x),
                                                         function(i) {
                                                           bin <- levels[[i]] %in% x[i]
                                                           if (length(bin) == 2) bin[1]
                                                           else bin
                                                         }))))
  dist(expand, method="binary")
}

