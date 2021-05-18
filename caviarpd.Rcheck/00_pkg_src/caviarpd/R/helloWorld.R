#' My First Function
#'
#' This function greets people.
#'
#' @param name A character vector giving the names of the people to greet.
#'
#' @returns A character vector of the greeting.
#'
#' @export
#'
helloWorld <- function(name) {
  paste0("Hello, ",name)
}
