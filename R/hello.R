# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


#' @title This is a title
#'
#' @description This is a description
#'
#' @param x The number
#'
#' @return The output from \code{\link{print}}
#' @export
#'
#' @examples
#' hello(5)
hello <- function(x) {
  x^2
}
