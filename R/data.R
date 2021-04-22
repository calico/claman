#' Growth-Limiting Metabolites
#' 
#' The NPLUG experiment has a full factorial design exploring relationship
#'   between nutrient limitation and growth rate in yeast.
#' 
#' @details
#' The experiment contains 25 primary experimental conditions:
#' \itemize{
#'   \item{5 limiting nutrients: nitrogen, phosphorous, leucine, uracil and
#'     glucose (i.e., NPLUG!)}
#'   \item{5 dilution rates (growth rate / ln(2)): from 0.05-0.3 1/hrs} 
#' }
#' 
#' Metabolomics of 
#'   limited yeast cultures growing at different rates. 
#' 
#' @return path to the NPLUG .mzrollDB dataset
#' 
#' @export
nplug_mzroll <- function() {
  
  path <- system.file("extdata", "nplug.mzrollDB", package = "claman")
  if (path == "" || !file.exists(path)) {
    stop("nplug_mzroll was not found")
  }

  return(path)
}