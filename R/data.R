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
#' @family nplug
#'
#' @export
nplug_mzroll <- function() {
  path <- system.file("extdata", "nplug.mzrollDB", package = "claman")
  if (path == "" || !file.exists(path)) {
    stop("nplug_mzroll was not found")
  }

  return(path)
}


#' NPLUG samples
#'
#' A table of meta-data for all NPLUG samples with variables:
#' \describe{
#'   \item{sample_name}{unique descriptor of each sample}
#'   \item{month}{month of sample generation}
#'   \item{replicate}{culture replicate}
#'   \item{DR}{dilution rate of culture - this is proportional to cells'
#'   growth rate}
#'   \item{limitation}{nutrient limiting growth
#'     \itemize{
#'       \item{GLU - carbon}
#'       \item{LEU - leucine (in a Leu4 auxotroph)}
#'       \item{NH4 - nitrogen}
#'       \item{PO4 - phosphorous}
#'       \item{URA - uracil (in a Ura2 auxotroph)}
#'     }
#'   }
#'   \item{exp_ref}{experimental or reference condition}
#'   \item{extraction}{filter- or pellet-based extraction}
#'   \item{condition}{Integer value for each unique condition. Here, that is
#'     a unique values of (month, limitation, DR, and extraction) for
#'     experimental samples, and unique values of (month and extraction) for
#'     reference samples.}
#'   \item{reference}{Condition number that this sample should be compared to.}
#'   }
#'
#' @family nplug
"nplug_samples"

#' NPLUG samples
#'
#' A table of meta-data for all NPLUG features with variables:
#' \describe{
#'   \item{compoundName}{compound name corresponding to the compoundName
#'     variable in the nplug_mzroll's features table.}
#'   \item{pathway}{A (rough) categorization of compounds into metabolic
#'     pathways.}
#'     }
#'
#' @family nplug
"nplug_compounds"

#' NPLUG MzRoll Augmented
#'
#' \link{nplug_mzroll} formatted with \link{process_mzroll} with sample
#'   (\link{nplug_samples}) and compound (\link{nplug_samples}) merged.
#'
#' @family nplug
"nplug_mzroll_augmented"

#' NPLUG MzRoll Normalized
#'
#' \link{nplug_mzroll_augmented} with injections collapsed using
#'   \link{collapse_injections}, followed by reference-sample normalization
#'   using \link{normalize_peaks}. Finally, reference samples and samples
#'   extracting using the pellet method were removed and sample names
#'   were cleaned up. These steps are described in the NPLUG vignette.
#'
#' @family nplug
"nplug_mzroll_normalized"
