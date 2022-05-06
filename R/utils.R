#' scale_rows
#'
#' Normalize (z-score) rows of a matrix
#'
#' @param x a matrix, possibly containing gene by samples
#' @return Row-normalized matrix
#'
#'
#' @examples
#' \donttest{
#' library(cinaR)
#' data(atac_seq_consensus_bm) # calls 'bed'
#' bed.row.normalized <- scale_rows(bed[,c(4:25)])
#' head(bed.row.normalized)
#' }
#' @export
scale_rows <- function (x){
  m <- apply(x, 1, mean, na.rm = TRUE)
  s <- apply(x, 1, stats::sd, na.rm = TRUE)
  return((x - m) / s)
}


#' show_comparisons
#'
#' returns the names of the created comparisons
#'
#' @param results output of the cinaR
#' @return comparisons created
#' @export
show_comparisons <- function (results) {
  # if enrichment not run!
  if(!is.null(results[["cp"]])){
    return(names(results$DA.peaks))
  } else { # if run
    return(names(results$DA.results$DA.peaks))
  }
}

#' verboseFn
#'
#' returns a printing function to be used with in the script
#' @param verbose boolean, determines whether the output going be printed or not
#' @return print function
#'
#' @export
verboseFn <- function(verbose){
  if(verbose){
    return(function(...) {
      message(...)
    })
  } else {
    return(function(...){
      return(invisible(NULL))
    })
  }
}


