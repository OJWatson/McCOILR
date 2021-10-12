#' Pipe operator
#'
#' See \code{\link[magrittr:pipe]{\%>\%}} for more details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL



#------------------------------------------------
#' Assign genotype calls based on a wsaf matrix
#'
#' This function triggers the c code for the categorical method
#' 
#' @param wsaf Within sample alle frequency matrix	An R data frame of SNP calling information. Row names are names of samples and column names are names of assays.
#' @param err  Error in sequencing for determining if is actually a 
#'   heterozygous call. Default = 0.1
#' 
#' @export
#' 

assignGTforREALMcCOIL <- function(wsaf, err = 0.1){
    
    GT <- matrix(NA, dim(wsaf)[1], dim(wsaf)[2])
    
    GT <- ifelse(wsaf > 1-err, 1,
                 ifelse(wsaf < 0+err, 0,
                        0.5))
    GT[is.na(GT)] <- -1
    return(GT)
    
}
