#' Returns the string and number of the surface boundary condition
#' used in \code{poel}
#' @details
#' "whole","confined","unconfined"
#' @export
#' @param bc numeric or character; the boundary condition; if this
#'             is missing, the full set of values is returned
#'
boundary.condition <- function(bc){
	#
	bcs <- data.frame(
		BC=c("whole","confined","unconfined"),
		BC.val=0:2)
	#
	ret <- if (missing(bc)){
		bcs
	} else {
		subset(bcs, BC==bc | BC.val==bc)
	}
	return(ret)
}
