#' Returns the string and number of the surface boundary condition
#' used in \code{poel}
#' @details
#' "whole","confined","unconfined"
# @export
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

#' @title Parameter sensitivity to Poisson's ratio
#' @description Functions to investigate sensitivity of normalized-parameters
#' @details
#' \code{\link{poel_lambda}} gives the Lame constant, normalized by
#' \eqn{\mu}, the elastic shear modulus; this has no dependence on the undrained
#' Poisson's ratio
#' 
#' \code{\link{poel_alpha}} gives the effective stress coefficient, normalized by
#' \eqn{1 / B} where B is Skempton's coefficient
#' 
#' \code{\link{poel_beta}} gives the bulk compressibility, normalized by
#' \eqn{1 / \mu / B^2}
#' 
#' \code{\link{poel_chi}} gives the Darcy conductivity, normalized by
#' \eqn{D / \mu / B^2} where D is the hydraulic diffusivity, which is proportional
#' to permeability
#' 
#' @param nu numeric; Poisson's ratio
#' @param nuu numeric; undrained Poisson's ratio
#' @name poel-params
#' @examples
#' a1 <- poel_alpha(0.25,0.33)
#' a2 <- poel_alpha(0.20,0.40)
#' a2/a1
#' 
#' # try multiple values
#' poel_alpha(0.25,seq(0,1,by=0.1))
NULL

#' @rdname poel-params
#' @export
poel_lambda <- function(nu){
  res <- 2*nu / (1 - 2*nu)
  return(res)
}
#' @rdname poel-params
#' @export
poel_alpha <- function(nu, nuu){
  res <- 3*(nuu - nu) / (1 - 2*nu) / (1 + nuu)
  res[nuu < nu] <- NA
  return(res)
}
#' @rdname poel-params
#' @export
poel_beta <- function(nu, nuu){
  res <- 9 * (1 - 2*nuu) * (nuu - nu) / 2 / (1 - 2*nu) / (1 + nuu)^2
  res[nuu < nu] <- NA
  return(res)
}
#' @rdname poel-params
#' @export
poel_chi <- function(nu, nuu){
  res <- 9 * (1 - nuu) * (nuu - nu) / 2 / (1 - nu) / (1 + nuu)^2
  res[nuu < nu] <- NA
  return(res)
}