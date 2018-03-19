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
#' @references 
#' H.-J. Kümpel; Poroelasticity: parameters reviewed. Geophys J Int 1991; 105 (3): 783-799. doi: 10.1111/j.1365-246X.1991.tb00813.x
#' 
#' @details
#' \code{\link{poel_lambda}} gives the Lame constant, normalized by
#' \eqn{\mu}, the elastic shear modulus; this has no dependence on the undrained
#' Poisson's ratio
#' 
#' \code{\link{poel_alpha}} gives the effective stress coefficient, normalized by
#' \eqn{1 / B} where \eqn{B} is Skempton's coefficient. Hence, to obtain Biot's pore
#' pressure coefficient (generally written as \eqn{\alpha}) take this value and
#' divide by \eqn{B}.
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
#' 
#' # Make a sensitivity matrix
#' # Poisson's ratios
#' nu <- nuu <- seq(0,0.5,by=0.01)
#' alpgrd <- outer(nu, nuu, "poel_alpha")
#' B <- 0.75 # Skempton's coefficient
#' biot <- alpgrd/B
#' 
#' fields::image.plot(nu, nuu, biot, asp=1)
#' contour(nu, nuu, biot, add=TRUE)
#' 
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

#' @title Convert vector into step-like vector
#' @description Take a vector of values (with some set of indices)
#' and create a step-like function by shifting indices.
#' @details This will be useful for creating source-time-functions
#' from realistic data, that are formatted appropriately for poel.
#' @param x,y numeric vectors giving the coordinates of the points to be
#' converted into a step-like function
#' @export
#' @examples 
#' 
#' SV <- stepvec(seq_len(10), rnorm(10))
#' SV
#' 
#' \dontrun{
#' plot(y ~ x, SV, type='b')
#' }
stepvec <- function(x, y){
  x <- as.vector(x)
  y <- as.vector(y)
  x.. <- rep(x, each=2)
  y.. <- rep(y, each=2)
  nx <- length(x..)
  xstep <- x..[-1]
  ystep <- y..[-nx]
  data.frame(n=seq_along(xstep), x=xstep, y=ystep)
}
