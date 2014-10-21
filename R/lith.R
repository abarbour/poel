#' Lith
#' @export
#' @examples
#' # just do something random to demonstrate
#' n <- 5 # five layers
#' deps <- seq_len(n)-1 # at these depths
#' rand <- runif(n)
#' (l <- lithology(deps, rand, rand, rand))
lithology <- function(Depth, 
                      ShearModulus=1.e9, Skempton=0.85, Diffusivity=0.1, 
                      Nu=0.25, Nuu=0.35, Number=NULL){
  if (is.null(Number)) Number <- seq_along(Depth)
  Depth <- as.numeric(Depth)
  stopifnot(all(Depth>=0))
  ord <- order(Depth)
  L <- cbind(Number, Depth, ShearModulus, Nu, Nuu, Skempton, Diffusivity)[ord, , drop=FALSE]
  class(L) <- 'lith'
  attr(L,"Units") <- c("","m","Pa","","","","m^2/s")
  return(L)
}

#' @rdname lithology
#' @export
plot.lith <- function(L, no.layout=TRUE, ...){
  Units <- attr(L,"Units")
  Deps <- L[,"Depth"]
  ldim <- dim(L)
  nms <- colnames(L)
  #
  nlay <- ldim[1]
  if (nlay==1){
    stop("Cannot plot this lithology -- only a single layer")
  }
  #
  layers <- seq_len(ldim[1])
  Deps.t <- Deps[layers[-1]]
  drng <- range(c(0,pretty(1.05*Deps[Deps>0])))
  #
  PLT <- function(n, ...){
    D <- L[,n]
    Du <- Units[n]
    Dnm <- nms[n]
    sf <- stepfun(Deps.t, D)
    plot(sf, lty=2, pch=NA, xaxs="i", xlab="Depth", main="", ylab=parse(text=Du), ...)
    lines(sf, vertical=FALSE, lwd=3, pch=NA, xaxs="i", ...)
    lims <- par("usr")
    text(lims[2]*0.98, lims[4]*0.95, Dnm, pos=2, font=2)
    return(sf)
  }
  toPlot <- c(3,4,5,6,7)
  if (!no.layout){
    layout(matrix(seq_along(toPlot),ncol=1))
  }
  op <- par(mar=c(3.8,4.1,1,1), mgp=c(2.5,1,0), las=1, lend="butt")
  on.exit(par(op))
  #
  invisible(lapply(toPlot, PLT, xlim=drng, ...))
}
