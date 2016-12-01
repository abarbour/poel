#' @title Read .inp files
#' @description Read-in and parse .inp files -- input files for poel -- for later use
#' @export
#' @param fi character; full-path to inp file
#' @param version the poel-version \code{fi} is formatted for; coerced to character
#' @seealso \code{\link{read_t}}
#' @examples
#'
#' # Example .inp #1
#' inpfi1 <- system.file('pe12_dk_s0_b1.inp', package='poel')
#' res1 <- read_inp(inpfi1)
#' print(res1)
#'
#' # Example .inp #2
#' inpfi2 <- system.file('pe12_dk_s0_b2.inp', package='poel')
#' res2 <- read_inp(inpfi2)
#' print(res2)
#'
read_inp <- function(fi, version=2012){

	# Setup some internal functions for parsing
	.splitter <- function(x., keep, TRANS=NULL, Stage=NULL){
	    if (!is.null(Stage)) message("parsing stage: ", Stage)
	    x. <- unlist(strsplit(trimws(x.),  " "))
		nc0 <- sapply(x., nchar) == 0
		x. <- x.[!nc0]
		if (missing(keep)) keep <- seq_along(x.)
		kept <- x.[keep]
		message(" + content: ", paste(kept, collapse="\t"))
		if (!missing(TRANS)){
			stopifnot(inherits(TRANS, 'function'))
			kept <- TRANS(kept)
		}
		return(kept)
	}
	N. <- as.numeric
	I. <- as.integer
	C. <- as.character
	CG. <- function(x) gsub("\"","", gsub("'","", as.character(x)))

	version <- match.arg(as.character(version), c('2012','2006'))
	if (version == '2006') stop('Parsing of version-',version," input files is currently unsupported.")

	# Get raw INP file
	.Inp. <- readr::read_lines(fi)

	# Identify commented lines
	comm.1 <- sapply(.Inp., function(x) grepl("^#", x), USE.NAMES=FALSE)

	# Get interlaced comments
	comments.interlaced <- .Inp.[comm.1] # comments -- not necessarily with footer comments

	# Get content to process -- may have footer comments (if applicable)
	Cont <- .Inp.[!comm.1]

	#
	# Parse raw content, section-by-section...
	#

	#       SOURCE PARAMETERS A: SOURCE GEOMETRY
	#       ====================================
	# 1. source top and bottom depth [m]
	# 2. source radius (> 0) [m]

	ind <- 1
	src <- .splitter(Cont[ind], 1:2, N., Stage='S-A')

	ind <- 2
	src_rad <- .splitter(Cont[ind], 1, N., Stage='S-A-r')

	#       SOURCE PARAMETERS B: SOURCE TYPE
	#       ================================
	# 1. selection of source type:
	#    0 = initial excess pore pressure within the source volume (initial value problem): IVP
	#    1 = injection within the source volume (boundary value problem): BVP

	ind <- 3
	src_type <- .splitter(Cont[ind], 1, I., Stage='S-B')
	is.ivp <- src_type==0
	src_type_lbl <- ifelse(is.ivp, "IVP-PP", "BVP-INJ")

	#       SOURCE PARAMETERS C: SOURCE TIME HISTORY
	#       ========================================
	# IF(initial pore pressure)THEN
	#   1. value of the initial pressure energy (pressure*volume) [Pa*m^3]
	#      Note: in this case, finite source volume is required.
	# ELSE IF(injection)THEN
	#   1. number of data lines describing the source time history
	#   2. listing of the injection rate time series

	ind <- 4
	first.data <- ind + 1
	stf <- if (is.ivp){
		last.data <- first.data
		data.frame(Init.ppvol = .splitter(Cont[ind], 1, I., Stage='S-C-ivp'))
	} else {
		ndata <- .splitter(Cont[ind], 1, I., Stage='S-C-bvp')
		last.data <- (first.data + ndata - 1)
		tmpstf <- plyr::ldply(Cont[first.data:last.data], .splitter, keep=1:3, TRANS=N.)
		names(tmpstf) <- c("Num", "Time.sec", "Injection.rate.mcps")
		tmpstf
	}

	#       RECEIVER PARAMETERS A: RECEIVER DEPTH SAMPLING
	#       ==============================================
	# 1. switch for equidistant steping (1/0 = yes/no)
	# 2. number of receiver depth samples (<= nzrmax defined in "peglobal.h")
	# 3. if equidistant, start depth [m], end depth [m]; else list of depths
	#    (all >= 0 and ordered from small to large!)

	ind <- last.data + 1
	rcv_type <- .splitter(Cont[ind], 1, I., Stage='R-A')
	is.eqd <- rcv_type==1

	ind <- ind + 1
	ndeps <- .splitter(Cont[ind], 1, I., Stage='R-A-ndep')

	ind <- ind + 1
	.keep.dep <- if (is.eqd){
		1:2
	} else {
		seq_len(ndeps)
	}
	rcv_depths <- .splitter(Cont[ind], .keep.dep, N., Stage='R-A-deps')

	#       RECEIVER PARAMETERS B: RECEIVER DISTANCE SAMPLING
	#       =================================================
	# 1. switch for equidistant steping (1/0 = yes/no)
	# 2. number of receiver distance samples (<= nrmax defined in "peglobal.h")
	# 3. if equidistant, start distance [m], end distance [m]; else list of
	#    distances (all >= 0 and ordered from small to large!)

	ind <- ind + 1
	rcv_type_dist <- .splitter(Cont[ind], 1, I., Stage='R-B')
	is.eqd_dist <- rcv_type_dist==1

	ind <- ind + 1
	ndists <- .splitter(Cont[ind], 1, I., Stage='R-B-ndist')

	ind <- ind + 1
	.keep.dist <- if (is.eqd_dist){
		1:2
	} else {
		seq_len(ndists)
	}
	rcv_dists <- .splitter(Cont[ind], .keep.dist, N., Stage='R-B-dist')

	#       RECEIVER PARAMETERS C: Time SAMPLING
	#       ====================================
	# 1. time window [s]
	# 2. number of time samples
	#    Note: the caracteristic diffusion time =
	#          max_receiver_distance^2 / diffusivity_of_source_layer

	ind <- ind + 1
	twin <- .splitter(Cont[ind], 1, N., Stage='T-A')

	ind <- ind + 1
	ntimes <- .splitter(Cont[ind], 1, I., Stage='T-A-nsamp')

	#       WAVENUMBER INTEGRATION PARAMETERS
	#       =================================
	# 1. relative accuracy (0.01 for 1% error) for numerical wavenumber integration;

	ind <- ind + 1
	int_accuracy <- .splitter(Cont[ind], 1, N., Stage='.-IntAcc')

	#       OUTPUTS A: DISPLACEMENT TIME SERIES
	#       ===================================
	# 1. select the 2 displacement time series (1/0 = yes/no)
	#    Note Ut = 0
	# 2. file names of these 2 time series

	ind <- ind + 1
	fiout_disps_write <- .splitter(Cont[ind], 1:2, I., Stage='OUT-A-d')

	ind <- ind + 1
	fiout_disps <- .splitter(Cont[ind], 1:2, CG., Stage='OUT-A-dnames')

	#       OUTPUTS B: STRAIN TENSOR & TILT TIME SERIES
	#       ===========================================
	# 1. select strain time series (1/0 = yes/no): Ezz, Err, Ett, Ezr (4 tensor
	#    components) and Tlt (= -dur/dz, the radial component of the vertical tilt).
	#    Note Ezt, Ert and Tlt (tangential tilt) = 0
	# 2. file names of these 5 time series

	ind <- ind + 1
	fiout_strns_write <- .splitter(Cont[ind], 1:5, I., Stage='OUT-B-s')

	ind <- ind + 1
	fiout_strns <- .splitter(Cont[ind], 1:5, CG., Stage='OUT-B-snames')

	#       OUTPUTS C: PORE PRESSURE & DARCY VELOCITY TIME SERIES
	#       =====================================================
	# 1. select pore pressure and Darcy velocity time series (1/0 = yes/no):
	#    Pp (excess pore pressure), Dvz, Dvr (2 Darcy velocity components)
	#    Note Dvt = 0
	# 2. file names of these 3 time series

	ind <- ind + 1
	fiout_pp_write <- .splitter(Cont[ind], 1:3, I., Stage='OUT-C-p')

	ind <- ind + 1
	fiout_pp <- .splitter(Cont[ind], 1:3, CG., Stage='OUT-C-pnames')

	#       OUTPUTS D: SNAPSHOTS OF ALL OBSERVABLES
	#       =======================================
	# 1. number of snapshots
	# 2. time[s] (within the time window, see above) and output filename of
	#    the 1. snapshot
	# 3. ...

	ind <- ind + 1
	nsnap <- .splitter(Cont[ind], 1, I., Stage='OUT-D-snap')

	first.snap <- ind + 1
	last.snap <- first.snap + nsnap - 1
	snaps <- plyr::ldply(Cont[first.snap:last.snap], .splitter, keep=1:2)
	names(snaps) <- c('Time.s','Snap')
	snaps$Snap <- CG.(snaps$Snap)

	#       GLOBAL MODEL PARAMETERS
	#       =======================
	# 1. switch for surface conditions:
	#    0 = without free surface (whole space),
	#    1 = unconfined free surface (p = 0),
	#    2 = confined free surface (dp/dz = 0).
	# 2. number of data lines of the layered model (<= lmax as defined in
	#    "peglobal.h") (see Note below)

	ind <- last.snap + 1
	surf_bc <- .splitter(Cont[ind], 1, I., Stage='MDL-surf-bc')

	surf_bc_lbl <- switch(paste0('BC',surf_bc), `BC0`='WHOLE-S', `BC1`='UNCONF-HS', `BC2`='CONF-HS')
	if (is.null(surf_bc_lbl)) warning('Bad surface boundary condition')

	ind <- ind + 1
	nlay <- .splitter(Cont[ind], 1, I., Stage='MDL-nlay')

	#       MULTILAYERED MODEL PARAMETERS
	#       =============================
	#
	#   Note: mu = shear modulus
	#         nu = Poisson ratio under drained condition
	#         nu_u = Poisson ratio under undrained condition (nu_u > nu)
	#         B = Skempton parameter (the change in pore pressure per unit change
	#             in confining pressure under undrained condition)
	#         D = hydraulic diffusivity
	#
	# no depth[m]      mu[Pa]    nu    nu_u       B     D[m^2/s]

	first.lay <- ind + 1
	last.lay <- first.lay + nlay - 1
	lays <- plyr::ldply(Cont[first.lay:last.lay], .splitter, keep=1:7, TRANS=N., Stage='MDL-lay')
	Lp <- ._lith_()
	names(lays) <- Lp$Params
	lays <- as.matrix(lays)
	class(lays) <- 'lith'
	attr(lays, "Units") <- Lp$Units
	attr(lays, 'Nlith') <- nlay

	stopifnot(poel::is.lith(lays))

	first.footer <- last.lay + 1
	last.footer <- max(first.footer, length(Cont))
	foot.inds <- first.footer:last.footer
	footer <- paste("#", Cont[foot.inds])

	Comments <- c(comments.interlaced, footer)

	Inp <- list(
		File = fi,
		Version = version,
		Raw = .Inp.,
		Comments = Comments,
		Content = Cont[-foot.inds],
		Parsed.Content = list(
			BC = c(Surface=surf_bc_lbl, Source=src_type_lbl),
			#Surface = c(BC=surf_bc, BC_type=surf_bc_lbl),
			Layers = lays,
			Source = list(Depth=src, Rad=src_rad, STF=stf),
			Receiver = list(
				DiscretizationType = c(Depth=rcv_type, Dist=rcv_type_dist, Time=1),
				Is.eqd = c(Depth=is.eqd, Dist=is.eqd_dist, Time=TRUE),
				N = c(Depth=ndeps, Dist=ndists, Time=ntimes),
				Loc = list(Depth=rcv_depths, Dist=rcv_dists, Time=c(0, twin))
			)
		)
	)
	class(Inp) <- c('poel-inp', version, src_type_lbl)
	return(Inp)
}

#' @title Methods for 'poel-inp' class
#' @description Methods for 'poel-inp' class
#' @rdname poel-inp-methods
#' @inheritParams read_inp
#' @export
#' @seealso \code{\link{read_inp}} \code{\link{read_t}}
`print.poel-inp` <- function(x){
    PC <- x[['Parsed.Content']]
    message("      POEL INP file: ", x[['File']])
    message(" Initial conditions: ", PC[['BC']][['Source']])
    message("Boundary conditions: ", PC[['BC']][['Surface']])
}

#' @rdname poel-inp-methods
#' @export
get_lith <- function(x, ...) UseMethod('get_lith')
#' @rdname poel-inp-methods
#' @export
`get_lith.poel-inp` <- function(x, ...){
    x[['Parsed.Content']][['Layers']]
}
