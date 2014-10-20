#' Construct an input file suitable for \code{poel}
#' @export
#' @param lith, the lithology xxx
#' @param ... additional parameters
#' @family Input-file
#'
inp <- function(lith, ...){ UseMethod("inp") }
#' @rdname inp
#' @export
inp.lith <- function(lith, ...){
  lith <- unclass(lith)
  inp(lith, ...)
}
#' @rdname inp
#' @export
inp.default <- function(lith,
                obs.depth=0,
                surface.bc=c("unconfined","confined","whole"),
                file="myinp",
                stf.params=list(file="stf/16dayave", src="PR-M", len=16*24),
                sim.params=list(time.win=1440, time.pts=time.win, accuracy=0.05),
                well.params=list(name="PR-M", st=0, en=100, rad=0.12),
                bsm.params=list(name=c("B089","B082"), dist=c(285,442), dep=c(NA,NA)),
                results.dir=".", ...){
  #
  lith <- as.matrix(lith)
  #
  #   # This is the input file of FORTRAN77 program "poel06" for modeling
  #   # coupled deformation-diffusion processes based on a multi-layered (half-
  #   # or full-space) poroelastic media induced by an injection (pump) of
  #   # from a borehole or by a (point) reservoir loading.
  #   #
  #   # by R. Wang,
  #   # GeoForschungsZentrum Potsdam
  #   # e-mail: wang@gfz-potsdam.de
  #   # phone 0049 331 2881209
  #   # fax 0049 331 2881204
  #   #
  #   # Last modified: Potsdam, Nov, 2006
  #   #
  #   ##############################################################
  #   ##                                                          ##
  #   ## Cylindrical coordinates (Z positive downwards!) are used ##
  #   ## If not others specified, SI Unit System is used overall! ##
  #   ##                                                          ##
  #   ## Tilt is positive when the upper end of a borehole tilt-  ##
  #   ## meter body moves away from the pumping well.             ##
  #   ##                                                          ##
  #   ##############################################################
  #   #
  #   #        SOURCE PARAMETERS
  #   #        =================
  #   # 1. start and end depth [m] of the injection (pump) screen, and well radius [m];
  #   # 2. number of input data lines for the source time function;
  #   # 3. listing of the source time series.
  #   #-------------------------------------------------------------------------------
  #   %(s_start_depth)g %(s_end_depth)g  %(s_radius)g                 |dble: s_start_depth, s_end_depth, s_radius;
  #   #-------------------------------------------------------------------------------
  #   2
  #   #-------------------------------------------------------------------------------
  #   # no    time    source_function
  #   # [-]   [hour]  [m^3/h for injection(+)/pump(-) rate]
  #   #-------------------------------------------------------------------------------
  #   %(source_function)s
  #   #-------------------------------------------------------------------------------
  #   #
  #   #        RECEIVER PARAMETERS
  #   #        ===================
  #   # 1. observation depth [m]
  #   # 2. switch for equidistant trace steping (1/0 = yes/no)
  #   # 3. number of distance samples [m] (<= nrmax defined in "peglobal.h")
  #   # 4. if equidistant, start trace distance, end trace distance; else list of
  #   #    trace distances (all >= 0 and ordered from small to large!)
  #   # 5. length of time window [h], number of time samples
  #   #-------------------------------------------------------------------------------
  #   %(receiver_depth)g              |dble: r_depth;
  #   %(sw_equidistant)i              |int: sw_equidistant;
  #   %(no_distances)i                |int: no_distances;
  #   %(str_distances)s               |dble: d_1,d_n; or d_1,d_2, ...;
  #   %(t_window)s %(no_t_samples)i   |dble: t_window; int: no_t_samples;
  #   #-------------------------------------------------------------------------------
  #   #
  #   #        WAVENUMBER INTEGRATION PARAMETERS
  #   #        =================================
  #   # 1. relative accuracy (0.01 for 1%% error) for numerical wavenumber integration;
  #   #-------------------------------------------------------------------------------
  #   %(accuracy)s                           |dble: accuracy;
  #   #-------------------------------------------------------------------------------
  #   #
  #   #        OUTPUTS A: DISPLACEMENT
  #   #        =======================
  #   # 1. select the 3 displacement time series (1/0 = yes/no)
  #   # 2. file names of these 3 time series
  #   #-------------------------------------------------------------------------------
  #   %(sw_t_files_1_3)s                                        |int: sw_t_files(1-3);
  #   %(t_files_1_3)s                                   |char: t_files(1-3);
  #   #-------------------------------------------------------------------------------
  #   #
  #   #        OUTPUTS B: STRAIN TENSOR & TILT
  #   #        ===============================
  #   # 1. select strain time series (1/0 = yes/no): Ezz, Err, Ett, Ezr, Ert, Etz
  #   #    (6 tensor components) and Tr (= -dur/dz, the radial component of the
  #   #    vertical tilt). Note Tt can be derived from Etz and Ut
  #   # 2. file names of these 7 time series
  #   #-------------------------------------------------------------------------------
  #   %(sw_t_files_4_10)s      |int: sw_t_files(4-10);
  #   %(t_files_4_10)s |char: t_files(4-10);
  #   #-------------------------------------------------------------------------------
  #   #
  #   #        OUTPUTS C: PORE PRESSURE & DARCY VELOCITY
  #   #        =========================================
  #   # 1. select pore pressure and Darcy velocity time series (1/0 = yes/no):
  #   #    P (excess pore pressure), Vz, Vr, Vt (3 Darcy velocity components)
  #   # 2. file names of these 4 time series
  #   #-------------------------------------------------------------------------------
  #   %(sw_t_files_11_14)s                              |int: sw_t_files(11-14);
  #   %(t_files_11_14)s                         |char: t_files(11-14);
  #   #-------------------------------------------------------------------------------
  #   #
  #   #        GLOBAL MODEL PARAMETERS
  #   #        =======================
  #   # 1. switch for surface conditions:
  surface.bc <- match.arg(surface.bc)
  #   #    0 = without free surface (whole space),
  #   #    1 = unconfined free surface (p = 0),
  #   #    2 = confined free surface (dp/dz = 0).
  isurfcon <- boundary.condition(surface.bc)$BC.val
  #switch(surface.bc, whole=0L, unconfined=1L, confined=2L)
  #   # 2. number of data lines of the layered model (<= lmax as defined in
  #   #    "peglobal.h") (see Note below)
  no_model_lines <- nrow(lith)
  #   #-------------------------------------------------------------------------------
  message(surface.bc)
  print(paste(isurfcon, inpstring(isurfcon=NA, cls="int")))
  #   %(isurfcon)i                   |int: isurfcon
  message(paste(no_model_lines,"model line(s)"))
  print(paste(no_model_lines, inpstring(no_model_lines=NA, cls="int")))
  #   %(no_model_lines)i             |int: no_model_lines;
  #   #-------------------------------------------------------------------------------
  #   #
  #   #        MULTILAYERED MODEL PARAMETERS
  #   #        =============================
  #   #
  #   # no depth[m] mu[Pa]    nu    nu_u   B     D[m^2/s]   Explanations
  #   #-------------------------------------------------------------------------------
  message(paste("model:"))
  #
  mdl.lines <- apply(lith, 1, paste, collapse=" ")
  print(paste(c("#", colnames(lith)), collapse=" "))
  print(as.data.frame(mdl.lines), row.names=FALSE)
  #for (l in mdl.lines) print(l)
  #   %(model)s
  #   #--------------------------end of all inputs------------------------------------
  #   
  #   Note for the model input format and the step-function approximation for models
  #   with depth gradient:
  #     
  #     The surface and the upper boundary of the half-space as well as the interfaces
  #   at which the poroelastic parameters are continuous, are all defined by a single
  #   data line; All other interfaces, at which the poroelastic parameters are
  #   discontinuous, are all defined by two data lines (upper-side and lower-side values). 
  #   This input format would also be needed for a graphic plot of the
  #   layered model. Layers which have different parameter values at top and bottom,
  #   will be treated as layers with a constant gradient, and will be discretised to a
  #   number of homogeneous sublayers. Errors due to the discretisation are limited
  #   within about 5%% (changeable, see peglobal.h).
}
#inp(1)

#' Generate an input file suitable to \code{poel}
#' @param x object
#' @param ... additional parameters
#' @export
#' @family Input-file
generator <- function(x, ...){ UseMethod("generator") }
#' @rdname generator
#' @export
generator.inp <- function(x, ...){
  #
}

#' Input line formatter
#' @details The user should not need to use this function.
#' @param ... objects to format
#' @param cls class of the input line
#' @param ns integer; the beginning index
#' @param ne integer; the ending index
#' @export
#' @family Input-file
#|dble: s_start_depth, s_end_depth, s_radius;
#|char: t_files(1-3);
#' @examples
#' inpstring(a=1)
#' inpstring(a=1,b=1)
inpstring <- function(..., cls=c("int", "dble", "char"), ns=NULL, ne=NULL){
      cls <- match.arg(cls)
      vars <- list(...)
      nms <- paste(names(vars), collapse=", ")
      nums <- if (!is.null(ns)){
	      sprintf("(%s);",paste(c(ns, ne), collapse="-"))
      } else {
	      ";"
      }
      Id <- sprintf("|%s: %s%s", cls, nms, nums)
      return(Id)
}
