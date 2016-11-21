#' @title Read .t files -- poel-result files
#' @description Read-in and maniputate .t files, the standard format used in poel output
#' @inheritParams read_inp
#' @seealso \code{\link{read_inp}}
#' @examples
#' tfi <- system.file('pp.t.gz', package='poel')
#' res <- read_t(tfi)
#' print(res)
read_t <- function(fi, version=2012){

    version <- match.arg(as.character(version), c('2012','2006'))
    if (version == '2006') stop('Parsing of version-',version," input files is currently unsupported.")

    # Try and find out what variable are being loaded
    readr::read_lines(fi, skip = 0, n_max = 1) %>% trimws -> oVars
    oVars %>% strsplit(., "  ") %>% unlist -> Vars
    if (length(Vars) == 1) oVars %>% strsplit(., " ") %>% unlist -> Vars

    # Load data
    readr::read_table(fi, col_names=TRUE) -> Dat

    # Get time values
    Time.var <- Vars[1]
    #   distance-depth values
    ZR.vars <- matrix(unlist(strsplit(Vars[-1], "_")), ncol=3, byrow=TRUE)
    #   and quanity in file
    Quant <- unique(ZR.vars[,1])

    Ts <- as.vector(Dat[,Time.var])
    Tinds <- seq_along(Ts)
    Zinds <- unique(as.numeric(gsub("^Z","",ZR.vars[,2])))
    Rinds <- unique(as.numeric(gsub("^R","",ZR.vars[,3])))

    # internal function to prepare the indices for numeric conversion
    .prep_index <- function(x){
        X <- as.numeric(gsub("^R", " ", gsub("^Z", "", unlist(strsplit(x, "_"))[2:3])))
        names(X) <- c("Z","R")
        return(X)
    }

    Inds <- list(Quant=Quant, Time=Ts, Ti=Tinds, Zi=Zinds, Ri=Rinds)

    # Tidy-up the data
    tidyr::gather(Dat, ZR, Value, -1) -> Dat.tdy

    # Go through each ZR and use .prep_index
    plyr::ldply(Dat.tdy$ZR, .prep_index) -> .ZR.
    .ZR.$Quantity <- Quant

    cbind(Dat.tdy, .ZR.) %>% tbl_df %>% dplyr::select(., 1, ZR, Z, R, Quantity, Value) -> Dat.tdy

    T.data <- list(File=fi,
                   Version = version,
                   Quantity = Quant,
                   Inds = Inds,
                   Dat = Dat.tdy)
    class(T.data) <- c("poel-t", version)
    return(T.data)
}

#' @title Methods for 'poel-t' class
#' @description Methods for 'poel-t' class
#' @rdname poel-t-methods
#' @inheritParams read_t
#' @param n integer; the maximum number of indices to print
#' @export
#' @seealso \code{\link{read_t}} \code{\link{read_inp}}
`print.poel-t` <- function(x, n=100){
    message("  poel results (.t file): ", x[['File']])
    message("      simulation indices: ")
    print(x[['Inds']], max = as.integer(n))
    message("      simulation values (", x[['Quantity']], "):")
    print(tbl_df(x[['Dat']]))
    return(invisible(x))
}
