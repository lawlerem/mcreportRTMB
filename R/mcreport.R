#' Make a quantity in an RTMB objective function available to mcreport
#' 
#' @param x 
#'     The quantity of which to take samples
#' 
#' @return 
#'     Used for its side effects when running mcreport
#' 
#' @rdname MCREPORT-macro
#' @export
MCREPORT<- function(x) internal_MCREPORT(x)
internal_MCREPORT<- function(x) invisible()

mcreporter<- function() {
    ans<- list()
    set<- function(nm, x) {
        ans[[nm]]<<- x
        NULL
    }
    report<- function(x) {
        # nm here gives the code that was inside the MCREPORT()
        #   that was called inside the objective function
        nm<- x |> 
            substitute() |>
            substitute() |>
            eval(parent.frame()) |>
            deparse()
        set(nm, x)
    }
    result<- function() ans
    clear<- function() ans<<- list()
    environment()
}

single_mcreport<- function(
        obj, 
        par = obj$env$last.par, 
        ...
    ) {
    package_env<- "mcreportRTMB" |> asNamespace()
    MCREPORT_ENV<- mcreporter()
    assignInMyNamespace(
        "internal_MCREPORT",
        MCREPORT_ENV$report
    )
    on.exit({
        rm(MCREPORT_ENV)
        assignInMyNamespace(
            "internal_MCREPORT",
            function(x) invisible()
        )
    })
    p<- obj$env$parList(par = par)
    attr(obj$env$data, "func")(p)
    MCREPORT_ENV$result()
}

#' Sample MCREPORTed quantities from the quasi-posterior of an RTMB objective 
#'     function
#' 
#' @param obj 
#'     The RTMB objective function.
#' @param replicates
#'     The number of samples.
#' @param sdr 
#'     The output of sdreport(obj, getJointPrecision = TRUE). If sdr is missing 
#'     or the supplied sdr does not have the joint precision matrix, sdreport 
#'     is called first.
#' @param parallel 
#'     The number of samples to run in parallel. Tries to use parallel::mclapply
#'     if parallal > 1, otherwise it uses lapply. Note that parallel::mclapply
#'     does not work on the Windows operating system.
#' @param silent 
#'     Should printing of progress be suppressed?
#' @param pivot
#'     See Matrix::chol.
#' 
#' @return 
#'     A list containing the samples of each MCREPORTed variable.
#'     If a variable is a vector, matrix, or array then it will be reported as
#'         an array with dimensions c(dim(var), replicates), else it will be
#'         returned as a list of length replicates
#'     The first replicate will be reported using the estimated parameters.
#'     The first replicate will be reported using the estimated parameters.
#' 
#' @export
mcreport<- function(
        obj,
        replicates = 100,
        sdr,
        parallel = 1,
        silent = FALSE,
        pivot = FALSE
    ) {
    need_sdr<- FALSE
    has_re<- length(obj$env$random) > 0
    if( missing(sdr) ) {
        need_sdr<- TRUE
    } else if( !("jointPrecision" %in% names(sdr)) & has_re ) {
        need_sdr<- TRUE
    }
    if( need_sdr ) {
        if( !silent ) {
            cat(
                paste0(
                    "Computing sdreport. "
                )
            )
        }
        sdr<- obj |> sdreport(getJointPrecision = TRUE)
    }
    if( !has_re ) sdr$jointPrecision<- sdr$cov.fixed |> solve()
    par_mean<- obj$env$last.par.best
    if( !silent ) {
        cat(
            paste0(
                "Computing precision cholesky decomposition. "
            )
        )
    }
    par_jpl<- sdr$jointPrecision |> Matrix::chol(pivot = pivot)
    if( !silent ) {
        cat(
            paste0(
                "Sampling parameter values.\n"
            )
        )
    }
    par_replicates<- (replicates * length(par_mean)) |> 
        rnorm() |>
        matrix(nrow = length(par_mean), ncol = replicates) |>
        Matrix::solve(
            a = par_jpl,
            b = _
        ) |>
        as.matrix()
    par_replicates<- (par_mean + par_replicates) |> t()
    par_replicates<- par_mean |> rbind(par_replicates)

    if( (parallel > 1) && requireNamespace("parallel", quietly = TRUE) ) {
        lapplyfn<- parallel::mclapply
    } else {
        lapplyfn<- lapply
    }
    mc_replicates<- seq(replicates + 1) |> 
        lapplyfn(
            function(i, ...) {
                if( !silent ) {
                    cat(
                        paste0(
                            "\rGetting mcreplicate: (",
                            i, " / ", replicates + 1, ")"
                        )
                    )
                    if( i == (replicates + 1) ) cat("\n")
                }
                return( obj |> single_mcreport(par_replicates[i, ]) )
            },
            mc.cores = parallel,
            mc.preschedule = FALSE
        )
    var_names<- mc_replicates[[1]] |> names()
    mc_replicates<- var_names |> 
        lapply(
            function(var) {
                x<- mc_replicates |> lapply(`[[`, var)
                try(
                    x<- x |> abind::abind(rev.along = 0),
                    silent = TRUE
                )
                return(x)
            }
        )
    names(mc_replicates)<- var_names
    mc_replicates<- c(
        par_replicates |> t() |> list(),
        mc_replicates
    )
    return( mc_replicates )
}
