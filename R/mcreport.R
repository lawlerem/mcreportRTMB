#' Make a quantity in an RTMB objective function available to mcreport
#' 
#' @param x The quantity to take samples of
#' 
#' @return Used for its side effects when running mcreport
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
        nm<- deparse(
            eval(
                substitute(
                    substitute(
                        x
                    )
                ),
                parent.frame()
            )
        )
        set(nm, x)
    }
    result<- function() ans
    clear<- function() ans<<- list()
    environment()
}

single_mcreport<- function(obj, par = obj$env$last.par, ...) {
    package_env<- asNamespace("mcreportRTMB")
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

#' Sample MCREPORTed quantities from the quasi-posterior of an RTMB objective function
#' 
#' @param obj The RTMB objective function.
#' @param n The number of samples.
#' @param sdr The output of sdreport(obj, getJointPrecision = TRUE). If sdr is missing or the supplied sdr does not have the joint precision matrix, sdreport is called first.
#' @param parallel The number of samples to run in parallel. Tries to use parallel::mclapply if parallal > 1, otherwise it uses lapply. Note that parallel::mclapply does not work on the Windows operating system.
#' @param trace Should progress be printed?
#' 
#' @return A list containing the samples of each MCREPORTed variable
#' 
#' @export
mcreport<- function(
        obj,
        n = 100,
        sdr,
        parallel = 1,
        trace = TRUE
    ) {
    need_sdr<- missing(sdr) | (
        !missing(sdr) && !("jointPrecision" %in% names(sdr))
    )
    if( need_sdr ) sdr<- sdreport(obj, getJointPrecision = TRUE)
    par_mean<- obj$env$last.par.best
    par_cov<- solve(sdr$jointPrecision)
    par_replicates<- MASS::mvrnorm(
        n = n,
        mu = par_mean,
        Sigma = par_cov
    )
    if( (parallel > 1) && requireNamespace("parallel", quietly = TRUE) ) {
        lapplyfn<- parallel::mclapply
    } else {
        lapplyfn<- lapply
    }
    mc_replicates<- lapplyfn(
        seq(n),
        function(i, ...) {
            if( trace ) {
                cat(
                    paste0(
                        "(",
                        i,
                        " / ",
                        n,
                        ")\n"
                    )
                )
            }
            single_mcreport(obj, par_replicates[i, ])
        },
        mc.cores = parallel,
        mc.preschedule = FALSE
    )
    var_names<- names(mc_replicates[[1]])
    mc_replicates<- lapply(
        var_names,
        function(var) {
            x<- abind::abind(
                lapply(mc_replicates, `[[`, var),
                rev.along = 0
            )
            return(x)
        }
    )
    names(mc_replicates)<- var_names
    return(mc_replicates)
}