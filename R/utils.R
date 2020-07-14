
#' Title
#'
#' @param param
#' @param prec
#' @param airt
#' @param area
#' @param obs
#' @param warmup
#' @param ep
#' @param incon
#' @param FUN_gof_args
#' @param as_hydromod
#' @param FUN_gof
#'
#' @return
#' @import TUWmodel
#' @examples
#' @keywords internal
hbv_single <-  function(prec,
                        airt,
                        ep,
                        area,
                        param,
                        obs,
                        warmup,
                        incon = NULL,
                        FUN_gof = hydroGOF::NSE,
                        FUN_gof_args = NULL,
                        as_hydromod = TRUE,
                        par_fixed = NULL,
                        disable_tr = FALSE) {

  # TODO: hbv_out as vector?
  # Add parameters excluded from optimization
  if (!is.null(par_fixed)) {
    param <- replace(par_fixed, is.na(par_fixed),param)
  }
  # if tr should not be used, set it to the value of ts
  if (disable_tr) {
    param[3] <- param[4]
  }
  # run the model
  hbv_out <- TUWmodel::TUWmodel(prec, airt, ep, area, param)
  # remove warmup period from simulated Q
  # TODO: This fails for multi-zone!! use tail instead?
  sim <- as.numeric(hbv_out$q)[-(1:warmup)]
  # if there are any NA in qsim, directly return NA as GoF
  if (any(is.na(sim)) && as_hydromod)
    return(list(GoF=NA,model.out=NA))
  FUN_gof_args = c(list(sim = sim, obs = obs), FUN_gof_args)
  gof <- (do.call(FUN_gof, FUN_gof_args))
  # combine sim and obs with any additional args to FUN_gof
  if (as_hydromod) {
        return(list(GoF=gof,model.out=sim))
  } else {
    # remove warmup from all components
    hbv_out$q <- sim
    idx <- c("qzones","q0","q1","q2","rain", "snow", "melt", "moist", "swe", "suz","slz")
    hbv_out[idx] <- lapply(hbv_out[idx], tail, n=length(sim))
    return(list(hbv_out=hbv_out,gof=gof))
  }
}



validate_input <- function(e) {

  if (is.null(e$param)) {
    e$param <- hbvPSO::tuwmodel_params_default
  }

    if (NCOL(e$obs) != 1) {
    stop("Observed discharge (\"obs\") must be univariate")
  }

  # param has to be 15 rows long

  nrow_par <- NROW(e$param) + sum(!is.na(e$par_fixed))
  ncol_par <- NCOL(e$param)
  if ((nrow_par != 15) || (ncol_par > 2 || ncol_par < 1)) {
    stop("Wrong number of parameters in param")
  }

  ts_names <- c("prec", "airt", "ep", "obs")
  not_numeric = names(Filter(function(x)
    !is.numeric(x), e[ts_names]))
  if (length(not_numeric) > 0) {
    not_numeric <- paste0(not_numeric, collapse = ", ")
    stop("The following arguments have the wrong class (must be vectors, matrices or zoo objects):",
         not_numeric)
  }
  # area fractions must sum up to one
  if (!isTRUE(all.equal(sum(e$area),1,tolerance=0.0001)))
    stop("The sum of \"area\" must be 1")

  wrong_dim <- names(Filter(function(x, n_area = length(e$area)) {
    ncols <- NCOL(x)
    return(ncols > 1 && ncols != n_area)
  }, e[ts_names]))

  if (length(wrong_dim) > 0) {
    wrong_dim <- paste0(wrong_dim, collapse = ", ")
    stop("The number of columns in the following time series do not match the \"area\" parameter:",
         wrong_dim)
  }

  if(!is.null(e$from)) {
    e$from <- tryCatch(as.Date(e$from), error = function(err) {NULL})
    if (is.null(e$from))
      stop("\from\" could not be cast to a date (using as.Date defaults)")
  }

  if(!is.null(e$to)) {
    e$to <- tryCatch(as.Date(e$to), error = function(err) {NULL})
    if (is.null(e$to))
      stop("\"to\" could not be cast to a date (using as.Date defaults)")
  }

  ts_as_zoo <- all(sapply(e[ts_names], zoo::is.zoo))

  if (ts_as_zoo) {
    no_start <- names(Filter(function(x) start(x) > e$from, e[ts_names]))
    if (length(no_start) > 0) {
      no_start <- paste0(no_start, collapse = ", ")
      stop("The following time series begin after \"from\":", no_start)
    }
    no_end <- names(Filter(function(x) end(x) < e$to, e[ts_names]))
    if (length(no_end) > 0) {
      no_end <- paste0(no_end, collapse = ", ")
      stop("The following time series end before \"to\":", no_end)
    }
    # cut timeseries to start/end
    e[ts_names] = lapply(e[ts_names], FUN = window,
                                 start = e$from, end = e$to)

  } else if ((!is.null(e$from) || !is.null(e$to))) {
    stop("If \"from\" and \"to\" are used, all timeseries (\"prec\", \"airt\", \"ep\", \"obs\") must be zoo objects")
  }

  if (!is.numeric(e$warmup)) {
    warmup <- tryCatch(as.Date(e$warmup,tz="UTZ"), error = function(err) {NULL})
    if (is.null(warmup)) {
      stop("Failed to convert \"warmup\" into date (must be numeric, Date or string which can be cast by \"as.Date\")")
    } else if (!ts_as_zoo){
      stop("If \"warmup\" is a date, all timeseries (\"prec\", \"airt\", \"ep\", \"obs\") must be zoo objects")
    } else if (is.null(e$from)) {
      stop("If \"warmup\" is a date, \"from\" must be set")
    } else if (warmup < e$from) {
      stop("If \"warmup\" is a date, it must be later than \"from\"")
    } else {
      e$warmup <- as.numeric(warmup - e$from)
    }
  }

  if (length(unique(sapply(e[ts_names], NROW))) != 1)
    stop("All timeseries (\"prec\", \"airt\", \"ep\", \"obs\") must be of equal length")

  # remove the warmup period from observed Q
  if (e$warmup >= NROW(e$obs))
    stop("Warmup is longer than the input timeseries")
  e$obs <- e$obs[-(1:e$warmup)]

  e$obs_zoo <- e$obs

  # convert any zoo ts to vector/matrix
  e[ts_names] <- lapply(e[ts_names], FUN = zoo::coredata)

  return(e)
}
