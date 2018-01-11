
#' Title
#'
#' @param param
#' @param prec
#' @param airt
#' @param area
#' @param obs
#' @param warmup
#' @param ep
#' @param elev_zones
#' @param pelev
#' @param telev
#' @param incon
#' @param FUN_gof_args
#' @param pcalt_applied
#' @param tcalt_applied
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
                        elev_zones,
                        obs,
                        warmup,
                        pelev = NULL,
                        telev = NULL,
                        incon = NULL,
                        FUN_gof = hydroGOF::NSE,
                        FUN_gof_args = NULL,
                        pcalt_applied = FALSE,
                        tcalt_applied = FALSE,
                        as_hydromod = TRUE) {

  # TODO: hbv_out as vector?
  if (tcalt_applied == FALSE)
    airt <- apply_tcalt(series = airt, elev_zones = elev_zones, telev = telev,
                        lapse_rate = param[16])
  if (pcalt_applied == FALSE)
    prec <- apply_pcalt(series = prec, elev_zones = elev_zones, pelev = pelev,
                        lapse_rate = param[17])
  # run the model
  hbv_out <- TUWmodel::TUWmodel(prec, airt, ep, area, param[1:15])
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
    +    return(list(GoF=gof,model.out=sim))
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

  # param has to be 15-17 rows long, tcalt/pcalt is optional and added as 0 if missing
  e$param <- as.matrix(e$param)
  nrow_par <- nrow(e$param)
  ncol_par <- ncol(e$param)
  if ((nrow_par > 17  || nrow_par < 15) || (ncol_par > 2 || ncol_par < 1)) {
    stop("Wrong number of parameters in param")
  } else if (nrow_par < 17) {
    nrow_missing <- 17-nrow_par
    names_missing <- list(c("tcalt","pcalt")[1:nrow_missing])
    e$param <- rbind(e$param, matrix(0,nrow_missing,ncol_par,dimnames=names_missing))
  }

  ts_names <- c("prec", "airt", "ep", "obs")
  not_numeric = names(Filter(function(x)
    !is.numeric(x), e[ts_names]))
  if (length(not_numeric) > 0) {
    not_numeric <- paste0(not_numeric, collapse = ", ")
    stop("The following arguments have the wrong class (must be vectors, matrices or zoo objects):",
         not_numeric)
  }
  # zones checking
  if (!is.null(e$elev_zones) &&  length(e$elev_zones) != length(e$area))
    stop("Elevation zone and area must have the same length")
  if (sum(e$area) != 1)
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
      stop("The following time series end after \"to\":", no_end)
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

  # apply {p,t}calt if not used for calibration. if set to zero,
  # the input series will be used as is
  tcalt <- unique(e$param[16, ])
  if (any(tcalt != 0) && (is.null(e$telev) || is.null(e$elev_zones)))
    stop("If tcalt is used (param[16,] != 0), Parameters \"telev\" and \"elev_zones\" must be set.")

  if (length(tcalt) == 1) {
    e$airt <- apply_tcalt( series = e$airt,
                                   elev_zones = e$elev_zones, telev = e$telev,
                                   lapse_rate = tcalt)
    e$tcalt_applied = TRUE
  }
  pcalt <- unique(e$param[17, ])
  if (any(pcalt != 0) && (is.null(e$pelev) || is.null(e$elev_zones)))
    stop("If pcalt is used (param[17,] != 0), Parameters \"pelev\" and \"elev_zones\" must be set.")

  if (length(pcalt) == 1) {
    e$prec <- apply_pcalt(series = e$prec,
                                  elev_zones = e$elev_zones, pelev = e$pelev,
                                  lapse_rate = pcalt)
    e$pcalt_applied = TRUE
  }

  return(e)
}



#' Title
#'
#' @param series
#' @param elev_zones
#' @param lapse_rate
#' @param pelev
#'
#' @return
#' @export
#' @keywords internal
#'
#' @examples
apply_pcalt <- function(series, elev_zones, pelev, lapse_rate) {
  # transform from percent to decimal proportion
  lapse_rate <- round((lapse_rate / 100),2)
  # if (lapse_rate < 1)
  #   warning("pcalt should be specified in decimal proportions, but the given values was above 1(",lapse_rate,"). Input was assumed to be in percent and divided by 100.")
  return(apply_lapse_rate(series, elev_zones, pelev, lapse_rate,
                          type = "rel", neg.tozero = TRUE))
}

#' Title
#'
#' @param series
#' @param elev_zones
#' @param lapse_rate
#' @param telev
#'
#' @return
#' @export
#' @keywords internal
#' @examples
apply_tcalt <- function(series, elev_zones, telev, lapse_rate) {
  # make sure tcalt is negative
  lapse_rate <- abs(lapse_rate) * -1
  return(apply_lapse_rate(series, elev_zones, telev, lapse_rate,
                          type = "abs", neg.tozero = FALSE))
}


#' Title
#'
#' @param series
#' @param elev_zones
#' @param elev_ref
#' @param lapse_rate
#' @param type
#' @param neg.tozero
#'
#' @return
#' @keywords internal
#' @examples
apply_lapse_rate <- function(series, elev_zones, elev_ref, lapse_rate,
                             type = "rel", neg.tozero = TRUE) {
  # TODO: does it make sense to allow lapse_rate for multi-columns series?
  if (lapse_rate == 0)
    return(series)
  if (NCOL(series) > 1)
    stop("If using pcalt/tcalt, the corresponding time series must be univariate")
  nelev_zones <- length(elev_zones)
  series <- matrix(rep(series, nelev_zones), ncol = nelev_zones)
  lapse_rates <- (elev_zones - elev_ref) / 100 * lapse_rate
  if (type == "rel") {
    lapse_rates <- lapse_rates + 1
    if (neg.tozero) {
      lapse_rates[lapse_rates < 0] <- 0
    }
    series_adj <- series %*% diag(lapse_rates)
  } else if (type == "abs") {
    series_adj <- sweep(series, 2, FUN = "+", lapse_rates)
    if (neg.tozero) {
      series_adj[series_adj < 0] <- 0
    }
  } else {
    stop("Type must be either \"rel\" (relative difference) or \"abs\" (absolute difference)")
  }
  return(round(series_adj, 5))
}
