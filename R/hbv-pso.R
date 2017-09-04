#' HBV Modelling with Particle Swarm Optimization
#'
#' Performs particle swarm optimization of the HBV hydrological model by wrapping \link[hydroPSO]{hydroPSO} and \link[TUWmodel]{TUWmodel}.
#' @param prec Precipitation input (mm/day) as zoo, matrix or numerical. If multivariate, each variable is the input for one zone.
#' @param airt Air Temperature input (degC/day) as zoo, matrix or numerical. If multivariate, each variable is the input for one zone.
#' @param ep Potential Evapotranspiration (mm/day) as zoo, matrix or numerical. If multivariate, each variable is the input for one zone.
#' @param area If input data is distributed into zones (multivariate zoo/matrix), a vector of the decimal proportion of area for each zone.
#' @param param Parameters as two-column (min,max) matrix or dataframe if optimization should be performed, and as vector otherwise.
#' \enumerate{
#' \item \code{SCF} snow correction factor (0.9-1.5);
#' \item \code{DDF} degree day factor (0.0-5.0 mm/degC/timestep);
#'\item \code{Tr} threshold temperature above which precipitation is rain (1.0-3.0 degC);
#' \item \code{Ts} threshold temperature below which precipitation is snow (-3.0-1.0 degC);
#' \item \code{Tm} threshold temperature above which melt starts (-2.0-2.0 degC);
#' \item \code{LPrat} parameter related to the limit for potential evaporation (0.0-1.0);
#' \item \code{FC} field capacity, i.e., max soil moisture storage (0-600 mm);
#' \item \code{BETA} the non linear parameter for runoff production (0.0-20.0);
#' \item \code{k0} storage coefficient for very fast response (0.0-2.0 timestep);
#' \item \code{k1} storage coefficient for fast response (2.0-30.0 timestep);
#' \item \code{k2} storage coefficient for slow response (30.0-250.0 timestep);
#' \item \code{lsuz} threshold storage state, i.e., the very fast response start if exceeded (1.0-100.0 mm);
#' \item \code{cperc} constant percolation rate (0.0-8.0 mm/timestep);
#' \item \code{bmax} maximum base at low flows (0.0-30.0 timestep);
#' \item \code{croute} free scaling parameter (0.0-50.0 timestep2/mm);
#' \item \code{tcalt} Lapse rate to adjust the temperature data by elevation zone (ºC/100m, decreasing with elevation)
#' \item \code{pcalt} Lapse rate to adjust the precipitation data by elevation zone (%/100m, increasing with elevation)
#' } The last two parameters are optional and used to transform the temperature/precipitation input (instead of being passed on to TUWmodel). To disable pcalt/tcalt, set them to zero or ommit from param.
#' See the example povided at \link[hbvPSO]{tuwmodel_params_default} for an example with the default ranges as specified in \link[TUWmodel]{TUWmodel}.
#' @param obs Observed Discharge (mm/day) as zoo or numerical
#' @param from Start of the modelling period (including warmup) as Date or string in standard date format. Requires input datasets to be zoo objects.
#' @param to End of the modelling period as Date or string in standard date format. Requires input datasets to be zoo objects.
#' @param warmup Warmup phase which is removed before calculating goodness of fit. Can be given as numeric (days removed from the model start date) or date (as Date object or string in default format which can be cast to Date by as.Date). If given as date, it marks the end of the warmup period.
#' @param telev Reference Elevation for the air temperature input, used to adjust temperature by tcalt.
#' @param pelev Reference Elevation for the precipitation input, used to adjust precipitation by pcalt.
#' @param elev_zones Vector of mean elevation for each zone. Only required if tcalt/pcalt is used.
#' @param incon vector/matrix of initial conditions for the model (\code{ncol} = number of zones):
#' \code{SSM0} soil moisture (mm);
#' \code{SWE0} snow water equivalent (mm);
#' \code{SUZ0} initial value for fast (upper zone) response storage (mm);
#' \code{SLZ0} initial value for slow (lower zone) response storage (mm)
#' @param outpath Path to the directory storing the output files as string. If not NULL, the following is set in in hydroPSO's control list: \code{drty.out=outpath} and \code{write2disk=TRUE}
#' @param hydroPSO_args Arguments passed on to \link[hydroPSO]{hydroPSO}
#' @param FUN_gof The function used to calculate goodness of fit. Must take sim and obs as first arguments.
#' @param FUN_gof_args Further arguments passed on to \code{FUN_gof}
#' @param plotting Toggles plotting of the results (with \link[hydroPSO]{plot_results} if optimization is performed, otherwise with \link[hydroPSO]{plot_out}). As alternative to \code{TRUE}, a list of arguments for the respective plotting functions can be provided.
#'
#' @return A list of the following items:
#' \enumerate{
#' \item \code{sim} simulated runoff (mm/day) of the best model run
#' \item \code{obs} observed runoff (mm/day)
#'\item \code{gof} goodness of fit of the best model run
#' \item \code{pso_out} \link[hydroPSO]{hydroPSO} output
#' \item \code{hbv_out} \link[TUWmodel]{TUWmodel} output from the best model run
#' }
#' @export
#' @importFrom hydroGOF NSE ggof
#' @import zoo
#'
#' @examples
hbv_pso <- function(prec = NULL,
                    airt = NULL,
                    ep = NULL,
                    area = 1,
                    param = hbvPSO::tuwmodel_params_default,
                    obs = NULL,
                    from = NULL,
                    to = NULL,
                    warmup = 0,
                    telev = NULL,
                    pelev = NULL,
                    elev_zones = NULL,
                    incon = NULL,
                    outpath=NULL,
                    hydroPSO_args = NULL,
                    FUN_gof = hydroGOF::NSE,
                    FUN_gof_args = NULL,
                    plotting=FALSE) {
  # FIXME: what if bestrun q returns NA
  # FIXME: obs/airt/ep/prec == NULL?
  # FIXME: why are my plot lines so thick
  # TODO: support specifying a validation period?
  # TODO: support parameter zones - convert pars to matrix if not vector of len15?
	# TODO: if mix of zoo/numeric for ts, somehow convert all to zoo with index of e.g. obs?
  # TODO: default control values?
  # TODO: no ts plot if obs!=zoo from plot_results (works for plot_out)
  if (is.null(hydroPSO_args))
    hydroPSO_args <- list()
  # simplify argument handling: plotting is a list if we should plot, FALSE otherwise
  if (isTRUE(plotting))
    plotting = list()
  args_list <- as.list(environment())
  if (!is.null(outpath) && !dir.exists(outpath)) {
      dir.create(outpath,recursive = TRUE)
  }
  # input checking
  if(is.list(plotting) && is.null(outpath))
    stop("If plotting is enabled, \"outpath\" must be set")

  if (NCOL(obs) != 1) {
    stop("Observed discharge (\"obs\") must be univariate")
  }

  # param has to be 15-17 rows long, tcalt/pcalt is optional and added as 0 if missing
  param <- as.matrix(param)
  nrow_par <- nrow(param)
  ncol_par <- ncol(param)
  if ((nrow_par > 17  || nrow_par < 15) || (ncol_par > 2 || ncol_par < 1)) {
    stop("Wrong number of parameters in param")
  } else if (nrow_par < 17) {
    nrow_missing <- 17-nrow_par
    names_missing <- list(c("tcalt","pcalt")[1:nrow_missing])
    param <- rbind(param, matrix(0,nrow_missing,ncol_par,dimnames=names_missing))
  }

  ts_names <- c("prec", "airt", "ep", "obs")
  not_numeric = names(Filter(function(x)
    !is.numeric(x), args_list[ts_names]))
  if (length(not_numeric) > 0) {
    not_numeric <- paste0(not_numeric, collapse = ", ")
    stop("The following arguments have the wrong class (must be vectors, matrices or zoo objects):",
      not_numeric)
  }
  # zones checking
  if (!is.null(elev_zones) &&  length(elev_zones) != length(area))
    stop("Elevation zone and area must have the same length")
  if (sum(area) != 1)
    stop("The sum of \"area\" must be 1")

  wrong_dim <- names(Filter(function(x, n_area = length(area)) {
    ncols <- NCOL(x)
    return(ncols > 1 && ncols != n_area)
    }, args_list[ts_names]))

  if (length(wrong_dim) > 0) {
    wrong_dim <- paste0(wrong_dim, collapse = ", ")
    stop("The number of columns in the following time series do not match the \"area\" parameter:",
      wrong_dim)
  }

  if(!is.null(from)) {
    from <- tryCatch(as.Date(from), error = function(err) {NULL})
    if (is.null(from))
      stop("\from\" could not be cast to a date (using as.Date defaults)")
  }

  if(!is.null(to)) {
    to <- tryCatch(as.Date(to), error = function(err) {NULL})
    if (is.null(to))
      stop("\to\" could not be cast to a date (using as.Date defaults)")
  }
  ts_as_zoo <- all(sapply(args_list[ts_names], zoo::is.zoo))

  if (ts_as_zoo) {
    no_start <- names(Filter(function(x) start(x) > from, args_list[ts_names]))
    if (length(no_start) > 0) {
      no_start <- paste0(no_start, collapse = ", ")
      stop("The following time series begin after \"from\":", no_start)
    }
    no_end <- names(Filter(function(x) end(x) < to, args_list[ts_names]))
    if (length(no_end) > 0) {
      no_end <- paste0(no_end, collapse = ", ")
      stop("The following time series end after \"to\":", no_end)
    }
    # cut timeseries to start/end
    args_list[ts_names] = lapply(args_list[ts_names], FUN = window,
                                 start = from, end = to)

  } else if ((!is.null(from) || !is.null(to))) {
      stop("If \"from\" and \"to\" are used, all timeseries (\"prec\", \"airt\", \"ep\", \"obs\") must be zoo objects")
  }

  if (!is.numeric(warmup)) {
  	warmup <- tryCatch(as.Date(warmup,tz="UTZ"), error = function(err) {NULL})
  	if (is.null(warmup)) {
  		stop("Failed to convert \"warmup\" into date (must be numeric, Date or string which can be cast by \"as.Date\")")
  	} else if (!ts_as_zoo){
  		stop("If \"warmup\" is a date, all timeseries (\"prec\", \"airt\", \"ep\", \"obs\") must be zoo objects")
  	} else if (is.null(from)) {
  		stop("If \"warmup\" is a date, \"from\" must be set")
  	} else if (warmup < from) {
  	  stop("If \"warmup\" is a date, it must be later than \"from\"")
  	} else {
  	  warmup <- as.numeric(warmup - from)
  	  args_list$warmup <- warmup
  	}
  }

  if (length(unique(sapply(args_list[ts_names], NROW))) != 1)
    stop("All timeseries (\"prec\", \"airt\", \"ep\", \"obs\") must be of equal length")

  # remove the warmup period from observed Q
  if (warmup >= NROW(args_list$obs))
    stop("Warmup is longer than the input timeseries")
  obs <- args_list$obs[-(1:warmup)]
  args_list$obs <- obs

  # convert any zoo ts to vector/matrix
  args_list[ts_names] <- lapply(args_list[ts_names], FUN = zoo::coredata)

  # apply {p,t}calt if not used for calibration. if set to zero,
  # the input series will be used as is
  tcalt <- unique(param[16, ])
  if (any(tcalt != 0) && (is.null(telev) || is.null(elev_zones)))
    stop("If tcalt is used (param[16,] != 0), Parameters \"telev\" and \"elev_zones\" must be set.")

  if (length(tcalt) == 1) {
    args_list$airt <- apply_tcalt( series = args_list$airt,
                                   elev_zones = elev_zones, telev = telev,
                                   lapse_rate = tcalt)
    args_list$tcalt_applied = TRUE
  }
  pcalt <- unique(param[17, ])
  if (any(pcalt != 0) && (is.null(pelev) || is.null(elev_zones)))
    stop("If pcalt is used (param[17,] != 0), Parameters \"pelev\" and \"elev_zones\" must be set.")

  if (length(pcalt) == 1) {
    args_list$prec <- apply_pcalt(series = args_list$prec,
                                  elev_zones = elev_zones, pelev = pelev,
                                  lapse_rate = pcalt)
    args_list$pcalt_applied = TRUE
  }
  # Clean up argument list before passing to other functions
  args_list <- args_list[!names(args_list) %in% c("hydroPSO_args", "param",
                                                  "outpath", "from", "to","plotting")]
  args_list <- args_list[!sapply(args_list, is.null)]

  do_optimize <- ifelse(ncol(param)==2,TRUE,FALSE)
  # Run hydroPSO if we have a parameter range
  if (do_optimize) {
    if (!is.null(outpath)) {
      if (!"control" %in% names(hydroPSO_args))
        hydroPSO_args$control <- list()
      hydroPSO_args$control$write2disk <- TRUE
      hydroPSO_args$control$drty.out <- outpath
    }
    hydroPSO_args_default <- list(lower=param[, 1], upper = param[, 2], fn="hbv_single")
    hydroPSO_args <- c(hydroPSO_args,hydroPSO_args_default[!(names(hydroPSO_args_default) %in% names(hydroPSO_args))])

    # run hydroPSO
    hydroPSO_args <- c(args_list, hydroPSO_args)
    bestpar <- do.call(hydroPSO::hydroPSO, hydroPSO_args)
    # Re-run hbv with best parameter set and save the output
    args_list$param <- bestpar$par
  } else {
    bestpar <- NA
    args_list$param <- param
  }

  args_list$gof_only <- FALSE
  bestrun <- do.call(hbv_single, args_list)

  if (!is.null(outpath)) {
    # TODO: rename modelout / keep bestmodelout for single run?
    sim_file_name <- ifelse(do_optimize, "BestModel_out.txt","Model_out.txt")
    sim_file <- file.path(outpath, sim_file_name)
    obs_file <- file.path(outpath, "Observations.txt")
    if (zoo::is.zoo(obs)) {
    	sim <- zoo::zoo(bestrun$hbv_out$q, order.by=zoo::index(obs))
      write.zoo(obs, obs_file, col.names=FALSE)
      write.zoo(sim, sim_file, col.names=FALSE)
    } else {
    	sim <- bestrun$hbv_out$q
      write.table(obs, obs_file, col.names=FALSE, quote=FALSE)
      write.table(sim, sim_file, col.names=FALSE, row.names=FALSE)
    }
  }
  if (is.list(plotting)) {
    if(do_optimize) {
      plot_args <- list(sim=zoo::coredata(sim),obs=obs,drty.out=outpath,do.png=TRUE,MinMax="max",
                          beh.thr=0.0, ftype="dm", legend.pos="right")
      FUN_plot <- hydroPSO::plot_results
    } else {
      plot_args <- list(sim=zoo::coredata(sim), obs=obs, ptype="ts", MinMax="max", ftype="dm",
                            do.png=TRUE, png.fname=file.path(outpath,"ModelOut_vs_Obs.png"))
      FUN_plot <- hydroPSO::plot_out
    }
    plot_args <- c(plotting,plot_args[!(names(plot_args) %in% names(plotting))])
    do.call(FUN_plot,plot_args)
  }
  return(list(sim = sim, obs = obs, gof = bestrun$gof, pso_out = bestpar,
              hbv_out = bestrun$hbv_out))
}



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
#' @param gof_only
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
                       gof_only = TRUE) {


  if (tcalt_applied == FALSE)
    airt <- apply_tcalt(series = airt, elev_zones = elev_zones, telev = telev,
                        lapse_rate = param[16])
  if (pcalt_applied == FALSE)
    prec <- apply_pcalt(series = prec, elev_zones = elev_zones, pelev = pelev,
                        lapse_rate = param[17])

  # TUWmodel does not model by zone if prec is a vector. Therefore prec must be
  # converted to a matrix of appropriate dimensions if zones are used
  nelev_zones <- length(area)
  if (nelev_zones > 1 && NCOL(prec) == 1)
    prec <- matrix(rep(prec, nelev_zones), ncol = nelev_zones)
  # run the model
  hbv_out <- TUWmodel::TUWmodel(prec, airt, ep, area, param[1:15])
  # remove warmup period from simulated Q
  sim <- as.numeric(hbv_out$q)[-(1:warmup)]
  # if there are any NA in qsim, directly return NA as GoF
  if (any(is.na(sim)) && gof_only)
    return(NA)
  FUN_gof_args = c(list(sim = sim, obs = obs), FUN_gof_args)
  gof <- (do.call(FUN_gof, FUN_gof_args))
  # combine sim and obs with any additional args to FUN_gof
  if (gof_only) {
    return(gof)
  } else {
    # remove warmup from all components
    hbv_out$q <- sim
    idx <- c("qzones","q0","q1","q1","rain", "melt", "moist", "swe", "suz","slz")
    hbv_out[idx] <- lapply(hbv_out[idx], tail, n=length(sim))
    return(list(hbv_out=hbv_out,gof=gof))
  }
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
