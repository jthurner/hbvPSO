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
#' See \link[hbvPSO]{tuwmodel_params_default} for an example with the default ranges as specified in \link[TUWmodel]{TUWmodel}.
#' @param obs Observed Discharge (mm/day) as zoo or numerical
#' @param from Start of the modelling period (including warmup) as Date or string in standard date format. Requires input datasets to be zoo objects.
#' @param to End of the modelling period as Date or string in standard date format. Requires input datasets to be zoo objects.
#' @param warmup Warmup phase which is removed before calculating goodness of fit. Can be given as numeric (days removed from the model start date) or date (as Date object or string in default format which can be cast to Date by as.Date). If given as date, it marks the start of the simulation period after warmup.
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
#' \dontrun{
#' # loading the example data from TUWmodel
#' data(example_TUWmodel, package="TUWmodel")
#' # extracting the input data for the lumped case (see TUWmodel examples):
#' # 1.) apply weighted means
#' prec <- apply(P_Vils, 1, weighted.mean, w=areas_Vils)
#' airt <- apply(T_Vils, 1, weighted.mean, w=areas_Vils)
#' ep <- apply(PET_Vils, 1, weighted.mean, w=areas_Vils)
#' # 2.) casting from named numeric vector to zoo
#' prec <- zoo(prec, order.by=as.Date(names(prec)))
#' airt <- zoo(airt, order.by=as.Date(names(airt)))
#' ep <- zoo(ep, order.by=as.Date(names(ep)))
#' obs <- zoo(Q_Vils, order.by=as.Date(names(Q_Vils)))
#' # setting up date range and warmup, limit max iterations to 500
#' control <- list(maxit=500)
#' from <- "1976-01-01"
#' to <- "1996-12-31"
#' warmup <- "1977-01-01"
#' # running hbv_PSO with default options, without plotting or parallel processing
#' res <- hbv_pso(prec=prec, airt=airt, ep=ep, obs=obs,hydroPSO_args = list(control=control))
#' }


hbv_pso <- function(prec = NULL,
                    airt = NULL,
                    ep = NULL,
                    area = 1,
                    param = hbvPSO::tuwmodel_params_default,
                    obs = NULL,
                    from = NULL,
                    to = NULL,
                    warmup = 0,
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

  # simplify argument handling: set up hydroPSO_args structure
  if (is.null(hydroPSO_args)) {
    hydroPSO_args <- list(control=list())
  } else if (is.null(hydroPSO_args$control)) {
    hydroPSO_args$control <- list()
  }
  if (is.null(outpath)) {
    hydroPSO_args$control$write2disk <- FALSE
  } else {
    hydroPSO_args$control$write2disk <- TRUE
    hydroPSO_args$control$drty.out <- outpath
  }
  # set up hydroPSO defaults
  hydroPSO_args$control$MinMax <- "max"

  # re-set defaults for FUN_gof and param if they were set to NULL
  if (is.null(FUN_gof)) {
    FUN_gof <- hydroGOF::NSE
    FUN_gof_args <- NULL
  }

  # force plotting to be a list if we should plot, FALSE otherwise
  if (isTRUE(plotting))
    plotting <- list()

  # capture arguments into list for processing & handing off to hydroPSO
  args_list <- as.list(environment())
  if (!is.null(outpath) && !dir.exists(outpath)) {
      dir.create(outpath,recursive = TRUE)
  }

  if(is.list(plotting) && is.null(outpath))
    stop("If plotting is enabled, \"outpath\" must be set")

  # input checking
  args_list <- validate_input(args_list)
  # obs == as.numeric(obs-warmup), obs_zoo == obs-warmup
  obs <- args_list$obs_zoo

  # Clean up argument list before passing to other functions
  param <- args_list$param
  args_list <- args_list[!names(args_list) %in% c("hydroPSO_args", "param",
                                                  "outpath", "from", "to",
                                                  "plotting", "obs_zoo")]
  args_list <- args_list[!sapply(args_list, is.null)]

  do_optimize <- ifelse(ncol(param)==2,TRUE,FALSE)
  # Run hydroPSO if we have a parameter range
  if (do_optimize) {
    hydroPSO_args_default <- list(lower=param[, 1], upper = param[, 2], fn="hbv_single")
    # merging user-provided arguments with default args
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

  args_list$as_hydromod <- FALSE
  bestrun <- do.call(hbv_single, args_list)

  if (zoo::is.zoo(obs)) {
    sim <- zoo::zoo(bestrun$hbv_out$q, order.by=zoo::index(obs))
  } else {
    sim <- bestrun$hbv_out$q
  }

  if (!is.null(outpath)) {
    # TODO: rename modelout / keep bestmodelout for single run?
    sim_file_name <- ifelse(do_optimize, "BestModel_out.txt","Model_out.txt")
    sim_file <- file.path(outpath, sim_file_name)
    obs_file <- file.path(outpath, "Observations.txt")
    if (zoo::is.zoo(obs)) {
      write.zoo(obs, obs_file, col.names=FALSE)
      write.zoo(sim, sim_file, col.names=FALSE)
    } else {
      write.table(obs, obs_file, col.names=FALSE, quote=FALSE)
      write.table(sim, sim_file, col.names=FALSE, row.names=FALSE)
    }
  }
  if (is.list(plotting)) {
    if(do_optimize) {
      param.names <- rownames(param)[param[,1] != param[,2]]
      # if gof.name is not defined, set it to the name of the objective function
      if (is.null(plotting$gof.name)) {
        plotting$gof.name <- gsub("^.*:","", deparse(substitute(FUN_gof)))
      }
      plot_args <- list(sim=zoo::coredata(sim), obs=obs, drty.out=outpath,
                        param.names = param.names, do.png=TRUE,MinMax="max",
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
