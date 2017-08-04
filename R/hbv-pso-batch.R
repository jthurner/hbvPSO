
#' Batch-run multiple Particle Swarm Optimizations of the HBV Model
#'
#' Batch-run multiple Particle Swarm Optimizations of the HBV Model using \link[hydroPSO]{hydroPSO} and \link[TUWmodel]{TUWmodel}.
#' @details
#' Runs multiple Particle Swarm Optimizations of HBV for different sets of input data (prec, airt, ep, obs, area, elevation_zones) and/or parameters.
#'
#' Each optimization run is given an identification string concatenated from the name of the input data folder, the name of the gof variable and an optional suffix. This identifier is used as folder name for the output files, in plot titles and file names, and in the result list returned by \code{hbv_pso_batch}.
#'
#' @section Directory layout:
#' \code{hbv_pso_batch} assumes you have set up your input data in the following way:
#' \itemize{
#' \item \code{base directory} A root directory containing the global parameter config and sub-directories for different sets of input data. The path to this directory is the first argument passed to \code{hbv_pso_batch}
#' \item \code{input data directories} Sub-directories to \code{base directory} containing sets of input data (prec,airt,ep,obs,area and elevation_zones).
#' \item \code{Data directory} Optional, sub-directory to \code{input data directories} containing a \code{Data} directory taken from HBV-Light (see "Input data definition" below).
#' }
#'
#' @section Parameter definition:
#' Parameters for the different optimization runs are defined by simple variable assignment, where the variable name corresponds to the parameter as used in \link[ittr]{hbv_pso}.
#' They can be set in the following locations, in ascending order of precedence (parameters assigned in any location will be overwritten by re-assignment in subsequent locations):
#' \enumerate{
#' \item \code{global_config.R} Default parameter definition for all optimization runs. Located at the root of the base directory
#' \item \code{arguments to hbv_pso_batch} Argument names represent the parameter name to be set (arguments will be unpacked into the environment). Setting individual parameters nested inside lists (e.g. hydroPSOA_args$control$maxit) is currently not possible, meaning any lists passed are not merged but replace the default value.
#' \item \code{localconf.R} Located inside any input data directory, parameters defined here overwrite all previous definitions and are local to the respective input data set.
#' }
#'
#' In addition to parameters described in \link[ittr]{hbv_pso}, the following parameters unique to \code{hbv_pso_batch} can be set:
#' \enumerate{
#' \item \code{suffix} String, optional. Will be appended to the optimization run identifier
#' \item \code{gof.name} String. Used as part of the optimization run identifier
#' \item \code{from_validation} Date or date string in default format, start of validation run with best parameter set found during optimization (optional)
#' \item \code{to_validation} Date or date string in default format, end of validation period
#' }
#' @section Input data definition:
#' For each data input folder, the following variables have to be defined and contain valid values as defined in \link[ittr]{hbv_pso}: prec, airt, ep, obs, area, elev_zones.
#' The input data can be assinged in the following ways (in descending order of precedence):
#' \enumerate{
#' \item assignment in localconf.R. Add code which assigns the input data to the corresponding variable name, e.g. \code{prec <- read.zoo(fpath, custom_options)}
#' \item text files inside input data directories. The files must be named as the input time series (with either .csv, .txt or no suffix) and produce a valid zoo object with index when read by \link[zoo]{read.zoo} with default arguments
#' \item import from HBV-light. Dropping a "Data" directory from HBV-light inside the input data directory will allow \code{hbv_pso_batch} to automatically import the data from the format used by HBV-light.
#' }
#' If a given input data set is not defined in localconf, the \code{hbv_pso_batch} will try to load it from a csv file, and then (if unsuccesful) from a HBV-light "Data" directory.
#'
#' @param basedir String, file path to the base directory containing glocal_config.R and directories with input data
#' @param ... Parameters set here will overwrite those set in glocal_config (but overwritten if also set in localconf)
#'
#' @return A list of the following items:
#' \enumerate{
#' \item \code{summary} data frame with id, gof, gov_valid, from, to, from_valid and to_valid for each optimization run
#' \item \code{results} list with hydro_pso output from all optimization runs, named by id
#' }
#' In addition, the output files from each optimization run are written inside the respective data input directory, in a sub-directory named by the identifier of the optimization run.
#'
#'
#'
#' @export
#'
#' @examples
hbv_pso_batch <- function(basedir,...) {
  start_time <- Sys.time()
  products <- list.dirs(basedir,recursive=FALSE, full.names=FALSE)
  all_runs <-lapply(products,hbv_pso_run,c(basedir=basedir,list(...)))
  # extract summary as dataframe and remove from results
  batch_summary <- as.data.frame(do.call(rbind,lapply(all_runs, `[[`, 3)))
  all_runs = lapply(all_runs,`[`,1:2)
  # flatten result list
  results <- unlist(all_runs,recursive = FALSE)
  run_time <- difftime(Sys.time(), start_time, units="secs")
  message("Execution time (hours:minutes:seconds): ",format(.POSIXct(run_time, tz="UTZ"), "%H:%M:%S"))
  return(list(summary=batch_summary,results=results))
}

#' Title
#'
#' @param product
#' @param basedir
#' @param ...
#'
#' @return
#' @keywords internal
#'
#' @examples
hbv_pso_run <- function(product,basedir,...){
  # TODO: read input data from .csv files in basedir
  # TODO: rename "prodcut" to something more generic
  # TODO: find a away to allow setting nested parameters through ...,
  #       e.g. maxit inside hydroGOF_args$control? list merges?
  # TODO: warn/stop if certain variabels exist in parent env, e.g. prec, airt etc
  # TODO: error handling if read.zoo fails/produces unexpected results?
  source(file.path(basedir,"global_config.R"),local=TRUE)
  # setting up paths and local config
  base_path <- file.path(basedir,product)
  localconf <- file.path(base_path,"localconf.R")
  if (file.exists(localconf)) {
    source(localconf,local=TRUE,chdir = TRUE)
  }
  list2env(list(...),env=environment())
  if (is.null(suffix)) {
    id <- paste(product,gof.name,sep="_")
  } else {
    id <- paste(product,gof.name,suffix,sep="_")
  }
  data_path <- file.path(base_path, "Data")
  output_path <- file.path(base_path,id)

  # set sim-obs plot titles and file names to id
  sim_obs_fname <- paste0(id,".png")
  if(isTRUE(plotting))
    plotting <- list(gof.name=gof.name, modelout.best.png.fname=sim_obs_fname,
                   main = id)

  # read in missing data from hbv-light and unpack them into the current environment
  var_names <- c("prec","airt","ep","obs","area","elev_zones")
  # missing_vars <- var_names[sapply(var_names,function(x) is.null(get0(x)))]
  missing_vars <- var_names[!sapply(var_names,exists,where=environment(),inherits=FALSE)]
  vars_from_csv <- sapply(missing_vars, function(x, base_path) {
    fp <- list.files(base_path, full.names = TRUE,
                     pattern = paste0("(?i)", x, "(\\.csv|\\.txt)?$"))[1]
    if (!is.na(fp)) {
      return(zoo::read.zoo(fp))
    }
    else {
      return(NULL)
    }
  }, base_path)
  # list2env(vars_from_csv, env = environment())
  # missing_vars <- lapply(var_names,function(x) is.null(get0(x)))
  list2env(Filter(Negate(is.null), vars_from_csv),env=environment())
  missing_vars <- sapply(var_names,exists,where=environment(),inherits=FALSE)
  missing_vars <- sapply(!missing_vars,list)
  if (dir.exists(data_path)) {
      ts_from_hbv_light <- do.call(parse_hbv_light,c(data_path,missing_vars))
      list2env(ts_from_hbv_light,envir=environment())
  }
  missing_vars <- var_names[!sapply(var_names,exists,where=environment(),inherits=FALSE)]
  if (length(missing_vars) > 0)
    stop("Could not find the input data for " ,id, ": ",paste(missing_vars))

  optimized <- hbv_pso(
    prec=prec,
    airt=airt,
    ep=ep,
    area=area,
    elev_zones=elev_zones,
    param=param,
    obs=obs,
    from=from,
    to=to,
    warmup=warmup,
    pelev = pelev,
    telev = telev,
    incon = incon,
    outpath = output_path,
    hydroPSO_args = hydroPSO_args,
    FUN_gof = FUN_gof,
    FUN_gof_args = FUN_gof_args,
    plotting = plotting
  )

  if (!is.null(from_validation)) {
    output_path <- file.path(output_path,"validation")
    sim_obs_fname <- file.path(output_path, paste0(id,"-validation",".png"))
    if(is.list(plotting))
      plotting <- list(png.fname=sim_obs_fname, main = id)

    optimized_validation = hbv_pso(
      prec=prec,
      airt=airt,
      ep=ep,
      area=area,
      elev_zones=elev_zones,
      param=optimized$pso_out$par[1:15],
      obs=obs,
      from=from,
      to=to_validation,
      warmup=from_validation,
      pelev = pelev,
      telev = telev,
      incon = incon,
      outpath = output_path,
      hydroPSO_args = list(control = control),
      FUN_gof = FUN_gof,
      FUN_gof_args = FUN_gof_args,
      plotting = plotting
    )
  }
  else {
    optimized_validation <- NULL
  }
  summary <- list(id = paste(product,suffix,sep="-"),
      gof_name = gof.name,
      gof = optimized$gof,
      gov_validation = ifelse(is.null(optimized_validation),NA,optimized_validation$gof),
      from = ifelse(is.null(from),NA, from),
      to = ifelse(is.null(to),NA, to),
      from_validation = ifelse(is.null(from_validation),NA, from_validation),
      to_validation = ifelse(is.null(to_validation),NA, to_validation)
  )
  return(setNames(list(optimized,optimized_validation,summary),c(id,paste0(id,"_valid"),"summary")))
}
