
#' (Batch-) Run Particle Swarm Optimizations of the HBV Model
#'
#' Run Particle Swarm Optimizations of the HBV Model using [`hydroPSO`]
#' and [`TUWmodel`] by importing configurations from R files. Supports
#' specifiying multiple configuration files to run optimizations in batches.
#'
#' @details Runs Particle Swarm Optimizations of HBV for given parameter sets
#'   defined in simple R scripts. Timeseries input data (prec, airt, ep, obs,
#'   area, elevation_zones) can be automatically imported from csv or HBV-light
#'   format.
#'
#'   Each optimization run is given an identification string concatenated from
#'   the name of the input data directory (i.e. the location of the
#'   configuration file), the name of the gof variable and an optional suffix.
#'   This identifier is used as directory name for the output files, in plot
#'   titles and file names, and in the result list returned by `hbv_pso_run`.
#'
#'
#' @section Parameter definition: Parameters for the different optimization runs
#'   are defined by simple variable assignment, where the variable name
#'   corresponds to the parameter as used in [hbv_pso()]. They are defined in
#'   simple R files passed to `hbv_pso_run` (see `configpath`). In order for the
#'   configuration files to be picked up when passing directory paths to
#'   `hbv_pso_run`, the configuration file names must start with "run".
#'
#'   Parameters can also be set by passing them as additional arguments to
#'   `hbv_pso_run`, which will overwrite corresponding values defined in
#'   configuration file(s).
#'   Be aware that when using this option, any lists passed are not merged but
#'   completely replace the values set in the configfile.
#'
#'   In addition to parameters described in [hbv_pso()], the following optional
#'   parameters unique to `hbv_pso_run` can also be set in configuration files:
#'   `suffix`, `gof.name`,`from_validation`, `to_validation` (see above).
#'
#' @section Input data definition: For each optimization run, the following
#'   variables have to be defined and contain valid values as defined in
#'   [hbv_pso()]: prec, airt, ep, obs, area, elev_zones. The input data
#'   can be assinged in the following ways (in descending order of precedence):
#'
#'   1. *assignment inside the configfile*: Add code which assigns the input
#`      data to the corresponding variable name, e.g. `prec <-
#`      read.zoo(fpath, custom_options)
#'   2. *text files inside input data directories*:
#'       The files must be named as the input time series (with either
#'       .csv, .txt or no suffix, e.g. `prec.csv`) and be located in the same
#'       directory as the configfile. They must produce a valid zoo object with
#'       index when read by [read.zoo] with default arguments
#'   3. *HBV-light input files*:
#'       Analogue to the text files above, input data can also be imported from
#'       files as used by HBV-light ("ptq.text" and "ep.txt")
#'
#' If a given input data set is not defined in the configuration file,
#' `hbv_pso_run` will try to load it from a corresponding file in the input data
#' directory - first from a csv file, and then from the corrseponding HBV-light
#' file (e.g. PTQ.txt for prec, airt or obs).
#'
#' @param configpath character vector of full path and/or file names. Any file
#'   given is sourced and an `hbv_pso` optimization run is started with the
#'   imported arguments. Any directory given is searched for R files starting
#'   with "run", and an optimization run is started for each found in the same
#'   way as for file paths.
#' @param recursive logical, should config files be searched in subdirectories?
#'   If true, `hbv_pso_run` will recurse once into any directory given in
#'   `configpath` before searching for configuration files

#' @param suffix String, will be appended to the optimization run identifier
#' @param gof.name String, Used as part of the optimization run identifier.
#' Shorthand for setting plotting$gof.name. If gof.name is not set directly,
#' the value from plotting$gof.name is used.
#' @param from_validation Date or date string in default format, start of
#' validation run with best parameter set found during optimization
#' @param to_validation Date or date string in default format, end of validation
#' period
#' @param ... Parameters passed to hbv_pso, will overwrite those set in in the
#'   configuration files.
#'
#' @return A list of the following items:
#'
#' 1. `summary` data
#'   frame with id, gof, gof_valid, from, to, from_valid and to_valid for each
#'   optimization run
#' 2. `results` list with hydro_pso output from all
#'   optimization runs, named by id
#'
#' In addition, the output files from each optimization run are written inside a
#' directory which is located in the same directory as the configuration file,
#' and named by the identifier of the respective optimization run.
#'
#' @md
#' @export
#'
#' @examples
hbv_pso_run <-
  function(configpath, recursive = FALSE, suffix = NULL, gof.name = NULL,
           from_validation = NULL, to_validation = NULL, ...) {
  # TODO: abstract common tasks (e.g. ts loading) into util function
  # TODO: plotting list<->TRUE conversion done in multiple places?
  start_time <- Sys.time()
  isdir <- sapply(configpath,function(x) file.info(x)$isdir)
  dirs <- configpath[isdir]
  if (recursive && length(dirs)>0) {
    dirs <- list.dirs(dirs,recursive=FALSE, full.names=TRUE)
  }
  configs <- list.files(dirs, full.names = TRUE,
                        pattern = "(?i)^run.*\\.R$")
  all_runs <-lapply(c(configs,configpath[!isdir]),hbv_pso_run_single,...)
  # build summary dataframe and remove from results
  batch_summary <- do.call(rbind,lapply(all_runs, `[[`, "summary"))
  batch_summary_par<- do.call(rbind,lapply(all_runs, `[[`, "summary_par"))
  all_runs = lapply(all_runs,`[[`,"results")
  # flatten result list
  results <- unlist(all_runs,recursive = FALSE)
  run_time <- difftime(Sys.time(), start_time, units="secs")
  print(batch_summary)
  message("Execution time (hours:minutes:seconds): ",format(.POSIXct(run_time, tz="UTZ"), "%H:%M:%S"))
  return(list(summary=batch_summary,summary_par = batch_summary_par, results=results))
}

#' Title
#'
#' @param configfile
#' @param ...
#'
#' @return
#' @keywords internal
#'
#' @examples
hbv_pso_run_single <- function(configfile,...){
  # TODO: find a away to allow setting nested parameters through ...,
  #       e.g. maxit inside hydroGOF_args$control? list merges?
  # TODO: error handling if read.zoo fails/produces unexpected results?
  # TODO: test everything with non-ts inputs

  print(paste("Running",configfile))
  config_env <- new.env()
  source(configfile, local = config_env, chdir = TRUE)
  list2env(list(...), envir = config_env)

  # easier handling of plotting arg: list if we should plot, FALSE otherwise
  plotting <- config_env$plotting
  if (is.null(plotting)) {
    plotting <- FALSE
  } else if (isTRUE(plotting)) {
    plotting <- list()
  }

  # ensure objective function defaults are correct if FUN_gof is NULL
  if (is.null(config_env$FUN_gof)) {
    config_env$FUN_gof <- hydroGOF::NSE
    config_env$FUN_gof_args <- NULL
    config_env$gof.name <- "NSE"
  }

  # if defined, gof.name sets plotting$gof.name. if not, get it from plotting
  gof.name <- config_env$gof.name
  if (is.list(plotting)) {
    if (is.null(gof.name)) {
      gof.name <- plotting$gof.name
    } else {
      plotting$gof.name <- gof.name
    }
  }

  suffix <- config_env$suffix

  config_dir <- dirname(configfile)
  id <- paste(c(basename(config_dir), gof.name, suffix), collapse = "_")
  if (is.null(config_env$outpath)) {
    output_path <- file.path(config_dir, "hbv", id)
  } else {
    output_path <- file.path(config_env$outpath, "hbv", id)
  }


  # set sim-obs plot titles and file names to id
  sim_obs_fname <- paste0(id,".png")
  if(is.list(plotting)) {
    plot_args <- list(modelout.best.png.fname=sim_obs_fname, main = id)
    if(is.list(plotting)) {
      plotting <- c(plotting,plot_args[!(names(plot_args) %in% names(plotting))])
    } else {
      plotting <- plot_args
    }
  }

  var_names <- c("prec","airt","ep","obs","area","elev_zones")
  # check for missing time series
  missing_vars <- var_names[sapply(var_names, function(x) is.null(config_env[[x]]))]
  if (length(missing_vars) > 0) {
    # read in missing data from csv files
    vars_from_csv <- sapply(missing_vars, function(x, config_dir) {
      fp <- list.files(config_dir, full.names = TRUE,
                       pattern = paste0("(?i)", x, "(\\.csv|\\.txt)?$"))[1]
      if (!is.na(fp)) {
        return(zoo::read.zoo(fp))
      }
      else {
        return(NULL)
      }
    }, simplify=FALSE, USE.NAMES=TRUE, config_dir)

    list2env(vars_from_csv,envir=config_env)
    # read in missing data from hbv-light files
    missing_vars <- sapply(var_names, function(x) is.null(config_env[[x]]))
    if (any(missing_vars)) {
        ts_from_hbv_light <- do.call(parse_hbv_light,c(hbv_light_dir=config_dir,as.list(missing_vars)))
        list2env(ts_from_hbv_light,envir=config_env)
    }
    missing_vars <- var_names[sapply(var_names, function(x) is.null(config_env[[x]]))]
    #TODO: area and elev_zone are optional!
    if (length(missing_vars) > 0)
      stop("Could not find the following input data for " ,id, ": ",paste(missing_vars,collapse=","))
  }


  # collecting arguments for hbv_pso
  hbv_pso_args <- list(
    prec = config_env$prec,
    airt = config_env$airt,
    ep = config_env$ep,
    area = config_env$area,
    elev_zones = config_env$elev_zones,
    param = config_env$param,
    obs = config_env$obs,
    from = config_env$from,
    to = config_env$to,
    warmup = config_env$warmup,
    pelev = config_env$pelev,
    telev = config_env$telev,
    incon = config_env$incon,
    outpath = output_path,
    hydroPSO_args = config_env$hydroPSO_args,
    FUN_gof = config_env$FUN_gof,
    FUN_gof_args = config_env$FUN_gof_args,
    plotting = plotting)

  # remove arguments with value NULL
  hbv_pso_args <- hbv_pso_args[!sapply(hbv_pso_args, is.null)]

  optimized <- do.call(hbv_pso,hbv_pso_args)

  if (!is.null(config_env$from_validation)) {
    hbv_pso_args$outpath <- file.path(output_path,"validation")
    hbv_pso_args$param <- optimized$pso_out$par
    hbv_pso_args$to <- config_env$to_validation
    hbv_pso_args$warmup <- config_env$from_validation

    if(is.list(plotting)) {
      sim_obs_fname <- file.path(output_path, paste0(id,"_ggof-validation",".png"))
      hbv_pso_args$plotting <- list(png.fname=sim_obs_fname, main = paste0(id,"-validation"))
    }

    optimized_validation <- do.call(hbv_pso,hbv_pso_args)
  } else {
    optimized_validation <- NULL
  }

  summary_base <- data.frame(
    id = paste(c(basename(config_dir),suffix),collapse = "-"),
    gof_name = ifelse(is.null(gof.name),NA, gof.name),
    gof = round(optimized$gof,3),
    stringsAsFactors = FALSE
  )
  summary_gof <- data.frame(
    gof_validation = ifelse(
      is.null(optimized_validation),NA,round(optimized_validation$gof,3)),
    from = ifelse(is.null(config_env$warmup),config_env$from, config_env$warmup),
    to = ifelse(is.null(config_env$to),NA, config_env$to),
    from_validation = ifelse(
      is.null(config_env$from_validation),NA, config_env$from_validation),
    to_validation = ifelse(is.null(config_env$to_validation),NA, config_env$to_validation),
    stringsAsFactors = FALSE
  )

  res <- list(results=list(),
              summary=cbind(summary_base,summary_gof),
              summary_par = cbind(summary_base,t(optimized$pso_out$par)))
  res$results[[id]] <- optimized
  res$results[[paste0(id,"_validation")]] <- optimized_validation
  return(res)
}
