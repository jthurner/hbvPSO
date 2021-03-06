#' Default parameter ranges for TUWmodel
#'
#' Input parameter ranges as given in the \link[TUWmodel]{TUWmodel} documentation for \code{param}.
#' Additionally includes \code{tcalt} and \code{pcalt} as used in \link[hbvPSO]{hbv_pso}, set to 0 (no pcalt/tcalt applied).
#'
#' @format A matrix with two columns, min and max, and the following parameters as rows:
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
#' }
"tuwmodel_params_default"

