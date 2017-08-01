#' Disaggregate Daily Means by Month-of-the-year into Daily Timeseries
#'
#' Wrapper around na.approx/na.spline which disaggregates monthly to daily values
#' by linear interpolation or cubic spline interpolation, respectively.
#'
#' The values for interpolation are mean daily data grouped by month (i.e. 12 values corresponding to
#' the calendar month, as produced by hydroTSM::monthlyfunction with FUN=mean).
#'
#' Monthly input data is fixed to the middle of each month for the interpolation, and values at
#' the left and right margins are interpolated by extending the monthly datapoints accordingly.
#' @param from Beginning of the interpolated time series, as Date object
#' @param to End of the interpolated time series, as Date object
#' @param values Vector of daily mean values for each month (length=12)
#' @param FUN interpolation function to use, must be either na.approx (linear) or na.spline (spline)
#' @param ... further arguments passed to approx/spline via na.approx/na.spline
#'
#' @return A daily timeseries (zoo) interpolated from daily mean values for each month of the year
#' @import zoo
#' @author Joschka Thurner, \email{joschka.thurner@th-koeln.de}
#' @examples
#' m2d <- monthly2daily(as.Date("2000-01-01"),as.Date("2001-12-31"),values=c(0:6,5:1),FUN=zoo::na.approx)
#' tail(m2d)
#' \dontrun{
#' # using na.approx directly for comparison
#' m <- zooreg(rep(c(0:6,5:1),2),as.yearmon("2000-01-01"),freq = 12)
#' sq <- seq(as.Date(start(m)), as.Date(end(m), frac = 1), by = "day")
#' zd <- na.approx(m, x = as.Date, xout = sq)
#'
#' # plot both results
#' # na.approx is missing values at the right margin, and sets the monthly value
#' # at the first day of the month
#' plot(m2d,col="green")
#' lines(zd,col="red")
#' legend('topright', c("na.approx","monthly2daily"), col = 2:3, lty = 1)
#' }
monthly2daily <-
  function(from, to, values, FUN = zoo::na.spline, ...) {
    #TODO: input checking - stop if from/to aren't dates? try to convert from string?
    if (!any(sapply(c(zoo::na.approx, zoo::na.spline), identical, FUN))) {
      stop("FUN must be either na.approx or na.spline")
    }
    # add one year on both ends to ensure from-to is always fully inside the data points
    start <- zoo::as.yearmon(from) - 1
    n_months <- round((zoo::as.yearmon(to) - start) * 12 + 13)
    # reorder monthly values if start != January
    firstmonth <- round(as.numeric(start %% 1 * 12))
    if (firstmonth > 0)
      values <- c(values[-c(1:firstmonth)], c(values[1:firstmonth]))
    # create monthly timeseries with values at middle of the month
    values <- rep(values, length.out = n_months)
    monthly = zoo::zooreg(values, start, freq = 12)
    index(monthly) <- zoo::as.Date(index(monthly), frac = 0.5)
    # interpolate to daily values
    idx <- seq(from, to, frac = 1, by = "day")
    daily <- FUN(monthly, x = as.Date, xout = idx, ...)
    return(round(daily,4))
  }
#
# library(xts)
# library(tempdisagg)
#
# #############
# # comparing na.spline, monthly2daily, tempdisagg
# from = as.Date("2000-01-01")
# to = as.Date("2001-12-31")
# values <- rep(c(0:6, 5:1), 2)
# # values <- rep(rnorm(12),2)
# # values <- rep(evap[,1],2)
# m <- zooreg(values, as.yearmon(from), freq = 12)
# idx_daily <- seq(from, to, frac = 1, by = "day")
#
# ################
# # na.spline
# da_spline <- na.spline(m, x = as.Date, xout = idx_daily)
#
# ################
# # tempdisagg
# m_d <- m
# index(m_d) <- zoo::as.Date(index(m_d))
# da_td <-
#   predict(
#     td(
#       as.xts(m_d) ~ 0,
#       lf.end = to,
#       hf = idx_daily,
#       conversion = "average",
#       method = "denton-cholette"
#     )
#   )
#
# ################
# # monthly2daily
# da_m2d <- monthly2daily(from, to, values, FUN = na.spline)
#
# ################
# # plot
# da <- as.xts(merge(da_spline, da_td, da_m2d))
# plot.xts(da,
#          screens = 1,
#          legend.loc = "top",
#          main = "Comparison of disaggregation functions - Values 0-1-2-3-4-5-6-5-4-3-2-1 for Jan-Dez")
