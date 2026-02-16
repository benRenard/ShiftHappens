#' Floods of the Rhone River at Beaucaire, France
#'
#' A data frame containing annual maximum stages (H, in meters) and discharges
#' (Q, in cubic meters per second) for the Rhone River at Beaucaire, France,
#' along with the associated uncertainties expressed as standard deviations (uH and uQ).
#' Details on the reconstruction of these long series can be found in the article
#' by Lucas et al. (2023) referenced below.
#' Note that years 1968, 1969 and 1970 are missing and are not included in the data frame.
#'
#' @format
#' \describe{
#'  \item{Year}{Year}
#'  \item{Date}{Date}
#'  \item{Time}{Time}
#'  \item{H}{Stage record}
#'  \item{uH}{Uncertainty of the stage expressed as standard deviation}
#'  \item{Q}{Discharge}
#'  \item{uQ}{Uncertainty of the discharge expressed as standard deviation}
#' }
#'
#' @source \doi{10.1016/j.jhydrol.2023.129840}
"RhoneRiverAMAX"

#' Gauging of the Ardèche River at Meyras, France, provided by UHPC Grand Delta
#'
#' A data frame containing stages (H, in meters) and discharges measurements
#' (Q, in cubic meters per second) for the Ardèche River at Meyras, France,
#' along with the associated uncertainties expressed as standard deviations (uQ) from 2001 until 2018
#'
#' @format
#' \describe{
#'   \item{Day}{Day}
#'   \item{Month}{Month}
#'   \item{Year}{Year}
#'   \item{Hour}{Hour}
#'   \item{Minute}{Minute}
#'   \item{Second}{Second}
#'   \item{Date}{Date}
#'   \item{H}{Stage measurement}
#'   \item{Q}{Discharge measurement}
#'   \item{uQ}{Discharge uncertainty expressed as a standard deviation}
#' }
#'
#' @source \doi{10.1029/2018WR023389}
"ArdecheRiverGaugings"

#' Stage record of the Ardèche River at Meyras, France, provided by UHPC Grand Delta
#'
#' A data frame containing stages (H, in meters) for the Ardèche River
#' at Meyras, France, from 07/11/2001 until 29/10/2018
#'
#' @format
#' \describe{
#'   \item{Date}{Date, in POSIXct format}
#'   \item{H}{Stage record}
#' }
#'
#' @source\url{https://theses.fr/2021GRALU006}
"ArdecheRiverStage"
