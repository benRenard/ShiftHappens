#' Linear regression
#'
#' Fit a linear regression (with possibly several predictors but a single predictand)
#'
#' @param x real matrix or data frame, predictors
#' @param y real vector, predictand
#' @param time vector (numeric or date), time
#' @inherit fittedModel return
#' @export
#' @importFrom stats lm sd coef
#' @examples
#' f=Fit_LinearRegression(x=RhoneRiverAMAX$H,y=RhoneRiverAMAX$Q,time=RhoneRiverAMAX$Date)
#' plot(f$data$x1,f$data$obs);lines(f$data$x1,f$data$sim,col='red')
#' plot(f$data$time,f$data$res)
Fit_LinearRegression <- function(x,y,time=1:NROW(y)){
  if(NROW(y)<=2){
    warning('NULL was returned because it not possible to perform linear regression with fewer than two points.')
    return(NULL)
  }
  # Linear regression
  mod <- stats::lm(as.matrix(y)~as.matrix(x))
  # Assemble results
  DF=as.data.frame(x)
  if(is.null(names(x))){names(DF) <- paste0('x',1:NCOL(DF))}
  DF=cbind(time=as.data.frame(time)[,1],DF,obs=as.matrix(y)[,1],sim=mod$fitted.values,
           res=mod$residuals,uRes=stats::sd(mod$residuals))
  pars=stats::coef(mod)
  out=fittedModel(data=DF,parameters=stats::coef(mod))
  return(out)
}
