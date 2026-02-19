#' Linear regression
#'
#' Fit a linear regression (with possibly several predictors but a single predictand)
#'
#' @param x real matrix or data frame, predictors
#' @param y real vector, predictand
#' @param time vector (numeric or date), time
#' @param uX real matrix or data frame, uncertainty in predictors (ignored)
#' @param uY real vector, uncertainty in predictand (ignored)
#' @inherit fittedModel return
#' @export
#' @importFrom stats lm sd coef
#' @examples
#' f=Fit_LinearRegression(x=RhoneRiverAMAX$H,y=RhoneRiverAMAX$Q,time=RhoneRiverAMAX$Date)
#' plot(f$data$x1,f$data$y);lines(f$data$x1,f$data$ysim,col='red')
#' plot(f$data$time,f$data$res)
Fit_LinearRegression <- function(x,y,time=1:NROW(y),uX=0*x,uY=0*y){
  if(NROW(y)<=2){
    warning('NULL was returned because it not possible to perform linear regression with fewer than two points.')
    return(NULL)
  }
  # Linear regression
  mod <- stats::lm(as.matrix(y)~as.matrix(x))
  # Assemble results
  uRes=stats::sd(mod$residuals)
  vYsim=uRes^2-uY^2;vYsim[vYsim<0]=0
  pars=stats::coef(mod)
  out=fittedModel(time=time,x=x,y=y,ysim=mod$fitted.values,res=mod$residuals,
                  uX=uX,uY=uY,uYsim=sqrt(vYsim),uRes=rep(uRes,length(y)),
                  parameters=as.numeric(pars))
  return(out)
}
