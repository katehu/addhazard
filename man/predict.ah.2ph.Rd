% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pred.2ph.R
\name{predict.ah.2ph}
\alias{predict.ah.2ph}
\title{Prediction Based on the Additive Hazards Model Fitted from Two-phase Sampling}
\usage{
\method{predict}{ah.2ph}(object, newdata, newtime, ...)
}
\arguments{
\item{object}{an object of class inhering from 'ah.2ph'.}

\item{newdata}{a dataframe of an individual's predictors.}

\item{newtime}{a given sequence of time points at which the prediction is performed.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
A dataframe including the given time points, predicted hazards, their
 standard errors, their variances, the phase I component of the variance for
 predicted hazards and the phase II component of the variance.
}
\description{
This function predicts a subject's overall hazard rates at given time points
based on this subject's covariate values. The prediction function is an object
from \code{\link{ah.2ph}}. The  estimating procedures follow Hu (2014).
}
\examples{
library(survival)
### fit an additive hazards model to  two-phase sampling data without calibration
fit1 <- ah.2ph(Surv(trel,relaps) ~ age + histol, data = nwts2ph,  R = in.ph2,  
               Pi = Pi, robust = FALSE, ties = 'break')

###  input the new data for prediction
newdata <- nwtsco[101,]
###  based on the fitted model fit1, perform prediction at time points t =3 and t= 5
predict(fit1, newdata, newtime = c(3,5))

### fit an additve hazards model to  two-phase sampling data with calibration
### The calibration variable is instit
fit2 <- ah.2ph(Surv(trel,relaps) ~ age + histol, data = nwts2ph, R = in.ph2, Pi = Pi,
                                   ties = 'break', robust = FALSE, calibration.variables = "instit")
### based on the fitted model fit2, perform prediction at time points t =3 and t= 5
predict(fit2, newdata, newtime = c(3,5))
}
\references{
Jie Hu (2014) A Z-estimation System for Two-phase Sampling with Applications
              to Additive Hazards Models and Epidemiologic Studies. Dissertation,
              University of Washington.
}
\seealso{
\code{\link{ah.2ph}} for fitting the additive hazards model with two-phase
}
