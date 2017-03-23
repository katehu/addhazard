#######################################################
# Author: Kate HU
# Date: July 27th, 2016
# Task: update finite sampling documentation
# Date: March 30th, 2017
#       clear problems in calling wts.pha2 and wts.cal
# To do: Interval censoring and add finite sampling argument
#        in the function
#######################################################
#' Fit Additive Hazards Regression Models to Two-phase Sampling
#'
#' The function fits a semiparametric additive hazards model
#'  \deqn{ \lambda(t|Z=z) = \lambda_0(t) + \beta'z.} to two-phase sampling data.
#'  The estimating procedures follow Hu (2014).
#'
#' @param formula a formula object for the regression model of the form
#'  response ~ predictors. The outcome is a survival object created by \code{\link[survival]{Surv}}.
#' @param data  a data frame. Input dataset.
#' @param ties a logical variable. FALSE if there are no ties in the censored failure times.
#' @param R  a phase II membership indicator. A vector of values of 0 and 1.
#'  The subject is selected to phase II if R = 1.
#' @param Pi  the  probability of a subject to be selected to the phase II subsample.
#' @param robust a logical variable.  Robust standard errors are provided if robust = TRUE.
#' @param calibration.variables  a vector of some column names of the data.
#'  These are the  variables available for every observation. They are used to
#'  calibrate the weight assigned to each subject in order to improve estimation efficiency.
#' @param ...	additional arguments to be passed to the low level regression fitting
#'  functions.
#' @return An object of class 'ah.2h' representing the fit.
#' @importFrom survival Surv
#' @export
#'
#' @note
#' This function estimates both model-based and robust standard errors. It can be
#' used to analyze case-cohort studies. It allows subsampling among cases. It can
#' incorporate the calibration procedure and analyze the combined dataset of phase
#' I and phase II samples.
#'
#' @seealso \code{\link{predict.ah.2ph}} for prediction based on fitted additive
#' hazards model with two-phase sampling and \code{\link{nwtsco}} for the description
#' of nwtsco dataset.
#'
#' @references
#' Jie Hu (2014) A Z-estimation System for Two-phase Sampling with Applications to
#' Additive Hazards Models and  Epidemiologic Studies. Dissertation, University of Washington.
#
#' @examples
#' library(survival)
#' ### fit an additive hazards model to two-phase sampling data without calibration
#' nwts2ph$trel <- nwts2ph$trel + runif(dim(nwts2ph)[1],0,1)*1e-8
#' fit1 <- ah.2ph(Surv(trel,relaps) ~ age + histol, ties = FALSE, data = nwts2ph, R = in.ph2, Pi = Pi,
#'  robust = FALSE,  calibration.variables = NULL)
#' summary(fit1)
#'
#' ### fit an additve hazards model with calibration on age
#' fit2 <- ah.2ph(Surv(trel,relaps) ~ age + histol, ties = FALSE, data = nwts2ph, R = in.ph2, Pi = Pi,
#'  robust = FALSE, calibration.variables = "age")
#' summary(fit2)
#'
#' ### calibrate on age square
#' ### note if users create a  calibration variable, then
#' ### the new variable should be added to the original data frame
#' nwts2ph$age2 <- nwts2ph$age^2
#' fit3 <- ah.2ph(Surv(trel,relaps) ~ age + histol, ties = FALSE, data = nwts2ph, R = in.ph2, Pi = Pi,
#'                                    robust = FALSE, calibration.variables = "age2")
#' summary(fit3)
#'
#' #############################################################################
#' ### When phase II samples are obtained by finite Sampling       #############
#' #############################################################################
#'
#' ### calculating the sample size for each straum
#' ### calculate the strata size
#' strt.size <- table(nwts2ph$strt)
#' ph2.strt.size <- table(subset(nwts2ph, in.ph2 == 1)$strt)
#' ### fit an additve hazards model with finite stratified sampling
#' ### calculate the sampling fractions
#' frac <- ph2.strt.size/strt.size
#' ### treating the problem as bernoulli sampling coupled with calibration on strata sizes
#' ### using frac as the sampling probilities
#' nwts2ph_by_FPSS <- nwts2ph
#' nwts2ph_by_FPSS$Pi <- NULL
#' for (i in 1:length(strt.size)){
#'   nwts2ph_by_FPSS$Pi[nwts2ph_by_FPSS$strt ==i] <- frac[i]
#' }
#'
#' ### create strt indicators, which become our calibration variables
#' for (i in 1:length(strt.size)){
#'    nwts2ph_by_FPSS$strt_ind <- as.numeric(nwts2ph_by_FPSS$strt ==i)
#'    names(nwts2ph_by_FPSS)[ncol(nwts2ph_by_FPSS)]= paste0("strt", i)
#' }
#' ### fit an additve hazards model with finate sampling
#' fit4 <- ah.2ph(Surv(trel,relaps) ~ age + histol,
#'                                    data = nwts2ph_by_FPSS, ties = FALSE,
#'                                    R = in.ph2, Pi = Pi,
#'                                    robust = FALSE,
#'                                    calibration.variables = c("strt1","strt2","strt3"))
#' summary(fit4)

ah.2ph <- function(formula, data, R, Pi, ties, robust = FALSE,
                   calibration.variables = NULL, ...) {

    # Z.pha2<-Z[R==1,] Y.pha2<-Y[R==1,] Creating ask
    Call <- match.call()
    R = data[, as.character(Call[["R"]])]
    Pi = data[, as.character(Call[["Pi"]])]
    calibration.variables = data[, calibration.variables]
    Pi.pha2 <- Pi[R == 1]
    data.pha2 <- data[R == 1, ]
    wts.pha2 <- as.numeric(1/Pi.pha2)
    data.pha2$wts.pha2 <- wts.pha2

    if (!length(calibration.variables)) {

        # Use the new weights and fit the model to the data
        # In ahaz, weights is extracted from data by calling the column name
        # Thus the varible name assigned to weights  has to to included in
        # the column name of data
        fit.A <- ah(formula, data = data.pha2, robust = robust,
            weights =  wts.pha2,  ties = ties)
        resid <- fit.A$resid
        temp <- resid * sqrt(1 - Pi.pha2)
        temp1 <- resid * sqrt(Pi.pha2)
    } else {
        aux <- as.matrix(calibration.variables)
        P <- t(aux) %*% (aux)
        aux.pha2 <- aux[R == 1, ]
        wts.cal <- cook.wts.cal(aux = aux, aux.pha2 = aux.pha2,
            P = P, wts.pha2 = wts.pha2)
        data.pha2$wts.cal <- wts.cal
        fit.A <- ah(formula, data = data.pha2, robust = robust,
            weights = wts.cal, ties = ties)
        resid <- fit.A$resid
        temp1 <- resid * sqrt(1/wts.cal)
        # multiplied by sqrt(wts.cal) because Qf= sum wts.cal*f
        Q <- t(aux.pha2 * sqrt(wts.cal)) %*% (resid/sqrt(wts.cal))
        # and resid already weighted by wts.cal
        resid.adj <- resid - (aux.pha2 %*% solve(P) %*% Q) * wts.cal
        temp <- resid.adj * sqrt((1 - Pi.pha2)/(Pi.pha2 * wts.cal))

    }

    var.pha1 <- fit.A$var
    iA <- fit.A$iA
    if (robust == TRUE) {
        var.pha1 <- iA %*% t(temp1) %*% temp1 %*% iA
    }

    var.pha2 <- iA %*% t(temp) %*% temp %*% iA
    var.tot <- var.pha1 + var.pha2

    fit <- NULL
    fit$coef <- fit.A$coef
    fit$var.pha1 <- var.pha1
    fit$var.pha2 <- var.pha2
    fit$var.tot <- var.pha1 + var.pha2
    fit$se <- sqrt(diag(var.tot))
    fit$Pi.pha2 <- Pi.pha2
    fit$wts.pha2 <- wts.pha2
    fit$calibration.variables <- calibration.variables
    fit$R <- R
    fit$call <- Call
    fit$fit.pha1 <- fit.A
    class(fit) <- "ah.2ph"
    fit
}




############# calculate the new weight #####################################

cook.wts.cal <- function(aux, aux.pha2, P, wts.pha2) {



    if (!is.matrix(aux.pha2))
        aux.pha2 <- as.matrix(aux.pha2)
    aux.tot <- apply(aux, 2, sum)
    aux.tot.pha2 <- apply(aux.pha2 * wts.pha2, 2, sum)
    ### phase I total, 1 x q

    L0 <- solve(P) %*% (aux.tot.pha2 - aux.tot)





    model.calibration <- function(L) {
        F <- NULL
        wts.fish <- as.vector(exp(-aux.pha2 %*% L) * wts.pha2)
        for (i in 1:dim(aux)[2]) {
            F[i] <- sum(wts.fish * aux.pha2[, i]) - aux.tot[i]
        }
        F
    }

    eval <- rootSolve::multiroot(model.calibration, start = L0)
    L <- eval$root
    # est.acc<-eval$est.acc
    wts.cal <- as.vector(exp(-aux.pha2 %*% L) * wts.pha2)
    return(wts.cal)
}
