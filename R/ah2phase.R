
#' Fit Additive Hazards Regression Models to Two-phase Sampling
#'
#' The function fits a semiparametric additive hazards model \deqn{ \lambda(t|Z=z) = \lambda_0(t) + \beta'z.} to 
#' two-phase sampling data. The estimating procedures follow Hu (2014). 
#'
#' @param formula a formula object for the regression model of the form response ~ predictors.
#'        The outcome is a survival object created by \code{\link[survival]{Surv}}.
#' @param data  a data frame. Input dataset.
#' @param R  a phase II membership indicator. A vector of values of 0 and 1. The subject is selected to phase II if R = 1.
#' @param Pi  the  probability of a subject to be selected to the phase II subsample. 
#' @param robust a logical variable.  Robust standard errors are provided if robust = TRUE.
#' @param calibration.variables  a vector of some column names of the data. These are the  variables available for every observation. 
#'      They are used to calibrate the weight assigned to each subject in order to improve estimation efficiency. 
#' @return An object of class 'ah.2h' representing the fit.
#' @importFrom survival Surv
#' @export
#'
#' @note 
#' This function estimates both model-based and robust standard errors. It can be used to analyze
#' case-cohort studies. It allows subsampling among cases. It can incoporate the calibration procedure 
#' and analyze the combined dataset of  phase I and phase II samples.
#'  
#'
#' @seealso \code{\link{predict.ah.2ph}} for prediction based on fitted additive hazards model with two-phase sampling
#' and \code{\link{nwtsco}} for the description of nwtsco dataset.
#'
#' @references
#' Jie Hu (2014) A Z-estimation System for Two-phase Sampling with Applications to Additive Hazards Models and 
#' Epidemiologic Studies. Dissertation, University of Washington.
# 
#' @examples
#' library(survival)
#' library(rootSolve)
#' ### load data 
#' nwts <- nwtsco[1:100,]
#'
#' ### create strata based on institutional histology and disease status
#' nwts$strt <- 1+nwts$instit
#' ### add a stratum containing all (relapsed) cases
#' nwts$strt[nwts$relaps==1] <- 3
#'
#' ######################################################################
#' ### When phase II samples are obtained by Bernoulli Sampling #########
#' ######################################################################
#'
#' ### assign phase II subsampling (Bernoulli) probabilities 
#' ### oversample unfavorable histology (instit =1) and cases
#' ### Pi = 0.5 for instit =0, Pi = 0.9 for instit =1 and relaps =1 
#' nwts$Pi <- 0.5 * (nwts$strt == 1) + 0.9 * (nwts$strt == 2) +
#'            0.9 * (nwts$strt == 3)
#'
#' ### generate phase II sampling indicators
#' N <- dim(nwts)[1]
#' nwts$in.ph2 <-  rbinom(N, 1, nwts$Pi)
#'
#' ### fit an additive hazards model to two-phase sampling data without calibration
#' fit1 <- ah.2ph(Surv(trel,relaps) ~ age + histol, data = nwts, R = in.ph2, Pi = Pi,
#'                                  robust = FALSE,  calibration.variables = NULL)
#' summary(fit1)
#'
#' 
#' ### fit an additve hazards model with calibration on age 
#' fit2 <- ah.2ph(Surv(trel,relaps) ~ age + histol, data = nwts, R = in.ph2, Pi = Pi, 
#'                                    robust = FALSE, calibration.variables = "age")
#' summary(fit2)
#'
#' ### calibrate on age square
#' ### note if users create a  calibration variable, then 
#' ### the new variable should be added to the original data frame
#' nwts$age2 <- nwts$age^2 
#' fit3 <- ah.2ph(Surv(trel,relaps) ~ age + histol, data = nwts, R = in.ph2, Pi = Pi, 
#'                                    robust = FALSE, calibration.variables = "age2")
#' summary(fit3)
#' 
#' #############################################################################
#' ### When phase II samples are obtained by finite Sampling       #############
#' #############################################################################
#'
#' ### calculating the sample size for each straum 
#' ### calculate the strata size
#' strt.size <- table(nwts$strt)
#' ### the subsample sizes are determined by the following probabilities 
#' ### 0.5 for nwts$strt == 1 , 0.9 for nwts$strt == 2 and nwts$strt == 3
#' strt.prob <- c(0.5, 0.9, 0.9)
#' ### determine the subsample size for each stratum
#' size.ph2 <- round(strt.size * strt.prob)
#' ### Generate phase II selection indicators 
#' nwts$in.ph2 <- rep(1, nrow(nwts))
#' for (i in 1:length(strt.size)){
#'    idx.selected <- sample(1:strt.size[i], size.ph2[i], replace = F)
#'    nwts$in.ph2[nwts$strt ==i][-idx.selected] <- 0
#' }
#' ### fit an additve hazards model with finite stratified sampling 
#' ### calculate the sampling fractions 
#' frac <- size.ph2/strt.size
#' ### treating the problem as bernoulli sampling coupled with calibration on strata sizes
#' ### using frac as the sampling probilities
#' nwts$Pi <- NULL
#' for (i in 1:length(strt.size)){
#'   nwts$Pi[nwts$strt ==i] <- frac[i]
#' }
#'
#' ### create strt indicators, which become our calibration variables
#' for (i in 1:length(strt.size)){
#'    nwts$strt_ind <- as.numeric(nwts$strt ==i)
#'    names(nwts)[ncol(nwts)]= paste0("strt", i)
#' }
#' ### fit an additve hazards model with finate sampling
#' fit4 <- ah.2ph(Surv(trel,relaps) ~ age + histol, data = nwts, R = in.ph2, Pi = Pi,
#'                                  robust = FALSE,  calibration.variables = c("strt1","strt2","strt3"))

#' summary(fit4)


ah.2ph <- function(formula, data, R, Pi, robust = FALSE, calibration.variables = NULL) {
    
    # Z.pha2<-Z[R==1,] Y.pha2<-Y[R==1,] Creating ask
    Call <- match.call()
    R = data[, as.character(Call[["R"]])]
    Pi = data[, as.character(Call[["Pi"]])]
    calibration.variables = data[, calibration.variables]
    Pi.pha2 <- Pi[R == 1]
    data.pha2 <- data[R == 1, ]
 
    data.pha2$wts.pha2 <- as.numeric(1/Pi.pha2)
    
    
    
    
    
    
    if (!length(calibration.variables)) {
        
        # Use the new weights and fit the model to the data
        fit.A <- ah(formula, data = data.pha2, robust = robust, 
            weights = wts.pha2)
        resid <- fit.A$resid
        
        temp <- resid * sqrt(1 - Pi.pha2)
        temp1 <- resid * sqrt(Pi.pha2)
        
        
    } else {
        
        aux <- as.matrix(calibration.variables)
        P <- t(aux) %*% (aux)
        aux.pha2 <- aux[R == 1, ]
        
        wts.cal <- cook.wts.cal(aux = aux, aux.pha2 = aux.pha2, 
            P = P, wts.pha2 = data.pha2$wts.pha2)
        data.pha2$wts.pha2 <- wts.cal
        
        fit.A <- ah(formula, data = data.pha2, robust = robust, 
            weights = wts.pha2)
        resid <- fit.A$resid
        temp1 <- resid * sqrt(1/wts.cal)
        
        
        Q <- t(aux.pha2 * sqrt(wts.cal)) %*% (resid/sqrt(wts.cal))  # multiplied by sqrt(wts.cal) because Qf= sum wts.cal*f
        # and resid already weighted by wts.cal
        resid.adj <- resid - (aux.pha2 %*% solve(P) %*% Q) * 
            wts.cal
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
    fit$wts.pha2 <- data.pha2$wts.pha2
    fit$calibration.variables <- calibration.variables
    fit$R <- R
    fit$call <- Call
    fit$fit.pha1 <- fit.A
    
    
    class(fit) <- "ah.2ph"
    fit
    
}




################################################################################################# calculate the new weight
################################################################################################# #####################################

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
    
    eval <- multiroot(model.calibration, start = L0)
    L <- eval$root
    # est.acc<-eval$est.acc
    wts.cal <- as.vector(exp(-aux.pha2 %*% L) * wts.pha2)
    return(wts.cal)
}




 
