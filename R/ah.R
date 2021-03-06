###############################################################################
# Author: Kate HU
# Date: July 27th, 2017
#       update documentation for breaking ties
# Date: March 20th, 2017
#       importFrom stats model.weight
# Date: July 1st, 2020
#       we break ties for users
# To do: Tie handling
###############################################################################
#' Fit Additive Hazards Regression Models
#'
#' Fit a semiparametric additive hazard model
#' '\deqn{ \lambda(t|Z=z) = \lambda_0(t) + \beta'z.}
#' The estimating procedures follow Lin & Ying (1994).
#'
#' @param formula a formula object for the regression model of the form
#'  response ~ predictors. The outcome is a survival object created by
#'  \code{\link[survival]{Surv}}.
#' @param data  a data frame. Input dataset.
#' @param robust a logical variable.  Robust standard errors are provided if
#'  robust == TRUE.
#' @param weights  a numeric vector. The weight of each observation.
#' @param ties a string. If there are ties in the survival time, when ties = 'break'
#'        a small random number is added to the survival time to break the ties.
#' @param seed an integer. Seed number used to generate random increment when 
#'        breaking ties. The default number is 20. 
#' @param ... additional arguments to be passed to the low level regression
#'  fitting functions.
#'
#' @return An object of class 'ah' representing the fit.
#'
#' @note
#'  The response variable is a survival object. The regression
#'  model can be univariate or multivariate. This function is built upon the function
#'  \code{\link[ahaz]{ahaz}} by Anders Gorst-Rasmussen.
#'
#'
#' @seealso \code{\link{predict.ah}} for prediction based on fitted
#'  \code{\link{ah}} model, \code{\link{nwtsco}} for the description of nwtsco dataset
#'
#' @references
#' Lin, D.Y. & Ying, Z. (1994). Semiparametric analysis of the additive risk model.
#'                              Biometrika; 81:61-71.
#'
#' @importFrom survival Surv
#' @importFrom stats residuals model.weights model.matrix model.extract terms runif
#' @export
#'
#' @examples
#' library(survival)
#' ### using the first 100 rows in nwtsco to build an additive hazards model
#' nwts<- nwtsco[1:100,]
#'
#' ### fit the additive hazards model to the data
#' ### the model-based standard errors are reported when setting robust = FALSE
#' fit1 <- ah(Surv(trel,relaps) ~ age + instit, data = nwts, robust = FALSE)
#' summary(fit1)
#'
#' ### fit the additive hazards model to the data with robust standard errors
#' fit2 <- ah(Surv(trel,relaps) ~ age + instit, data = nwts, robust = TRUE)
#' summary(fit2)
#'
#' ### when there are ties, break the ties by setting ties = 'break'
#' nwts_all <- nwtsco
#' fit3 <- ah(Surv(trel,relaps) ~ age + instit, ties = 'break', data = nwts_all, robust = TRUE)
#' summary(fit3)
#' ### users could break the ties on their own by
#' nwts_all$trel <- nwtsco$trel + runif(dim(nwts_all)[1],0,1)*1e-10
#' fit3 <- ah(Surv(trel,relaps) ~ age + instit, data = nwts_all, robust = TRUE)
#' summary(fit3)

ah <- function(formula, data, robust, weights, ties, seed = 20, ...) {
  Call <- match.call()
  indx <- match(c("formula", "data", "weights"), names(Call), nomatch = 0)
  if (indx[1] == 0) 
    stop("A formula argument is required")
  temp <- Call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")
  temp$formula <- if (missing(data)) 
    terms(formula) else terms(formula, data = data)
  
  mf <- eval(temp, parent.frame())
  if (nrow(mf) == 0) 
    stop("No (non-missing) observations")
  Terms <- terms(mf)
  Y <- model.extract(mf, "response")
  if (!inherits(Y, "Surv")) 
    stop("Response must be a survival object")
  type <- attr(Y, "type")
  if (type != "right" && type != "counting") 
    stop(paste("additive hazard model doesn't support \"", type, "\" survival data", sep = ""))
  weights <- model.weights(mf)
  nobs <- nrow(Y)
  X <- model.matrix(Terms, mf)[, -1]
  
  if (missing(robust) || is.null(robust)) 
    robust <- FALSE
  if (missing(weights) || is.null(weights)) 
    weights <- rep(1, nobs)
  
  if (any(duplicated(Y[Y[, "status"] == 1, "time"]))) {
    if (missing(ties) || ties != "break" || is.null(ties)) {
      stop(paste("There are ties in the survival time. Consider adding
                       an argument ties = 'break' in the function call. 
                       The program will break ties by adding a small random number 
                       to the survival time. The default seed number is 20."))
    } else if (ties == "break") {
      set.seed(seed)
      Y[, "time"] = runif(dim(Y)[1], 0, 1) * 1e-10 + Y[, "time"]
    }
  }
  
  fit <- ahaz::ahaz(Y, X, weights = weights, robust = robust)
  fit$coef <- summary(fit)$coef[, 1]
  fit$se <- summary(fit)$coef[, 2]
  
  fit$iA <- solve(fit$D)
  fit$var <- fit$iA %*% fit$B %*% fit$iA
  fit$resid <- residuals(fit)
  fit$npar <- fit$nvar
  
  # } else{ fit<-ah.fit(X, Y, wts= weights,robust=robust)
  # fit$resid<-residuals(fit,type='pseudoscore')*weights }
  
  # if (!is.null(fit$coef) && any(is.na(fit$coef))) { vars <- (1:length(fit$coef))[is.na(fit$coefs)] msg
  # <- paste('X matrix deemed to be singular; variable', paste(vars, collapse = ' ')) if (singular.ok)
  # warning(msg) else stop(msg) }
  
  fit$weights <- weights
  fit$robust = robust
  fit$terms <- Terms
  fit$formula <- formula(Terms)
  fit$model <- mf
  fit$call <- Call
  fit$nevent <- sum(Y[, ncol(Y)])
  fit$nobs <- nobs
  fit$data <- data
  
  fit$x <- X
  fit$y <- Y
  
  #  fit$xlevels <- xlevels 
  #  fit$assign <- assign if (robust) { 
  #  fit$naive.var <- fit$var
  #  temp<-residuals.ah(fit,type='pseudoscore')
  #  B.rob <- t(temp*weights) %*% (temp*weights)
  #  fit$var<-fit$iA%*%B.rob %*%fit$iA 
  #  fit$B.rob<-B.rob } 
  #  fit2 <- c(fit, list(x = X, y = Y, weights = weights)) 
  #  temp <- residuals.coxph(fit2, type = 'dfbeta', weighted = TRUE) 
  #  fit2$linear.predictors <- 0 * fit$linear.predictors 
  #  temp0 <- residuals.coxph(fit2, type = 'score', weighted = TRUE)
  #} 
  # fit$var <- t(temp) %*% temp u <- apply(as.matrix(temp0), 2, sum) 
  # fit$rscore <- coxph.wtest(t(temp0) %*% temp0, u, control$toler.chol)$test } 
  # if (length(fit$coefficients) && is.null(fit$wald.test)) {
  #   nabeta <- !is.na(fit$coefficients)
  #   if (is.null(init)) 
  #     temp <- fit$coefficients[nabeta] 
  #   else 
  #     temp <- (fit$coefficients - init[1:length(fit$coefficients)])[nabeta] 
  #   fit$wald.test <- coxph.wtest(fit$var[nabeta, nabeta], temp, control$toler.chol)$test
  # } 
  #   na.action <- attr(mf,'na.action')
  #   if (length(na.action)) fit$na.action <- na.action
  class(fit) <- "ah"
  fit
}
