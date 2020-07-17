#' Prediction Based on the Additive Hazards Model Fitted from Two-phase Sampling
#'
#' This function predicts a subject's overall hazard rates at given time points
#' based on this subject's covariate values. The prediction function is an object
#' from \code{\link{ah.2ph}}. The  estimating procedures follow Hu (2014).
#'
#' @param object an object of class inhering from 'ah.2ph'.
#' @param newdata a dataframe of an individual's predictors.
#' @param newtime  a given sequence of time points at which the prediction is performed.
#' @param ...  further arguments passed to or from other methods.
#'
#' @return A dataframe including the given time points, predicted hazards, their
#'  standard errors, their variances, the phase I component of the variance for
#'  predicted hazards and the phase II component of the variance.
#'
#' @seealso \code{\link{ah.2ph}} for fitting the additive hazards model with two-phase
#   sampling and \code{\link{nwtsco}} for the description of nwtsco dataset
#'
#' @importFrom survival Surv
#' @importFrom stats delete.response model.matrix terms
#' @export
#'
#' @references
#' Jie Hu (2014) A Z-estimation System for Two-phase Sampling with Applications
#'               to Additive Hazards Models and Epidemiologic Studies. Dissertation,
#'               University of Washington.
#'
#' @examples
#' library(survival)
#' ### fit an additive hazards model to  two-phase sampling data without calibration
#' fit1 <- ah.2ph(Surv(trel,relaps) ~ age + histol, data = nwts2ph,  R = in.ph2,  
#'                Pi = Pi, robust = FALSE, ties = 'break')
#'
#' ###  input the new data for prediction
#' newdata <- nwtsco[101,]
#' ###  based on the fitted model fit1, perform prediction at time points t =3 and t= 5
#' predict(fit1, newdata, newtime = c(3,5))
#'
#' ### fit an additve hazards model to  two-phase sampling data with calibration
#' ### The calibration variable is instit
#' fit2 <- ah.2ph(Surv(trel,relaps) ~ age + histol, data = nwts2ph, R = in.ph2, Pi = Pi,
#'                                    ties = 'break', robust = FALSE, calibration.variables = "instit")
#' ### based on the fitted model fit2, perform prediction at time points t =3 and t= 5
#' predict(fit2, newdata, newtime = c(3,5))

predict.ah.2ph <- function(object, newdata, newtime, ...) {
  
  # phase I complete data object
  object1 <- object$fit.pha2
  ######################### creat model matrix for the new data ##########
  tt <- terms(object1)
  if (!inherits(object, "ah.2ph")) {
    warning("calling predict.2ph(<fake-ah-object>) ...")
  }
  Terms <- delete.response(tt)
  # identify factor variable
  factor.var.list <- names(object1$model)[sapply(object1$model,
                                                 is.factor) == T]
  # assign levels to newdata variable according to the original data
  if (length(factor.var.list) != 0) {
    for (i in 1:length(factor.var.list)) {
      factor.var <- factor.var.list[i]
      newdata[[factor.var]] <- factor(newdata[[factor.var]],
                                      levels = levels(object1$model[[factor.var]]))
    }
  }
  # retrieve the new model matrix based on the new data
  pred <- model.matrix(Terms, newdata)
  # delete the intercept
  pred <- as.numeric(pred[, -1])
  ###########################################################################
  Pi.pha2 <- object$Pi.pha2
  if (object1$robust) {
    warning("robust is set to be FALSE, which means the additive hazards model you input is assumed misspecified. Prediction based on this model may not be valid")
  }
  ## No calibration
  if (!length(object$calibration.variables)) {
    
    pred.result <- predict.ah(object1, newdata, newtime,
                              level = 0.95)
    L <- pred.result$L
    L.var.pha1 <- pred.result$L.se^2
    
    ### calculate the  variance due to phase II sampling
    ingr <- ingredients(object1)
    
    ############################ find the positions of a list of given newtime on the
    ############################ original timeline######################
    t.pos <- NULL
    t.unique <- ingr$t.unique
    
    for (i in 1:length(newtime)) {
      # newtime s is less than t_1, then we report
      if (newtime[i] < t.unique[1]) {
        print(paste("The data used for building the ah model does not have
                               enough information for predicting such a small t, t=",
                    newtime[i]))
      } else {
        # find the position of newtime[i] on the original timeline
        c <- rank(c(newtime[i], t.unique), ties.method = "max")[1]
        # when newtime s equals one of the t.unique t_k, the above
        # rank function returns k+1 on the original timelines when
        # newtime s in between t_k and t_{k+1}, it returns k+1 on
        # the original timelines when newtime s = t_1, it returns 1
        # thus
        t.pos[i] <- ifelse(c == 1, 1, c - 1)
      }
    }
    
    ### predicted outcomes are calculated based on fitting an ah
    ### model to only phase II data using assigned weights
    
    # Calculate the phase II variance
    L.resid <- do.call("cook.L.resid", c(ingr, list(time.pos = ingr$match.event)))
    ingr$wts <- sqrt(1 - Pi.pha2)
    L.var.pha2 <- do.call("L.rvar", c(ingr, object1,
                                      list(t.pos = t.pos, pred = pred, newtime = newtime,
                                           L.resid = L.resid)))
    
  } else {
    # Calibration
    pred.result <- predict.ah(object1, newdata, newtime,
                              level = 0.95)
    L <- pred.result$L
    
    ### Calculate the 1st component of the variance
    L.var.pha1 <- pred.result$L.se^2
    wts.cal <- object1$weights
    ingr <- ingredients(object1)
    # retrieve calibration variables
    aux<-as.matrix(object$calibration.variables)
    aux.pha2<-aux[object$R==1,]
    P<-t(aux)%*%(aux)
    t.pos <- NULL
    t.unique <- ingr$t.unique
    
    for (i in 1:length(newtime)) {
      # newtime s is less than t_1, then we report
      if (newtime[i] < t.unique[1]) {
        print(paste("The data used for building the ah model does not
                               have enough information for predicting such a small t, t=",
                    newtime[i]))
      } else {
        # find the position of newtime[i] on the original timeline
        c <- rank(c(newtime[i], t.unique), ties.method = "max")[1]
        # when newtime s equals one of the t.unique t_k, the above
        # rank function returns k+1 on the original timelines when
        # newtime s in between t_k and t_{k+1}, it returns k+1 on
        # the original timelines when newtime s = t_1, it returns 1
        # thus
        t.pos[i] <- ifelse(c == 1, 1, c - 1)
      }
    }
    
    ###Calculate the variance of the predicted value See page 109 ####
    #### of Jie Hu' thesis for the detailed formula #######
    ## retrieve the martingale integral, i.e.  part of \psi in
    ## the formula on page 109 of Jie Hu' thesis
    
    L.resid <- do.call("cook.L.resid", c(ingr, 
                                         list(time.pos = ingr$match.event)))
    
    L.var.pha2 <- do.call("L.rvar.calibration", 
                          c(ingr, 
                            object1, 
                            list(t.pos = t.pos, 
                                 pred = pred, 
                                 newtime = newtime, 
                                 L.resid = L.resid, 
                                 new.wts = sqrt((1 - Pi.pha2)/(wts.cal *Pi.pha2)), 
                                 aux.pha2 = aux.pha2, 
                                 P = P, 
                                 wts.cal = wts.cal)))
  }

# L.var.pha2<-L.rvar(match.event=match.event,new.wts=sqrt(1-Pi.pha2),L.ncol=L.ncol,
# eta.cum=eta.cum, resid=resid, L.resid=L.resid,iA=iA, B=B)

  L.var = L.var.pha1 + L.var.pha2
  L.se = sqrt(L.var)
  return(data.frame(L = L, L.se = L.se, L.var = L.var, L.var.pha1 = L.var.pha1,
                  L.var.pha2 = L.var.pha2))
}
