#######################################################
# Author: Kate HU
# Date: July 27th, 2016
# Task: write function when assuming there may be ties
# To do: interval censoring
#######################################################

require("ahaz")
require("matrix")

ah.fit<-function(Y,Z,wts,robust){

  n.obs <- dim(Z)[1]
  n.par <- dim(Z)[2]
  if (missing(wts) || is.null(wts)) wts <- rep(1, n.obs) #?

  t.fail <- Y[,1]
  censor <- Y[,2]


  # Centralized the covariates
  # Z <- scale( Z,center= rep(1,n.par),scale=F)

  # Divide the timeline to small intervals according to  events or censoring

  t.sorted <- sort(t.fail)
  t.unique <- unique(t.sorted)
  #  The length of eta
  eta.ncol <- length(t.unique)
  t.diff <- t.unique - c(0,t.unique[-eta.ncol])

  # Section 1: calculate eta(t) ------------------------------------------------
  # eta(t):  average Z for patients still at risk at time point t
  # eta(t) only changes when a subject fails or censored, i.e., t.unique


  ### We calculate the numerator and the denomiators of eta separtedly
  ### We create the spaces for the numerator and the denomiators of eta
  ### Z_i is a n.par x 1 vector , therefore


  ### We calcuate the numerator and the denomiators of eta
  ### At each unique failure time, we calcuate for each covariate(age,sex,...)
  ### Then we calculate for all unique failure times
  z.death <- eta.num <- eta.den <- matrix(data=NA,n row=n.par, ncol = eta.ncol)
  n.death<-NULL
  eta.den.vec<-NULL
  for ( i in seq_along(eta.ncol)){
  	 idx1 <- t.fail >= t.unique[i]
     idx2 <- which(t.fail == t.unique[i])

     n.death[i] <- sum((censor*wts)[idx2])
     eta.den.vec[i]<-sum(wts[idx1])

  	 for (j in 1:n.par){
	   eta.num[j,i] <-sum((Z*wts)[idx1,j])

         ########################## IF we calculate Lambda ###############
    ### calculate the cumulative eta
    ### The sum of Z of the subjests who have died at the jth dicrete time

     z.death[j,i]<-sum((Z[,j]*censor*wts)[idx2])
  		}
   eta.den[,i]<-eta.den.vec[i]
  }



  ### Calculate eta, which is a  n.par X eta.col matrix
  dLambda1<-n.death/eta.den.vec
  #print(dLambda1)
  eta<-eta.num/eta.den
  ###############################################################################

  ###############################################################################
  ###Calculate theta. It doens't dependent on t

  ### Calculate the numerator of theta



  ###Generate a new eta matrix(vector) corresponding to each observation such that
  ### eta.each.person_i is the average of Z for patients still at risk at
  ###the failure time of  individual i

  ### create the space fo this  n.obs x n.par matrix (vector)

  #### Individuals


  ### find which
  match.eta<-match(t, t.unique)
  eta.i<-matrix(as.vector(eta[,match.eta]),ncol=n.par,byrow=TRUE)


  ###Calculate theta.numerator


  	   # z.cen.1<-(Z-eta.i)*wts*censor
   A=fit$D




    #theta.denominator<-Reduce("+",theta.denominator.timepoint)

  ###Calculate theta.denominator


 # A<-matrix(0,n.par,n.par)

  #for ( i in 1:eta.ncol){

  #  center<-Z[failure.time>= time.interval[i]]-		  #eta.each.person[failure.time>= time.interval[i]]

   #   idx1<-t>=  t.unique[i]
    #  rep.times<-sum(idx1)
     # eta.stack<-matrix(rep(as.vector(eta[,i]),rep.times), byrow=TRUE,ncol=n.par)
     # temp<-(Z[idx1,]-eta.stack)*sqrt(wts[idx1])

    	#A<-
    	#A+t(temp)%*%temp *t.diff[i]

 #}





   # theta.denominator<-Reduce("+",theta.denominator.timepoint)

  ### Calculate theta
    iA<-solve(A)
    theta.num<-fit$d
    theta<-iA%*%theta.num



  ###Calculate the variance
  ###B part


  B<-fit$B
  var<-iA%*% B %*% iA


 eta.cum<-as.matrix(apply(t(eta)*t.diff,2,cumsum))



 #for (k in 1: n.obs){
  #       int.upper<-match(t[k],t.unique)
   #     for (l in 1: int.upper){
    #        temp<-(Z[k,]-eta[,l])
     #       A.i[,,k] <-A.i[,,k]+temp%*%t(temp)*t.diff[l]*wts[k]
      #      z.cen.2[k,]<-z.cen.2[k,]+temp*dLambda1[l]*wts[k]

       # }

 # cumhaz<-function(fit){
#    eta.cum<-matrix(data=NA,nrow=fit$n.par,ncol=fit$eta.ncol)
#     Lambda1<-rep(0,fit$eta.ncol)
#     eta.cum[,1]<-fit$eta[,1]*fit$t.diff[1]
#     Lambda1[1]<-fit$dLambda1[1]
#     for ( i in 2:fit$eta.ncol){
#         eta.cum[,i]<-eta.cum[,i-1]+fit$eta[,i]*fit$t.diff[i]
#         Lambda1[i]<-Lambda1[i-1]+ fit$dLambda1[i]
#     }
#     print(Lambda1)
#     return(Lambda1-t(eta.cum)%*%fit$coef)
#
# }

effects<-as.vector(Z%*%theta)

  fit<-list(coef=theta,
              var=var,
              A=A,
               iA=iA,
                eta=eta,
                eta.cum=eta.cum,
                 z.cen.1=z.cen.1,
                 match.eta=match.eta,
                 eta.i=eta.i,
                 eta.den.vec=eta.den.vec,
                 dLambda1=dLambda1,
                 n.par=n.par,
                 n.obs=n.obs,
                 z.death=z.death,
                 n.death=n.death,
                 eta.ncol=eta.ncol,
                 censor=censor,
                 time=t,
                 t.diff=t.diff,
                 t.unique=t.unique,
                 effects=effects,
                  robust=robust)


 lambda0<-dLambda1-as.vector(t(eta)%*%theta*t.diff)
 Lambda0<-cumsum(lambda0)



 M<-(censor-Lambda0[match.eta]-effects*t)



z.resid<-Z*M
eta.resid2<-apply(t(eta)*lambda0,2,cumsum)


eta.resid3<-eta.cum[match.eta,]*effects


eta.resid<-eta.i*censor-eta.resid2[match.eta,]-eta.resid3

 ###nobs x n.par
   z.cen.resid<-(z.resid-eta.resid)

    #A.i<-array(rep(0,n.par*n.par*n.obs),dim=c(n.par,n.par,n.obs))
   #for (k in 1: n.obs){
    #     int.upper<-match(t[k],t.unique)
     #   for (l in 1: int.upper){
      #      temp<-(Z[k,]-eta[,l])
       #     A.i[,,k] <-A.i[,,k]+temp%*%t(temp)*t.diff[l]*wts[k]
        #    z.cen.2[k,]<-z.cen.2[k,]+temp*dLambda1[l]*wts[k]

        #}
    #}


      fit$Lambda0<-Lambda0
      fit$lambda0<-lambda0

      fit$M<-M
     fit$resid<-z.cen.resid

      #z.cen.3<-as.matrix(apply(A.i,3,function(s){s%*%fit$coef}))
      #print(A.i[,,1]%*%theta)
      #print(A.i[,,2]%*%theta)
      #print(dim(z.cen.3))
      #print(dim(z.cen.2))
      #fit$z.cen<-(z.cen.1-z.cen.2-z.cen.3)

      #print(z.cen.2)


 return(fit)
}


residuals.ah<-function (object, type = c("martingale", "pseudoscore","pseudodfbeta","martingale.increment")){

    type <- match.arg(type)
    if (type=="martingale") return(object$M)
    if (type=="pseudoscore") return(object$resid)
    if (type=="pseudodfbeta") return(object$iA%*%object$resid)
    #if (type=="martingale.increment"){
     # dM<-matrix(0,nrow=object$n.obs,ncol=object$eta.ncol)
      #for  ( i in 1: object$eta.ncol){
       #    dM[object$match.eta<=i,] <- (object$censor-object$lambda0[i]-object$effects*object$t.diff[i])[object$match.eta<=i]
        #  }
      #return(dM)
    #}
}
