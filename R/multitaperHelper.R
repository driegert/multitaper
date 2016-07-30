##     The multitaper R package
##     Multitaper and spectral analysis package for R
##     Copyright (C) 2011 Karim Rahim 
##
##     Written by Karim Rahim and Wesley Burr.
## 
##     This file is part of the multitaper package for R.
##     http://cran.r-project.org/web/packages/multitaper/index.html
## 
##     The multitaper package is free software: you can redistribute it and 
##     or modify it under the terms of the GNU General Public License as 
##     published by the Free Software Foundation, either version 2 of the 
##     License, or any later version.
##
##     The multitaper package is distributed in the hope that it will be 
##     useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
##     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.
##
##     You should have received a copy of the GNU General Public License
##     along with multitaper.  If not, see <http://www.gnu.org/licenses/>.
##
##     If you wish to report bugs please contact the author:
## 
##     Karim Rahim
##     karim.rahim@gmail.com


.mweave <-  function (x,dw,swz,ndata,nord,ssqswz,dt_) {
    
    out <- .Fortran("mweave", as.double(x), as.double(dw),
                    as.double(swz), as.integer(ndata),
                    as.integer(nord), as.double(ssqswz),
                    cntr=double(1), as.double(dt_),
                    spz=double(1), varc=double(1),
                    PACKAGE='multitaper')
    return(list(cntr=out$cntr, spz=out$spz, varc=out$varc))
}

.HF4mp1 <- function(cft, swz, nord, ssqswz) {

    ## ######################################
    ## The notation and function names were chosen
    ## to map to original fortran (f77) code.
    ## Note to obtain:  swz <- apply(dw, 2, sum)
    ## swz is the zeroth frequency Fourier transform of the
    ## Slepian sequences. It is H_k(0) from P ercival and Walden (1993)
    ## pages 497--399.
    ## (just to define dw) dw <- dpssIN$v*sqrt(deltaT)    
    ## Vectorized from original F77 code
    ## Equation (13.5) of:
    ##   Thomson, D.J. Spectrum Estimation and Harmonic Analysis,
    ##   Proceedings of the IEEE, 1982.

    cmv <- (cft %*% swz) /ssqswz
    ssqave <-  (Mod(cmv)^2)*ssqswz
    swz <- as.matrix(swz)
    
    ssqres <- apply( Mod(cft - (cmv %*% t(swz)))^2,
                    1, sum)
    F_<- (nord-1)*ssqave/ssqres
    
    return(list(Ftest=F_,cmv=cmv))
}


.mw2wta <- function(sa, nfreq, nord,
                    var, dt_, ev, evp=(1-ev),
                    tol=.03, maxadaptiveiteration=100) {

    ## this is equation (5.3) and (5.4) form 
    ##   Thomson, D.J. Spectrum Estimation and Harmonic Analysis,
    ##   Proceedings of the IEEE, 1982.
    ## note that the weights are squared, they are |d_k(f)^2 from equation
    ## (5.4)

    out <- .Fortran("mw2wta", as.double(sa),
                    wt=matrix(as.double(0), nfreq, nord),
                    as.integer(nfreq), as.integer(nord),
                    s=double(nfreq), as.double(ev), as.double(evp),
                    dofs=double(nfreq), dofav=double(1),
                    as.double(var), as.double(dt_),
                    as.double(tol),
                    as.integer(maxadaptiveiteration),
                    mxiter=integer(1), aviter=double(1),
                    PACKAGE='multitaper')
    
    return(list(s=out$s, wt=out$wt, dofs=out$dofs, dofav=out$dofav,
                mxiter=out$mxiter, aviter=out$aviter))
}

.mw2jkw <-  function(sa, nfreq, nord, var, dt_, ev,
                     evp=(1-ev), tol=.03,
                     maxadaptiveiteration=100) {
    
    nordP2 <-  nord+2
    out <- .Fortran("mw2jkw", as.double(sa),
                    wt=matrix(as.double(0), nfreq, nord),
                    as.integer(nfreq), as.integer(nord),
                    s=double(nfreq), as.double(ev),
                    as.double(evp), dofs=double(nfreq),
                    dofav=double(1), as.double(var),
                    as.double(dt_), as.double(tol),
                    sjk=double(nordP2), varjk=double(nfreq),
                    bcjk=double(nfreq),
                    matrix(as.double(0), nord, nordP2),
                    double(nordP2), double(nord),
                    as.integer(maxadaptiveiteration),
                    mxiter=integer(1),
                    PACKAGE='multitaper')
                   
    
    return(list(s=out$s, wt=out$wt, dofs=out$dofs, dofav=out$dofav,
                mxiter=out$mxiter,
                varjk=out$varjk, bcjk=out$bcjk, sjk=out$sjk))
}

.qsF <- function(nFreqs,nFFT,k,cft,useAdapt,kadapt) {

    out <- .Fortran("quickSineF", as.integer(nFreqs),
                    as.integer(nFFT), as.integer(k),
                    cft=cft, as.logical(useAdapt),
                    kadapt=matrix(data=as.double(kadapt),nrow=nFreqs,ncol=1),
                    spec=matrix(data=double(nFreqs),nrow=nFreqs,ncol=1),
                    PACKAGE='multitaper')

    return(list(spec=out$spec))
}

.cF <- function(n,v) {

  out <- .Fortran("curbF",as.integer(n),as.double(v),PACKAGE='multitaper')
  opt <- out[[2]]
  return(list(opt=opt))
}

.nF <- function(n,i1,i2,s) {

   out <- .Fortran("northF",as.integer(n),as.integer(i1),
                   as.integer(i2),sx=matrix(data=as.double(s),nrow=n,ncol=1),
                   ds=double(1), dds=double(1),
                   PACKAGE='multitaper')

   return(list(ds=out$ds, dds=out$dds))
}

.adaptSine <- function(ntimes, k, nFreqs, sx, nFFT, cft, df, fact) {

    out <- .Fortran("adapt",as.integer(ntimes),as.integer(k),
                    as.integer(nFreqs),sx=matrix(data=as.double(sx),nrow=nFreqs,ncol=1),
                    as.integer(nFFT), cft=cft, as.double(df),
                    kopt=double(nFreqs),fact=as.double(fact),PACKAGE='multitaper')

    return(list(spec=out$sx,kadapt=out$kopt))

}

# things that are currently hardcoded:
## returnZeroFreq=TRUE, Ftest=FALSE, jkCIProb=0.95, 
.determineBandwidth <- function(timeSeries,
                                nFFT,
                                dpssIN,
                                adaptiveWeighting,
                                maxAdaptiveIterations,
                                n,
                                deltaT,
                                sigma2,
                                series,
                                dtUnits,
                                ...){
  # initial values:
  nw <- seq(2, 30)
  w <- nw / n
  k <- floor(2*nw)
  
  ### NEED TO FIX: How to determine break points in an unsupervised manner
  # if (is.complex(timeSeries)){
  #   probs <- seq(0, 1, 0.1/2) ## this
  # } else {
  #   probs <- seq(0, 1, 0.1) ## and this..
  # }
  
  mse <- rep(0, length(nw))
  bias2 <- rep(0, length(nw))
  Evar <- rep(0, length(nw))
  # Mjkvar <- rep(0, length(nw))
  for (i in 1:length(nw)){
    ## This is here from testing
    # mtm.obj <- spec.mtm(timeSeries, nw = nw[i], k = k[i], nFFT = nFFT, jackknife = TRUE, plot=F
    # , deltat = deltaT, adaptiveWeighting = TRUE)
    mtm.obj <- .spec.mtm.dpss(timeSeries=timeSeries,
                              nw=nw[i], k=k[i], nFFT=nFFT,
                              dpssIN=FALSE, returnZeroFreq=TRUE, #dpssIN=FALSE may not work...
                              Ftest=FALSE, jackknife=TRUE, jkCIProb=0.95,
                              adaptiveWeighting = adaptiveWeighting,
                              maxAdaptiveIterations=maxAdaptiveIterations,
                              returnInternals=FALSE, # might not need these?
                              n=n, deltaT=deltaT, sigma2=sigma2, series=series,
                              dtUnits=dtUnits, ...)
    
    # setup the spline basis and derivatives:
    ####
    # subsample the frequencies - every 2w
    deltaf <- mtm.obj$freq[2] - mtm.obj$freq[1]
    
    # original code:
    # freqSubIdx <- seq(1, length(mtm.obj$freq), by = round(2*w[i] / (deltaT*deltaf))) # not sure if rounding will cause issues... probably not
    # spline.freq <- mtm.obj$freq[freqSubIdx]
    
    # Try all the frequencies:
    freqSubIdx <- 1:length(mtm.obj$freq)
    spline.freq <- mtm.obj$freq

    # expected jackknife variance - Haley, Eq (43)
    expJkVar <- ( ((k[i]-1)^3)/(k[i] * (k[i] - 0.5)) ) * ( (2 / (k[i]-2)^2) + 0.5*(trigamma((k[i]-1)/2) - trigamma((k[i]-2)/2)) )
    Evar[i] <- expJkVar
    
    #################### Choice 1 - b-splines - 5th order ###############
    # breaks <- .chooseKnots(mtm.obj$spec[freqSubIdx], spline.freq, spar = length(freqSubIdx) * expJkVar)
    # 
    # # fit the splines to the spectrum
    # splineModel <- lm(mtm.obj$spec[freqSubIdx] ~ bsplineS(x = spline.freq, norder = 6
    #                                           , breaks = breaks, nderiv = 0) - 1) # no intercept
    # # 2nd derivative of spline fit
    # splineModel.d2 <- bsplineS(x = spline.freq, norder = 6
    #                            , breaks = breaks, nderiv = 2) %*%
    #   matrix(coefficients(splineModel), ncol = 1)
    #####################################################################
    
    ### choice 2 - penalized smoothing splines in R ###
    splineModel <- smooth.Pspline(spline.freq, mtm.obj$spec[freqSubIdx], norder = 3, method = 4)
    splineModel.d2 <- predict(splineModel, xarg = spline.freq, nderiv = 2)
    #########################################
    
    ### Choice 3 - cubic smoothing splines
    # splineModel <- smooth.spline(spline.freq, mtm.obj$spec[freqSubIdx])
    # splineModel.d2 <- predict(splineModel, x = spline.freq, deriv = 2)$y
    
    ### REMOVE THIS AFTER ###
    bias2[i] <- mean( log( (mtm.obj$spec[freqSubIdx] + abs(w[i]^3 * splineModel.d2 / 3)) / mtm.obj$spec[freqSubIdx] )^2 )
    logbias2 <- log( (mtm.obj$spec[freqSubIdx] + abs(w[i]^3 * splineModel.d2 / 3)) / mtm.obj$spec[freqSubIdx] )^2
    #########################
    mse[i] <- mean(logbias2 + mtm.obj$mtm$jk$varjk[freqSubIdx])
    # mse[i] <- mean(logbias2 + expJkVar)
  }
  
  nw[which(mse == min(mse))]
}

# choose the correct number of knots to use in the bsplines
## based on the argument 's' from: 
# http://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.UnivariateSpline.html
.chooseKnots <- function(spec, freqSamp, spar = NULL){
  nknots <- 0
  if (is.null(spar)){
    spar <- 3.2 #### Make this Eq (43) in Haley (expected jk variance)
  }
  
  probs <- seq(0, 1, length.out = nknots + 2)
  breaks <- quantile(freqSamp, probs = probs)
  
  splineModel <- lm(spec ~ bsplineS(x = freqSamp, norder = 6
                                    , breaks = breaks, nderiv = 0) - 1)
  
  while (sum((fitted.values(splineModel) - spec)^2) > spar){
    if (nknots + 6 >= length(spec)){
      warning("Could not get the correct number of knots based on smoothing parameter.
              All points used.")
      
      break
    }
    
    nknots <- nknots + 1
    probs <- seq(0, 1, length.out = nknots + 2)
    breaks <- quantile(freqSamp, probs = probs)
    
    splineModel <- lm(spec ~ bsplineS(x = freqSamp, norder = 6
                                      , breaks = breaks, nderiv = 0) - 1)
  }
  
  breaks
}