source("CFO_tox_utils.R")
source("ACFO_tox_utils.R")
# dlt.time: time to event, if 0, no event
# delta.time, time interval of patient arrive


# Below functions are to impute missing y 
#------------------------------------------------------------------------------------------
fracImpute <- function(enter.times, dlt.times, current.time, tau){ 
  
  #args:
  # enter.times: The enter times of the patients, a vector 
  # dlt.times: The DLT times of the patients, if no DLT, 0, a vector
  # current.time: Current time point: a value
  # tau: Observing window size
  
  #return:
  # ym: Imputed y for the patient with no DLT and follow-up time < tau
  
  assesstime = enter.times+tau;	
  dlt.times[dlt.times==0]= tau+1;
  yo = (dlt.times<=tau)*(assesstime<=current.time)+(dlt.times<=(current.time-enter.times))*(current.time<assesstime);		
  No.impute <- FALSE
  if (sum(yo)==0)	{
    No.impute <- TRUE
    ym <- yo
    #stop("fraction design takes effect when at least one DLT has been observed")
  }
  if (sum(yo)!=0){			
    otime = yo*dlt.times+(1-yo)*((current.time-enter.times)*(current.time<assesstime)+tau*(assesstime<=current.time))			
    kmfit = survfit(Surv(otime,yo)~1)	
    #plot(kmfit,xlab="Time (days)",ylab="Survival probability",
    #     main="Kaplan-Meier estimate",cex.lab=1.3,cex.main=1.3,cex.sub=1.3)	
    ym = yo
    
    for (i in 1:length(yo)){
      if (current.time<assesstime[i] & yo[i]==0){
        ym[i]=(kmfit$surv[tail(which(kmfit$time<=(current.time-assesstime[i]+tau+0.00001)),n=1)]- kmfit$surv[tail(which(kmfit$time<=tau),n=1)])/
          kmfit$surv[tail(which(kmfit$time<=(current.time-assesstime[i]+tau+0.00001)),n=1)]
      }
    }
    
  }
  obsIdxs <- current.time >= assesstime
  obsIdxs[yo==1] <- TRUE
  #print(No.impute)
  #print(ym)
  #print(yo)
  
  res <- list(y.impute=ym, y.raw=yo, obsIdxs=obsIdxs, No.impute=No.impute)
  res
}

TITEImpute2.one <- function(followup.times, tau, y, n, prior.paras){
  #args:
  #   followup.times: The follow-up times of the pending patients at the dose level
  #   tau: Observing window size
  #   y: Num of Observed DLT at the dose level 
  #   n: Num of patients with observed results at the dose level
  #   prior.paras: a vector of 2, prior when estimating ptilde
  
  #return: 
  #  ym: imputed y 
  
  p.tilde <- (y+prior.paras[1])/(n+sum(prior.paras))
  #ym <- p.tilde * (1-followup.times/tau)
  ym <- p.tilde * (1-followup.times/tau) /((1-p.tilde)+p.tilde * (1-followup.times/tau))
  #    ym <- p.tilde * (1-followup.times/tau) /(1-p.tilde)
  #    ym[ym >1] <- 1 # trunc the value
  ym
}

TITEImpute2 <- function(enter.times, dlt.times, current.time, tau, dose.levels, ndose, prior.paras){
  #args:
  # enter.times: The enter times of the patients, a vector 
  # dlt.times: The DLT times of the patients, if no DLT before tau, 0, a vector
  # current.time: Current time point: a value
  # tau: Observing window size
  # dose.levels: dose level for each subject
  # ndose: num of dose levels
  # prior.paras: a vector of 2, prior when estimating ptilde
  
  #return:
  # ym: Imputed y for the patient with no DLT and follow-up time < tau
  
  assesstime = enter.times + tau;	
  followup.times <- current.time - enter.times
  dlt.times[dlt.times==0]= tau+1;
  yo <- (dlt.times<=tau)*(assesstime<=current.time)+(dlt.times<=followup.times)*(current.time<assesstime);		
  obsIdxs <- current.time >= assesstime
  obsIdxs[yo==1] <- TRUE
  ym <- yo
  for (i in 1:ndose){
    doseIdxs <- dose.levels == i
    if (sum(1-obsIdxs[doseIdxs]!=0)){
      y <- sum(yo[doseIdxs])
      n <- sum(doseIdxs[obsIdxs==1])
      kpidxs <- doseIdxs & (obsIdxs!=1)
      ym.part <- TITEImpute2.one(followup.times[kpidxs], tau, y, n, prior.paras)
      ym[kpidxs] <- ym.part
    }
  }
  res <- list(y.impute=ym, y.raw=yo, obsIdxs=obsIdxs)
  res
}

#------------------------------------------------------------------------------------------
# The function is to obtain the DLT results (with TITE) for each subject
#------------------------------------------------------------------------------------------
gen.tite<-function(dist=1, n, pi, tau=1, alpha=0.5){
  #args:
  #   dist: TITE distribution, 1-uniform, 2-weibull, 3-log-log
  #   n: Num of subjects to generate
  #   pi: Target DLT rate, pi=Pr(T<=tau)
  #   tau: Maximal window size
  #   alpha: Parameter for generate time
  #Return:
  #   if no DLT, tox.t=0
  ############ subroutines ############
  weib<-function(n, pi, pihalft)
  {
    ## solve parameters for Weibull given pi=1-S(T) and phalft=1-S(T/2)
    alpha = log(log(1-pi)/log(1-pihalft))/log(2);
    lambda = -log(1-pi)/(tau^alpha);
    t = (-log(runif(n))/lambda)^(1/alpha);
    return(t);
  }
  
  llogit<-function(n, pi, pihalft)
  {
    ## solve parameters for log-logistic given pi=1-S(T) and phalft=1-S(T/2)
    alpha = log((1/(1-pi)-1)/(1/(1-pihalft)-1))/log(2);
    lambda = (1/(1-pi)-1)/(tau^alpha);
    t = ((1/runif(n)-1)/lambda)^(1/alpha);
    return(t);
  }  
  ############ end of subroutines ############
  
  
  tox = rep(0, n);
  t.tox = rep(0, n);
  
  #### uniform
  if(dist==1) {  # 50% event in (0, 1/2T)
    tox = rbinom(n, 1, pi);
    ntox.st = sum(tox);
    t.tox[tox==1]=runif(ntox.st, 0, tau);
    t.tox[tox==0]=0;
  }
  #### Weibull
  if(dist==2)
  {
    pihalft = alpha*pi;  # alpha*100% event in (0, 1/2T)
    t.tox = weib(n, pi, pihalft);
    tox[t.tox<=tau]=1;
    ntox.st = sum(tox);
    t.tox[tox==0]=0;
  }
  #### log-logistic
  if(dist==3)
  {
    pihalft = alpha*pi;  # alpha*100% event in (0, 1/2T)
    t.tox = llogit(n, pi, pihalft);
    tox[t.tox<=tau]=1;
    ntox.st = sum(tox);
    t.tox[tox==0]=0;
  }
  return(list(tox=tox, t.tox=t.tox, ntox.st=ntox.st));
}

#------------------------------------------------------------------------------------------
# Propose the next dose level with CFO design
#------------------------------------------------------------------------------------------
CFOlateonset.next.dose <- function(curDose, phi, tau, impute.method, 
                                   enter.times, dlt.times, current.t, doses, tover.doses, simu=FALSE, add.args=list()){
  #args:
  #   curDose: the current dose level
  #   phi: Target DLT rate
  #   tau: maximal observing win size
  #   impute.method: impute method: frac: fractional design, TITE: TITE way
  #   enter.times: enter.times of each subject
  #   dlt.times: dlt.times of each subject, let 0 if no DLT
  #   current.t: The current time 
  #   doses: Dose level for each subject
  #   tover.doses: over dose index for each dose, 1 over-toxic, 0 safe
  #   simu: Whether simulation or not, if TRUE, return tover.dose also
  #return:
  #   dose: Recommend dose, 0 early stopping
  
  ndose <- length(tover.doses)
  ## Obtain effective results
  if (impute.method == "frac"){
    impute.res <- fracImpute(enter.times, dlt.times, current.t, tau)
    y.raw <- impute.res$y.raw
    y.impute <- impute.res$y.impute
    if (impute.res$No.impute){
      cn <- sum(doses==curDose)
      # should be cn <- sum(doses==curDose & impute.res$obsIdxs)
      cy <- sum(y.raw[doses==curDose])
    }
    else{
      cn <- sum(doses==curDose)
      cy <- sum(y.impute[doses==curDose])
    }
  }else if(impute.method == "TITE"){
    impute.res <-  TITEImpute2(enter.times, dlt.times, current.t, tau, doses, ndose, c(phi/2, 1-phi/2))
    y.raw <- impute.res$y.raw
    y.impute <- impute.res$y.impute
    
    cy <- sum(y.impute[doses==curDose])
    cn <- sum(doses==curDose)
  }
  
  ### early stop
  add.args <- c(list(y=cy, n=cn), add.args)
  if (overdose.fn(phi, add.args)){
    tover.doses[curDose:ndose] <- 1
  }
  
  if (tover.doses[1] == 1){
    dose <- 0
  }else{
    ### the results for current 3 dose levels
    if (curDose==1){
      cys <- c(NA, cy, sum(y.impute[doses==2]))
      cns <- c(NA, cn, sum(doses==2))
      cover.doses <- c(NA, tover.doses[1:2])
    }else if (curDose==ndose){
      cys <- c(sum(y.impute[doses==(curDose-1)]), cy, NA)
      cns <- c(sum(doses==(curDose-1)), cn, NA)
      cover.doses <- c(tover.doses[(curDose-1):curDose], NA)
    }else {
      cys <- c(sum(y.impute[doses==(curDose-1)]), cy, sum(y.impute[doses==(curDose+1)]))
      cns <- c(sum(doses==(curDose-1)), cn, sum(doses==(curDose+1)))
      cover.doses <- tover.doses[(curDose-1):(curDose+1)]
    }
    
    dose.chg <- make.decision.CFO.fn(phi, cys, cns, add.args$alp.prior, add.args$bet.prior, cover.doses) - 2
    dose <- dose.chg + curDose
  }
  
  if (simu){
    res <- list(dose=dose, tover.doses=tover.doses)
    return(res)
  }else{
    res <- list(dose=dose)
    return(res)
  }
}
#------------------------------------------------------------------------------------------
# Propose the next dose level with aCFO design
#------------------------------------------------------------------------------------------
aCFOlateonset.next.dose <- function(curDose, phi, tau, impute.method, 
                                   enter.times, dlt.times, current.t, doses, tover.doses, simu=FALSE, add.args=list()){
  #args:
  #   curDose: the current dose level
  #   phi: Target DLT rate
  #   tau: maximal observing win size
  #   impute.method: impute method: frac: fractional design, TITE: TITE way
  #   enter.times: enter.times of each subject
  #   dlt.times: dlt.times of each subject, let 0 if no DLT
  #   current.t: The current time 
  #   doses: Dose level for each subject
  #   tover.doses: over dose index for each dose, 1 over-toxic, 0 safe
  #   simu: Whether simulation or not, if TRUE, return tover.dose also
  #return:
  #   dose: Recommend dose, 0 early stopping
  
  ndose <- length(tover.doses)
  tys = NULL
  tns = NULL
  ## Obtain effective results
  if (impute.method == "frac"){
    impute.res <- fracImpute(enter.times, dlt.times, current.t, tau)
    y.raw <- impute.res$y.raw
    y.impute <- impute.res$y.impute
    if (impute.res$No.impute){
      for (i in 1:ndose){
        tys <- c(tys, sum(y.raw[doses==i]))
        tns <- c(tns, sum(doses==i))
      } 
      # should be cn <- sum(doses==curDose & impute.res$obsIdxs)
    }else{
      for (i in 1:ndose){
        tys <- c(tys, sum(y.impute[doses==i]))
        tns <- c(tns, sum(doses==i))
      } 
    }
  }else if(impute.method == "TITE"){
    impute.res <-  TITEImpute2(enter.times, dlt.times, current.t, tau, doses, ndose, c(phi/2, 1-phi/2))
    y.raw <- impute.res$y.raw
    y.impute <- impute.res$y.impute
    for (i in 1:ndose){
      tys <- c(tys, sum(y.impute[doses==i]))
      tns <- c(tns, sum(doses==i))
    } 
  }else if(impute.method == "No"){
    assesstime = enter.times+tau;	
    dlt.times[dlt.times==0]= tau+1;
    y <- (dlt.times<=tau)*(assesstime<=current.t)
    for (i in 1:ndose){
      tys <- c(tys, sum(y[doses==i]))
      tns <- c(tns, sum(doses==i))
    }
  }
  
  ### early stop
  add.args <- c(list(y=tys[curDose], n=tns[curDose]), add.args)
  if (overdose.fn(phi, add.args)){
    tover.doses[curDose:ndose] <- 1
  }
  
  if (tover.doses[1] == 1){
    dose <- 0
  }else{
    ### the results for all dose levels
    dose.chg <- make.decision.ACFO.fn(phi, tys, tns, add.args$alp.prior, add.args$bet.prior, tover.doses, curDose, diag=FALSE) - 2
    dose <- dose.chg + curDose
  }
  
  if (simu){
    res <- list(dose=dose, tover.doses=tover.doses)
    return(res)
  }else{
    res <- list(dose=dose)
    return(res)
  }
}


#------------------------------------------------------------------------------------------
# main function to RUN simulation
#------------------------------------------------------------------------------------------

Simu.Fn <- function(phi, ps, tau, cohortsize, ncohort,                                       
                    accrual, tite.dist, accrual.dist, design=1, impute.method="frac",
                    init.dose=1, add.args=list()){
  #args:
  #   phi: Target DLT rate
  #   ps: True dlt rates for each dose level      
  #   tau: maximal observing win size
  #   cohortsize: cohort size
  #   ncohort: num of cohort
  #   accrual: the accrual rate, i.e., the number of patients accrued in tau time 
  #   tite.dist: Time to event distibution, 1-uniform, 2-weibull, 3-log-log
  #   accrual.dist: accrual distribution  0: fixed, 1, uniform, 2: exponential,  
  #   design: the phase I design, 1: fCFO or TITE-CFO, 2: f-aCFO or TITE-aCFO, 3: benchmark-aCFO
  #   impute.method: impute method: frac: fractional design, TITE: TITE way
  #   init.dose: Initial dose level, 1 by default
  #   add.args: additional parameters
  #       alp.prior, bet.prior: prior parameter for TITE
  #       CV: critical value of early stopping for CFO and TITE-BOIN
  #       crmCI.CV: coverage prob of CI for CRM
  #       ps.prior: prior skeleton for CRM
  #       boinEps: The parameter set the boundary for BOIN, by default, boinEps = 0.4
  
  ndose <- length(ps)
  if (is.null(add.args$alp.prior)){
    add.args <- c(add.args, list(alp.prior=phi, bet.prior=1-phi))
  }
  
  earlystop <- 0
  enter.times <- NULL # enter time of each subject
  dlt.times <- NULL # dlt time of each subject
  dlts <- NULL # dlt event for each subject
  doses <- NULL # dose level for each subject
  current.t<- 0
  curDose <- init.dose  #current dose level
  
  tover.doses <- rep(0, ndose)
  
  for (i in 1:ncohort){
    curP <- ps[curDose]    
    
    if (accrual.dist==0){
      delta.times <- rep(0, cohortsize)
    }else if (accrual.dist == 1){
      delta.times <- cumsum(c(0, runif(cohortsize-1, 0,  2*tau/accrual)))
    }else if (accrual.dist == 2){
      delta.times <- cumsum(c(0, rexp(cohortsize-1, rate=accrual/tau)))
    }
    enter.times <- c(enter.times, current.t+delta.times)
    
    # obtain the results of the patients
    obscohort <- gen.tite(tite.dist, cohortsize, curP, alpha=0.5, tau=tau);
    dlt.times <- c(dlt.times, obscohort$t.tox);
    dlts <- c(dlts, obscohort$tox);
    doses <- c(doses, rep(curDose, cohortsize));
    # print(dlt.times)
    # print(dlts)
    # print(curP)
    
    # Move to next cohort 
    if (i != ncohort){
      if (accrual.dist==0){
        delta.time <- tau*cohortsize/accrual
      }else if (accrual.dist == 1){
        delta.time <- runif(1, 0, 2*tau/accrual)
      }else if (accrual.dist == 2){
        delta.time <- rexp(1, rate=accrual/tau)
      }
    }else{
      delta.time <- tau
    }
    current.t<- enter.times[length(enter.times)] + delta.time;
    
    if (design == 3){
      current.t <- enter.times[length(enter.times)] + tau
      res <- aCFOlateonset.next.dose(curDose, phi, tau, impute.method, enter.times, 
                                     dlt.times, current.t, doses, tover.doses, TRUE, add.args)
      tover.doses <- res$tover.doses
      current.t <- current.t + delta.time
    }else if (design == 1){
      #fCFO or TITE-CFO 
      res <- CFOlateonset.next.dose(curDose, phi, tau, impute.method, enter.times, 
                                    dlt.times, current.t, doses, tover.doses, TRUE, add.args)
      tover.doses <- res$tover.doses
    }else if(design == 2){
      #f-aCFO or TITE-aCFO
      res <- aCFOlateonset.next.dose(curDose, phi, tau, impute.method, enter.times, 
                                    dlt.times, current.t, doses, tover.doses, TRUE, add.args)
      tover.doses <- res$tover.doses
    }else{
      res <- list()
      res$dose <- curDose
    }
    
    if (res$dose==0){
      earlystop <- 1
      break()
    }else{
      curDose <- res$dose
    }
  }
  
  tns <- NULL
  tys <- NULL
  assess.t <- enter.times + tau
  y.raw <- (dlt.times!=0)*1
  for (j in 1:ndose){
    tns <- c(tns, sum(doses==j))
    tys <- c(tys, sum(y.raw[doses==j]))
  }
  if (earlystop==0){
    MTD <- select.mtd(phi, tns, tys)$MTD
  }else{
    MTD <- 99
  }
  res <- list(MTD=MTD, dose.ns=tns, DLT.ns=tys, ps=ps, target=phi, total.time=assess.t[length(assess.t)])
  res
  }