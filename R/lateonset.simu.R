#' Find the maximum tolerated dose (MTD) for a single trial using the CFO-type and aCFO-type designs with late-onset toxicities.
#'
#' Use this function to find the maximum tolerated dose (MTD) for the CFO-type and aCFO-type designs with late-onset toxicities, 
#' specifically, including Time-to-event CFO (TITE-CFO) design, fractional CFO (fCFO), benchmark CFO, TITE-aCFO design, 
#' the f-aCFO design and benchmark aCFO design.
#'
#' @usage lateonset.simu(phi, p.true, tau, cohortsize, ncohort, accrual, tite.dist, accrual.dist, 
#'        design, init.dose=1, add.args=list(alp.prior=phi, bet.prior=1-phi), seed=NULL)
#'
#' @param phi the target DLT rate.
#' @param p.true the true DLT rates under the different dose levels.
#' @param tau maximal assessment window size.
#' @param cohortsize the sample size in each cohort.
#' @param ncohort the total number of cohorts.
#' @param accrual the accrual rate, i.e., the number of patients accrued in tau time.
#' @param tite.dist the distribution of the time to DLT events. \code{tite.dist=1} corresponds to a uniform distribution, 
#'                  \code{tite.dist=2} corresponds to a Weibull distribution, and \code{tite.dist=3} corresponds to a 
#'                  log-logistic distribution.
#' @param accrual.dist the distribution of the arrival times of patients. When \code{tite.dist=1}, it corresponds to all 
#'                     patients in each cohort arriving simultaneously at a given accrual rate. When \code{tite.dist=2}, 
#'                     it corresponds to a uniform distribution, and when \code{tite.dist=3}, it corresponds to an 
#'                     exponential distribution.
#' @param design option for selecting different designs, including TITE-CFO, TITE-aCFO, fCFO, f-aCFO, bCFO, and
#'               b-aCFO. Specifically, bCFO refers to the benchmark CFO, and b-aCFO denotes the benchmark aCFO.
#' @param init.dose the dose level assigned to the first cohort. The default value \code{init.level} is 1.
#' @param add.args additional parameters, usually set as list(alp.prior=phi, bet.prior=1-phi) by default. \code{alp.prior} 
#'                 and \code{bet.prior} represent the parameters of the prior distribution for the true DLT rate at 
#'                 any dose level. This prior distribution is specified as Beta( \code{alpha.prior}, \code{beta.prior}).
#' @param seed an integer to set as the seed of the random number generator for reproducible results, the default is set to NULL.
#' @details The \code{lateonset.simu()} function is developed for determining the MTD in a single trial using CFO-type and 
#'          aCFO-type designs with late-onset toxicities. Specifically, it includes the TITE-CFO design, fCFO design, benchmark 
#'          CFO design, TITE-aCFO design, and f-aCFO design. The \code{design} parameter allows for selecting different designs. 
#'          Please note that the results returned by \code{lateonset.simu()} represent a single trial outcome. In practical 
#'          applications, multiple simulations are required to ascertain the MTD, and this can be achieved using the 
#'          \code{CFO.oc()} function.
#'
#' @return The \code{lateonset.simu()} function returns a list object comprising the following components: the target DLT 
#'         rate ($target), the actual DLT rates under different dose levels ($p.true), the selected MTD ($MTD),
#'         the list that includes the dose level assigned to each cohort($dose.list), the total number of DLTs and 
#'         patients for all dose levels ($DLT.ns and $dose.ns), and the duration of the trial in months ($total.time).
#'         
#' @author Jialu Fang
#' 
#' @references Jin, H., & Yin, G. (2022). CFO: Calibration-free odds design for phase I/II clinical trials. 
#'             \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066. \cr
#'             Jin, H., & Yin, G. (2023). Time‐to‐event calibration‐free odds design: A robust efficient design for 
#'             phase I trials with late‐onset outcomes. \emph{Pharmaceutical Statistics}. \cr
#'             Yin, G., Zheng, S., & Xu, J. (2013). Fractional dose-finding methods with late-onset toxicity in 
#'             phase I clinical trials. \emph{Journal of Biopharmaceutical Statistics}, 23(4), 856-870.
#' @export
#'
#' @examples
#' phi <- 0.2; ncohort <- 12; cohortsize <- 3
#' p.true <-c(0.01, 0.05, 0.10, 0.14, 0.20, 0.26, 0.34)
#'  tau <- 3; accrual <- 6; tite.dist <- 2; accrual.dist <- 1
#' ## find the MTD for a single TITE-CFO trial
#' lateonset.simu (phi, p.true, tau, cohortsize, ncohort, accrual, tite.dist, accrual.dist, 
#'                 design='TITE-CFO', init.dose=1, add.args=list(alp.prior=phi, bet.prior=1-phi))
#' ## find the MTD for a single TITE-aCFO trial
#' lateonset.simu (phi, p.true, tau, cohortsize, ncohort, accrual, tite.dist, accrual.dist, 
#'                 design='TITE-aCFO', init.dose=1, add.args=list(alp.prior=phi, bet.prior=1-phi))
#' ## find the MTD for a single fCFO trial
#' lateonset.simu (phi, p.true, tau, cohortsize, ncohort, accrual, tite.dist, accrual.dist, 
#'                 design='fCFO', init.dose=1, add.args=list(alp.prior=phi, bet.prior=1-phi))
#' ## find the MTD for a single f-aCFO trial
#' lateonset.simu (phi, p.true, tau, cohortsize, ncohort, accrual, tite.dist, accrual.dist, 
#'                 design='f-aCFO', init.dose=1, add.args=list(alp.prior=phi, bet.prior=1-phi))
#' ## find the MTD for a single benchmark CFO trial
#' lateonset.simu (phi, p.true, tau, cohortsize, ncohort, accrual, tite.dist, accrual.dist, 
#'                 design='bCFO', init.dose=1, add.args=list(alp.prior=phi, bet.prior=1-phi))
#' ## find the MTD for a single benchmark aCFO trial
#' lateonset.simu (phi, p.true, tau, cohortsize, ncohort, accrual, tite.dist, accrual.dist, 
#'                 design='b-aCFO', init.dose=1, add.args=list(alp.prior=phi, bet.prior=1-phi))
lateonset.simu <- function(phi, p.true, tau, cohortsize, ncohort, accrual, tite.dist, accrual.dist, 
                    design, init.dose=1, add.args=list(alp.prior=phi, bet.prior=1-phi), seed=NULL){
  
  ###############################################################################
  ###############define the functions used for main function#####################
  ###############################################################################
  
  # The function is to obtain the DLT results (with TITE) for each subject
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
  
  ###############################################################################
  ############################MAIN DUNCTION######################################
  ############################################################################### 
  
  if (design == 'TITE-CFO'){accumulation = FALSE; impute.method = "TITE"
  }else if (design == 'fCFO'){accumulation = FALSE; impute.method = "frac"
  }else if (design == 'bCFO'){accumulation = FALSE; impute.method = "No"
  }else if (design == 'TITE-aCFO'){accumulation = TRUE; impute.method = "TITE"
  }else if (design == 'f-aCFO'){accumulation = TRUE; impute.method = "frac"
  }else if (design == 'b-aCFO'){accumulation = TRUE; impute.method = "No"}
  
  set.seed(seed)
  ndose <- length(p.true)
  doselist <- rep(0, ncohort)
  
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
    curP <- p.true[curDose]
    doselist[i] <- curDose
    
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
    
    if (design == 'bCFO' || design == 'b-aCFO'){
      current.t <- enter.times[length(enter.times)] + tau
      lateonset.simu 
      res <- lateonset.next(curDose, phi, tau, impute.method, enter.times, dlt.times, current.t, accumulation,
                            doses, tover.doses, TRUE, add.args)
      tover.doses <- res$tover.doses
      current.t <- current.t + delta.time
    }else{
      res <- lateonset.next(curDose, phi, tau, impute.method, enter.times, dlt.times, current.t, accumulation,
                                doses, tover.doses, TRUE, add.args)
      tover.doses <- res$tover.doses
    }
    
    if (res$curDose==0){
      earlystop <- 1
      break()
    }else{
      curDose <- res$curDose
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
  out <- list(MTD=MTD, dose.list=doselist, dose.ns=tns, DLT.ns=tys, p.true=p.true, 
              target=phi, total.time=assess.t[length(assess.t)])
  class(out) <- "cfo"
  return(out)
}
