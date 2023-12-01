#' Determination of the dose level for next cohort in the CFO-type and aCFO-type designs with late-onset toxicities
#' 
#' Propose the next dose level in the CFO-type and aCFO-type designs with late-onset toxicities, specifically, including 
#' Time-to-event CFO (TITE-CFO) design, fractional CFO (fCFO) design, benchmark CFO design, TITE-aCFO design, f-aCFO 
#' design and benchmark aCFO design.
#' 
#' @usage lateonset.next(curDose, phi, tau, impute.method, enter.times, dlt.times,
#'        current.t, accumulation, doses, tover.doses=c(), simu=FALSE,
#'        add.args=list(alp.prior=phi, bet.prior=1-phi), seed=NULL)
#'
#' @param curDose the current dose level.
#' @param phi the target DLT rate.
#' @param tau maximal assessment window size
#' @param impute.method the imputing method for handling pending DLT data. \code{impute.method = 'frac'} corresponds to 
#'                      fractional framework. \code{impute.method = 'TITE'} corresponds to time-to-event framework.
#'                      \code{impute.method = 'No'} implies no use of any imputing method, corresponding to the 
#'                      benchmark CFO and benchmark aCFO designs.
#' @param enter.times the time at which each subject existing in the trial enters the trial.
#' @param dlt.times the time to DLT for each subject existing in the trial.  If no DLT occurs for a certain subject, 
#'                  \code{dlt.times} is set to 0.
#' @param current.t the current time.
#' @param accumulation set \code{accumulation=FALSE} to conduct the CFO-type design; set \code{accumulation=TRUE} to 
#'                     conduct the aCFO-type design.
#' @param doses the dose level for each subject existing in the trial.
#' @param tover.doses whether the dose level (from the first to last dose level) is over-toxic or not. 
#'                    The value is set as 1 if the dose level is overly toxicity; otherwise, it is set to 0.
#'                    It should be predetermined. If not predetermined, it is set to NaN and is assigned values of 0 
#'                    for each dose level.
#' @param simu whether simulation or not, if \code{simu=TRUE}, \code{lateonset.next()} also return \code{tover.doses}.
#' @param add.args additional parameters, usually set as list(alp.prior=phi, bet.prior=1-phi) by default. \code{alp.prior} 
#'                 and \code{bet.prior} represent the parameters of the prior distribution for the true DLT rate at 
#'                 any dose level. This prior distribution is specified as Beta( \code{alpha.prior}, \code{beta.prior}).
#' @param seed an integer to set as the seed of the random number generator for reproducible results.
#'
#' @details Late-onset outcomes commonly occur in phase I trials involving targeted agents or immunotherapies. As a 
#'          result, the TITE framework and fractional framework serve as two imputation methods to handle pending data 
#'          related to late-onset outcomes. This approach extends the original designs to integrate time information 
#'          for delayed outcomes, leading to the development of TITE-CFO, fCFO, TITE-aCFO, and f-aCFO designs. \cr
#'          In the TITE framework context, an assumption about the distribution of time to DLT must be pre-specified, 
#'          whereas the fractional framework does not require justification for a specific distribution of the time to 
#'          DLT. Consequently, fCFO and f-aCFO adapt to a more diverse range of scenarios.\cr
#'          The function \code{lateonset.next()} also provides the option to set \code{impute.method = "No"} to execute 
#'          the benchmark CFO and aCFO design. These two methods await complete observation of toxicity outcomes for 
#'          the previous cohorts before determining the next dose assignment. This enhances precision but comes at the 
#'          expense of a prolonged trial duration.
#'          
#' @note   The \code{tover.doses()} should be predetermined. It can be calculated based on the number of patients enrolled 
#'          ($tns) and the number of DLTs observed ($tys), following the overdose control rule. If not predetermined, it is 
#'          set to NaN and is assigned default values of 0 for each dose level.
#' 
#' @return The \code{lateonset.next()} function returns the recommended dose level for treating the next cohort of 
#'         patients ($curDose). if \code{simu=TRUE}, \code{lateonset.next()} also return a vector indicating whether the 
#'         dose level (from the first to last dose level) is over-toxic or not ($tover.doses).
#' 
#' @author Jialu Fang 
#' 
#' @references Jin, H., & Yin, G. (2022). CFO: Calibration-free odds design for phase I/II clinical trials. 
#'             \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066. \cr
#'             Jin, H., & Yin, G. (2023). Time‐to‐event calibration‐free odds design: A robust efficient design for 
#'             phase I trials with late‐onset outcomes. \emph{Pharmaceutical Statistics}. \cr
#'             Yin, G., Zheng, S., & Xu, J. (2013). Fractional dose-finding methods with late-onset toxicity in 
#'             phase I clinical trials. \emph{Journal of Biopharmaceutical Statistics}, 23(4), 856-870.
#' @import survival
#' @importFrom utils tail
#' @export
#'
#' @examples
#' ## Given the parameters for the function, the unit for time-related parameters is in months. 
#' ## The parameter generation follows the following guidelines: 
#' ## 1) the assessment window ([0, $tau]) covers three months
#' ## 2) the accrual rate is set at two patients per month. 
#' ## 3) the arrival times of patients are distributed uniformly.
#' ## 4) the time to DLT events is simulated using a Weibull distribution, with 50% of these events 
#' ##    occurring in the first half of the assessment window.
#' phi<-0.2
#' add.args=list(alp.prior=phi, bet.prior=1-phi)
#' enter.times<-c(0,0.082,0.343,0.554,1.18,1.88,2.15,2.68,2.74,3.30,3.99,
#'                4.52,5.51,5.93,6.18,6.68,7.05,7.92)
#' dlt.times<-c(0,0,0,0,0,0,0,0,0,0.934,0,0,0,0,0,1.5,0.6962,0)
#' current.t<-8.134413
#' doses<-c(1,1,1,2,2,2,3,3,3,4,4,4,3,3,3,4,4,4)
#' tover.doses<-c(0,0,0,0,0,0,0)
#' ## determine the dose level for the next cohort using the TITE-CFO design
#' lateonset.next(curDose=4, phi, tau=3, impute.method="TITE", enter.times, dlt.times, current.t, 
#'                accumulation = FALSE, doses, tover.doses, simu=FALSE, add.args)
#' ## determine the dose level for the next cohort using the TITE-aCFO design
#' lateonset.next(curDose=4, phi, tau=3, impute.method="TITE", enter.times, dlt.times, current.t, 
#'                accumulation = TRUE, doses, tover.doses, simu=FALSE, add.args)
#' ## determine the dose level for the next cohort using the f-CFO design
#' lateonset.next(curDose=4, phi, tau=3, impute.method="frac", enter.times, dlt.times, current.t, 
#'                accumulation = FALSE, doses, tover.doses, simu=FALSE, add.args)
#' ## determine the dose level for the next cohort using the f-aCFO design
#' lateonset.next(curDose=4, phi, tau=3, impute.method="frac", enter.times, dlt.times, current.t, 
#'                accumulation = TRUE, doses, tover.doses, simu=FALSE, add.args)
#' ## determine the dose level for the next cohort using the benchmark CFO design
#' lateonset.next(curDose=4, phi, tau=3, impute.method="No", enter.times, dlt.times, current.t, 
#'                accumulation = FALSE, doses, tover.doses, simu=FALSE, add.args)
#' ## determine the dose level for the next cohort using the benchmark aCFO design
#' lateonset.next(curDose=4, phi, tau=3, impute.method="No", enter.times, dlt.times, current.t, 
#'                accumulation = TRUE, doses, tover.doses, simu=FALSE, add.args)
#' 
lateonset.next <- function(curDose, phi, tau, impute.method, enter.times, dlt.times, current.t, 
                           accumulation, doses, tover.doses=c(), simu=FALSE, 
                           add.args=list(alp.prior=phi, bet.prior=1-phi), seed=NULL){
  ###############################################################################
  ###############define the functions used for main function#####################
  ###############################################################################
  
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
      kmfit = survival::survfit(survival::Surv(otime,yo)~1)	
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
    
    
    res <- list(y.impute=ym, y.raw=yo, obsIdxs=obsIdxs, No.impute=No.impute)
    res
  }
  
  TITEImpute.one <- function(followup.times, tau, y, n, prior.paras){
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
  
  TITEImpute <- function(enter.times, dlt.times, current.time, tau, dose.levels, ndose, prior.paras){
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
        ym.part <- TITEImpute.one(followup.times[kpidxs], tau, y, n, prior.paras)
        ym[kpidxs] <- ym.part
      }
    }
    res <- list(y.impute=ym, y.raw=yo, obsIdxs=obsIdxs)
    res
  }
  
  # posterior probability of pj >= phi given data
  post.prob.fn <- function(phi, y, n, alp.prior=0.1, bet.prior=0.9){
    alp <- alp.prior + y 
    bet <- bet.prior + n - y
    1 - pbeta(phi, alp, bet)
  }
  
  overdose.fn <- function(phi, add.args=list()){
    y <- add.args$y
    n <- add.args$n
    alp.prior <- add.args$alp.prior
    bet.prior <- add.args$bet.prior
    pp <- post.prob.fn(phi, y, n, alp.prior, bet.prior)
    # print(data.frame("prob of overdose" = pp))
    if ((pp >= 0.95) & (add.args$n>=3)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
  
  ###############################################################################
  ############################MAIN DUNCTION######################################
  ############################################################################### 
  set.seed(seed)
  ndose <- length(tover.doses)
  if (is.null(tover.doses)){
    tover.doses <- rep(0,ndose)
    warning("tover.doses is set to NaN and is assigned default values of ", 
            paste(tover.doses, collapse = ","))
  }
  if (is.null(add.args$alp.prior)){
    add.args <- c(add.args, list(alp.prior=phi, bet.prior=1-phi))
  }
  tys = NULL
  tns = NULL
  ## Obtain effective results
  if (impute.method == "frac"){
    impute.res <- fracImpute(enter.times, dlt.times, current.t, tau)
    y.raw <- impute.res$y.raw
    y.impute <- impute.res$y.impute
    if (impute.res$No.impute){
      if (accumulation == FALSE){
        cn <- sum(doses==curDose)
        cy <- sum(y.raw[doses==curDose])
      }else{
        for (i in 1:ndose){
          tys <- c(tys, sum(y.raw[doses==i]))
          tns <- c(tns, sum(doses==i))
        } 
      }
    }
    else{
      if (accumulation == FALSE){
        cn <- sum(doses==curDose)
        cy <- sum(y.impute[doses==curDose])
      }else{
        for (i in 1:ndose){
          tys <- c(tys, sum(y.impute[doses==i]))
          tns <- c(tns, sum(doses==i))
        }
      }
    }
  }else if(impute.method == "TITE"){
    impute.res <-  TITEImpute(enter.times, dlt.times, current.t, tau, doses, ndose, c(phi/2, 1-phi/2))
    y.raw <- impute.res$y.raw
    y.impute <- impute.res$y.impute
    if (accumulation == FALSE){
      cy <- sum(y.impute[doses==curDose])
      cn <- sum(doses==curDose)
    }else{
      for (i in 1:ndose){
        tys <- c(tys, sum(y.impute[doses==i]))
        tns <- c(tns, sum(doses==i))
      } 
    }
  }else if(impute.method == "No"){
    assesstime = enter.times+tau;	
    dlt.times[dlt.times==0]= tau+1;
    y.impute <- (dlt.times<=tau)*(assesstime<=current.t)
    if (accumulation == FALSE){
      cy <- sum(y.impute[doses==curDose])
      cn <- sum(doses==curDose)
    }else{
      for (i in 1:ndose){
        tys <- c(tys, sum(y.impute[doses==i]))
        tns <- c(tns, sum(doses==i))
      }
    }
  }
  
  ### early stop
  if (accumulation == FALSE) add.args <- c(list(y=cy, n=cn), add.args) else add.args <- c(list(y=tys[curDose], n=tns[curDose]), add.args)
  if (overdose.fn(phi, add.args)){
    tover.doses[curDose:ndose] <- 1
  }
  if (tover.doses[1] == 1){
    dose <- 0
  }else{
    if (accumulation == FALSE){
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
      curDose <- CFO.next(phi, cys, cns, cover.doses, curDose, add.args)$curDose
    }else{
      curDose <- aCFO.next (phi, tys, tns, tover.doses, curDose, add.args)$curDose
    }
  }
  
  if (simu){
    res <- list(curDose=curDose, tover.doses=tover.doses)
    return(res)
  }else{
    res <- list(curDose=curDose)
    return(res)
  }
}

