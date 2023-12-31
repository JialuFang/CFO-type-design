#' Find the maximum tolerated dose (MTD) for a single Calibration-Free Odds (CFO) or accumulative CFO (aCFO) trial.
#' 
#' Use this function to find the maximum tolerated dose (MTD) for a single Calibration-Free Odds (CFO) or accumulative CFO (aCFO) trial.
#'
#' @usage CFO.simu(phi, p.true, ncohort, init.level=1, cohortsize=3, 
#'        add.args=list(alp.prior=phi, bet.prior=1-phi), accumulation = FALSE, seed=NULL)
#'
#' @param phi the target DLT rate.
#' @param p.true the true DLT rates under the different dose levels.
#' @param ncohort the total number of cohorts.
#' @param init.level the dose level assigned to the first cohort. The default value \code{init.level} is 1.
#' @param cohortsize the sample size in each cohort. The default value \code{cohortsize} is 3.
#' @param accumulation set \code{accumulation=FALSE} to conduct the CFO design; set \code{accumulation=TRUE} to conduct the 
#'                     aCFO design.
#' @param add.args additional parameters, usually set as list(alp.prior=phi, bet.prior=1-phi) by default. \code{alp.prior} 
#'                 and \code{bet.prior} represent the parameters of the prior distribution for the true DLT rate at 
#'                 any dose level. This prior distribution is specified as Beta( \code{alpha.prior}, \code{beta.prior}).
#' @param seed an integer to set as the seed of the random number generator for reproducible results, the default is set to NULL.
#'                            
#' @details The \code{CFO.simu()} function is designed to determine the Maximum Tolerated Dose (MTD) for a single CFO or aCFO 
#'          trial. If \code{accumulation = FALSE}, this trial corresponds to the CFO design. If \code{accumulation = TRUE}, it
#'          corresponds to the aCFO design. \cr
#'          Given the toxicity outcomes from previous cohorts, each cohort is sequentially assigned to the most suitable dose 
#'          level based on the CFO or aCFO decision rule. Early stopping criteria are incorporated into the CFO or aCFO design 
#'          to ensure patient safety and benefit. If there is substantial evidence indicating that the current dose level 
#'          exhibits excessive toxicity (\eqn{\Pr(p_C > \phi|x_C, m_C \geq 3) > 0.95}), we exclude the current dose level as 
#'          well as higher dose levels from the trial. Upon the predefined maximum sample size is reached or the lowest dose 
#'          level is overly toxicity, the experiment is concluded, and the MTD is determined using isotonic regression.
#'
#' @return The \code{CFO.simu} function returns a list object comprising the following components: 
#'         \itemize{
#'         \item{target: }{the target DLT rate.}
#'         \item{p.true: }{the actual DLT rates under different dose levels.}
#'         \item{MTD: }{the selected MTD.}
#'         \item{dose.list: }{the list that includes the dose level assigned to each cohort.}
#'         \item{DLT.ns: }{the total number of DLTs for all dose levels.}
#'         \item{dose.ns: }{the total number of patients for all dose levels}
#'         \item{tover.doses: }{the over-toxicity status for all dose levels (1 for yes).}
#'         \item{DLT.list: }{The DLT outcome observed at each subject.}
#'         }
#' 
#' @author Jialu Fang
#' 
#' @references Jin, H., & Yin, G. (2022). CFO: Calibration-free odds design for phase I/II clinical trials. 
#'             \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066.
#'
#' @examples
#' phi <- 0.2; ncohort <- 12; cohortsize <- 3
#' p.true <-c(0.01, 0.05, 0.10, 0.14, 0.20, 0.26, 0.34)
#' ## find the MTD for a single CFO trial
#' CFO.simu(phi, p.true, ncohort, init.level=1, cohortsize=3,
#'          add.args=list(alp.prior=phi, bet.prior=1-phi), accumulation = FALSE)
#' ## find the MTD for a single aCFO trial
#' CFO.simu(phi, p.true, ncohort, init.level=1, cohortsize=3,
#'          add.args=list(alp.prior=phi, bet.prior=1-phi), accumulation = TRUE)
#' @import BOIN
#' @export
CFO.simu <- function(phi, p.true, ncohort, init.level=1, cohortsize=3,
                      add.args=list(alp.prior=phi, bet.prior=1-phi), accumulation = FALSE, seed=NULL){
  ###############################################################################
  ###############define the functions used for main function#####################
  ###############################################################################
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
  earlystop <- 0
  ndose <- length(p.true)
  doselist <- rep(0, ncohort)
  curDose <- init.level
  
  tys <- rep(0, ndose) # number of responses for different doses.
  tns <- rep(0, ndose) # number of subject for different doses.
  tover.doses <- rep(0, ndose) # Whether each dose is overdosed or not, 1 yes
  DLTlist <- c()
  
  for (i in 1:ncohort){
    pc <- p.true[curDose]
    doselist[i] <- curDose
    
    # sample from current dose
    cres <- rbinom(cohortsize, 1, pc)
    DLTlist <- c(DLTlist, cres)
    
    # update results
    tys[curDose] <- tys[curDose] + sum(cres)
    tns[curDose] <- tns[curDose] + cohortsize
    
    cy <- tys[curDose]
    cn <- tns[curDose]
    
    add.args <- c(list(y=cy, n=cn, tys=tys, tns=tns, curDose=curDose), add.args)
    
    if (overdose.fn(phi, add.args)){
      tover.doses[curDose:ndose] <- 1
    }
    
    if (tover.doses[1] == 1){
      earlystop <- 1
      break()
    }
    
    if (accumulation == FALSE){
      # the results for current 3 dose levels
      if (curDose!=1){
        cys <- tys[(curDose-1):(curDose+1)]
        cns <- tns[(curDose-1):(curDose+1)]
        cover.doses <- tover.doses[(curDose-1):(curDose+1)]
        #cover.doses <- c(0, 0, 0) # No elimination rule
      }else{
        cys <- c(NA, tys[1:(curDose+1)])
        cns <- c(NA, tns[1:(curDose+1)])
        cover.doses <- c(NA, tover.doses[1:(curDose+1)])
        #cover.doses <- c(NA, 0, 0) # No elimination rule
      }
      curDose <- CFO.next(phi, cys, cns, cover.doses, curDose, add.args)$curDose
    } else if (accumulation == TRUE){
      curDose <- aCFO.next (phi, tys, tns, tover.doses, curDose, add.args)$curDose
    }
  }
  
  
  if (earlystop==0){
    MTD <- select.mtd(phi, tns, tys)$MTD
  }else{
    MTD <- 99
  }
  out<-list(target=phi, p.true=p.true, MTD=MTD, dose.list=doselist, dose.ns=tns,
            DLT.ns=tys, tover.doses=tover.doses, DLT.list=DLTlist)
  class(out) <- "cfo"
  return(out)
}

