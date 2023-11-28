#' 
#' Determination of the dose level for next cohort
#' 
#' Use the function to determine the dose movement based on the toxicity outcomes of the enrolled cohorts.
#'
#' @usage CFO.next(phi, cys, cns, cover.doses=c(0,0,0), curDose,
#'        add.args=list(alp.prior=phi,bet.prior=1-phi))
#'
#' @param phi the target DLT rate
#' @param cys the current number of DLTs observed in patients for the left, current, and right dose levels.
#' @param cns the current number of patients for the left, current, and right dose levels.
#' @param alp.prior,bet.prior the parameters of the prior distribution for the true DLT rate at any dose level.
#'                            This prior distribution is set to Beta( \code{alpha.prior}, \code{beta.prior}). 
#'                            The default value is \code{phi} and \code{1-phi}.
#' @param curDose the current dose level.
#' @param cover.doses whether the dose level (left, current and right) is over-toxic or not. 
#'                    The value is set as 1 if the dose level is overly toxicity; otherwise, it is set to 0.
#'
#' @details The CFO design determines the dose level for the next cohort by assessing evidence from the current 
#'          dose level and its adjacent levels. This evaluation is based on odds ratios denoted as \eqn{O_k}, where 
#'          k = L, C, R represents left, current, and right dose levels. Additionally, we define \eqn{\overline{O}_k = 1/O_k}. 
#'          The ratio \eqn{O_C / \overline{O}_{L}} indicates the inclination for de-escalation, while \eqn{\overline{O}_C / O_R} 
#'          quantifies the tendency for escalation. Threshold values \eqn{\gamma_L} and \eqn{\gamma_R} are chosen to 
#'          minimize the probability of making incorrect decisions.The decision process is summarized in Table 1
#'          of Jin and Yin (2022).
#'          An overdose control rule is implemented to ensure patient safety. If the data suggest excessive 
#'          toxicity at the current dose level, we exclude that level and those higher levels. Two scenarios 
#'          lead to a decision on one side only: when the current dose is at the boundary (the first or last dose level) 
#'          or when higher dose levels have been eliminated.
#'          
#' @note    The \code{cover.doses()} should be predetermined. It can be calculated based on the number of patients enrolled 
#'          ($cns) and the number of DLTs observed ($cys), following the overdose control rule. If not predetermined, it is 
#'          assigned default values of 0 for each dose level.
#'          
#' @return The \code{CFO.next()} function returns a list object comprising the following elements: the target DLT 
#'         rate ($target), the current counts of DLTs and patients for the left, current, and right dose levels ($cys and $cns), 
#'         the decision for whether to move to the left or right dose level for the next cohort ($decision), and the current 
#'         dose after movement ($curDose).
#'         
#' @author Jialu Fang and Wenliang Wang
#' 
#' @references Jin, H., & Yin, G. (2022). CFO: Calibration-free odds design for phase I/II clinical trials. 
#'             \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066.
#' 
#' @examples
#' ## determine the dose level for the next cohort of new patients
#' cys <- c(0,1,0); cns <- c(3,6,0)
#' CFO.next(phi=0.2, cys=cys, cns=cns, cover.doses=c(0,0,0), curDose=3,
#'          add.args=list(alp.prior=0.2, bet.prior=0.8))
#' 
#' @import stats
#' @export
CFO.next <- function(phi, cys, cns, cover.doses=c(0,0,0), curDose,
                     add.args=list(alp.prior=phi, bet.prior=1-phi)){
  ###############################################################################
  ###############define the functions used for main function#####################
  ###############################################################################
  # posterior probability of pj >= phi given data
  post.prob.fn <- function(phi, y, n, alp.prior=0.1, bet.prior=0.1){
    alp <- alp.prior + y 
    bet <- bet.prior + n - y
    1 - pbeta(phi, alp, bet)
  }
  
  prob.int <- function(phi, y1, n1, y2, n2, alp.prior, bet.prior){
    alp1 <- alp.prior + y1
    alp2 <- alp.prior + y2
    bet1 <- bet.prior + n1 - y1
    bet2 <- bet.prior + n2 - y2
    fn.min <- function(x){
      dbeta(x, alp1, bet1)*(1-pbeta(x, alp2, bet2)) 
    }
    fn.max <- function(x){
      pbeta(x, alp1, bet1)*dbeta(x, alp2, bet2)
    }
    const.min <- integrate(fn.min, lower=0, upper=1)$value
    const.max <- integrate(fn.max, lower=0, upper=1)$value
    p1 <- integrate(fn.min, lower=0, upper=phi)$value/const.min
    p2 <- integrate(fn.max, lower=0, upper=phi)$value/const.max
    
    list(p1=p1, p2=p2)
  }
  
  
  OR.values <- function(phi, y1, n1, y2, n2, alp.prior, bet.prior, type){
    ps <- prob.int(phi, y1, n1, y2, n2, alp.prior, bet.prior)
    if (type=="L"){
      pC <- 1 - ps$p2
      pL <- 1 - ps$p1
      oddsC <- pC/(1-pC)
      oddsL <- pL/(1-pL)
      OR <- oddsC*oddsL
      
    }else if (type=="R"){
      pC <- 1 - ps$p1
      pR <- 1 - ps$p2
      oddsC <- pC/(1-pC)
      oddsR <- pR/(1-pR)
      OR <- (1/oddsC)/oddsR
    }
    return(OR)
  }
  
  All.OR.table <- function(phi, n1, n2, type, alp.prior, bet.prior){
    ret.mat <- matrix(rep(0, (n1+1)*(n2+1)), nrow=n1+1)
    for (y1cur in 0:n1){
      for (y2cur in 0:n2){
        ret.mat[y1cur+1, y2cur+1] <- OR.values(phi, y1cur, n1, y2cur, n2, alp.prior, bet.prior, type)
      }
    }
    ret.mat
  }
  
  # compute the marginal prob when lower < phiL/phiC/phiR < upper
  # i.e., Pr(Y=y|lower<phi<upper)
  margin.phi <- function(y, n, lower, upper){
    C <- 1/(upper-lower)
    fn <- function(phi) {
      dbinom(y, n, phi)*C
    }
    integrate(fn, lower=lower, upper=upper)$value
  }
  
  # Obtain the table of marginal distribution of (y1, y2) 
  # after intergrate out (phi1, phi2)
  # under H0 and H1
  # H0: phi1=phi, phi < phi2 < 2phi
  # H1: phi2=phi, 0   < phi1 < phi
  margin.ys.table <- function(n1, n2, phi, hyperthesis){
    if (hyperthesis=="H0"){
      p.y1s <- dbinom(0:n1, n1, phi)
      p.y2s <- sapply(0:n2, margin.phi, n=n2, lower=phi, upper=2*phi)
    }else if (hyperthesis=="H1"){
      p.y1s <- sapply(0:n1, margin.phi, n=n1, lower=0, upper=phi)
      p.y2s <- dbinom(0:n2, n2, phi)
    }
    p.y1s.mat <- matrix(rep(p.y1s, n2+1), nrow=n1+1)
    p.y2s.mat <- matrix(rep(p.y2s, n1+1), nrow=n1+1, byrow=TRUE)
    margin.ys <- p.y1s.mat * p.y2s.mat
    margin.ys
  }
  
  
  # Obtain the optimal gamma for the hypothesis test
  optim.gamma.fn <- function(n1, n2, phi, type, alp.prior, bet.prior){
    OR.table <- All.OR.table(phi, n1, n2, type, alp.prior, bet.prior) 
    ys.table.H0 <- margin.ys.table(n1, n2, phi, "H0")
    ys.table.H1 <- margin.ys.table(n1, n2, phi, "H1")
    
    argidx <- order(OR.table)
    sort.OR.table <- OR.table[argidx]
    sort.ys.table.H0 <- ys.table.H0[argidx]
    sort.ys.table.H1 <- ys.table.H1[argidx]
    n.tol <- length(sort.OR.table)
    
    if (type=="L"){
      errs <- rep(0, n.tol-1)
      for (i in 1:(n.tol-1)){
        err1 <- sum(sort.ys.table.H0[1:i])
        err2 <- sum(sort.ys.table.H1[(i+1):n.tol])
        err <- err1 + err2
        errs[i] <- err
      }
      min.err <- min(errs)
      if (min.err > 1){
        gam <- 0
        min.err <- 1
      }else {
        minidx <- which.min(errs)
        gam <- sort.OR.table[minidx]
      }
    }else if (type=='R'){
      errs <- rep(0, n.tol-1)
      for (i in 1:(n.tol-1)){
        err1 <- sum(sort.ys.table.H1[1:i])
        err2 <- sum(sort.ys.table.H0[(i+1):n.tol])
        err <- err1 + err2
        errs[i] <- err
      }
      min.err <- min(errs)
      if (min.err > 1){
        gam <- 0
        min.err <- 1
      }else {
        minidx <- which.min(errs)
        gam <- sort.OR.table[minidx]
      }
      
    }
    list(gamma=gam, min.err=min.err)
  }
  
  ###############################################################################
  ############################MAIN DUNCTION######################################
  ###############################################################################
  if (is.null(add.args$alp.prior)){
    add.args <- c(add.args, list(alp.prior=phi, bet.prior=1-phi))
  }
  alp.prior <- add.args$alp.prior
  bet.prior <- add.args$bet.prior
  
  if (cover.doses[2] == 1){
    index <- -1
    decision <- "de-escalation"
  }
  else{
    if (is.na(cys[1]) & (cover.doses[3]==1)){
      index <- 0
      decision <- "stay"
    }
    else  if (is.na(cys[1]) & (!(cover.doses[3]==1))){
      gam2 <- optim.gamma.fn(cns[2], cns[3], phi, "R", alp.prior, bet.prior)$gamma 
      OR.v2 <- OR.values(phi, cys[2], cns[2], cys[3], cns[3], alp.prior, bet.prior, type="R")
      if (OR.v2>gam2){
        index <- 1
        decision <- "escalation"
      }else{
        index <- 0
        decision <- "stay"
      }
    }
    else  if (is.na(cys[3]) | (cover.doses[3]==1)){
      gam1 <- optim.gamma.fn(cns[1], cns[2], phi, "L", alp.prior, bet.prior)$gamma 
      OR.v1 <- OR.values(phi, cys[1], cns[1], cys[2], cns[2], alp.prior, bet.prior, type="L")
      if (OR.v1>gam1){
        index <- -1
        decision <- "de-escalation"
      }else{
        index <- 0
        decision <- "stay"
      }
    }
    else  if (!(is.na(cys[1]) | is.na(cys[3]) | cover.doses[3]==1)){
      gam1 <- optim.gamma.fn(cns[1], cns[2], phi, "L", alp.prior, bet.prior)$gamma 
      gam2 <- optim.gamma.fn(cns[2], cns[3], phi, "R", alp.prior, bet.prior)$gamma 
      OR.v1 <- OR.values(phi, cys[1], cns[1], cys[2], cns[2], alp.prior, bet.prior, type="L")
      OR.v2 <- OR.values(phi, cys[2], cns[2], cys[3], cns[3], alp.prior, bet.prior, type="R")
      v1 <- OR.v1 > gam1
      v2 <- OR.v2 > gam2
      if (v1 & !v2){
        index <- -1
        decision <- "de-escalation"
      }else if (!v1 & v2){
        index <- 1
        decision <- "escalation"
      }else{
        index <- 0
        decision <- "stay"
      }
    }
  }
  curDose <- curDose+index
  out <- list(target=phi, cys=cys, cns=cns, decision=decision, curDose = curDose)
  class(out) <- "CFO"
  return(out)
}
