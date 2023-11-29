#' 
#' Determination of the dose level for next cohort in aCFO design
#' 
#' In aCFO design, use the function to determine the dose movement based on the toxicity outcomes of the enrolled cohorts.
#'
#' @usage aCFO.next(phi, tys, tns, tover.doses=c(), curDose, 
#'        add.args=list(alp.prior=phi, bet.prior=1-phi))
#'
#' @param phi the target DLT rate.
#' @param tys the current number of DLTs observed in patients for all dose levels.
#' @param tns the current number of patients for all dose levels.
#' @param tover.doses whether the dose level (from the first to last dose level) is over-toxic or not. 
#'                    The value is set as 1 if the dose level is overly toxicity; otherwise, it is set to 0.
#'                    It should be predetermined. If not predetermined, it is set to NaN and is assigned values of 0 
#'                    for each dose level.
#' @param curDose the current dose level.
#' @param add.args additional parameters, usually set as list(alp.prior=phi, bet.prior=1-phi) by default. \code{alp.prior} 
#'                 and \code{bet.prior} represent the parameters of the prior distribution for the true DLT rate at 
#'                 any dose level. This prior distribution is specified as Beta( \code{alpha.prior}, \code{beta.prior}).
#'
#' @details The aCFO design design incorporate the dose information of all positions (from the lowest to the 
#'          highest dose levels) into the trial decision-making. Prior to assigning dose levels for new patient 
#'          cohorts, aCFO compares the evidence from the current dose level with all doses to its left and right. 
#'          This design is rooted in the odds ratio, specifically \eqn{O_C / \overline{O}_{J}} and 
#'          \eqn{\overline{O}_C / O_R} from the CFO design. By aggregating odds ratios from the left and right sides, 
#'          it forms two collective statistics: \eqn{ {\rm OR}_L =\sum_{i=1}^{J} O_C/ \overline{O}_{L_i} } for dose 
#'          de-escalation (movement to the left) and \eqn{ {\rm OR}_R = \sum_{i=1}^{H} \overline{O}_C / O_{R_i} } 
#'          for dose escalation (movement to the right), where J and H represent the counts of doses on the left and 
#'          right sides of the current dose, respectively. For the new statistic \eqn{ {\rm OR}_L } and \eqn{ {\rm OR}_R }, 
#'          their corresponding thresholds are derived by summing up its individual thresholds \eqn{\gamma_{L_i}} and 
#'          \eqn{\gamma_{R_i}}, i.e., \eqn{\sum_{i=1}^{J}\gamma_{L_i}} and \eqn{\sum_{i=1}^{H}\gamma_{R_i}}. \cr
#'          Besides，The aCFO design retains the same early stopping criteria as the CFO design. Overall，while preserving 
#'          the nature of the CFO design (model-free and calibration-free), the aCFO designs enhance the efﬁciency by 
#'          incorporating more dose information. 
#'          
#' @note    The \code{tover.doses()} should be predetermined. It can be calculated based on the number of patients enrolled 
#'          ($tns) and the number of DLTs observed ($tys), following the overdose control rule. If not predetermined, it is 
#'          set to NaN and is assigned default values of 0 for each dose level.
#' 
#' @return The \code{aCFO.next()} function returns a list object comprising the following elements: the target DLT 
#'         rate ($target), the current number of DLTs and patients for the left, current, and right dose levels ($cys and $cns), 
#'         the decision in the aCFO design of whether to move to the left or right dose level for the next cohort ($decision), 
#'         and the current dose after movement ($curDose).
#'         
#' @author Jialu Fang 
#' 
#' @references Jin, H., & Yin, G. (2022). CFO: Calibration-free odds design for phase I/II clinical trials. 
#'             \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066.
#'
#' @examples
#' ## determine the dose level for the next cohort of new patients
#' tys <- c(0,0,1,0,0,0,0); tns <- c(3,3,6,0,0,0,0)
#' aCFO.next(phi=0.2, tys=tys, tns=tns, tover.doses=c(0,0,0,0,0,0,0), curDose=3,
#'          add.args=list(alp.prior=0.2, bet.prior=0.8))
#' 
#' @import stats
#' @export
aCFO.next <- function(phi, tys, tns, tover.doses=c(), curDose, add.args=list(alp.prior=phi, bet.prior=1-phi)){
  
###############################################################################
###############define the functions used for main function#####################
###############################################################################
  
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
  
  OR.union.values <- function(phi, ctns, ctys, alp.prior, bet.prior, type){
    ndose <- length(ctys)
    if (type=="L"){
      OR.list <- rep(0, ndose-1)
      for (i in 1:(ndose-1)){
        OR.list[i] <- OR.values(phi, ctys[i], ctns[i], ctys[ndose], ctns[ndose], alp.prior, bet.prior, type)
      }
    }else if (type=="R"){
      OR.list <- rep(0, ndose-1)
      for (i in 2:ndose){
        OR.list[i-1] <- OR.values(phi, ctys[1], ctns[1], ctys[i], ctns[i], alp.prior, bet.prior, type)
      }
    }
    return(sum(OR.list))
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
  
  optim.gamma.union.fn <- function(ctns, phi, type, alp.prior, bet.prior){
    ndose <- length(ctns)
    if (type == "L"){
      gamma.list <- rep(0, ndose-1)
      for (i in 1:(ndose-1)){
        gamma.list[i] <- optim.gamma.fn(ctns[i], ctns[ndose], phi, type, alp.prior, bet.prior)$gamma
      }
    }else if (type == "R"){
      gamma.list <- rep(0, ndose-1)
      for (i in 2:ndose){
        gamma.list[i-1] <- optim.gamma.fn(ctns[1], ctns[i], phi, type, alp.prior, bet.prior)$gamma
      }
    }
    return(sum(gamma.list))
  }
  
  ###############################################################################
  ############################MAIN DUNCTION######################################
  ###############################################################################
  ndose <- length(tys)
  if (is.null(tover.doses)){
    tover.doses <- rep(0,ndose)
    warning("tover.doses is set to NaN and is assigned default values of ", 
            paste(tover.doses, collapse = ","))
  }
  
  if (is.null(add.args$alp.prior)){
    add.args <- c(add.args, list(alp.prior=phi, bet.prior=1-phi))
  }
  alp.prior <- add.args$alp.prior
  bet.prior <- add.args$bet.prior
  
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
      OR.v2 <- OR.union.values(phi, tns[curDose:ndose], tys[curDose:ndose], alp.prior, bet.prior, type="R")
      gam2 <- optim.gamma.union.fn(tns[curDose:ndose], phi, "R", alp.prior, bet.prior)
      if (OR.v2>gam2){
        index <- 1
        decision <- "escalation"
      }else{
        index <- 0
        decision <- "stay"
      }
    }
    else  if (is.na(cys[3]) | (cover.doses[3]==1)){
      gam1 <- optim.gamma.union.fn(tns[1:curDose], phi, "L", alp.prior, bet.prior)
      OR.v1 <- OR.union.values(phi, tns[1:curDose], tys[1:curDose], alp.prior, bet.prior, type="L")
      if (OR.v1>gam1){
        index <- -1
        decision <- "de-escalation"
      }else{
        index <- 0
        decision <- "stay"
      }
    }
    else  if (!(is.na(cys[1]) | is.na(cys[3]) | cover.doses[3]==1)){
      gam1 <- optim.gamma.union.fn(tns[1:curDose], phi, "L", alp.prior, bet.prior)
      gam2 <- optim.gamma.union.fn(tns[curDose:ndose], phi, "R", alp.prior, bet.prior)
      OR.v1 <- OR.union.values(phi, tns[1:curDose], tys[1:curDose], alp.prior, bet.prior, type="L")
      OR.v2 <- OR.union.values(phi, tns[curDose:ndose], tys[curDose:ndose], alp.prior, bet.prior, type="R")
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
  out <- list(target=phi, tys=tys, tns=tns, decision=decision, curDose = curDose)
  return(out)
}

