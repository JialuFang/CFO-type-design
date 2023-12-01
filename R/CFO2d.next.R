
#' CFO2d.next function
#'
#' This function is used to determine the next dose level in a two-dimensional 
#' dose-finding clinical trial setting. It is based on the Continual Feedback On (CFO) method.
#'
#' @param phi The efficacy threshold.
#' @param tys A matrix of the total number of successes for each dose combination.
#' @param tns A matrix of the total number of patients for each dose combination.
#' @param cur A vector of the current dose indices in the horizontal and vertical direction.
#' @param alp.prior Parameter alpha for the prior beta distribution.
#' @param bet.prior Parameter beta for the prior beta distribution.
#' @param seed Optional; an integer to set as the seed of the random number generator for reproducible results.
#' 
#' @return A vector with the next dose indices in the horizontal and vertical directions.
#' @export
#' @examples
#' tns <- matrix(c(3, 3, 3, 0,
#'                 0, 0, 6, 0,
#'                 0, 0, 0, 0), 
#'               nrow = 3, ncol = 4, byrow = TRUE)
#' tys <- matrix(c(0, 0, 1, 0,
#'                 0, 0, 2, 0,
#'                 0, 0, 0, 0), 
#'               nrow = 3, ncol = 4, byrow = TRUE)
#' cur <- c(2,3)
#' CFO2d.next(0.3, tys, tns, cur, seed = 1)
#' 

CFO2d.next <- function(phi, tys, tns, cur, alp.prior=0.1, bet.prior=0.9, seed=NULL){
  
  # posterior probability of pj >= phi given data
  post.prob.fn <- function(phi, y, n, alp.prior, bet.prior){
    alp <- alp.prior + y 
    bet <- bet.prior + n - y
    1 - pbeta(phi, alp, bet)
  }
  
  
  prob.int <- function(phi, y1, n1, y2, n2, alp.prior, bet.prior){
    alp1 <- alp.prior + y1
    alp2 <- alp.prior + y2
    bet1 <- alp.prior + n1 - y1
    bet2 <- alp.prior + n2 - y2
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
    }else if (type=="D"){
      pC <- 1 - ps$p2
      pD <- 1 - ps$p1
      oddsC <- pC/(1-pC)
      oddsD <- pD/(1-pD)
      OR <- oddsC*oddsD
    }else if (type=="U"){
      pC <- 1 - ps$p1
      pU <- 1 - ps$p2
      oddsC <- pC/(1-pC)
      oddsU <- pU/(1-pU)
      OR <- (1/oddsC)/oddsU
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
  
  make.decision.1dCFO.fn <- function(phi, cys, cns, alp.prior, bet.prior, cover.doses, diag=FALSE){
    if (cover.doses[2] == 1){
      return(1)
    }else{
      if (is.na(cys[1]) & (cover.doses[3]==1)){
        return(2)
      }else  if (is.na(cys[1]) & (!(cover.doses[3]==1))){
        gam2 <- optim.gamma.fn(cns[2], cns[3], phi, "R", alp.prior, bet.prior)$gamma
        OR.v2 <- OR.values(phi, cys[2], cns[2], cys[3], cns[3], alp.prior, bet.prior, type="R")
        if (OR.v2>gam2){
          return(3)
        }else{
          return(2)
        }
      }else  if (is.na(cys[3]) | (cover.doses[3]==1)){
        gam1 <- optim.gamma.fn(cns[1], cns[2], phi, "L", alp.prior, bet.prior)$gamma
        OR.v1 <- OR.values(phi, cys[1], cns[1], cys[2], cns[2], alp.prior, bet.prior, type="L")
        if (OR.v1>gam1){
          return(1)
        }else{
          return(2)
        }
        
      }else  if (!(is.na(cys[1]) | is.na(cys[3]) | cover.doses[3]==1)){
        gam1 <- optim.gamma.fn(cns[1], cns[2], phi, "L", alp.prior, bet.prior)$gamma
        gam2 <- optim.gamma.fn(cns[2], cns[3], phi, "R", alp.prior, bet.prior)$gamma
        OR.v1 <- OR.values(phi, cys[1], cns[1], cys[2], cns[2], alp.prior, bet.prior, type="L")
        OR.v2 <- OR.values(phi, cys[2], cns[2], cys[3], cns[3], alp.prior, bet.prior, type="R")
        v1 <- OR.v1 > gam1
        v2 <- OR.v2 > gam2
        if (v1 & !v2){
          return(1)
        }else if (!v1 & v2){
          return(3)
        }else{
          return(2)
        }
      }
    }
  }

  make.decision.2dCFO.fn <- function(phi, cys, cns, alp.prior, bet.prior, cover.doses){
    cidx.A <- 0
    cidx.B <- 0
    # horizontal direction
    idx.chg.A <- make.decision.1dCFO.fn(phi, cys[2,], cns[2,], alp.prior, bet.prior, cover.doses[2,]) - 2
    # vertical direction
    idx.chg.B <- make.decision.1dCFO.fn(phi, cys[,2], cns[,2], alp.prior, bet.prior, cover.doses[,2]) - 2
    
    if (idx.chg.A == 1 & idx.chg.B == 1){
      ### horizontal and vertical only
      
      OR.R <- OR.values(phi, cys[2,2], cns[2,2], cys[2,3], cns[2,3], alp.prior, bet.prior, type="R")
      OR.U <- OR.values(phi, cys[2,2], cns[2,2], cys[3,2], cns[3,2], alp.prior, bet.prior, type="R")
      
      message(paste('OR.R: ',OR.R))
      message(paste('OR.U: ',OR.U))
      
      if (OR.R == OR.U){
        rand <- rbinom(1,1,0.5)
        if(rand == 0){
          cidx.A <- 1
        } else {
          cidx.B <- 1
        }
      } else if (OR.R > OR.U){
        cidx.B <- 1
      } else {
        cidx.A <- 1
      }
      
    } else if (idx.chg.A == -1 & idx.chg.B == -1){
      if (is.na(cys[2,1]) & is.na(cys[1,2])){
        cidx.A <- 0
        cidx.B <- 0
      } else if (is.na(cys[2,1])){
        cidx.A <- -1
      } else if (is.na(cys[1,2])){
        cidx.B <- -1
      } else {
        OR.L <- OR.values(phi, cys[2,2], cns[2,2], cys[2,1], cns[2,1], alp.prior, bet.prior, type="L")
        OR.D <- OR.values(phi, cys[2,2], cns[2,2], cys[1,2], cns[1,2], alp.prior, bet.prior, type="L")
        
        message(paste('OR.L: ',OR.L))
        message(paste('OR.D: ',OR.D))
        
        if (OR.L == OR.D){
          rand <- rbinom(1,1,0.5)
          if(rand == 0){
            cidx.A <- -1
          } else {
            cidx.B <- -1
          }
        } else if (OR.L > OR.D){
          cidx.B <- -1
        } else {
          cidx.A <- -1
        }
      }
    } else if (idx.chg.A == 1 & idx.chg.B == -1){
      DCR <- make.decision.1dCFO.fn(phi, c(cys[1,2],cys[2,2],cys[2,3]), c(cns[1,2],cns[2,2],cns[2,3]), alp.prior, 
                                    bet.prior, c(cover.doses[1,2],cover.doses[2,2],cover.doses[2,3])) - 2
      if (DCR == 1){
        cidx.B <- 1
      } else if (DCR == -1){
        cidx.A <- -1
      }
    } else if (idx.chg.A == -1 & idx.chg.B == 1){
      LCU <- make.decision.1dCFO.fn(phi, c(cys[2,1],cys[2,2],cys[3,2]), c(cns[2,1],cns[2,2],cns[3,2]), alp.prior, 
                                    bet.prior, c(cover.doses[2,1],cover.doses[2,2],cover.doses[3,2])) - 2
      if (LCU == 1){
        cidx.A <- 1
      } else if (LCU == -1){
        cidx.B <- -1
      }
    } else if (idx.chg.A == 1 & idx.chg.B == 0){
      cidx.B <- 1
    } else if (idx.chg.A == 0 & idx.chg.B == 1){
      cidx.A <- 1
    } else if (idx.chg.A == -1 & idx.chg.B == 0){
      cidx.B <- -1
    } else if (idx.chg.A == 0 & idx.chg.B == -1){
      cidx.A <- -1
    }
    
    return (c(cidx.A, cidx.B))
  }
  
  
  overdose.fn <- function(phi, add.args=list()){
    y <- add.args$y
    n <- add.args$n
    alp.prior <- add.args$alp.prior
    bet.prior <- add.args$bet.prior
    pp <- post.prob.fn(phi, y, n, alp.prior, bet.prior)
    # message(pp)
    if ((pp >= 0.95) & (add.args$n>=3)){
      return(TRUE)
    }else{
      return(FALSE)
    }
    
  }
  

  set.seed(seed)
  ndose.A <- length(tns[,1])
  ndose.B <- length(tns[1,])
  cidx.A <- cur[1]
  cidx.B <- cur[2]
  
  tover.doses <- matrix(0, ndose.A, ndose.B) # Whether each dose is overdosed or not, 1 yes
  
  cy <- tys[cidx.A, cidx.B]
  cn <- tns[cidx.A, cidx.B]
  
  add.args <- list(y=cy, n=cn, tys=tys, tns=tns, cidx.A=cidx.A, cidx.B=cidx.B)
  
  if (tover.doses[1,1] == 1){
    message("The lowest dose is overly toxic.")
    return(c(99, 99))
  }
  
  if (cidx.A!=1 & cidx.B!=1 & cidx.A!=ndose.A & cidx.B!=ndose.B){
    # no boundary
    cys <- tys[(cidx.A-1):(cidx.A+1), (cidx.B-1):(cidx.B+1)]
    cns <- tns[(cidx.A-1):(cidx.A+1), (cidx.B-1):(cidx.B+1)]
    cover.doses <- tover.doses[(cidx.A-1):(cidx.A+1), (cidx.B-1):(cidx.B+1)]
  } else if (cidx.A==1 & cidx.B==1){
    # (1, 1)
    cys <- rbind(c(NA,NA,NA),cbind(c(NA,NA),tys[1:2,1:2]))
    cns <- rbind(c(NA,NA,NA),cbind(c(NA,NA),tns[1:2,1:2]))
    cover.doses <- rbind(c(NA,NA,NA),cbind(c(NA,NA),tover.doses[1:2,1:2]))
  } else if (cidx.A==ndose.A & cidx.B==ndose.B){
    # (nA, nB)
    cys <- rbind(cbind(tys[(cidx.A-1):cidx.A,(cidx.B-1):cidx.B],c(NA,NA)), c(NA,NA,NA))
    cns <- rbind(cbind(tns[(cidx.A-1):cidx.A,(cidx.B-1):cidx.B],c(NA,NA)), c(NA,NA,NA))
    cover.doses <- rbind(cbind(tover.doses[(cidx.A-1):cidx.A,(cidx.B-1):cidx.B],c(NA,NA)), c(NA,NA,NA))
  } else if (cidx.A==1 & cidx.B==ndose.B){
    # (1, nB) 
    cys <- rbind(c(NA,NA,NA),cbind(tys[1:2,(cidx.B-1):cidx.B],c(NA,NA)))
    cns <- rbind(c(NA,NA,NA),cbind(tns[1:2,(cidx.B-1):cidx.B],c(NA,NA)))
    cover.doses <- rbind(c(NA,NA,NA),cbind(tover.doses[1:2,(cidx.B-1):cidx.B],c(NA,NA)))
  } else if (cidx.A==ndose.A & cidx.B==1){
    # (nA, 1) 
    cys <- rbind(cbind(c(NA,NA), tys[(cidx.A-1):cidx.A,1:2]),c(NA,NA,NA))
    cns <- rbind(cbind(c(NA,NA), tns[(cidx.A-1):cidx.A,1:2]),c(NA,NA,NA))
    cover.doses <- rbind(cbind(c(NA,NA), tover.doses[(cidx.A-1):cidx.A,1:2]),c(NA,NA,NA))
  } else if (cidx.A==1 & cidx.B!=1){
    # (1, 2:(nB-1))
    cys <- rbind(c(NA,NA,NA), tys[1:2, (cidx.B-1):(cidx.B+1)])
    cns <- rbind(c(NA,NA,NA), tns[1:2, (cidx.B-1):(cidx.B+1)])
    cover.doses <- rbind(c(NA,NA,NA), tover.doses[1:2, (cidx.B-1):(cidx.B+1)])
  } else if (cidx.A!=1 & cidx.B==1){
    # (2:(nA-1), 1)
    cys <- cbind(c(NA,NA,NA), tys[(cidx.A-1):(cidx.A+1), 1:2])
    cns <- cbind(c(NA,NA,NA), tns[(cidx.A-1):(cidx.A+1), 1:2])
    cover.doses <- cbind(c(NA,NA,NA), tover.doses[(cidx.A-1):(cidx.A+1), 1:2])
  } else if (cidx.A==ndose.A & cidx.B!=ndose.B){
    # (nA, 2:(nB-1))
    cys <- rbind(tys[(ndose.A-1):ndose.A, (cidx.B-1):(cidx.B+1)], c(NA,NA,NA))
    cns <- rbind(tns[(ndose.A-1):ndose.A, (cidx.B-1):(cidx.B+1)], c(NA,NA,NA))
    cover.doses <- rbind(tover.doses[(ndose.A-1):ndose.A, (cidx.B-1):(cidx.B+1)], c(NA,NA,NA))
  } else if (cidx.A!=ndose.A & cidx.B==ndose.B){
    # (2:(nA-1), nB)
    cys <- cbind(tys[(cidx.A-1):(cidx.A+1), (cidx.B-1):cidx.B], c(NA,NA,NA))
    cns <- cbind(tns[(cidx.A-1):(cidx.A+1), (cidx.B-1):cidx.B], c(NA,NA,NA))
    cover.doses <- cbind(tover.doses[(cidx.A-1):(cidx.A+1), (cidx.B-1):cidx.B], c(NA,NA,NA))
  } else {
    message('no such case')
  }
  
  idx.chg <- make.decision.2dCFO.fn(phi, cys, cns, 0.1, 0.9, cover.doses)

  return(idx.chg+cur)
}



tns <- matrix(c(3, 3, 3, 0,
                0, 0, 6, 0,
                0, 0, 0, 0), 
              nrow = 3, ncol = 4, byrow = TRUE)

tys <- matrix(c(0, 0, 1, 0,
                0, 0, 2, 0,
                0, 0, 0, 0), 
              nrow = 3, ncol = 4, byrow = TRUE)
cur <- c(2,3)
CFO2d.next(0.3, tys, tns, cur)
