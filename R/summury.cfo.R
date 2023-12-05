#' Generate descriptive summary for objects returned by other functions
#'
#' Generate descriptive summary for objects returned by other functions.
#'
#' @param object the object returned by other functions.
#' @param ... ignored arguments
#'
#'
#' @details \code{summary()} prints the objects returned by other functions. Additionally, in the example, 
#'          we set \code{nsimu=100} for testing time considerations. In reality, \code{nsimu} is typically 
#'          set to 5000 to ensure the accuracy of the results.
#'
#' @return \code{summary()} prints the objects returned by other functions.
#'
#' @author Jialu Fang
#'
#' @examples
#' ## setting
#' nsimu <- 100; phi <- 0.2; ncohort <- 12; cohortsize <- 3; init.level <- 1
#' p.true <-c(0.01, 0.05, 0.10, 0.14, 0.20, 0.26, 0.34)
#' add.args=list(alp.prior=phi, bet.prior=1-phi)
#' tau <- 3; accrual <- 6; tite.dist <- 2; accrual.dist <- 1
#' 
#' ## summarize the object returned by CFO.next()
#' decision <- CFO.next(phi=0.2, cys=c(0,1,0), cns=c(3,6,0), cover.doses=c(0,0,0), 
#'                      curDose=3, add.args)
#' summary(decision)
#' 
#' ## summarize the object returned by CFO.simu()
#' aCFOtrial <- CFO.simu(phi, p.true, ncohort, init.level, cohortsize=3,add.args, 
#'                     accumulation = TRUE)
#' summary(aCFOtrial)
#' 
#' ## summarize the object returned by lateonset.simu()
#' faCFOtrial <- lateonset.simu (phi, p.true, tau, cohortsize, ncohort, accrual, tite.dist, 
#'                 accrual.dist, design='f-aCFO', init.level, add.args)
#' summary(faCFOtrial)
#' 
#' ## summarize the object returned by CFO.oc()
#' faCFOsimu <- CFO.oc (nsimu, design='f-aCFO', phi, p.true, ncohort, init.level, cohortsize,
#'                       tau, accrual, tite.dist, accrual.dist, add.args)
#' summary(faCFOsimu)
#' 
summary.cfo<- function (object, ...)
{
  if (!is.null(object$simu.oc)) {
    if (object$errStop == 0){
      cat("No instance of early stopping was observed in",
          object$simu.setup$nsimu, "simulations. \n")
    }else{
      cat("In", object$simu.setup$nsimu, "simulations, early stopping occurred",
          object$errStop, "times. \n")
      cat("Among simulations where early stopping did not occur:")
    }
    
    cat("selection percentage at each dose level:\n")
    cat(formatC(object$selPercent, digits = 3, format = "f"),
        sep = "  ", "\n")
    cat("number of patients treated at each dose level:\n")
    cat(formatC(object$nPatients, format = "d"),
        sep = "  ", "\n")
    cat("number of toxicity observed at each dose level:\n")
    cat(formatC(object$nTox,  format = "d"),
        sep = "  ", "\n")
    cat("percentage of correct selection of the MTD:", 
        formatC(object$simu.oc$MTDSel, digits = 3, format = "f"), "\n")
    cat("percentage of patients allocated to the MTD:", 
        formatC(object$simu.oc$MTDAllo, digits = 3, format = "f"), "\n")
    cat("percentage of selecting a dose above the MTD:",
        formatC(object$simu.oc$overSel, digits = 3, format = "f")," \n")
    cat("percentage of allocating patients at dose levels above the MTD:",
        formatC(object$simu.oc$overAllo, digits = 3, format = "f")," \n")
    cat("percentage of the patients suffering DLT:",
        formatC(object$simu.oc$averDLT, digits = 3, format = "f")," \n")
    
    if (!is.null(object$simu.oc$averDur)){
      cat("Average trial duration if trials with late-onset toxicities:",
          formatC(object$simu.oc$averDur, digits = 3, format = "f"),"% \n")
    }
  }
  
  if(!is.null(object$MTD)){
    if (length(object$MTD) == 1) {
      if (object$MTD == 99) {
        cat("All tested doses are overly toxic. No MTD should be selected! \n\n")
      }
      else {
        cat("The selected MTD is dose level ", object$MTD, "\n\n")
      }
      cat("For",length(object$dose.list),"cohorts, the dose level assigned to each cohort is: \n")
      cat(formatC(object$dose.list, format = "d"), sep = "  ", "\n")
      cat("number of toxicity observed at each dose level:\n")
      cat(formatC(object$DLT.ns, format = "d"), sep = "  ", "\n")
      cat("number of patients treated at each dose level:\n")
      cat(formatC(object$dose.ns, format = "d"), sep = "  ", "\n")
      if (!is.null(object$total.time)){
        cat("the duration of the trial in months:",
            formatC(object$total.time, digits = 3, format = "f"),"% \n")
      }
    }
  }
  
  if(!is.null(object$decision)){
    cat("The decision regarding the direction of movement is", object$decision, "\n")
    cat("The next cohort will be assigned to dose level", object$curDose, "\n")
  }
}

