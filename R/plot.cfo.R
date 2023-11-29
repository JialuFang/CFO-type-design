#' Plot the simulation results for CFO-type and aCFO-type designs
#'
#' Plot the objects returned by other functions, including (1) dose allocation of a single trial;
#' (2) operating characteristics of multiple simulations, including selesction percentage and the
#' number of patients treated at each dose
#' 
#' @param x the object returned by other functions
#' @param ... ignored arguments
#' @param name the name of the object to be plotted.
#'             User doesn't need to input this parameter.
#'
#' @return \code{plot()} returns a figure or a series of figures depending on the object entered
#' 
#' @author Jialu Fang
#' 
#' @importFrom grDevices dev.flush dev.hold devAskNewPage
#' @importFrom graphics axis barplot mtext par plot
#' @export
#'
#' @examples
#' ############design without late-onset outcomes################
#' nsimu <- 100; phi <- 0.2; ncohort <- 12; cohortsize <- 3; init.level <- 1
#' p.true <-c(0.01, 0.05, 0.10, 0.14, 0.20, 0.26, 0.34)
#' add.args=list(alp.prior=phi, bet.prior=1-phi)
#' ## CFO design
#' CFOtrial <- aCFO.simu(phi, p.true, ncohort, init.level, cohortsize=3,add.args, accumulation = FALSE)
#' plot(CFOtrial)
#' CFOsimu <- CFO.oc (nsimu, design='CFO', phi, p.true, ncohort, init.level, cohortsize,
#'                     tau=NaN, accrual=NaN, tite.dist=NaN, accrual.dist=NaN, add.args)
#' plot(CFOsimu)
#' ## aCFO design
#' aCFOtrial <- aCFO.simu(phi, p.true, ncohort, init.level, cohortsize=3,add.args, 
#'                     accumulation = TRUE)
#' plot(aCFOtrial)
#' aCFOsimu <- CFO.oc (nsimu, design='aCFO', phi, p.true, ncohort, init.level, cohortsize,
#'                     tau=NaN, accrual=NaN, tite.dist=NaN, accrual.dist=NaN, add.args)
#' plot(aCFOsimu)
#' 
#' 
#' ##############design with late-onset outcomes################
#' tau <- 3; accrual <- 6; tite.dist <- 2; accrual.dist <- 1
#' ## TITE-CFO design
#' TITECFOtrial <- lateonset.simu (phi, p.true, tau, cohortsize, ncohort, accrual, tite.dist, 
#'                 accrual.dist, design='TITE-CFO', init.level, add.args)
#' plot(TITECFOtrial)
#' TITECFOsimu <- CFO.oc (nsimu, design='TITE-CFO', phi, p.true, ncohort, init.level, cohortsize,
#'                       tau, accrual, tite.dist, accrual.dist, add.args)
#' plot(TITECFOsimu)
#' ## TITE-aCFO design
#' TITEaCFOtrial <- lateonset.simu (phi, p.true, tau, cohortsize, ncohort, accrual, tite.dist, 
#'                  accrual.dist,design='TITE-aCFO', init.level, add.args)
#' plot(TITEaCFOtrial)
#' TITEaCFOsimu <- CFO.oc (nsimu, design='TITE-aCFO', phi, p.true, ncohort, init.level, cohortsize,
#'                       tau, accrual, tite.dist, accrual.dist, add.args)
#' plot(TITEaCFOsimu)
#' ## fCFO design
#' fCFOtrial <- lateonset.simu (phi, p.true, tau, cohortsize, ncohort, accrual, tite.dist,  
#'                 accrual.dist, design='fCFO', init.level, add.args)
#' plot(fCFOtrial)
#' fCFOsimu <- CFO.oc (nsimu, design='fCFO', phi, p.true, ncohort, init.level, cohortsize,
#'                       tau, accrual, tite.dist, accrual.dist, add.args)
#' plot(fCFOsimu)                      
#' ## f-aCFO design
#' faCFOtrial <- lateonset.simu (phi, p.true, tau, cohortsize, ncohort, accrual, tite.dist, 
#'                 accrual.dist, design='f-aCFO', init.level, add.args)
#' plot(faCFOtrial)
#' fCFOsimu <- CFO.oc (nsimu, design='f-aCFO', phi, p.true, ncohort, init.level, cohortsize,
#'                       tau, accrual, tite.dist, accrual.dist, add.args)
#' plot(faCFOsimu)
plot.cfo<- function (x,..., name = deparse(substitute(x)))
{
  new.obj = unlist(strsplit(name, split = "\\$"))
  strpattern = "none"
  if (length(new.obj) >= 2) {
    strpattern = new.obj[2]
  }
  assign("objectPlot", get(new.obj[1]))
  
  if (!is.element(strpattern, c("none", names(objectPlot)))) {
    warning("Please double check and specify the variable to be plotted...\n")
  }
  else {
    if (!is.null(objectPlot$simu.oc)) { #oc for one-dim multiple simulations
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
      dev.flush()
      dev.hold()
      par(mar = c(5, 6, 4, 2))
      bplot = barplot(objectPlot$selPercent*100, ylab = "selection percentage (%)",
                      ylim = c(0, 100), cex.names = 1, xaxt = "n",
                      cex.lab = 1.3)
      axis(1, at = bplot, labels = seq(1, length(objectPlot$selPercent)))
      mtext("Selection percentage", 3, line = 0, cex = 1.3)
      mtext("Dose level", 1, line = 2, cex = 1)
      dev.flush()
      dev.hold()
      bplot = barplot(objectPlot$nPatients, ylab = "number of patients",
                      ylim = c(0, max(objectPlot$nPatients))*1.3, cex.names = 1,
                      beside = FALSE, xaxt = "n", cex.lab = 1.3)
      axis(1, at = bplot, labels = seq(1, length(objectPlot$nPatients)))
      mtext("Patient allocation", 3, line = 0, cex = 1.3)
      mtext("Dose level", 1, line = 2, cex = 1)
      dev.flush()
      dev.hold()
      bplot = barplot(objectPlot$nTox, ylab = "number of toxicities",
                      ylim = c(0, max(objectPlot$nTox)*1.3), cex.names = 1,
                      beside = FALSE, xaxt = "n", cex.lab = 1.3)
      axis(1, at = bplot, labels = seq(1, length(objectPlot$nTox)))
      mtext("Observed toxicity", 3, line = 0, cex = 1.3)
      mtext("Dose level", 1, line = 2, cex = 1)
    } else if (!is.null(objectPlot$dose.list)) {
      par(mar = c(5, 6, 4, 2))
      x <- seq(1, length(objectPlot$dose.list))
      stepplot = plot(x, objectPlot$dose.list, type = "s", lwd = 2, col = "black", xaxt = "n",
                  main = "Step line for dose allocation", xlab = "Cohort index", ylab = "Dose level")
      axis(1, at = x, labels = x)
    }
  }
}
