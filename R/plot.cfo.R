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
#' @details \code{plot()} returns a figure or a series of figures depending on the object entered.
#'           Additionally, in the example, we set \code{nsimu=100} for testing time considerations. 
#'           In reality, \code{nsimu} is typically set to 5000 to ensure the accuracy of the results.
#'
#' @return \code{plot()} returns a figure or a series of figures depending on the object entered
#' 
#' @author Jialu Fang
#' 
#' @importFrom grDevices dev.flush dev.hold devAskNewPage
#' @importFrom graphics axis barplot mtext par plot
#' @import ggplot2
#' @export
#'
#' @examples
#' ############design without late-onset outcomes################
#' nsimu <- 100; phi <- 0.2; ncohort <- 12; cohortsize <- 3; init.level <- 1
#' p.true <-c(0.01, 0.05, 0.10, 0.14, 0.20, 0.26, 0.34)
#' add.args=list(alp.prior=phi, bet.prior=1-phi)
#' ## CFO design
#' CFOtrial <- CFO.simu(phi, p.true, ncohort, init.level, cohortsize=3,add.args, accumulation = FALSE)
#' plot(CFOtrial)
#' CFOsimu <- CFO.oc (nsimu, design='CFO', phi, p.true, ncohort, init.level, cohortsize,
#'                     tau=NaN, accrual=NaN, tite.dist=NaN, accrual.dist=NaN, add.args)
#' plot(CFOsimu)
#' ## aCFO design
#' aCFOtrial <- CFO.simu(phi, p.true, ncohort, init.level, cohortsize=3,add.args,
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
#' faCFOsimu <- CFO.oc (nsimu, design='f-aCFO', phi, p.true, ncohort, init.level, cohortsize,
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
      if (!is.null(objectPlot$total.time)){
        dose <- objectPlot$dose.list
        DLT <- objectPlot$DLT.list
        ncohort <- length(objectPlot$dose.list)
        cohortsize <- sum(objectPlot$dose.ns)/ncohort
        
        # Generate y_labels
        y_labels <- seq(1, max(dose))
        
        # Generate sequences for each patient
        sequences <- objectPlot$enter.times
        
        # Generate dose_levels for each patient
        dose_levels <- rep(dose, each = cohortsize)
        
        # Generate DLT_observed for each patient
        DLT_observed <- matrix(DLT, nrow = cohortsize, ncol = ncohort)
        
        new_seq <- ifelse(objectPlot$DLT.times!=0, sequences+objectPlot$DLT.times, NaN)
        new_y <- ifelse(objectPlot$DLT.times!=0, dose_levels+0.3, NaN)
        
        df <- data.frame(sequence = sequences, dose_levels = dose_levels, DLT_observed = DLT_observed)
        dfnew <- data.frame(sequence = sequences, dose_levels = dose_levels, new_seq = new_seq, new_y = new_y)
        dfnew <- na.omit(dfnew)
        
        # Create the plot
        p <- ggplot(df, aes(x = sequence, y = dose_levels)) +
          geom_point(aes(shape = factor(DLT_observed,levels=c(0,1,2))), color = 'black', size = 2.5) +
          geom_step(direction = 'hv', color = 'black') +
          scale_y_continuous(breaks = 1:length(y_labels), labels = y_labels) +
          labs(x = "Sequence of patients treated", 
               y = "Combined dose level",
               fill = 'DLT observed') +
          theme_minimal() +
          theme(text = element_text(size = 12), legend.title=element_blank(), legend.position = c(1, 0), legend.justification = c(1, 0)) +
          scale_shape_manual(values = c(1, 16, 4), labels = c('DLT not observed', 'DLT observed',"DLT time"), drop = FALSE)
        
        for (row in 1:(nrow(dfnew))){
          xuse=c(dfnew[row,"sequence"],dfnew[row,"new_seq"])
          yuse=c(dfnew[row,"dose_levels"],dfnew[row,"new_y"])
          dfuse <-data.frame(xuse=xuse, yuse=yuse)
          p <- p + geom_point(aes(x = xuse[2], y = yuse[2]), shape = 4,size = 2.5, data = dfuse)+
            geom_step(aes(x = xuse, y = yuse), data = dfuse,direction = 'vh',
                      linetype = 2)
        }
        print(p)
      }
      else{
        dose <- objectPlot$dose.list
        DLT <- objectPlot$DLT.list
        ncohort <- length(objectPlot$dose.list)
        cohortsize <- sum(objectPlot$dose.ns)/ncohort
        
        # Generate y_labels
        y_labels <- seq(1, max(dose))
        
        # Generate sequences for each patient
        sequences <- 1:(ncohort * cohortsize)
        
        # Generate dose_levels for each patient
        dose_levels <- rep(dose, each = cohortsize)
        
        # Generate DLT_observed for each patient
        DLT_observed <- matrix(DLT, nrow = cohortsize, ncol = ncohort)
        
        df <- data.frame(sequence = sequences, dose_levels = dose_levels, DLT_observed = DLT_observed)
        
        # Create the plot
        p <- ggplot(df, aes(x = sequence, y = dose_levels)) +
          geom_point(aes(fill = as.factor(DLT_observed)), color = 'black', shape = 21, size = 2.5) +
          geom_step(direction = 'hv', color = 'black') +
          scale_y_continuous(breaks = 1:length(y_labels), labels = y_labels) +
          labs(x = "Sequence of patients treated", 
               y = "Combined dose level",
               fill = 'DLT observed') +
          theme_minimal() +
          theme(text = element_text(size = 12), legend.title=element_blank(), legend.position = c(1, 0), legend.justification = c(1, 0)) +
          scale_fill_manual(values = c('white', 'black'), labels = c('DLT not observed', 'DLT observed'))
        
        # Display the plot
        print(p)
      }
    }
  }
}
