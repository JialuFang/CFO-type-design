source("utilities.R")
source("./setting/Sce_5dose_phi33.R")

MTDs <- c()
res.ls <- list()
ndose <- length(p.trues[[1]])
nsimu <- 5000

for (flag in 1:length(p.trues)){
  file.name <- paste0("./use_results/","MTDSimu_",ndose, "dose_phi_", target*100, "_fix_",  flag, ".Rdata")
  load(file.name)
  #res.ls[[flag]] <-  post.process.random(results)[c("fcfo", "facfo", "titecfo", "titeacfo","benchacfo"), ]
  res.ls[[flag]] <-  post.process.random(results)[c("cfo","acfo","crm","boin"), ]
  print(res.ls[[flag]])
  flag <- flag + 1
}
Ress <- list()
for (nam in names(res.ls[[1]])){
  cRes <- do.call(rbind, lapply(res.ls, function(x){x[[nam]]}))
  colnames(cRes) <- row.names(res.ls[[1]])
  Ress[[nam]] <- cRes
}

Ress$MTD.Allo[6,] <- c(0,0,0,0)
Ress$Over.Allo[6,] <- c(0,0,0,0)


singlePlot <- function(cRes, main="", ylab="Percentage (%)"){
  par(cex.axis = 1.4)
  plot(cRes[, 1]*100, type = "b", col=1, lwd=2, lty=1, pch=1, ylim=c(min(cRes)-0.02, max(cRes)+0.02)*100, 
       xlab="Senarios", main=main, cex.lab = 1.4, cex.main = 1.5,
       xaxt="n", ylab=ylab)
  lines(cRes[, 3]*100, type = "b", col=3, lwd=2, lty=3, pch=3)
  lines(cRes[, 4]*100, type = "b", col=4, lwd=2, lty=4, pch=4)
  lines(cRes[, 2]*100, type = "b", col=2, lwd=2, lty=2, pch=2)
  #lines(cRes[, 5]*100, type = "b", col=alpha(5, 0.8), lwd=1, lty=5, pch=5)
  axis(1, at=1:length(p.trues), labels=1:length(p.trues))
  #legend("topleft", toupper(colnames(cRes)), col=1:3, lty=1:3, pch=1:3, lwd=2)
}


fig.name <- paste0("./use/", "SimuLate_", ndose, "dose_phi_", 100*target,  "_fix",  ".jpeg")
png(filename=fig.name, unit="in", height=6, width=10, res=300)
par(mfrow=c(2, 3), oma = c(3,1,1,1))
singlePlot(Ress[[1]], main="(A) MTD Selection")
singlePlot(Ress[[2]], main="(B) MTD Allocation")
singlePlot(Ress[[3]], main="(C) Overdose Selection")
singlePlot(Ress[[4]], main="(D) Overdose Allocation")
singlePlot(Ress[[6]], main="(E) Average DLT Rate")
#singlePlot(Ress[[7]]/100, main="(F) Average Trial Duration", ylab="Time (in months)")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')

legend("bottom", c("CFO","aCFO","CRM","BOIN"), col=c(1:4), lty=1:4, pch=1:4, lwd=2, 
      xpd=TRUE, horiz = TRUE, cex = 1.5, seg.len=1)

#legend("bottom", c("fCFO","f-aCFO","TITE-CFO", "TITE-aCFO","Benchmark aCFO"), col=c(1:5), lty=1:5, pch=1:5, lwd=2, 
#      xpd=TRUE, horiz = TRUE, cex = 1.5, seg.len=1)
par(mfrow=c(1, 1))
dev.off()

