# quick script to compute percentiles and plot distributions of quality scores and depths

setwd("TestDPdistr/")

for (fin in c("scaffold1.qc", ""scaffold2.qc", ""scaffold3.qc", ""scaffold4.qc", ""scaffold5.qc") {
  
  #fin <- "scaffold1.qc"
  cat("", file=paste(fin,".info",sep="",collapse=""))

  namesF <- read.table("../All_bams_names.txt")
  nF <- namesF$V1

  pdf(paste(fin,".pdf",sep="",collapse=""))
  par(mfrow=c(3,1))

## sample depth
  deps <- t(read.table(paste(fin,".depthSample", sep="",collapse=""),head=F, stringsAsFactors=F))
## per sample
  max_x <- 2 * nrow(deps)/ncol(deps)
  for (i in 1:ncol(deps)) {
    barplot(height=deps[1:max_x,i], names.arg=seq(1,max_x)-1, xlab="Sample Depth", ylab="Counts", main=nF[i])
  }

invisible(dev.off())

M <- c()
for (i in 1:ncol(deps)) {
  n <- c()
  for (j in 1:length(deps[,i])) {
    n <- c(n, j*deps[,i][j])
  }
  M <- c(M, sum(n)/sum(deps[,i]))
}

MeanCov <- cbind(M, namesF)[,1:2]

#MeanCov$MaxDP <- MeanCov$M*3
MeanCov$MaxDP <- ifelse(MeanCov$M >20.0, MeanCov$M*2, MeanCov$M*3)
MeanCov$MinDP <- MeanCov$M/2
MeanCov$MinDP[MeanCov$MinDP < 4] = 4
MeanCov$MinDP  <- round(MeanCov$MinDP)
MeanCov$MaxDP  <- round(MeanCov$MaxDP)
print(MeanCov)
}

##MeanCov has been used to create MAxDpPerSample.txt file manually
