citation <- "Mignan, A., J. Woessner (2012). Estimating the magnitude of completeness for earthquake catalogs, Community Online Resource for Statistical Seismicity Analysis, doi:10.5078/corssa-00180805. Available at http://www.corssa.org."
print(citation)

data <- read.csv('Rcatalogo.csv')
data <- data.frame(data)
bin <- 0.1
cat <- data
mag <- cat$mag

library(nlstools)
library(dplyr)

fmd <- function(mag,mbin){
  mi <- seq(min(round(mag/mbin)*mbin), max(round(mag/mbin)*mbin), mbin)
  nbm <- length(mi)
  cumnbmag <- numeric(nbm)
  nbmag <- numeric(nbm)
  for(i in 1:nbm) cumnbmag[i] <- length(which(mag > mi[i]-mbin/2))
  cumnbmagtmp <- c(cumnbmag,0)
  nbmag <- abs(diff(cumnbmagtmp))
  res <- list(m=mi, cum=cumnbmag, noncum=nbmag)
  return(res)
}

maxc <- function(mag,mbin){
  FMD <- fmd(mag,mbin)
  Mc <- FMD$m[which(FMD$noncum == max(FMD$noncum))[1]]
  return(list(Mc=Mc))
}

maxc.var <- maxc(mag,bin)

plot(FMD$m, FMD$cum, 
     log = "y", 
     xlab = "Magnitude", 
     ylab = "Log (N / Ntot)", 
     col = "red", pch = 1)

points(FMD$m, FMD$noncum, 
       pch = 2, col = "blue")

abline(v = maxc.var$Mc, lwd = 1)

legend("topright", c("Cumulative", "Incremental", "Mc: 1.8 (MAXC)"), 
       cex = 0.8, pch = c(1, 2, NA), 
       col = c("red", "blue", "black"),
       lwd = c(NA, NA, 1))

bootstrap_estimation <- function(data, bin, nbsample) {
  set.seed(123)
  
  Mc_bootstrap <- numeric(nbsample)
  
  for (i in 1:nbsample) {
    resampled_data <- sample(data, replace = TRUE)
    maxc_result <- maxc(resampled_data, bin)
    Mc_bootstrap[i] <- maxc_result$Mc
    
  }
  
  result_table <- data.frame(
    Parameter = c('Mc'),
    Mean = c(mean(Mc_bootstrap)),
    Dev.std = c(sd(Mc_bootstrap))
  )
  
  print(result_table)
  
  return(
    list(
      Mc_bootstrap = Mc_bootstrap
    )
  )
}


nbsample <- 500
result <- bootstrap_estimation(mag, bin, nbsample)

