citation <- "Mignan, A., J. Woessner (2012). Estimating the magnitude of completeness for earthquake catalogs, Community Online Resource for Statistical Seismicity Analysis, doi:10.5078/corssa-00180805. Available at http://www.corssa.org."
print(citation)

bin <- 0.1

data_GK <- read.csv('GK_declustered.csv')
data_GK <- data.frame(data_GK)
mag_GK <- data_GK$magnitude

data_GRU <- read.csv('GRU_declustered.csv')
data_GRU <- data.frame(data_GRU)
mag_GRU <- data_GRU$magnitude

data_URH <- read.csv('URH_declustered.csv')
data_URH <- data.frame(data_URH)
mag_URH <- data_URH$magnitude


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

fmd_gk <- fmd(mag_GK,bin)
fmd_gru <- fmd(mag_GRU, bin)
fmd_urh <- fmd(mag_URH, bin)
fmd_obs <- fmd(mag, bin)

grid()

plot(fmd_obs$m, fmd_obs$cum, log = "y", xlab = "Magnitude", ylab = "Number of Events", 
     type = "n", mgp = c(2, 0.8, 0), cex.lab=0.8, cex.axis=0.8)

points(fmd_obs$m, 
       fmd_obs$noncum, 
       col = "#32CD32", 
       pch = 17, 
       cex = 0.7, 
       type = 'o')

points(fmd_gk$m, 
       fmd_gk$noncum, 
       col = "blue", 
       pch = 17, 
       cex = 0.7, 
       type = 'o')

points(fmd_gru$m, 
       fmd_gru$noncum, 
       col = "red", 
       pch = 17, 
       cex = 0.7, 
       type = 'o')

points(fmd_urh$m, 
       fmd_urh$noncum, 
       col = "#DAA520", 
       pch = 17, 
       cex = 0.7, 
       type = 'o')


# Add a legend
legend("topright", 
       legend = c("Observed", "GardnerKnopoff", "Grünthal", "Uhrhammer"),
       col = c("#32CD32", "blue", "red", "#DAA520"), 
       pch = 17, 
       cex = 0.8,
       title = "Catalogs non-cum FMD", 
       lty = c(1,1,1,1))

# Add gridlines


maxc <- function(mag, mbin){
  FMD <- fmd(mag, mbin)
  Mc <- FMD$m[which(FMD$noncum == max(FMD$noncum))[1]]
  return(list(Mc = Mc))
}


maxc_GK <- maxc(mag_GK, bin)
maxc_GRU <- maxc(mag_GRU, bin)
maxc_URH <- maxc(mag_URH, bin)

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
result <- bootstrap_estimation(mag_GK, bin, nbsample)

nbsample <- 500
result <- bootstrap_estimation(mag_GRU, bin, nbsample)

nbsample <- 500
result <- bootstrap_estimation(mag_URH, bin, nbsample)










