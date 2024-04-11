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

maxc <- function(mag,mbin){
  FMD <- fmd(mag,mbin)
  Mc <- FMD$m[which(FMD$noncum == max(FMD$noncum))[1]]
  return(list(Mc=Mc))
}

# Modification to Mco initial range values
mbs <- function(mag, mbin) {
  McBound <- maxc(mag, mbin)$Mc
  Mco <- McBound - 0.7 + (seq(40) - 1) / 10  # Modify the sequence range for Mco
  bi <- numeric(40)
  unc <- numeric(40)
  
  for (i in 1:40) {
    indmag <- which(mag > Mco[i] - mbin / 2)
    nbev <- length(indmag)
    bi[i] <- log10(exp(1)) / (mean(mag[indmag]) - (Mco[i] - mbin / 2))
    unc[i] <- 2.3 * bi[i]^2 * sqrt(sum((mag[indmag] - mean(mag[indmag]))^2) / (nbev * (nbev - 1)))
  }
  
  bave <- numeric(40)  # Adjusted range for bave
  for (i in 1:40) bave[i] <- mean(bi[i:(i + 5)])
  dbi_old <- abs(diff(bi))
  indMBS_old <- which(dbi_old <= 0.03)
  dbi <- abs(bave[1:40] - bi[1:40])  # Adjusted range for dbi
  indMBS <- which(dbi <= unc[1:40])  # Adjusted range for unc
  Mc <- Mco[indMBS[1]]
  
  return(list(Mc = Mc, Mco = Mco, 
              bi = bi, unc = unc, 
              bave = bave, indMBS = indMBS))
}

mbs_GK <- mbs(mag_GK, bin)
mbs_GRU <- mbs(mag_GRU, bin)
mbs_URH <- mbs(mag_URH, bin)


bootstrap_estimation <- function(data, bin, nbsample) {
  set.seed(123)
  
  Mc_bootstrap <- numeric(nbsample)
  
  for (i in 1:nbsample) {
    resampled_data <- sample(data, replace = TRUE)
    mbs_var <- mbs(resampled_data, bin)
    Mc_bootstrap[i] <- mbs_var$Mc
  }
  
  result_table <- data.frame(
    Parameter = c('Mc', 'b', 'a'),
    Mean = c(mean(Mc_bootstrap, na.rm=TRUE)),
    Dev.std = c(sd(Mc_bootstrap,na.rm=TRUE))
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







