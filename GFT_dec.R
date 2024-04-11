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

#Goodness-of-fit test (GFT) [Wiemer & Wyss, 2000]
gft <- function(mag,mbin){
  FMD <- fmd(mag,mbin)
  McBound <- maxc(mag,mbin)$Mc
  Mco <- McBound-0.4+(seq(15)-1)/10
  R <- numeric(15)
  for(i in 1:15){
    indmag <- which(mag >= Mco[i]-mbin/2)
    b <- log10(exp(1))/(mean(mag[indmag])-(Mco[i]-mbin/2))
    a <- log10(length(indmag))+b*Mco[i]
    FMDcum_model <- 10^(a-b*FMD$m)
    indmi <- which(FMD$m >= Mco[i])
    R[i] <- sum(abs(FMD$cum[indmi]-FMDcum_model[indmi]))/sum(FMD$cum[indmi])*100
    #in Wiemer&Wyss [2000]: 100-R
  }
  indGFT <- which(R <= 5) #95% confidence
  if(length(indGFT) != 0){
    Mc <- Mco[indGFT[1]]
    best <- "95%"
  } else{
    indGFT <- which(R <= 10) #90% confidence
    if(length(indGFT) != 0){
      Mc <- Mco[indGFT[1]]
      best <- "90%"
    } else{
      Mc <- McBound
      best <- "MAXC"
    }
  }
  
  return(list(Mc=Mc, best=best, Mco=Mco, R=R))
}


gft_GK <- gft(mag_GK, bin)
gft_GRU <- gft(mag_GRU, bin)
gft_URH <- gft(mag_URH, bin)


bootstrap_estimation <- function(data, bin, nbsample) {
  set.seed(123)
  
  Mc_bootstrap <- numeric(nbsample)
  
  for (i in 1:nbsample) {
    resampled_data <- sample(data, replace = TRUE)
    gft_result <- gft(resampled_data, bin)
    Mc_bootstrap[i] <- gft_result$Mc
    
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

# Example usage

nbsample <- 500
result <- bootstrap_estimation(mag_GK, bin, nbsample)

nbsample <- 500
result <- bootstrap_estimation(mag_GRU, bin, nbsample)

nbsample <- 500
result <- bootstrap_estimation(mag_URH, bin, nbsample)






