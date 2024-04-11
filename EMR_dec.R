citation <- "Mignan, A., J. Woessner (2012). Estimating the magnitude of completeness for earthquake catalogs, Community Online Resource for Statistical Seismicity Analysis, doi:10.5078/corssa-00180805. Available at http://www.corssa.org."
print(citation)

library(nlstools)
library(dplyr)

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


emr <- function(mag,mbin){
  FMD <- fmd(mag,mbin)
  nbm <- length(FMD$m)
  McMAXC <- maxc(mag,mbin)$Mc
  mu <- abs(McMAXC/4); 
  sig <- abs(McMAXC/2)
  if(mu > 1)mu <- abs(McMAXC/10); sig <- abs(McMAXC/20)
  print(mu)
  print(sig)
  McBound <- McMAXC
  Mco <- McBound-0.5+(seq(9)-1)/10
  params <- numeric(9*4); dim(params) <- c(9,4) #a, b, mu, sigma
  prob <- numeric(9)
  savedmodel <- numeric(9*nbm); dim(savedmodel) <- c(9,nbm)
  for(i in 1:9){
    indmag <- which(mag > Mco[i]-mbin/2)
    nbev <- length(indmag)
    b <- log10(exp(1))/(mean(mag[indmag])-(Mco[i]-mbin/2))
    a <- log10(length(indmag))+b*Mco[i]
    cumN <- 10^(a-b*FMD$m)
    params[i,1] <- a; params[i,2] <- b
    cumNtmp <- 10^(a-b*(max(FMD$m)+mbin))
    cumNtmp <- c(cumN, cumNtmp)
    N <- abs(diff(cumNtmp))
    data <- data.frame(N=N, m=FMD$m, Nd=FMD$noncum)
    indLow <- which(FMD$m < Mco[i]); indHigh <- which(FMD$m >= Mco[i])
    dataTest <- data.frame(N=data$N[indLow], m=data$m[indLow], Nd=data$Nd[indLow])
    dataTmp <- data.frame(N=data$N[indHigh], m=data$m[indHigh], Nd=data$Nd[indHigh])
    checkNo0 <- which(dataTest$Nd != 0)
    dataTest <- data.frame(N=dataTest$N[checkNo0], m=dataTest$m[checkNo0],
                           Nd=dataTest$Nd[checkNo0])
    #Nmax <- max(dataTmp$Nd)
    Nmax <- max(dataTest$Nd)
    #Nmax <- dataTest$Nd[length(dataTest$Nd)]
    Mmintmp <- min(dataTest$m)
    dataTest$Nd <- dataTest$Nd/Nmax
    dataTest$m <- dataTest$m-Mmintmp
    data4fit <- data.frame(N=dataTest$Nd, m=dataTest$m)
    #non-linear least squares fit
    nlsfit <- nls(N~pnorm(m, mean=mean, sd=sd), data=data4fit,
                  start=list(mean=mu, sd=sig), control=list(maxiter=100,
                                                            warnOnly = TRUE))
    params[i,3] <- coef(nlsfit)["mean"]; params[i,4] <- coef(nlsfit)["sd"]
    dataTest$N <- pnorm(dataTest$m, mean=coef(nlsfit)["mean"],
                        sd=coef(nlsfit)["sd"])*Nmax
    dataTest$m <- dataTest$m+Mmintmp
    dataTest$Nd <- dataTest$Nd*Nmax
    dataPred <- data.frame(N=c(dataTest$N, dataTmp$N), m=c(dataTest$m, dataTmp$m),
                           Nd=c(dataTest$Nd, dataTmp$Nd))
    dataPred$N <- round(dataPred$N)
    savedmodel[i,c(checkNo0,indHigh)] <- dataPred$N
    #CHECK EMR METHOD#
    #pdf(paste(wd,"plot_NonCumModel_",Mco[i],".pdf", sep=""))
    #plot(dataPred$m, dataPred$Nd, pch=18, xlab="Magnitude",ylab="Cumulative Number")#, log="y")
    #points(dataPred$m, dataPred$N, pch=1)
    #abline(v=Mco[i], lty="dashed")
    #legend("topright", c("Data","EMR model"), cex=0.8, lty=c(0,0), pch=c(18,1))
    #dev.off()
    #write.table(dataPred, file=paste("file_NonCumModel_",Mco[i],".txt", sep=""))
    #Logarithm to the basis of 10 of Poisson probability density
    probtmp <- numeric(nbm)
    CheckNo0 <- which(dataPred$N != 0)
    Pmodel <- dataPred$N[CheckNo0]; Pdata <- dataPred$Nd[CheckNo0]
    probtmp[CheckNo0] <- 1/log(10)*(-Pmodel+Pdata*log(Pmodel)-lgamma(Pdata+1))
    prob[i] <- -sum(probtmp)
  }
  indbestfit <- which(prob == min(prob, na.rm=TRUE))
  res <- list(Mc=Mco[indbestfit], a=params[indbestfit,1], b=params[indbestfit,2],
              mu=params[indbestfit,3], sigma=params[indbestfit,4],
              model=savedmodel[indbestfit,], Mco=Mco, prob=prob)
  return(res)
}

emr_GK <- emr(mag_GK, bin)
emr_GRU <- emr(mag_GRU, bin)
emr_URH <- emr(mag_URH, bin)



bootstrap_estimation <- function(data, bin, nbsample) {
  set.seed(123)
  
  Mc_bootstrap <- numeric(nbsample)
  b_bootstrap <- numeric(nbsample)
  a_bootstrap <- numeric(nbsample)
  mu_bootstrap <- numeric(nbsample)
  sigma_bootstrap <- numeric(nbsample)
  
  for (i in 1:nbsample) {
    resampled_data <- sample(data, replace = TRUE)
    emr_result <- emr(resampled_data, bin)
    Mc_bootstrap[i] <- emr_result$Mc
    a_bootstrap[i] <- emr_result$a
    b_bootstrap[i] <- emr_result$b
    mu_bootstrap[i] <- emr_result$mu
    sigma_bootstrap[i] <- emr_result$sigma
  }
  
  result_table <- data.frame(
    Parameter = c('Mc', 'b', 'a', 'mu', 'sigma'),
    Mean = c(mean(Mc_bootstrap,na.rm=TRUE), 
             mean(b_bootstrap,na.rm=TRUE), 
             mean(a_bootstrap,na.rm=TRUE), 
             mean(mu_bootstrap,na.rm=TRUE), 
             mean(sigma_bootstrap,na.rm=TRUE)),
    Dev.std = c(sd(Mc_bootstrap,na.rm=TRUE), 
                sd(b_bootstrap,na.rm=TRUE), 
                sd(a_bootstrap,na.rm=TRUE),
                sd(mu_bootstrap,na.rm=TRUE),
                sd(sigma_bootstrap,na.rm=TRUE))
  )
  
  print(result_table)
  
  return(
    list(
      Mc_bootstrap = Mc_bootstrap,
      b_bootstrap = b_bootstrap,
      a_bootstrap = a_bootstrap,
      mu_bootstrap = mu_bootstrap,
      sigma_bootstrap = sigma_bootstrap
    )
  )
}

nbsample <- 500
result <- bootstrap_estimation(mag_GK, bin, nbsample)

nbsample <- 500
result <- bootstrap_estimation(mag_GRU, bin, nbsample)

nbsample <- 500
result <- bootstrap_estimation(mag_URH, bin, nbsample)




