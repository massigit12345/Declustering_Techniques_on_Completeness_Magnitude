library(nlstools)
library(dplyr)


data <- read.csv('Rcatalogo.csv')
data <- data.frame(data)
bin <- 0.1
cat <- data
mag <- cat$mag


citation <- "Mignan, A., J. Woessner (2012). Estimating the magnitude of completeness for earthquake catalogs, Community Online Resource for Statistical Seismicity Analysis, doi:10.5078/corssa-00180805. Available at http://www.corssa.org."
print(citation)


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
  mu <- abs(McMAXC/2); sig <- abs(McMAXC/4)
  if(mu > 1)mu <- abs(McMAXC/10); sig <- abs(McMAXC/20)
  McBound <- McMAXC
  Mco <- McBound-0.3+(seq(9)-1)/10
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
                  start=list(mean=mu, sd=sig), control=list(maxiter=100, warnOnly = TRUE))
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
              model=savedmodel[indbestfit,], Mco=Mco, prob=prob, indbestfit=indbestfit)
  return(res)
}

emr.variable <- emr(mag,bin)

FMD <- fmd(mag,bin)

# Plot the FMD
plot(FMD$m, 
     emr.variable$model, 
     log="y", 
     xlab="Magnitude", 
     ylab="Number of events", 
     pch=16,
     mgp = c(2, 0.8, 0), 
     cex.lab=1, 
     cex.axis=1, 
     col='red', 
     type = "n")

points(FMD$m, 
       emr.variable$model, 
       col='red', 
       type = 'o')

points(FMD$m, 
       FMD$noncum, 
       type='o', 
       pch=16)


legend("right", legend = c("EMR Model", "Observed FMD"), 
       lty = c(1, 1), 
       pch = c(1, 16), 
       col = c("red", "black"), 
       cex=1)

legend("topright", c(paste("a =", round(emr.variable$a, digit=2)),
                     paste("b =", round(emr.variable$b, digit=2)), 
                     paste("mu =", round(emr.variable$mu, digit=2)), 
                     paste("sigma =", round(emr.variable$sigma, digit=2)),
                     paste("Mc =", round(emr.variable$Mc, digit=2))),
                     cex=1.2, 
                     title = "EMR model parameter")

abline(v=emr.variable$Mc)



# Create Vector Mco
Mco <- emr.variable$Mco

# Create Vector Prob
prob <- round(emr.variable$prob, 
              digits = 4)

# Plotting the MLE vs Mco
plot(Mco, prob, type = "b", 
     xlab = "Magnitude cut-off", 
     ylab = "MLE",
     mgp = c(2, 0.8, 0), 
     cex.lab=0.8, 
     cex.axis=0.8)

grid()

# Adding a line at Mc
abline(v = emr.variable$Mc, 
       col = "red", 
       lty = 2, 
       bty = "n")

# Adding labels for the magnitude bins and probabilities
text(Mco, prob, 
     labels = round(prob, digit=2), 
     pos = 3, 
     offset = c(1,1), 
     cex=0.6 )




bootstrap_estimation <- function(data, bin, nbsample) {
  # set.seed(123)
  
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
    Mean = c(mean(Mc_bootstrap), 
             mean(b_bootstrap), 
             mean(a_bootstrap), 
             mean(mu_bootstrap), 
             mean(sigma_bootstrap)),
    
    Dev.std = c(sd(Mc_bootstrap), 
                sd(b_bootstrap), 
                sd(a_bootstrap), 
                sd(mu_bootstrap), 
                sd(sigma_bootstrap))
  )
  
  # print(result_table)
  
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
result <- bootstrap_estimation(mag, bin, nbsample)










