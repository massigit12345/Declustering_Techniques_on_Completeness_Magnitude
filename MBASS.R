citation <- "Amorèse D., 2007, Applying a Change-Point Detection Method on Frequency-Magnitude Distributions, Bull. Seism. Soc. Am. 97, no 5, 1742-1749. October 2007Bulletin of the Seismological Society of America 97(5):1742-1749 DOI:10.1785/0120060181" 
print(citation)

library(bootstrap)

data <- read.csv('Rcatalogo.csv')
data <- data.frame(data)
bin <- 0.1
cat <- data
mag <- cat$mag

data_GK <- read.csv('GK_declustered.csv')
data_GK <- data.frame(data_GK)
mag_GK <- data_GK$magnitude

data_GRU <- read.csv('GRU_declustered.csv')
data_GRU <- data.frame(data_GRU)
mag_GRU <- data_GRU$magnitude

data_URH <- read.csv('URH_declustered.csv')
data_URH <- data.frame(data_URH)
mag_URH <- data_URH$magnitude

fmbass <- function(a, delta = 0.1, plot = TRUE, alldisc = FALSE) {
  if (plot) {
    par(mfrow = c(1, 1))
  }
  tau <- numeric()
  pva <- numeric()
  minmag <- min(a, na.rm = TRUE)
  g_r <- hist(a, plot = FALSE, breaks = seq((minmag - delta/2), (max(a, na.rm = TRUE) + delta/2), delta))
  n <- length(g_r$counts)
  xc <- seq(minmag, max(a, na.rm = TRUE), delta)[1:(n-1)]
  log_nc <- log10((1/delta) * (length(a) - cumsum(g_r$counts)[1:(n-1)]) * delta)
  x <- seq(minmag, max(a, na.rm = TRUE), delta)
  log_n <- log10((1/delta) * g_r$counts * delta)
  x <- x[is.finite(log_n)]
  log_n <- log_n[is.finite(log_n)]
  sl <- diff(log_n) / diff(x)
  xsl <- x[2:length(x)]
  if (plot) {
    plot(xc, 10^log_nc, type = "p", ylim = c(1, length(a)), log = "y", xlab = "Magnitude", ylab = "Number of events", pch = 1)
    points(x, 10^log_n, pch = 2)
  }
  niter <- 3
  N <- length(sl)
  j <- 0
  k <- 0
  SA <- vector(length = N)
  while (j < niter) {
    for (i in seq(1, N, 1)) SA[i] <- abs(2 * sum(rank(sl)[1:i]) - i * (N + 1))
    n1 <- which(SA == SA[order(SA)[length(order(SA))]])
    xn1 <- sl[1:n1[1]]
    xn2 <- sl[-(1:n1[1])]
    if ((n1[1] > 2) && (n1[1] <= (N-2)) && (wilcox.test(xn1, xn2, exact = FALSE, correct = TRUE)[3] < 0.05)) {
      k <- k + 1
      pva[k] <- wilcox.test(xn1, xn2, exact = FALSE, correct = TRUE)[3]
      tau[k] <- n1[1]
      if (k > 1) {
        medsl1 <- median(sl[1:n0])
        medsl2 <- median(sl[-(1:n0)])
        sl[1:n0] <- sl[1:n0] + medsl1
        sl[(n0+1):length(sl)] <- sl[(n0+1):length(sl)] + medsl2
      }
      medsl1 <- median(sl[1:n1[1]])
      medsl2 <- median(sl[-(1:n1[1])])
      sl[1:n1[1]] <- sl[1:n1[1]] - medsl1
      sl[(n1[1]+1):length(sl)] <- sl[(n1[1]+1):length(sl)] - medsl2
      n0 <- n1[1]
    }
    j <- j + 1
  }
  v_pva <- as.vector(pva, mode = "numeric")
  ip <- order(v_pva)
  m0 <- c(signif(xsl[tau[ip[1]]]), signif(xsl[tau[ip[2]]]))
  if (alldisc) {
    return((list(discmag = xsl[tau], p = v_pva, m0 = m0)))
  }
  return(m0)
}

mbass <- function(a, delta = 0.1, plot = TRUE, alldisc = TRUE, bs = 0) {
  mba <- function(x) {
    fmbass(x, delta, plot, alldisc)
  }
  if (bs == 0) {
    res <- mba(a)  # Actual FMD analyzed
    return(res)
  } else {
    results <- bootstrap(a, bs, mba)  # Run bootstrap with bs replicates
    m0_values <- numeric(bs)  # Store m0 values for each bootstrap replicate
    
    for (i in 1:bs) {
      # Extract m0 values from the bootstrap results
      m0_values[i] <- results[["thetastar"]][[i]][["m0"]][1]
    }
    # Calculate and print mean and standard deviation of m0 values
    mean_m0 <- mean(m0_values)
    sd_m0 <- sd(m0_values)
    
    cat("Mean of m0: ", mean_m0, "\n")
    cat("Standard Deviation of m0: ", sd_m0, "\n")
    
    return(list(mean_m0 = mean_m0, sd_m0 = sd_m0, m0s = m0_values))
  }
}


fmbass.var <- fmbass(mag, alldisc = FALSE, plot=FALSE)
mbass.mc <- fmbass(mag, alldisc= FALSE, plot=FALSE)[1]

mbass.var <- mbass(mag, delta = 0.1, plot = TRUE, alldisc = TRUE, bs = 500)


# Plot occurrences of main discontinuity identified as Mc 
m0_freq <- table(mbass.var$m0s)

frequency_df <- as.data.frame(m0_freq)
colnames(frequency_df) <- c("Values", "Occurrences")

frequency_df <- frequency_df[order(frequency_df$Values), ]

barplot(frequency_df$Occurrences, names.arg = frequency_df$Values, xlab = "Main Discontinuity", ylab = "Occurrences", main = "Distribution of m0 Values")

# Repeat the same procedure for declustered catalog