citation <- "Earthquake Frequency-Magnitude Distribution & Network Statistics"
author <- "A.Mignan"
link <- "https://github.com/amignan/rseismNet"
print(citation)
print(author)
print(link)

library(devtools)
library(rseismNet)

mus <- c(1, 2, 3)  # Define different values for mu
sigmas <- seq(0.1, 3.0, by = 0.1)  # Define a sequence for sigma values

list_params_curved <- list()

for (mu in mus) {
  for (sigma in sigmas) {
    list_params_curved[[length(list_params_curved) + 1]] <- list(beta = log(10), mu = mu, sigma = sigma)
  }
}

# print(list_params_curved)

# list_params_curved <- list(list(beta = log(10), mu = 1.06, sigma = 0.33))

run_simulations_for_curved_sets <- function(n_eq, list_params_curved, legend_cex = 1.3, n_sample = 500, namePlotDir = "dir") {
  
  plots_dir <- namePlotDir
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir)
  }
  
  for (i in seq_along(list_params_curved)) {
    params_curved <- list_params_curved[[i]]
    
    m_curved <- rseismNet::bfmd.sim(n_eq, params_curved)
    mdistr_curved <- rseismNet::fmd(m_curved)
    
    mc_mode_curved <- rseismNet::mc.val(m_curved, "mode")
    mc_mbass_curved <- rseismNet::mc.val(m_curved, "mbass")
    mc_gft_curved <- rseismNet::mc.val(m_curved, "gft")
    
    beta_mode_curved <- rseismNet::beta.mle(m_curved, mc_mode_curved)
    beta_mbass_curved <- rseismNet::beta.mle(m_curved, mc_mbass_curved)
    beta_gft_curved <- rseismNet::beta.mle(m_curved, mc_gft_curved)
    
    file_name <- sprintf("%s/simulation_plot_%d.png", plots_dir, i)
    
    png(filename = file_name, width = 1500, height = 800)
    
    par(mfrow = c(1, 2), oma = c(5, 5, 5, 5), mar = c(5, 5, 4, 2) + 0.1, bg = "white")
    
    # title_curved <- paste("Curved Set", i, 'n.event', n_eq)
    plot(mdistr_curved$mi, mdistr_curved$Ni, log = "y", type="p",
         # main = title_curved, 
         xlab = "Magnitude", ylab = "Log N",
         pch = 21, col =  "#FC4E07", cex.main = 2, cex.lab = 1.8, cex.axis = 1.5)
    points(mdistr_curved$mi, mdistr_curved$ni, type = "p", pch = 24, col =  "#0073C2FF")
    
    abline(v = mc_mode_curved, col = "#FF00CC")
    abline(v = mc_mbass_curved, col = "#999500", lty = "dashed")
    abline(v = mc_gft_curved, col = "#FF9900", lty = "dotdash")
    
    safe_plot_line <- function(mc_value, beta_value, mdistr_curved, col, lty) {
      idx <- which(mdistr_curved$mi >= mc_value)
      if (length(idx) > 0 && mdistr_curved$ni[idx[1]] > 0 && !is.na(beta_value) && !is.infinite(beta_value)) {
        intercept <- log10(mdistr_curved$ni[idx[1]]) + beta_value / log(10) * mc_value
        slope <- -beta_value / log(10)
        if (is.finite(intercept) && is.finite(slope)) {
          abline(a = intercept, b = slope, col = col, lty = lty)
        }
      }
    }
    
    safe_plot_line(mc_mode_curved, beta_mode_curved, mdistr_curved, "#FF00CC", 1)
    safe_plot_line(mc_mbass_curved, beta_mbass_curved, mdistr_curved, "#999500", "dashed")
    safe_plot_line(mc_gft_curved, beta_gft_curved, mdistr_curved, "#FF9900", "dotdash")
    
    legend_info_curved <- c(paste("n. events :", n_eq),
                            paste("\u03B2 :", format(params_curved$beta, digits=4, nsmall=2)),
                            paste("\u03BC :", format(params_curved$mu, digits=4, nsmall=2)),
                            paste("\u03C3 :", format(params_curved$sigma, digits=4, nsmall=2)),
                            paste("Mode:", format(mc_mode_curved, digits=4, nsmall=2)),
                            paste("Mbass:", format(mc_mbass_curved, digits=4, nsmall=2)),
                            paste("Gft:", format(mc_gft_curved, digits=4, nsmall=2)))
    
    legend_colors_curved <- c("black", "black", "black", "#FF00CC", "#999500", "#FF9900")
    legend_lty_curved <- c(NA, NA, NA, 1, 2, 3)
    
    legend("topright", legend=legend_info_curved, col=legend_colors_curved, lty=legend_lty_curved, cex=1.3, bg='white')
    
    mc_mode_bootstrap <- sapply(1:n_sample, function(i)
      rseismNet::mc.val(sample(m_curved, replace = T), "mode"))
    mean(mc_mode_bootstrap, na.rm = T)
    sd(mc_mode_bootstrap, na.rm = T)
    
    mc_mbass_bootstrap <- sapply(1:n_sample, function(i) 
      rseismNet::mc.val(sample(m_curved, replace = T), "mbass"))
    mean(mc_mbass_bootstrap, na.rm = T)
    sd(mc_mbass_bootstrap, na.rm = T)
    
    mc_gft_bootstrap <- sapply(1:n_sample, function(i) 
      rseismNet::mc.val(sample(m_curved, replace = T), "gft"))
    mean(mc_gft_bootstrap, na.rm = T)
    sd(mc_gft_bootstrap, na.rm = T)
    
    stddev <- seq(0,3,0.1)
    
    mci_mode <- round(mean(mc_mode_bootstrap, na.rm = T) + stddev * sd(mc_mode_bootstrap, na.rm = T), digits = 1)
    beta_var_mode <- sapply(1:length(stddev), function(i) rseismNet::beta.mle(m_curved, mci_mode[i]))
    
    mci_mbass <- round(mean(mc_mbass_bootstrap, na.rm = T) + stddev * sd(mc_mbass_bootstrap, na.rm = T), digits = 1)
    beta_var_mbass <- sapply(1:length(stddev), function(i) rseismNet::beta.mle(m_curved, mci_mbass[i]))
    
    mci_gft <- round(mean(mc_gft_bootstrap, na.rm = T) + stddev * sd(mc_gft_bootstrap, na.rm = T), digits = 1)
    beta_var_gft <- sapply(1:length(stddev), function(i) rseismNet::beta.mle(m_curved, mci_gft[i]))
    
    # title_betasense <- paste(legend_info_curved[2],legend_info_curved[3], legend_info_curved[4], legend_info_curved[5], legend_info_curved[6])
    plot(stddev, beta_var_mbass, type = "l", col = "#999500", xlab = "\u03C3", ylab = "\u03B2", ylim = c(log(10) - 0.5, log(10) + 0.3),# main = title_betasense,
         cex.main = 2, cex.lab = 1.8, cex.axis = 1.5)
    lines(stddev, beta_var_gft, col = "#FF9900")
    abline(h = log(10))
    legend_info_second <- c("Mbass", "Gft")
    legend_colors_second <- c("#999500", "#FF9900")
    legend("topright", legend=legend_info_second, col=legend_colors_second, lty=c(1,1), cex=1.3, bg='white')
    
    dev.off() 
  }
  
  par(mfrow = c(1, 1)) 
}


n_eq5000 <- 0.5e4
n_eq10000 <- 1e4
n_eq20000 <- 2e4
n_eq25000 <- 2.5e4

# Execute the simulations with the defined parameter sets
run_simulations_for_curved_sets(n_eq5000, list_params_curved)
run_simulations_for_curved_sets(n_eq10000, list_params_curved)
run_simulations_for_curved_sets(n_eq20000, list_params_curved)
run_simulations_for_curved_sets(n_eq25000, list_params_curved)

























