# load data
data = read.csv(file = "Project-1/1_snow_particles.csv")
data$endpoint[52] = 3.25 # s.t. all the bins are of comparable size

# parameters obtained by optimizing binned log likelihood
mu1_opt = -2.019334
mu2_opt = -0.4550672
sigma1_opt = 0.6051886
sigma2_opt = 0.3024675
tau_opt = 0.6510166

# starting values for EM
mu1 <- -3
mu2 <- -0.5
sigma1 <- 0.3
sigma2 <- 0.4
tau <- 0.7

num_inbin = rep(0, length(data$X)) # store the number of datapoints in each bin here

for (i in 1:length(data$X)) {
  #rounded number of data in the i-th bin, need to jitter the data in the next step
  num_inbin[i] = round(data$retained....[i] / 100 * data$particles.detected[1]) 
}

# function to sample rom a bi-log-normal distribution
rbilnorm <- function(N, mu1, mu2, sigma1, sigma2, tau){
  ind <- I(runif(N) > tau)
  X <- rep(0,N)
  X[ind] <- rlnorm(sum(ind), mu1, sigma1)
  X[!ind] <- rlnorm(sum(!ind), mu2, sigma2)
  return(X)
}

# pdf for the bi-log-mormal distribution
dbilnorm <- function(x, mu1, mu2, sigma1, sigma2, tau){
  y <- (1-tau)*dlnorm(x, mu1, sigma1) + tau * dlnorm(x, mu2, sigma2)
  return(y)
}

# EM algorithm for the bi-log-normal model
EM = function(x, mu1, mu2, sigma1, sigma2, tau) {
  n = length(x)
  lold = 0
  lnew = 10
  while (abs(lold-lnew) > 10^-2) {
    #calculate log likelihood for convergence criteria
    lold = sum(log((1 - tau)*dlnorm(x, mean = mu1, sd = sigma1) + tau * dlnorm(x, mean = mu2, sd = sigma2)))
    p = dlnorm(x, mean = mu2, sd = sigma2) * tau / dbilnorm(x, mu1, mu2, sigma1, sigma2, tau)
    #update the parameters according to ML estimates
    tau = 1/n * sum(p)
    mu1 = sum(log(x) * (1 - p)) / sum(1 - p)
    mu2 = sum(log(x) * p) / sum(p)
    sigma1 = sqrt(sum((1 - p) * ((log(x) - mu1)^2)) / sum(1 - p))
    sigma2 = sqrt(sum(p * ((log(x) - mu2)^2)) / sum(p))
    #calculate log likelihood for convergence criteria
    lnew = sum(log((1 - tau)*dlnorm(x, mean = mu1, sd = sigma1) + tau * dlnorm(x, mean = mu2, sd = sigma2)))
  }
  return(c(mu1, mu2, sigma1, sigma2, tau))
}

# function calculating the negative log likelihood of the binned data (which we want to minimize)
negloglikelihood = function(par, data, num_inbin) {
  logli = 0 # start at 0 and sum op at each datapoint
  mu1 = par[1]
  mu2 = par[2]
  sigma1 = par[3]
  sigma2 = par[4]
  tau = par[5]
  for (i in 1:length(data$X)) {
    logli = logli + (num_inbin[i]) * log((cdfbilognorm(data$endpoint[i], mu1, mu2, sigma1, sigma2, tau) - 
                                            cdfbilognorm(data$startpoint[i], mu1, mu2, sigma1, sigma2, tau)))
  }
  return(-logli) #return the negative log likelihood
}

# function calculating the cdf of a bilognornal distribution<
cdfbilognorm = function(x, mu1, mu2, sigma1, sigma2, tau) {
  y <- (1 - tau) * plnorm(x, mu1, sigma1) + tau * plnorm(x, mu2, sigma2)
  return(y)
}

# initialize parametric bootstrap
set.seed(3) # for reproducability
n_bootstrap = 100 # number of bootstrap samples
ndata_boot = 705044 # number of data in bootstrap sample (same as original data)

# generate bootstrap datasets
bootstrap_data = matrix(nrow = ndata_boot, ncol = n_bootstrap) 
for (i in 1:n_bootstrap) {
  bootstrap_data[, i] = rbilnorm(ndata_boot, mu1_opt, mu2_opt, sigma1_opt, sigma2_opt, tau_opt)
}

#bin the bootstrap data
binned_bootstrap_data = matrix(nrow = length(data$X), ncol = n_bootstrap)
for (i in 1:n_bootstrap) {
  for (j in 1:length(data$X)) {
    counter = 0
    for (k in 1:ndata_boot) {
      if (data$startpoint[j] <= bootstrap_data[k, i] && data$endpoint[j] > bootstrap_data[k, i]) {
        counter = counter + 1
      }
    }
    binned_bootstrap_data[j, i] = counter # number of data points in j-th bin
  }
}

# jitter the binned bootstrap samples
jittered_boot_data = matrix(nrow = ndata_boot, ncol = n_bootstrap)
for (k in 1:n_bootstrap) {
  counter = 0
  for (i in 1:length(data$X)) {
    interval_length = data$endpoint[i] - data$startpoint[i]
    for(j in 1:binned_bootstrap_data[i, k]) {
      if (binned_bootstrap_data[i, k] == 0) {
        next
      }
      jittered_boot_data[j + counter, k] = data$startpoint[i] + 
        runif(1, min = 0, max = interval_length)
    }
    counter = counter + binned_bootstrap_data[i, k]
  }
}

# estimate the parameters for each bootstrap sample
estim_boot = matrix(nrow = 5, ncol = n_bootstrap)
for (i in 1:n_bootstrap) {
  estim_boot[, i] = EM(jittered_boot_data[, i], mu1, mu2, sigma1, sigma2, tau)
}

# use optimizetion for the binned bootstrap data
estim_boot_opt = matrix(nrow = 5, ncol = n_bootstrap)
for (i in 1:n_bootstrap) {
  out = optim(par = estim_boot[, i], fn = negloglikelihood, method="Nelder-Mead",
              data = data, num_inbin = binned_bootstrap_data[, i])
  estim_boot_opt[, i] = as.numeric(unlist(out[1]))
}

# Calculate the kolmogorov smirnoff statisitc
# To calculate it we create data in such a way to create ecdf for binned data
# That means we will put each datapoint at the end of its bin
counter = 0
n = sum(num_inbin)
ecdf_data = rep(0, times = n)
for (i in 1:length(num_inbin)) {
  for(j in 1:num_inbin[i]) {
    if (num_inbin[i] == 0) {
      next
    }
    ecdf_data[j + counter] = data$endpoint[i]
  }
  counter = counter + num_inbin[i]
}
ecdf_binned = ecdf(ecdf_data) # now calculate ecdf using the prepared data
grid = seq(from = 0, to = 3.25, by = 0.01) # grid to compare distribution functions on (for KS statistic)
kol_smir = 0 
for (i in 1:length(grid)) {
  a = abs(ecdf_binned(grid[i]) - 
            cdfbilognorm(grid[i], mu1_opt, mu2_opt, sigma1_opt, sigma2_opt, tau_opt))
  if (a > kol_smir) {
    kol_smir = a
  }
}

# calculate T_b^* (similar to KS statistic)
T_b = rep(0, n_bootstrap)
for (k in 1:n_bootstrap) {
  counter = 0
  ecdf_boot_data = rep(0, ndata_boot)
  for (i in 1:length(data$X)) {
    for(j in 1:binned_bootstrap_data[i, k]) {
      if (binned_bootstrap_data[i, k] == 0) {
        next
      }
      ecdf_boot_data[j + counter] = data$endpoint[i] 
    }
    counter = counter + binned_bootstrap_data[i, k]
  }
  ecdf_boot = ecdf(ecdf_boot_data)
  for (j in 1:length(grid)) {
    a = abs(ecdf_boot(grid[j]) - 
              cdfbilognorm(grid[j], estim_boot_opt[1, k], estim_boot_opt[2, k], 
                           estim_boot_opt[3, k], estim_boot_opt[4, k], estim_boot_opt[5, k]))
    if (a > T_b[k]) {
      T_b[k] = a
    }
  }
}

# calculate the estimated p value
a = 0
for (k in 1:n_bootstrap) {
  if (T_b[k] > kol_smir) {
    a = a + 1
  }
}
p = 1 / (n_bootstrap + 1) * (a + 1)
print(p)