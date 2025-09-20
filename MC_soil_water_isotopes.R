rm(list = ls())
set.seed(42)
nsyth = 1e6 # number of iterations
# load packages
pacman::p_load(tidyverse, readxl)

# load data 
dat1 = read_xlsx("data/Triple_Clumped_Combined_Stats.xlsx", sheet = "triple_clumped")
dat2 = read_xlsx("data/Triple_Clumped_Combined_Stats.xlsx", sheet = "clumped")

# calculating soil water d18O and Dp17O
sims1 = dat1[,1]
for (i in 1:nrow(dat1)) {
  d47 = rnorm(nsyth, dat1$D47[i], dat1$D47_sd[i])
  T47 = sqrt(3.91e-2 * 1e6 / (d47 - .154)) - 273.15
  alpha18 = exp((1.61e4 / (T47 + 273.15) - 24.6) / 1e3)
  theta = -1.39 / (273.15 + T47) + 0.5305
  alpha17 = alpha18 ^ theta 
  d18_c = rnorm(nsyth, dat1$d18c[i], dat1$d18c_sd[i])
  d17_c = rnorm(nsyth, dat1$d17c[i], dat1$d17c_sd[i])
  d18_sw = (d18_c + 1e3) / alpha18 - 1e3
  d17_sw = (d17_c + 1e3) / alpha17 - 1e3
  dp18_sw = 1e3 * log(d18_sw / 1e3 + 1)
  dp17_sw = 1e3 * log(d17_sw / 1e3 + 1)
  Dp17_sw = (dp17_sw - 0.528 * dp18_sw) * 1e3
  sims1$T47[i] = mean(T47)
  sims1$T47_sd[i] = sd(T47)
  sims1$d18sw[i] = mean(d18_sw)
  sims1$d18sw_sd[i] = sd(d18_sw)
  sims1$dp18sw[i] = mean(dp18_sw)
  sims1$dp18sw_sd[i] = sd(dp18_sw)
  sims1$Dp17sw[i] = mean(Dp17_sw)
  sims1$Dp17sw_sd[i] = sd(Dp17_sw)
}

# calculating soil water d18O
sims2 = dat2[,2]
for (i in 1:nrow(dat2)) {
  d47 = rnorm(nsyth, dat2$D47[i], dat2$D47_sd[i])
  T47 = sqrt(3.91e-2 * 1e6 / (d47 - .154)) - 273.15
  alpha18 = exp((1.61e4 / (T47 + 273.15) - 24.6) / 1e3)
  d18c_vpdb = rnorm(nsyth, dat2$d18c[i], dat2$d18c_sd[i])
  d18c_smow = 1.03091 * d18c_vpdb + 30.91
  d18_sw = (d18c_smow + 1e3) / alpha18 - 1e3
  sims2$T47[i] = mean(T47)
  sims2$T47_sd[i] = sd(T47)
  sims2$d18sw[i] = mean(d18_sw)
  sims2$d18sw_sd[i] = sd(d18_sw)
}

write_csv(sims1, "out/Guadix_Dp17sw.csv")
write_csv(sims2, "out/Guadix_d18sw.csv")
