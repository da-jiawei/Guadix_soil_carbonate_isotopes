# Set up time series ----
ai = function(){
  source("code/constructors.R")
  ## Read data
  td = read.csv("data/BBNP_data.csv")
  md = read.csv("data/Shatsky_data.csv")
  
  ## Protect against overly optimistic uncertainties
  td$d18O.stdev = sqrt(0.1^2 + td$d18O.stdev^2)
  td$d13C.stdev = sqrt(0.1^2 + td$d13C.stdev^2)
  
  ## Remove PETM and earliest terrestrial data, make ages negative
  md = md[md$age <= 55.741 | md$age >= 55.942,]
  td = td[td$Age <= 59 & td$Age >= 53, ]
  md$age = -md$age
  td$Age = -td$Age
  
  ## Parse data into series
  ### Marine
  d18Of = na.exclude(md[c("age", "d18O")])
  d18Of$d18O.stdev = rep(0.05)
  d13Cf = na.exclude(md[c("age", "d13C")])
  d13Cf$d13C.stdev = rep(0.05)
  mgcaf = na.exclude(md[c("age", "MgCa")])
  mgcaf$MgCa.stdev = rep(0.1)
  d11BGrub = d11BTsac = na.exclude(md[c("age", "d11B", "d11Bse", "species")])
  d11BGrub = d11BGrub[d11BGrub$species == "Grub", -4]
  d11BTsac = d11BTsac[d11BTsac$species == "Tsac", -4]
  
  ### Terrestrial
  d13Cc = na.exclude(td[c("Age", "d13C", "d13C.stdev")])
  d18Oc = na.exclude(td[c("Age", "d18O", "d18O.stdev")])
  D47c = na.exclude(td[c("Age", "D47", "D47.stderr")])
  
  ages = ts(d18Of$age, d13Cf$age, mgcaf$age, d11BGrub$age, d11BTsac$age,
            d13Cc$Age, d18Oc$Age, D47c$Age,
            c(d18Of$age, d13Cf$age, mgcaf$age, d11BGrub$age, d11BTsac$age),
            c(d13Cc$Age, d18Oc$Age, D47c$Age),
            seq(-59, -53, by = 0.2))
}


plot.jpi = function(x, y, n = 500, ylab = deparse(substitute(y)), ylim = range(y), ...){
  # x = timeseries ages, y = simslist matrix
  inds = sample(seq_along(y[, 1]), n)
  plot(x, y[inds[1],], type = "l", xlab = "Age (Ma)", ylab = ylab, ylim = ylim, 
       col = rgb(0, 0, 0, 0.05), ...)
  for(i in inds[-1]){
    lines(x, y[i,], col = rgb(0, 0, 0, 0.1))
  }
  
  m = apply(y, 2, median)
  points(x, m, pch = 21, cex = 1.2, lwd = 2, col = "red", bg = "white")
  # points(ai, rep(min(y), length(ai)))
}

get.ind = function(obs, ages){
  ai = sapply(obs, function(x){max(which(x >= ages))})
}

# Add time series to plot w 2 prob density envelopes
tsdens = function(d, base = "black"){
  # Check dimensions of d
  if(ncol(d) != 6){stop("d cols should be should be time, 5%, 25%, 50%, 75%, 95% CI")}
  
  base.rgb = col2rgb(base)
  cols = c(rgb(base.rgb[1]/255, base.rgb[2]/255, base.rgb[3]/255, alpha = 0.25), 
           rgb(base.rgb[1]/255, base.rgb[2]/255, base.rgb[3]/255, alpha = 0.25),
           rgb(base.rgb[1]/255, base.rgb[2]/255, base.rgb[3]/255, alpha = 1))
  
  polygon(c(d[, 1], rev(d[, 1])), c(d[, 2], rev(d[, 6])), 
          col = cols[1], border = NA)
  polygon(c(d[, 1], rev(d[, 1])), c(d[, 3], rev(d[, 5])), 
          col = cols[2], border = NA)
  lines(d[, 1], d[, 4], col = cols[3], lwd = 2)
}

# Plot gamma prior on precision as sd

plot.pre = function(a, g){
  pre = rgamma(1e6, a, g)
  v = 1 / pre
  sd = sqrt(v)
  plot(density(sd))
}
