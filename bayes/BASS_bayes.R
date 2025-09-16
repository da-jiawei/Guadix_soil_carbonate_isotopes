model{
  # Data model ----
  for(i in 1:length(ages)) {
    dp18sw.obs[i, 1] ~ dnorm(dp18sw[i], dp18sw.pre[i])
    dp18sw.pre[i] = 1 / .2 ^ 2
  }
  
  for (i in 1:length(ages)) {
    Dp17sw.obs[i, 1] ~ dnorm(Dp17sw[i], Dp17sw.pre[i])
    Dp17sw.pre[i] = 1 / 15 ^ 2
  }
  
  # steady state soil water isotope model ----
  for (i in 1:length(ages)) {
    R18p[i] = (d18p[i] / 1000 + 1) * R18smow
    dp18p[i] = log(R18p[i] / R18smow) * 1000
    dp17p[i] = Dp17p / 1000 + 0.528 * dp18p[i]
    R17p[i] = exp(dp17p[i] / 1000) * R17smow
    dT[i] = MSST - MAST
    airT[i] = MAST + dT[i] * sin(2 * pi * tsc)
    airT.K[i] = airT[i] + 273.15
    Tsoil[i] = MAST + (dT[i] * sin(2 * pi * tsc - depth[i] / d)) / exp(depth[i] / d) 
    Tsoil.K[i] = Tsoil[i] + 273.15
    
    # equilibrium fractionation factor for vapor to liquid - (Horita and Wesolowski 1994)
    alpha18.eq[i] = 1 / exp(((1.137e6 / (airT.K[i] ^ 2) - 0.4156e3/airT.K[i] - 2.0667) /1000))
    alpha17.eq[i] = alpha18.eq[i] ^ theta.eq
    
    # water vapor in air assumed to be in equilibrium with reservoir water
    R18a[i] = R18p[i] * alpha18.eq[i]
    R17a[i] = R17p[i] * alpha17.eq[i]
    
    es[i] = (0.611 * exp(17.502 * airT[i] / (airT[i] + 240.97))) * 1000 # saturated water vapor pressure from Tetens formula
    N.sat[i] = 0.01802 * es[i] / (Rgas * airT.K[i]) # saturated water vapor concentration at a given temperature
    z.bar[i] = N.sat[i] * Dsoil.eff / (evap[i] * rho) # penetration depth (m)
    z.ef1[i] = (1 - RH[i]) * z.bar[i] # the thickness of the water vapor phase region (m)
    
    # liquid water diffusivity (m2/s) (Easteal 1984)
    Dlo[i] = exp(1.6766 + 1.6817 * (1000 / airT.K[i]) - 0.5773 * (1000 / airT.K[i])^2) * 10^-9 
    Dlo.eff[i] = Dlo[i] * pore * tort # effective diffusivity of liquid water (m2/s)
    z.hat[i] = Dlo.eff[i] / evap[i] # the decay length (mean penetration depth)
    
    # Evaporation front
    h.ef[i] = RH[i] + z.ef[i] / z.bar[i] # humidity at the evaporation front
    R18.ef[i] = (alpha18.diff * R18p[i] * (z.ef[i] / z.bar[i]) + RH[i] * R18a[i]) / (h.ef[i] * alpha18.eq[i])
    R17.ef[i] = (alpha17.diff * R17p[i] * (z.ef[i] / z.bar[i]) + RH[i] * R17a[i]) / (h.ef[i] * alpha17.eq[i])
    
    z[i] = depth[i] / 100 # in meter unit
    hs[i] = min(RH[i] + z[i] / z.bar[i], 1)
    z.ef[i] = max(z.ef1[i], 1e-10)
    z.f[i] = (theta.mean / a.theta) * log(z[i] / z.ef[i]) # the modified depth function
    
    R18sw[i] = ifelse(z[i] <= z.ef[i], 
                     (alpha18.diff * R18p[i] * z[i] / z.bar[i] + RH[i] * R18a[i]) / 
                       (hs[i] * alpha18.eq[i]),
                     (R18.ef[i] - R18p[i]) * exp(-z.f[i] / z.hat[i]) + R18p[i])
    R17sw[i] = ifelse(z[i] <= z.ef[i], 
                     (alpha17.diff * R17p[i] * z[i] / z.bar[i] + RH[i] * R17a[i]) / 
                       (hs[i] * alpha17.eq[i]),
                     (R17.ef[i] - R17p[i]) * exp(-z.f[i] / z.hat[i]) + R17p[i])
    
    # d18sw = ((R18sw / R18smow) - 1) * 1000
    dp18sw[i] = log(R18sw[i] / R18smow) * 1000
    dp17sw[i] = log(R17sw[i] / R17smow) * 1000
    Dp17sw[i] = 1000 * (dp17sw[i] - 0.528 * dp18sw[i])
  }
  
  # priors of environmental parameters ----
  for (i in 1:length(ages)) {
    RH[i] ~ dunif(.1, .9)
    d18p[i] ~ dunif(-50, -5)
    # Dp17p[i] ~ dunif(0, 40)
    # MAST[i] ~ dunif(10, 20)
    # MSST[i] ~ dunif(25, 35)
    # tsc[i] ~ dunisf(0, 1)
    depth[i] ~ dunif(0, 100)
    evap[i] ~ dunif(1e-10, 1e-9)
  }
  
  # constants
  MAST = 20
  MSST = 30
  tsc = .25
  # RH = .8
  Dp17p = 10
  # evap = 2e-10
  pi = 3.1415926
  d = sqrt((2 * 0.0007) / ((2 * 3.1415 / 3.154e7) * 0.3)) # decay depth
  # evap = 5e-10
  R18smow = 0.0020052
  R17smow = 0.0003799
  R18vpdb = 0.0020672
  R17vpdb = 0.0003860
  Rgas = 8.314462 # gas constant
  theta.o = 0.05 # disconnected water content
  rho = 1000 # liquid water density (kg/m3)
  Dair = 2.44E-05 # water vapor diffusivity in air (m2/s) (Merlivat, 1978)
  a.theta = 0.05 # rate of increase of water content with depth (m-1)
  pore = 0.45 # for loess
  tort = 0.7
  theta.mean = 0.065
  theta.eq = 0.529
  # diffusion fractionation factor (Merlivat 1978)
  alpha18.diff = 1.028489
  theta.diff = 0.5185
  alpha17.diff = alpha18.diff ^ theta.diff
  Dsoil.eff = Dair * tort * (pore - theta.o) # The effective diffusivity of water vapor (m2/s)
}






