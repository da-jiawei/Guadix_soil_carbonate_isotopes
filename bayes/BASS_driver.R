rm(list = ls())
# load libraries
libs = c("R2jags", "tidyverse", "readxl", "ggpubr")
invisible(lapply(libs, library, character.only = TRUE))
# theme for ggplot
theme = theme(axis.ticks.length=unit(0.15, "cm"),
              axis.ticks = element_line(colour = "black"),
              text = element_text(color = "black", size = 10),
              axis.title = element_text(size = 12), 
              axis.text = element_text(color = "black", size = 10),
              legend.text = element_text(size = 10),
              legend.title = element_text(size = 10),
              panel.grid = element_blank())
# load functions
source('bayes/helpers.R')
cat("\014")

## Read and groom data ----
guadix = read_xlsx("data/Triple_Clumped_Combined.xlsx", sheet = 1)
guadix_triple = guadix[1:9, c(1:3, 14:16)]
names(guadix_triple) = c("age", "ID", "T47", "dp18sw", "dp17sw", "Dp17sw")
guadix_triple = guadix_triple |>
  mutate(across(c("age", "T47"), as.numeric))

## Parse data into series
dp18sw = guadix_triple[c("dp18sw")]
Dp17sw = guadix_triple[c("Dp17sw")]
ages = guadix_triple$age

## MCMC ----
# priors
d = list(ages = ages, 
         dp18sw.obs = dp18sw, 
         Dp17sw.obs = Dp17sw)

parms = c("d18p", "depth", "RH", "evap")

system.time({post.clp = jags.parallel(d, NULL, parms, "bayes/BASS_bayes.R",
                                      n.iter = 2e5, n.chains = 3, n.burnin = 5e4)})

View(post.clp$BUGSoutput$summary)
for (i in 1:length(parms)) {
  param = parms[i]
  plot.jpi(ages, post.clp$BUGSoutput$sims.list[[param]], n = 100, ylab = param)
}

post_data = guadix_triple |>
  select(age, T47, dp18sw, Dp17sw)
for (i in 1:length(parms)) {
  post_data$Rhat = post.clp$BUGSoutput$summary[grep(parms[i], rownames(post.clp$BUGSoutput$summary)), "Rhat"]
  post_data$n.eff = post.clp$BUGSoutput$summary[grep(parms[i], rownames(post.clp$BUGSoutput$summary)), "n.eff"]
  post_data[, paste0("post_", parms[i])] = post.clp$BUGSoutput$mean[parms[i]]
  post_data[, paste0("post_", parms[i], "_sd")] = post.clp$BUGSoutput$sd[parms[i]]
  for (p in 1:nrow(post_data)) {
    if (post_data$Rhat[p] <= 1.05 & post_data$n.eff[p] >= 200){
      post_data[p, paste0(parms[i], "_eff")] = "positive"
    } else {
      post_data[p, paste0(parms[i], "_eff")] = "negative"
    }
  }
}
write.csv(post_data, file = "bayes/BASS_bayes_vary_evap.csv")

## plot ----
post_data = read.csv("bayes/BASS_bayes_vary_evap.csv")
ggplot(post_data, aes(x = T47, y = post_d18p)) +
  geom_errorbar(aes(ymin = post_d18p - post_d18p_sd, ymax = post_d18p + post_d18p_sd),
                linewidth = .2, width = 0, color = "grey80") +
  geom_point(aes(alpha = d18p_eff),
             size = 3, shape = 21) +
  scale_alpha_manual(values = c("positive" = 1, "negative" = .3)) +
  theme_bw() + theme +
  guides(alpha = "none") +
  labs(x = expression(paste("T"[Delta][47]*" (",degree, "C)")),
       y = expression(delta^"'18"*"O"[p]*" (\u2030, VSMOW)"),
       fill = "")

ggplot(post_data, aes(x = T47, y = post_RH * 1e2)) +
  geom_errorbar(aes(ymin = (post_RH - post_RH_sd) * 1e2, ymax = (post_RH + post_RH_sd) * 1e2),
                linewidth = .2, width = 0, color = "grey80") +
  geom_point(aes(alpha = RH_eff),
             size = 3, shape = 21) +
  scale_alpha_manual(values = c("positive" = 1, "negative" = .3)) +
  theme_bw() + theme +
  guides(alpha = "none") +
  labs(x = expression(paste("T"[Delta][47]*" (",degree, "C)")),
       y = "RH (%)",
       fill = "")
ggarrange(p1, p2, nrow = 1, ncol = 2, align = "hv",
          common.legend = TRUE, legend = "right")
ggsave("figures/Bayes_BASS.jpg", width = 7, height = 3.4, dpi = 500)

ggplot(post_data, aes(x = temp, y = post_evap)) +
  geom_errorbar(aes(xmin = temp - temp.se, xmax = temp + temp.se),
                linewidth = .2, width = 0, color = "grey80") +
  geom_errorbar(aes(ymin = post_evap - post_evap_sd, ymax = post_evap + post_evap_sd),
                linewidth = .2, width = 0, color = "grey80") +
  geom_point(aes(fill = section, alpha = evap_eff),
             size = 3, shape = 21) +
  scale_fill_brewer(palette = "Paired") +
  scale_alpha_manual(values = c("positive" = 1, "negative" = .3)) +
  theme_bw() + theme +
  guides(alpha = "none") +
  labs(x = expression(paste("T"[Delta][47]*" (",degree, "C)")),
       y = expression(delta^"'18"*"O"[sw]*" (\u2030, VSMOW)"),
       fill = "")

