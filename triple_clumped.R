rm(list = ls())
pacman::p_load(tidyverse, readxl, patchwork)
theme = theme(
  panel.grid = element_blank(),
  axis.ticks.length = unit(2, units = "mm"),
  axis.title = element_text(size = 12),
  axis.text = element_text(size = 10, color = "black"),
  legend.title = element_text(size = 10),
  legend.background = element_blank(),
  legend.position = c(.2 , .3)
)


aron = read_xlsx("data/Aron2023.xlsx")
kelson = read_xlsx("data/Kelson2023.xlsx", sheet = "data")
guadix = read_xlsx("data/Triple_Clumped_Combined.xlsx", sheet = 1)
guadix_triple = guadix[1:9, c(1:3, 14:16)]
names(guadix_triple) = c("age", "ID", "T47", "dp18sw", "dp17sw", "Dp17sw")
guadix_triple = guadix_triple |>
  mutate(across(c("age", "T47"), as.numeric))
guadix_clumped = guadix[13:28, c(1:3, 8)]
names(guadix_clumped) = guadix_clumped[1,]
guadix_clumped = guadix_clumped[-1,]
guadix_clumped = guadix_clumped |>
  mutate(across(c("Age", "T", "d18Osw"), as.numeric))
bayes = read_csv("bayes/BASS_bayes_vary_evap.csv")
bayes = bayes |>
  mutate(d18p_low = post_d18p - post_d18p_sd,
         d18p_high = post_d18p + post_d18p_sd)

# dp18-DP17 ----
p1 = ggplot(guadix_triple) +
  geom_point(aes(x = dp18sw, y = Dp17sw, fill = T47),
             shape = 21, size = 3) +
  scale_fill_distiller(palette = "RdBu") +
  theme_bw() + theme +
  labs(x = expression(delta^"'18"*"O"[sw]*" (\u2030)"),
       y = expression(Delta^"'17"*"O"[sw]*" (per meg)"),
       fill = expression(paste("T"[Delta][47]*" (", degree, "C)")))

p2 = ggplot(guadix_triple) +
  geom_point(data = aron, aes(x = dp18p, y = Dp17p),
             shape = 22, size = 2, color = "gray90") +
  geom_point(aes(x = dp18sw, y = Dp17sw, fill = T47),
             shape = 21, size = 3) +
  scale_fill_distiller(palette = "RdBu") +
  ggtitle("vs. meteoric water") +
  theme_bw() + theme +
  labs(x = expression(delta^"'18"*"O"[sw]*" (\u2030)"),
       y = expression(Delta^"'17"*"O"[sw]*" (per meg)"),
       fill = expression(paste("T"[Delta][47]*" (", degree, "C)")))

p3 = ggplot(guadix_triple) +
  geom_point(data = kelson, aes(x = dp18sw, y = Dp17sw),
             shape = 23, size = 2, color = "gray90") +
  geom_point(aes(x = dp18sw, y = Dp17sw, fill = age),
             shape = 21, size = 3) +
  scale_fill_viridis_c(option = "mako") +
  ggtitle("vs. Holocene soil water") +
  theme_bw() + theme +
  labs(x = expression(delta^"'18"*"O"[sw]*" (\u2030)"),
       y = expression(Delta^"'17"*"O"[sw]*" (per meg)"),
       fill = "Age (Ma)") 

p4 = ggplot(bayes, aes(x = T47, y = post_d18p, fill = age)) +
  geom_errorbar(aes(ymin = d18p_low,
                    ymax = d18p_high),
                width = 0, linewidth = .1) +
  geom_point(shape = 21, size = 3) +
  scale_fill_viridis_c(option = "mako") +
  theme_bw() + theme +
  labs(x = expression(paste("T"[Delta][47]*" (", degree, "C)")),
       y = expression(delta^"18"*"O"[p]*" (\u2030)"),
       fill = "Age (Ma)") 

p1 + p4 + p2 + p3 + 
  plot_layout(ncol = 2)
ggsave("figure/guadix_dp18_Dp17.png", width = 6.7, height = 7.4, dpi = 300)

# time series ----
pal = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99")

png("figure/timeseries.png", 4, 5.3, units = "in", res = 300)
par(mar = c(4, 4, 1, 4))
plot(-1, 0, xlim = c(2.5, 4), ylim = c(0, 4), axes = FALSE,
     xlab = "", ylab = "")

yext = range(bayes$d18p_low, bayes$d18p_high)
tix = seq(-13, -5, 2)
fit = loess(post_d18p ~ age, data = bayes, span = 0.6)
x.new = seq(min(bayes$age), max(bayes$age), length.out = 200)
y.pred = predict(fit, newdata = data.frame(age = x.new))
loess.rs = cbind(x.new,
                 3 + (y.pred - min(tix)) / diff(range(tix)))
lines(loess.rs[,1], loess.rs[,2], col = pal[3], lwd = 2)
d18p.rs = cbind(bayes$age,
                3 + (bayes$post_d18p - min(tix)) / diff(range(tix)),
                3 + (bayes$d18p_low - min(tix)) / diff(range(tix)),
                3 + (bayes$d18p_high - min(tix)) / diff(range(tix)))
arrows(d18p.rs[, 1], d18p.rs[, 3], 
       d18p.rs[, 1], d18p.rs[, 4], 
       col = "black", angle=90, length=0, code = 0)
points(d18p.rs[,1], d18p.rs[,2], pch = 21, col = pal[4], bg = "white")
axis(4, 3 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(delta^"18"*"O"[p]*" (\u2030)"), 4, line = 2.5, at = 3.5)


yext = range(guadix_triple$Dp17sw)
tix = seq(100, -50, -20)
fit = loess(Dp17sw ~ age, data = guadix_triple, span = 0.6)
x.new = seq(min(guadix_triple$age), max(guadix_triple$age), length.out = 200)
y.pred = predict(fit, newdata = data.frame(age = x.new))
loess.rs = cbind(x.new,
                 3 - (y.pred - min(tix)) / diff(range(tix)))
lines(loess.rs[,1], loess.rs[,2], col = pal[1], lwd = 2)
Dp17sw.rs = cbind(guadix_triple$age,
                  3 - (guadix_triple$Dp17sw - min(tix)) / diff(range(tix)))
points(Dp17sw.rs[,1], Dp17sw.rs[,2], pch = 21, col = pal[2], bg = "white")
axis(2, 3 - (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(Delta^"'17"*"O"[sw]*" (per meg)"), 2, line = 2.5, at = 2.5)

yext = range(guadix_clumped$T)
tix = seq(-8, 24, 4)
fit = loess(T ~ Age, data = guadix_clumped, span = 0.5)
x.new = seq(min(guadix_clumped$Age), max(guadix_clumped$Age), length.out = 200)
y.pred = predict(fit, newdata = data.frame(Age = x.new))
loess.rs = cbind(x.new,
                 1 + (y.pred - min(tix)) / diff(range(tix)))
lines(loess.rs[,1], loess.rs[,2], col = pal[3], lwd = 2)
clumped.rs = cbind(guadix_clumped$Age,
                   1 + (guadix_clumped$T - min(tix)) / diff(range(tix)))
points(clumped.rs[,1], clumped.rs[,2], pch = 21, col = pal[4], bg = "white")
axis(4, 1 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(paste("T"[Delta][47]*" (", degree, "C)")), 4, line = 2.5, at = 1.5)

yext = range(guadix_clumped$d18Osw)
tix = seq(-11, -5, 2)
fit = loess(d18Osw ~ Age, data = guadix_clumped, span = 0.5)
x.new = seq(min(guadix_clumped$Age), max(guadix_clumped$Age), length.out = 200)
y.pred = predict(fit, newdata = data.frame(Age = x.new))
loess.rs = cbind(x.new,
                 0 + (y.pred - min(tix)) / diff(range(tix)))
lines(loess.rs[,1], loess.rs[,2], col = pal[1], lwd = 2)
d18sw.rs = cbind(guadix_clumped$Age,
                   0 + (guadix_clumped$d18Osw - min(tix)) / diff(range(tix)))
points(d18sw.rs[,1], d18sw.rs[,2], pch = 21, col = pal[2], bg = "white")
axis(2, 0 + (tix - min(tix)) / diff(range(tix)), tix)
mtext(expression(delta^"'18"*"O"[sw]*" (\u2030)"), 2, line = 2.5, at = 0.5)

text(2.6, 3.5, "a", cex = 1, col = "black", font = 2)
text(2.6, 2.5, "b", cex = 1, col = "black", font = 2)
text(2.6, 1.5, "c", cex = 1, col = "black", font = 2)
text(2.6, 0.5, "d", cex = 1, col = "black", font = 2)
axis(1, at = seq(2.5, 4, .1))
mtext("Age (Ma)", 1, line = 2)

dev.off()
