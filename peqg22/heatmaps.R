library(ggplot2)

ss1=read.csv("~/gsccs-data/replicates/sxs/ss1.csv", header = F)
ss2=read.csv("~/gsccs-data/replicates/sxs/ss2.csv", header = F)

s =read.csv("~/gsccs-data/replicates/Lxs/s.csv",  header = F)
L1=read.csv("~/gsccs-data/replicates/Lxs/L1.csv", header = F)
L2=read.csv("~/gsccs-data/replicates/Lxs/L2.csv", header = F)
L3=read.csv("~/gsccs-data/replicates/Lxs/L3.csv", header = F)
L4=read.csv("~/gsccs-data/replicates/Lxs/L4.csv", header = F)
L5=read.csv("~/gsccs-data/replicates/Lxs/L5.csv", header = F)

# bin the Ls

cslfrac.s = read.csv("~/gsccs-data/replicates/sxs/cslfrac.csv", header = F)
cslfrac.L = read.csv("~/gsccs-data/replicates/Lxs/cslfrac.csv", header = F)

cslfrac.s$V1 = rnorm(length(cslfrac.s$V1))
cslfrac.L$V1 = rnorm(length(cslfrac.L$V1))

sxs.df = data.frame(x=ss1$V1,y=ss2$V1,value=cslfrac.s$V1)

ggplot(sxs.df, aes(x = x, y = y, fill = value)) +
  scale_x_log10() +
  scale_y_log10() +
  geom_tile() +
  xlab("sₕ") +
  ylab("sₚ") +
  theme_minimal()

Lxs.df = data.frame(x=L1$V1,y=s$V1,value=cslfrac.L$V1)

ggplot(Lxs.df, aes(x = x, y = y, fill = value)) +
  scale_y_log10() +
  geom_tile() +
  xlab("Lₚ") +
  ylab("sₚ") +
  theme_minimal()
