library(ggplot2)

ss1=read.csv("~/gsccs-data/replicates/sxs/ss1.csv", header = F)
ss2=read.csv("~/gsccs-data/replicates/sxs/ss2.csv", header = F)

s=read.csv("~/gsccs-data/replicates/Lxs/s.csv", header = F)
L1=read.csv("~/gsccs-data/replicates/Lxs/L1.csv", header = F)
L2=read.csv("~/gsccs-data/replicates/Lxs/L2.csv", header = F)
L3=read.csv("~/gsccs-data/replicates/Lxs/L3.csv", header = F)
L4=read.csv("~/gsccs-data/replicates/Lxs/L4.csv", header = F)
L5=read.csv("~/gsccs-data/replicates/Lxs/L5.csv", header = F)

# bin the Ls

cslfrac.s = read.csv("~/gsccs-data/replicates/sxs/cslfrac.csv", header = F)
cslfrac.L = read.csv("~/gsccs-data/replicates/Lxs/cslfrac.csv", header = F)

df = data.frame(x=ss1$V1,y=ss2$V1,value=s$V1)

ggplot(df, aes(x = x, y = y, fill = value)) +
  geom_tile() 
