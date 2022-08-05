require(ggplot2)
require(gridExtra)
require(ggdark)

params = read.csv("~/gsccs-data/params.csv")
tsdf = read.csv("~/gsccs-data/time-series.csv")

colors <- c("Parasite" = "purple", "Host" = "orange")

abn = ggplot(tsdf) + 
  geom_line(aes(x=t,y=Nh,color="Host")) +
  geom_line(aes(x=t,y=Np,color="Parasite")) +
  labs(x = "Years",        
       y = "Global Census Size",        
       color = "Legend") +   
  scale_color_manual(values = colors) +
  theme_minimal(lege)

trt = ggplot(tsdf) + 
  geom_line(aes(x=t,y=zh,color="Host")) +
  geom_line(aes(x=t,y=zp,color="Parasite")) +
  labs(x = "Years",
       y = "Global Trait Means",        
       color = "Legend") +   
  scale_color_manual(values = colors, guide = "none") +
  theme_minimal()

pol = ggplot(tsdf) + 
  geom_line(aes(x=t,y=muh,color="Host")) +
  geom_line(aes(x=t,y=mup,color="Parasite")) +
  labs(x = "Time",
       y = "Num Polymorphic Loci") +   
  scale_color_manual(values = colors) +
  theme_minimal()

tcr = ggplot(tsdf) + 
  geom_line(aes(x=t,y=pntcorr)) +
  labs(x = "Time",        
       y = "Spatial Correlation of Traits",        
       color = "Legend") +   
  scale_color_manual(values = colors) +
  ylim(c(-1,1)) +
  theme_minimal()

grid.arrange(abn,trt,tcr,nrow=1)
