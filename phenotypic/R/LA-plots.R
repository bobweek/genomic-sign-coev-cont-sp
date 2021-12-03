require(ggplot2)
require(gridExtra)
require(scales)

ISL = read.csv("~/gits/genomic-sign-coev-cont-sp/phenotypic/julia/ISL.csv")
CLS = read.csv("~/gits/genomic-sign-coev-cont-sp/phenotypic/julia/CLS.csv")
LMD = read.csv("~/gits/genomic-sign-coev-cont-sp/phenotypic/julia/LMD.csv")
LMD.B = read.csv("~/gits/genomic-sign-coev-cont-sp/phenotypic/julia/LMD_B.csv")
CLS.B = read.csv("~/gits/genomic-sign-coev-cont-sp/phenotypic/julia/CLS_B.csv")
CLS.sml = subset(CLS,σₕ<=10 & σₚ<=10)
LMD.sml = subset(LMD,σₕ<=10 & σₚ<=10)
CLS.med = subset(CLS,σₕ<=100 & σₚ<=100)
LMD.med = subset(LMD,σₕ<=100 & σₚ<=100)
CLS.big = subset(CLS,σₕ>=100 & σₚ>=100)
LMD.big = subset(LMD,σₕ>=100 & σₚ>=100)

# parasite LA under island model
ISL.pl.p = ggplot(ISL,aes(x=σₕ,y=σₚ,z=LAₚ)) + 
  scale_x_log10() + scale_y_log10() + 
  ggtitle("Parasite Local Adpataion") +
  geom_contour_filled(bins=100, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
ISL.pl.h = ggplot(ISL,aes(x=σₕ,y=σₚ,z=LAₕ)) +
  scale_x_log10() + scale_y_log10() + 
  ggtitle("Host Local Adpataion") +
  geom_contour_filled(bins=100, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

grid.arrange(ISL.pl.p,ISL.pl.h,nrow=1)

# making plot of just classical index
CLS.sml.pl.p = ggplot(CLS.sml,aes(x=σₕ,y=σₚ,z=LAₚ)) + 
  scale_x_log10() + scale_y_log10() + 
  ylab("σₚ") +
  ggtitle("Parasite Local Adpataion") +
  geom_contour_filled(bins=100, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
CLS.sml.pl.h = ggplot(CLS.sml,aes(x=σₕ,y=σₚ,z=LAₕ)) + 
  scale_x_log10() + scale_y_log10() + 
  ylab("σₚ") +
  ggtitle("Host Local Adpataion") +
  geom_contour_filled(bins=100, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
grid.arrange(CLS.sml.pl.p,CLS.sml.pl.h,nrow=1)

sml_la = which(abs(CLS.sml$`LAₕ`)==min(abs(CLS.sml$`LAₕ`)))
CLS.sml$`σₕ`[sml_la]/CLS.sml$`σₚ`[sml_la]

ggplot(CLS.sml,aes(x=σₕ/σₚ)) + 
  scale_x_log10() + 
  geom_smooth(aes(y=LAₚ),color="red") +
  geom_smooth(aes(y=LAₕ),color="black") +
  annotate("text",label="Parasite",x=7.5,y=5e-6,color="red",size=3) +
  annotate("text",label="Host",x=7.5,y=-5e-6, size=3) +
  geom_vline(xintercept=1,linetype="dotted") +
  ylab("Local Adaptation") +
  xlab("Relative Dispersal Distance (Host/Parasite)") +
  theme_minimal() + 
  theme(text = element_text(size=10))

plot(CLS.sml$`σₕ`/CLS.sml$`σₚ`,CLS.sml$`LAₕ`,xlab="σₕ/σₚ",ylab="Local Adaptation",cex=0.3)
points(CLS.sml$`σₕ`/CLS.sml$`σₚ`,CLS.sml$`LAₚ`,cex=0.3,col="red")
text(9,2.5e-6,"Parasite",col="red")
text(9,-2.5e-6,"Host")
abline(v=1)
text(1.2,-1.1e-5,"1")

# small scale parasite LA
CLS.sml.pl.p = ggplot(CLS.sml,aes(x=σₕ,y=σₚ,z=LAₚ)) + 
  scale_x_log10() + scale_y_log10() + 
  ylab("Parasite Local Adpataion\nσₚ") +
  ggtitle("Classical Index") +
  geom_contour_filled(bins=100, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
LMD.sml.pl.p = ggplot(LMD.sml,aes(x=σₕ,y=σₚ,z=LAₚ)) + scale_x_log10() + scale_y_log10() + 
  ggtitle("Modified Index") +
  geom_contour_filled(bins=100, show.legend = F) +
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# small scale host LA
CLS.sml.pl.h = ggplot(CLS.sml,aes(x=σₕ,y=σₚ,z=LAₕ)) + 
  scale_x_log10() + scale_y_log10() + 
  ylab("Host Local Adpataion\nσₚ") +
  geom_contour_filled(bins=100, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal()
LMD.sml.pl.h = ggplot(LMD.sml,aes(x=σₕ,y=σₚ,z=LAₕ)) + scale_x_log10() + scale_y_log10() + 
  geom_contour_filled(bins=100, show.legend = F) +
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal()

grid.arrange(CLS.sml.pl.p, LMD.sml.pl.p, CLS.sml.pl.h, LMD.sml.pl.h, nrow=2)

# med scale parasite LA
CLS.med.pl.p = ggplot(CLS.med,aes(x=σₕ,y=σₚ,z=LAₚ)) + scale_x_log10() + scale_y_log10() + 
  ylab("Parasite Local Adpataion\nσₚ") +
  ggtitle("Classical Index") +
  geom_contour_filled(bins=100, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
LMD.med.pl.p = ggplot(LMD.med,aes(x=σₕ,y=σₚ,z=LAₚ)) + scale_x_log10() + scale_y_log10() + 
  ggtitle("Modified Index") +
  geom_contour_filled(bins=100, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

# med scale host LA
CLS.med.pl.h = ggplot(CLS.med,aes(x=σₕ,y=σₚ,z=LAₕ)) + 
  scale_x_log10() + scale_y_log10() + 
  ylab("Host Local Adpataion\nσₚ") +
  geom_contour_filled(bins=100, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal()
LMD.med.pl.h = ggplot(LMD.med,aes(x=σₕ,y=σₚ,z=LAₕ)) + scale_x_log10() + scale_y_log10() + 
  geom_contour_filled(bins=100, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal()

grid.arrange(CLS.med.pl.p, LMD.med.pl.p, CLS.med.pl.h, LMD.med.pl.h,nrow=2)

# large scale parasite LA
CLS.big.pl.p = ggplot(CLS.big,aes(x=σₕ,y=σₚ,z=LAₚ)) + scale_x_log10() + scale_y_log10() + 
  ylab("Parasite Local Adpataion\nσₚ") +
  ggtitle("Classical Index") +
  geom_contour_filled(bins=100, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
LMD.big.pl.p = ggplot(LMD.big,aes(x=σₕ,y=σₚ,z=LAₚ)) + scale_x_log10() + scale_y_log10() + 
  ggtitle("Modified Index") +
  geom_contour_filled(bins=100, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

# large scale host LA
CLS.big.pl.h = ggplot(CLS.big,aes(x=σₕ,y=σₚ,z=LAₕ)) + 
  scale_x_log10() + scale_y_log10() + 
  ylab("Host Local Adpataion\nσₚ") +
  geom_contour_filled(bins=100, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal()
LMD.big.pl.h = ggplot(LMD.big,aes(x=σₕ,y=σₚ,z=LAₕ)) + scale_x_log10() + scale_y_log10() + 
  geom_contour_filled(bins=100, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal()

grid.arrange(CLS.big.pl.p, LMD.big.pl.p, CLS.big.pl.h, LMD.big.pl.h,nrow=2)

# large scale parasite LA (BxB)
CLS.B.p = ggplot(CLS.B,aes(x=Bₕ,y=Bₚ,z=LAₚ)) + #scale_x_log10() + scale_y_log10() + 
  ylab("Parasite Local Adpataion\nBₚ") +
  ggtitle("Classical Index") +
  geom_contour_filled(bins=100, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
LMD.B.p = ggplot(LMD.B,aes(x=Bₕ,y=Bₚ,z=LAₚ)) + #scale_x_log10() + scale_y_log10() + 
  ggtitle("Modified Index") +
  geom_contour_filled(bins=100, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

# large scale host LA (BxB)
CLS.B.h = ggplot(CLS.B,aes(x=Bₕ,y=Bₚ,z=LAₕ)) + 
  #scale_x_log10() + scale_y_log10() + 
  ylab("Host Local Adpataion\nBₚ") +
  geom_contour_filled(bins=100, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal()
LMD.B.h = ggplot(LMD.B,aes(x=Bₕ,y=Bₚ,z=LAₕ)) + #scale_x_log10() + scale_y_log10() + 
  geom_contour_filled(bins=100, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal()

grid.arrange(CLS.B.p, LMD.B.p, CLS.B.h, LMD.B.h,nrow=2)

ggtitle("σₕ=σₚ=10")
