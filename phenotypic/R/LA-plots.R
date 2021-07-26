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
