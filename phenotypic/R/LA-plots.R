require(ggplot2)
require(gridExtra)
require(scales)

ISL = read.csv("~/gits/genomic-sign-coev-cont-sp/phenotypic/julia/ISL.csv")
CLS = read.csv("~/gits/genomic-sign-coev-cont-sp/phenotypic/julia/CLS.csv")
LMD = read.csv("~/gits/genomic-sign-coev-cont-sp/phenotypic/julia/LMD.csv")
CLS.sml = subset(CLS,σₕ<=10 & σₚ<=10)
LMD.sml = subset(LMD,σₕ<=10 & σₚ<=10)
CLS.big = subset(CLS,σₕ>=10 & σₚ>=10)
LMD.big = subset(LMD,σₕ>=10 & σₚ>=10)

# small scale parasite LA
ISL.pl.p = ggplot(ISL,aes(x=σₕ,y=σₚ,z=LAₚ)) + 
  scale_x_log10() + scale_y_log10() + 
  ggtitle("Parasite Local Adpataion") +
  geom_contour_filled(bins=1000, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
ISL.pl.h = ggplot(ISL,aes(x=σₕ,y=σₚ,z=LAₕ)) + 
  scale_x_log10() + scale_y_log10() + 
  ggtitle("Host Local Adpataion") +
  geom_contour_filled(bins=1000, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

grid.arrange(ISL.pl.p,ISL.pl.h,nrow=1)

# small scale parasite LA
CLS.sml.pl.p = ggplot(CLS.sml,aes(x=σₕ,y=σₚ,z=LAₚ)) + 
  scale_x_log10() + scale_y_log10() + 
  ylab("Parasite Local Adpataion\nσₚ") +
  ggtitle("Classical Index") +
  geom_contour_filled(bins=1000, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
LMD.sml.pl.p = ggplot(LMD.sml,aes(x=σₕ,y=σₚ,z=LAₚ)) + scale_x_log10() + scale_y_log10() + 
  ggtitle("Modified Index") +
  geom_contour_filled(bins=1000, show.legend = F) +
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# small scale host LA
CLS.sml.pl.h = ggplot(CLS.sml,aes(x=σₕ,y=σₚ,z=LAₕ)) + 
  scale_x_log10() + scale_y_log10() + 
  ylab("Host Local Adpataion\nσₚ") +
  geom_contour_filled(bins=1000, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal()
LMD.sml.pl.h = ggplot(LMD.sml,aes(x=σₕ,y=σₚ,z=LAₕ)) + scale_x_log10() + scale_y_log10() + 
  geom_contour_filled(bins=1000, show.legend = F) +
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal()

grid.arrange(CLS.sml.pl.p, LMD.sml.pl.p, CLS.sml.pl.h, LMD.sml.pl.h, nrow=2)

# large scale parasite LA
CLS.big.pl.p = ggplot(CLS.big,aes(x=σₕ,y=σₚ,z=LAₚ)) + scale_x_log10() + scale_y_log10() + 
  ylab("Parasite Local Adpataion\nσₚ") +
  ggtitle("Classical Index") +
  geom_contour_filled(bins=1000, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
LMD.big.pl.p = ggplot(LMD.big,aes(x=σₕ,y=σₚ,z=LAₚ)) + scale_x_log10() + scale_y_log10() + 
  ggtitle("Modified Index") +
  geom_contour_filled(bins=1000, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

# large scale host LA
CLS.big.pl.h = ggplot(CLS.big,aes(x=σₕ,y=σₚ,z=LAₕ)) + 
  scale_x_log10() + scale_y_log10() + 
  ylab("Host Local Adpataion\nσₚ") +
  geom_contour_filled(bins=1000, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal()
LMD.big.pl.h = ggplot(LMD.big,aes(x=σₕ,y=σₚ,z=LAₕ)) + scale_x_log10() + scale_y_log10() + 
  geom_contour_filled(bins=1000, show.legend = F) + 
  geom_contour(aes(colour = after_stat(level))) + 
  scale_colour_viridis_c(name="LA",labels=scientific_format()) + 
  theme_minimal()

grid.arrange(CLS.big.pl.p, LMD.big.pl.p, CLS.big.pl.h, LMD.big.pl.h,nrow=2)

