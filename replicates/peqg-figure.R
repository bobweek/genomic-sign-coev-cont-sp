require(ggplot2)
require(gridExtra)
library(latex2exp)

trait.diff = seq(-5,5,by=0.01)
pr.inf = function(dz){
  exp(-dz^2/2)
}

inf.pr = data.frame(pr.inf(trait.diff),trait.diff)
colnames(inf.pr) = c("pr","dz")

inf.pr$`Infection Probability`

xlb = TeX("Host Trait - Parasite Trait ($\\Delta z$)")
ylb = TeX("Infection Probability ($\\alpha$)")

phint = ggplot(inf.pr) + 
  geom_line(aes(x=dz,y=pr),color="black") +
  annotate("text", x=-3, y=0.8, 
           label= TeX("$\\alpha(\\Delta z)=\\exp\\left(-\\frac{\\gamma}{2}\\Delta z^2\\right)$"),
           color="black") + 
  xlab(xlb) +
  ylab(ylb) +
  # ggtitle("Phenotypic Interface = Trait Matching") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
