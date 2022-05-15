require(ggplot2)
require(viridis)
require(ggdark)

mmdf = read.csv("~/gsccs-data/minmax.csv")

ndatpts = length(list.files(path = "~/gsccs-data/ind-data/", pattern = "csv"))

for(i in 1:ndatpts){
  num = sprintf("%04d", i)
  datfile = paste("~/gsccs-data/ind-data/indData",num,".csv",sep="")
  inddf = read.csv(datfile)
  
  pz <- ggplot(inddf) +
    geom_point(aes(color=z,x=x,y=y),alpha=0.5,stroke=0.75,size=0.5) +
    facet_grid(.~spp) +
    scale_colour_viridis() +
    # scale_colour_viridis(limits=c(mmdf$min,mmdf$max)) +
    dark_theme_gray()
  pzf = paste("~/gsccs-data/ind-data/z",num,".png",sep="")
  ggsave(pzf,pz,width=16,height=8,bg="black")
  
}
