require(ggplot2)
require(viridis)
require(ggdark)

ndatpts = length(list.files(path = "~/gsccs-data/ind-data/", pattern = "csv"))

mmdf = read.csv("~/gsccs-data/minmax.csv")
mnn = min(mmdf$min)
mxx = max(mmdf$max)

species_names <- list(
  '1'="Host",
  '2'="Parasite"
)
species_labeller <- function(variable,value){
  return(species_names[value])
}

pars = read.csv("~/gsccs-data/params.csv")

for(i in 1:ndatpts){
  num = sprintf("%04d", i)
  datfile = paste("~/gsccs-data/ind-data/indData",num,".csv",sep="")
  inddf = read.csv(datfile)
  
  ttl = paste("Generation:",(i-1)*pars$tempres)

  pz <- ggplot(inddf) +
    geom_point(aes(color=z,x=x,y=y),alpha=0.5,stroke=0.75,size=0.5) +
    facet_grid(.~spp, labeller=species_labeller) +
    scale_colour_viridis(limits=c(mnn,mxx),name = "Trait\nValue") +
    ggtitle(ttl) +
    xlab("Latitude") +
    ylab("Longitude") +
    dark_theme_gray()
  
  pzf = paste("~/gsccs-data/ind-data/z",num,".png",sep="")
  ggsave(pzf,pz,width=16,height=8,bg="black")
  
}
