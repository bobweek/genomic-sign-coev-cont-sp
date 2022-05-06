require(ggplot2)
require(gridExtra)
require(viridis)

mmdf = read.csv("~/gsccs-data/rast-data/minmax.csv")

nums = seq(1,4551,by=10)

for(i in nums){
  num = sprintf("%04d", i)
  datfile = paste("~/gsccs-data/rast-data/rast",num,".csv",sep="")
  rastdf = read.csv(datfile)
  nlocs = length(rastdf$x)/2
  rastdf$spp = c(rep(0,nlocs),rep(1,nlocs))
  
  ggplot(rastdf) +
    geom_raster(aes(fill=N,x=x,y=y)) +
    facet_grid(.~spp) +
    scale_fill_viridis(limits=c(mmdf$minN,mmdf$maxN)) +
    theme_minimal()
  
  ggplot(rastdf) +
    geom_raster(aes(fill=z,x=x,y=y)) +
    facet_grid(.~spp) +
    scale_fill_viridis(limits=c(mmdf$minz,mmdf$maxz)) +
    theme_minimal()
  
  ggplot(rastdf) +
    geom_raster(aes(fill=v,x=x,y=y)) +
    facet_grid(.~spp) +
    scale_fill_viridis(limits=c(mmdf$minv,mmdf$maxv)) +
    theme_minimal()
}