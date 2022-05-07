# ffmpeg -r 30 -i N%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p N.mp4
# ffmpeg -r 30 -i z%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p z.mp4
# ffmpeg -r 30 -i v%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p v.mp4

require(ggplot2)
require(gridExtra)
require(viridis)

mmdf = read.csv("~/gsccs-data/rast-data/minmax.csv")

nums = seq(1,4551,by=10)

for(i in 1:length(nums)){
  num = sprintf("%04d", nums[i])
  datfile = paste("~/gsccs-data/rast-data/rast",num,".csv",sep="")
  rastdf = read.csv(datfile)
  nlocs = length(rastdf$x)/2
  rastdf$spp = c(rep(0,nlocs),rep(1,nlocs))
  
  pnum = sprintf("%04d",i)
  
  pN = ggplot(rastdf) +
    geom_raster(aes(fill=N,x=x,y=y)) +
    facet_grid(.~spp) +
    scale_fill_viridis(limits=c(mmdf$minN,mmdf$maxN)) +
    theme_minimal()
  
  pNf = paste("~/gsccs-data/rast-data/N",pnum,".png",sep="")
  ggsave(pNf,pN,width=8,height=4,bg="white")
  
  pz = ggplot(rastdf) +
    geom_raster(aes(fill=z,x=x,y=y)) +
    facet_grid(.~spp) +
    scale_fill_viridis(limits=c(mmdf$minz,mmdf$maxz)) +
    theme_minimal()
  
  pzf = paste("~/gsccs-data/rast-data/z",pnum,".png",sep="")
  ggsave(pzf,pz,width=8,height=4,bg="white")
  
  pv = ggplot(rastdf) +
    geom_raster(aes(fill=v,x=x,y=y)) +
    facet_grid(.~spp) +
    scale_fill_viridis(limits=c(mmdf$minv,mmdf$maxv)) +
    theme_minimal()

  pvf = paste("~/gsccs-data/rast-data/v",pnum,".png",sep="")
  ggsave(pvf,pv,width=8,height=4,bg="white")
  
}