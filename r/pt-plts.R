# ffmpeg -r 30 -i N%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p N.mp4
# ffmpeg -r 30 -i z%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p z.mp4
# ffmpeg -r 30 -i v%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p v.mp4

# check if kennel needs vax's

require(ggplot2)
require(gridExtra)
require(viridis)
require(ggdark)

mmdf = read.csv("~/gsccs-data/ind-data/minmax.csv")

for(i in 0:599){
  num = sprintf("%04d", i)
  datfile = paste("~/gsccs-data/ind-data/indData",num,".csv",sep="")
  inddf = read.csv(datfile)
  
  pz <- ggplot(inddf) +
    geom_point(aes(color=z,x=x,y=y),alpha=0.5,stroke=0.75,size=0.5) +
    facet_grid(.~spp) +
    scale_colour_viridis() +
    dark_theme_gray()
  pzf = paste("~/gsccs-data/ind-data/z",num,".png",sep="")
  ggsave(pzf,pz,width=16,height=8,bg="black")
  
}
