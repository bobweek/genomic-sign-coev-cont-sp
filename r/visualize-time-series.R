require(ggplot2)
require(gridExtra)

tsdf = read.csv("time-series.csv")
tsdf$X = seq(1,10*length(tsdf$X),by=10)
colnames(tsdf)[1]="t"

glb_abun_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=Nh),col="orange") +
  geom_line(aes(x=t,y=Np),col="purple") +
  xlab("Time") +
  ylab("Global Census Size") +
  theme_minimal()

Ntcorr = cor(tsdf$Nh_m[4:length(tsdf$t)],tsdf$Np_m[4:length(tsdf$t)])
lcl_abun_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=Nh_m),col="orange") +
  geom_line(aes(x=t,y=Np_m),col="purple") +
  # xlim(c(100,1000)) +
  xlab("Time") +
  ylab("Mean Local Census Size") +
  theme_minimal()

Ncorr_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=Ncorr)) +
  # xlim(c(100,1000)) +
  xlab("Time") +
  ylab("Interspecific Spatial Correlation of Abundance") +
  ylim(c(0,1)) +
  theme_minimal()

lcl_m_trt_m_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=zh_m),col="orange") +
  geom_line(aes(x=t,y=zp_m),col="purple") +
  # xlim(c(100,1000)) +
  xlab("Time") +
  ylab("Spatial Avg Local Mean Trait") +
  theme_minimal()

lcl_m_trt_v_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=zh_stdv),col="orange") +
  geom_line(aes(x=t,y=zp_stdv),col="purple") +
  # xlim(c(100,1000)) +
  xlab("Time") +
  ylab("Stdev of Local Mean Traits") +
  theme_minimal()

zcorr_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=zcorr)) +
  # xlim(c(100,1000)) +
  xlab("Time") +
  ylab("Spatial Correlation of Mean Traits") +
  ylim(c(-1,1)) +
  theme_minimal()

lcl_trt_v_m_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=vh_m),col="orange") +
  geom_line(aes(x=t,y=vp_m),col="purple") +
  # xlim(c(100,1000)) +
  xlab("Time") +
  ylab("Spatial Avg Local Trait Variance") +
  theme_minimal()

lcl_m_trt_v_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=vh_stdv),col="orange") +
  geom_line(aes(x=t,y=vp_stdv),col="purple") +
  # xlim(c(100,1000)) +
  xlab("Time") +
  ylab("Stdev of Local Trait Variance") +
  theme_minimal()

vcorr_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=vcorr)) +
  # xlim(c(100,1000)) +
  xlab("Time") +
  ylab("Spatial Correlation of Trait Variances") +
  ylim(c(-1,1)) +
  theme_minimal()

rhrp_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=rh_rp)) +
  # xlim(c(100,1000)) +
  xlab("Time") +
  ylab("Interspecific Spatial Correlation of Growth Rates") +
  ylim(c(0,1)) +
  theme_minimal()

DBzp_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=DBzh_m),col="orange") +
  geom_line(aes(x=t,y=DBzp_m),col="purple") +
  # xlim(c(100,1000)) +
  xlab("Time") +
  ylab("Phenotypic Response to Biotic Selection") +
  theme_minimal()

un_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=pr_unparasitized),col="orange") +
  geom_line(aes(x=t,y=pr_unhosted),col="purple") +
  # xlim(c(100,1000)) +
  xlab("Time") +
  ylab("Unparasitized/Unhosted") +
  ylim(c(0,1)) +
  theme_minimal()

NhDBzp_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=Nh_DBzp)) +
  # xlim(c(100,1000)) +
  xlab("Time") +
  ylab("Correlation of Host Abund & Para Biotic Response") +
  ylim(c(-1,1)) +
  theme_minimal()

#
#
#
colnames(tsdf)