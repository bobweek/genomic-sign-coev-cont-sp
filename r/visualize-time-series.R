require(ggplot2)
require(gridExtra)

tsdf = read.csv("~/gsccs-data/time-series.csv")
tsdf$X = seq(1,10*length(tsdf$X),by=10)
colnames(tsdf)[1]="t"

colors <- c("Parasite" = "purple", "Host" = "orange")

glb_abun_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=Nh,color="Host")) +
  geom_line(aes(x=t,y=Np,color="Parasite")) +
  labs(x = "Time",        
       y = "Global Census Size",        
       color = "Legend") +   
  scale_color_manual(values = colors) +
  theme_minimal()

Ntcorr = cor(tsdf$Nh_m[4:length(tsdf$t)],tsdf$Np_m[4:length(tsdf$t)])
lcl_abun_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=Nh_m,color="Host")) +
  geom_line(aes(x=t,y=Np_m,color="Parasite")) +
  # xlim(c(100,1000)) +
  labs(x = "Time",        
       y = "Mean Local Census Size",        
       color = "Legend") +   
  scale_color_manual(values = colors) +
  theme_minimal()

Ncorr_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=Ncorr)) +
  # xlim(c(100,1000)) +
  labs(x = "Time",        y = "Spatial Avg Local Mean Trait",        color = "Legend") +   scale_color_manual(values = colors) +
  ylab("Interspecific Spatial Correlation of Abundance") +
  ylim(c(0,1)) +
  theme_minimal()

lcl_m_trt_m_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=zh_m,color="Host")) +
  geom_line(aes(x=t,y=zp_m,color="Parasite")) +
  labs(x = "Time",
       y = "Spatial Avg Local Mean Trait",
       color = "Legend") +
  scale_color_manual(values = colors) +
  theme_minimal()

lcl_m_trt_v_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=zh_stdv,color="Host")) +
  geom_line(aes(x=t,y=zp_stdv,color="Parasite")) +
  labs(x = "Time",        
       y = "Stdev of Local Mean Traits",        
       color = "Legend") +   
  scale_color_manual(values = colors) +
  theme_minimal()

zcorr_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=zcorr)) +
  labs(x = "Time",        
       y = "Spatial Correlation of Mean Traits",        
       color = "Legend") +   
  scale_color_manual(values = colors) +
  ylim(c(-1,1)) +
  theme_minimal()

lcl_trt_v_m_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=sqrt(vh_m),color="Host")) +
  geom_line(aes(x=t,y=sqrt(vp_m),color="Parasite")) +
  labs(x = "Time",        
       y = "Spatial Avg Local Trait Stddev",        
       color = "Legend") +   
  scale_color_manual(values = colors) +
  theme_minimal()

lcl_m_trt_v_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=vh_stdv,color="Host")) +
  geom_line(aes(x=t,y=vp_stdv,color="Parasite")) +
  labs(x = "Time",        
       y = "Stdev of Local Trait Variance",        
       color = "Legend") +   
  scale_color_manual(values = colors) +
  theme_minimal()

vcorr_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=vcorr)) +
  labs(x = "Time",        
       y = "Spatial Correlation of Trait Variances",        
       color = "Legend") +   
  scale_color_manual(values = colors) +
  ylim(c(-1,1)) +
  theme_minimal()

rhrp_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=rh_rp)) +
  labs(x = "Time",        
       y = "Interspecific Spatial Correlation of Growth Rates",        
       color = "Legend") +   
  scale_color_manual(values = colors) +
  ylim(c(0,1)) +
  theme_minimal()

DBzp_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=DBzh_m,color="Host")) +
  geom_line(aes(x=t,y=DBzp_m,color="Parasite")) +
  labs(x = "Time",        
       y = "Phenotypic Response to Biotic Selection",        
       color = "Legend") +   
  scale_color_manual(values = colors) +
  theme_minimal()

un_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=pr_unparasitized,color="Host")) +
  geom_line(aes(x=t,y=pr_unhosted,color="Parasite")) +
  labs(x = "Time",        
       y = "Unparasitized/Unhosted",        
       color = "Legend") +   
  scale_color_manual(values = colors) +
  ylim(c(0,1)) +
  theme_minimal()

NhDBzp_p = ggplot(tsdf) + 
  geom_line(aes(x=t,y=Nh_DBzp)) +
  labs(x = "Time",        
       y = "Correlation of Host Abund & Para Biotic Response",        
       color = "Legend") +   
  scale_color_manual(values = colors) +
  ylim(c(-1,1)) +
  theme_minimal()

#
#
#
colnames(tsdf)