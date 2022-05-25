require(ggplot2)
require(ggdark)

iscaf_cor = read.csv("~/gsccs-data/iscaf-cor.csv")
iscaf_cov = read.csv("~/gsccs-data/iscaf-cov.csv")
csl_iscaf_cor = read.csv("~/gsccs-data/csl-iscaf-cor.csv")
csl_iscaf_cov = read.csv("~/gsccs-data/csl-iscaf-cov.csv")

iscaf_cor_flat = data.frame(unname(as.vector(data.matrix(iscaf_cor))))
iscaf_cov_flat = data.frame(unname(as.vector(data.matrix(iscaf_cov))))
colnames(iscaf_cor_flat) <- c("iscaf")
colnames(iscaf_cov_flat) <- c("iscaf")

ncsl = dim(csl_iscaf_cov)[1]*dim(csl_iscaf_cov)[2]
csl_iscaf_cor_flat = data.frame(unname(as.vector(data.matrix(csl_iscaf_cor))),rep(0,ncsl))
csl_iscaf_cov_flat = data.frame(unname(as.vector(data.matrix(csl_iscaf_cov))),rep(0,ncsl))
colnames(csl_iscaf_cor_flat) <- c("iscaf","y")
colnames(csl_iscaf_cov_flat) <- c("iscaf","y")

lim_cor = max(abs(csl_iscaf_cor_flat$iscaf))
lim_cov = max(abs(csl_iscaf_cov_flat$iscaf))

ild_cor = ggplot(iscaf_cor_flat, aes(x=iscaf)) + 
  geom_density() +
  geom_point(data=csl_iscaf_cor_flat, aes(x=iscaf,y=y),color="red") +
  xlab("Interspecific Spatial Correlations") +
  ylab("Density") +
  xlim(c(-lim_cor,lim_cor)) +
  theme_minimal()

ild_cov = ggplot(iscaf_cov_flat, aes(x=iscaf)) + 
  geom_density() +
  geom_point(data=csl_iscaf_cov_flat, aes(x=iscaf,y=y),color="red") +
  xlab("Interspecific Spatial Covariances") +
  ylab("Density") +
  xlim(c(-lim_cov,lim_cov)) +
  theme_minimal()

ggsave("~/gsccs-data/")

# h_ntl_snps = read.csv("~/gsccs-data/h_ntl_snps.csv")
# h_csl_snps = read.csv("~/gsccs-data/h_csl_snps.csv")
# colnames(h_ntl_snps) <- c("Host Neutral Snps")
# colnames(h_csl_snps) <- c("Host Causal Snps")
# 
# p_ntl_snps = read.csv("~/gsccs-data/p_ntl_snps.csv")
# p_csl_snps = read.csv("~/gsccs-data/p_csl_snps.csv")
# colnames(p_ntl_snps) <- c("Parasite Neutral Snps")
# colnames(p_csl_snps) <- c("Parasite Causal Snps")
# 
# xcoord = seq(from=0,to=4,length.out=1e8)
# hycoord = rep(1,1e8)
# pycoord = rep(0,1e8)
# 
# pts = data.frame(rep(xcoord,2),c(hycoord,pycoord))
