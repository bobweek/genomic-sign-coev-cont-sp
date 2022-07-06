require(ggplot2)
require(ggdark)

nsp_ild = read.csv("~/gsccs-data/spacecomp/nonspatial/ild-flat.csv", header = F)
sp_ild = read.csv("~/gsccs-data/spacecomp/spatial/ild-flat.csv", header = F)

nsp_ild = nsp_ild$V1[nsp_ild$V1 != 0]
sp_ild  = sp_ild$V1[sp_ild$V1 != 0]

ild = data.frame(c(sp_ild,nsp_ild),c(rep("spatial",length(sp_ild)),rep("nonspatial",length(nsp_ild))))
colnames(ild) = c("ild","model")

ild$ild <- factor(ild$model, as.character(ild$model))

ggplot(ild) + 
  # geom_histogram(aes(x=ild,y=..density..,fill=model),
  #                bins=100,alpha=0.5,position = "identity") +
  geom_density(aes(x=ild,y=log(..density..),fill=factor(model,levels = c("spatial","nonspatial"))),alpha=0.25,color="black") +
  xlim(c(0,2e-3)) +
  # ylim(c(0,25e3)) +
  xlab("Interspecific Linkage Disequilibrium") +
  ylab("log(Density)") +
  theme_minimal()+ guides(fill=guide_legend(title="Model"))

nsp_ild = read.csv("~/gsccs-data/spacecomp/nonspatial/ild-flat-cor.csv", header = F)
sp_ild = read.csv("~/gsccs-data/spacecomp/spatial/ild-flat-cor.csv", header = F)

nsp_ild = abs(nsp_ild$V1[nsp_ild$V1 != 0])
sp_ild  = abs(sp_ild$V1[sp_ild$V1 != 0])

ild = data.frame(c(sp_ild,nsp_ild),c(rep("spatial",length(sp_ild)),rep("nonspatial",length(nsp_ild))))
colnames(ild) = c("ild","model")

ild$ild <- factor(ild$model, as.character(ild$model))

ggplot(ild) + 
  # geom_histogram(aes(x=ild,y=log(..density..),fill=model),
                 # bins=100,alpha=0.5,position = "identity") +
  geom_density(aes(x=ild,y=log(..density..),fill=factor(model,levels = c("spatial","nonspatial"))),alpha=0.25,color="black") +
  xlim(c(0,0.03)) +
  # ylim(c(0,25e3)) +
  xlab("Interspecific Linkage Disequilibrium") +
  ylab("Density") +
  theme_minimal()+ guides(fill=guide_legend(title="Model"))


csl_iscaf_cor = read.csv("~/gsccs-data/csl-iscaf-cor-flat.csv", header = F)
csl_iscaf_cov = read.csv("~/gsccs-data/csl-iscaf-cov-flat.csv", header = F)

kprs = list(allcor=1,allcov=2,cslcor=3,cslcov=4)
dfs = list(allcor=iscaf_cor,allcov=iscaf_cov,cslcor=csl_iscaf_cor,cslcov=csl_iscaf_cov)

for(i in 1:2){
  dfs[[i]]$V1 = abs(dfs[[i]]$V1)
  qtl = quantile(dfs[[i]]$V1,0.99)
  kprs[[i]] = which(dfs[[i]]$V1>qtl)
}
for(i in 3:4){
  dfs[[i]]$V1 = abs(dfs[[i]]$V1)
  qtl = quantile(dfs[[i-2]]$V1,0.99)
  kprs[[i]] = which(dfs[[i]]$V1>qtl)
}

lumpnms = c(rep("all",length(kprs$allcor)),rep("csl",length(kprs$cslcor)))
lumpcrs = c(dfs$allcor$V1[kprs$allcor],dfs$cslcor$V1[kprs$cslcor])
cor_df = data.frame(lumpcrs,lumpnms)
colnames(cor_df) = c("iscaf","set")

lumpnms = c(rep("all",length(kprs$allcov)),rep("csl",length(kprs$cslcov)))
lumpcvs = c(dfs$allcov$V1[kprs$allcov],dfs$cslcov$V1[kprs$cslcov])
cov_df = data.frame(lumpcvs,lumpnms)
colnames(cov_df) = c("iscaf","set")

ggplot(cor_df) + 
  geom_histogram(aes(x=iscaf,y=..density..,fill=set),
                 bins=50,alpha=0.25,position = "identity") +
  xlab("99th Perc Spatial Correlations") +
  ylab("Density") +
  dark_theme_gray()

ggplot(cov_df) + 
  geom_histogram(aes(x=iscaf,y=..density..,fill=set),
                 bins=50,alpha=0.25,position = "identity") +
  xlab("99th Perc Spatial Covariances") +
  ylab("Density") +
  dark_theme_gray()
