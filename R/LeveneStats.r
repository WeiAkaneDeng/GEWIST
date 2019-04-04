  Levene <- function(obs, zs, zvars, verbose=F){
 N <- rowSums(obs, na.rm=T)
 Z <- rowSums(zs*obs/N,na.rm=T)
 
 df1 <- 2 - rowSums(zs==0)
 df2 <- rowSums(obs)-(3-rowSums(zs==0))
 
 pval <- 1-pf(N*rowSums((zs-Z)^2*obs, na.rm=T)/(3-1)/rowSums(zvars*(obs-1),na.rm=T),df1=df1,df2=df2)  
 pval[rowSums(zvars > 0, na.rm=TRUE) <= rep(2,length(N))] <- NA
if (verbose == F)return("Levene_pvalue"=pval)
			else return("LeveneStats"=N*rowSums((zs-Z)^2*obs, na.rm=T)/(3-1)/rowSums(zvars*(obs-1),na.rm=T))
 }