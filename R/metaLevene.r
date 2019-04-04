
metaLevene <- function(geno_count, z_avg, z_var, verbose=F) {
   
   if (length(dim(geno_count)) != 2) 
		stop("geno_count is not an array")
   if (length(dim(z_avg)) != 2) 
		stop("z_avg is not an array")
   if (length(dim(z_var)) != 2) 
		stop("z_var is not an array")

   if (dim(z_avg)[1] != dim(z_var)[1] | dim(geno_count)[1] != dim(z_var)[1]) 
		stop("Please make sure geno_count, z_avg and z_var all have equal number of rows")	
	
	N <- sum(geno_count, na.rm=T)
	Zs <- colSums(z_avg*geno_count,na.rm=T)/colSums(geno_count,na.rm=T)
	Z <- sum(Zs*colSums(geno_count)/N,na.rm=T)
 
	
	meta_stats_n <- (N-3)*sum(colSums(geno_count)*(Zs-Z)^2)
	
	meta_stats_d <- (3-1)*(sum(z_var*geno_count+z_avg^2*geno_count)
			-sum((colSums(z_avg*geno_count)/colSums(geno_count))^2*colSums(geno_count)))
	
	 meta_levene <- 1-pf(meta_stats_n/meta_stats_d,df1=2,df2=N-3)
 

		if (verbose == F)return("meta_levene"=meta_levene)
			else return("LeveneStats"=meta_stats_n/meta_stats_d)

}