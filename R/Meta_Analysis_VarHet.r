#!/bin/bash


#nohup R CMD BATCH --no-save '--args Input_Directory="/Users/user/infiles/" OUTDIR="/Users/user/outfiles/" TRAIT = "BMI"' Meta_Analysis_VarHet.R&



args=(commandArgs(TRUE))

##args is now a list of character vectors
if(length(args)==0){
    print("No arguments supplied.")
}else{
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}


########################################################
##
## Step0: ### read in data  
##
########################################################

input_file_names <- list.files(paste(Input_Directory))
INPUT_FILES <- list()
for (i in 1:length(input_file_names)){
	INPUT_FILES[[i]] <- read.table(paste(Input_Directory, input_file_names[i], sep=""), head=T, sep="\t", na.strings=".")
}


########################################################
##
## Step1: Source all functions
##
########################################################

source("LeveneStats.r")
source("metaLevene.r")

report_meta <- function(outputs=data_output){

if (class(outputs) != "list") stop(" the given input is not a list")

mega_stats <- list()
record_ind_pval <- list()

	for (i in 1:length(outputs)){
	analyzed_var  <-  outputs[[i]]

	analyzed_var$Levene <- Levene(obs=cbind(analyzed_var$N_0, analyzed_var$N_1, analyzed_var$N_2),
						          zs=cbind(analyzed_var$ABSLT_MEAN_PHENO_0*analyzed_var$RESID_SD_PHENO_0, 
								        analyzed_var$ABSLT_MEAN_PHENO_1*analyzed_var$RESID_SD_PHENO_1, 
									    analyzed_var$ABSLT_MEAN_PHENO_2*analyzed_var$RESID_SD_PHENO_2), 
										zvars=cbind((analyzed_var$ABSLT_SD_PHENO_0*analyzed_var$RESID_SD_PHENO_0)^2,
										(analyzed_var$ABSLT_SD_PHENO_1*analyzed_var$RESID_SD_PHENO_1)^2,
										(analyzed_var$ABSLT_SD_PHENO_2*analyzed_var$RESID_SD_PHENO_2)^2), verbose=F) 

	
	subset_analyzed_var <- subset(analyzed_var, analyzed_var$CHR >0 & analyzed_var$CHR < 23 & !is.na(analyzed_var$LRT)& !is.na(analyzed_var$Levene) )

	mega_stats[[i]] <- data.frame("MarkerName"=subset_analyzed_var$MarkerName, 
				 "N_0" = subset_analyzed_var$N_0, "N_1" = subset_analyzed_var$N_1, "N_2" = subset_analyzed_var$N_2,
			     "Z_0" = subset_analyzed_var$ABSLT_MEAN_PHENO_0*subset_analyzed_var$RESID_SD_PHENO_0, 
			     "Z_1" = subset_analyzed_var$ABSLT_MEAN_PHENO_1*subset_analyzed_var$RESID_SD_PHENO_1, 
			     "Z_2" = subset_analyzed_var$ABSLT_MEAN_PHENO_2*subset_analyzed_var$RESID_SD_PHENO_2, 
			     "Z_var0" = (subset_analyzed_var$ABSLT_SD_PHENO_0*subset_analyzed_var$RESID_SD_PHENO_0)^2,
			     "Z_var1" = (subset_analyzed_var$ABSLT_SD_PHENO_1*subset_analyzed_var$RESID_SD_PHENO_1)^2,
			     "Z_var2" = (subset_analyzed_var$ABSLT_SD_PHENO_2*subset_analyzed_var$RESID_SD_PHENO_2)^2,
				 "MAF" = subset_analyzed_var$MAF)
	names(mega_stats[[i]]) <- c("MarkerName", paste(names(mega_stats[[i]])[2:11], "_", i, sep=""))		     
	
	record_ind_pval[[i]] <- data.frame("MarkerName"=subset_analyzed_var$MarkerName, 
									"Levene" = subset_analyzed_var$Levene) 
									#"Levene_Ef" = subset_analyzed_var$effsize)
	names(record_ind_pval[[i]]) <- c("MarkerName", paste(names(record_ind_pval[[i]])[2:3], "_", i, sep=""))	
	return(list=c(mega_stats, record_ind_pval)
	} 

 
 
metaLevene <- function(geno_count, z_avg, z_var, w_defined) {

	meta_stats_n <- (sum(geno_count)-3)*sum((colSums(z_avg*geno_count)/colSums(geno_count)
			-sum(z_avg*geno_count)/(sum(geno_count)))^2*colSums(geno_count))

	meta_stats_d <- (3-1)*(sum(z_var*(geno_count-1)+z_avg^2*geno_count)-sum((colSums(z_avg*geno_count)/colSums(geno_count))^2*colSums(geno_count)))
	
	 meta_levene <- 1-pf(meta_stats_n/meta_stats_d,df1=2,df2=sum(geno_count)-3)
	 
#### only study 1
	
	### dim(geno_count) == dim(w_define)
	geno_count_defined <- geno_count*w_defined
	N <- sum(geno_count_defined)

	w <- matrix(c(outer(geno_count_defined,colSums(geno_count_defined),"/")[,,1][,1],
				  outer(geno_count_defined,colSums(geno_count_defined),"/")[,,2][,2],
				  outer(geno_count_defined,colSums(geno_count_defined),"/")[,,3][,3]), 3,3, byrow=F)
	
	v <- colSums(geno_count_defined)/sum(geno_count_defined)

	Z_avg_meta <- colSums(w*z_avg)	
	Z_grand_avg_meta <- sum(v*Z_avg_meta)
	
	weight_stats_n <- N/2*(sum(Z_avg_meta^2*v) - (Z_grand_avg_meta)^2)
		
	weight_stats_d <- N/(N-3)*sum(v*(colSums(w*z_var) - colSums(z_var)/v/N + colSums(w*z_avg^2) - colSums(w*z_avg)^2))

	adj_levene <- 1-pf(weight_stats_n/weight_stats_d,df1=2,df2=N-3)
	
		return(c("meta_pval" = meta_levene, "adjusted_meta_pval" = adj_levene))

}


	combine_stats <- function(merge_data, w_defined,v_defined){
		merge_data <- as.numeric(merge_data)
		if (length(merge_data) < 9) stop("input statistics per genotype of length >= 9 are expected")	
		nb_studies <- round((length(merge_data))/9)
		geno_count <- matrix(rep(-999, 3*nb_studies), nrow=nb_studies, ncol=3, byrow=T)
		z_avg <- matrix(rep(-999, 3*nb_studies), nrow=nb_studies, ncol=3, byrow=T) 
		z_var <- matrix(rep(-999, 3*nb_studies), nrow=nb_studies, ncol=3, byrow=T)
			for (j in 1:nb_studies){
				geno_count[j,] <- as.numeric(merge_data[(1+10*(j-1)):(1+10*(j-1)+2)])
				z_avg[j,] <- as.numeric(merge_data[(4+10*(j-1)):(4+10*(j-1)+2)])
				z_var[j,] <- as.numeric(merge_data[(7+10*(j-1)):(7+10*(j-1)+2)])
								}						
		meta_pval <- metaLevene("geno_count"=geno_count, "z_avg" = z_avg, "z_var" = z_var,w_defined = w_defined)							

			#if (verbose == F)
			return(meta_pval)
			#	else return(list("GENO_COUNT"=geno_count, "Z_AVG" = z_avg, "Z_VAR" = z_var ))	 
			}					
 

########################################################
##
## Step3: Meta analyze data
##
########################################################      
  
  
    resultsCombo <- report_meta(outputs = INPUT_FILES)
    rm(INPUT_FILES)
    
    merged_STATS <- Reduce(function(x,y) {merge(x,y, by="MarkerName")}, resultsCombo$mega_stats)  
    individual_pvals <- Reduce(function(x,y) {merge(x,y, by="MarkerName")}, resultsCombo$record_ind_pval)

	#write.table(merged_STATS, paste(OUTDIR, trait, "_merged_output.txt", sep=""),  col.names=T, quote=F, row.names=F, na= "NA", sep="\t")
			
	pval_data <- apply(merged_STATS[,-1], 1, combine_stats, w_defined = matrix(c(1,1,1,1,1,1,1,1,1),3,3))[1,]
		
	pval_data_mat <- as.data.frame(pval_data)
	
	write.table(data.frame("MarkerName" = merged_STATS[,1], "Pval" = pval_data[,1]), file = paste(OUTDIR, TRAIT, "merged_Pvals_Equal_Weight.txt"), quote=F, sep="\t", col.names=T, row.names=F)
		

	