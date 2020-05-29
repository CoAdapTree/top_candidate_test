##Identifying top candidate genes for each environment
##Pooja Singh
##pooja.singh09@gmail.com
##April2019
##CoAdapTree
##works for as many environments as your heart desires
##binom expectation error has been fixed
##removing genes with no outliers from the pout calculation
##per environment
##including dbinom in calculation for plotting in mahanttan plots


library(dplyr)
library(tidyr)


## Read in the gene annotation dataset
gene_anno <- read.table("JP_Pita2pangenome.filtered.genes.final.gff", sep = "\t")

# Filter out extra information
GA <- data.frame("CHROM" = gene_anno$V1,
                 "Start_Pos" = gene_anno$V4,
                 "End_Pos" = gene_anno$V5,
                 "Gene" = gene_anno$V9)

GA <- arrange(GA, CHROM, as.numeric(Start_Pos))



##read in correlations data
rho1 <- read.table("pooja", header = TRUE)
rho <- separate(rho1, snp,  c("CHROM", "POS"), sep="-")



# get unique envs
envs <- unique(rho$env)



## top candidate test
top_candidate_test<-function(data, gene_list, binomial_cutoff, threshold){
  res_name<-c("Environment", "Chromosome","Gene_Start", "Gene_End", "SNPs_in_Gene","Observed_Outliers", "Gene_Name", "Pout")
  
  
  #Loop by environment
  Res<-list()
  Count <- 0

  for (i in envs) {
    Count <- Count + 1
    sub1 <- data[data[, "env"] == unique(data$env)[Count],]
    genes <- gene_list
    rho_threshold<-quantile(sub1$pvalue, threshold) # Calculate the pvalue at the quantile specified by the threshold parameter
    
    results <- matrix(NA, nrow=nrow(genes), ncol=8)

    #loop by gene
    for(t in 1:nrow(genes)){
      chrom <- genes[t,1]
      genemax <- genes[t,3]
      genemin <- genes[t,2]
      gname <- genes[t,4]
      sub2 <- sub1[sub1[,"CHROM"] == chrom & sub1[,"POS"] <= genemax & sub1[,"POS"] >= genemin,]
      SNP_in_Gene <- nrow(sub2)
      
      results[t,1] <-i #enviroment
      results[t,2] <-as.character(chrom) #scaffold
      results[t,3] <-genemin #gene start position (bp) 
      results[t,4] <-genemax #gene end position (bp)
      results[t,5] <-SNP_in_Gene #number of SNPs in gene
      results[t,6] <-sum(sub2[,"pvalue"] < rho_threshold) # the number of outliers in the gene
      results[t,7] <-as.character(gname) #Gene name
    }
    colnames(results)<-res_name
    results1 <- data.frame(results)
    as.numeric.factor <- function(x) {as.numeric(levels(x))[x]} # convert factors to numeric function
    results2 <- results1[as.numeric.factor(results1$Observed_Outliers) > 0,]
    results2$Pout <- sum(as.numeric.factor(results2$Observed_Outliers)) / sum(as.numeric.factor(results2$SNPs_in_Gene)) #probability of outliers Pout    
    Res[[i]]<- results2
  }
  Res1 <- do.call(rbind, Res)
  Res2 <- Res1[complete.cases(Res1), ]

  # calculate expected number of outliers using bionomial distribution
  Res2$Expected_Outliers <- qbinom(binomial_cutoff, as.numeric.factor(Res2$SNPs_in_Gene), Res2$Pout)
  
  return(Res2)
}



# Run top candidate test
Cands <- top_candidate_test(rho, GA, 0.9999, 0.01)


write.table(Cands, "top_candidates_pooja", sep="\t", col.names = T, row.names = F, quote = F)
