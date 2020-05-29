### PSingh JAN 2020 ####
### pooja.singh09@gmail.com ###
### CoAdapTree #####



### this script plots top candidates using the R graphics package ggplot2. It takes the output of the top_candidates.R script and plots the number of SNPs against the observed outliers for each gene. A single point represents a gene. A line representing the expected number of outliers predicted from a binomial distribution is overlayed.

options('stringsAsFactors'=FALSE)

library(dplyr)
library(ggplot2)

Cands <- read.table("top_candidates_split_snp_env_spearmans_rho_all.txt", header=T, sep="\t")


### plot and save to one file

plot_list = list()

for (i in unique(Cands$Environment)){
	sub <- subset(Cands, Environment == i)

	p = ggplot(data = sub %>% arrange(SNPs_in_Gene))+
  	geom_point(aes(SNPs_in_Gene,Observed_Outliers, colour = Observed_Outliers > Expected_Outliers), 
             size = 0.1)+
  	geom_line(aes(SNPs_in_Gene, Expected_Outliers), colour="firebrick")+
  	labs(x="Number of SNPs per gene", y="Number of outliers per gene")+
	ggtitle(paste0("Jackpine ",i, " B = 0.9999 Q = 0.01"))+
  	scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "grey"))+
  	theme_bw()+
    	ylim(c(0, max(as.numeric(Cands$Observed_Outliers))))+
    	xlim(c(0, max(as.numeric(Cands$SNPs_in_Gene))))+
    	theme(
      	panel.grid.major = element_blank(),
      	panel.grid.minor = element_blank(),
      	legend.position = "none"
   	)
	plot_list[[i]] = p
}


pdf("top_candidate_plots_rho.pdf")
for (i in unique(Cands$Environment)) {
    print(plot_list[[i]])
}
dev.off()

# Extract top candidates
Top_Cands <- Cands[Cands$Observed_Outliers > Cands$Expected_Outliers, c(1:9)]

# Save data
write.table(Top_Cands, file = "top_candidates_split_snp_env_spearmans_rho_all_final.txt", sep = "\t", row.names=FALSE, quote=FALSE)
