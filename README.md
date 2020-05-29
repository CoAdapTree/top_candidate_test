###Pooja Singh
###pooja.singh09@gmail.com
###CoAdapTree


Top_candidate_test as described in Yeaman, Hodgins et al 2016, Convergent local adaptation to climate in distantly related conifers)

top_candidates3_2env.R takes a gff annotation file and the output from the GEA spearmans rho analysis implemented in snp_env_association_spearmans3_parallelise.R for example:


snp     env     spearmansrho    pvalue
super4-15421    AHM     0.14909892997935        0.358498794832594
super4-47333    AHM     0.410221513046743       0.00856271980030036
super4-47418    AHM     -0.376483846516514      0.0166462386728531
super4-47461    AHM     -0.416326350480616      0.00753769112244098
super4-47513    AHM     0.278796986042015       0.0814941669997785
super4-47541    AHM     -0.271186302205992      0.0905260371101551

But it can be also be used with the output of a GWAS analysis which has columsn with scaffold-pos and some significance value.

top_candidates3_2env.R is run per environment and then the output is concatenated and used as input for top_candidates_select_plot.R. This produces the final list of singificant genes and top candidate style plots.


