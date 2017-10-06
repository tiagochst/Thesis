
# For each differently methylated probes we will get the 
# 20 nearby genes (10 downstream and 10 upstream)
nearGenes <- GetNearGenes(data = mae, 
                          probes =  diff.probes$probe, 
                          numFlankingGenes = 20,
                          cores = 1)

# This step is the most time consuming. Depending on the size of the groups
# and the number of probes found previously it migh take hours
Hypo.pair <- get.pair(data = mae,
                      nearGenes = nearGenes,
                      group.col = "definition",
                      group1 = "Primary solid Tumor",
                      group2 = "Solid Tissue Normal",
                      permu.dir = "Results_hypo/permu",
                      permu.size = 10000, 
                      mode = "unsupervised",
                      minSubgroupFrac = 0.4, # 40% of samples to create U and M
                      raw.pvalue = 0.001,   
                      Pe = 0.001, 
                      filter.probes = TRUE,
                      filter.percentage = 0.05,
                      filter.portion = 0.3,
                      dir.out = "Results_hypo",
                      cores = 1,
                      label = "hypo")
# Number of pairs: 2950 