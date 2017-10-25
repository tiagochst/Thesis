
# For each differently methylated probes we will get the
# 20 nearby genes (10 downstream and 10 upstream)
nearGenes <- GetNearGenes(data = mae,
                          probes =  diff.probes$probe,
                          numFlankingGenes = 20,
                          cores = 1)

Hypo.pair <- get.pair(data = mae,
                      nearGenes = nearGenes,
                      group.col = "definition",
                      group1 = "Primary solid Tumor",
                      group2 = "Solid Tissue Normal",
                      permu.dir = "Results_hypo/permu",
                      permu.size = 10000,
                      mode = "unsupervised",
                      minSubgroupFrac = 0.4,
                      raw.pvalue = 0.001,
                      Pe = 0.001,
                      filter.probes = TRUE,
                      filter.percentage = 0.05,
                      filter.portion = 0.3,
                      dir.out = "Results_hypo",
                      cores = 1,
                      label = "hypo")
