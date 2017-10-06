
TF <- get.TFs(data = mae, 
              group.col = "definition",
              group1 = "Primary solid Tumor",
              group2 = "Solid Tissue Normal",
              minSubgroupFrac = 0.4, # Set to 1 if supervised mode
              enriched.motif = enriched.motif,
              dir.out = "Results_hypo", 
              cores = 1, 
              label = "hypo")

paste(sort(unique(TF$top.potential.TF.family)),collapse = ",")
# "EMX1,ESR1,FOXA1,GATA3,HOMEZ,LMX1B,MYB,MZF1,NR2F6,OVOL1,PBX1,RARA,SPDEF,ZKSCAN1,ZSCAN16"
paste(sort(unique(TF$top.potential.TF.subfamily)),collapse = ",")
# "AR,EMX1,FOXA1,FOXD2,GATA3,HOMEZ,LMX1B,MYB,NR2E3,PBX1,ZKSCAN1,ZSCAN16"