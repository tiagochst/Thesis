
metBoxPlot(mae,
           group.col = "definition",
           group1 = "Primary solid Tumor",
           group2 = "Solid Tissue Normal",
           diff.dir = "hypo",
           probe = "cg14058239",
           minSubgroupFrac = 0.2)
           
# Highlight subgroups in the colData matrix
# you can verify possible inputs with colnames(colData(mae))
metBoxPlot(mae,
           group.col = "definition",
           group1 = "Primary solid Tumor",
           group2 = "Solid Tissue Normal",
           diff.dir = "hypo",
           legend.col = "subtype_PAM50.mRNA",
           probe = "cg14058239",
           minSubgroupFrac = 0.2)             