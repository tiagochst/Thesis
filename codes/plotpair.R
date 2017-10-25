# by region and without detail about DNA methylation
schematic.plot(data = mae,
               pair = Hypo.pair,
               byCoordinate = list(chr="chr1", start = 40000429, end = 42070429))
# by probe and with detail about DNA methylation
schematic.plot(data = mae,
               group.col = "definition",
               group1 = "Primary solid Tumor",
               group2 = "Solid Tissue Normal",
               pair = Hypo.pair,
               statehub.tracks = "hg38/ENCODE/mcf-7.16mark.segmentation.bed",
               byProbe = "cg04723436")

# by gene ID and adding annotation tracks from StateHub
schematic.plot(data = mae,
               group.col = "definition",
               group1 = "Primary solid Tumor",
               group2 = "Solid Tissue Normal",
               pair = Hypo.pair,
               byGeneID = "ENSG00000107485", # GATA3
               statehub.tracks = "hg38/ENCODE/mcf-7.16mark.segmentation.bed")
