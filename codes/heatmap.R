# Get molecular subytpes/purity from cell paper
# https://doi.org/10.1016/j.cell.2015.09.033
file <- "http://ars.els-cdn.com/content/image/1-s2.0-S0092867415011952-mmc2.xlsx"
downloader::download(file, basename(file))
subtypes <- readxl::read_excel(basename(file), skip = 2)

subtypes$sample <- substr(subtypes$Methylation,1,16)
meta.data <- merge(colData(mae),subtypes,by = "sample",all.x = T)
meta.data <- meta.data[match(colData(mae)$sample,meta.data$sample),]
meta.data <- S4Vectors::DataFrame(meta.data)
rownames(meta.data) <- meta.data$sample
stopifnot(all(meta.data$patient == colData(mae)$patient))
colData(mae) <- meta.data

heatmapPairs(data = mae, 
             group.col = "definition",
             group1 = "Primary solid Tumor",
             group2 = "Solid Tissue Normal",
             annotation.col = c("TumorPurity",
                                "Final.Pathology",
                                "PR.IHC",
                                "HER2.IHC",
                                "ER.IHC",
                                "PAM50"),
             pairs = Hypo.pair,
             filename = "heatmap.pdf")