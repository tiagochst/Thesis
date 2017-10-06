
# Install funcivar package
if(!"funciVar" %in% installed.packages()[,1]) devtools::install_github("Simon-Coetzee/funcivar")
library(funciVar)
library(GenomicRanges)
library(dplyr)
# Load DNA methylation platform 450K manifest (hg38) and select only probes paired
distal.probes <- ELMER::get.feature.probe(feature = NULL,genome = "hg38", met.platform = "450K") 
getPair_hypo_pairs_significant <- readr::read_csv("Results_hypo/getPair.hypo.pairs.significant.csv")
paired.probes <- unique(getPair_hypo_pairs_significant$Probe)
paired.probes <- distal.probes[names(distal.probes) %in% paired.probes]

# Download state for breast cancer cell line (mcf-7)
base <- "http://s3-us-west-2.amazonaws.com/statehub-trackhub/tracks/5813b67f46e0fb06b493ceb0/hg38/ENCODE/"
# download tracks (search used: "encode hg38 h3k27ac h3k4me1 h3k4me3 ctcf")
state <- c("mcf-7.16mark.segmentation.bed",
           "bipolar_spindle_neuron.8mark.segmentation.bed",
           "cardiac_muscle_cell.9mark.segmentation.bed",
           "cd14-positive_monocyte.9mark.segmentation.bed",
           "dohh2.8mark.segmentation.bed",
           "fibroblast_of_dermis.8mark.segmentation.bed",
           "fibroblast_of_lung.13mark.segmentation.bed",
           "gm12878.12mark.segmentation.bed",
           "hct116.12mark.segmentation.bed",
           "hela-s3.13mark.segmentation.bed",
           "hepatocyte.9mark.segmentation.bed",
           "induced_pluripotent_stem_cell.7mark.segmentation.bed",
           "k562.19mark.segmentation.bed",
           "mcf-7.16mark.segmentation.bed",
           "neutrophil.8mark.segmentation.bed")

bed <- paste0(base,state)
dir.create("state_tracks", showWarnings = FALSE)
for( i in bed) if(!file.exists(file.path("state_tracks",basename(i)))) downloader::download(i,file.path("state_tracks",basename(i)))

esegs <- GetSegmentations(files =  dir("state_tracks",full.names = T)) %>% unlist

enrichmet <- CalculateEnrichment(variants = list(bg = distal.probes, fg = paired.probes),
                    features = esegs,
                    feature.type = "segmentations",
                    prior = c(a=0.5, b=0.5))

PlotEnrichment(variant.enrichment = enrichmet, 
               value = "difference", 
               block1 = "state", 
               color.by = "sample", 
               ncol = 3)