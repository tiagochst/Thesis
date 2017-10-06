#--------------------------    DESCRIPTION    --------------------------
# Function to evaluate the overlaps between two genomic Interaction Ranges
# with start1/end1 and start2/end2
#--------------------------       STEPS       --------------------------
# => Step 1 - Compare gene-probe pairs found to chiapet
# We will download chiapet file and compare to the gene probe pairs found
# We check if there is an overlap between the range start1/end1 of both files
# and also start1/end1 to the start2/end2
# => Step 2 - Random analysis: compare random pairs to chiapet file
# Create a gene probe pair with the same number of size and pairs
# and do the same comparison
#-----------------------------------------------------------------------
library(GenomicRanges)
library(readr)
library(SummarizedExperiment)
library(ELMER)

# Auxiliary function
getOverlaps <- function(pairs,chiapet, genome = "hg19"){
  df <- NULL
  # We will need to create an ID for each interaction, as the liftover is one to many
  chiapet$ID <- rownames(chiapet)
  chiapet.left <- makeGRangesFromDataFrame(chiapet,start.field = "start1",end.field = "end1",seqnames.field = "chrom1",keep.extra.columns = T)
  chiapet.right <- makeGRangesFromDataFrame(chiapet,start.field = "start2",end.field = "end2",seqnames.field = "chrom2",keep.extra.columns = T)  
  if(genome == "hg38"){
    file <- "http://hgdownload.cse.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz"
    if(!file.exists(gsub("\\.gz","",basename(file)))) {
      downloader::download(file, basename(file))
      R.utils::gunzip(basename(file)) 
    }
    chain <- rtracklayer::import.chain(gsub("\\.gz","",basename(file)))
    chiapet.left <- unlist(rtracklayer::liftOver(chiapet.left,chain))
    chiapet.right <- unlist(rtracklayer::liftOver(chiapet.right,chain))
  }
  probes.gr <- makeGRangesFromDataFrame(pairs,start.field = "probe_start",end.field = "probe_end",seqnames.field = "probe_seqnames")
  genes.gr <- makeGRangesFromDataFrame(pairs,start.field = "gene_start_position",end.field = "gene_end_position",seqnames.field = "gene_seqnames")
  
  message("Find overlaps")
  # Left side probe, right side gene
  probe.overlap <- findOverlaps(probes.gr, chiapet.left)
  gene.overlap <- findOverlaps(genes.gr, chiapet.right)
  m1 <- chiapet.left[subjectHits(probe.overlap)]
  m2 <- chiapet.right[subjectHits(gene.overlap)] 
  merged.left <- merge(as.data.frame(m1), as.data.frame(m2), by = "ID")

  # Nb of overlaps btw probes and regions
  probe.overlap <- findOverlaps(probes.gr, chiapet.right)
  gene.overlap <- findOverlaps(genes.gr, chiapet.left)
  m1 <- chiapet.right[subjectHits(probe.overlap)]
  m2 <- chiapet.left[subjectHits(gene.overlap)]
  merged.right <- merge(as.data.frame(m1), as.data.frame(m2), by = "ID")
  
  merged <- rbind(merged.left, merged.right)
  merged <- merged[!duplicated(merged$ID),]
  
  nb.pairs <- nrow(unique(as.data.frame(pairs)[,c("Probe","GeneID")]))
  
  df <- data.frame(
    "Number pairs:" = nb.pairs,
    "Interaction in paper" = length(chiapet.left),
    "probe-gene that matched paper Interaction" = length(unique(merged$ID)), # we will consider only the gene instead of transcript 
    "percentage probe-gene that matched paper Interaction" = length(unique(merged$ID))/nb.pairs, 
    "Percentage of paper pairs infered" = length(unique(merged$ID))/(length(chiapet.left))
  )
  return(df)
}

#--------------------    MCF7 CHIAPET DOWNLOAD ----------------
# ChIA-PET Libraries MCF7_saturated - all interaction (enhancer-promoter, enhancer-enhancer, promoter-promoter)
# more info in http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S0092867411015170/1-s2.0-S0092867411015170-mmc1.pdf/272196/html/S0092867411015170/9e4e1c089f5568ede37f1a69a0597172/mmc1.pdf
file <- "http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S0092867411015170/1-s2.0-S0092867411015170-mmc2.xls/272196/html/S0092867411015170/15d30775b695c55d1f9fb145284c6876/mmc2.xls"
if(!file.exists(basename(file))) downloader::download(file, basename(file))
mmc2 <- readxl::read_excel(basename(file), sheet = "Table S3h. MCF7_saturated IXN", col_names = TRUE, skip = 3)
colnames(mmc2) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2","pet count", 
                    "p-value", "FDR", "overlap with structural variation from DNA-PET")
mcf7.sat.ixn <- mmc2

#--------------------   READ PAIRS  ------------------------
# Get pairs infered and add genomic positions
pairs <- read_csv("Results_hypo/getPair.hypo.pairs.statistic.with.empirical.pvalue.csv")
pairs <- pairs[pairs$Pe < 0.001 & pairs$Raw.p < 0.001 ,]

# ---  Map gene to its transcripts positions and probes to its ranges ---
# We will consider all transcripts (multiple for one gene)
tss <- as.data.frame(getTSS(genome = "hg38"))[,c("seqnames","strand","start_position","end_position","ensembl_gene_id")]
colnames(tss) <- paste0("gene_",colnames(tss))

data("hm450.manifest.hg38")
probes.ranges <- as.data.frame(hm450.manifest.hg38)[,1:5]
colnames(probes.ranges) <- paste0("probe_",colnames(probes.ranges))
probes.ranges$Probe <- rownames(probes.ranges)
pairs <- merge(pairs,tss,by.x = "GeneID", by.y = "gene_ensembl_gene_id")
pairs <- merge(pairs,probes.ranges,by = "Probe")

#---------------------------  Verify overlap  ----------------------------
df <- getOverlaps(pairs, mcf7.sat.ixn, "hg38")
print(df)
readr::write_csv(as.data.frame(df),path = "eval_near_genes_all_probes.csv")

# ------------------ Random comparison -------------------------------
# Create the random set for validation
distal.enhancer <- get.feature.probe(genome = "hg38", 
                                     met.platform = "450K", 
                                     feature = NULL) # get distal probes
distal.enhancer <- distal.enhancer[!names(distal.enhancer) %in% pairs$Probe,] # Select probes were not used
nb.pairs <- nrow(pairs)
nb.probes <- length(unique(pairs$Probe))
genes <- TCGAbiolinks:::get.GRCh.bioMart("hg38",as.granges = TRUE)
df.random <- NULL
for(i in 1:100){
  # We will get the double of random probes, because some will not be used in case it does not matches the position
  # Example: real probe + gene R10 and random probe does not have R10. Discart and get next random
  random.probes <- distal.enhancer[sample(1:length(distal.enhancer), nb.probes * 2),]
  near.genes <- GetNearGenes(TRange = random.probes, 
                             geneAnnot = genes, 
                             numFlankingGenes = 30)
  near.genes.df <- data.table::rbindlist(near.genes)
  # Now we should get the exactly same genes positions
  # if probe 1 was linked to R4 and L10 the random probe 1 will also be linked to its R4 and L10
  near.genes.linked <- NULL
  eval <- 1
  for(p in 1:length(near.genes)){
    side <- unique(pairs[pairs$Probe == unique(pairs$Probe)[eval],"Sides"])
    aux <-  near.genes.df[near.genes.df$Target == names(near.genes)[p],]
    same <- aux[aux$Side %in% side,]
    # If I do not have the same nearby positions
    if(length(side) != nrow(same)) next
    near.genes.linked <- rbind(near.genes.linked, same)
    eval <- eval+ 1
    if(length(unique(near.genes.linked$Target)) == nb.probes) break
  }
  print(length(unique(near.genes.linked$Probe)))
  colnames(near.genes.linked)[1] <- "Probe" 
  near.genes.linked <- merge(as.data.frame(near.genes.linked),tss,by.x = "GeneID", by.y = "gene_ensembl_gene_id")
  near.genes.linked <- merge(near.genes.linked,probes.ranges,by ="Probe")
  ret <- getOverlaps(near.genes.linked, mcf7.sat.ixn, "hg38")
  df.random <- rbind(df.random,ret)
  print(df.random)
}
save(df,df.random, file = "results_mcf7.rda")

# Plotting output
library(plyr)
library(ggthemes)
library(ggplot2)

data.random <- data.frame(mean(df.random$percentage.probe.gene.that.matched.paper.Interaction * df.random$Number.pairs.), sd(df.random$percentage.probe.gene.that.matched.paper.Interaction  * df.random$Number.pairs.), type = "random")
colnames(data.random) <- c("mean","sd","type")
data.real <- data.frame(mean(df$percentage.probe.gene.that.matched.paper.Interaction * df$Number.pairs.), sd(df$percentage.probe.gene.that.matched.paper.Interaction * df$Number.pairs.), type = "putative pair")
colnames(data.real) <- c("mean","sd","type")
data.real$sd <- NA
data <- rbind(data.random,data.real)

limits <- aes(ymax = mean + sd, ymin = mean - sd)
dodge <- position_dodge(width=0.9)
p <- ggplot(data, aes( fill = type, y=mean, x=type)) +
  geom_bar(position=dodge, stat="identity") +
  geom_errorbar(position=dodge,limits, width=0.25) +
  theme_solarized() +
  scale_colour_solarized("blue") + guides(fill=FALSE) +
  annotate("text", x=2, y=data$mean[2] + 5, label= paste0(data$mean[2],"(",round(data$mean[2]*100/unique(df.random$Number.pairs.),2),"%)")) + 
  annotate("text", x = 1, y=data$mean[1] + 30, label= paste0(data$mean[1],"(",round(data$mean[1]*100/unique(df.random$Number.pairs.),2),"%)")) +
  ylab("# of probe-gene pairs overlapping with MCF7 ChIA-PET interaction") +  xlab("") 
ggsave(p, filename = "validation.png")