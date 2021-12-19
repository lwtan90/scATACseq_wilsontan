# Analysis of scATAC-seq for Human Fetal Hearts  
Author: Wilson Tan
Date: Dec 19, 2021

## Load required libraries  
```
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
require(data.table)
set.seed(1234)
```  


## Read count matrix  
```
matrixfile = "GSE165837_CARE_ATAC_merged_matrix.mtx.gz" 
featurefile = "GSE165837_CARE_ATAC_merged_features.txt.gz"
barcodefile = "GSE165837_CARE_ATAC_merged_barcodes.txt.gz"
fragmentfile = "cleaned_sorted_fragment.bed.gz"
genome = "hg38"
counts = ReadMtx(mtx=matrixfile,features=featurefile,cells=barcodefile,feature.column=1)
```  
 
## Forming Seurat object  
```
chrom_assay <- CreateChromatinAssay(
   counts = counts,
   sep = c(":", "-"),
   genome = genome,
   fragments = fragmentfile,
   min.cells = 10,
   min.features = 200
)  

heart <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)
```  

## Perform annotation of Peaks  
```
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(heart) <- annotations

save(heart,file="step1.RData")
```  

# compute nucleosome signal score per cell
```
heart <- NucleosomeSignal(object = heart)

# compute TSS enrichment score per cell
heart <- TSSEnrichment(object = heart, fast = FALSE)
save(heart,file="step2.RData")

# Count total fragment
total_fragments <- CountFragments("cleaned_sorted_fragment.bed.gz")
ind = match(colnames(heart),total_fragments$CB)
heart$fragments = total_fragments$frequency_count[ind]
heart <- FRiP(object=heart,assay="peaks",total.fragments="fragments")
save(heart,file="step3.RData")

# Fraction of peaks in blacklisted region
heart$blacklist_fraction <- FractionCountsInRegion(
  object = heart, 
  assay = 'peaks',
  regions = blacklist_hg38
)

# add blacklist ratio and fraction of reads in peaks
heart$pct_reads_in_peaks <- heart$FRiP
heart$blacklist_ratio <- heart$blacklist_fraction
heart$high.tss <- ifelse(heart$TSS.enrichment > 1.5, 'High', 'Low')

heart@meta.data$newvar = rep("Y")
heart@meta.data$newvar[ is.na(heart@meta.data$fragments) == TRUE] = "U"
heart2 = subset(x = heart, subset = newvar == "Y")
 
heart2 <- NucleosomeSignal(object = heart2)
 
options(bitmapType='cairo')
png("TSSplot.png",width=4000,height=2500,res=300)
TSSPlot(heart2, group.by = 'high.tss') + NoLegend()
dev.off()

heart$nucleosome_group <- ifelse(heart$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
png("Fragment_histogram.png",width=2500,height=2500,res=300)
FragmentHistogram(object = heart2, group.by = 'nucleosome_group')
dev.off()

png("QC_violinplot.png",width=4000,height=2500,res=300)
VlnPlot(
  object = heart2,
  features = c('pct_reads_in_peaks', 'fragments', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

save(heart,file="preQC_heart.RData")

heart <- subset(
  x = heart,
  subset = fragments > 2000 &
    fragments < 100000 &
    FRiP > 0.5 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 5 &
    TSS.enrichment > 1.5
)
heart
```  

### Normalization  
```
heart <- RunTFIDF(heart)
heart <- FindTopFeatures(heart, min.cutoff = 'q0')
heart <- RunSVD(heart)

pdf("DepthCor.pdf",width=5,height=5)
DepthCor(heart)
dev.off()
```  

### Run UMAP  
```
heart <- RunUMAP(object = heart, reduction = 'lsi', dims = 2:30)
heart <- FindNeighbors(object = heart, reduction = 'lsi', dims = 2:30)
heart <- FindClusters(object = heart, verbose = FALSE, algorithm = 3)

options(bitmapType='cairo')
png("UMAP.png",width=2500,height=2300,res=300)
DimPlot(object = heart, label = TRUE) + NoLegend()
dev.off()

save(heart,file="UMAP.RData")
```  


### Gene Activities  
```
gene.activities <- GeneActivity(heart)
heart[['RNA']] <- CreateAssayObject(counts = gene.activities)
heart <- NormalizeData(
  object = heart,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(heart$nCount_RNA)
)

DefaultAssay(heart) <- 'RNA'

options(bitmapType='cairo')
png("scATACseq_featureplot.png",width=4000,height=4000,res=300)
FeaturePlot(
  object = heart,
  features = c('TNNT2', 'TTN', 'MYL2', 'CD34', 'CD14', 'MYH6', 'COL1A1','MYH7', 'MYL7'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
dev.off()

save(heart,file="gene.activities.RData")
```  

## Merge scRNA-seq and scATAC-seq
```
#### Load RNA first (and change to heart.rna)
heart.rna = readRDS("heart.rds")
heart.rna <- RenameIdents(
object = heart.rna,
"0"	= "ENDO",
"1"	= "FB",
"2"	= "FB",
"3"	= "Pericyte",
"4"	= "CM_MYL7",
"5"	= "ENDO",
"6"	= "CM",
"7"	= "FB",
"8" = "CM",
"9"	= "Macrophage",
"10" = "CM",
"11" = "CM_MYL7",
"12" =	"FB",
"13" = "CM_MYL7",
"14" = "MyoFB",
"15" = "CM",
"16" = "SMC",
"17" = "Neuron", 
"18" = "Adipose",
"19" = "Epicardial",
"20" = "Lymphocyte",
"21" = "Adipose"

)

#### Load ATAC-seq
load("heart.rna.RData")
load("gene.activities.RData")
load("heart.gene.RData")

heart[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(heart) <- "ACTIVITY"
heart <- NormalizeData(heart)
heart <- ScaleData(heart, features = rownames(heart))

transfer.anchors <- FindTransferAnchors(reference = heart.rna, query = heart, features = VariableFeatures(object = heart.rna), reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
save(transfer.anchors,file="transfer.anchors.RData")

heart.rna$celltype = Idents(heart.rna)
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = heart.rna$celltype, weight.reduction = heart[["lsi"]], dims = 2:30)

save(celltype.predictions,file="celltype.pred.RData")

heart <- AddMetaData(heart, metadata = celltype.predictions)


options(bitmapType='cairo')
png("lifedOver_UMAP.png",width=5000,height=2500,res=300)
p1 <- DimPlot(heart, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
    NoLegend() + scale_colour_hue(drop = FALSE)
p2 <- DimPlot(heart.rna, group.by = "celltype", label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + 
    NoLegend()
p1 | p2
dev.off()



png("prediction.score.png",width=2000,height=2000,res=300)
hist(heart$prediction.score.max)
abline(v = 0.5, col = "red")
dev.off()



heart.filtered <- subset(heart, subset = prediction.score.max > 0.5)
heart.filtered$predicted.id <- factor(heart.filtered$predicted.id, levels = levels(pbmc.rna))  # to make the colors match

png("lifedOver_UMAP_cleaned.png",width=5000,height=2500,res=300)
p1 <- DimPlot(heart.filtered, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + NoLegend() + scale_colour_hue(drop = FALSE)
p2 <- DimPlot(heart.rna, group.by = "celltype", label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + NoLegend()
p1 + p2
dev.off()


save(heart.filtered,file="final.ATACseq.RData")

heart = heart.filtered
```  


### Running Cicero  
```
library(monocle3)
require(cicero)
require(SeuratWrappers)
####remotes::install_github('satijalab/seurat-wrappers')

heart.cds <- as.cell_data_set(x = heart)
heart.cds = cluster_cells(cds = heart.cds,reduction_method="UMAP")
heart.cicero <- make_cicero_cds(heart.cds, reduced_coordinates = reducedDims(heart.cds)$UMAP)
save(heart.cicero,file="heart.cicero.RData")

genome = seqlengths(heart)
genome.df <- data.frame("chr" = names(genome), "length" = genome)
conns <- run_cicero(heart.cicero, genomic_coords = genome.df, sample_num = 100)

save(conns,file="conns.Rdata")

conns.filtered = conns[!is.na(conns$coaccess),]
dim(conns.filtered)

ccans <- generate_ccans(conns)
save(ccans,file="ccans.RData")

links <- ConnectionsToLinks(conns = conns, ccans = ccans)
DefaultAssay(heart) = "peaks"
Links(heart)=links
linked = data.frame(chr=seqnames(Links(heart)),start=start(ranges(Links(heart))),end=end(ranges(Links(heart))))
linked = cbind(linked,mcols(Links(heart)))


save(heart,file="link.added.heart.RData")
save(linked ,file="linked.RData")
```  


#### Call peaks
```
peaks <- CallPeaks(object = heart,group.by = "predicted.id",macs2.path = "macs2")
```

### Visualizing Peaks  
```
require(rtracklayer)
gene_anno = readGFF("gencode.v38.annotation.gtf")
gene_anno$chromosome = gene_anno$seqid
gene_anno$gene = gene_anno$gene_id
gene_anno$transcript = gene_anno$transcript_id
gene_anno$symbol = gene_anno$gene_name

options(bitmapType='cairo')
png("connection_TTN.png",width=3000,height=1500,res=300)
plot_connections(conns, "chr2", 178244554, 179088858,gene_model = gene_anno, coaccess_cutoff = 0.25, connection_width = 0.5, collapseTranscripts = "longest" )
dev.off()
```  



pos <- subset(gene_anno, strand == "+")
pos <- pos[order(pos$start),] 
# remove all but the first exons per transcript
pos <- pos[!duplicated(pos$transcript),] 
# make a 1 base pair marker of the TSS
pos$end <- pos$start + 1 

neg <- subset(gene_anno, strand == "-")
neg <- neg[order(neg$start, decreasing = TRUE),] 
# remove all but the first exons per transcript
neg <- neg[!duplicated(neg$transcript),] 
neg$start <- neg$end - 1

gene_annotation_sub <- rbind(pos, neg)

# Make a subset of the TSS annotation columns containing just the coordinates 
# and the gene name
gene_annotation_sub <- gene_annotation_sub[,c("chromosome", "start", "end", "symbol")]

# Rename the gene symbol column to "gene"
names(gene_annotation_sub)[4] <- "gene"

heart.cds <- annotate_cds_by_site(heart.cds, gene_annotation_sub)
heart.cds <- estimate_size_factors(heart.cds)
heart.cds <- preprocess_cds(heart.cds, method = "LSI")
heart.cds <- reduce_dimension(heart.cds, reduction_method = 'UMAP', preprocess_method = "LSI")
heart.cds <- cluster_cells(heart.cds)
heart.cds <- learn_graph(heart.cds)
# cell ordering can be done interactively by leaving out "root_cells"
heart.cds <- order_cells(heart.cds)





################### Differential ATAC-seq between clusters ###################
load("link.added.heart.RData")
Idents(heart) = heart$predicted.id
DefaultAssay(heart) <- 'peaks'
heart.markers <- FindAllMarkers(heart, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2)



###FeaturePlot
options(bitmapType='cairo')
plot1 <- VlnPlot(
  object = heart,
  features = "chr3-69828757-69829257",
  pt.size = 0.1
)
plot2 <- FeaturePlot(
  object = heart,
  features = "chr3-69828757-69829257",
  pt.size = 0.1
)

png("chr3-69828757-69829257_featureplot.png",width=2500,height=1500,res=300)
plot1 | plot2
dev.off()

# set plotting order
##levels(heart) <- c("CM","FB","ENDO","SMC")

options(bitmapType='cairo')
png("chr5-16539211-16539711_Region.png",width=5000,height=3000,res=300)
CoveragePlot(
  object = heart,
  region = "chr5-16539211-16539711",
  extend.upstream = 40000,
  extend.downstream = 40000
)
dev.off()
















### Rename Clusters 

heart <- RenameIdents(
 object = heart,
"0"="UNK2",
"1"="FB",
"2"="CM1",
"3"="CM1",
"4"="ENDO",
"5"="CM1",
"6"="UNK2",
"7"="FB",
"8"="CM3",
"9"="CM1",
"10"="CM3",
"11"="CM5",
"12"="ENDO",
"13"="CM6"
)

da_peaks_0_6 <- FindMarkers(
  object = heart,
  ident.1 = "0",
  ident.2 = "6",
  test.use = 'LR'
)
peaks_0_6 <- da_peaks_0_6[abs(da_peaks_0_6$avg_log2FC) > 0.3, ]
peaks_0_6 = peaks_0_6[ grep("KI",rownames(peaks_0_6), invert=TRUE),]
closest_genes_peaks_0_6 <- ClosestFeature(heart, regions = rownames(peaks_0_6) )
peaks_0_6 = cbind(peaks_0_6,closest_genes_peaks_0_6[ match(rownames(peaks_0_6),closest_genes_peaks_0_6$query_region), ])
peaks_0_6$comp = rep("0_6")
write.table(peaks_0_6,file="annotated_peaks_0_6.txt",sep="\t",quote=FALSE)
head(peaks_0_6)



da_peaks_0_4 <- FindMarkers(
  object = heart,
  ident.1 = "0",
  ident.2 = "4",
  test.use = 'LR'
)
peaks_0_4 <- da_peaks_0_4[abs(da_peaks_0_4$avg_log2FC) > 0.3, ]
peaks_0_4 = peaks_0_4[ grep("KI",rownames(peaks_0_4), invert=TRUE),]
closest_genes_peaks_0_4 <- ClosestFeature(heart, regions = rownames(peaks_0_4) )
peaks_0_4 = cbind(peaks_0_4,closest_genes_peaks_0_4[ match(rownames(peaks_0_4),closest_genes_peaks_0_4$query_region), ])
peaks_0_4$comp = rep("0_4")
write.table(peaks_0_4,file="annotated_peaks_0_4.txt",sep="\t",quote=FALSE)
head(peaks_0_4)


da_peaks_0_13 <- FindMarkers(
  object = heart,
  ident.1 = "0",
  ident.2 = "13",
  test.use = 'LR'
)
peaks_0_13 <- da_peaks_0_13[abs(da_peaks_0_13$avg_log2FC) > 0.3, ]
peaks_0_13 = peaks_0_13[ grep("KI",rownames(peaks_0_13), invert=TRUE),]
closest_genes_peaks_0_13 <- ClosestFeature(heart, regions = rownames(peaks_0_13) )
peaks_0_13 = cbind(peaks_0_13,closest_genes_peaks_0_13[ match(rownames(peaks_0_13),closest_genes_peaks_0_13$query_region), ])
peaks_0_13$comp = rep("0_13")
write.table(peaks_0_13,file="annotated_peaks_0_13.txt",sep="\t",quote=FALSE)
head(peaks_0_13)


da_peaks_0_1 <- FindMarkers(
  object = heart,
  ident.1 = "0",
  ident.2 = "1",
  test.use = 'LR'
)
peaks_0_1 <- da_peaks_0_1[abs(da_peaks_0_1$avg_log2FC) > 0.3, ]
peaks_0_1 = peaks_0_1[ grep("KI",rownames(peaks_0_1), invert=TRUE),]
peaks_0_1 = peaks_0_1[ grep("GL",rownames(peaks_0_1), invert=TRUE),]
closest_genes_peaks_0_1 <- ClosestFeature(heart, regions = rownames(peaks_0_1) )
peaks_0_1 = cbind(peaks_0_1,closest_genes_peaks_0_1[ match(rownames(peaks_0_1),closest_genes_peaks_0_1$query_region), ])
peaks_0_1$comp = rep("0_1")
write.table(peaks_0_1,file="annotated_peaks_0_1.txt",sep="\t",quote=FALSE)
head(peaks_0_1)




da_peaks_6_1 <- FindMarkers(
  object = heart,
  ident.1 = "6",
  ident.2 = "1",
  test.use = 'LR'
)
peaks_6_1 <- da_peaks_6_1[abs(da_peaks_6_1$avg_log2FC) > 0.3, ]
peaks_6_1 = peaks_6_1[ grep("KI",rownames(peaks_6_1), invert=TRUE),]
peaks_6_1 = peaks_6_1[ grep("GL",rownames(peaks_6_1), invert=TRUE),]
closest_genes_peaks_6_1 <- ClosestFeature(heart, regions = rownames(peaks_6_1) )
peaks_6_1 = cbind(peaks_6_1,closest_genes_peaks_6_1[ match(rownames(peaks_6_1),closest_genes_peaks_6_1$query_region), ])
peaks_6_1$comp = rep("6_1")
write.table(peaks_6_1,file="annotated_peaks_6_1.txt",sep="\t",quote=FALSE)
head(peaks_6_1)

da_peaks_2_1 <- FindMarkers(
  object = heart,
  ident.1 = "2",
  ident.2 = "1",
  test.use = 'LR'
)
peaks_2_1 <- da_peaks_2_1[abs(da_peaks_2_1$avg_log2FC) > 0.3, ]
peaks_2_1 = peaks_2_1[ grep("KI",rownames(peaks_2_1), invert=TRUE),]
peaks_2_1 = peaks_2_1[ grep("GL",rownames(peaks_2_1), invert=TRUE),]
closest_genes_peaks_2_1 <- ClosestFeature(heart, regions = rownames(peaks_2_1) )
peaks_2_1 = cbind(peaks_2_1,closest_genes_peaks_2_1[ match(rownames(peaks_2_1),closest_genes_peaks_2_1$query_region), ])
peaks_2_1$comp = rep("2_1")
write.table(peaks_2_1,file="annotated_peaks_2_1.txt",sep="\t",quote=FALSE)
head(peaks_2_1)

da_peaks_2_0 <- FindMarkers(
  object = heart,
  ident.1 = "2",
  ident.2 = "0",
  test.use = 'LR'
)
peaks_2_0 <- da_peaks_2_0[abs(da_peaks_2_0$avg_log2FC) > 0.3, ]
peaks_2_0 = peaks_2_0[ grep("KI",rownames(peaks_2_0), invert=TRUE),]
peaks_2_0 = peaks_2_0[ grep("GL",rownames(peaks_2_0), invert=TRUE),]
closest_genes_peaks_2_0 <- ClosestFeature(heart, regions = rownames(peaks_2_0) )
peaks_2_0 = cbind(peaks_2_0,closest_genes_peaks_2_0[ match(rownames(peaks_2_0),closest_genes_peaks_2_0$query_region), ])
peaks_2_0$comp = rep("2_0")
write.table(peaks_2_0,file="annotated_peaks_2_0.txt",sep="\t",quote=FALSE)
head(peaks_2_0)

da_peaks_3_0 <- FindMarkers(
  object = heart,
  ident.1 = "3",
  ident.2 = "0",
  test.use = 'LR'
)
peaks_3_0 <- da_peaks_3_0[abs(da_peaks_3_0$avg_log2FC) > 0.3, ]
peaks_3_0 = peaks_3_0[ grep("KI",rownames(peaks_3_0), invert=TRUE),]
peaks_3_0 = peaks_3_0[ grep("GL",rownames(peaks_3_0), invert=TRUE),]
closest_genes_peaks_3_0 <- ClosestFeature(heart, regions = rownames(peaks_3_0) )
peaks_3_0 = cbind(peaks_3_0,closest_genes_peaks_3_0[ match(rownames(peaks_3_0),closest_genes_peaks_3_0$query_region), ])
peaks_3_0$comp = rep("3_0")
write.table(peaks_3_0,file="annotated_peaks_3_0.txt",sep="\t",quote=FALSE)
head(peaks_3_0)

da_peaks_8_0 <- FindMarkers(
  object = heart,
  ident.1 = "8",
  ident.2 = "0",
  test.use = 'LR'
)
peaks_8_0 <- da_peaks_8_0[abs(da_peaks_8_0$avg_log2FC) > 0.3, ]
peaks_8_0 = peaks_8_0[ grep("KI",rownames(peaks_8_0), invert=TRUE),]
peaks_8_0 = peaks_8_0[ grep("GL",rownames(peaks_8_0), invert=TRUE),]
closest_genes_peaks_8_0 <- ClosestFeature(heart, regions = rownames(peaks_8_0) )
peaks_8_0 = cbind(peaks_8_0,closest_genes_peaks_8_0[ match(rownames(peaks_8_0),closest_genes_peaks_8_0$query_region), ])
peaks_8_0$comp = rep("8_0")
write.table(peaks_8_0,file="annotated_peaks_8_0.txt",sep="\t",quote=FALSE)
head(peaks_8_0)

da_peaks_15_0 <- FindMarkers(
  object = heart,
  ident.1 = "15",
  ident.2 = "0",
  test.use = 'LR'
)
peaks_15_0 <- da_peaks_15_0[abs(da_peaks_15_0$avg_log2FC) > 0.3, ]
peaks_15_0 = peaks_15_0[ grep("KI",rownames(peaks_15_0), invert=TRUE),]
peaks_15_0 = peaks_15_0[ grep("GL",rownames(peaks_15_0), invert=TRUE),]
closest_genes_peaks_15_0 <- ClosestFeature(heart, regions = rownames(peaks_15_0) )
peaks_15_0 = cbind(peaks_15_0,closest_genes_peaks_15_0[ match(rownames(peaks_15_0),closest_genes_peaks_15_0$query_region), ])
peaks_15_0$comp = rep("15_0")
write.table(peaks_15_0,file="annotated_peaks_15_0.txt",sep="\t",quote=FALSE)
head(peaks_15_0)


da_peaks_12_0 <- FindMarkers(
  object = heart,
  ident.1 = "12",
  ident.2 = "0",
  test.use = 'LR'
)
peaks_12_0 <- da_peaks_12_0[abs(da_peaks_12_0$avg_log2FC) > 0.3, ]
peaks_12_0 = peaks_12_0[ grep("KI",rownames(peaks_12_0), invert=TRUE),]
peaks_12_0 = peaks_12_0[ grep("GL",rownames(peaks_12_0), invert=TRUE),]
closest_genes_peaks_12_0 <- ClosestFeature(heart, regions = rownames(peaks_12_0) )
peaks_12_0 = cbind(peaks_12_0,closest_genes_peaks_12_0[ match(rownames(peaks_12_0),closest_genes_peaks_12_0$query_region), ])
peaks_12_0$comp = rep("12_0")
write.table(peaks_12_0,file="annotated_peaks_12_0.txt",sep="\t",quote=FALSE)
head(peaks_12_0)


da_peaks_3_4 <- FindMarkers(
  object = heart,
  ident.1 = "3",
  ident.2 = "4",
  test.use = 'LR'
)
peaks_3_4 <- da_peaks_3_4[abs(da_peaks_3_4$avg_log2FC) > 0.3, ]
peaks_3_4 = peaks_3_4[ grep("KI",rownames(peaks_3_4), invert=TRUE),]
peaks_3_4 = peaks_3_4[ grep("GL",rownames(peaks_3_4), invert=TRUE),]
closest_genes_peaks_3_4 <- ClosestFeature(heart, regions = rownames(peaks_3_4) )
peaks_3_4 = cbind(peaks_3_4,closest_genes_peaks_3_4[ match(rownames(peaks_3_4),closest_genes_peaks_3_4$query_region), ])
peaks_3_4$comp = rep("3_4")
write.table(peaks_3_4,file="annotated_peaks_3_4.txt",sep="\t",quote=FALSE)
head(peaks_3_4)

da_peaks_3_4 <- FindMarkers(object = heart,ident.1 = "3",ident.2 = "4",test.use = 'LR')
peaks_3_4 <- da_peaks_3_4[abs(da_peaks_3_4$avg_log2FC) > 0.3, ]
peaks_3_4 = peaks_3_4[ grep("KI",rownames(peaks_3_4), invert=TRUE),]
peaks_3_4 = peaks_3_4[ grep("GL",rownames(peaks_3_4), invert=TRUE),]
closest_genes_peaks_3_4 <- ClosestFeature(heart, regions = rownames(peaks_3_4) )
peaks_3_4 = cbind(peaks_3_4,closest_genes_peaks_3_4[ match(rownames(peaks_3_4),closest_genes_peaks_3_4$query_region), ])
peaks_3_4$comp = rep("3_4")
write.table(peaks_3_4,file="annotated_peaks_3_4.txt",sep="\t",quote=FALSE)
head(peaks_3_4)


DAframe = data.frame(p_val=c(),avg_log2FC=c(),pct.1=c(),pct.2=c(),p_val_adj=c(),tx_id=c(),gene_name=c(),gene_id=c(),gene_biotype=c(),type=c(),closest_region=c(),query_region=c(),distance=c(),comp=c())
for(i in 7:18)
{
	for(j in (i+1):18)
	{
		print(paste(i,j))
		da = FindMarkers(object = heart,ident.1 = i,ident.2 = j,test.use = 'LR')
		print(nrow(da))
		
		da = da[ grep("KI",rownames(da), invert=TRUE),]
		da = da[ grep("GL",rownames(da), invert=TRUE),]
		if(nrow(da)<2){
			next
		}
		print(dim(da))
		closest_genes_da <- ClosestFeature(heart, regions = rownames(da) )
		da = cbind(da,closest_genes_da[ match(rownames(da),closest_genes_da$query_region), ])
		print(dim(da))
		da=da[da$avg_log2FC>0,]
		if(nrow(da)<2){
			next
		}
		da$comp = rep(i)
		
		DAframe = rbind(DAframe,da)
	}
}

write.table(DAframe,file="DA_frame.txt",sep="\t",quote=FALSE,row.names=F)

DA = rbind(peaks_0_6,peaks_0_4,peaks_0_13,peaks_0_1,peaks_6_1,peaks_2_1,peaks_2_0,peaks_8_0,peaks_15_0,peaks_3_0,peaks_12_0,peaks_3_4)
write.table(DA,file="DA.txt",sep="\t",quote=FALSE,row.names=F)


da_peaks_CM1_FB <- FindMarkers(
  object = heart,
  ident.1 = "CM1",
  ident.2 = "FB",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'fragments'
)
peaks_CM1_FB <- da_peaks_CM1_FB[abs(da_peaks_CM1_FB$avg_log2FC) > 0.3, ]
closest_genes_peaks_CM1_FB <- ClosestFeature(heart, regions = rownames(peaks_CM1_FB) )
peaks_CM1_FB = cbind(peaks_CM1_FB,closest_genes_peaks_CM1_FB[ match(rownames(peaks_CM1_FB),closest_genes_peaks_CM1_FB$query_region), ])
peaks_CM1_FB$comp = rep("CM1_FB")
write.table(peaks_CM1_FB,file="annotated_peaks_CM1_FB.txt",sep="\t",quote=FALSE)
head(peaks_CM1_FB)



da_peaks_CM1_ENDO <- FindMarkers(
  object = heart,
  ident.1 = "CM1",
  ident.2 = "ENDO",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'fragments'
)
peaks_CM1_ENDO <- da_peaks_CM1_ENDO[abs(da_peaks_CM1_ENDO$avg_log2FC) > 0.3, ]
closest_genes_peaks_CM1_ENDO <- ClosestFeature(heart, regions = rownames(peaks_CM1_ENDO) )
peaks_CM1_ENDO = cbind(peaks_CM1_ENDO,closest_genes_peaks_CM1_ENDO[ match(rownames(peaks_CM1_ENDO),closest_genes_peaks_CM1_ENDO$query_region), ])
peaks_CM1_ENDO$comp = rep("CM1_ENDO")
write.table(peaks_CM1_ENDO,file="annotated_peaks_CM1_ENDO.txt",sep="\t",quote=FALSE)
head(peaks_CM1_ENDO)


da_peaks_FB_ENDO <- FindMarkers(
  object = heart,
  ident.1 = "FB",
  ident.2 = "ENDO",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'fragments'
)
peaks_FB_ENDO <- da_peaks_FB_ENDO[abs(da_peaks_FB_ENDO$avg_log2FC) > 0.3, ]
closest_genes_peaks_FB_ENDO <- ClosestFeature(heart, regions = rownames(peaks_FB_ENDO) )
peaks_FB_ENDO = cbind(peaks_FB_ENDO,closest_genes_peaks_FB_ENDO[ match(rownames(peaks_FB_ENDO),closest_genes_peaks_FB_ENDO$query_region), ])
peaks_FB_ENDO$comp = rep("FB_ENDO")
write.table(peaks_FB_ENDO,file="annotated_peaks_FB_ENDO.txt",sep="\t",quote=FALSE)
head(peaks_FB_ENDO)



da_peaks_CM1_UNK2 <- FindMarkers(
  object = heart,
  ident.1 = "CM1",
  ident.2 = "UNK2",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'fragments'
)
peaks_CM1_UNK2 <- da_peaks_CM1_UNK2[abs(da_peaks_CM1_UNK2$avg_log2FC) > 0.3, ]
closest_genes_peaks_CM1_UNK2 <- ClosestFeature(heart, regions = rownames(peaks_CM1_UNK2) )
peaks_CM1_UNK2 = cbind(peaks_CM1_UNK2,closest_genes_peaks_CM1_UNK2[ match(rownames(peaks_CM1_UNK2),closest_genes_peaks_CM1_UNK2$query_region), ])
peaks_CM1_UNK2$comp = rep("CM1_UNK2")
write.table(peaks_CM1_UNK2,file="annotated_peaks_CM1_UNK2.txt",sep="\t",quote=FALSE)
head(peaks_CM1_UNK2)



da_peaks_CM1_CM5 <- FindMarkers(
  object = heart,
  ident.1 = "CM1",
  ident.2 = "CM5",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'fragments'
)
peaks_CM1_CM5 <- da_peaks_CM1_CM5[abs(da_peaks_CM1_CM5$avg_log2FC) > 0.3, ]
closest_genes_peaks_CM1_CM5 <- ClosestFeature(heart, regions = rownames(peaks_CM1_CM5) )
peaks_CM1_CM5 = cbind(peaks_CM1_CM5,closest_genes_peaks_CM1_CM5[ match(rownames(peaks_CM1_CM5),closest_genes_peaks_CM1_CM5$query_region), ])
peaks_CM1_CM5$comp = rep("CM1_CM5")
write.table(peaks_CM1_CM5,file="annotated_peaks_CM1_CM5.txt",sep="\t",quote=FALSE)
head(peaks_CM1_CM5)


da_peaks_CM1_CM6 <- FindMarkers(
  object = heart,
  ident.1 = "CM1",
  ident.2 = "CM6",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'fragments'
)
peaks_CM1_CM6 <- da_peaks_CM1_CM6[abs(da_peaks_CM1_CM6$avg_log2FC) > 0.3, ]
closest_genes_peaks_CM1_CM6 <- ClosestFeature(heart, regions = rownames(peaks_CM1_CM6) )
peaks_CM1_CM6 = cbind(peaks_CM1_CM6,closest_genes_peaks_CM1_CM6[ match(rownames(peaks_CM1_CM6),closest_genes_peaks_CM1_CM6$query_region), ])
peaks_CM1_CM6$comp = rep("CM1_CM6")
write.table(peaks_CM1_CM6,file="annotated_peaks_CM1_CM6.txt",sep="\t",quote=FALSE)
head(peaks_CM1_CM6)



da_peaks_ENDO_CM5 <- FindMarkers(
  object = heart,
  ident.1 = "ENDO",
  ident.2 = "CM5",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'fragments'
)
peaks_ENDO_CM5 <- da_peaks_ENDO_CM5[abs(da_peaks_ENDO_CM5$avg_log2FC) > 0.3, ]
closest_genes_peaks_ENDO_CM5 <- ClosestFeature(heart, regions = rownames(peaks_ENDO_CM5) )
peaks_ENDO_CM5 = cbind(peaks_ENDO_CM5,closest_genes_peaks_ENDO_CM5[ match(rownames(peaks_ENDO_CM5),closest_genes_peaks_ENDO_CM5$query_region), ])
peaks_ENDO_CM5$comp = rep("ENDO_CM5")
write.table(peaks_ENDO_CM5,file="annotated_peaks_ENDO_CM5.txt",sep="\t",quote=FALSE)
head(peaks_ENDO_CM5)


da_peaks_ENDO_CM6 <- FindMarkers(
  object = heart,
  ident.1 = "ENDO",
  ident.2 = "CM6",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'fragments'
)
peaks_ENDO_CM6 <- da_peaks_ENDO_CM6[abs(da_peaks_ENDO_CM6$avg_log2FC) > 0.3, ]
closest_genes_peaks_ENDO_CM6 <- ClosestFeature(heart, regions = rownames(peaks_ENDO_CM6) )
peaks_ENDO_CM6 = cbind(peaks_ENDO_CM6,closest_genes_peaks_ENDO_CM6[ match(rownames(peaks_ENDO_CM6),closest_genes_peaks_ENDO_CM6$query_region), ])
peaks_ENDO_CM6$comp = rep("ENDO_CM6")
write.table(peaks_ENDO_CM6,file="annotated_peaks_ENDO_CM6.txt",sep="\t",quote=FALSE)
head(peaks_ENDO_CM6)



da_peaks_FB_CM5 <- FindMarkers(
  object = heart,
  ident.1 = "FB",
  ident.2 = "CM5",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'fragments'
)
peaks_FB_CM5 <- da_peaks_FB_CM5[abs(da_peaks_FB_CM5$avg_log2FC) > 0.3, ]
closest_genes_peaks_FB_CM5 <- ClosestFeature(heart, regions = rownames(peaks_FB_CM5) )
peaks_FB_CM5 = cbind(peaks_FB_CM5,closest_genes_peaks_FB_CM5[ match(rownames(peaks_FB_CM5),closest_genes_peaks_FB_CM5$query_region), ])
peaks_FB_CM5$comp = rep("FB_CM5")
write.table(peaks_FB_CM5,file="annotated_peaks_FB_CM5.txt",sep="\t",quote=FALSE)
head(peaks_FB_CM5)


da_peaks_FB_CM6 <- FindMarkers(
  object = heart,
  ident.1 = "FB",
  ident.2 = "CM6",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'fragments'
)
peaks_FB_CM6 <- da_peaks_FB_CM6[abs(da_peaks_FB_CM6$avg_log2FC) > 0.3, ]
closest_genes_peaks_FB_CM6 <- ClosestFeature(heart, regions = rownames(peaks_FB_CM6) )
peaks_FB_CM6 = cbind(peaks_FB_CM6,closest_genes_peaks_FB_CM6[ match(rownames(peaks_FB_CM6),closest_genes_peaks_FB_CM6$query_region), ])
peaks_FB_CM6$comp = rep("FB_CM6")
write.table(peaks_FB_CM6,file="annotated_peaks_FB_CM6.txt",sep="\t",quote=FALSE)
head(peaks_FB_CM6)

da_peaks_CM5_CM6 <- FindMarkers(
  object = heart,
  ident.1 = "CM5",
  ident.2 = "CM6",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'fragments'
)
peaks_CM5_CM6 <- da_peaks_CM5_CM6[abs(da_peaks_CM5_CM6$avg_log2FC) > 0.3, ]
closest_genes_peaks_CM5_CM6 <- ClosestFeature(heart, regions = rownames(peaks_CM5_CM6) )
peaks_CM5_CM6 = cbind(peaks_CM5_CM6,closest_genes_peaks_CM5_CM6[ match(rownames(peaks_CM5_CM6),closest_genes_peaks_CM5_CM6$query_region), ])
peaks_CM5_CM6$comp = rep("CM5_CM6")
write.table(peaks_CM5_CM6,file="annotated_peaks_CM5_CM6.txt",sep="\t",quote=FALSE)
head(peaks_CM5_CM6)

finalregion = c( rownames(peaks_CM1_CM3),rownames(peaks_CM1_FB),rownames(peaks_CM1_ENDO),rownames(peaks_FB_ENDO),rownames(da_peaks_CM1_UNK2), rownames(da_peaks_CM1_CM5), rownames(da_peaks_CM1_CM6), rownames(da_peaks_ENDO_CM5), rownames(da_peaks_ENDO_CM6), rownames(da_peaks_FB_CM5), rownames(da_peaks_FB_CM6), rownames(peaks_CM5_CM6))
finalregion = unique(finalregion)
 

### reperform UMAP
subset.matrix <- heart@assays$peaks@counts[finalregion, ] # Pull the raw expression matrix from the original Seurat object containing only the genes of interest
heart2 <- CreateSeuratObject(subset.matrix) # Create a new Seurat object with just the genes of interest
orig.ident <- Idents(heart) # Pull the identities from the original Seurat object as a data.frame
heart2 <- AddMetaData(object = heart2, metadata=heart@meta.data) # Add the idents to the meta.data slot
Idents(heart2) = orig.ident


heart2 <- RunTFIDF(heart2)
heart2 <- FindTopFeatures(heart2, min.cutoff = 'q0')
heart2 <- RunSVD(heart2)

heart2 <- RunUMAP(object = heart2,reduction = 'lsi',dims = 1:30)
heart2 <- FindNeighbors(object = heart2, reduction = 'lsi', dims = 2:30)
heart2 <- FindClusters(object = heart2, verbose = FALSE, algorithm = 3)


options(bitmapType='cairo')
png("UMAP.png",width=2500,height=2300,res=300)
DimPlot(object = heart, label = TRUE) + NoLegend()
dev.off()

options(bitmapType='cairo')
png("UMAP2.png",width=2500,height=2300,res=300)
DimPlot(object = heart2, label = TRUE) + NoLegend()
dev.off()

save(heart2,file="UMAP2.RData")


### cluster using only CM1 markers
finalregion = c( rownames(peaks_CM1_FB)[ abs(peaks_CM1_FB$avg_log2FC) > 0.2 ], rownames(peaks_CM1_ENDO)[ abs(peaks_CM1_ENDO$avg_log2FC) > 0.2 ])
subset.matrix <- heart@assays$peaks@counts[finalregion, ] # Pull the raw expression matrix from the original Seurat object containing only the genes of interest
heart2 <- CreateSeuratObject(subset.matrix) # Create a new Seurat object with just the genes of interest
orig.ident <- Idents(heart) # Pull the identities from the original Seurat object as a data.frame
heart2 <- AddMetaData(object = heart2, metadata=heart@meta.data) # Add the idents to the meta.data slot
Idents(heart2) = orig.ident


heart2 <- RunTFIDF(heart2)
heart2 <- FindTopFeatures(heart2, min.cutoff = 'q0')
heart2 <- RunSVD(heart2)

heart2 <- RunUMAP(object = heart2,reduction = 'lsi', dims = 1:30)
heart2 <- FindNeighbors(object = heart2, reduction = 'lsi', dims = 2:30)
heart2 <- FindClusters(object = heart2, verbose = FALSE, algorithm = 3)
Idents(heart2) = orig.ident

options(bitmapType='cairo')
png("UMAP2.png",width=2500,height=2300,res=300)
DimPlot(object = heart2, label = TRUE) + NoLegend()
dev.off()

save(heart2,file="UMAP2_CM1_FB_ENDO.RData")












heart2.new <- UpdateSeuratObject(object = heart2)

gene.activities <- GeneActivity(heart2)
gene.activities <- GeneActivity(heart2)
heart2[['RNA']] <- CreateAssayObject(counts = gene.activities)
heart2 <- NormalizeData(
  object = heart2,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(heart2$nCount_RNA)
)

DefaultAssay(heart) <- 'RNA'

options(bitmapType='cairo')
png("scATACseq_featureplot_umap2.png",width=4000,height=4000,res=300)
FeaturePlot(
  object = heart2,
  features = c('TNNT2', 'TTN', 'MYL2', 'CD34', 'CD14', 'MYH6', 'COL1A1','MYH7', 'MYL7'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
dev.off()


### extracting data
DApeak = heart@assays$peaks@data[ finalregion, ]
z = (DApeak-rowMeans(DApeak))/apply(DApeak,1,sd)

require(reshape)
require(ggplot2)


plotdata = melt(as.matrix(z))
names(plotdata)=c("peak","sample","z")
metadata = heart@meta.data
metadata$CELLTYPE = Idents(heart)
plotdata$CELLTYPE = metadata$CELLTYPE[match(plotdata$sample,rownames(metadata))]

options(bitmapType='cairo')
png("DE_heatmap.png",width=4000,height=4000,res=300)
p1<-ggplot(plotdata,aes(x=sample,y=peak))+geom_tile(aes(fill=z))+theme_bw()+theme(panel.grid=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.spacing=unit(0,"lines"))+facet_grid(de.group~cluster,space="free",scale="free")
p1<-p1+scale_fill_gradientn(values=rescale(c(-2,-1,-0.5,0,1,2)),colors=c("midnightblue","blue3","blue","turquoise","white","orange","red"))+xlab("")+ylab("")
print(p1)
dev.off()

 

#### Integration
load(file="gene.activity.RData") 
atac.heart = heart
DefaultAssay(atac.heart) <- 'RNA'
atac.heart <- RenameIdents(
object = atac.heart,
"0"="UNK2",
"1"="FB",
"2"="CM",
"3"="CM",
"4"="ENDO",
"5"="CM",
"6"="UNK2",
"7"="FB",
"8"="CM3",
"9"="CM",
"10"="CM3",
"11"="CM5",
"12"="ENDO",
"13"="CM6"
)


 rnaseq.heart = readRDS("../../../../scRNAseq/fetalHeart/heart_tutorial.rds")
 rnaseq.heart <- RenameIdents(
 object = rnaseq.heart,
"0"="CM",
"1"="FB2",
"2"="CM",	
"3"="CM",	
"4"="FB",	
"5"="CM",	
"6"="ENDO",	
"7"="ENDO",	
"8"="CM",	
"9"="MACRO",	
"10"="CM",	
"11"="PERI",	
"12"="MC",	
"13"="FB2",	
"14"="ENDO3",	
"15"="FB",	
"16"="NEU",	
"17"="HBF",	
"18"="ADI",	
"19"="MACRO",	
"20"="FB"

)

rnaseq.heart$celltype = Idents(rnaseq.heart)
options(bitmapType='cairo')
png("scRNAseq_UMAP.png",width=3000,height=3000,res=300)
plot1 <- DimPlot(
  object =rnaseq.heart,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot1 
dev.off()



rnaseq.heart$celltype = Idents(rnaseq.heart)

transfer.anchors <- FindTransferAnchors(
  reference = rnaseq.heart,
  query = atac.heart,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = rnaseq.heart$celltype,
  weight.reduction = atac.heart[['lsi']],
  dims = 2:30
)


atac.heart <- AddMetaData(object = atac.heart, metadata = predicted.labels)
save(atac.heart,transfer.anchors,predicted.labels,file="RNASEQ_ATACSEQ.RDATA")

options(bitmapType='cairo')
png("scRNAseq_atacseq_UMAP.png",width=3000,height=2000,res=300)
plot1 <- DimPlot(
  object =rnaseq.heart,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = atac.heart,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 + plot2
dev.off()














### Perform differential analysis  
DefaultAssay(heart) <- 'peaks'

da_peaks_CM_FB <- FindMarkers(
  object = heart,
  ident.1 = "CM",
  ident.2 = "FB",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'fragments'
)

plot1 <- VlnPlot(
  object = heart,
  features = rownames(da_peaks_CM_FB)[1],
  pt.size = 0.1
)
plot2 <- FeaturePlot(
  object = heart,
  features = rownames(da_peaks_CM_FB)[1],
  pt.size = 0.1
)

pdf("da_peaks_CM_FB.pdf",width=15,height=8)
plot1 | plot2
dev.off()

da_peaks_CM_FB.de <- rownames(da_peaks_CM_FB)
da_peaks_CM_FB.closest_genes <- ClosestFeature(heart, regions = da_peaks_CM_FB.de)

da_peaks_8_5 <- FindMarkers(
  object = heart,
  ident.1 = "5",
  ident.2 = "8",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'fragments'
)
peaks_8_5 <- da_peaks_8_5[abs(da_peaks_8_5$avg_log2FC) > 0.3, ]
closest_genes_peaks_8_5 <- ClosestFeature(heart, regions = rownames(peaks_8_5) )
peaks_8_5 = cbind(peaks_8_5,closest_genes_peaks_8_5[ match(rownames(peaks_8_5),closet_genes_peaks_8_5$query_region),] )
peaks_8_5 = cbind(peaks_8_5,closest_genes_peaks_8_5[ match(rownames(peaks_8_5),closest_genes_peaks_8_5$query_region), ])
write.table(peaks_8_5,file="annotated_peaks_10-5.txt",sep="\t",quote=FALSE)
head(peaks_8_5)

da_peaks_10_5 <- FindMarkers(
  object = heart,
  ident.1 = "5",
  ident.2 = "10",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'fragments'
)
peaks_10_5 <- da_peaks_10_5[abs(da_peaks_10_5$avg_log2FC) > 0.3, ]
closest_genes_peaks_10_5 <- ClosestFeature(heart, regions = rownames(peaks_10_5) )
peaks_10_5 = cbind(peaks_10_5,closest_genes_peaks_10_5[ match(rownames(peaks_10_5),closet_genes_peaks_10_5$query_region),] )
peaks_10_5 = cbind(peaks_10_5,closest_genes_peaks_10_5[ match(rownames(peaks_10_5),closest_genes_peaks_10_5$query_region), ])
write.table(peaks_10_5,file="annotated_peaks_10-5.txt",sep="\t",quote=FALSE)
head(peaks_10_5)

da_peaks_10_8 <- FindMarkers(
  object = heart,
  ident.1 = "8",
  ident.2 = "10",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'fragments'
)
options(bitmapType='cairo')
plot1 <- VlnPlot(
  object = heart,
  features = rownames(da_peaks_8_5)[1],
  pt.size = 0.1
)
plot2 <- FeaturePlot(
  object = heart,
  features = rownames(da_peaks_8_5)[1],
  pt.size = 0.1
)

png("da_peaks_8_5.png",width=2500,height=1500,res=300)
plot1 | plot2
dev.off()

da_peaks_ENDO_FB.de <- rownames(da_peaks_ENDO_FB)
da_peaks_ENDO_FB.closest_genes <- ClosestFeature(heart, regions = da_peaks_ENDO_FB.de)


da_peaks_ENDO_CM <- FindMarkers(
  object = heart,
  ident.1 = "ENDO",
  ident.2 = "CM",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'fragments'
)

plot1 <- VlnPlot(
  object = heart,
  features = rownames(da_peaks_ENDO_CM)[1],
  pt.size = 0.1
)
plot2 <- FeaturePlot(
  object = heart,
  features = rownames(da_peaks_ENDO_CM)[1],
  pt.size = 0.1
)

pdf("da_peaks_ENDO_CM.pdf",width=15,height=8)
plot1 | plot2
dev.off()


da_peaks_ENDO_CM.de <- rownames(da_peaks_ENDO_CM)
da_peaks_ENDO_CM.closest_genes <- ClosestFeature(heart, regions = da_peaks_ENDO_CM.de)
 


### Plotting UCSC Genome Browser  

# set plotting order
levels(heart) <- c("CM","FB","ENDO","SMC")

pdf("Example_Region.pdf",width=5,height=4)
CoveragePlot(
  object = heart,
  region = rownames(da_peaks_ENDO_CM)[1],
  extend.upstream = 100000,
  extend.downstream = 100000
)
dev.off()
