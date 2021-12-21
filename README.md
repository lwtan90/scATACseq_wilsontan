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

## Link to Rdata  
Link to R temp data(https://drive.google.com/drive/folders/1p_sj0EcSj3jvbKIALuS0LtYoFBP1eh1L?usp=sharing).  

## Read count matrix  
```
setwd("C:/Users/wlwtan/Dropbox/My PC (CVI-CZL8ZC3)/Documents/tutorial_RData/")
matrixfile = "GSE165837_CARE_ATAC_merged_matrix.mtx.gz" 
featurefile = "GSE165837_CARE_ATAC_merged_features.txt.gz"
barcodefile = "GSE165837_CARE_ATAC_merged_barcodes.txt.gz"
fragmentfile = "cleaned_sorted_fragment.bed.gz"
genome = "hg38"
counts = ReadMtx(mtx=matrixfile,features=featurefile,cells=barcodefile,feature.column=1)

class(counts)
dim(counts)
```  
  
Output:  
```
class(counts)
[1] "dgCMatrix"
attr(,"package")
[1] "Matrix"

dim(counts)
[1] 501805  80083
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
  assay = "peaks"
)
```  
  
Output:  
```
> heart

An object of class Seurat 
501346 features across 79475 samples within 1 assay 
Active assay: peaks (501346 features, 0 variable features)
```  


## Adding annotation to chrom_assay  
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
  
Output:  
```
> head(Annotation(heart))
GRanges object with 6 ranges and 5 metadata columns:
                  seqnames        ranges strand |           tx_id   gene_name
                     <Rle>     <IRanges>  <Rle> |     <character> <character>
  ENSE00001489430     chrX 276322-276394      + | ENST00000399012      PLCXD1
  ENSE00001536003     chrX 276324-276394      + | ENST00000484611      PLCXD1
  ENSE00002160563     chrX 276353-276394      + | ENST00000430923      PLCXD1
  ENSE00001750899     chrX 281055-281121      + | ENST00000445062      PLCXD1
  ENSE00001489388     chrX 281192-281684      + | ENST00000381657      PLCXD1
  ENSE00001719251     chrX 281194-281256      + | ENST00000429181      PLCXD1
                          gene_id   gene_biotype     type
                      <character>    <character> <factor>
  ENSE00001489430 ENSG00000182378 protein_coding     exon
  ENSE00001536003 ENSG00000182378 protein_coding     exon
  ENSE00002160563 ENSG00000182378 protein_coding     exon
  ENSE00001750899 ENSG00000182378 protein_coding     exon
  ENSE00001489388 ENSG00000182378 protein_coding     exon
  ENSE00001719251 ENSG00000182378 protein_coding     exon
  -------
  seqinfo: 25 sequences from hg38 genome
> length(Annotation(heart))
[1] 3021151

```  
 

## compute nucleosome signal score per cell
```
heart <- NucleosomeSignal(object = heart)
```  
Output:  
![alt text](https://github.com/lwtan90/scATACseq_wilsontan/blob/main/img/nucleosome_pattern%20copy.png)  



## compute TSS enrichment score per cell  
```
heart <- TSSEnrichment(object = heart, fast = FALSE)
save(heart,file="step2.RData")
```  
Output:
![alt text](https://github.com/lwtan90/scATACseq_wilsontan/blob/main/img/TSSenrichment.png)  
  
## other QCs  
```
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
```  
Output:  
![alt text](https://github.com/lwtan90/scATACseq_wilsontan/blob/main/img/QC_violinplot.png)

## FIltering low quality cells  
```
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

png("DepthCor.png",width=2000,height=2000,res=300)
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

Output:  
![alt text](https://github.com/lwtan90/scATACseq_wilsontan/blob/main/img/UMAP.png)


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
Note: we can use scRNA-seq to infer cell type, and subsequently perform label transfer to annotate scATAC-seq.  
This RNA-seq data has been processed using seurat, and stored in heart.rds.  

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
save(heart.rna,file="heart.rna.RData")
```  

## Load scATAC-seq Data  
```
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
![alt text](https://github.com/lwtan90/scATACseq_wilsontan/blob/main/img/lifedOver_UMAP_cleaned.png)

## Running Cicero  

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



## Differential ATAC-seq between clusters  
```
load("link.added.heart.RData")
Idents(heart) = heart$predicted.id
DefaultAssay(heart) <- 'peaks'

## This is one quick way to get all markers for all cell type
heart.markers <- FindAllMarkers(heart, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2)

## This is another way to get pairwise comparison
da_peaks_CM_FB <- FindMarkers(
  object = heart,
  ident.1 = "CM",
  ident.2 = "FB",
  test.use = 'LR'
)
##closest_genes_peaks_CM_FB <- ClosestFeature(heart, regions = rownames(da_peaks_CM_FB) )
##da_peaks_CM_FB = cbind(da_peaks_CM_FB,closest_genes_peaks_CM_FB[ match(rownames(da_peaks_CM_FB),closest_genes_peaks_CM_FB$query_region), ])
##da_peaks_CM_FB$comp = rep("0_4")
##write.table(da_peaks_CM_FB,file="annotated_peaks_CM_FB.txt",sep="\t",quote=FALSE)
```  


## FeaturePlot  
```
options(bitmapType='cairo')  
plot1 <- VlnPlot(
  object = heart,
  features = "chr1-16017879-16018379",
  pt.size = 0.1
)
plot2 <- FeaturePlot(
  object = heart,
  features = "chr1-16012798-16013298",
  pt.size = 0.1
)

png("chr3-69828757-69829257_featureplot.png",width=2500,height=1500,res=300)
plot1 | plot2
dev.off()


options(bitmapType='cairo')
png("chr3-69828757-69829257_Region.png",width=5000,height=3000,res=300)
CoveragePlot(
  object = heart,
  region = "chr3-69828757-69829257",
  extend.upstream = 40000,
  extend.downstream = 40000
)
dev.off()
```  
![alt text](https://github.com/lwtan90/scATACseq_wilsontan/blob/main/img/chr3-69828757-69829257_Region.png)  
![alt text](https://github.com/lwtan90/scATACseq_wilsontan/blob/main/img/chr3-69828757-69829257_featureplot.png)  


## Motif Analysis  
```
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
require(motifmatchr)

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
heart <- AddMotifs(
  object = heart,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks_CM_FB[da_peaks_CM_FB$p_val < 0.005, ])
enriched.motifs <- FindMotifs(
  object = heart,
  features = top.da.peak
)

MotifPlot(object = heart,motifs = head(rownames(enriched.motifs)))

```  


## Run ChroomVAR  
```
heart <- RunChromVAR(
  object = mouse_brain,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(heart) <- 'chromvar'

# look at the activity of Mef2c
differential.activity <- FindMarkers(
  object = heart,
  ident.1 = 'CM',
  ident.2 = 'FB',
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

MotifPlot(
  object = heart,
  motifs = head(rownames(differential.activity)),
  assay = 'peaks'
)
```
