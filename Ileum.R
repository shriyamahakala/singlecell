library(Seurat)
library(dplyr)
library(ggplot2)
library(devtools)
library(dittoSeq)
library(ggpubr)
library(ggtree)
library(tidyverse)
library(ggdendro)
library(cowplot)
library(patchwork) 
library(hdf5r)
library(DropletUtils)
library(rhdf5)

# Load the PBMC dataset
GCA40_I.data <- Read10X(data.dir = "C:/Users/shriy/Documents/Research/GCA40_I/filtered_feature_bc_matrix/")

#need access to with_introns 
#GCA40_I.data <- Read10X(data.dir = "~/OneDrive - Emory University/Hemsley_renewal/batch3/with_introns/GCA40_I/filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data).
GCA40_I <- CreateSeuratObject(counts = GCA40_I.data, project = "GCA40_I")
GCA40_I
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
GCA40_I[["percent.mt"]] <- PercentageFeatureSet(GCA40_I, pattern = "^MT-")
GCA40_I[["percent.ribo"]] <- PercentageFeatureSet(GCA40_I, pattern = "^RP[SL]")
# Visualize QC metrics as a violin plot
VlnPlot(GCA40_I, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)

# Load the PBMC dataset
GCA42_I.data <- Read10X(data.dir = "C:/Users/shriy/Documents/Research/GCA42_I/filtered_feature_bc_matrix/")
#GCA42_I.data <- Read10X(data.dir = "~/OneDrive - Emory University/Hemsley_renewal/batch3/with_introns/GCA42_I/filtered_feature_bc_matrix/")

# Initialize the Seurat object with the raw (non-normalized data).
GCA42_I <- CreateSeuratObject(counts = GCA42_I.data, project = "GCA42_I")
GCA42_I
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
GCA42_I[["percent.mt"]] <- PercentageFeatureSet(GCA42_I, pattern = "^MT-")
GCA42_I[["percent.ribo"]] <- PercentageFeatureSet(GCA42_I, pattern = "^RP[SL]")
# Visualize QC metrics as a violin plot
VlnPlot(GCA42_I, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)

# Load the PBMC dataset
GCA43_I.data <- Read10X(data.dir = "C:/Users/shriy/Documents/Research/GCA43_I/filtered_feature_bc_matrix/")
#GCA43_I.data <- Read10X(data.dir = "~/OneDrive - Emory University/Hemsley_renewal/batch3/with_introns/GCA43_I/filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data).
GCA43_I <- CreateSeuratObject(counts = GCA43_I.data, project = "GCA43_I")
GCA43_I
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
GCA43_I[["percent.mt"]] <- PercentageFeatureSet(GCA43_I, pattern = "^MT-")
GCA43_I[["percent.ribo"]] <- PercentageFeatureSet(GCA43_I, pattern = "^RP[SL]")
# Visualize QC metrics as a violin plot
VlnPlot(GCA43_I, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)

# Load the PBMC dataset
GCA44_I.data <- Read10X(data.dir = "C:/Users/shriy/Documents/Research/GCA44_I/filtered_feature_bc_matrix/")
#GCA44_I.data <- Read10X(data.dir = "~/OneDrive - Emory University/Hemsley_renewal/batch3/with_introns/GCA44_I/filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data).
GCA44_I <- CreateSeuratObject(counts = GCA44_I.data, project = "GCA44_I")
GCA44_I
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
GCA44_I[["percent.mt"]] <- PercentageFeatureSet(GCA44_I, pattern = "^MT-")
GCA44_I[["percent.ribo"]] <- PercentageFeatureSet(GCA44_I, pattern = "^RP[SL]")
# Visualize QC metrics as a violin plot
VlnPlot(GCA44_I, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)

#---------------------------------------------------------------------------
#Run CCA integration--> after merging 8 datasets
#---------------------------------------------------------------------------
prot.list<-list(GCA40_I,GCA42_I,GCA43_I, GCA44_I)
#prot.list <- SplitObject(prot, split.by = "orig.ident")
prot.list <- lapply(X = prot.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = prot.list)
prot.list <- lapply(X = prot.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
data.anchors <- FindIntegrationAnchors(object.list = prot.list, anchor.features = features, reduction = "rpca")
#this command creates an 'integrated' data assay

gc()

data.combined <- IntegrateData(anchorset = data.anchors)
# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(data.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
data.combined <- ScaleData(data.combined, verbose = FALSE)
data.combined <- RunPCA(data.combined, npcs = 30, verbose = FALSE)
data.combined <- RunUMAP(data.combined, reduction = "pca", dims = 1:30)
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:30)
data.combined <- FindClusters(data.combined, resolution = 0.5)
# Visualization
DimPlot(data.combined, reduction = "umap", label = TRUE)
DimPlot(data.combined, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE)
#Add disease condition to the metadata
DefaultAssay(data.combined) <- "RNA"


#---------------------------------------------------------------------------
#Run merge #no integration
#---------------------------------------------------------------------------
prot_merge <- merge(GCA40_I, y = c(GCA42_I, GCA43_I,GCA44_I), add.cell.ids = c("GCA40_I", "GCA42_I", "GCA43_I","GCA44_I"), project = "prot_merge")
prot_merge
pbmc=prot_merge
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), npcs = 50)
DimPlot(pbmc, reduction = "pca")

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100, dims = 50)
pbmc <- ScoreJackStraw(pbmc, dims = 1:50)
JackStrawPlot(pbmc, dims = 1:50)
ElbowPlot(pbmc,ndims = 50)

pbmc <- FindNeighbors(pbmc, dims = 1:50)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:50)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap", label = TRUE)


# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc=data.combined
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(pbmc.markers, file="~/OneDrive - Emory University/Hemsley_renewal/markers/not_integrated_Ileum_batch3.csv")
write.csv(pbmc.markers, file="~/OneDrive - Emory University/Hemsley_renewal/markers/integrated_Ileum_batch3.csv")

write.csv(pbmc.markers, file="~/OneDrive - Emory University/Hemsley_renewal/markers/With_introns_not_integrated_Ileum_batch3.csv")
write.csv(pbmc.markers, file="~/OneDrive - Emory University/Hemsley_renewal/markers/With_introns_integrated_Ileum_batch3.csv")

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
write.csv(top10, file="~/OneDrive - Emory University/Hemsley_renewal/markers/top10_With_introns_not_integrated_Ileum_batch3.csv")
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

saveRDS(pbmc, file = "~/OneDrive - Emory University/Hemsley_renewal/seurat/not_integrated_Ileum_batch3.rds")
saveRDS(pbmc, file = "~/OneDrive - Emory University/Hemsley_renewal/seurat/integrated_Ileum_batch3.rds")

saveRDS(pbmc, file = "~/OneDrive - Emory University/Hemsley_renewal/seurat/With_introns_not_integrated_Ileum_batch3.rds")
saveRDS(pbmc, file = "~/OneDrive - Emory University/Hemsley_renewal/seurat/With_introns_integrated_Ileum_batch3.rds")



#-------------------------------------------------------------
#UMAP based mt%
#-------------------------------------------------------------
prot3=data.combined
temp1=as.data.frame(Cells(prot3))
temp1$mtf=prot3$percent.mt
#head(temp1)
#setwd("F:/Hemsley/Rough_data/")
#write.csv(temp1, file="temp1.csv")
#temp1=read.csv(file="temp1.csv", header = TRUE)
#head(temp1)
#temp<-temp5 %>% mutate(group=if_else(mtf < 5, "less5p","g5p"))
temp<-temp1 %>% mutate(group=case_when(mtf <= 30~ "lessthan30",mtf >=30 & mtf <= 50 ~ "30to50", TRUE ~  "greaterthan50" ))
#temp<-temp5 %>% mutate(group=case_when(mtf < 5~ "less5p",mtf >=5 & mtf < 20 ~ "g5p", mtf >=20 & mtf < 40 ~ "g2040p",mtf >=40 & mtf < 60 ~ "g4060p",mtf >=60 & mtf < 80 ~ "g6080p",   TRUE ~ "high" ))
prot31=prot3
Idents(prot31)
Idents(prot31)=temp$group
DimPlot(prot31, group.by = "ident")


#-------------------------------------------------------------
#read rds objects Ileum_batch3 --> downstream analysis of gene lists 
#-------------------------------------------------------------

Without_introns_not_integrated<-readRDS("~/OneDrive - Emory University/Hemsley_renewal/seurat/not_integrated_Ileum_batch3.rds")
Without_introns_integrated<-readRDS("~/OneDrive - Emory University/Hemsley_renewal/seurat/integrated_Ileum_batch3.rds")

With_introns_not_integrated<-readRDS("~/OneDrive - Emory University/Hemsley_renewal/seurat/With_introns_not_integrated_Ileum_batch3.rds")
With_introns_integrated<-readRDS("~/OneDrive - Emory University/Hemsley_renewal/seurat/With_introns_integrated_Ileum_batch3.rds")

#high UMI
highUMI<-subset(Without_introns_integrated, subset=nCount_RNA>170000)
highUMI
counts<-as.data.frame(highUMI@assays$RNA@counts)
counts<-counts[rowSums(counts[])>0,]
dim(counts)
#write.csv(counts,"~/OneDrive - Emory University/Hemsley_renewal/testing/original/counts_highUMI_Without_introns_integrated.csv")
cell1<-subset(counts, select="GAGACCCAGATGGTAT-1_1")
cell2<-subset(counts, select="AGAAGTAAGACTTGTC-1_3")
cell3<-subset(counts, select="TACCGAACAAGCTGCC-1_3")
cell4<-subset(counts, select="TCACTCGGTATTGAGA-1_3")
cell5<-subset(counts, select="TCAGCCTCACTAGGTT-1_3")
cell6<-subset(counts, select="TCATGCCTCATCGCAA-1_3")
cell7<-subset(counts, select="TCGCACTTCCCAGTGG-1_3")
cell8<-subset(counts, select="TTACGCCGTCAGTCCG-1_3")
cell9<-subset(counts, select="TCAAGTGAGTGCTAGG-1_4")
cell10<-subset(counts, select="TGCATGACAATGACCT-1_4")
cell11<-subset(counts, select="TGTTGGAAGGCTGAAC-1_4")
cell11<-filter(cell11, cell11$`TGTTGGAAGGCTGAAC-1_4`>0)
write.csv(cell11,"~/OneDrive - Emory University/Hemsley_renewal/testing/original/cell11_highUMI.csv")

#high UMI
lowUMI<-subset(Without_introns_integrated, subset=nCount_RNA<501)
lowUMI
counts<-as.data.frame(lowUMI@assays$RNA@counts)
counts<-counts[rowSums(counts[])>0,]
dim(counts)
#write.csv(counts,"~/OneDrive - Emory University/Hemsley_renewal/testing/original/counts_lowUMI_Without_introns_integrated.csv")

#epi_PTPRC 
epi_PTPRC<-subset(Without_introns_integrated, subset=EPCAM>0 & PTPRC>0)#dp+ cells
epi_PTPRC
counts_epi_PTPRC<-as.data.frame(epi_PTPRC@assays$RNA@counts)
counts_epi_PTPRC<-counts_epi_PTPRC[1:10]#first 10 dp+ cells
counts_epi_PTPRC<-counts_epi_PTPRC[rowSums(counts_epi_PTPRC[])>0,]#10 randon dp+ cells with counts>0
dim(counts_epi_PTPRC)
row_dp<-as.data.frame(rownames(counts_epi_PTPRC))
#write.csv(counts_epi_PTPRC,"~/OneDrive - Emory University/Hemsley_renewal/testing/original/first_10_cells_counts_epi_PTPRC_Without_introns_integrated.csv")

#epi
epi<-subset(Without_introns_integrated, subset=EPCAM>0 & PTPRC<=0)
epi
counts_epi1<-as.data.frame(epi@assays$RNA@counts)
counts_epi1<-counts_epi1[1:10]#first 10 epi+ cells
counts_epi1<-counts_epi1[rowSums(counts_epi1[])>0,]#10 randon dp+ cells with counts>0
dim(counts_epi1)
row_epi<-as.data.frame(rownames(counts_epi1))
#write.csv(counts_epi1,"~/OneDrive - Emory University/Hemsley_renewal/testing/original/first_10_cells_counts_epi_Without_introns_integrated.csv")

#PTPRC 
PTPRC<-subset(Without_introns_integrated, subset=EPCAM<=0 & PTPRC>0)
PTPRC
counts_PTPRC1<-as.data.frame(PTPRC@assays$RNA@counts)
counts_PTPRC1<-counts_PTPRC1[1:10]#10 randon dp+ cells
counts_PTPRC1<-counts_PTPRC1[rowSums(counts_PTPRC1[])>0,]#10 randon dp+ cells with counts>0
dim(counts_PTPRC1)
#write.csv(counts_PTPRC1,"~/OneDrive - Emory University/Hemsley_renewal/testing/original/first_10_cells_counts_PTPRC_Without_introns_integrated.csv")

#compare the counts1 gene list to counts_epi1 gene list
t<-data.frame(rownames(counts_epi1))
t1<-data.frame(rownames(counts_epi_PTPRC))
t2<-as.data.frame(rownames(counts_PTPRC1))
colnames(t)[1]<-"genes"
colnames(t1)[1]<-"genes"
colnames(t2)[1]<-"genes"

tf<-inner_join(t1,t)#7051
tf2<-inner_join(t1,t2)
tf<-data.frame(intersect(t$genes, t1$genes))#7051
tf2<-as.data.frame(intersect(t1$`rownames(counts_epi_PTPRC)`, t2$`rownames(counts_PTPRC1)`))#7050

uni<-merge(t,t1)
uni<-as.data.frame(t$`rownames(counts_epi1)`%notin%t1$`rownames(counts1)`)
uni<-as.data.frame(uni)
uni1<-uni[is.na(t1$`rownames(counts1)`, )]


counts<-as.data.frame(epi_PTPRC@assays$RNA@counts)
write.csv(counts,file = "~/OneDrive - Emory University/Hemsley_renewal/markers/epi_counts_Without_introns_integrated.csv" )

molecule_info<-"~/OneDrive - Emory University/Hemsley_renewal/batch3/with_introns/GCA40_I/molecule_info.h5"
filtered_feature_bc_matrix<-"~/OneDrive - Emory University/Hemsley_renewal/batch3/without_introns/GCA40_I/filtered_feature_bc_matrix.h5"
molecule_info<- read10xMolInfo(molecule_info)
filtered_feature_bc_matrix<- h5ls(filtered_feature_bc_matrix)
molecule_info


table(Without_introns_integrated$nCount_RNA>0)



draw.pairwise.venn(area1 = 8073,                        # Create pairwise venn diagram
                   area2 = 8770,
                   cross.area = 5936,
                   fill = c("violetred", "slateblue2"),
                   lty = "blank",
                   cex=1)

t<-as.vector(rownames(counts_epi1))
t1<-as.vector(rownames(counts1))


rn<-union(t,t1)
rn<-as.data.frame(rn)
temp<-temp1 %>% mutate(group=case_when(mtf <= 30~ "lessthan30",mtf >=30 & mtf <= 50 ~ "30to50", TRUE ~  "greaterthan50" ))

rn1<-rn %>% mutate(group=case_when(rn %in% t~"t",rn %in% t1 ~ "t1"))
rn2<-ifelse(rn1$group%in% t,1,0)
rn2<-as.data.frame(rn2)
rn1$value<-rn2$rn2
t<-c("t")
t1<-c("t1")
rn2<-rn1%>% mutate(group=case_when(rn1$group %in% t ~"1",rn1$group %in% t1 ~ "2"))


dat <- data.frame(t = as.integer(rn %in% t),
                  t1 = as.integer(rn %in% t1)
)

row.names(dat) <- rn

dat
library(tidyr)
library(dplyr)
library(ggplot2)
dat$genes<-rownames(dat)
dat %>% 
  pivot_longer(-1) %>%
  mutate(genes = factor(genes, paste0("genes", 1:12228))) %>%
  ggplot(aes(dat, genes, fill = value)) + geom_tile(color = "black") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))
dat1<-as.matrix(dat)
dat2<-t(dat1)
dat2<-dat2[-3,]
dat2<-as.data.frame(dat2)
c<-as.data.frame(colnames(dat2))
r<-as.data.frame(rownames(dat2))
ggplot(rn1,aes(x=group,y=rn,))

write.csv(rn1, "~/Desktop/rn1.csv")
book2<-read.csv("~/Desktop/Book2.csv")
rn3<-rn1[1:30,]
ggplot(book2,aes(x=group,y=rn, fill=value))+ geom_tile(color = "black") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

ggplot(rn1,aes(x=group,y=rn, fill=value))+ geom_tile(color = "black") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))



