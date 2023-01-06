library(dplyr)
library(Seurat)
library(patchwork)

#Feature = gene/row   Cell = column

data <- Read10X(data.dir = "C:/Users/shriy/Documents/Research/filtered_feature_bc_matrix")

sobj<- CreateSeuratObject(counts = data, project = "GCA43_R", min.cells = 3, min.features = 200) 

sobj

sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")

#VlnPlot(sobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#plot1 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(sobj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2



#filter out cell with genes over 2500 and less than 200 and mitochondrial counts over 5%
sobj <- subset(sobj, subset = nFeature_RNA < 10000)


sobj

sobj <- NormalizeData(sobj, normalization.method = "LogNormalize", scale.factor = 10000)

sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 9000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sobj), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sobj)
plot2 <- LabelPoints(plot = plot1, points = top10)
plot2


all.genes <- rownames(sobj)
sobj <- ScaleData(sobj, features = all.genes)

sobj <- RunPCA(sobj, features = VariableFeatures(object = sobj))

#DimHeatmap(sobj, dims = 1:18, cells = 500, balanced = TRUE)

sobj <- JackStraw(sobj, num.replicate = 100)

sobj1 <- ScoreJackStraw(sobj, dims = 1:20)

JackStrawPlot(sobj1, dims = 1:20)


ElbowPlot(sobj)

sobj <- FindNeighbors(sobj, dims = 1:15) 
sobj <- FindClusters(sobj, resolution = 0.5)


sobj <- RunUMAP(sobj, dims = 1:15)
DimPlot(sobj, reduction = "umap")

FeaturePlot(sobj, features = c("BEST4", "CA4", "OTOP2", "CA7"))#BEST4+ enterocyte

FeaturePlot(sobj, features = c("CHGA", "CHGB", "CPE", "NEUROD1", "PYY","NTS"))#enteroendocrine cells

FeaturePlot(sobj, features = c("ZG16", "CLCA1", "FFAR4", "TFF3", "SPINK4","MUC2"))#goblet

FeaturePlot(sobj, features = c("LYZ", "CA7", "SPIB", "CA4", "FKBP1A","PLA2G2A","DEFA6","PRSS2","DEFA5", "REG3A"))#Paneth

FeaturePlot(sobj, features = c("MKI67", "PCNA", "TOP2A", "CCNA2", "MCM5"))#TA

FeaturePlot(sobj, features = c("LGR5", "RGMB", "SMOC2", "ASCL2","OLFM4"))#stem

#-------------------------------------------------------------------------------
#Immune Markers
#-------------------------------------------------------------------------------
FeaturePlot(sobj, features = c("PTPRC"))#Immune cells

FeaturePlot(sobj, features = c("VPREB3","VPREB1", "CD79A","CD79B"))#pre-Bcells

FeaturePlot(sobj, features = c("SDC1","JCHAIN","IGHA1","IGHM","MZB1","IGHG1"))#Plasma 

FeaturePlot(sobj, features = c("CD3D"))#T cells

FeaturePlot(sobj, features = c("KLRzD1", "NCAM1", "GNLY", "NKG7","GZMA","GZMK"))#NK

FeaturePlot(sobj, features = c("CD11C", "CD163", "CD68", "C1QA", "C1QB","SEPP1","C1QC"))#macrophages

FeaturePlot(sobj, features = c("CXCL8", "CCL3","IL1B","S100A8"))#Inf.macros

FeaturePlot(sobj, features = c("HBB", "HBA1","HBA2"))#erythrocytes

#-------------------------------------------------------------------------------
#stromal Markers
#-------------------------------------------------------------------------------

FeaturePlot(sobj, features = c("ACTA2","COL1A2","COL3A1","ADAMDEC1",""))#stromal




sobj.markers <- FindAllMarkers(sobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sobj.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

write.csv(sobj.markers,"C:\\Users\\shriy\\Documents\\Research\\Markers.csv", row.names = TRUE)
