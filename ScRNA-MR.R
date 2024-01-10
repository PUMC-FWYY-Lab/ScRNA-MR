library(Seurat)
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
library(ggplot2)
library(reshape2)

setwd("~/AF_atlas/原始数据/人心房数据")
counts <- Read10X(data.dir = "AF1")
Object <- CreateSeuratObject(counts = counts,project = "AF1", min.features=200)
Object <- RenameCells(object = Object, add.cell.id = "AF1")
Object[["percent.mt"]] <- PercentageFeatureSet(object = Object, pattern = "^MT-")
dim(Object) # 33694  4155
VlnPlot(object = Object , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
Object_QC<-subset(x = Object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & 
                    nCount_RNA < 10^(mean(log10(Object$nCount_RNA))+2*sd(log10(Object$nCount_RNA))) & 
                    nCount_RNA > 10^(mean(log10(Object$nCount_RNA))-2*sd(log10(Object$nCount_RNA))) & 
                    percent.mt < 25)
dim(Object_QC) # 33694  3706
Object_QC[["condition"]] <- "AF1"

AF1<-Object_QC

#-----------------------------------------

setwd("~/AF_atlas/原始数据/人心房数据")
counts <- Read10X(data.dir = "AF2")
Object <- CreateSeuratObject(counts = counts,project = "AF2", min.features=200)
Object <- RenameCells(object = Object, add.cell.id = "AF2")
Object[["percent.mt"]] <- PercentageFeatureSet(object = Object, pattern = "^MT-")
dim(Object) # 33694   166
VlnPlot(object = Object , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
Object_QC<-subset(x = Object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & 
                    nCount_RNA < 10^(mean(log10(Object$nCount_RNA))+2*sd(log10(Object$nCount_RNA))) & 
                    nCount_RNA > 10^(mean(log10(Object$nCount_RNA))-2*sd(log10(Object$nCount_RNA))) & 
                    percent.mt < 20)
dim(Object_QC) # 33694   133
Object_QC[["condition"]] <- "AF2"

AF2<-Object_QC

#-----------------------------------------

setwd("~/AF_atlas/原始数据/人心房数据")
counts <- Read10X(data.dir = "AF3")
Object <- CreateSeuratObject(counts = counts,project = "AF3", min.features=200)
Object <- RenameCells(object = Object, add.cell.id = "AF3")
Object[["percent.mt"]] <- PercentageFeatureSet(object = Object, pattern = "^MT-")
dim(Object) # 33694  1403
VlnPlot(object = Object , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
Object_QC<-subset(x = Object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & 
                    nCount_RNA < 10^(mean(log10(Object$nCount_RNA))+2*sd(log10(Object$nCount_RNA))) & 
                    nCount_RNA > 10^(mean(log10(Object$nCount_RNA))-2*sd(log10(Object$nCount_RNA))) & 
                    percent.mt < 25)
dim(Object_QC) # 33694  1186
Object_QC[["condition"]] <- "AF3"

AF3<-Object_QC

#-----------------------------------------

setwd("~/AF_atlas/原始数据/人心房数据")
counts <- Read10X(data.dir = "AF4")
Object <- CreateSeuratObject(counts = counts,project = "AF4", min.features=200)
Object <- RenameCells(object = Object, add.cell.id = "AF4")
Object[["percent.mt"]] <- PercentageFeatureSet(object = Object, pattern = "^MT-")
dim(Object) # 33694   955
VlnPlot(object = Object , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
Object_QC<-subset(x = Object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & 
                    nCount_RNA < 10^(mean(log10(Object$nCount_RNA))+2*sd(log10(Object$nCount_RNA))) & 
                    nCount_RNA > 10^(mean(log10(Object$nCount_RNA))-2*sd(log10(Object$nCount_RNA))) & 
                    percent.mt < 25)
dim(Object_QC) # 33694   820
Object_QC[["condition"]] <- "AF4"

AF4<-Object_QC

#-----------------------------------------

setwd("~/AF_atlas/原始数据/人心房数据")
counts <- Read10X(data.dir = "AF5")
Object <- CreateSeuratObject(counts = counts,project = "AF5", min.features=200)
Object <- RenameCells(object = Object, add.cell.id = "AF5")
Object[["percent.mt"]] <- PercentageFeatureSet(object = Object, pattern = "^MT-")
dim(Object) # 33694  7674
VlnPlot(object = Object , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
Object_QC<-subset(x = Object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & 
                    nCount_RNA < 10^(mean(log10(Object$nCount_RNA))+2*sd(log10(Object$nCount_RNA))) & 
                    nCount_RNA > 10^(mean(log10(Object$nCount_RNA))-2*sd(log10(Object$nCount_RNA))) & 
                    percent.mt < 25)
dim(Object_QC) # 33694  6446
Object_QC[["condition"]] <- "AF5"

AF5<-Object_QC

#-----------------------------------------

setwd("~/AF_atlas/原始数据/人心房数据")
counts <- Read10X(data.dir = "AF6")
Object <- CreateSeuratObject(counts = counts,project = "AF6", min.features=200)
Object <- RenameCells(object = Object, add.cell.id = "AF6")
Object[["percent.mt"]] <- PercentageFeatureSet(object = Object, pattern = "^MT-")
dim(Object) # 36601  4183
VlnPlot(object = Object , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
Object_QC<-subset(x = Object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & 
                    nCount_RNA < 10^(mean(log10(Object$nCount_RNA))+2*sd(log10(Object$nCount_RNA))) & 
                    nCount_RNA > 10^(mean(log10(Object$nCount_RNA))-2*sd(log10(Object$nCount_RNA))) & 
                    percent.mt < 25)
dim(Object_QC) #  36601  3678
Object_QC[["condition"]] <- "AF6"

AF6<-Object_QC

#-----------------------------------------

setwd("~/AF_atlas/原始数据/人心房数据")
counts <- Read10X(data.dir = "AF7")
Object <- CreateSeuratObject(counts = counts,project = "AF7", min.features=200)
Object <- RenameCells(object = Object, add.cell.id = "AF7")
Object[["percent.mt"]] <- PercentageFeatureSet(object = Object, pattern = "^MT-")
dim(Object) # 36601  9138
VlnPlot(object = Object , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
Object_QC<-subset(x = Object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & 
                    nCount_RNA < 10^(mean(log10(Object$nCount_RNA))+2*sd(log10(Object$nCount_RNA))) & 
                    nCount_RNA > 10^(mean(log10(Object$nCount_RNA))-2*sd(log10(Object$nCount_RNA))) & 
                    percent.mt < 25)
dim(Object_QC) #  36601  8105
Object_QC[["condition"]] <- "AF7"

AF7<-Object_QC

#-----------------------------------------

setwd("~/AF_atlas/原始数据/人心房数据")
counts <- Read10X(data.dir = "C1")
Object <- CreateSeuratObject(counts = counts,project = "C1", min.features=200)
Object <- RenameCells(object = Object, add.cell.id = "C1")
Object[["percent.mt"]] <- PercentageFeatureSet(object = Object, pattern = "^MT-")
dim(Object) # 33694  5414
VlnPlot(object = Object , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
Object_QC<-subset(x = Object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & 
                    nCount_RNA < 10^(mean(log10(Object$nCount_RNA))+2*sd(log10(Object$nCount_RNA))) & 
                    nCount_RNA > 10^(mean(log10(Object$nCount_RNA))-2*sd(log10(Object$nCount_RNA))) & 
                    percent.mt < 25)
dim(Object_QC) #  33694  4914
Object_QC[["condition"]] <- "C1"

C1<-Object_QC

#-----------------------------------------

setwd("~/AF_atlas/原始数据/人心房数据")
counts <- Read10X(data.dir = "C2")
Object <- CreateSeuratObject(counts = counts,project = "C2", min.features=200)
Object <- RenameCells(object = Object, add.cell.id = "C2")
Object[["percent.mt"]] <- PercentageFeatureSet(object = Object, pattern = "^MT-")
dim(Object) #  33694  2752
VlnPlot(object = Object , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
Object_QC<-subset(x = Object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & 
                    nCount_RNA < 10^(mean(log10(Object$nCount_RNA))+2*sd(log10(Object$nCount_RNA))) & 
                    nCount_RNA > 10^(mean(log10(Object$nCount_RNA))-2*sd(log10(Object$nCount_RNA))) & 
                    percent.mt < 25)
dim(Object_QC) #  33694  2492
Object_QC[["condition"]] <- "C2"

C2<-Object_QC

#-----------------------------------------

setwd("~/AF_atlas/原始数据/人心房数据")
counts <- Read10X(data.dir = "C3")
Object <- CreateSeuratObject(counts = counts,project = "C3", min.features=200)
Object <- RenameCells(object = Object, add.cell.id = "C3")
Object[["percent.mt"]] <- PercentageFeatureSet(object = Object, pattern = "^MT-")
dim(Object) #  36601  3019
VlnPlot(object = Object , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
Object_QC<-subset(x = Object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & 
                    nCount_RNA < 10^(mean(log10(Object$nCount_RNA))+2*sd(log10(Object$nCount_RNA))) & 
                    nCount_RNA > 10^(mean(log10(Object$nCount_RNA))-2*sd(log10(Object$nCount_RNA))) & 
                    percent.mt < 25)
dim(Object_QC) #  36601  2695
Object_QC[["condition"]] <- "C3"

C3<-Object_QC

#-----------------------------------------

setwd("~/AF_atlas/原始数据/人心房数据")
counts <- Read10X(data.dir = "C4")
Object <- CreateSeuratObject(counts = counts,project = "C4", min.features=200)
Object <- RenameCells(object = Object, add.cell.id = "C4")
Object[["percent.mt"]] <- PercentageFeatureSet(object = Object, pattern = "^MT-")
dim(Object) #  36601  2065
VlnPlot(object = Object , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
Object_QC<-subset(x = Object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & 
                    nCount_RNA < 10^(mean(log10(Object$nCount_RNA))+2*sd(log10(Object$nCount_RNA))) & 
                    nCount_RNA > 10^(mean(log10(Object$nCount_RNA))-2*sd(log10(Object$nCount_RNA))) & 
                    percent.mt < 25)
dim(Object_QC) #  36601  1862
Object_QC[["condition"]] <- "C4"

C4<-Object_QC

#-----------------------------------------

setwd("~/AF_atlas/原始数据/人心房数据")
counts <- Read10X(data.dir = "C5")
Object <- CreateSeuratObject(counts = counts,project = "C5", min.features=200)
Object <- RenameCells(object = Object, add.cell.id = "C5")
Object[["percent.mt"]] <- PercentageFeatureSet(object = Object, pattern = "^MT-")
dim(Object) #  36601  3829
VlnPlot(object = Object , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
Object_QC<-subset(x = Object, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & 
                    nCount_RNA < 10^(mean(log10(Object$nCount_RNA))+2*sd(log10(Object$nCount_RNA))) & 
                    nCount_RNA > 10^(mean(log10(Object$nCount_RNA))-2*sd(log10(Object$nCount_RNA))) & 
                    percent.mt < 25)
dim(Object_QC) #  36601  1862
Object_QC[["condition"]] <- "C5"

C5<-Object_QC


#---------------------------------------------------------------
AF_combine<-merge(AF1,c(AF2,AF3,AF4,AF5,AF6,AF7,C1,C2,C3,C4,C5))
dim(AF_combine)  #   48361 39466

Idents(AF_combine) <- AF_combine$condition
patient_ids <- unique(Idents(AF_combine))

# 定义质控函数
single_cell_qc <- function(patient_id, kidney_map) {
  patient_data <- subset(kidney_map, idents = patient_id)
  
  # 质控
  patient_data_QC <- patient_data
  
  # 去除MT和ERCC
  mt_ercc_genes <- grep(pattern = "^(MT-|ERCC)", rownames(patient_data_QC), value = TRUE)
  patient_data_QC_no_ERCC_MT <- patient_data_QC[!rownames(patient_data_QC) %in% mt_ercc_genes, ]
  
  # 去除双重细胞
  patient_data_SCT <- SCTransform(patient_data_QC_no_ERCC_MT)
  patient_data_SCT <- RunPCA(patient_data_SCT)
  patient_data_SCT <- RunUMAP(patient_data_SCT, dims = 1:10)
  
  sweep.res.list <- paramSweep_v3(patient_data_SCT, PCs = 1:10, sct = T)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  nExp_poi <- round(0.075 * nrow(patient_data_SCT@meta.data))
  
  patient_data_filtered <- doubletFinder_v3(patient_data_SCT, PCs = 1:10, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  
  # 获取DF.classifications的列名
  df_classification_colname <- grep("^DF.classifications", colnames(patient_data_filtered@meta.data), value = TRUE)
  patient_data_filtered <- patient_data_filtered[, patient_data_filtered@meta.data[, df_classification_colname] == "Singlet"]
  
  # 将处理后的数据存储在一个列表中并返回
  processed_data <- list(
    QC_data = patient_data_QC,
    QC_no_ERCC_MT_data = patient_data_QC_no_ERCC_MT,
    filtered_data = patient_data_filtered
  )
  
  return(processed_data)
}


# 在循环开始之前，创建空列表以记录每个处理步骤的细胞数目
cell_counts_before_qc <- c()
cell_counts_after_qc <- c()
cell_counts_no_MT_ERCC <- c()
cell_counts_no_doublets <- c()

# 创建一个空的 Seurat 对象用于存储合并后的数据
AF_QC <- NULL

for (patient_id in patient_ids) {
  patient_data <- subset(AF_combine, idents = patient_id)
  
  cell_counts_before_qc[[patient_id]] <- dim(patient_data)[2]
  
  # 获取各处理步骤的细胞数目
  patient_data_filtered <- single_cell_qc(patient_id, AF_combine)
  
  cell_counts_after_qc[[patient_id]] <- dim(patient_data_filtered$QC_data)[2]
  cell_counts_no_MT_ERCC[[patient_id]] <- dim(patient_data_filtered$QC_no_ERCC_MT_data)[2]
  cell_counts_no_doublets[[patient_id]] <- dim(patient_data_filtered$filtered_data)[2]
  
  # 合并处理后的患者数据
  if (is.null(AF_QC)) {
    AF_QC <- patient_data_filtered$filtered_data
  } else {
    AF_QC <- merge(AF_QC, patient_data_filtered$filtered_data)
  }
}

# 数据存储
saveRDS(AF_QC, file = "~/AF_atlas/原始数据/人心房数据/AF_QC.RDS")

dim(AF_QC)

# 循环结束后，将记录的细胞数目整合成一个数据框
patient_ids <- names(cell_counts_before_qc)
cell_counts_df <- data.frame(
  patient_id = patient_ids,
  before_qc = unlist(cell_counts_before_qc),
  after_qc = unlist(cell_counts_after_qc),
  no_MT_ERCC = unlist(cell_counts_no_MT_ERCC),
  no_doublets = unlist(cell_counts_no_doublets)
)

# 使用 ggplot2 创建堆叠条形图以展示每个处理步骤的细胞数目变化
cell_counts_long <- melt(cell_counts_df, id.vars = "patient_id", variable.name = "step", value.name = "cell_count")
ggplot(cell_counts_long, aes(x = patient_id, y = cell_count, fill = step)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(x = "Patient ID", y = "Cell count", title = "Cell count changes per patient after each step")

table(Idents(AF_QC))

AF_QC<-subset(AF_QC,idents = c("AF1","AF3","AF4","AF5","AF6","AF7",
                               "C1","C2","C3","C4","C5"))

dim(AF_QC) # 24101 36382
table(Idents(AF_QC))

saveRDS(AF_QC, file = "~/AF_atlas/原始数据/人心房数据/AF_QC.RDS")


# 定义变量
input_seurat_object <- AF_QC
split_attr <- "condition"

# 使用 SplitObject 函数分割数据
Object.list <- SplitObject(input_seurat_object, split.by = split_attr)

# 在 SplitObject 后，合并细胞数量低于500的患者
cell_threshold <- 500
small_patients <- names(Object.list)[sapply(Object.list, function(x) ncol(x) < cell_threshold)]

if (length(small_patients) > 0) {
  if (length(small_patients) == 1) {
    # 找到细胞数量最接近500的患者
    closest_patient_id <- names(Object.list)[which.min(abs(sapply(Object.list, function(x) ncol(x)) - cell_threshold))]
    cat("Merging patient ", small_patients[1], " with ", closest_patient_id, "\n")
    combined_patients <- merge(x = Object.list[[small_patients[1]]], y = Object.list[[closest_patient_id]])
    
    # 使用"_"连接患者编号作为新的患者ID
    combined_patient_id <- paste(c(small_patients[1], closest_patient_id), collapse = "_")
    cat("Assigning new patient ID: ", combined_patient_id, "\n")
    combined_patients[[split_attr]] <- factor(combined_patients[[split_attr]], levels = unique(combined_patients[[split_attr]]), labels = combined_patient_id)
    
    # 从 Object.list 中移除已合并的患者
    cat("Removing merged patients from Object.list...\n")
    Object.list <- Object.list[!(names(Object.list) %in% c(small_patients[1], closest_patient_id))]
    
    # 将合并后的患者添加到 Object.list
    cat("Adding combined patient to Object.list...\n")
    Object.list[[combined_patient_id]] <- combined_patients
    
  } else {
    # 确保患者存在于 Object.list 中
    small_patients <- small_patients[small_patients %in% names(Object.list)]
    
    # 合并细胞数量较少的患者
    combined_small_patients <- Object.list[[small_patients[1]]]
    for (patient_id in small_patients[2:length(small_patients)]) {
      cat("Merging patient ", patient_id, "\n")
      combined_small_patients <- merge(x = combined_small_patients, y = Object.list[[patient_id]])
    }
    
    # 使用"_"连接患者编号作为新的患者ID
    combined_patient_id <- paste(small_patients, collapse = "_")
    cat("Assigning new patient ID: ", combined_patient_id, "\n")
    combined_small_patients[[split_attr]] <- factor(combined_small_patients[[split_attr]], levels = unique(combined_small_patients[[split_attr]]), labels = combined_patient_id)
    
    # 从 Object.list 中移除已合并的患者
    cat("Removing small patients from Object.list...\n")
    Object.list <- Object.list[!(names(Object.list) %in% small_patients)]
    
    # 将合并后的患者添加到 Object.list
    cat("Adding combined patient to Object.list...\n")
    Object.list[[combined_patient_id]] <- combined_small_patients
  }
}


for (patient_id in names(Object.list)) {
  cat(patient_id, " - Number of cells: ", ncol(Object.list[[patient_id]]), "\n")
}

# 获取可用基因
available_genes <- rownames(input_seurat_object)

available_cc_genes <- list(
  s.genes = intersect(cc.genes$s.genes, available_genes),
  g2m.genes = intersect(cc.genes$g2m.genes, available_genes)
)

process_patient_data <- function(patient_data) {
  # SCTransform
  patient_data <- SCTransform(patient_data, verbose = T)
  cat("CellCycleScoring ", patient_data$condition[1], "\n")
  
  # 尝试使用24个bins进行CellCycleScoring，如有错误则自动减少bins
  num_bins <- 24
  success <- FALSE
  while (num_bins > 0) {
    # 尝试 CellCycleScoring
    success <- tryCatch(
      {
        patient_data <- CellCycleScoring(patient_data, s.features = available_cc_genes$s.genes, g2m.features = available_cc_genes$g2m.genes, set.ident = TRUE, bins = num_bins)
        TRUE
      },
      error = function(e) {
        cat("Error with", num_bins, "bins. Retrying with", num_bins - 1, "bins.\n")
        FALSE
      }
    )
    if (success) break
    num_bins <- num_bins - 1
  }
  
  # 如果成功进行细胞周期校正，则执行以下操作
  if (success) {
    # 计算细胞周期差异
    patient_data$CC.Difference <- patient_data$S.Score - patient_data$G2M.Score
    
    # 再次执行 SCTransform 并回归 CC.Difference
    patient_data <- SCTransform(patient_data, vars.to.regress = "CC.Difference", verbose = T)
  }
  
  return(patient_data)
}

# 处理每个患者的数据
for (patient_id in names(Object.list)) {
  Object.list[[patient_id]] <- process_patient_data(Object.list[[patient_id]])
}

# 后续分析
features <- SelectIntegrationFeatures(object.list = Object.list, nfeatures = 3000)
Object.list <- PrepSCTIntegration(object.list = Object.list, anchor.features = features)
Object.anchors <- FindIntegrationAnchors(object.list = Object.list, normalization.method = 'SCT', anchor.features = features)

saveRDS(Object.anchors, file = "~/AF_atlas/原始数据/人心房数据/Object.anchors.RDS")

Object.combined.sct <- IntegrateData(anchorset = Object.anchors, normalization.method = 'SCT')

saveRDS(Object.combined.sct, file = "~/AF_atlas/原始数据/人心房数据/Object.combined.sct.RDS")

Object.combined.sct <- RunPCA(Object.combined.sct, verbose = FALSE)
Object.combined.sct <- RunUMAP(Object.combined.sct, reduction = "pca", dims = 1:30)
Object.combined.sct <- FindNeighbors(Object.combined.sct, dims = 1:30)
Object.combined.sct <- FindClusters(Object.combined.sct, resolution = 0.5)

DefaultAssay(Object.combined.sct) <- "SCT"
Object.combined.sct <- PrepSCTFindMarkers(object = Object.combined.sct, verbose = T)

saveRDS(Object.combined.sct, file = "~/AF_atlas/原始数据/人心房数据/Object.combined.sct.RDS")

DimPlot(Object.combined.sct,label = T)
FeaturePlot(Object.combined.sct,features = c("CDH5","PECAM1","TEK")) # EC
FeaturePlot(Object.combined.sct,features = c("CD68","CD14","CSF1R","ITGAM","ITGAX")) # MP
VlnPlot(Object.combined.sct,features = c("CD68","CD14","CSF1R","ITGAM","ITGAX"))
FeaturePlot(Object.combined.sct,features = c("ACTA2","CNN1","MYH11","CSPG4","RGS5"))# SMC
FeaturePlot(Object.combined.sct,features = c("COL1A1","COL1A2","COL3A1","DCN","PDGFRA"))# FB
FeaturePlot(Object.combined.sct,features = c("CD1C","CLEC10A","FCER1A")) # DC
FeaturePlot(Object.combined.sct,features = c("FCGR3B","PTGS2","S100A8")) # Neutrophil
VlnPlot(Object.combined.sct,features = c("FCGR3B","PTGS2","S100A8"))
FeaturePlot(Object.combined.sct,features = c("CD3E","CD8A","GZMB")) # T cell
FeaturePlot(Object.combined.sct,features = c("CD79A","CD79B","MS4A1")) # B cell
FeaturePlot(Object.combined.sct,features = c("KIT","TPSAB1","TPSB2")) # Mast cell

pal <- c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", 
         "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", 
         "#E8601C", "#DC050C")

pal_30<-c("#EBA6B0","#6085BF","#8E82B0","#E4756A","#E66A6C",
          "#DCB9CE","#8BC078","#E56B7F","#F4B151","#F5BC53",
          "#65C0E8","#35AC9A","#EB9671","#31A5D9","#DE5758",
          "#C597B9","#D64448","#C2197E","#6AB8B1","#D579A5",
          "#C8B3C9","#C28339","#BEBFD7","#31A5D9","#7580B5",
          "#F49F2E","#7064A0","#87C9D3","#3A6CAF","#375BA4")

Object.combined.sct<-RenameIdents(Object.combined.sct,
                                   "10"="EC",
                                   "17"="EC",
                                   "13"="EC",
                                   "19"="EC",
                                   "9"="DC",
                                   "12"="SMC",
                                   "14"="SMC",
                                   "23"="SMC",
                                   "8"="FB",
                                   "18"="B_cell",
                                   "2"="T_cell",
                                   "3"="T_cell",
                                   "4"="T_cell",
                                   "5"="T_cell",
                                   "0"="T_cell",
                                   "6"="T_cell",
                                   "15"="T_cell",
                                   "16"="T_cell",
                                   "20"="T_cell",
                                   "21"="T_cell",
                                   "1"="MP",
                                   "11"="MP",
                                   "23"="MP",
                                   "7"="Neutrophil")

table(Idents(Object.combined.sct))

Object.combined.sct$cell_type<-Idents(Object.combined.sct)

list<-levels(Object.combined.sct)
cluster<- setdiff(list, "22")
dim(Object.combined.sct) # 23740 36382
Object.combined.sct<- subset(Object.combined.sct, idents = cluster)
dim(Object.combined.sct) # 23740 36327
table(Object.combined.sct)

AF_renamed<-Object.combined.sct
DimPlot(AF_renamed,cols = pal)
table(AF_renamed$cell_type)

table(AF_renamed$condition)

AF_renamed$condition <- sapply(AF_renamed$condition, function(x) {
  ifelse(grepl("^AF", x), "AF", 
         ifelse(grepl("^C", x), "Ctrl", x))
})

# 检查修改后的表
table(AF_renamed$condition)

saveRDS(AF_renamed, file = "~/AF_atlas/Data/单细胞分析/AF_renamed.RDS")

#---------------------------------- 数据转化为h5ad
library(SeuratDisk)

DefaultAssay(AF_renamed)

SaveH5Seurat(AF_renamed, filename = "AF_renamed.h5Seurat",overwrite = T)

Convert("AF_renamed.h5Seurat", dest = "h5ad",assay = "RNA",overwrite = T,verbose = T)

#--------------------------------------

matrix<-AF_renamed@assays$RNA@counts
dim(matrix)
meta<-AF_renamed@meta.data
dim(meta)

library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

genes <- rownames(matrix) 
gene.mapping <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                      filters = 'hgnc_symbol',
                      values = genes,
                      mart = ensembl)

# 创建一个基因名称到Ensembl ID的映射向量，仅包括映射表中存在的基因
gene.mapping.vector <- setNames(gene.mapping$ensembl_gene_id, gene.mapping$hgnc_symbol)

# 确定哪些基因名称已经在映射表中找到对应的Ensembl ID
matched_genes <- genes %in% gene.mapping$hgnc_symbol
# 仅保留那些能够匹配到Ensembl ID的基因
matrix_matched <- matrix[matched_genes, ]

# 更新行名为Ensembl ID，仅用已匹配的基因的Ensembl ID
rownames(matrix_matched) <- gene.mapping.vector[rownames(matrix_matched)]


head(matrix_matched)

AF_renamed_ensembl<-CreateSeuratObject(counts = matrix_matched,meta.data = meta)
dim(AF_renamed_ensembl)
table(Idents(AF_renamed_ensembl))
table(AF_renamed_ensembl$cell_type)

saveRDS(AF_renamed_ensembl, file = "~/AF_atlas/Data/单细胞分析/AF_renamed_ensembl.RDS")

SaveH5Seurat(AF_renamed_ensembl, filename = "AF_renamed_ensembl.h5Seurat",overwrite = T)

Convert("AF_renamed_ensembl.h5Seurat", dest = "h5ad",assay = "RNA",overwrite = T,verbose = T)


#-------------------------------------- 后续分析房颤数据-----------------------------------
DimPlot(AF_renamed,cols = pal,label = T)

table(AF_renamed$condition,AF_renamed$cell_type)
table(AF_renamed$cell_type)
VlnPlot(AF_renamed,features = "FAM13B")

#-------------------------------------- 孟德尔随机化---------------------------------------
setwd("~/AF_atlas/Data/in silico genes")

#af_data <- read.csv("AF.csv")
library(dplyr)
#sig_genes<-af_data %>% filter(Sig==1)
#candidate_genes<-unique(sig_genes$Ensembl_ID) # 167
setwd("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据")
candidate_genes_有意义靶点基因补充基因名 <- readRDS("~/AF_atlas/Data/in silico genes/candidate_genes_有意义靶点基因补充基因名.RDS")
candidate_genes<-rownames(candidate_genes_有意义靶点基因补充基因名)
writeLines(candidate_genes, "All_candidate_genes.txt")

# zcat Heart_Atrial_Appendage.tsv.gz | head -n 1 > All_candidate_genes.tsv
# zcat Heart_Atrial_Appendage.tsv.gz | awk 'NR==FNR{genes[$1]; next} FNR!=1 && $7 in genes' All_candidate_genes.txt - >> All_candidate_genes.tsv

# 本地循环clump
library(TwoSampleMR)
library(MRInstruments)
#devtools::install_github("explodecomputer/plinkbinr")
library(plinkbinr)
#get_plink_exe() #######此行代码运行之后会显示--[1]"plink exe的路径"
#devtools::install_github("mrcieu/ieugwasr")
library(ieugwasr)

# 心脏组织eQTL
#zcat Heart_Atrial_Appendage.tsv.gz | head -n 2

heart_exp_dat <- read_exposure_data(
  filename = "~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/All_candidate_genes.tsv",
  sep = "\t",
  phenotype_col = "gene_id",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "maf",
  pval_col = "pvalue"
)

head(heart_exp_dat)

heart_exp<- heart_exp_dat %>% filter(pval.exposure <5e-8)
length(unique(heart_exp$exposure)) # 101

saveRDS(heart_exp, file = "~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/heart_exp_for_all.RDS")

#---------------------------------------------- 心脏本地clump---------------------------------------------- 
heart_exp <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/heart_exp_for_all.RDS")

# 创建一个空的数据框，用于存储 heart_eqtl_local_clumped
heart_eqtl_local_clumped <- data.frame()

# 获取所有不同的 id.exposure 值
unique_ids <- unique(heart_exp$id.exposure)

# 循环遍历不同的 id.exposure 值
for (id in unique_ids) {
  # 从 lung_exp_dat 中获取当前 id.exposure 的子集
  subset_data <- filter(heart_exp, id.exposure == id)
  
  # 使用 tryCatch 处理可能的错误
  clump_result <- tryCatch({
    ld_clump(
      clump_kb = 10000,
      clump_r2 = 0.001,
      clump_p = 0.99,
      pop = "EUR",
      dat = dplyr::tibble(rsid = subset_data$SNP, pval = subset_data$pval.exposure, id = subset_data$id.exposure),
      plink_bin = "/usr/local/lib/R/site-library/plinkbinr/bin/plink_Linux",
      bfile = "/home/wangxihe/Lung_Atlas/单细胞孟德尔/文章部分二：孟德尔随机化/GWAS原始数据/1kg.v3/EUR"
    )
  }, error = function(e) {
    message("Error in ld_clump for id.exposure ", id, ":", conditionMessage(e))
    NULL
  })
  
  # 如果 clump_result 不为空，将其添加到数据框中
  if (!is.null(clump_result)) {
    clump_result$id.exposure <- id  # 添加 id.exposure 列
    
    # 创建 heart_eqtl_local_clumped 子集并将其添加到 heart_eqtl_local_clumped 中
    subset_local_clumped <- subset(heart_exp, SNP %in% clump_result$rsid)
    heart_eqtl_local_clumped <- rbind(heart_eqtl_local_clumped, subset_local_clumped)
  }
}

saveRDS(heart_eqtl_local_clumped, file = "~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/heart_eqtl_local_clumped_for_all.RDS")

#-------------------------------------------------------心脏eQTL孟德尔-----------------------------------------
heart_eqtl_local_clumped_for_all <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/heart_eqtl_local_clumped_for_all.RDS")
filtered_exp_dat<-heart_eqtl_local_clumped_for_all

setwd("~/AF_atlas/Data/孟德尔随机化分析数据/GWAS数据")
filename="GCST90038611"
data <- read.table(paste0(filename, ".tsv"), header = TRUE, sep = "\t")

head(data)
# 存储为 CSV 文件
write.csv(data, file = paste0(filename, ".csv"), row.names = FALSE)

filename_CSV<-paste0(filename, ".csv")

outcome_dat <- read_outcome_data(
  snps = filtered_exp_dat$SNP,
  filename = filename_CSV,
  sep = ",",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value"
)

harmonise_data <- harmonise_data(exposure_dat = filtered_exp_dat, outcome_dat = outcome_dat)
res <- mr(harmonise_data)
num=length(unique(res$exposure))
num
adj_P=0.05/num

sig_targets<-res %>% filter(pval<adj_P)
length(unique(sig_targets$exposure))

#-------------------------------------------------------心脏eQTL 共定位-----------------------------------------
library(coloc)
data_eQTL <- read.table("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/All_candidate_genes.tsv", header = TRUE, sep = "\t")
length(unique(data_eQTL$gene_id))
head(data_eQTL)

sig_genes<-sig_targets$exposure

coloc_combine_results <- data.frame()

folder_path2 <- "~/AF_atlas/Data/孟德尔随机化分析数据/GWAS数据/"
filename="GCST90038611"
data_outcome<- read.csv(file.path(folder_path2, paste0(filename, ".csv")))
head(data_outcome)

for (gene in sig_genes) {
  # 从 data_eQTL 中提取数据
  sub_data_eQTL <- data_eQTL[data_eQTL$gene_id == gene, ]
  data1 <- data.frame(
    SNP = sub_data_eQTL$rsid,
    CHR = sub_data_eQTL$chromosome,
    BP = sub_data_eQTL$position,
    A1 = sub_data_eQTL$alt,
    A2 = sub_data_eQTL$ref,
    BETA = sub_data_eQTL$beta,
    SE = sub_data_eQTL$se,
    MAF = sub_data_eQTL$maf,
    N = 372,
    pvalue = sub_data_eQTL$pvalue
  )
  
  # 从 data_outcome 中提取数据
  
  data2<-data.frame(SNP=data_outcome$variant_id,CHR=data_outcome$chromosome,
                    BP=data_outcome$base_pair_location,A1=data_outcome$effect_allele,A2=data_outcome$other_allele,
                    BETA=data_outcome$beta,SE=data_outcome$standard_error,pvalue=data_outcome$p_value)
  
  # 合并数据并去重
  data = merge(data1, data2, by = "SNP")
  data = data[!duplicated(data$SNP),]
  
  # 对齐效应等位基因
  data = data %>% filter((A1.x == A1.y & A2.x == A2.y) | (A1.x == A2.y & A2.x == A1.y))
  data = data %>% mutate(BETA.y = ifelse(A1.x == A1.y, BETA.y, -BETA.y))
  
  # 计算方差VAR
  data$VAR.x = data$SE.x^2
  data$VAR.y = data$SE.y^2
  data = data[data$VAR.x != 0 & data$VAR.y != 0,]
  
  # 拆分整理
  data1 = data[, c("BETA.x", "VAR.x", "SNP", "MAF", "N")]
  data2 = data[, c("BETA.y", "VAR.y", "SNP")]
  colnames(data1) = c("beta", "varbeta", "snp", "MAF", "N")
  colnames(data2) = c("beta", "varbeta", "snp")
  data1 = as.list(data1)
  data2 = as.list(data2)
  # 声明表型类型,二分类表型"cc",连续型表型"quant" 
  data1$type = "quant"
  data2$type = "cc"
  
  #运行 coloc.abf 函数
  coloc_res = coloc.abf(data2, data1)
  
}


for (gene in sig_genes) {
  # 从 data_eQTL 中提取数据
  sub_data_eQTL <- data_eQTL[data_eQTL$gene_id == gene, ]
  data1 <- data.frame(
    SNP = sub_data_eQTL$rsid,
    CHR = sub_data_eQTL$chromosome,
    BP = sub_data_eQTL$position,
    A1 = sub_data_eQTL$ref,
    A2 = sub_data_eQTL$alt,
    BETA = sub_data_eQTL$beta,
    SE = sub_data_eQTL$se,
    MAF = sub_data_eQTL$maf,
    N = 515,
    pvalue = sub_data_eQTL$pvalue
  )
  
  # 从 data_outcome 中提取数据
  
  data2<-data.frame(SNP=data_outcome$rsids,CHR=data_outcome$X.chrom,
                    BP=data_outcome$pos,A1=data_outcome$ref,A2=data_outcome$alt,
                    BETA=data_outcome$beta,SE=data_outcome$sebeta,pvalue=data_outcome$pval)
  
  # 合并数据并去重
  data = merge(data1, data2, by = "SNP")
  data = data[!duplicated(data$SNP),]
  
  # 对齐效应等位基因
  data = data %>% filter((A1.x == A1.y & A2.x == A2.y) | (A1.x == A2.y & A2.x == A1.y))
  data = data %>% mutate(BETA.y = ifelse(A1.x == A1.y, BETA.y, -BETA.y))
  
  # 计算方差VAR
  data$VAR.x = data$SE.x^2
  data$VAR.y = data$SE.y^2
  data = data[data$VAR.x != 0 & data$VAR.y != 0,]
  
  # 拆分整理
  data1 = data[, c("BETA.x", "VAR.x", "SNP", "MAF", "N")]
  data2 = data[, c("BETA.y", "VAR.y", "SNP")]
  colnames(data1) = c("beta", "varbeta", "snp", "MAF", "N")
  colnames(data2) = c("beta", "varbeta", "snp")
  data1 = as.list(data1)
  data2 = as.list(data2)
  # 声明表型类型,二分类表型"cc",连续型表型"quant" 
  data1$type = "quant"
  data2$type = "cc"
  
  #运行 coloc.abf 函数
  coloc_res = coloc.abf(data2, data1)
  PP.H4.abf <- coloc_res$summary["PP.H4.abf"]
  
  # 从 candidate_genes 中提取方向信息
  candidate_genes[gene,]
  res_sig[res_sig$exposure==gene,]
  snp_gene<-subset(filtered_exp_dat,exposure==gene)
  
  coloc_results<-data.frame(SNP=snp_gene$SNP,
                            exposure=res_sig[res_sig$exposure==gene,]$exposure,
                            nsnp=res_sig[res_sig$exposure==gene,]$nsnp,
                            b=res_sig[res_sig$exposure==gene,]$b,
                            se=res_sig[res_sig$exposure==gene,]$se,
                            pval=res_sig[res_sig$exposure==gene,]$pval,
                            PP.H4.abf=PP.H4.abf)
  
  coloc_combine_results <- rbind(coloc_combine_results, coloc_results)
  
}

# 打印最终结果
candidate_genes_direction<-candidate_genes[rownames(candidate_genes) %in% coloc_combine_results$exposure,]
candidate_genes_direction$exposure<-rownames(candidate_genes_direction)

merged_data <- data.frame(left_join(coloc_combine_results,candidate_genes_direction, by = "exposure"))
merged_data$num_p_adj<-num
merged_data$adj_p_level<-0.05/num
colnames(merged_data)
part1<-merged_data[,c("exposure","method","nsnp","b","se","pval","PP.H4.abf","SNP","num_p_adj","adj_p_level")]
part2<-merged_data[,22:72]
merged_data<-cbind(part1,part2)

assign(paste0(filename, "_coloc_merged_data"), merged_data)
folder_path <- "~/Lung_Atlas/单细胞孟德尔/文章部分二：孟德尔随机化/肺eQTL/"
rds_file_path3 <- file.path(folder_path, filename, paste0(filename, "_coloc_merged_data.RDS"))

saveRDS(get(paste0(filename, "_coloc_merged_data")), file = rds_file_path3)

#---------------------------------------------- 血eQTL---------------------------------------------- 
Blood_eqtl_local_clumped <- readRDS("~/Lung_Atlas/单细胞孟德尔/文章部分二：孟德尔随机化/全血/Blood_eqtl_local_clumped.RDS")
head(Blood_eqtl_local_clumped)

setwd("~/AF_atlas/Data/in silico genes")

af_data <- read.csv("AF.csv")
library(dplyr)
sig_genes<-af_data %>% filter(Sig==1)
length(unique(sig_genes$Gene)) # 167

candidate_genes<-intersect(sig_genes$Ensembl_ID,Blood_eqtl_local_clumped$exposure)

filtered_exp_dat<-Blood_eqtl_local_clumped %>% filter(exposure %in% candidate_genes)

setwd("~/AF_atlas/Data/孟德尔随机化分析数据/GWAS数据")
filename="GCST006414"
data <- read.table(paste0(filename, ".tsv"), header = TRUE, sep = "\t")

head(data,20)
# 存储为 CSV 文件
write.csv(data, file = paste0(filename, ".csv"), row.names = FALSE)

filename_CSV<-paste0(filename, ".csv")

outcome_dat <- read_outcome_data(
  snps = filtered_exp_dat$SNP,
  filename = filename_CSV,
  sep = ",",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value"
)

harmonise_data <- harmonise_data(exposure_dat = filtered_exp_dat, outcome_dat = outcome_dat)
res <- mr(harmonise_data)
length(unique(res$exposure))
adj_P=0.05/50

sig_targets<-res %>% filter(pval<adj_P)

length(unique(sig_targets$exposure))














