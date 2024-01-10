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


#----------------------------------------------------------------------------------------
library(Seurat)
#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
library(ggplot2)
library(reshape2)


#-------------------------------------- 心脏孟德尔随机化-----------------------------------
heart_exp_dat <- read_exposure_data(
  filename = "~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/Heart_Atrial_Appendage.tsv",
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
length(unique(heart_exp$exposure)) # 27

saveRDS(heart_exp, file = "~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/heart_exp_for_all_heart.RDS")

#---------------------------------------------- 心脏本地clump---------------------------------------------- 
heart_exp <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/heart_exp_for_all_heart.RDS")
length(unique(heart_exp$exposure))

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

length(unique(heart_eqtl_local_clumped$exposure))

saveRDS(heart_eqtl_local_clumped, file = "~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/heart_eqtl_local_clumped_for_all_heart.RDS")

#-------------------------------------------------------心脏eQTL孟德尔-----------------------------------------
heart_eqtl_local_clumped_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/heart_eqtl_local_clumped_for_all_heart.RDS")
filtered_exp_dat<-heart_eqtl_local_clumped_for_all_heart

setwd("~/AF_atlas/Data/孟德尔随机化分析数据/GWAS数据")
filename="GCST006414"
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
saveRDS(res, file = "~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/heart_eqtl_res_for_all_heart.RDS")
heart_eqtl_res_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/heart_eqtl_res_for_all_heart.RDS")
num=length(unique(heart_eqtl_res_for_all_heart $exposure))
num
adj_P=0.05/num

sig_targets<-heart_eqtl_res_for_all_heart  %>% filter(pval<adj_P)
length(unique(sig_targets$exposure)) #132

#-------------------------------------------------------心脏eQTL 共定位-----------------------------------------
candidate_genes<-unique(sig_targets$exposure) # 132

setwd("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据")
writeLines(candidate_genes, "Candidate_genes_for_all_heart_eQTLs.txt")

# zcat Heart_Atrial_Appendage.tsv.gz | head -n 1 > Candidate_genes_for_all_heart_eQTLs.tsv
# zcat Heart_Atrial_Appendage.tsv.gz | awk 'NR==FNR{genes[$1]; next} FNR!=1 && $7 in genes' Candidate_genes_for_all_heart_eQTLs.txt - >> Candidate_genes_for_all_heart_eQTLs.tsv

library(coloc)
data_eQTL <- read.table("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/Candidate_genes_for_all_heart_eQTLs.tsv", header = TRUE, sep = "\t")
length(unique(data_eQTL$gene_id))
head(data_eQTL)

{
  # 假设data_eQTL是您的数据框，rsid是其中的一个列
  # 去除空白
  data_eQTL$rsid <- trimws(data_eQTL$rsid)
  
  # 检查格式 - 假设所有的variant_id都应以"rs"开头
  valid_format <- grepl("^rs", data_eQTL$rsid)
  
  # 排除NA值或空字符串
  valid_entries <- !is.na(data_eQTL$rsid) & data_eQTL$rsid != ""
  
  # 合并所有的检查
  valid_variant_ids <- valid_format & valid_entries
  
  # 提取有效的variant_id
  valid_data_eQTL <- data_eQTL[valid_variant_ids, ]
  
  # 检查无效的variant_id，这可以帮助您识别出问题所在
  invalid_data_eQTL <- data_eQTL[!valid_variant_ids, ]
  
}

setwd("~/AF_atlas/Data/孟德尔随机化分析数据/GWAS数据")
candidate_snp<-unique(valid_data_eQTL$rsid)
writeLines(candidate_snp, "Candidate_snps_for_all_heart_eQTLs.txt")

# zcat GCST006414.tsv.gz | head -n 1 > Candidate_snps_for_all_heart_eQTLs.tsv
# zcat GCST006414.tsv.gz | awk 'NR==FNR{snps[$1]; next} FNR>1 && $2 in snps' Candidate_snps_for_all_heart_eQTLs.txt - >> Candidate_snps_for_all_heart_eQTLs.tsv

{
  heart_eqtl_res_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/heart_eqtl_res_for_all_heart.RDS")
  num=length(unique(heart_eqtl_res_for_all_heart $exposure))
  num
  adj_P=0.05/num
  sig_targets<-heart_eqtl_res_for_all_heart  %>% filter(pval<adj_P)
  length(unique(sig_targets$exposure)) #132
  
  sig_genes<-unique(sig_targets$exposure)
}

coloc_combine_results <- data.frame()
coloc_combine_snps <- data.frame()

data_outcome<- read.table("~/AF_atlas/Data/孟德尔随机化分析数据/GWAS数据/Candidate_snps_for_all_heart_eQTLs.tsv", 
                          header = TRUE, sep = "\t")
head(data_outcome)

{
  heart_eqtl_local_clumped_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/heart_eqtl_local_clumped_for_all_heart.RDS")
  filtered_exp_dat<-heart_eqtl_local_clumped_for_all_heart
  heart_eqtl_res_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/heart_eqtl_res_for_all_heart.RDS")
  res_sig<-heart_eqtl_res_for_all_heart 
  heart_eqtl_res_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/heart_eqtl_res_for_all_heart.RDS")
  num=length(unique(heart_eqtl_res_for_all_heart $exposure))
  num
  adj_P=0.05/num
}

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
  
  PP.H4.abf <- coloc_res$summary["PP.H4.abf"]
  coloc_results <-res_sig[res_sig$exposure==gene,]
  coloc_results <-coloc_results[,c("exposure","method","nsnp","b","se","pval")]
  coloc_results$PP.H4.abf=PP.H4.abf
  coloc_results$num_p_adj<-num
  coloc_results$adj_p_level<-0.05/num
  
  coloc_combine_results <- rbind(coloc_combine_results, coloc_results)
  
  coloc_snps <-subset(filtered_exp_dat,exposure==gene)
  coloc_snps <-coloc_snps[,c("pval.exposure","SNP","exposure")]
  coloc_combine_snps<-rbind(coloc_combine_snps, coloc_snps)
  
}

saveRDS(coloc_combine_results, file = "~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/coloc_combine_results.RDS")
saveRDS(coloc_combine_snps, file = "~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/coloc_combine_snps.RDS")


#---------------------------------------------心脏eQTL在芬兰数据库中复现-----------------------------------------
{
  heart_eqtl_res_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/heart_eqtl_res_for_all_heart.RDS")
  num=length(unique(heart_eqtl_res_for_all_heart $exposure))
  num
  adj_P=0.05/num
  sig_targets<-heart_eqtl_res_for_all_heart  %>% filter(pval<adj_P)
  length(unique(sig_targets$exposure)) #132
  
  sig_genes<-unique(sig_targets$exposure)
}

#从所有心脏eQTL数据中获取发现队列的有意义基因对应的SNP
heart_eqtl_local_clumped_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/heart_eqtl_local_clumped_for_all_heart.RDS")
head(heart_eqtl_local_clumped_for_all_heart)

filtered_exp_dat<-heart_eqtl_local_clumped_for_all_heart %>% filter(exposure %in% sig_genes)
length(unique(filtered_exp_dat$exposure))

file_content <- readLines("~/AF_atlas/Data/孟德尔随机化分析数据/GWAS数据/finngen_R9_I9_AF")
finngen <- read.table("~/AF_atlas/Data/孟德尔随机化分析数据/GWAS数据/finngen_R9_I9_AF", 
                                               header = TRUE, sep = "\t", col.names = unlist(strsplit(file_content[1], "\t")))
head(finngen,3)
# 将DataFrame保存为CSV文件
write.csv(finngen, file = "~/AF_atlas/Data/孟德尔随机化分析数据/GWAS数据/finngen.csv", row.names = FALSE)

filename="~/AF_atlas/Data/孟德尔随机化分析数据/GWAS数据/finngen.csv"

outcome_duplicate <- read_outcome_data(
  snps = filtered_exp_dat$SNP,
  filename = filename,
  sep = ",",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt",
  pval_col = "pval"
)

harmonise_data <- harmonise_data(exposure_dat = filtered_exp_dat, outcome_dat = outcome_duplicate)
res_duplicated <- mr(harmonise_data)
length(unique(res_duplicated$exposure))
saveRDS(res_duplicated, file = "~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/res_duplicated.RDS")

sig_duplicated<-res_duplicated %>% filter(pval<0.05)

length(unique(sig_duplicated$exposure))

#---------------------------------------------- 血eQTL孟德尔---------------------------------------------- 
Blood_eqtl_local_clumped <- readRDS("~/Lung_Atlas/单细胞孟德尔/文章部分二：孟德尔随机化/全血/Blood_eqtl_local_clumped.RDS")
head(Blood_eqtl_local_clumped)

filtered_exp_dat<-Blood_eqtl_local_clumped

setwd("~/AF_atlas/Data/孟德尔随机化分析数据/GWAS数据")
filename="GCST006414"
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
saveRDS(res, file = "~/AF_atlas/Data/孟德尔随机化分析数据/血eQTL数据/Blood_eqtl_res_for_all_heart.RDS")
Blood_eqtl_res_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/血eQTL数据/Blood_eqtl_res_for_all_heart.RDS")
num=length(unique(Blood_eqtl_res_for_all_heart$exposure))
num
adj_P=0.05/num

sig_targets<-Blood_eqtl_res_for_all_heart %>% filter(pval<adj_P)
length(unique(sig_targets$exposure))

#-------------------------------------------------------血eQTL 共定位-----------------------------------------
{
  Blood_eqtl_res_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/血eQTL数据/Blood_eqtl_res_for_all_heart.RDS")
  num=length(unique(Blood_eqtl_res_for_all_heart$exposure))
  num
  adj_P=0.05/num
  sig_targets<-Blood_eqtl_res_for_all_heart %>% filter(pval<adj_P)
  length(unique(sig_targets$exposure))
  
  sig_genes<-unique(sig_targets$exposure)
}


candidate_genes<-unique(sig_targets$exposure) # 151

setwd("~/Lung_Atlas/单细胞孟德尔/文章部分二：孟德尔随机化/GTEx")
writeLines(candidate_genes, "Candidate_genes_for_all_blood_eQTLs.txt")

# head Whole_Blood.tsv -n 1 > Candidate_genes_for_all_blood_eQTLs.tsv
# awk 'NR==FNR{genes[$1]; next} FNR!=1 && $7 in genes' Candidate_genes_for_all_blood_eQTLs.txt Whole_Blood.tsv >> Candidate_genes_for_all_blood_eQTLs.tsv

library(coloc)
data_eQTL <- read.table("~/AF_atlas/Data/孟德尔随机化分析数据/血eQTL数据/Candidate_genes_for_all_blood_eQTLs.tsv", header = TRUE, sep = "\t")
length(unique(data_eQTL$gene_id))
head(data_eQTL)

{
  # 假设data_eQTL是您的数据框，rsid是其中的一个列
  # 去除空白
  data_eQTL$rsid <- trimws(data_eQTL$rsid)
  
  # 检查格式 - 假设所有的variant_id都应以"rs"开头
  valid_format <- grepl("^rs", data_eQTL$rsid)
  
  # 排除NA值或空字符串
  valid_entries <- !is.na(data_eQTL$rsid) & data_eQTL$rsid != ""
  
  # 合并所有的检查
  valid_variant_ids <- valid_format & valid_entries
  
  # 提取有效的variant_id
  valid_data_eQTL <- data_eQTL[valid_variant_ids, ]
  
  # 检查无效的variant_id，这可以帮助您识别出问题所在
  invalid_data_eQTL <- data_eQTL[!valid_variant_ids, ]
  
}

setwd("~/AF_atlas/Data/孟德尔随机化分析数据/GWAS数据")
candidate_snp<-unique(valid_data_eQTL$rsid)
writeLines(candidate_snp, "Candidate_snps_for_all_blood_eQTLs.txt")

# zcat GCST006414.tsv.gz | head -n 1 > Candidate_snps_for_all_blood_eQTLs.tsv
# zcat GCST006414.tsv.gz | awk 'NR==FNR{snps[$1]; next} FNR>1 && $2 in snps' Candidate_snps_for_all_blood_eQTLs.txt - >> Candidate_snps_for_all_blood_eQTLs.tsv


coloc_combine_results <- data.frame()
coloc_combine_snps <- data.frame()

data_outcome<- read.table("~/AF_atlas/Data/孟德尔随机化分析数据/GWAS数据/Candidate_snps_for_all_blood_eQTLs.tsv", 
                          header = TRUE, sep = "\t")
head(data_outcome)

{
  Blood_eqtl_local_clumped <- readRDS("~/Lung_Atlas/单细胞孟德尔/文章部分二：孟德尔随机化/全血/Blood_eqtl_local_clumped.RDS")
  filtered_exp_dat<-Blood_eqtl_local_clumped
  Blood_eqtl_res_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/血eQTL数据/Blood_eqtl_res_for_all_heart.RDS")
  res_sig<-Blood_eqtl_res_for_all_heart
  Blood_eqtl_res_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/血eQTL数据/Blood_eqtl_res_for_all_heart.RDS")
  num=length(unique(Blood_eqtl_res_for_all_heart$exposure))
  num
  adj_P=0.05/num
}

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
    N = 838,
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
  
  PP.H4.abf <- coloc_res$summary["PP.H4.abf"]
  coloc_results <-res_sig[res_sig$exposure==gene,]
  coloc_results <-coloc_results[,c("exposure","method","nsnp","b","se","pval")]
  coloc_results$PP.H4.abf=PP.H4.abf
  coloc_results$num_p_adj<-num
  coloc_results$adj_p_level<-0.05/num
  
  coloc_combine_results <- rbind(coloc_combine_results, coloc_results)
  
  coloc_snps <-subset(filtered_exp_dat,exposure==gene)
  coloc_snps <-coloc_snps[,c("pval.exposure","SNP","exposure")]
  coloc_combine_snps<-rbind(coloc_combine_snps, coloc_snps)
  
}

saveRDS(coloc_combine_results, file = "~/AF_atlas/Data/孟德尔随机化分析数据/血eQTL数据/coloc_combine_results.RDS")
saveRDS(coloc_combine_snps, file = "~/AF_atlas/Data/孟德尔随机化分析数据/血eQTL数据/coloc_combine_snps.RDS")

#---------------------------------------------血eQTL在芬兰数据库中复现-----------------------------------------
{
  Blood_eqtl_res_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/血eQTL数据/Blood_eqtl_res_for_all_heart.RDS")
  num=length(unique(Blood_eqtl_res_for_all_heart$exposure))
  num
  adj_P=0.05/num
  sig_targets<-Blood_eqtl_res_for_all_heart %>% filter(pval<adj_P)
  length(unique(sig_targets$exposure)) # 151
  
  sig_genes<-unique(sig_targets$exposure)
}

#从所有血eQTL数据中获取发现队列的有意义基因对应的SNP
Blood_eqtl_local_clumped <- readRDS("~/Lung_Atlas/单细胞孟德尔/文章部分二：孟德尔随机化/全血/Blood_eqtl_local_clumped.RDS")
head(Blood_eqtl_local_clumped)

filtered_exp_dat<-Blood_eqtl_local_clumped %>% filter(exposure %in% sig_genes)
length(unique(filtered_exp_dat$exposure)) # 151

file_content <- readLines("~/AF_atlas/Data/孟德尔随机化分析数据/GWAS数据/finngen_R9_I9_AF")
finngen <- read.table("~/AF_atlas/Data/孟德尔随机化分析数据/GWAS数据/finngen_R9_I9_AF", 
                      header = TRUE, sep = "\t", col.names = unlist(strsplit(file_content[1], "\t")))
head(finngen,3)
# 将DataFrame保存为CSV文件
write.csv(finngen, file = "~/AF_atlas/Data/孟德尔随机化分析数据/GWAS数据/finngen.csv", row.names = FALSE)

filename="~/AF_atlas/Data/孟德尔随机化分析数据/GWAS数据/finngen.csv"

outcome_duplicate <- read_outcome_data(
  snps = filtered_exp_dat$SNP,
  filename = filename,
  sep = ",",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt",
  pval_col = "pval"
)

harmonise_data <- harmonise_data(exposure_dat = filtered_exp_dat, outcome_dat = outcome_duplicate)
res_duplicated <- mr(harmonise_data)
length(unique(res_duplicated$exposure))
saveRDS(res_duplicated, file = "~/AF_atlas/Data/孟德尔随机化分析数据/血eQTL数据/res_duplicated.RDS")

sig_duplicated<-res_duplicated %>% filter(pval<0.05)

length(unique(sig_duplicated$exposure))

#-------------------------------------- 心室孟德尔随机化-----------------------------------
heart_exp_dat <- read_exposure_data(
  filename = "~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/Heart_Left_Ventricle.tsv",
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
length(unique(heart_exp$exposure)) # 27

saveRDS(heart_exp, file = "~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/heart_exp_for_all_heart.RDS")

#---------------------------------------------- 心室本地clump---------------------------------------------- 
heart_exp <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/heart_exp_for_all_heart.RDS")
length(unique(heart_exp$exposure))

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

length(unique(heart_eqtl_local_clumped$exposure))

saveRDS(heart_eqtl_local_clumped, file = "~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/heart_eqtl_local_clumped_for_all_heart.RDS")

#-------------------------------------- 心室孟德尔随机化-----------------------------------
heart_exp_dat <- read_exposure_data(
  filename = "~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/Heart_Left_Ventricle.tsv",
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
length(unique(heart_exp$exposure)) # 27

saveRDS(heart_exp, file = "~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/heart_exp_for_all_heart.RDS")

#---------------------------------------------- 心室本地clump---------------------------------------------- 
heart_exp <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/heart_exp_for_all_heart.RDS")
length(unique(heart_exp$exposure))

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

length(unique(heart_eqtl_local_clumped$exposure))

saveRDS(heart_eqtl_local_clumped, file = "~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/heart_eqtl_local_clumped_for_all_heart.RDS")

#-------------------------------------------------------心室eQTL孟德尔-----------------------------------------
heart_eqtl_local_clumped_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/heart_eqtl_local_clumped_for_all_heart.RDS")
filtered_exp_dat<-heart_eqtl_local_clumped_for_all_heart

setwd("~/AF_atlas/Data/孟德尔随机化分析数据/GWAS数据")
filename="GCST006414"
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
saveRDS(res, file = "~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/heart_eqtl_res_for_all_heart.RDS")
heart_eqtl_res_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/heart_eqtl_res_for_all_heart.RDS")
num=length(unique(heart_eqtl_res_for_all_heart$exposure))
num
adj_P=0.05/num

sig_targets<-heart_eqtl_res_for_all_heart  %>% filter(pval<adj_P)
length(unique(sig_targets$exposure)) #114

#-------------------------------------------------------心室eQTL 共定位-----------------------------------------
candidate_genes<-unique(sig_targets$exposure) # 114

setwd("~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据")
writeLines(candidate_genes, "Candidate_genes_for_all_heart_eQTLs.txt")

# zcat Heart_Left_Ventricle.tsv.gz | head -n 1 > Candidate_genes_for_all_heart_eQTLs.tsv
# zcat Heart_Left_Ventricle.tsv.gz | awk 'NR==FNR{genes[$1]; next} FNR!=1 && $7 in genes' Candidate_genes_for_all_heart_eQTLs.txt - >> Candidate_genes_for_all_heart_eQTLs.tsv

library(coloc)
data_eQTL <- read.table("~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/Candidate_genes_for_all_heart_eQTLs.tsv", header = TRUE, sep = "\t")
length(unique(data_eQTL$gene_id))
head(data_eQTL)

{
  # 假设data_eQTL是您的数据框，rsid是其中的一个列
  # 去除空白
  data_eQTL$rsid <- trimws(data_eQTL$rsid)
  
  # 检查格式 - 假设所有的variant_id都应以"rs"开头
  valid_format <- grepl("^rs", data_eQTL$rsid)
  
  # 排除NA值或空字符串
  valid_entries <- !is.na(data_eQTL$rsid) & data_eQTL$rsid != ""
  
  # 合并所有的检查
  valid_variant_ids <- valid_format & valid_entries
  
  # 提取有效的variant_id
  valid_data_eQTL <- data_eQTL[valid_variant_ids, ]
  
  # 检查无效的variant_id，这可以帮助您识别出问题所在
  invalid_data_eQTL <- data_eQTL[!valid_variant_ids, ]
  
}

setwd("~/AF_atlas/Data/孟德尔随机化分析数据/GWAS数据")
candidate_snp<-unique(valid_data_eQTL$rsid)
writeLines(candidate_snp, "Candidate_vent_snps_for_all_heart_eQTLs.txt")

# zcat GCST006414.tsv.gz | head -n 1 > Candidate_vent_snps_for_all_heart_eQTLs.tsv
# zcat GCST006414.tsv.gz | awk 'NR==FNR{snps[$1]; next} FNR>1 && $2 in snps' Candidate_vent_snps_for_all_heart_eQTLs.txt - >> Candidate_vent_snps_for_all_heart_eQTLs.tsv

{
  heart_eqtl_res_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/heart_eqtl_res_for_all_heart.RDS")
  num=length(unique(heart_eqtl_res_for_all_heart $exposure))
  num
  adj_P=0.05/num
  sig_targets<-heart_eqtl_res_for_all_heart  %>% filter(pval<adj_P)
  length(unique(sig_targets$exposure)) #114
  
  sig_genes<-unique(sig_targets$exposure)
}

coloc_combine_results <- data.frame()
coloc_combine_snps <- data.frame()

data_outcome<- read.table("~/AF_atlas/Data/孟德尔随机化分析数据/GWAS数据/Candidate_vent_snps_for_all_heart_eQTLs.tsv", 
                          header = TRUE, sep = "\t")
head(data_outcome)

{
  heart_eqtl_local_clumped_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/heart_eqtl_local_clumped_for_all_heart.RDS")
  filtered_exp_dat<-heart_eqtl_local_clumped_for_all_heart
  heart_eqtl_res_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/heart_eqtl_res_for_all_heart.RDS")
  res_sig<-heart_eqtl_res_for_all_heart 
  heart_eqtl_res_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/heart_eqtl_res_for_all_heart.RDS")
  num=length(unique(heart_eqtl_res_for_all_heart $exposure))
  num
  adj_P=0.05/num
}

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
    N = 386,
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
  
  PP.H4.abf <- coloc_res$summary["PP.H4.abf"]
  coloc_results <-res_sig[res_sig$exposure==gene,]
  coloc_results <-coloc_results[,c("exposure","method","nsnp","b","se","pval")]
  coloc_results$PP.H4.abf=PP.H4.abf
  coloc_results$num_p_adj<-num
  coloc_results$adj_p_level<-0.05/num
  
  coloc_combine_results <- rbind(coloc_combine_results, coloc_results)
  
  coloc_snps <-subset(filtered_exp_dat,exposure==gene)
  coloc_snps <-coloc_snps[,c("pval.exposure","SNP","exposure")]
  coloc_combine_snps<-rbind(coloc_combine_snps, coloc_snps)
  
}

saveRDS(coloc_combine_results, file = "~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/coloc_combine_results.RDS")
saveRDS(coloc_combine_snps, file = "~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/coloc_combine_snps.RDS")

#---------------------------------------------心室eQTL在芬兰数据库中复现-----------------------------------------
{
  heart_eqtl_res_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/heart_eqtl_res_for_all_heart.RDS")
  num=length(unique(heart_eqtl_res_for_all_heart $exposure))
  num
  adj_P=0.05/num
  sig_targets<-heart_eqtl_res_for_all_heart  %>% filter(pval<adj_P)
  length(unique(sig_targets$exposure)) # 114
  
  sig_genes<-unique(sig_targets$exposure)
}

#从所有心室eQTL数据中获取发现队列的有意义基因对应的SNP
heart_eqtl_local_clumped_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/heart_eqtl_local_clumped_for_all_heart.RDS")
head(heart_eqtl_local_clumped_for_all_heart)

filtered_exp_dat<-heart_eqtl_local_clumped_for_all_heart %>% filter(exposure %in% sig_genes)
length(unique(filtered_exp_dat$exposure))

file_content <- readLines("~/AF_atlas/Data/孟德尔随机化分析数据/GWAS数据/finngen_R9_I9_AF")
finngen <- read.table("~/AF_atlas/Data/孟德尔随机化分析数据/GWAS数据/finngen_R9_I9_AF", 
                      header = TRUE, sep = "\t", col.names = unlist(strsplit(file_content[1], "\t")))
head(finngen,3)
# 将DataFrame保存为CSV文件
write.csv(finngen, file = "~/AF_atlas/Data/孟德尔随机化分析数据/GWAS数据/finngen.csv", row.names = FALSE)

filename="~/AF_atlas/Data/孟德尔随机化分析数据/GWAS数据/finngen.csv"

outcome_duplicate <- read_outcome_data(
  snps = filtered_exp_dat$SNP,
  filename = filename,
  sep = ",",
  snp_col = "rsids",
  beta_col = "beta",
  se_col = "sebeta",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "af_alt",
  pval_col = "pval"
)

harmonise_data <- harmonise_data(exposure_dat = filtered_exp_dat, outcome_dat = outcome_duplicate)
res_duplicated <- mr(harmonise_data)
length(unique(res_duplicated$exposure))
saveRDS(res_duplicated, file = "~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/res_duplicated.RDS")

sig_duplicated<-res_duplicated %>% filter(pval<0.05)

length(unique(sig_duplicated$exposure)) # 93

#---------------------- 心耳、心室和血有意义基因展示---------------------------------
adjust_b_values <- function(df, b_col, exposure_col) {
  # 对每个 exposure 分组
  df_grouped <- df %>% group_by(!!sym(exposure_col))
  
  # 对每个组应用函数来调整 b 值
  df_adjusted <- df_grouped %>% mutate(
    b_adjusted = case_when(
      sum(.data[[b_col]] > 0) > sum(.data[[b_col]] < 0) ~ abs(.data[[b_col]]),
      sum(.data[[b_col]] < 0) > sum(.data[[b_col]] > 0) ~ -abs(.data[[b_col]]),
      TRUE ~ .data[[b_col]] # 如果正负相等，则保持不变
    )
  ) %>% ungroup()
  
  # 检测正负相等的情况
  equal_counts <- df_adjusted %>% 
    group_by(!!sym(exposure_col)) %>% 
    summarise(
      Positive = sum(.data[[b_col]] > 0),
      Negative = sum(.data[[b_col]] < 0)
    ) %>% 
    filter(Positive == Negative)
  
  if (nrow(equal_counts) > 0) {
    cat("存在正负相等基因：", toString(equal_counts[[exposure_col]]), "\n")
  }
  
  return(df_adjusted)
}


genes_num_combine<-data.frame()

coloc_combine_results <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/coloc_combine_results.RDS")
res_duplicated <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/res_duplicated.RDS")

coloc_combine_results <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/coloc_combine_results.RDS")
res_duplicated <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/res_duplicated.RDS")

coloc_combine_results <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/血eQTL数据/coloc_combine_results.RDS")
res_duplicated <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/血eQTL数据/res_duplicated.RDS")

{
  Atrial_genes<-coloc_combine_results
  Atrial_genes <- Atrial_genes %>%
    group_by(exposure) %>%
    # 使用filter和if_else组合来选择每个组中的行
    filter(if(n() == 1 && method == "Wald ratio") {
      TRUE  # 如果只有一行且为Wald ratio，则保留
    } else if(n() >= 1 && any(method == "Inverse variance weighted")) {
      method == "Inverse variance weighted"  # 如果有多行，只保留Inverse variance weighted
    } else {
      FALSE  # 如果没有Wald ratio或Inverse variance weighted，则不保留
    }) %>%
    ungroup()  # 取消分组
  Atrial_genes_clean <- na.omit(Atrial_genes)
  Atrial_genes_adjusted <- adjust_b_values(Atrial_genes_clean, "b","exposure")
  all_sig_genes<-Atrial_genes_adjusted %>% filter(b_adjusted>0)
  all_sig_genes_num<-length(unique(all_sig_genes$exposure))
  
  pph4_genes<-filter(all_sig_genes,PP.H4.abf>0.8)
  pph4_genes_num<-length(unique(pph4_genes$exposure))
  
  res_duplicated <- res_duplicated %>%
    group_by(exposure) %>%
    # 使用filter和if_else组合来选择每个组中的行
    filter(if(n() == 1 && method == "Wald ratio") {
      TRUE  # 如果只有一行且为Wald ratio，则保留
    } else if(n() >= 1 && any(method == "Inverse variance weighted")) {
      method == "Inverse variance weighted"  # 如果有多行，只保留Inverse variance weighted
    } else {
      FALSE  # 如果没有Wald ratio或Inverse variance weighted，则不保留
    }) %>%
    ungroup()  # 取消分组
  res_duplicated_clean<-na.omit(res_duplicated)
  sig_duplicated <- adjust_b_values(res_duplicated_clean, "b","exposure")
  sig_duplicated_adjusted<-sig_duplicated %>% filter(pval<0.05 & b_adjusted>0)
  sig_duplicated_num<-length(unique(sig_duplicated_adjusted$exposure))
  
  sig_duplicated_pph4<-sig_duplicated_adjusted %>% filter(exposure %in% pph4_genes$exposure)
  sig_duplicated_pph4_num<-length(unique(sig_duplicated_pph4$exposure))
  
  genes_num<-data.frame(all_sig_genes_num=all_sig_genes_num,
                        pph4_genes_num=pph4_genes_num,
                        sig_duplicated_num,
                        sig_duplicated_pph4_num)
  
  genes_num_combine<-rbind(genes_num_combine,genes_num)
}

rowname<-c("Atrial","Ventricle","Blood")

rownames(genes_num_combine)<-rowname

genes_num_combine


genes_num_combine<-data.frame()

coloc_combine_results <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/coloc_combine_results.RDS")
res_duplicated <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/res_duplicated.RDS")

coloc_combine_results <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/coloc_combine_results.RDS")
res_duplicated <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/res_duplicated.RDS")

coloc_combine_results <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/血eQTL数据/coloc_combine_results.RDS")
res_duplicated <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/血eQTL数据/res_duplicated.RDS")


{
  Atrial_genes<-coloc_combine_results
  Atrial_genes <- Atrial_genes %>%
    group_by(exposure) %>%
    # 使用filter和if_else组合来选择每个组中的行
    filter(if(n() == 1 && method == "Wald ratio") {
      TRUE  # 如果只有一行且为Wald ratio，则保留
    } else if(n() >= 1 && any(method == "Inverse variance weighted")) {
      method == "Inverse variance weighted"  # 如果有多行，只保留Inverse variance weighted
    } else {
      FALSE  # 如果没有Wald ratio或Inverse variance weighted，则不保留
    }) %>%
    ungroup()  # 取消分组
  Atrial_genes_clean <- na.omit(Atrial_genes)
  Atrial_genes_adjusted <- adjust_b_values(Atrial_genes_clean, "b","exposure")
  all_sig_genes<-Atrial_genes_adjusted %>% filter(b_adjusted<0)
  all_sig_genes_num<-length(unique(all_sig_genes$exposure))
  
  pph4_genes<-filter(all_sig_genes,PP.H4.abf>0.8)
  pph4_genes_num<-length(unique(pph4_genes$exposure))
  
  res_duplicated <- res_duplicated %>%
    group_by(exposure) %>%
    # 使用filter和if_else组合来选择每个组中的行
    filter(if(n() == 1 && method == "Wald ratio") {
      TRUE  # 如果只有一行且为Wald ratio，则保留
    } else if(n() >= 1 && any(method == "Inverse variance weighted")) {
      method == "Inverse variance weighted"  # 如果有多行，只保留Inverse variance weighted
    } else {
      FALSE  # 如果没有Wald ratio或Inverse variance weighted，则不保留
    }) %>%
    ungroup()  # 取消分组
  res_duplicated_clean<-na.omit(res_duplicated)
  sig_duplicated <- adjust_b_values(res_duplicated_clean, "b","exposure")
  sig_duplicated_adjusted<-sig_duplicated %>% filter(pval<0.05 & b_adjusted<0)
  sig_duplicated_num<-length(unique(sig_duplicated_adjusted$exposure))
  
  sig_duplicated_pph4<-sig_duplicated_adjusted %>% filter(exposure %in% pph4_genes$exposure)
  sig_duplicated_pph4_num<-length(unique(sig_duplicated_pph4$exposure))
  
  genes_num<-data.frame(all_sig_genes_num=all_sig_genes_num,
                        pph4_genes_num=pph4_genes_num,
                        sig_duplicated_num,
                        sig_duplicated_pph4_num)
  
  genes_num_combine<-rbind(genes_num_combine,genes_num)
}

rowname<-c("Atrial","Ventricle","Blood")

rownames(genes_num_combine)<-rowname

genes_num_combine


#---------------------- 靶点基因森林图，只展示IVW方法和Ward方法---------------------------------
coloc_combine_results <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/coloc_combine_results.RDS")
res_duplicated <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/res_duplicated.RDS")

coloc_combine_results <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/coloc_combine_results.RDS")
res_duplicated <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/res_duplicated.RDS")

coloc_combine_results <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/血eQTL数据/coloc_combine_results.RDS")
res_duplicated <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/血eQTL数据/res_duplicated.RDS")

{
  Atrial_genes<-coloc_combine_results
  Atrial_genes <- Atrial_genes %>%
    group_by(exposure) %>%
    # 使用filter和if_else组合来选择每个组中的行
    filter(if(n() == 1 && method == "Wald ratio") {
      TRUE  # 如果只有一行且为Wald ratio，则保留
    } else if(n() >= 1 && any(method == "Inverse variance weighted")) {
      method == "Inverse variance weighted"  # 如果有多行，只保留Inverse variance weighted
    } else {
      FALSE  # 如果没有Wald ratio或Inverse variance weighted，则不保留
    }) %>%
    ungroup()  # 取消分组
  Atrial_genes_clean <- na.omit(Atrial_genes)
  Atrial_genes_adjusted <- adjust_b_values(Atrial_genes_clean, "b","exposure")
  all_sig_genes<-Atrial_genes_adjusted
  all_sig_genes_num<-length(unique(all_sig_genes$exposure))
  
  pph4_genes<-filter(all_sig_genes,PP.H4.abf>0.8)
  pph4_genes_num<-length(unique(pph4_genes$exposure))
  
  res_duplicated <- res_duplicated %>%
    group_by(exposure) %>%
    # 使用filter和if_else组合来选择每个组中的行
    filter(if(n() == 1 && method == "Wald ratio") {
      TRUE  # 如果只有一行且为Wald ratio，则保留
    } else if(n() >= 1 && any(method == "Inverse variance weighted")) {
      method == "Inverse variance weighted"  # 如果有多行，只保留Inverse variance weighted
    } else {
      FALSE  # 如果没有Wald ratio或Inverse variance weighted，则不保留
    }) %>%
    ungroup()  # 取消分组
  res_duplicated_clean<-na.omit(res_duplicated)
  sig_duplicated <- adjust_b_values(res_duplicated_clean, "b","exposure")
  sig_duplicated_adjusted<-sig_duplicated %>% filter(pval<0.05)
  sig_duplicated_pph4<-sig_duplicated_adjusted %>% filter(exposure %in% pph4_genes$exposure)
}
length(unique(sig_duplicated_pph4$exposure))

heart_eqtl_res_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/heart_eqtl_res_for_all_heart.RDS")
heart_eqtl_res_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/heart_eqtl_res_for_all_heart.RDS")
heart_eqtl_res_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/血eQTL数据/Blood_eqtl_res_for_all_heart.RDS")

subset_show_genes<-filter(heart_eqtl_res_for_all_heart,exposure %in% unique(sig_duplicated_pph4$exposure))
subset_show_genes <- subset_show_genes %>%
  group_by(exposure) %>%
  # 使用filter和if_else组合来选择每个组中的行
  filter(if(n() == 1 && method == "Wald ratio") {
    TRUE  # 如果只有一行且为Wald ratio，则保留
  } else if(n() >= 1 && any(method == "Inverse variance weighted")) {
    method == "Inverse variance weighted"  # 如果有多行，只保留Inverse variance weighted
  } else {
    FALSE  # 如果没有Wald ratio或Inverse variance weighted，则不保留
  }) %>%
  ungroup()  # 取消分组
subset_show_genes <- na.omit(subset_show_genes)

subset_show_genes<-generate_odds_ratios(subset_show_genes)
head(subset_show_genes)
length(unique(subset_show_genes$exposure))


library(clusterProfiler)
library(org.Hs.eg.db)

# 假设subset_show_genes$exposure包含了ENSEMBL IDs
genes <- subset_show_genes$exposure

genes_mapped <-mapIds(org.Hs.eg.db,
                      keys = genes,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")

# 将genes_mapped转换为数据框
genes_mapped_df <- data.frame(exposure = names(genes_mapped), symbol = genes_mapped, stringsAsFactors = FALSE)

# 合并genes_mapped_df和subset_show_genes
subset_show_genes <- merge(subset_show_genes, genes_mapped_df, by = "exposure", all.x = TRUE)

subset_show_genes <- subset_show_genes %>%
  mutate(symbol = ifelse(is.na(symbol), exposure, symbol))

subset_show_genes

library(dplyr)

# 替换特定的symbol值
subset_show_genes <- subset_show_genes %>%
  mutate(symbol = ifelse(exposure == "ENSG00000243225", "RPLP1(pseudo)", symbol))

subset_show_genes <- subset_show_genes %>%
  mutate(symbol = ifelse(exposure == "ENSG00000246323", "FAM13B-AS1", symbol))

subset_show_genes <- subset_show_genes %>%
  mutate(symbol = ifelse(exposure == "ENSG00000250260", "Lnc-WNT8A-4", symbol))

subset_show_genes <- subset_show_genes %>%
  mutate(symbol = ifelse(exposure == "ENSG00000246323", "FAM13B-AS1", symbol))

subset_show_genes <- subset_show_genes %>%
  mutate(symbol = ifelse(exposure == "ENSG00000250091", "DNAH10OS", symbol))

subset_show_genes <- subset_show_genes %>%
  mutate(symbol = ifelse(exposure == "ENSG00000273151", "Lnc-ADAP1-1", symbol))

subset_show_genes

ggplot(subset_show_genes, aes(x = reorder(symbol, or), y = or, color = or)) +
  geom_point() +
  geom_errorbar(aes(ymin = or_lci95, ymax = or_uci95), width = 0) +
  scale_color_gradient(low = "blue", high = "red") +
  coord_flip() +  # 翻转坐标轴，使得exposure在y轴
  theme_minimal() +
  labs(title = "Odds Ratios with 95% CI for Different Exposures",
       x = "Odds Ratio (OR)",
       y = "Exposure") +
  theme(legend.position = "right") +
  # 当or > 1时，箭头位于误差线的上端
  geom_segment(data = subset(subset_show_genes, or > 1),
               aes(x = reorder(symbol, or), xend = reorder(symbol, or), 
                   y = or_uci95, yend = or_uci95 + 0.01),
               arrow = arrow(type = "closed", length = unit(0.1, "cm"))) +
  # 当or < 1时，箭头位于误差线的下端
  geom_segment(data = subset(subset_show_genes, or < 1),
               aes(x = reorder(symbol, or), xend = reorder(symbol, or), 
                   y = or_lci95, yend = or_lci95 - 0.01),
               arrow = arrow(type = "closed", length = unit(0.1, "cm")))


#-------------------------------不同数据集气泡图展示---------------------------------------------------
library(tidyverse)
library(ggrepel)
library(ggfun)
library(grid)

####----load data----####

heart_eqtl_res_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/heart_eqtl_res_for_all_heart.RDS")
coloc_combine_results <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/coloc_combine_results.RDS")

heart_eqtl_res_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/heart_eqtl_res_for_all_heart.RDS")
coloc_combine_results <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/coloc_combine_results.RDS")

heart_eqtl_res_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/血eQTL数据/Blood_eqtl_res_for_all_heart.RDS")
coloc_combine_results <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/血eQTL数据/coloc_combine_results.RDS")


subset_show_genes<-heart_eqtl_res_for_all_heart
subset_show_genes <- subset_show_genes %>%
  group_by(exposure) %>%
  # 使用filter和if_else组合来选择每个组中的行
  filter(if(n() == 1 && method == "Wald ratio") {
    TRUE  # 如果只有一行且为Wald ratio，则保留
  } else if(n() >= 1 && any(method == "Inverse variance weighted")) {
    method == "Inverse variance weighted"  # 如果有多行，只保留Inverse variance weighted
  } else {
    FALSE  # 如果没有Wald ratio或Inverse variance weighted，则不保留
  }) %>%
  ungroup()  # 取消分组
subset_show_genes <- na.omit(subset_show_genes)

head(subset_show_genes)

subset_show_genes<-generate_odds_ratios(subset_show_genes)

df<-data.frame(SYMBOL=subset_show_genes$exposure,
               log2FoldChange=log2(subset_show_genes$or),
               padj= subset_show_genes$pval)

head(df)
df$change <- ifelse(df$log2FoldChange > 0, "Up", 
                    ifelse(df$log2FoldChange < 0, "Down", "Normal"))

# 假设subset_show_genes$exposure包含了ENSEMBL IDs
df$exposure<-df$SYMBOL
head(df)
genes<-df$SYMBOL

genes_mapped <-mapIds(org.Hs.eg.db,
                      keys = genes,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")

# 将genes_mapped转换为数据框
genes_mapped_df <- data.frame(exposure = names(genes_mapped), symbol = genes_mapped, stringsAsFactors = FALSE)

# 合并genes_mapped_df和subset_show_genes
subset_show_genes <- merge(df, genes_mapped_df, by = "exposure", all.x = TRUE)
head(subset_show_genes)
subset_show_genes <- subset_show_genes %>%
  mutate(symbol = ifelse(is.na(symbol), exposure, symbol))

df<-subset_show_genes

adj_p_level<-unique(coloc_combine_results$adj_p_level)

####----plot----####
ggplot(data = df) + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj), 
                 color = log2FoldChange,
                 size = -log10(padj))) + 
  geom_point(data =  df %>%
               tidyr::drop_na() %>%
               #dplyr::filter(change != "Normal") %>%
               dplyr::arrange(desc(-log10(padj))) %>%
               dplyr::slice(1:20),
             aes(x = log2FoldChange, y = -log10(padj),
                 # fill = log2FoldChange,
                 size = -log10(padj)),
             shape = 21, show.legend = F, color = "#000000") +
  geom_text_repel(data =  df %>% 
                    tidyr::drop_na() %>% 
                    dplyr::filter(change != "Normal") %>%
                    dplyr::arrange(desc(-log10(padj))) %>%
                    dplyr::slice(1:15) %>%
                    dplyr::filter(change == "Up"),
                  aes(x = log2FoldChange, y = -log10(padj), label = symbol),
                  box.padding = 0.5,
                  nudge_x = 0.5,
                  nudge_y = 0.2,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  # segment.angle = 10,
                  direction = "y", 
                  hjust = "left",
                  max.overlaps = Inf
  ) + 
  geom_text_repel(data =  df %>% 
                    tidyr::drop_na() %>% 
                    #dplyr::filter(change != "Normal") %>%
                    dplyr::arrange(desc(-log10(padj))) %>%
                    dplyr::slice(1:15) %>%
                    dplyr::filter(change == "Down"),
                  aes(x = log2FoldChange, y = -log10(padj), label = symbol),
                  box.padding = 0.5,
                  nudge_x = -0.2,
                  nudge_y = 0.2,
                  segment.curvature = -0.1,
                  segment.ncp = 3,
                  segment.angle = 20,
                  direction = "y", 
                  hjust = "right",
                  max.overlaps = Inf
  ) + 
  scale_color_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                        values = seq(0, 1, 0.2)) +
  scale_fill_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                       values = seq(0, 1, 0.2)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = -log10(adj_p_level), linetype = 4) + 
  scale_size(range = c(1,7)) + 
  xlim(c(-3, 3)) + 
  ylim(c(-1, 90)) + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        legend.background = element_roundrect(color = "#808080", linetype = 1),
        axis.text = element_text(size = 13, color = "#000000"),
        axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
  ) + 
  annotate(geom = "text", x = 2.5, y = 0.25, label = "p = 0.05", size = 5) + 
  coord_cartesian(clip = "off") + 
  annotation_custom(
    grob = grid::segmentsGrob(
      y0 = unit(-10, "pt"),
      y1 = unit(-10, "pt"),
      arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "first"),
      gp = grid::gpar(lwd = 3, col = "#74add1")
    ), 
    xmin = -3, 
    xmax = -1,
    ymin = 90,
    ymax = 90
  ) +
  annotation_custom(
    grob = grid::textGrob(
      label = "Down",
      gp = grid::gpar(col = "#74add1")
    ),
    xmin = -3, 
    xmax = -1,
    ymin = 90,
    ymax = 90
  ) +
  annotation_custom(
    grob = grid::segmentsGrob(
      y0 = unit(-10, "pt"),
      y1 = unit(-10, "pt"),
      arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "last"),
      gp = grid::gpar(lwd = 3, col = "#d73027")
    ), 
    xmin = 3, 
    xmax = 1,
    ymin = 90,
    ymax = 90
  ) +
  annotation_custom(
    grob = grid::textGrob(
      label = "Up",
      gp = grid::gpar(col = "#d73027")
    ),
    xmin = 3, 
    xmax = 1,
    ymin = 90,
    ymax = 90
  ) 

#-------------------------------心脏心室和血基因交集可视化-----------------------------------------

adjust_b_values <- function(df, b_col, exposure_col) {
  # 对每个 exposure 分组
  df_grouped <- df %>% group_by(!!sym(exposure_col))
  
  # 对每个组应用函数来调整 b 值
  df_adjusted <- df_grouped %>% mutate(
    b_adjusted = case_when(
      sum(.data[[b_col]] > 0) > sum(.data[[b_col]] < 0) ~ abs(.data[[b_col]]),
      sum(.data[[b_col]] < 0) > sum(.data[[b_col]] > 0) ~ -abs(.data[[b_col]]),
      TRUE ~ .data[[b_col]] # 如果正负相等，则保持不变
    )
  ) %>% ungroup()
  
  # 检测正负相等的情况
  equal_counts <- df_adjusted %>% 
    group_by(!!sym(exposure_col)) %>% 
    summarise(
      Positive = sum(.data[[b_col]] > 0),
      Negative = sum(.data[[b_col]] < 0)
    ) %>% 
    filter(Positive == Negative)
  
  if (nrow(equal_counts) > 0) {
    cat("存在正负相等基因：", toString(equal_counts[[exposure_col]]), "\n")
  }
  
  return(df_adjusted)
}


coloc_combine_results <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/coloc_combine_results.RDS")
res_duplicated <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/res_duplicated.RDS")

coloc_combine_results <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/coloc_combine_results.RDS")
res_duplicated <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/res_duplicated.RDS")

coloc_combine_results <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/血eQTL数据/coloc_combine_results.RDS")
res_duplicated <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/血eQTL数据/res_duplicated.RDS")

{
  Atrial_genes<-coloc_combine_results
  Atrial_genes <- Atrial_genes %>%
    group_by(exposure) %>%
    # 使用filter和if_else组合来选择每个组中的行
    filter(if(n() == 1 && method == "Wald ratio") {
      TRUE  # 如果只有一行且为Wald ratio，则保留
    } else if(n() >= 1 && any(method == "Inverse variance weighted")) {
      method == "Inverse variance weighted"  # 如果有多行，只保留Inverse variance weighted
    } else {
      FALSE  # 如果没有Wald ratio或Inverse variance weighted，则不保留
    }) %>%
    ungroup()  # 取消分组
  Atrial_genes_clean <- na.omit(Atrial_genes)
  Atrial_genes_adjusted <- adjust_b_values(Atrial_genes_clean, "b","exposure")
  all_sig_genes<-Atrial_genes_adjusted
  all_sig_genes_num<-length(unique(all_sig_genes$exposure))
  
  pph4_genes<-filter(all_sig_genes,PP.H4.abf>0.8)
  pph4_genes_num<-length(unique(pph4_genes$exposure))
  
  res_duplicated <- res_duplicated %>%
    group_by(exposure) %>%
    # 使用filter和if_else组合来选择每个组中的行
    filter(if(n() == 1 && method == "Wald ratio") {
      TRUE  # 如果只有一行且为Wald ratio，则保留
    } else if(n() >= 1 && any(method == "Inverse variance weighted")) {
      method == "Inverse variance weighted"  # 如果有多行，只保留Inverse variance weighted
    } else {
      FALSE  # 如果没有Wald ratio或Inverse variance weighted，则不保留
    }) %>%
    ungroup()  # 取消分组
  res_duplicated_clean<-na.omit(res_duplicated)
  sig_duplicated <- adjust_b_values(res_duplicated_clean, "b","exposure")
  sig_duplicated_adjusted<-sig_duplicated %>% filter(pval<0.05)
  sig_duplicated_pph4<-sig_duplicated_adjusted %>% filter(exposure %in% pph4_genes$exposure)
}
length(unique(sig_duplicated_pph4$exposure))

heart_eqtl_res_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心脏eQTL数据/heart_eqtl_res_for_all_heart.RDS")
heart_eqtl_res_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/心室eQTL数据/heart_eqtl_res_for_all_heart.RDS")
heart_eqtl_res_for_all_heart <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/血eQTL数据/Blood_eqtl_res_for_all_heart.RDS")

subset_show_genes<-filter(heart_eqtl_res_for_all_heart,exposure %in% unique(sig_duplicated_pph4$exposure))
subset_show_genes <- subset_show_genes %>%
  group_by(exposure) %>%
  # 使用filter和if_else组合来选择每个组中的行
  filter(if(n() == 1 && method == "Wald ratio") {
    TRUE  # 如果只有一行且为Wald ratio，则保留
  } else if(n() >= 1 && any(method == "Inverse variance weighted")) {
    method == "Inverse variance weighted"  # 如果有多行，只保留Inverse variance weighted
  } else {
    FALSE  # 如果没有Wald ratio或Inverse variance weighted，则不保留
  }) %>%
  ungroup()  # 取消分组
subset_show_genes <- na.omit(subset_show_genes)

subset_show_genes<-generate_odds_ratios(subset_show_genes)
head(subset_show_genes)
length(unique(subset_show_genes$exposure))


library(clusterProfiler)
library(org.Hs.eg.db)

# 假设subset_show_genes$exposure包含了ENSEMBL IDs
genes <- subset_show_genes$exposure

genes_mapped <-mapIds(org.Hs.eg.db,
                      keys = genes,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")

# 将genes_mapped转换为数据框
genes_mapped_df <- data.frame(exposure = names(genes_mapped), symbol = genes_mapped, stringsAsFactors = FALSE)

# 合并genes_mapped_df和subset_show_genes
subset_show_genes <- merge(subset_show_genes, genes_mapped_df, by = "exposure", all.x = TRUE)

subset_show_genes <- subset_show_genes %>%
  mutate(symbol = ifelse(is.na(symbol), exposure, symbol))

subset_show_genes

library(dplyr)

# 替换特定的symbol值
subset_show_genes <- subset_show_genes %>%
  mutate(symbol = ifelse(exposure == "ENSG00000243225", "RPLP1(pseudo)", symbol))

subset_show_genes <- subset_show_genes %>%
  mutate(symbol = ifelse(exposure == "ENSG00000246323", "FAM13B-AS1", symbol))

subset_show_genes <- subset_show_genes %>%
  mutate(symbol = ifelse(exposure == "ENSG00000250260", "Lnc-WNT8A-4", symbol))

subset_show_genes <- subset_show_genes %>%
  mutate(symbol = ifelse(exposure == "ENSG00000246323", "FAM13B-AS1", symbol))

subset_show_genes <- subset_show_genes %>%
  mutate(symbol = ifelse(exposure == "ENSG00000250091", "DNAH10OS", symbol))

subset_show_genes <- subset_show_genes %>%
  mutate(symbol = ifelse(exposure == "ENSG00000273151", "Lnc-ADAP1-1", symbol))

subset_show_genes

ven1<-data.frame(symbol=subset_show_genes$symbol,or=subset_show_genes$or,source="Atrial",exposure=subset_show_genes$exposure)
head(ven1)
ven2<-data.frame(symbol=subset_show_genes$symbol,or=subset_show_genes$or,source="Ventricle",exposure=subset_show_genes$exposure)
head(ven2)
ven3<-data.frame(symbol=subset_show_genes$symbol,or=subset_show_genes$or,source="Blood",exposure=subset_show_genes$exposure)
head(ven3)

ven_combine<-rbind(ven1,ven2,ven3)
ven_combine

ven_combine$category <- ifelse(ven_combine$or > 1, "positive", "negative")

ven_combine_unique <- ven_combine %>%  distinct()
ven_combine_unique
saveRDS(ven_combine_unique, file = "~/AF_atlas/Data/孟德尔随机化分析数据/ven_combine_unique.RDS")

ven_combine_unique <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/ven_combine_unique.RDS")

d1=data.frame(from="origin", to=paste("group", seq(1,6), sep=""))
d2=data.frame(from=rep(d1$to[1], each=19), to=paste("subgroup", seq(1,19), sep="_"))
d3=data.frame(from=rep(d1$to[2], each=16), to=paste("subgroup", seq(20,35), sep="_"))
d4=data.frame(from=rep(d1$to[3], each=11), to=paste("subgroup", seq(36,46), sep="_"))
d5=data.frame(from=rep(d1$to[4], each=11), to=paste("subgroup", seq(47,57), sep="_"))
d6=data.frame(from=rep(d1$to[5], each=8), to=paste("subgroup", seq(58,65), sep="_"))
d7=data.frame(from=rep(d1$to[6], each=1), to=paste("subgroup", seq(66,66), sep="_"))

edges=rbind(d1, d2,d3,d4,d5,d6,d7)
num=1+6+19+16+11+11+8+1

vertices = data.frame(
  name = unique(c(as.character(edges$from), as.character(edges$to))) , 
  value = runif(num)
) 

vertices$group = edges$from[ match( vertices$name, edges$to ) ]
head(vertices)

used_data<-data.frame(symbol=ven_combine_unique$symbol,
                      or=ven_combine_unique$or,
                      source=ven_combine_unique$source,
                      category=ven_combine_unique$category)


used_data$source <- factor(used_data$source, levels = c('Atrial', 'Ventricle', 'Blood'))
used_data$category <- factor(used_data$category, levels = c('positive', 'negative'))

# Order the data by 'source', then by 'category', and finally by 'or'
used_data <- used_data[order(used_data$source, used_data$category, used_data$or), ]

# If you need to reset row names to be sequential after sorting
rownames(used_data) <- seq(nrow(used_data))

# View the sorted data
print(used_data)

subgroup_rows <- grep("subgroup_", vertices$name)
num_rows <- length(subgroup_rows)

# 检查used_data的行数是否足够
if (nrow(used_data) >= num_rows) {
  vertices$value[subgroup_rows] <- used_data$or[1:num_rows]
  vertices$name[subgroup_rows] <- used_data$symbol[1:num_rows]
} else {
  stop("used_data中的行数不足以替换所有subgroup_的值。")
}

# 查看替换后的vertices
print(vertices)

length(vertices$name)

# 找出重复的顶点名称
dup_names <- vertices$name[duplicated(vertices$name)]

# 为每个重复的顶点名称添加后缀以使其唯一
for (name in dup_names) {
  # 找到所有重复名称的行
  dup_rows <- which(vertices$name == name)
  
  # 为这些名称添加后缀（除了第一个）
  vertices$name[dup_rows] <- paste0(name, "_", seq_along(dup_rows))
}

# 确保现在所有的顶点名称都是唯一的
if (!any(duplicated(vertices$name))) {
  print("所有顶点名称现在都是唯一的。")
}

length(vertices$name)
vertices
print(edges)

# 找到vertices中以group或origin开头的行的索引
group_origin_start_index <- which(grepl("^(group|origin)", vertices$name))

# 查看找到的索引
print(group_origin_start_index)

# 确保存在足够的subgroup用于替换
if (nrow(vertices) - length(group_origin_start_index) < 66) {
  stop("vertices中的subgroup数量不足以替换edges中所有subgroup的值。")
}

# 替换edges中的subgroup
for (i in 1:66) {
  subgroup_name <- paste0("subgroup_", i)
  subgroup_rows <- which(edges$to == subgroup_name)
  if (length(subgroup_rows) > 0) {
    edges$to[subgroup_rows] <- vertices$name[length(group_origin_start_index) + i]
  }
}

# 查看替换后的edges
print(edges)
vertices

vertices$id=NA
myleaves=which(is.na( match(vertices$name, edges$from) ))
nleaves=length(myleaves)
vertices$id[ myleaves ] = seq(1:nleaves)
vertices$angle= 90 - 360 * vertices$id / nleaves

vertices$hjust<-ifelse( vertices$angle < -90, 1, 0)

# flip angle BY to make them readable
vertices$angle<-ifelse(vertices$angle < -90, vertices$angle+180, vertices$angle)

# 查看测试数据
head(edges)

head(vertices)

# Create a graph object
mygraph <- graph_from_data_frame( edges, vertices=vertices )

# Make the plot
ggraph(mygraph, layout = 'dendrogram', circular = TRUE) + 
  geom_edge_diagonal(colour="grey") + #设置节点边的颜色
  # 设置节点的标签，字体大小，文本注释信息
  geom_node_text(aes(x = x*1.15, y=y*1.15, filter = leaf, label=name, angle = angle, hjust=hjust*1, colour=group), size=2.5, alpha=1) +
  # 设置节点的大小，颜色和透明度
  geom_node_point(aes(filter = leaf, x = x*1.07, y=y*1.07, colour=group, size=value, alpha=1)) +
  # 设置颜色的画板
  scale_colour_manual(values= rep( brewer.pal(8,"Dark2") , 3)) +
  # 设置节点大小的范围
  scale_size_continuous( range = c(1,10) ) +
  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(c(0,0,0,0),"cm"),
  ) +
  expand_limits(x = c(-1.6, 1.6), y = c(-1.6, 1.6))


#--------------------------------------------绘制韦恩图---------------------------------------

ven_combine_unique <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/ven_combine_unique.RDS")

data_show<-data.frame(symbol=ven_combine_unique$symbol,
                      or=ven_combine_unique$or,
                      source=ven_combine_unique$source,
                      sub=ven_combine_unique$category)

#devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
data<-data_show

atrial <- subset(data, source == "Atrial")
ventricle <- subset(data, source == "Ventricle")
blood <- subset(data, source == "Blood")

# 用list构建集合
list_data <- list(
  Atrial = atrial$symbol,
  Ventricle = ventricle$symbol,
  Blood = blood$symbol
)

# 绘制韦恩图
ggvenn(list_data,fill_color = c("#add8e6", "#90ee90", "#dda0dd"),stroke_size=0.3)

#-------------------------------有意义基因在房颤队列中的表达-----------------------------------------
ven_combine_unique <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/ven_combine_unique.RDS")
ven_combine_unique

setwd("~/AF_atlas/Data/GSE")
gene_data <- read.table("GSE128188.txt", header = TRUE, sep = "\t", check.names = FALSE)
head(gene_data)

library(edgeR)
library(ggplot2)
library(ggpubr)

plot_expression_with_stats <- function(data, expression_column, group_column, title="Expression Across Groups") {
  # 检查每个组的正态性
  groups <- unique(data[[group_column]])
  normality_tests <- sapply(groups, function(group) {
    group_data <- data[data[[group_column]] == group, expression_column]
    
    # 确保 group_data 是数值型
    if (is.numeric(group_data)) {
      return(shapiro.test(group_data)$p.value)
    } else {
      return(NA)  # 非数值数据返回 NA
    }
  })
  
  # 判断是否使用 t 检验（排除 NA 值）
  use_t_test <- all(normality_tests > 0.05, na.rm = TRUE)
  
  # 选择统计方法
  method <- if(use_t_test) "t.test" else "wilcox"
  
  # 绘制箱线图
  p <- ggplot(data, aes(x = !!sym(group_column), y = !!sym(expression_column), fill = !!sym(group_column))) +
    geom_boxplot() +
    geom_jitter(width = 1, size = 1.5, alpha = 0.6) +
    scale_fill_manual(values = c("SR_LA" = "#F3E7EF", "AF_LA" = "#CBB8D4", "SR_RA" = "#DAE5F0", "AF_RA" = "#B9C9E3")) +
    labs(title = title, x = "Group", y = "Expression") +
    theme_minimal() +
    stat_compare_means(method = method, comparisons = list(c("SR_LA", "AF_LA"), c("SR_RA", "AF_RA")), label = "p.signif")
  
  return(p)
}


row.names(gene_data) <- gene_data$Gene
gene_data <- gene_data[, -1]
head(gene_data)
boxplot(gene_data)

# 创建 sample_mapping 数据框
sample_mapping <- data.frame(
  short_name = c("7-RA", "7-LA", "12b-RA", "12b-LA", "20-RA", "20-LA", "27-RA", "27-LA", "30-RA", "30-LA", "5-RA", "5-LA", "14-RA", "14-LA", "25-RA", "25-LA", "28-RA", "28-LA", "31-RA", "31-LA"),
  full_name = c("SR_RA_7", "SR_LA_7", "SR_RA_12b", "SR_LA_12b", "SR_RA_20", "SR_LA_20", "SR_RA_27", "SR_LA_27", "SR_RA_30", "SR_LA_30", "AF_RA_5", "AF_LA_5", "AF_RA_14", "AF_LA_14", "AF_RA_25", "AF_LA_25", "AF_RA_28", "AF_LA_28", "AF_RA_31", "AF_LA_31")
)

colnames(gene_data) <- sample_mapping$full_name[match(colnames(gene_data), sample_mapping$short_name)]
head(gene_data)

group <- sapply(colnames(gene_data), function(col_name) {
  rhythm <- ifelse(grepl("SR", col_name), "SR", "AF")
  chamber <- ifelse(grepl("LA", col_name), "LA", "RA")
  paste(rhythm, chamber, sep = "_")
})

print(group)

data <- DGEList(counts = gene_data, group = group)

keep <- rowSums(cpm(data) > 1) >= 2
table(keep)
data <- data[keep, , keep.lib.sizes = FALSE]
dim(data)
data$samples$lib.size <- colSums(data$counts)
data <- calcNormFactors(data)

data<-cpm(data)
boxplot(data)

data<-log(data+1)
boxplot(data)

head(data)

synpo2l_data <- data.frame(SYNPO2L = data["GMCL1", ], group = group)
long_data <- pivot_longer(synpo2l_data, cols = -group, names_to = "sample", values_to = "expression")

# 定义分组的顺序
long_data$group <- factor(long_data$group, levels = c("SR_LA", "AF_LA", "SR_RA", "AF_RA"))

plot_expression_with_stats(long_data, "expression", "group")

library(patchwork)

plot_gene_expression <- function(gene) {
  synpo2l_data <- data.frame(Expression = data[gene, ], group = group)
  long_data <- pivot_longer(synpo2l_data, cols = -group, names_to = "sample", values_to = "expression")
  long_data$group <- factor(long_data$group, levels = c("SR_LA", "AF_LA", "SR_RA", "AF_RA"))
  plot_expression_with_stats(long_data, "expression", "group", title = paste("Expression of", gene, "Across Groups"))
}

ven_combine_unique <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/ven_combine_unique.RDS")
atrial_genes<-filter(ven_combine_unique,source=="Atrial") %>% arrange(-or)
Ventricle_genes<-filter(ven_combine_unique,source=="Ventricle")%>% arrange(-or)
Blood_genes<-filter(ven_combine_unique,source=="Blood")%>% arrange(-or)

# 为ven_combine_unique$symbol中的每个基因生成一个图
gene_plots <- lapply(atrial_genes$symbol[atrial_genes$symbol %in% rownames(data)][1:27], plot_gene_expression)
gene_plots <- lapply(Ventricle_genes$symbol[Ventricle_genes$symbol %in% rownames(data)][1:18], plot_gene_expression)
gene_plots <- lapply(Blood_genes$symbol[Blood_genes$symbol %in% rownames(data)][1:7], plot_gene_expression)

# 使用patchwork将所有图组合成一个图
combined_plot <- wrap_plots(gene_plots, ncol = 5)

# 显示组合图
combined_plot

genes_chosen<-c("SYNPO2L","SLC27A6","GTF2I","THRB")
genes_chosen<-c("GYPC","MTSS1","ASAH1","SYNE2")
genes_chosen<-c("ZBTB38","GSDMB","ORMDL3","COG5")

gene_plots <- lapply(genes_chosen, plot_gene_expression)

# 使用patchwork将所有图组合成一个图
combined_plot <- wrap_plots(gene_plots, ncol = 4)

# 显示组合图
combined_plot

plot_gene_expression("SYNPO2L")

#-------------------------------有意义基因富集分析-----------------------------------------
library(clusterProfiler)

ven_combine_unique <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/ven_combine_unique.RDS")
atrial_genes<-filter(ven_combine_unique,source=="Atrial") 

keys <- mapIds(org.Hs.eg.db,
               keys = atrial_genes$exposure,
               column = "ENTREZID",
               keytype = "ENSEMBL",
               multiVals = "first")

go <- enrichGO(
  gene = keys,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.5,
  readable=F) #是否将gene ID转换为gene symbol

dotplot(go, split="ONTOLOGY",showCategory = 10,label_format=50)+ 
  facet_grid(ONTOLOGY~.,scale="free")+
  theme(panel.grid = element_blank())+#修改主题
  theme(axis.title =element_text(size = 12, color = 'black'),
        axis.text.y =element_text(size = 12),
        legend.title=element_text(size=12))+
  scale_color_gradient(high="#FC8D62",low="#4b5cc4")#设置颜色


#-------------------------- 展示细胞密度图-----------------------------
AF_renamed <- readRDS("~/AF_atlas/Data/单细胞分析/AF_renamed.RDS")
Idents(AF_renamed)<-AF_renamed$cell_type
DimPlot(AF_renamed,label = T,cols = pal_30)
DimPlot(AF_renamed,cols = pal_30)
table(AF_renamed$condition)
Idents(AF_renamed)<-AF_renamed$condition

Ctrl<-subset(AF_renamed,idents="Ctrl")
AF<-subset(AF_renamed,idents="AF")

object<-Ctrl
object<-AF
{
  library(viridis)
  coord = Embeddings(object = object, reduction = "umap")
  coord = coord[,c(1,2)]
  colnames(coord) = c("UMAP_1", "UMAP_2")
  coord = data.frame(ID = rownames(coord), coord)
  meta = object@meta.data
  meta = data.frame(ID = rownames(meta), meta,stringsAsFactors = F)
  meta = left_join(meta, coord, by = 'ID')
  
  theme_black <- function(base_size = 12, base_family = "") {
    theme_grey(base_size = base_size, base_family = base_family) %+replace%
      theme(
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        panel.background = element_rect(fill = "black", color  =  NA),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.background = element_rect(fill = "grey30", color = "grey10"),
        strip.text.x = element_text(size = base_size*0.8, color = "white"),
        strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),
        plot.background = element_rect(color = "black", fill = "black"),
        plot.title = element_text(size = base_size*1.2, color = "white"),
        plot.margin = unit(rep(0, 4), "lines")
      )
  }
  
  ggplot(data = coord, mapping = aes(x = UMAP_1, y = UMAP_2)) +
    stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
    geom_point(color = 'white', size = 0.05) +
    scale_fill_viridis(option = "magma") +
    theme_black()
}



#-------------------------------有意义基因在不同细胞的表达-----------------------------------------
AF_renamed <- readRDS("~/AF_atlas/Data/单细胞分析/AF_renamed.RDS")

ven_combine_unique <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/ven_combine_unique.RDS")
atrial_genes<-filter(ven_combine_unique,source=="Atrial") %>% arrange(-or)
Ventricle_genes<-filter(ven_combine_unique,source=="Ventricle")%>% arrange(-or)

table(Idents(Atrial))
#devtools::install_github('junjunlab/scRNAtoolVis')
library(scRNAtoolVis)
AverageHeatmap(object = AF_renamed,markerGene = Ventricle_genes$symbol,
               assays = "SCT",cluster_columns = T,
               cluster_rows = T,htCol = c("#1f78b4","#a6cee3","#fc8d62"))

object<-AF_renamed
Idents(object)<-object$condition

AF<-subset(object,idents="AF")
Ctrl<-subset(object,idents="Ctrl")
Idents(AF)<-AF$cell_type
Idents(Ctrl)<-Ctrl$cell_type

p1<-AverageHeatmap(object = Ctrl,markerGene = atrial_genes$symbol,
                   assays = "SCT")
p2<-AverageHeatmap(object = AF,markerGene = atrial_genes$symbol,
                   assays = "SCT")

p3<-AverageHeatmap(object = Ctrl,markerGene = Ventricle_genes$symbol,
                   assays = "SCT")
p4<-AverageHeatmap(object = AF,markerGene = Ventricle_genes$symbol,
                   assays = "SCT")

p1+p2

p3+p4

#----------------------------------------------在房颤中计算基因分数----------------------------------------

ven_combine_unique <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/ven_combine_unique.RDS")
data_used<-AF_renamed
data_used$condition <- factor(data_used$condition, levels = c("Ctrl", "AF"))
col<-c(pal_30[2],pal_30[1])

atrial_genes<-filter(ven_combine_unique,source=="Atrial") %>% arrange(or)%>% filter(category=="positive")
atrial_genes_positive<-list(atrial_genes$symbol)
data_used<-AddModuleScore(data_used,features = atrial_genes_positive,name = "atrial_genes_positive")

VlnPlot(data_used,features = "atrial_genes_positive1",pt.size = 0,sort = "increasing",
        split.by = "condition",split.plot = T)+ 
  scale_fill_manual(values = col)+ geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.shape = NA)

atrial_genes<-filter(ven_combine_unique,source=="Atrial") %>% arrange(or)%>% filter(category=="negative")
atrial_genes_negative<-list(atrial_genes$symbol)
data_used<-AddModuleScore(data_used,features = atrial_genes_negative,name = "atrial_genes_negative")

VlnPlot(data_used,features = "atrial_genes_negative1",pt.size = 0,sort = "increasing",
        split.by = "condition",split.plot = T)+ 
  scale_fill_manual(values = col)+ geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.shape = NA)

Ventricle_genes<-filter(ven_combine_unique,source=="Ventricle")%>% arrange(or)%>% filter(category=="positive")
Ventricle_genes_positive<-list(Ventricle_genes$symbol)
data_used<-AddModuleScore(data_used,features = Ventricle_genes_positive,name = "Ventricle_genes_positive")

VlnPlot(data_used,features = "Ventricle_genes_positive1",pt.size = 0,sort = "increasing",
        split.by = "condition",split.plot = T)+ 
  scale_fill_manual(values = col)+ geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.shape = NA)

Ventricle_genes<-filter(ven_combine_unique,source=="Ventricle")%>% arrange(or)%>% filter(category=="negative")
Ventricle_genes_negative<-list(Ventricle_genes$symbol)
data_used<-AddModuleScore(data_used,features = Ventricle_genes_negative,name = "Ventricle_genes_negative")

VlnPlot(data_used,features = "Ventricle_genes_negative1",pt.size = 0,sort = "increasing",
        split.by = "condition",split.plot = T)+ 
  scale_fill_manual(values = col)+ geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.shape = NA)

#---------------------------------------------画图展示in silico genes------------------------------------

library(ggplot2)
library(dplyr)
library(readr)
library(readr)

# 读取数据
deletion <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/T_cell/T_cell.csv", header = TRUE, sep = ",")
overexpression <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/T_cell/T_cell_over.csv", header = TRUE, sep = ",")

deletion <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/FB/FB.csv", header = TRUE, sep = ",")
overexpression <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/FB/FB_over.csv", header = TRUE, sep = ",")

deletion <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/EC/EC.csv", header = TRUE, sep = ",")
overexpression <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/EC/EC_over.csv", header = TRUE, sep = ",")

deletion <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/MP/MP.csv", header = TRUE, sep = ",")
overexpression <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/MP/MP_over.csv", header = TRUE, sep = ",")

deletion <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/DC/DC.csv", header = TRUE, sep = ",")
overexpression <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/DC/DC_over.csv", header = TRUE, sep = ",")

deletion <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/SMC/SMC.csv", header = TRUE, sep = ",")
overexpression <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/SMC/SMC_over.csv", header = TRUE, sep = ",")

deletion <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/Neutrophil/Neutrophil.csv", header = TRUE, sep = ",")
overexpression <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/Neutrophil/Neutrophil_over.csv", header = TRUE, sep = ",")


{
  deletion_filtered <- filter(deletion, Sig == 1)
  overexpression_filtered <- filter(overexpression, Sig == 1)
  deletion_filtered$Category <- 'deletion'
  overexpression_filtered$Category <- 'overexpression'
  common_genes <- intersect(deletion_filtered$Ensembl_ID, overexpression_filtered$Ensembl_ID)
  length(unique(common_genes))
  
  deletion_subset <- deletion_filtered[deletion_filtered$Ensembl_ID %in% common_genes, ]
  overexpression_subset <- overexpression_filtered[overexpression_filtered$Ensembl_ID%in% common_genes, ]
  
  head(deletion_subset)
  head(overexpression_subset)
  
  common_ids <- merge(deletion_subset, overexpression_subset, by = "Ensembl_ID")
  
  common_ids <- common_ids[sign(common_ids$Shift_to_goal_end.x) == sign(common_ids$Shift_to_goal_end.y), ]
  
  common_ids
  dim(common_ids)
  
  common_ids$source <- ifelse(abs(common_ids$Shift_to_goal_end.x) > abs(common_ids$Shift_to_goal_end.y), 
                              "deletion", 
                              "overexpression")
  
  # 2. 创建最终的子集
  final_subset <- rbind(
    deletion_subset[deletion_subset$Ensembl_ID %in% common_ids$Ensembl_ID[common_ids$source == "deletion"], ],
    overexpression_subset[overexpression_subset$Ensembl_ID %in% common_ids$Ensembl_ID[common_ids$source == "overexpression"], ]
  )
  
  # 查看结果
  head(final_subset)
  dim(final_subset)
  final_subset
  
  # 找到Category为deletion的Ensembl_ID
  deletion_ids <- final_subset %>% 
    filter(Category == "deletion") %>% 
    .$Ensembl_ID
  
  # 找到Category为overexpression的Ensembl_ID
  overexpression_ids <- final_subset %>% 
    filter(Category == "overexpression") %>% 
    .$Ensembl_ID
  
  # 更新overexpression_filtered：去除Category为deletion的Ensembl_ID对应的行
  overexpression_final <- overexpression_filtered %>% 
    filter(!Ensembl_ID %in% deletion_ids)
  
  # 更新deletion_filtered：去除Category为overexpression的Ensembl_ID对应的行
  deletion_final <- deletion_filtered %>% 
    filter(!Ensembl_ID %in% overexpression_ids)
  
  # 查看更新后的数据集
  head(deletion_final)
  head(overexpression_final)
  
  merged_data <- rbind(deletion_final, overexpression_final)
  head(merged_data)
  dim(merged_data)
  
  # 排序基因，以便绘图时有序
  merged_data <- merged_data %>%
    arrange(desc(Shift_to_goal_end))
  head(merged_data)
  
  # 为绘图生成行号（排名）
  merged_data$Rank <- seq_along(merged_data$Shift_to_goal_end)
  
  # 选择颜色方案
  color_palette <- brewer.pal(n = 3, name = "Set1")
  #color_palette
  #color_palette <- c("#377EB8" ,"#E41A1C","#4DAF4A" )
  
  # 绘制图形
  ggplot(merged_data, aes(x = Rank, y = Shift_to_goal_end, group = Category, color = Category)) +
    geom_line() +
    geom_point(aes(color = Category), size = 2) +
    scale_color_manual(values = color_palette) +
    theme_minimal() +
    labs(x = "Rank", y = "Shift to normal", title = "In silico treatment") +
    geom_text_repel(data = merged_data %>% 
                      tidyr::drop_na() %>% 
                      dplyr::arrange(desc(Shift_to_goal_end)) %>%
                      dplyr::slice(1:5),
                    aes(x = Rank, y = Shift_to_goal_end, label = Gene_name),
                    box.padding = 0.3,
                    nudge_x = 100,
                    nudge_y = 0.0001,
                    segment.curvature = -0.1,
                    segment.ncp = 3,
                    direction = "y", 
                    hjust = "left",
                    max.overlaps = Inf
    )+
    geom_text_repel(data = merged_data %>% 
                      tidyr::drop_na() %>% 
                      dplyr::arrange(desc(-Shift_to_goal_end)) %>%
                      dplyr::slice(1:5),
                    aes(x = Rank, y = Shift_to_goal_end, label = Gene_name),
                    box.padding = 0.3,
                    nudge_x = -150,
                    nudge_y = 0.0001,
                    segment.curvature = 0.1,
                    segment.ncp = 3,
                    direction = "y", 
                    hjust = "left",
                    max.overlaps = Inf
    )
}

length(unique(deletion_final$Ensembl_ID))
length(unique(overexpression_final$Ensembl_ID))


#---------------------------------------雷达图展示in silico 基因-----------------------------
# T cell 164/509 
# FB 145/1623
# EC 63/1348
# MP 100/290
# DC 71/525
# SMC 61/294
# Neutrophil 84/236
devtools::install_github("ricardo-bion/ggradar")
library(ggradar)

data <- data.frame(
  Group = c("deletion", "overexpression"),
  T_cell = c(164, 509),
  FB = c(145, 1623),
  EC = c(63, 1348),
  MP = c(100, 290),
  DC = c(71, 525),
  SMC = c(61, 294),
  Neutrophil = c(84, 236)
)
data
ggradar(data,grid.min = 50,grid.mid = 100,grid.max = 1700,centre.y = 10,
        values.radar = c(20, 90, 180))

data <- data.frame(
  Group = c("deletion"),
  T_cell = c(164),
  FB = c(145),
  EC = c(63),
  MP = c(100),
  DC = c(71),
  SMC = c(61),
  Neutrophil = c(84)
)

ggradar(data,grid.min = 20,grid.mid = 90,grid.max = 180,centre.y = 10,
        values.radar = c(20, 90, 180))

data <- data.frame(
  Group = c("overexpression"),
  T_cell = c(509),
  FB = c(1623),
  EC = c(1348),
  MP = c(290),
  DC = c(525),
  SMC = c(294),
  Neutrophil = c(236)
)
data
ggradar(data,grid.min = 200,grid.mid = 1000,grid.max = 1700,centre.y = 0,
        values.radar = c(200, 500, 1700))

data <- data.frame(
  Group = c("deletion", "overexpression"),
  T_cell = c(164, 509),
  FB = c(145, 600),
  EC = c(63, 600),
  MP = c(100, 290),
  DC = c(71, 525),
  SMC = c(61, 294),
  Neutrophil = c(84, 236)
)
data
ggradar(data,grid.min = 50,grid.mid = 300,grid.max = 600,centre.y = 10,
        plot.legend = F,group.point.size = 4,
        group.line.width = 1)

#---------------------------------------整合所有in silico 基因-----------------------------

# 读取数据
deletion <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/T_cell/T_cell.csv", header = TRUE, sep = ",")
overexpression <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/T_cell/T_cell_over.csv", header = TRUE, sep = ",")

deletion <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/FB/FB.csv", header = TRUE, sep = ",")
overexpression <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/FB/FB_over.csv", header = TRUE, sep = ",")

deletion <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/EC/EC.csv", header = TRUE, sep = ",")
overexpression <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/EC/EC_over.csv", header = TRUE, sep = ",")

deletion <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/MP/MP.csv", header = TRUE, sep = ",")
overexpression <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/MP/MP_over.csv", header = TRUE, sep = ",")

deletion <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/DC/DC.csv", header = TRUE, sep = ",")
overexpression <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/DC/DC_over.csv", header = TRUE, sep = ",")

deletion <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/SMC/SMC.csv", header = TRUE, sep = ",")
overexpression <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/SMC/SMC_over.csv", header = TRUE, sep = ",")

deletion <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/Neutrophil/Neutrophil.csv", header = TRUE, sep = ",")
overexpression <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/Neutrophil/Neutrophil_over.csv", header = TRUE, sep = ",")

{
  deletion_filtered <- filter(deletion, Sig == 1)
  overexpression_filtered <- filter(overexpression, Sig == 1)
  deletion_filtered$Category <- 'deletion'
  overexpression_filtered$Category <- 'overexpression'
  common_genes <- intersect(deletion_filtered$Ensembl_ID, overexpression_filtered$Ensembl_ID)
  length(unique(common_genes))
  
  deletion_subset <- deletion_filtered[deletion_filtered$Ensembl_ID %in% common_genes, ]
  overexpression_subset <- overexpression_filtered[overexpression_filtered$Ensembl_ID%in% common_genes, ]
  
  head(deletion_subset)
  head(overexpression_subset)
  
  common_ids <- merge(deletion_subset, overexpression_subset, by = "Ensembl_ID")
  
  common_ids <- common_ids[sign(common_ids$Shift_to_goal_end.x) == sign(common_ids$Shift_to_goal_end.y), ]
  
  common_ids
  dim(common_ids)
  
  common_ids$source <- ifelse(abs(common_ids$Shift_to_goal_end.x) > abs(common_ids$Shift_to_goal_end.y), 
                              "deletion", 
                              "overexpression")
  
  # 2. 创建最终的子集
  final_subset <- rbind(
    deletion_subset[deletion_subset$Ensembl_ID %in% common_ids$Ensembl_ID[common_ids$source == "deletion"], ],
    overexpression_subset[overexpression_subset$Ensembl_ID %in% common_ids$Ensembl_ID[common_ids$source == "overexpression"], ]
  )
  
  # 查看结果
  head(final_subset)
  dim(final_subset)
  final_subset
  
  # 找到Category为deletion的Ensembl_ID
  deletion_ids <- final_subset %>% 
    filter(Category == "deletion") %>% 
    .$Ensembl_ID
  
  # 找到Category为overexpression的Ensembl_ID
  overexpression_ids <- final_subset %>% 
    filter(Category == "overexpression") %>% 
    .$Ensembl_ID
  
  # 更新overexpression_filtered：去除Category为deletion的Ensembl_ID对应的行
  overexpression_final <- overexpression_filtered %>% 
    filter(!Ensembl_ID %in% deletion_ids)
  
  # 更新deletion_filtered：去除Category为overexpression的Ensembl_ID对应的行
  deletion_final <- deletion_filtered %>% 
    filter(!Ensembl_ID %in% overexpression_ids)
  
  # 查看更新后的数据集
  head(deletion_final)
  head(overexpression_final)
  
  merged_data <- rbind(deletion_final, overexpression_final)
}

head(merged_data)
names(merged_data) <- paste0(names(merged_data), "_T_cell")
names(merged_data) <- paste0(names(merged_data), "_FB")
names(merged_data) <- paste0(names(merged_data), "_EC")
names(merged_data) <- paste0(names(merged_data), "_MP")
names(merged_data) <- paste0(names(merged_data), "_DC")
names(merged_data) <- paste0(names(merged_data), "_SMC")
names(merged_data) <- paste0(names(merged_data), "_Neutrophil")

head(merged_data)

T_cell<-merged_data
dim(T_cell)
FB<-merged_data
dim(FB)
EC<-merged_data
dim(EC)
MP<-merged_data
dim(MP)
DC<-merged_data
dim(DC)
SMC<-merged_data
dim(SMC)
Neutrophil<-merged_data
dim(Neutrophil)

combine_data<-merge(T_cell,FB,by.x="Ensembl_ID_T_cell",by.y="Ensembl_ID_FB",all=T)
combine_data<-merge(combine_data,EC,by.x="Ensembl_ID_T_cell",by.y="Ensembl_ID_EC",all=T)
combine_data<-merge(combine_data,MP,by.x="Ensembl_ID_T_cell",by.y="Ensembl_ID_MP",all=T)
combine_data<-merge(combine_data,DC,by.x="Ensembl_ID_T_cell",by.y="Ensembl_ID_DC",all=T)
combine_data<-merge(combine_data,SMC,by.x="Ensembl_ID_T_cell",by.y="Ensembl_ID_SMC",all=T)
combine_data<-merge(combine_data,Neutrophil,by.x="Ensembl_ID_T_cell",by.y="Ensembl_ID_Neutrophil",all=T)

head(combine_data)
dim(combine_data)

#saveRDS(combine_data, file = "~/AF_atlas/Data/in silico genes/combined_data_所有疾病靶点基因.RDS")

combined_data_CMs <- readRDS("~/AF_atlas/Data/in silico genes/combined_data_CMs.RDS")
head(combined_data_CMs)

combined_data_所有疾病靶点基因 <- readRDS("~/AF_atlas/Data/in silico genes/combined_data_所有疾病靶点基因.RDS")
combine_data_all_cells<-merge(combined_data_所有疾病靶点基因,combined_data_CMs,
                              by.x="Ensembl_ID_T_cell",by.y="Ensembl_ID_CM",all=T)
head(combine_data_all_cells)
table(combine_data_all_cells$Category_T_cell)
saveRDS(combine_data_all_cells, file = "~/AF_atlas/Data/in silico genes/combined_data_所有疾病靶点基因.RDS")

#-----------------deletion
T_cell_deletion<-T_cell %>% filter(Category_T_cell=="deletion")
FB_deletion<-FB %>% filter(Category_FB=="deletion")
EC_deletion<-EC %>% filter(Category_EC=="deletion")
MP_deletion<-MP %>% filter(Category_MP=="deletion")
DC_deletion<-DC %>% filter(Category_DC=="deletion")
SMC_deletion<-SMC %>% filter(Category_SMC=="deletion")
Neutrophil_deletion<-Neutrophil %>% filter(Category_Neutrophil=="deletion")

deletion_subset<-merge(T_cell_deletion,FB_deletion,by.x="Ensembl_ID_T_cell",by.y="Ensembl_ID_FB",all=T)
deletion_subset<-merge(deletion_subset,EC_deletion,by.x="Ensembl_ID_T_cell",by.y="Ensembl_ID_EC",all=T)
deletion_subset<-merge(deletion_subset,MP_deletion,by.x="Ensembl_ID_T_cell",by.y="Ensembl_ID_MP",all=T)
deletion_subset<-merge(deletion_subset,DC_deletion,by.x="Ensembl_ID_T_cell",by.y="Ensembl_ID_DC",all=T)
deletion_subset<-merge(deletion_subset,SMC_deletion,by.x="Ensembl_ID_T_cell",by.y="Ensembl_ID_SMC",all=T)
deletion_subset<-merge(deletion_subset,Neutrophil_deletion,by.x="Ensembl_ID_T_cell",by.y="Ensembl_ID_Neutrophil",all=T)

saveRDS(deletion_subset, file = "~/AF_atlas/Data/in silico genes/combined_deletion_所有deletion疾病靶点基因.RDS")

combined_deletion_所有deletion疾病靶点基因 <- readRDS("~/AF_atlas/Data/in silico genes/combined_deletion_所有deletion疾病靶点基因.RDS")
deletion_subset<-combined_deletion_所有deletion疾病靶点基因
exps_columns <- grep("^Shift_to_goal_end_|Ensembl_ID_T_cell", colnames(deletion_subset))
exps_in_silico<-deletion_subset[,exps_columns]
exps_in_silico_unique <- unique(exps_in_silico)
rownames(exps_in_silico_unique)<-exps_in_silico_unique$Ensembl_ID_T_cell

#-----------------overexpression
T_cell_over<-T_cell %>% filter(Category_T_cell=="overexpression")
FB_over<-FB %>% filter(Category_FB=="overexpression")
EC_over<-EC %>% filter(Category_EC=="overexpression")
MP_over<-MP %>% filter(Category_MP=="overexpression")
DC_over<-DC %>% filter(Category_DC=="overexpression")
SMC_over<-SMC %>% filter(Category_SMC=="overexpression")
Neutrophil_over<-Neutrophil %>% filter(Category_Neutrophil=="overexpression")

overexpression_subset<-merge(T_cell_over,FB_over,by.x="Ensembl_ID_T_cell",by.y="Ensembl_ID_FB",all=T)
overexpression_subset<-merge(overexpression_subset,EC_over,by.x="Ensembl_ID_T_cell",by.y="Ensembl_ID_EC",all=T)
overexpression_subset<-merge(overexpression_subset,MP_over,by.x="Ensembl_ID_T_cell",by.y="Ensembl_ID_MP",all=T)
overexpression_subset<-merge(overexpression_subset,DC_over,by.x="Ensembl_ID_T_cell",by.y="Ensembl_ID_DC",all=T)
overexpression_subset<-merge(overexpression_subset,SMC_over,by.x="Ensembl_ID_T_cell",by.y="Ensembl_ID_SMC",all=T)
overexpression_subset<-merge(overexpression_subset,Neutrophil_over,by.x="Ensembl_ID_T_cell",by.y="Ensembl_ID_Neutrophil",all=T)

saveRDS(overexpression_subset, file = "~/AF_atlas/Data/in silico genes/combined_over_所有over疾病靶点基因.RDS")

combined_over_所有over疾病靶点基因 <- readRDS("~/AF_atlas/Data/in silico genes/combined_over_所有over疾病靶点基因.RDS")
overexpression_subset<-combined_over_所有over疾病靶点基因
exps_columns <- grep("^Shift_to_goal_end_|Ensembl_ID_T_cell", colnames(overexpression_subset))
exps_in_silico<-overexpression_subset[,exps_columns]
head(exps_in_silico)
exps_in_silico_unique <- unique(exps_in_silico)
rownames(exps_in_silico_unique)<-exps_in_silico_unique$Ensembl_ID_T_cell

head(exps_in_silico_unique)
exps_in_silico_unique<-exps_in_silico_unique[,-1]
head(exps_in_silico_unique)

library(ClusterGVis)
any(is.na(exps_in_silico_unique))
exps_in_silico_unique[is.na(exps_in_silico_unique)] <- 0

symbols <- mapIds(org.Hs.eg.db, 
                  keys = rownames(exps_in_silico_unique), 
                  column = "SYMBOL", 
                  keytype = "ENSEMBL", 
                  multiVals = "first")

valid_indices <- !is.na(symbols)
dim(exps_in_silico_unique)
exps_in_silico_unique <- exps_in_silico_unique[valid_indices, ]
dim(exps_in_silico_unique)
symbols_no_na <- symbols[valid_indices]

# 将有效的基因符号设置为行名
rownames(exps_in_silico_unique) <- symbols_no_na

scaled_exps_in_silico <- scale(exps_in_silico_unique)

ck <- clusterData(exp = scaled_exps_in_silico,
                  cluster.method = "mfuzz",
                  cluster.num = 7,
                  seed = 42,
                  scaleData = F)

#获取要展示的基因
derta<-0.12
{
  percentage_to_extract <- 0.10 *derta
  
  # 计算每个簇中元素的数量
  cluster_sizes <- table(ck$wide.res$cluster)
  
  # 计算每个簇应抽取的元素数量
  num_to_extract <- round(cluster_sizes * percentage_to_extract)
  
  # 从每个簇中随机抽取索引
  set.seed(123)  # 设置随机种子以获得可重现的结果
  extracted_indices <- unlist(mapply(function(size, n) sample(size, n), size = cluster_sizes, n = num_to_extract, SIMPLIFY = FALSE))
  
  # 获取对应的基因名称
  markGenes <- rownames(exps_in_silico_unique)[extracted_indices]
}

head(ck$wide.res)
cluster_gene_df <- ck$wide.res[, c("gene", "cluster")]
# 查看新的数据框
head(cluster_gene_df)

# 创建一个空的数据框来存储结果
termBP <- data.frame(id = character(0), term = character(0), pval = numeric(0))

# 获取不同的 cluster 值
clusters <- unique(cluster_gene_df$cluster)

# 遍历每个 cluster
for (cluster_value in clusters) {
  # 从 cluster_gene_df 中选择特定 cluster 的子集
  cluster_subset <- subset(cluster_gene_df, cluster == cluster_value)
  
  # 获取基因的映射
  mapped_entrez_ids <- mapIds(org.Hs.eg.db,
                              keys = cluster_subset$gene,
                              column = "ENTREZID",
                              keytype = "SYMBOL",
                              multiVals = "first")
  
  # 进行 GO 富集分析
  go_BP_result <- enrichGO(gene = mapped_entrez_ids,
                           OrgDb = org.Hs.eg.db,
                           ont = "BP",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)
  
  # 提取 GO 富集结果的 Description 和 qvalue 列，并添加 cluster 的 id
  cluster_BP <- go_BP_result@result[, c("Description", "qvalue")]
  cluster_BP$id <- cluster_value
  
  # 将结果追加到 termanno2 数据框中
  termBP <- rbind(termBP, cluster_BP)
}

# 查看结果
head(termBP)

BP <- data.frame(
  id = termBP$id, 
  term = termBP$Description, 
  pval = termBP$qvalue)

BP <- BP[BP$pval < 0.05, ]

cluster_labels <- paste0("C", 1:length(unique(termBP$id)))
# 将数值标签映射为标签
BP$id <- cluster_labels[match(BP$id, unique(termBP$id))]

saveRDS(BP, file = "~/AF_atlas/Data/in silico genes/BP_用于画图的deletion的BP富集分析.RDS")

library(dplyr)
# 根据id分组，然后按照pval升序排列，取前三个
top_terms <- BP %>%
  group_by(id) %>%
  arrange(pval) %>%
  slice_head(n = 5)

setwd("~/AF_atlas/Data/in silico genes")

pal <- c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", 
         "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", 
         "#E8601C", "#DC050C")

pal_30<-c("#EBA6B0","#6085BF","#8E82B0","#E4756A","#E66A6C",
          "#DCB9CE","#8BC078","#E56B7F","#F4B151","#F5BC53",
          "#65C0E8","#35AC9A","#EB9671","#31A5D9","#DE5758",
          "#C597B9","#D64448","#C2197E","#6AB8B1","#D579A5",
          "#C8B3C9","#C28339","#BEBFD7","#31A5D9","#7580B5",
          "#F49F2E","#7064A0","#87C9D3","#3A6CAF","#375BA4")

pdf('termlf.pdf',height = 10,width = 12,onefile = F)

cluster_num<-length(unique(ck$wide.res$cluster))
sample_num<-ncol(ck$wide.res)-3
length(top_terms$term)
visCluster(object = ck,
           plot.type = "both",
           column_names_rot = 45,
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = top_terms,
           line.side = "left",
           add.box = T,
           sample.col = pal_30[1:sample_num],
           ctAnno.col = pal[1:cluster_num],
           go.col = "Black"
)
dev.off()

#----------------------------------------------非心肌细胞与silico基因取交集----------------------------------------

combined_data_所有疾病靶点基因 <- readRDS("~/AF_atlas/Data/in silico genes/combined_data_所有疾病靶点基因.RDS")
head(combined_data_所有疾病靶点基因)
ven_combine_unique <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/ven_combine_unique.RDS")
head(ven_combine_unique)

inter_genes<-intersect(ven_combine_unique$exposure,unique(combined_data_所有疾病靶点基因$Ensembl_ID_T_cell))

inter_in_eQTLs_to_silico<-ven_combine_unique %>% filter(exposure %in% inter_genes)
inter_in_silico_to_eQTLs<-combined_data_所有疾病靶点基因 %>% filter(Ensembl_ID_T_cell %in% inter_genes)


#------------------------------------------心肌细胞与silico基因取交集-----------------------------------------------

deletion <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/CM/CM.csv", header = TRUE, sep = ",")
overexpression <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/CM/CM_over.csv", header = TRUE, sep = ",")

{
  deletion_filtered <- filter(deletion, Sig == 1)
  overexpression_filtered <- filter(overexpression, Sig == 1)
  deletion_filtered$Category <- 'deletion'
  overexpression_filtered$Category <- 'overexpression'
  common_genes <- intersect(deletion_filtered$Ensembl_ID, overexpression_filtered$Ensembl_ID)
  length(unique(common_genes))
  
  deletion_subset <- deletion_filtered[deletion_filtered$Ensembl_ID %in% common_genes, ]
  overexpression_subset <- overexpression_filtered[overexpression_filtered$Ensembl_ID%in% common_genes, ]
  
  head(deletion_subset)
  head(overexpression_subset)
  
  common_ids <- merge(deletion_subset, overexpression_subset, by = "Ensembl_ID")
  
  common_ids <- common_ids[sign(common_ids$Shift_to_goal_end.x) == sign(common_ids$Shift_to_goal_end.y), ]
  
  common_ids
  dim(common_ids)
  
  common_ids$source <- ifelse(abs(common_ids$Shift_to_goal_end.x) > abs(common_ids$Shift_to_goal_end.y), 
                              "deletion", 
                              "overexpression")
  
  # 2. 创建最终的子集
  final_subset <- rbind(
    deletion_subset[deletion_subset$Ensembl_ID %in% common_ids$Ensembl_ID[common_ids$source == "deletion"], ],
    overexpression_subset[overexpression_subset$Ensembl_ID %in% common_ids$Ensembl_ID[common_ids$source == "overexpression"], ]
  )
  
  # 查看结果
  head(final_subset)
  dim(final_subset)
  final_subset
  
  # 找到Category为deletion的Ensembl_ID
  deletion_ids <- final_subset %>% 
    filter(Category == "deletion") %>% 
    .$Ensembl_ID
  
  # 找到Category为overexpression的Ensembl_ID
  overexpression_ids <- final_subset %>% 
    filter(Category == "overexpression") %>% 
    .$Ensembl_ID
  
  # 更新overexpression_filtered：去除Category为deletion的Ensembl_ID对应的行
  overexpression_final <- overexpression_filtered %>% 
    filter(!Ensembl_ID %in% deletion_ids)
  
  # 更新deletion_filtered：去除Category为overexpression的Ensembl_ID对应的行
  deletion_final <- deletion_filtered %>% 
    filter(!Ensembl_ID %in% overexpression_ids)
  
  # 查看更新后的数据集
  head(deletion_final)
  head(overexpression_final)
  
  merged_data <- rbind(deletion_final, overexpression_final)
}

head(merged_data)
names(merged_data) <- paste0(names(merged_data), "_CM")

saveRDS(merged_data, file = "~/AF_atlas/Data/in silico genes/combined_data_CMs.RDS")

combined_data_CMs <- readRDS("~/AF_atlas/Data/in silico genes/combined_data_CMs.RDS")
head(combined_data_CMs)
ven_combine_unique <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/ven_combine_unique.RDS")
head(ven_combine_unique)

inter_CMs_genes<-intersect(ven_combine_unique$exposure,unique(combined_data_CMs$Ensembl_ID_CM))

inter_silico_to_CM_eQTL<-ven_combine_unique %>% filter(exposure %in% inter_CMs_genes)
inter_silico_to_CM_eQTL

combined_data_CMs_subset<-filter(combined_data_CMs,Ensembl_ID_CM %in% inter_CMs_genes)
combined_data_CMs_subset

#----------------------------------------------心肌细胞silico基因曲线图展示----------------------------------------
deletion <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/CM/CM.csv", header = TRUE, sep = ",")
overexpression <- read.csv("~/AF_atlas/Data/in silico genes/每种细胞/CM/CM_over.csv", header = TRUE, sep = ",")

{
  deletion_filtered <- filter(deletion, Sig == 1)
  overexpression_filtered <- filter(overexpression, Sig == 1)
  deletion_filtered$Category <- 'deletion'
  overexpression_filtered$Category <- 'overexpression'
  common_genes <- intersect(deletion_filtered$Ensembl_ID, overexpression_filtered$Ensembl_ID)
  length(unique(common_genes))
  
  deletion_subset <- deletion_filtered[deletion_filtered$Ensembl_ID %in% common_genes, ]
  overexpression_subset <- overexpression_filtered[overexpression_filtered$Ensembl_ID%in% common_genes, ]
  
  head(deletion_subset)
  head(overexpression_subset)
  
  common_ids <- merge(deletion_subset, overexpression_subset, by = "Ensembl_ID")
  
  common_ids <- common_ids[sign(common_ids$Shift_to_goal_end.x) == sign(common_ids$Shift_to_goal_end.y), ]
  
  common_ids
  dim(common_ids)
  
  common_ids$source <- ifelse(abs(common_ids$Shift_to_goal_end.x) > abs(common_ids$Shift_to_goal_end.y), 
                              "deletion", 
                              "overexpression")
  
  # 2. 创建最终的子集
  final_subset <- rbind(
    deletion_subset[deletion_subset$Ensembl_ID %in% common_ids$Ensembl_ID[common_ids$source == "deletion"], ],
    overexpression_subset[overexpression_subset$Ensembl_ID %in% common_ids$Ensembl_ID[common_ids$source == "overexpression"], ]
  )
  
  # 查看结果
  head(final_subset)
  dim(final_subset)
  final_subset
  
  # 找到Category为deletion的Ensembl_ID
  deletion_ids <- final_subset %>% 
    filter(Category == "deletion") %>% 
    .$Ensembl_ID
  
  # 找到Category为overexpression的Ensembl_ID
  overexpression_ids <- final_subset %>% 
    filter(Category == "overexpression") %>% 
    .$Ensembl_ID
  
  # 更新overexpression_filtered：去除Category为deletion的Ensembl_ID对应的行
  overexpression_final <- overexpression_filtered %>% 
    filter(!Ensembl_ID %in% deletion_ids)
  
  # 更新deletion_filtered：去除Category为overexpression的Ensembl_ID对应的行
  deletion_final <- deletion_filtered %>% 
    filter(!Ensembl_ID %in% overexpression_ids)
  
  # 查看更新后的数据集
  head(deletion_final)
  head(overexpression_final)
  
  merged_data <- rbind(deletion_final, overexpression_final)
  head(merged_data)
  dim(merged_data)
  
  # 排序基因，以便绘图时有序
  merged_data <- merged_data %>%
    arrange(desc(Shift_to_goal_end))
  head(merged_data)
  
  # 为绘图生成行号（排名）
  merged_data$Rank <- seq_along(merged_data$Shift_to_goal_end)
  
  # 选择颜色方案
  color_palette <- brewer.pal(n = 3, name = "Set1")
  #color_palette
  #color_palette <- c("#377EB8" ,"#E41A1C","#4DAF4A" )
  
  # 绘制图形
  ggplot(merged_data, aes(x = Rank, y = Shift_to_goal_end, group = Category, color = Category)) +
    geom_line() +
    geom_point(aes(color = Category), size = 2) +
    scale_color_manual(values = color_palette) +
    theme_minimal() +
    labs(x = "Rank", y = "Shift to normal", title = "In silico treatment") +
    geom_text_repel(data = merged_data %>% 
                      tidyr::drop_na() %>% 
                      dplyr::arrange(desc(Shift_to_goal_end)) %>%
                      dplyr::slice(1:5),
                    aes(x = Rank, y = Shift_to_goal_end, label = Gene_name),
                    box.padding = 0.3,
                    nudge_x = 400,
                    nudge_y = 0.0001,
                    segment.curvature = -0.1,
                    segment.ncp = 3,
                    direction = "y", 
                    hjust = "left",
                    max.overlaps = Inf
    )+
    geom_text_repel(data = merged_data %>% 
                      tidyr::drop_na() %>% 
                      dplyr::arrange(desc(-Shift_to_goal_end)) %>%
                      dplyr::slice(1:5),
                    aes(x = Rank, y = Shift_to_goal_end, label = Gene_name),
                    box.padding = 0.3,
                    nudge_x = -500,
                    nudge_y = 0.0001,
                    segment.curvature = 0.1,
                    segment.ncp = 3,
                    direction = "y", 
                    hjust = "left",
                    max.overlaps = Inf
    )
}

overexpression_subset<-filter(merged_data,Category=="overexpression")
dim(overexpression_subset)
length(unique(overexpression_subset$Ensembl_ID)) # 1033

deletion_subset<-filter(merged_data,Category=="deletion")
dim(deletion_subset)
length(unique(deletion_subset$Ensembl_ID)) # 134

#----------------------------------------------心肌细胞silico基因富集分析----------------------------------------
combined_data_CMs <- readRDS("~/AF_atlas/Data/in silico genes/combined_data_CMs.RDS")
head(combined_data_CMs)

CM_deletion<-filter(combined_data_CMs,Category_CM=="deletion")
CM_overexpression<-filter(combined_data_CMs,Category_CM=="overexpression")

keys <- mapIds(org.Hs.eg.db,
               keys = CM_deletion$Ensembl_ID_CM,
               column = "ENTREZID",
               keytype = "ENSEMBL",
               multiVals = "first")

go <- enrichGO(
  gene = keys,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.5,
  readable=F) #是否将gene ID转换为gene symbol

dotplot(go, split="ONTOLOGY",showCategory = 5,label_format=50)+ 
  facet_grid(ONTOLOGY~.,scale="free")+
  theme(panel.grid = element_blank())+#修改主题
  theme(axis.title =element_text(size = 12, color = 'black'),
        axis.text.y =element_text(size = 12),
        legend.title=element_text(size=12))+
  scale_color_gradient(high="#FC8D62",low="#4b5cc4")#设置颜色

keys <- mapIds(org.Hs.eg.db,
               keys = CM_overexpression$Ensembl_ID_CM,
               column = "ENTREZID",
               keytype = "ENSEMBL",
               multiVals = "first")

go <- enrichGO(
  gene = keys,
  OrgDb = org.Hs.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.5,
  readable=F) #是否将gene ID转换为gene symbol

dotplot(go, split="ONTOLOGY",showCategory = 5,label_format=50)+ 
  facet_grid(ONTOLOGY~.,scale="free")+
  theme(panel.grid = element_blank())+#修改主题
  theme(axis.title =element_text(size = 12, color = 'black'),
        axis.text.y =element_text(size = 12),
        legend.title=element_text(size=12))+
  scale_color_gradient(high="#FC8D62",low="#4b5cc4")#设置颜色


#----------------------------------------------所有基因upset分析----------------------------------------
library(UpSetR)

ven_combine_unique <- readRDS("~/AF_atlas/Data/孟德尔随机化分析数据/ven_combine_unique.RDS")
head(ven_combine_unique)
table(ven_combine_unique$source,ven_combine_unique$category)

ven_subset <- ven_combine_unique %>%
  mutate(
    Atrial_negative = as.integer(source == "Atrial" & category == "negative"),
    Atrial_positive = as.integer(source == "Atrial" & category == "positive"),
    Blood_negative = as.integer(source == "Blood" & category == "negative"),
    Blood_positive = as.integer(source == "Blood" & category == "positive"),
    Ventricle_negative = as.integer(source == "Ventricle" & category == "negative"),
    Ventricle_positive = as.integer(source == "Ventricle" & category == "positive")
  ) %>%
  select(exposure, Atrial_negative, Atrial_positive, Blood_negative, Blood_positive, Ventricle_negative, Ventricle_positive)

head(ven_subset)

ven_subset_combined <- ven_subset %>%
  group_by(exposure) %>%
  summarise(
    Atrial_negative = max(Atrial_negative),
    Atrial_positive = max(Atrial_positive),
    Blood_negative = max(Blood_negative),
    Blood_positive = max(Blood_positive),
    Ventricle_negative = max(Ventricle_negative),
    Ventricle_positive = max(Ventricle_positive)
  )
dim(ven_subset_combined)

ven_subset_combined_for_upset<-data.frame(ven_subset_combined)

any(duplicated(ven_subset_combined_for_upset$exposure))

upset(ven_subset_combined_for_upset,nsets = 100,nintersects = 100,
      shade.alpha = 0.1,
      mb.ratio = c(0.618, 0.382),
      matrix.color = '#225EA8',
      set_size.show = T,
      show.numbers = F)


combined_data_所有疾病靶点基因 <- readRDS("~/AF_atlas/Data/in silico genes/combined_data_所有疾病靶点基因.RDS")
head(combined_data_所有疾病靶点基因)

process_shift_col <- function(shift_col_name, data) {
  cat_col_name <- paste0("Category_", gsub("Shift_to_goal_end_", "", shift_col_name))
  
  towards_SR_deletion_col <- paste0(shift_col_name, "_towards_SR_deletion")
  towards_SR_overexpression_col <- paste0(shift_col_name, "_towards_SR_overexpression")
  away_SR_deletion_col <- paste0(shift_col_name, "_away_SR_deletion")
  away_SR_overexpression_col <- paste0(shift_col_name, "_away_SR_overexpression")
  
  data %>%
    mutate(
      !!towards_SR_deletion_col := as.integer((!!sym(cat_col_name) == "deletion") & (!is.na(!!sym(shift_col_name))) & (!!sym(shift_col_name) > 0)),
      !!towards_SR_overexpression_col := as.integer((!!sym(cat_col_name) == "overexpression") & (!is.na(!!sym(shift_col_name))) & (!!sym(shift_col_name) > 0)),
      !!away_SR_deletion_col := as.integer((!!sym(cat_col_name) == "deletion") & (!is.na(!!sym(shift_col_name))) & (!!sym(shift_col_name) < 0)),
      !!away_SR_overexpression_col := as.integer((!!sym(cat_col_name) == "overexpression") & (!is.na(!!sym(shift_col_name))) & (!!sym(shift_col_name) < 0))
    ) %>%
    mutate(across(starts_with(shift_col_name), ~ replace_na(., 0)))
}

# 应用处理函数到所有Shift_to_goal_end列
shift_cols <- grep("Shift_to_goal_end", names(combined_data_所有疾病靶点基因), value = TRUE)
result_data <- combined_data_所有疾病靶点基因

for (col in shift_cols) {
  result_data <- process_shift_col(col, result_data)
}

# 查看结果
head(result_data)

# 选择第一列和最后32列
selected_columns <- c(1, (ncol(result_data)-31):ncol(result_data))
new_dataset <- result_data[, selected_columns]

# 查看新数据集
head(new_dataset)
any(duplicated(new_dataset$Ensembl_ID_T_cell))

new_dataset_combined <- new_dataset %>%
  group_by(Ensembl_ID_T_cell) %>%
  summarise(across(everything(), ~ max(., na.rm = TRUE)))

new_dataset_combined_for_upset<-data.frame(new_dataset_combined)

any(duplicated(new_dataset_combined_for_upset$Ensembl_ID_T_cell))

upset_data<-merge(new_dataset_combined_for_upset,ven_subset_combined_for_upset,
                              by.x="Ensembl_ID_T_cell",by.y="exposure",all=T)
head(upset_data)
any(duplicated(upset_data$Ensembl_ID_T_cell))

#数据验证
upset_data[is.na(upset_data)] <- 0
head(upset_data)
table(upset_data$Atrial_negative)
table(upset_data$Atrial_positive)
table(upset_data$Ventricle_negative)
table(upset_data$Ventricle_positive)
table(upset_data$Blood_negative)
table(upset_data$Blood_positive)
table(upset_data$Shift_to_goal_end_T_cell_towards_SR_deletion)
table(upset_data$Shift_to_goal_end_T_cell_towards_SR_overexpression)
table(upset_data$Shift_to_goal_end_T_cell_away_SR_deletion)
table(upset_data$Shift_to_goal_end_T_cell_away_SR_overexpression)

saveRDS(upset_data, file = "~/AF_atlas/Data/in silico genes/upset_data用于画图.RDS")

combined_data_所有疾病靶点基因 <- readRDS("~/AF_atlas/Data/in silico genes/combined_data_所有疾病靶点基因.RDS")
table(combined_data_所有疾病靶点基因$Category_T_cell)

#+++++++++++++++++++++++++++++++++++++++++++ 画图+++++++++++++++++++++++++++++
upset_data <- readRDS("~/AF_atlas/Data/in silico genes/upset_data用于画图.RDS")

selected_columns <- c(grep("T_cell", names(upset_data), value = TRUE),"Ensembl_ID_T_cell",
                      "Atrial_negative", "Atrial_positive", "Blood_negative",
                      "Blood_positive", "Ventricle_negative", "Ventricle_positive")

# 从upset_data中选择这些列
upset_subset <- upset_data %>% select(all_of(selected_columns))

extended_pal_30 <- rep(pal_30, length.out = 10)
extended_pal_300 <- rep(pal_30, length.out = 17)

upset(upset_subset,nsets = 100,nintersects = 100,
      shade.alpha = 0.1,
      main.bar.color = extended_pal_300, 
      sets.bar.color = extended_pal_30, 
      mb.ratio = c(0.618, 0.382),
      matrix.color = '#225EA8')

inter_matrix<-upset_subset[,c("Ensembl_ID_T_cell","Ventricle_negative","Shift_to_goal_end_T_cell_towards_SR_overexpression")]
row_sums <- rowSums(inter_matrix[,-1])
selected_rows <- inter_matrix[row_sums > 1, ]
selected_rows

inter_matrix<-upset_subset[,c("Ensembl_ID_T_cell","Atrial_positive","Shift_to_goal_end_T_cell_towards_SR_overexpression")]
row_sums <- rowSums(inter_matrix[,-1])
selected_rows <- inter_matrix[row_sums > 1, ]
selected_rows

#------------------------------------------------
selected_columns <- c("Ensembl_ID_T_cell",grep("FB", names(upset_data), value = TRUE),
                      "Atrial_negative", "Atrial_positive", "Blood_negative",
                      "Blood_positive", "Ventricle_negative", "Ventricle_positive")

# 从upset_data中选择这些列
upset_subset <- upset_data %>% select(all_of(selected_columns))

extended_pal_30 <- rep(pal_30, length.out = 10)
extended_pal_300 <- rep(pal_30, length.out = 20)

upset(upset_subset,nsets = 100,nintersects = 100,
      shade.alpha = 0.1,
      main.bar.color = extended_pal_300, 
      sets.bar.color = extended_pal_30, 
      mb.ratio = c(0.618, 0.382),
      matrix.color = '#225EA8')

inter_matrix<-upset_subset[,c("Ensembl_ID_T_cell","Atrial_negative","Shift_to_goal_end_FB_away_SR_overexpression")]
row_sums <- rowSums(inter_matrix[,-1])
selected_rows <- inter_matrix[row_sums > 1, ]
selected_rows

inter_matrix<-upset_subset[,c("Ensembl_ID_T_cell","Blood_positive","Shift_to_goal_end_FB_away_SR_overexpression")]
row_sums <- rowSums(inter_matrix[,-1])
selected_rows <- inter_matrix[row_sums > 1, ]
selected_rows

inter_matrix<-upset_subset[,c("Ensembl_ID_T_cell","Blood_positive","Shift_to_goal_end_FB_towards_SR_overexpression")]
row_sums <- rowSums(inter_matrix[,-1])
selected_rows <- inter_matrix[row_sums > 1, ]
selected_rows

inter_matrix<-upset_subset[,c("Ensembl_ID_T_cell","Ventricle_positive","Shift_to_goal_end_FB_towards_SR_overexpression")]
row_sums <- rowSums(inter_matrix[,-1])
selected_rows <- inter_matrix[row_sums > 1, ]
selected_rows

inter_matrix<-upset_subset[,c("Ensembl_ID_T_cell","Ventricle_positive","Atrial_positive","Shift_to_goal_end_FB_towards_SR_overexpression")]
row_sums <- rowSums(inter_matrix[,-1])
selected_rows <- inter_matrix[row_sums > 1, ]
selected_rows

#------------------------------------------------
selected_columns <- c("Ensembl_ID_T_cell",grep("EC", names(upset_data), value = TRUE),
                      "Atrial_negative", "Atrial_positive", "Blood_negative",
                      "Blood_positive", "Ventricle_negative", "Ventricle_positive")

# 从upset_data中选择这些列
upset_subset <- upset_data %>% select(all_of(selected_columns))

extended_pal_30 <- rep(pal_30, length.out = 10)
extended_pal_300 <- rep(pal_30, length.out = 18)

upset(upset_subset,nsets = 100,nintersects = 100,
      shade.alpha = 0.1,
      main.bar.color = extended_pal_300, 
      sets.bar.color = extended_pal_30, 
      mb.ratio = c(0.618, 0.382),
      matrix.color = '#225EA8')

inter_matrix<-upset_subset[,c("Ensembl_ID_T_cell","Ventricle_negative","Shift_to_goal_end_EC_away_SR_overexpression")]
row_sums <- rowSums(inter_matrix[,-1])
selected_rows <- inter_matrix[row_sums > 1, ]
selected_rows

inter_matrix<-upset_subset[,c("Ensembl_ID_T_cell","Atrial_negative","Shift_to_goal_end_EC_away_SR_overexpression")]
row_sums <- rowSums(inter_matrix[,-1])
selected_rows <- inter_matrix[row_sums > 1, ]
selected_rows

inter_matrix<-upset_subset[,c("Ensembl_ID_T_cell","Ventricle_positive","Atrial_positive","Shift_to_goal_end_EC_away_SR_overexpression")]
row_sums <- rowSums(inter_matrix[,-1])
selected_rows <- inter_matrix[row_sums > 1, ]
selected_rows


#------------------------------------------------
selected_columns <- c("Ensembl_ID_T_cell",grep("MP", names(upset_data), value = TRUE),
                      "Atrial_negative", "Atrial_positive", "Blood_negative",
                      "Blood_positive", "Ventricle_negative", "Ventricle_positive")

# 从upset_data中选择这些列
upset_subset <- upset_data %>% select(all_of(selected_columns))

extended_pal_30 <- rep(pal_30, length.out = 10)
extended_pal_300 <- rep(pal_30, length.out = 16)

upset(upset_subset,nsets = 100,nintersects = 100,
      shade.alpha = 0.1,
      main.bar.color = extended_pal_300, 
      sets.bar.color = extended_pal_30, 
      mb.ratio = c(0.618, 0.382),
      matrix.color = '#225EA8')


inter_matrix<-upset_subset[,c("Ensembl_ID_T_cell","Ventricle_negative","Atrial_negative","Shift_to_goal_end_MP_towards_SR_overexpression")]
row_sums <- rowSums(inter_matrix[,-1])
selected_rows <- inter_matrix[row_sums > 1, ]
selected_rows



#------------------------------------------------
selected_columns <- c("Ensembl_ID_T_cell",grep("DC", names(upset_data), value = TRUE),
                      "Atrial_negative", "Atrial_positive", "Blood_negative",
                      "Blood_positive", "Ventricle_negative", "Ventricle_positive")

# 从upset_data中选择这些列
upset_subset <- upset_data %>% select(all_of(selected_columns))

extended_pal_30 <- rep(pal_30, length.out = 10)
extended_pal_300 <- rep(pal_30, length.out = 16)

upset(upset_subset,nsets = 100,nintersects = 100,
      shade.alpha = 0.1,
      main.bar.color = extended_pal_300, 
      sets.bar.color = extended_pal_30, 
      mb.ratio = c(0.618, 0.382),
      matrix.color = '#225EA8')

inter_matrix<-upset_subset[,c("Ensembl_ID_T_cell","Atrial_positive","Shift_to_goal_end_DC_towards_SR_overexpression")]
row_sums <- rowSums(inter_matrix[,-1])
selected_rows <- inter_matrix[row_sums > 1, ]
selected_rows



#------------------------------------------------
selected_columns <- c("Ensembl_ID_T_cell",grep("SMC", names(upset_data), value = TRUE),
                      "Atrial_negative", "Atrial_positive", "Blood_negative",
                      "Blood_positive", "Ventricle_negative", "Ventricle_positive")

# 从upset_data中选择这些列
upset_subset <- upset_data %>% select(all_of(selected_columns))

extended_pal_30 <- rep(pal_30, length.out = 10)
extended_pal_300 <- rep(pal_30, length.out = 16)

upset(upset_subset,nsets = 100,nintersects = 100,
      shade.alpha = 0.1,
      main.bar.color = extended_pal_300, 
      sets.bar.color = extended_pal_30, 
      mb.ratio = c(0.618, 0.382),
      matrix.color = '#225EA8')


inter_matrix<-upset_subset[,c("Ensembl_ID_T_cell","Atrial_negative","Shift_to_goal_end_SMC_towards_SR_deletion")]
row_sums <- rowSums(inter_matrix[,-1])
selected_rows <- inter_matrix[row_sums > 1, ]
selected_rows



#------------------------------------------------
selected_columns <- c("Ensembl_ID_T_cell",grep("Neutrophil", names(upset_data), value = TRUE),
                      "Atrial_negative", "Atrial_positive", "Blood_negative",
                      "Blood_positive", "Ventricle_negative", "Ventricle_positive")

# 从upset_data中选择这些列
upset_subset <- upset_data %>% select(all_of(selected_columns))

extended_pal_30 <- rep(pal_30, length.out = 10)
extended_pal_300 <- rep(pal_30, length.out = 15)

upset(upset_subset,nsets = 100,nintersects = 100,
      shade.alpha = 0.1,
      main.bar.color = extended_pal_300, 
      sets.bar.color = extended_pal_30, 
      mb.ratio = c(0.618, 0.382),
      matrix.color = '#225EA8')





#------------------------------------------------
selected_columns <- c("Ensembl_ID_T_cell",grep("CM", names(upset_data), value = TRUE),
                      "Atrial_negative", "Atrial_positive", "Blood_negative",
                      "Blood_positive", "Ventricle_negative", "Ventricle_positive")

# 从upset_data中选择这些列
upset_subset <- upset_data %>% select(all_of(selected_columns))

extended_pal_30 <- rep(pal_30, length.out = 10)
extended_pal_300 <- rep(pal_30, length.out = 20)

upset(upset_subset,nsets = 100,nintersects = 100,
      shade.alpha = 0.1,
      main.bar.color = extended_pal_300, 
      sets.bar.color = extended_pal_30, 
      mb.ratio = c(0.618, 0.382),
      matrix.color = '#225EA8')


inter_matrix<-upset_subset[,c("Ensembl_ID_T_cell","Atrial_positive","Shift_to_goal_end_CM_towards_SR_overexpression")]
row_sums <- rowSums(inter_matrix[,-1])
selected_rows <- inter_matrix[row_sums > 1, ]
selected_rows

inter_matrix<-upset_subset[,c("Ensembl_ID_T_cell","Ventricle_positive","Shift_to_goal_end_CM_away_SR_overexpression")]
row_sums <- rowSums(inter_matrix[,-1])
selected_rows <- inter_matrix[row_sums > 1, ]
selected_rows

inter_matrix<-upset_subset[,c("Ensembl_ID_T_cell","Ventricle_negative","Shift_to_goal_end_CM_towards_SR_overexpression")]
row_sums <- rowSums(inter_matrix[,-1])
selected_rows <- inter_matrix[row_sums > 1, ]
selected_rows

inter_matrix<-upset_subset[,c("Ensembl_ID_T_cell","Atrial_negative","Shift_to_goal_end_CM_towards_SR_overexpression")]
row_sums <- rowSums(inter_matrix[,-1])
selected_rows <- inter_matrix[row_sums > 1, ]
selected_rows

inter_matrix<-upset_subset[,c("Ensembl_ID_T_cell","Ventricle_negative","Atrial_negative","Shift_to_goal_end_CM_towards_SR_deletion")]
row_sums <- rowSums(inter_matrix[,-1])
selected_rows <- inter_matrix[row_sums > 1, ]
selected_rows


























