library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(ComplexHeatmap)
library(Seurat)
library(SingleR)

#read data
fn <- c(list.files("./", pattern = "*.h5", full.names = T))

#create seurat object
data_gem2 <- list()
for (i in 1:3) {
  data_gem2[[i]] <- Read10X_h5(fn[i])
  data_gem2[[i]] <- CreateSeuratObject(counts = data_gem2[[i]], 
                                      project = file_tab$name[i],
                                      min.cells = 10)
  data_gem2[[i]]$percent.mt <- PercentageFeatureSet(data_gem2[[i]], pattern = "^MT-")
}

#merge data
data_gem2 <- merge(data_gem2[[1]], y = c(data_gem2[[2]], data_gem2[[3]]), data_gem2[[4]]),
                  add.cell.ids = file_tab$name)

#quality control
data_gem2$cell_id <- names(data_gem2$orig.ident)
data_gem_sub <- subset(data_gem2, cell_id %in% data_vdj_single_tab$cell_id)
VlnPlot(data_gem_sub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(data_gem2, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(data_gem2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
data_gem_sub <- subset(data_gem_sub, nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)

#integration
data_gem_sub2 <- SplitObject(data_gem_sub, split.by = "ident")
data_gem_sub2 <- lapply(data_gem_sub2, function(x){
  x <- NormalizeData(x) #normalization
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000) #identify variable features
})

features <- SelectIntegrationFeatures(data_gem_sub2)
data_anchors <- FindIntegrationAnchors(data_gem_sub2, anchor.features = features)
data_gem_comb2 <- IntegrateData(data_anchors)
DefaultAssay(data_gem_comb2) <- "integrated"

#scale
data_gem_comb2 <- ScaleData(data_gem_comb2, verbose = F)

#PCA & cluster
data_gem_comb2 <- RunPCA(data_gem_comb2, npcs = 30, verbose = F)
data_gem_comb2 <- RunUMAP(data_gem_comb2, reduction = "pca", dims = 1:30)
data_gem_comb2 <- FindNeighbors(data_gem_comb2, reduction = "pca", dims = 1:30)
data_gem_comb2 <- FindClusters(data_gem_comb2, resolution = 1)
data_gem_comb2 <- RunTSNE(data_gem_comb2, reduction = "pca", dims = 1:30)

DimPlot(data_gem_comb2, reduction = "umap", label = F, label.size = 4, pt.size = .5) +
  scale_color_discrete(labels = paste0("C", 0:8)) + 
  labs(title = "") +
  NoLegend()

cluster_marker <- FindAllMarkers(data_gem_comb2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#cell type annotation
monaco_ref <- MonacoImmuneData()
data_gem_comb_single2 <- as.SingleCellExperiment(data_gem_comb2)
pred_gem_main <- SingleR(test = data_gem_comb_single2, ref = monaco_ref, assay.type.test = 1,
                           labels = monaco_ref$label.main)
pred_gem_fine <- SingleR(test = data_gem_comb_single2, ref = monaco_ref, assay.type.test = 1,
                           labels = monaco_ref$label.fine, method = "wilcox")
pred_gem <- merge(data.frame(cell_id = rownames(pred_gem_main), main_label = pred_gem_main$pruned.labels, stringsAsFactors = F),
                   data.frame(cell_id = rownames(pred_gem_fine), fine_label = pred_gem_fine$pruned.labels, stringsAsFactors = F), 
                   all = T, by = c("cell_id"), sort = F)
pred_gem$label <- pred_gem$main_label
pred_gem$label[which(pred_gem$label == "B cells")] <- pred_gem$fine_label[which(pred_gem$main_label == "B cells")]

cell_type_tab <- data.frame(main = monaco_ref$label.main, 
                            fine = monaco_ref$label.fine,
                            stringsAsFactors = F)
cell_type_tab$filter <- paste(cell_type_tab$main, cell_type_tab$fine, sep = "/")

#remove unmatched main label and fine label
pred_gem$filter <- paste(pred_gem$main_label, pred_gem$fine_label, sep = "/")
pos <- which(pred_gem$filter %in% cell_type_tab$filter)
pred_gem$label[-pos] <- "Unclear"
pred_gem <- pred_gem[, -5]

#add singleR label to seurat object
data_gem_comb2$singler_label <- multi_gsub(data_gem_comb2$seurat_clusters, 
                                          rownames(pred_gem_fine), pred_gem_fine$labels)
										  
#modify label
data_gem_comb2$isotype <- multi_gsub(names(data_gem_comb2$orig.ident), 
                                    data_vdj_labeled$cell_id, data_vdj_labeled$c_gene_H)
data_gem_comb2$isotype[which(data_gem_comb2$isotype == "")] <- "Unk"
data_gem_comb2$isotype <- factor(data_gem_comb2$isotype, 
                                levels = c("IGHD", "IGHM", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2", "Unk"))
		
data_gem_comb2$label_modf <- data_gem_comb2$singler_label
data_gem_comb2$label_modf[which(data_gem_comb2$isotype %in% isotype[3:8] & data_gem_comb2$singler_label == "Non-switched memory B cells")] <- "Switched memory B cells"
data_gem_comb2$label_modf[which(data_gem_comb2$isotype %in% isotype[1:2] & data_gem_comb2$singler_label == "Switched memory B cells")] <- "Non-switched memory B cells"
data_gem_comb2$label_modf <- factor(data_gem_comb2$label_modf, 
                                   levels = unique(data_gem_comb2$label_modf))	

DimPlot(data_gem_comb2, reduction = "umap", label = F, label.size = 4, pt.size = .5,
        group.by = "label_modf") +
  scale_color_discrete(labels = c("Naïve B cells", "Non-switched MBCs",
                                  "Switched MBCs", "Plasmablasts")) +
  guides(color = "none") +
  labs(title = "")		

#selected & validated mAbs
h <- data_vdj_single_tab
h$aaseq_h <- paste0(h$fwr1_H, h$cdr1_H, h$fwr2_H, h$cdr2_H, h$fwr3_H, h$cdr3_H, h$fwr4_H)
h$aaseq_l <- paste0(h$fwr1_L, h$cdr1_L, h$fwr2_L, h$cdr2_H, h$fwr3_L, h$cdr3_L, h$fwr4_L)
h <- h[, c(44,3,5,45,17,50,24,26,46,38,51)]
h$cdr3len_h <- nchar(h$cdr3_H) - 2
names(h) <- names(ab_selc)[-c(5,11,14,15,17)]
nt <- list()
for (i in 1:nrow(ab_selc)) {
  nt[[i]] <- h[which(h$HCDR3 == ab_selc$HCDR3[i] & h$LCDR3 == ab_selc$LCDR3[i]),]
  if(nrow(nt[[i]]) > 0){
    nt[[i]] <- rbind(ab_selc[i, -c(5,11,14,15,17)], nt[[i]])
    nt[[i]]$cl <- paste("cl.", i)
  }
}
df <- do.call(rbind, nt)
df$group <- "sc"
df$group[which(df$mAb %in% ab_selc$mAb[which(ab_selc$group == "TRUE")])] <- "Delta binding"
df$group[which(df$mAb %in% ab_bind_tab$mAb)] <- "Cross-variant\nbinding"
n <- data_gem_comb2
n$validated <- "All"
n$validated[which(n$integrated_snn_res.1 %in% c(6))] <- "C6"
m <- df[,c(1,13,14)]
t <- unique(m[which(m$mAb %in% ab_selc$mAb), 2:3])
m$validated <- multi_gsub(m$cl, t$cl, t$group)
m <- m[which(m$mAb %in% data_vdj_single_tab$cell_id),]
n$validated <- factor(n$validated, levels = unique(n$validated)[c(4,2,3,1)])
n@active.ident <- factor(n$validated, levels = unique(n$validated)[c(1,3,2,4)])

DimPlot(n, reduction = "umap", label = F, pt.size = .6, 
        order = rev(levels(Idents(n))),
        group.by = "validated", cols = rev(c("red", "grey60", "#FDCDAC", "grey90"))) +
  labs(title = "", x = "UMAP 1", y = "UMAP 2") +
  NoLegend()  
  
m <- n@meta.data
m <- m[which(m$validated %in% c("Delta binding", "Cross-variant\nbinding")),]
m$label <- NA
m$label[which(m$cell_id == "Y9_TGACTAGGTGTCAATC-1")] <- "YB9-120"
m$label[which(m$cell_id == "Y12_CCATGTCTCGGCCGAT-1")] <- "YB12-197"
m$label[which(m$cell_id == "Y13_AGCGGTCCACGTGAGA-1")] <- "YB13-208"
m$label[which(m$cell_id == "Y9_TTATGCTGTAAATGAC-1")] <- "YB9-258"
m$label[which(m$cell_id == "Y13_GAACGGAGTCTCGTTC-1")] <- "YB13-292"
m$col <- "black"
m$col[which(is.na(m$label))] <- "white"
m$integrated_snn_res.0.75 <- paste0("C", m$integrated_snn_res.0.75)
m$validated <- factor(m$validated, levels = unique(m$validated))
m <- arrange(m, desc(col))

ggplot(m, aes(x = integrated_snn_res.0.75, y = validated, label = label)) +
  geom_jitter(aes(fill = validated, color = col), width = 0.25, height = 0.25, 
              shape = 21, size = 2) +
  scale_fill_manual(values = rev(c("red", "grey60"))) +
  scale_color_manual(values = c("black", "white")) +
  scale_y_discrete(labels = c("Delta binding", "Cross-variant\nbinding")) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "grey", linetype = "dashed", size = 0.5),
        panel.spacing = unit(.3, "lines"),
        panel.border = element_rect(colour = "black", size = 0.5, fill = NA),
        axis.line.x.bottom = element_line(colour = "black", size = 0.5),
        axis.line.y.left = element_line(colour = "black", size = 0.5),
        axis.title.x.bottom = element_text(size = 12, family = "Arial", face = "bold"),
        axis.text.y.left = element_text(size = 12, family = "Arial", color = "black", hjust = 0.5),
        axis.text.x.bottom = element_text(size = 10, family = "Arial", hjust = 0.5, color = "black")) +
  labs(x = "Cluster", y = NULL) +
  guides(color = "none", fill = "none")

#DEG
marker_signif <- subset(cluster_marker, p_val_adj < 0.05)
aver_exp_cl <- AverageExpression(data_gem_comb2)
m <- aver_exp_cl$integrated %>% as.data.frame()
m <- apply(m, 1, function(x){z.score(x)}) %>% t()

nt <- list()
for (i in 1:10) {
  nt[[i]] <- head(marker_signif[which(marker_signif$cluster == (i-1)),6:7], 5)
}
t <- do.call(rbind, nt)
mat <- m[which(rownames(m) %in% t$gene),]
colnames(mat) <- paste0("C", 0:8)
mat <- t(mat)

col_fun <- colorRamp2(c(-2, 0, 2), c("#1B9E77", "white", "#D95F02"))
Heatmap(mat, show_row_dend = F, show_column_dend = F, show_heatmap_legend = T, 
        show_row_names = T, show_column_names = T,
        col = col_fun, 
        column_order = unique(t$gene), 
        row_order = rownames(mat),
        row_names_side = "left", 
        column_names_side = "bottom", column_names_rot = 45, column_names_centered = F,
        border = "black",
        row_names_gp = gpar(family = "Arial"),
        column_names_gp = gpar(family = "Arial"),
        heatmap_legend_param = list(title = "Average\nexpression",
                                    family = "Arial", border = "black"))
