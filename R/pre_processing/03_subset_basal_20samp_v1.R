# geo ----
setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_basal_20samp_v1";dir.create(odir)
library(future)
plan("multiprocess", workers = 20)
options(future.globals.maxSize = Inf)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "03_subset_basal_v2"
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_basal_20samp_v1/"

# core ----
setwd("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_basal_20samp_v1"; dir.create(odir)
library(future)
plan("multiprocess", workers = 50)
options(future.globals.maxSize = Inf)
setwd(paste0("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "03_subset_basal_v2"
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_basal_20samp_v1/"
# stress_sig = read.delim("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_fibroblasts_final1/stromal_stress_signature_032720.txt")
# stress_sig = as.character(stress_sig$stromal.stress.sig)


# read original object --------------
tum = read_rds(paste0(wd, "03_subset_basal_v2/dim_12k_20_tum.integrated.ns_integrated_v4_intgeratedassay_integrated_snn_res.0.3_vst_method_9.rds"))
tum$sheet_patient_id[tum$sheet_ts_sampleid %in% c("hbca25_c")] = "pt37"

png(paste0(save_path, "vlnplot_nCount_RNA_sample_id.png"), width=14, height=3, units = "in", res = 300)
VlnPlot(tum, features = "nCount_RNA", group.by = "final_sample_id", pt.size = 0 )
dev.off()

# subset and filter ----------
Idents(tum) <- "batch"
tum_sub = subset(tum, idents = c("b2"), invert = F) #endo_vas
Idents(tum_sub) <- "sheet_exp_proc"
tum_sub = subset(tum_sub, idents = c("3hr_digestion"), invert = F) #endo_vas
Idents(tum_sub) <- "sheet_figure1_grp_updated"
tum_sub = subset(tum_sub, idents = c("Reduction Mammoplasty"), invert = F) #endo_vas

write_rds(tum_sub, "epi_basal_sub_3hr_b2_reduction_v1.rds")
tum_sub = read_rds("epi_basal_sub_3hr_b2_reduction_v1.rds")

# *** recluster only ---------
seu_obj = tum_sub
#seu_obj = NormalizeData(seu_object) %>% FindVariableFeatures(nfeatures = 5000) 
seu_obj <- SCTransform(seu_obj)
all.genes <- rownames(seu_obj)
#seu_obj <- ScaleData(seu_obj, features = all.genes)
seu_obj <- RunPCA(seu_obj, features = VariableFeatures(object = seu_obj))
seu_obj <- FindNeighbors(seu_obj, dims = 1:18)
seu_obj <- FindClusters(seu_obj, resolution = c(0.5, 0.1, 0.05, 0.075, 0.2, 0.15, 0.3))
seu_obj <- RunUMAP(seu_obj, dims = 1:18)
write_rds(seu_obj, "seu_obj_noscaling.rds")
write_rds(seu_obj, "seu_obj_SCTscaling.rds")

seu_obj = read_rds("seu_obj_noscaling.rds")
DimPlot(seu_obj, reduction = "umap", group.by = "SCT_snn_res.0.2")#RNA_
DimPlot(seu_obj, reduction = "umap", group.by = "sheet_patient_id")
DimPlot(seu_obj, reduction = "umap", group.by = "RNA_snn_res.0.2", split.by = "sheet_patient_id")
#seu_obj$RNA_snn_res.0.2 = factor(seu_obj$RNA_snn_res.0.2, levels = c("0", "1", "2", "3", "4", "5", "6"))
png("umap_BASAL_0.3_v1.png", width=4, height=3, units = "in", res = 300)
DimPlot(seu_obj, reduction = "umap", group.by = "RNA_snn_res.0.3", cols = colors_dark)
dev.off()
png("umap_BASAL_0.3_splitptv1.png", width=12, height=12, units = "in", res = 300)
DimPlot(seu_obj, reduction = "umap", group.by = "RNA_snn_res.0.3", split.by = "sheet_patient_id", pt.size = 1,cols = colors_dark, ncol = 4)
dev.off()


Idents(seu_obj) <- "RNA_snn_res.0.2"

mks = FindAllMarkers(seu_obj,assay = "RNA")
mks %>%  
  group_by(cluster) %>%
  slice_max(avg_logFC, n = 10) -> top10;write.csv(mks, "mks_0.2.csv")
mks = read.csv("mks_0.3.csv")
png("heatmap_FULL_all_labelled_BASAL_0.2_v1.png", width=12, height=9, units = "in", res = 300)
DoHeatmap(seu_obj, features = as.character(top10$gene), group.colors = colors_dark, assay = "RNA") + NoLegend()
dev.off()
tavg = AverageExpression(seu_obj, features = unique(c("EPCAM", "ACTA2","MYLK","TAGLN", "KRT5", as.character(top10$gene))), return.seurat = T)
DoHeatmap(tavg, features = unique(c("EPCAM", "ACTA2","MYLK","TAGLN", "KRT5", as.character(top10$gene))), assay = "RNA", label = TRUE,
          slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = colors_dark, raster = F)  +
  scale_fill_gradient2(low = "#29A7B8", mid = "white", high = "#E85041", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=9)) -> tcell_heatmap
png("heatmap_all_labelled_BASAL_0.2_v1.png", width=7, height=12, units = "in", res = 300)
tcell_heatmap
dev.off()

# **** \\\ check my clusters ------
tk_obj = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumhr_v2/clus_out.rds")
kevin_obj = load("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumhr_v10/hbca_lumhr.RObj")
kevin_obj = hbca.lumhr

tk_obj_df = data.frame(cellname = tk_obj$cellname, tk_states = tk_obj$cellstate_tk2)
clus_out_meta = seu_obj@meta.data
clus_out_meta = left_join(clus_out_meta, tk_obj_df)

seu_obj@meta.data = clus_out_meta
seu_obj$RNA_snn_res.0.2 = as.character(seu_obj$RNA_snn_res.0.2)
seu_obj$RNA_snn_res.0.2[(is.na(seu_obj$RNA_snn_res.0.2))] = "unkn"
rownames(seu_obj@meta.data) = seu_obj$cellname

kevin_obj_df = data.frame(cellname = kevin_obj$cellname, kevin_states = kevin_obj$lumhr_states)
clus_out_meta = seu_obj@meta.data
clus_out_meta = left_join(clus_out_meta, kevin_obj_df)

seu_obj@meta.data = clus_out_meta
seu_obj$kevin_states = as.character(seu_obj$kevin_states)
seu_obj$kevin_states[(is.na(seu_obj$kevin_states))] = "unkn"
rownames(seu_obj@meta.data) = seu_obj$cellname


p00 = DimPlot(tk_obj, reduction  = "umap", group.by = "cellstate_tk2", ncol = 1, label = T, cols = colors_dark, pt.size = 0.0001)
p0 = DimPlot(seu_obj, reduction  = "umap", group.by = "kevin_states", ncol = 1, label = T, cols = colors_dark, pt.size = 0.0001)
p1 = DimPlot(seu_obj, reduction  = "umap", group.by = "tk_states", ncol = 1, label = T, cols = colors_dark, pt.size = 0.0001)
p2 = DimPlot(seu_obj, reduction  = "umap", group.by = "RNA_snn_res.0.2", ncol = 1, label = T, cols = colors_dark, pt.size = 0.0001)
p3 = DimPlot(seu_obj, reduction  = "umap", group.by = "tk_states", split.by = "final_sample_id", ncol = 6, label = F, cols = colors_dark, pt.size = 1)
png("umap_tk_v2.png", width=20, height=4, units = "in", res = 300) #
cowplot::plot_grid(p00, p0, p1, p2, ncol = 4)
dev.off()

png("umap_tk_split_by_pt_v2.png", width=11, height=11, units = "in", res = 300) #
p3
dev.off()

png("umap_v2.png", width=4, height=4, units = "in", res = 300) #
p2
dev.off()

png("umap_pt_v2.png", width=6, height=4, units = "in", res = 300) #
DimPlot(seu_obj, reduction  = "umap", group.by = "sheet_patient_id", ncol = 1, label = T, cols = colors_dark, pt.size = 0.0001)
dev.off()

# *** recluster with scaling umi count effecrs ---------
seu_object = read_rds(glue("{save_path}/epi_basal_sub_3hr_b2_reduction_v1.rds"))
DefaultAssay(seu_object) <- "RNA"
seu_obj = NormalizeData(seu_object) %>% FindVariableFeatures(nfeatures = 5000) 
all.genes <- rownames(seu_obj)
seu_obj <- ScaleData(seu_obj, features = all.genes, vars.to.regress = "nCount_RNA")
seu_obj <- RunPCA(seu_obj, features = VariableFeatures(object = seu_obj))
seu_obj <- FindNeighbors(seu_obj, dims = 1:35)
seu_obj <- FindClusters(seu_obj, resolution = c(0.5, 0.1, 0.05, 0.075, 0.2, 0.15, 0.3))
seu_obj <- RunUMAP(seu_obj, dims = 1:35)
DimPlot(seu_obj, reduction = "umap", group.by = "RNA_snn_res.0.2")
DimPlot(seu_obj, reduction = "umap", group.by = "sheet_patient_id")
DimPlot(seu_obj, reduction = "umap", group.by = "RNA_snn_res.0.2", split.by = "sheet_patient_id")
seu_obj$RNA_snn_res.0.2 = factor(seu_obj$RNA_snn_res.0.2, levels = c("0", "1", "2", "3", "4", "5", "6"))
Idents(seu_obj) <- "RNA_snn_res.0.2"

mks = FindAllMarkers(seu_obj,assay = "RNA")
mks %>%  
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_logFC) -> top10
DoHeatmap(seu_obj, features = top10$gene, group.colors = colors_dark) + NoLegend()
tavg = AverageExpression(seu_obj, features = unique(c("EPCAM", "ESR1","PGR","AR", "TOP2A", top10$gene)), return.seurat = T)
DoHeatmap(tavg, features = unique(c("EPCAM", "ESR1","PGR","AR", "TOP2A",top10$gene)), assay = "RNA", label = TRUE,
          slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = colors_dark, raster = F)  +
  scale_fill_gradient2(low = "#29A7B8", mid = "white", high = "#E85041", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5)) -> tcell_heatmap
png("heatmap_all_labelled_lumHR_0.2_v2.png", width=5, height=16, units = "in", res = 300)
tcell_heatmap
dev.off()

# **** \\\ check my clusters ------
tk_obj = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumhr_v2/clus_out.rds")
kevin_obj = load("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumhr_v10/hbca_lumhr.RObj")
kevin_obj = hbca.lumhr

tk_obj_df = data.frame(cellname = tk_obj$cellname, tk_states = tk_obj$cellstate_tk2)
clus_out_meta = seu_obj@meta.data
clus_out_meta = left_join(clus_out_meta, tk_obj_df)

seu_obj@meta.data = clus_out_meta
seu_obj$RNA_snn_res.0.2 = as.character(seu_obj$RNA_snn_res.0.2)
seu_obj$RNA_snn_res.0.2[(is.na(seu_obj$RNA_snn_res.0.2))] = "unkn"
rownames(seu_obj@meta.data) = seu_obj$cellname

kevin_obj_df = data.frame(cellname = kevin_obj$cellname, kevin_states = kevin_obj$lumhr_states)
clus_out_meta = seu_obj@meta.data
clus_out_meta = left_join(clus_out_meta, kevin_obj_df)

seu_obj@meta.data = clus_out_meta
seu_obj$kevin_states = as.character(seu_obj$kevin_states)
seu_obj$kevin_states[(is.na(seu_obj$kevin_states))] = "unkn"
rownames(seu_obj@meta.data) = seu_obj$cellname


p00 = DimPlot(tk_obj, reduction  = "umap", group.by = "cellstate_tk2", ncol = 1, label = T, cols = colors_dark, pt.size = 0.0001)
p0 = DimPlot(seu_obj, reduction  = "umap", group.by = "kevin_states", ncol = 1, label = T, cols = colors_dark, pt.size = 0.0001)
p1 = DimPlot(seu_obj, reduction  = "umap", group.by = "tk_states", ncol = 1, label = T, cols = colors_dark, pt.size = 0.0001)
p2 = DimPlot(seu_obj, reduction  = "umap", group.by = "RNA_snn_res.0.2", ncol = 1, label = T, cols = colors_dark, pt.size = 0.0001)
p3 = DimPlot(seu_obj, reduction  = "umap", group.by = "tk_states", split.by = "final_sample_id", ncol = 6, label = F, cols = colors_dark, pt.size = 1)
png("umap_tk_v3.png", width=20, height=4, units = "in", res = 300) #
cowplot::plot_grid(p00, p0, p1, p2, ncol = 4)
dev.off()

png("umap_tk_split_by_pt_v3.png", width=11, height=11, units = "in", res = 300) #
p3
dev.off()

stress_sig = list(c("CXCL8","CXCL1","HMOX1","SAT1","SOD2","CXCL3","AKR1C1","GAPDH","TXN","CXCL2","EDNRB","MEDAG","CD59",    
                    "EIF5A","CLEC2B","TXNRD1","TPI1","SERPINE1", "SLC43A3","PRDX1","PGAM1","WTAP","ADRM1","EIF3I","PRELID1","RAD23A",  
                    "TUBB","PSMB2","ATP1A1","MPZL1","NEAT1","PIM3","AKR1C2","HSPA5"))

# *** recluster with scaling patient effecrs ---------
seu_object = read_rds(glue("{save_path}/epi_basal_sub_3hr_b2_reduction_v1.rds"))
DefaultAssay(seu_object) <- "RNA"
seu_obj = NormalizeData(seu_object) %>% FindVariableFeatures(nfeatures = 5000) 
seu_obj@assays$RNA@var.features = seu_obj@assays$RNA@var.features[!seu_obj@assays$RNA@var.features %in% stress_sig[[1]]]
all.genes <- rownames(seu_obj)
seu_obj <- ScaleData(seu_obj, features = all.genes, vars.to.regress = "sheet_patient_id")
seu_obj <- RunPCA(seu_obj, features = VariableFeatures(object = seu_obj))
seu_obj <- FindNeighbors(seu_obj, dims = 1:18)
seu_obj <- FindClusters(seu_obj, resolution = c(0.1, 0.05, 0.075, 0.2))
seu_obj <- RunUMAP(seu_obj, dims = 1:18)
seu_obj = read_rds("seu_obj_stressscaled.rds")
DimPlot(seu_obj, reduction = "umap", group.by = "RNA_snn_res.0.2")
DimPlot(seu_obj, reduction = "umap", group.by = "sheet_patient_id")
DimPlot(seu_obj, reduction = "umap", group.by = "RNA_snn_res.0.2", split.by = "sheet_patient_id")

png("umap_tk_bypt_0.3.png", width=4, height=3, units = "in", res = 300) #
DimPlot(seu_obj, reduction = "umap", group.by = "RNA_snn_res.0.3")
dev.off()
png("umap_tk_bypt_0.2pt.png", width=4, height=3, units = "in", res = 300) #
DimPlot(seu_obj, reduction = "umap", group.by = "sheet_patient_id")
dev.off()

Idents(seu_obj) <- "RNA_snn_res.0.2"

mks = FindAllMarkers(seu_obj,assay = "RNA");write.csv(mks, "mks_stress_0.2.csv")
mks %>%  
  group_by(cluster) %>%
  slice_max(avg_logFC, n = 20) -> top10
DoHeatmap(seu_obj, features = top10$gene, group.colors = colors_dark) + NoLegend()
tavg = AverageExpression(seu_obj, features = unique(c("EPCAM", "ACTA2","TOP2A","MYLK", "SLPI", top10$gene)), return.seurat = T)
DoHeatmap(tavg, features = unique(c("EPCAM", "ACTA2","TOP2A","MYLK", "SLPI",top10$gene)), assay = "RNA", label = TRUE,
          slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = colors_dark, raster = F)  +
  scale_fill_gradient2(low = "#29A7B8", mid = "white", high = "#E85041", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5)) -> tcell_heatmap
png("heatmap_all_labelled_BASAL_0.2_v3.png", width=5, height=16, units = "in", res = 300)
tcell_heatmap
dev.off()

png("featureplot_BASAL_0.2_v1.png", width=16, height=16, units = "in", res = 300)
FeaturePlot(seu_obj, c("rna_EPCAM", "rna_KRT5", "rna_KRT14", "rna_ACTA2", "rna_MYLK", "rna_TAGLN", "rna_TOP2A", "rna_SLPI"))
FeaturePlot(seu_obj, c("TNFRSF12A", "TPM4", "HES1", "EGR1", "AREG", "TAGLN", "ACTA2", "OSMR", "BCL3"), order = T)
dev.off()

# *** recluster with SCT scaling patient effecrs ---------
seu_object = read_rds(glue("{save_path}/epi_basal_sub_3hr_b2_reduction_v1.rds"))
DefaultAssay(seu_object) <- "RNA"
seu_obj = SCTransform(seu_object, vars.to.regress = "sheet_patient_id")
seu_obj@assays$SCT@var.features = seu_obj@assays$SCT@var.features[!seu_obj@assays$SCT@var.features %in% stress_sig[[1]]]
all.genes <- rownames(seu_obj)
#seu_obj <- ScaleData(seu_obj, features = all.genes, vars.to.regress = "sheet_patient_id")
seu_obj <- RunPCA(seu_obj, features = VariableFeatures(object = seu_obj))
seu_obj <- FindNeighbors(seu_obj, dims = 1:18)
seu_obj <- FindClusters(seu_obj, resolution = c(0.1, 0.05, 0.075, 0.2))
seu_obj <- RunUMAP(seu_obj, dims = 1:18)
write_rds(seu_obj, "seu_obj_SCT_stressscaled.rds")
seu_obj = read_rds("seu_obj_SCT_stressscaled.rds")
DimPlot(seu_obj, reduction = "umap", group.by = "SCT_snn_res.0.2")
DimPlot(seu_obj, reduction = "umap", group.by = "sheet_patient_id")
DimPlot(seu_obj, reduction = "umap", group.by = "RNA_snn_res.0.2", split.by = "sheet_patient_id")

png("umap_tk_bypt_0.3.png", width=4, height=3, units = "in", res = 300) #
DimPlot(seu_obj, reduction = "umap", group.by = "RNA_snn_res.0.3")
dev.off()
png("umap_tk_bypt_0.2pt.png", width=4, height=3, units = "in", res = 300) #
DimPlot(seu_obj, reduction = "umap", group.by = "sheet_patient_id")
dev.off()

Idents(seu_obj) <- "SCT_snn_res.0.2"

mks = FindAllMarkers(seu_obj,assay = "RNA");write.csv(mks, "mks_stress_SCT_0.2.csv")
mks %>%  
  group_by(cluster) %>%
  slice_max(avg_logFC, n = 12) -> top10
DoHeatmap(seu_obj, features = top10$gene, group.colors = colors_dark) + NoLegend()
tavg = AverageExpression(seu_obj, features = unique(c("EPCAM", "ACTA2","TOP2A","MYLK", "SLPI", top10$gene)), return.seurat = T)
DoHeatmap(tavg, features = unique(c("EPCAM", "ACTA2","TOP2A","MYLK", "SLPI",top10$gene)), assay = "RNA", label = TRUE,
          slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = colors_dark, raster = F)  +
  scale_fill_gradient2(low = "#29A7B8", mid = "white", high = "#E85041", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5)) -> tcell_heatmap
png("heatmap_all_labelled_BASAL_0.2_v3.png", width=5, height=16, units = "in", res = 300)
tcell_heatmap
dev.off()

png("featureplot_BASAL_0.2__SCT_v1.png", width=16, height=16, units = "in", res = 300)
FeaturePlot(seu_obj, c("TNFRSF12A", "TPM4", "HES1", "EGR1", "AREG", "TAGLN", "ACTA2", "OSMR", "BCL3"), order = T)
dev.off()


# TEST NMF ----------
# remotes::install_github("jbergenstrahle/STUtility")
# library(NMF)
mat_nmf <- as.matrix(seu_obj@assays$SCT@scale.data)
sc_nonneg <- function(x){x-min(x)}
data = apply(mat_nmf, 1, sc_nonneg) %>% t() # need transform since apply transformed the matrix
data = data[which(rowSums(data)!=0),]
data = data[,which(colSums(data)!=0)]
data_nmf <- nmf(data, rank = 2:4, nrun=10,  seed=21,.opt='vP30')
data_nmf = read_rds("data_nmf.rds")

# score
pdf("metricsnmf.pdf")
plot(data_nmf)
dev.off()

pdf("consensusnmf.pdf")
consensusmap(data_nmf, labCol=NA, labRow=NA)

rank2 = data_nmf$fit$'2'
dim(rank2@consensus)  
sl_rank2 = silhouette(rank2)
df_sl_rank2 = sl_rank2[, 1:3] %>% data.frame() %>% as_tibble()
rownames(df_sl_rank2) = 

# *** \\\ cluster single patient --------------
seu_object = read_rds(glue("{save_path}/epi_basal_sub_3hr_b2_reduction_v1.rds"))

Idents(seu_object) <- "sheet_patient_id"
pt54 = subset(seu_object, idents = "pt54")
DefaultAssay(pt54) <- "RNA"
pt54 = NormalizeData(pt54) %>% FindVariableFeatures(nfeatures = 5000) 
all.genes <- rownames(pt54)
pt54 <- ScaleData(pt54, features = all.genes)
pt54 <- RunPCA(pt54, features = VariableFeatures(object = pt54))
pt54 <- FindNeighbors(pt54, dims = 1:12)
pt54 <- FindClusters(pt54, resolution = c(0.5, 0.1, 0.05, 0.075, 0.2, 0.15, 0.3))
pt54 <- RunUMAP(pt54, dims = 1:12)
DimPlot(pt54, reduction = "umap", sgroup.by = "RNA_snn_res.0.2")
Idents(pt54) <- "RNA_snn_res.0.2"
mks = FindAllMarkers(pt54,assay = "RNA");write.csv(mks, "mks_pt54_0.2.csv")
mks %>%  
  group_by(cluster) %>%
  slice_max(avg_logFC, n = 12) -> top10
DoHeatmap(pt54, features = top10$gene, group.colors = colors_dark) + NoLegend()

Idents(seu_object) <- "sheet_patient_id"
pt13 = subset(seu_object, idents = "pt13")
DefaultAssay(pt13) <- "RNA"
pt13 = NormalizeData(pt13) %>% FindVariableFeatures(nfeatures = 5000) 
all.genes <- rownames(pt13)
pt13 <- ScaleData(pt13, features = all.genes)
pt13 <- RunPCA(pt13, features = VariableFeatures(object = pt13))
pt13 <- FindNeighbors(pt13, dims = 1:12)
pt13 <- FindClusters(pt13, resolution = c(0.5, 0.1, 0.05, 0.075, 0.2, 0.15, 0.3))
pt13 <- RunUMAP(pt13, dims = 1:12)
DimPlot(pt13, reduction = "umap", group.by = "RNA_snn_res.0.2")
Idents(pt13) <- "RNA_snn_res.0.2"
mks = FindAllMarkers(pt13,assay = "RNA");write.csv(mks, "mks_pt13_0.2.csv")
mks %>%  
  group_by(cluster) %>%
  slice_max(avg_logFC, n = 12) -> top10
DoHeatmap(pt13, features = top10$gene, group.colors = colors_dark) + NoLegend()

Idents(seu_object) <- "sheet_patient_id"
pt18 = subset(seu_object, idents = "pt18")
DefaultAssay(pt18) <- "RNA"
pt18 = NormalizeData(pt18) %>% FindVariableFeatures(nfeatures = 5000) 
all.genes <- rownames(pt18)
pt18 <- ScaleData(pt18, features = all.genes)
pt18 <- RunPCA(pt18, features = VariableFeatures(object = pt18))
pt18 <- FindNeighbors(pt18, dims = 1:12)
pt18 <- FindClusters(pt18, resolution = c(0.5, 0.1, 0.05, 0.075, 0.2, 0.15, 0.3))
pt18 <- RunUMAP(pt18, dims = 1:12)
DimPlot(pt18, reduction = "umap", group.by = "RNA_snn_res.0.2")
Idents(pt18) <- "RNA_snn_res.0.2"
mks = FindAllMarkers(pt18,assay = "RNA");write.csv(mks, "mks_pt18_0.2.csv")
mks %>%  
  group_by(cluster) %>%
  slice_max(avg_logFC, n = 12) -> top10
DoHeatmap(pt18, features = top10$gene, group.colors = colors_dark) + NoLegend()

Idents(seu_object) <- "sheet_patient_id"
pt57 = subset(seu_object, idents = "pt57")
DefaultAssay(pt57) <- "RNA"
pt57 = NormalizeData(pt57) %>% FindVariableFeatures(nfeatures = 5000) 
all.genes <- rownames(pt57)
pt57 <- ScaleData(pt57, features = all.genes)
pt57 <- RunPCA(pt57, features = VariableFeatures(object = pt57))
pt57 <- FindNeighbors(pt57, dims = 1:12)
pt57 <- FindClusters(pt57, resolution = c(0.5, 0.1, 0.05, 0.075, 0.2, 0.15, 0.3))
pt57 <- RunUMAP(pt57, dims = 1:12)
DimPlot(pt57, reduction = "umap", group.by = "RNA_snn_res.0.2")
Idents(pt57) <- "RNA_snn_res.0.2"
mks = FindAllMarkers(pt57,assay = "RNA");write.csv(mks, "mks_pt57_0.2.csv")
mks %>%  
  group_by(cluster) %>%
  slice_max(avg_logFC, n = 12) -> top10
DoHeatmap(pt57, features = top10$gene) + NoLegend()

stress_sig = c("JUNB", "JUN", "FOS","FOSB", "DNAJB1", "HSPA1A", "HSPA1B", "CXCL8","CXCL1","HMOX1","SAT1","SOD2","CXCL3","AKR1C1","GAPDH","TXN","CXCL2","EDNRB","MEDAG","CD59",    
                    "EIF5A","CLEC2B","TXNRD1","TPI1","SERPINE1", "SLC43A3","PRDX1","PGAM1","WTAP","ADRM1","EIF3I","PRELID1","RAD23A",  
                    "TUBB","PSMB2","ATP1A1","MPZL1","NEAT1","PIM3","AKR1C2","HSPA5")
Idents(seu_object) <- "sheet_patient_id"
pt57 = subset(seu_object, idents = "pt57")
DefaultAssay(pt57) <- "RNA"
pt57 = NormalizeData(pt57) %>% FindVariableFeatures(nfeatures = 5000)
pt57@assays$RNA@var.features = pt57@assays$RNA@var.features[!pt57@assays$RNA@var.features %in% stress_sig]
all.genes <- rownames(pt57)
pt57 <- ScaleData(pt57)
pt57 <- RunPCA(pt57, features = VariableFeatures(object = pt57))
pt57 <- FindNeighbors(pt57, dims = 1:12)
pt57 <- FindClusters(pt57, resolution = c(0.5, 0.1, 0.05, 0.075, 0.2, 0.15, 0.3))
pt57 <- RunUMAP(pt57, dims = 1:12)
DimPlot(pt57, reduction = "umap", group.by = "RNA_snn_res.0.2")
Idents(pt57) <- "RNA_snn_res.0.2"
mks = FindAllMarkers(pt57,assay = "RNA");write.csv(mks, "mks_pt57_0.2.csv")
mks %>%  
  group_by(cluster) %>%
  slice_max(avg_logFC, n = 12) -> top10
DoHeatmap(pt57, features = top10$gene) + NoLegend()

Idents(tum) <- "sheet_patient_id"
pt56 = subset(tum, idents = "pt56")
DefaultAssay(pt56) <- "RNA"
pt56 = NormalizeData(pt56) %>% FindVariableFeatures(nfeatures = 5000) 
all.genes <- rownames(pt56)
pt56 <- ScaleData(pt56, features = all.genes)
pt56 <- RunPCA(pt56, features = VariableFeatures(object = pt56))
pt56 <- FindNeighbors(pt56, dims = 1:12)
pt56 <- FindClusters(pt56, resolution = c(0.5, 0.1, 0.05, 0.075, 0.2, 0.15, 0.3))
pt56 <- RunUMAP(pt56, dims = 1:12)
DimPlot(pt56, reduction = "umap", group.by = "RNA_snn_res.0.2")
Idents(pt56) <- "RNA_snn_res.0.2"
mks = FindAllMarkers(pt56,assay = "RNA");write.csv(mks, "mks_pt56_0.2.csv")
mks %>%  
  group_by(cluster) %>%
  slice_max(avg_logFC, n = 12) -> top10
DoHeatmap(pt56, features = top10$gene) + NoLegend()

stress_sig = c("JUNB", "JUN", "FOS","FOSB", "DNAJB1", "HSPA1A", "HSPA1B", "CXCL8","CXCL1","HMOX1","SAT1","SOD2","CXCL3","AKR1C1","GAPDH","TXN","CXCL2","EDNRB","MEDAG","CD59",    
               "EIF5A","CLEC2B","TXNRD1","TPI1","SERPINE1", "SLC43A3","PRDX1","PGAM1","WTAP","ADRM1","EIF3I","PRELID1","RAD23A",  
               "TUBB","PSMB2","ATP1A1","MPZL1","NEAT1","PIM3","AKR1C2","HSPA5")
Idents(tum) <- "sheet_patient_id"
pt56 = subset(tum, idents = "pt56")
DefaultAssay(pt56) <- "RNA"
pt56 = NormalizeData(pt56) %>% FindVariableFeatures(nfeatures = 5000)
pt56@assays$RNA@var.features = pt56@assays$RNA@var.features[!pt56@assays$RNA@var.features %in% stress_sig]
all.genes <- rownames(pt56)
pt56 <- ScaleData(pt56)
pt56 <- RunPCA(pt56, features = VariableFeatures(object = pt56))
pt56 <- FindNeighbors(pt56, dims = 1:12)
pt56 <- FindClusters(pt56, resolution = c(0.5, 0.1, 0.05, 0.075, 0.2, 0.15, 0.3))
pt56 <- RunUMAP(pt56, dims = 1:12)
DimPlot(pt56, reduction = "umap", group.by = "RNA_snn_res.0.2")
Idents(pt56) <- "RNA_snn_res.0.2"
mks = FindAllMarkers(pt56,assay = "RNA");write.csv(mks, "mks_pt56_0.2.csv")
mks %>%  
  group_by(cluster) %>%
  slice_max(avg_logFC, n = 12) -> top10
DoHeatmap(pt56, features = top10$gene) + NoLegend()


# **** \\\ check my clusters ------
tk_obj = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumhr_v2/clus_out.rds")
kevin_obj = load("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumhr_v10/hbca_lumhr.RObj")
kevin_obj = hbca.lumhr

tk_obj_df = data.frame(cellname = tk_obj$cellname, tk_states = tk_obj$cellstate_tk2)
clus_out_meta = seu_obj@meta.data
clus_out_meta = left_join(clus_out_meta, tk_obj_df)

seu_obj@meta.data = clus_out_meta
seu_obj$RNA_snn_res.0.2 = as.character(seu_obj$RNA_snn_res.0.2)
seu_obj$RNA_snn_res.0.2[(is.na(seu_obj$RNA_snn_res.0.2))] = "unkn"
rownames(seu_obj@meta.data) = seu_obj$cellname

kevin_obj_df = data.frame(cellname = kevin_obj$cellname, kevin_states = kevin_obj$lumhr_states)
clus_out_meta = seu_obj@meta.data
clus_out_meta = left_join(clus_out_meta, kevin_obj_df)

seu_obj@meta.data = clus_out_meta
seu_obj$kevin_states = as.character(seu_obj$kevin_states)
seu_obj$kevin_states[(is.na(seu_obj$kevin_states))] = "unkn"
rownames(seu_obj@meta.data) = seu_obj$cellname


p00 = DimPlot(tk_obj, reduction  = "umap", group.by = "cellstate_tk2", ncol = 1, label = T, cols = colors_dark, pt.size = 0.0001)
p0 = DimPlot(seu_obj, reduction  = "umap", group.by = "kevin_states", ncol = 1, label = T, cols = colors_dark, pt.size = 0.0001)
p1 = DimPlot(seu_obj, reduction  = "umap", group.by = "tk_states", ncol = 1, label = T, cols = colors_dark, pt.size = 0.0001)
p2 = DimPlot(seu_obj, reduction  = "umap", group.by = "RNA_snn_res.0.2", ncol = 1, label = T, cols = colors_dark, pt.size = 0.0001)
p3 = DimPlot(seu_obj, reduction  = "umap", group.by = "tk_states", split.by = "final_sample_id", ncol = 6, label = F, cols = colors_dark, pt.size = 1)
png("umap_tk_v4.png", width=20, height=4, units = "in", res = 300) #
cowplot::plot_grid(p00, p0, p1, p2, ncol = 4)
dev.off()

png("umap_tk_split_by_pt_v4.png", width=11, height=11, units = "in", res = 300) #
p3
dev.off()

png("umap_tk_split_by_pt_v4_check.png", width=11, height=11, units = "in", res = 300) #
DimPlot(seu_obj, reduction  = "umap", group.by = "RNA_snn_res.0.2", split.by = "final_sample_id", ncol = 6, label = F, cols = colors_dark, pt.size = 1)
dev.off()

VlnPlot(seu_obj, "nCount_RNA", group.by = "RNA_snn_res.0.2")

# END HERE -----------

# method 9 ----------
# Log T integration with no scaling
library(future)
plan("multiprocess", workers = 85)
#plan("sequential")

# ** Clus function --------
# run clustering function
#seu_object = read_rds(glue("{save_path}/t_subset_v1.rds"))
seu_object = read_rds(glue("{save_path}/epi_lumhr_sub_3hr_b2_reduction_v1.rds"))
DefaultAssay(seu_object) <- "RNA"
#seu_object$new_pt = 
dims = 35; k.param = 20; cols = colors_dark; var_split = "sheet_patient_id"; var_scale = "nCount_RNA"
fts1 = c("nFeature_RNA","nCount_RNA", "percent.mito","percent.rb", "percent.stress")
fts2 = c("rna_EPCAM","rna_PTPRC","rna_PECAM1","rna_DCN","rna_KRT5",  "rna_CD86", "rna_COL1A1", "rna_PIP")
fts3 = c("rna_S100A8", "rna_SERPINB4", "rna_FBLN5", "rna_CCL2", "rna_LTF","rna_SAA1", "rna_TNFAIP2", "rna_KYNU", "rna_NDRG1")
fts4 = c("rna_KRT15", "rna_KRT14", "rna_SCGB1D2", "rna_PTN", "rna_CITED4","rna_MYC", "rna_LEF3", "rna_CALML5", "rna_RGS2",  "rna_PCNA", "rna_MKI67", "rna_LTF")


ref_set = NULL

clus_out = func_cluster_vst(tum_int = seu_object, var_scale = var_scale, var_split = var_split, dims = dims, k.param = k.param, 
                            markers1 = fts1, markers2 = fts2, markers3 = fts3, markers4 = fts4, save_path = save_path,
                            cols = colors_dark, ref_set = ref_set)

char_ident = "integrated_snn_res.0.3"

p1 = DimPlot(clus_out, split.by = "sheet_patient_id", group.by = "cellstate_tk1", cols = colors_darj, label = F, shuffle = T, pt.size = 0.6, ncol = 10)
png("vln_allpts_v1.png", width=24, height=20, units = "in", res = 300)
p1
dev.off()

p1 = DimPlot(clus_out, split.by = "sheet_patient_id", group.by = "integrated_snn_res.0.2", 
             cols = colors_dark, label = F, shuffle = F, pt.size = 0.6, ncol = 10)
png("dim_allpts_v1.png", width=24, height=20, units = "in", res = 300)
p1
dev.off()

p2 = DimPlot(clus_out, group.by = "integrated_snn_res.0.4", cols = colors_dark, label = T)
png("umap0.2_v1.png", width=4, height=3, units = "in", res = 300)
p2
dev.off()
FeaturePlot(clus_out, features = c("rna_SCGB"))
VlnPlot(clus_out, group.by = "integrated_snn_res.0.4", features = c("nCount_RNA", "percent.mito"), pt.size = 0)
VlnPlot(clus_out, group.by = "integrated_snn_res.0.4", features = c("EPCAM", "SLPI", "LTF", "KRT5", 
                                                                    "KRT6B", "KRT14", "KRT15"), pt.size = 0, cols = colors_dark, assay = "RNA")


# ** Markers function (9)--------
dims = 35; k.param = 20;
char_ident = "integrated_snn_res.0.2"  #integrated_snn_res.0.1
ht_fts1 = c("EPCAM", "KRT17", "AR")
#clus_out = read_rds(file = glue("{save_path}/dim_35k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")) 
mark_out = func_markers_vst(clus_out = clus_out, char_ident = char_ident, fts1 = ht_fts1, dims = dims, k.param = k.param, cols = colors_dark, 
                            save_path = save_path)

clus_out = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_v2/dim_35k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")
markers = read.csv("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_v2/dims_35k_20_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.2_vst_method_9.csv")
# clus_out = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v1/dim_30k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")
# markers = read.csv("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v1/dims_30k_20_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.4_vst_method_9.csv")
tum = clus_out

markers_v2 = markers %>% dplyr::filter(gene %in% markers$gene)

top_markers = markers_v2 %>%
  dplyr::group_by(cluster) %>%
  top_n(12, avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()

#Idents(tum.integrated.ns)
#tum.integrated.ns$clusters1 <- factor(tum.integrated.ns$clusters1, levels = levels(Idents(tum.integrated.ns)))
message("Making Heatmaps")
DefaultAssay(tum) <- "integrated"
Idents(tum) <- "integrated_snn_res.0.2"
p1 = DoHeatmap(tum, features = c("EPCAM", top_markers), size = 2,raster = FALSE,
               group.by = char_ident, slot = "scale.data", disp.min = -1.5, disp.max = 1.5, group.colors = colors_dark) + 
  scale_color_manual(values = colors_dark) +
  scale_fill_gradient2(low="steelblue1", high = "tomato3", mid = "white", midpoint = 0) + 
  theme(axis.text.y.left = element_text(size=5))

png("heatmap_all_labelled_res2_cells_v1.png", width=12, height=8, units = "in", res = 200)
p1
dev.off()

tavg = AverageExpression(tum, features = unique(c("EPCAM", "KIT", "TOP2A", top_markers)), return.seurat = T)
DoHeatmap(tavg, features = unique(c("EPCAM", "KIT", "TOP2A", top_markers)), assay = "integrated", label = TRUE,
          slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = colors_dark, raster = F)  +
  scale_fill_gradient2(low = "#29A7B8", mid = "white", high = "#E85041", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5)) -> tcell_heatmap
png("heatmap_all_labelled_lumsec_0.2_v1.png", width=12, height=12, units = "in", res = 300)
tcell_heatmap
dev.off()


# add module score ------
basal_sig = list(c("KRT14","KRT17","KRT5","ACTG2","TUBB2B","COL17A1","LAMA3","SERPINB5"))
lumhr_sig = list(c("AREG","ANKRD30A","AGR2","AGR3","TMC5","DNAJC12","RASEF","MLPH","CA12","C3orf52","ACADSB","MAPK8"))
lumsec_sig = list(c("SLPI","LTF","KRT15","MMP7","CCL28","PIGR","RARRES1","ALDH1A3","LCN2","SNORC","S100A1","GABRP"))
clus_out = AddModuleScore(object = clus_out, features = basal_sig, assay = "RNA", name = "basal_sig")
clus_out = AddModuleScore(object = clus_out, features = lumhr_sig, assay = "RNA", name = "lumhr_sig")
clus_out = AddModuleScore(object = clus_out, features = lumsec_sig, assay = "RNA", name = "lumsec_sig")
FeaturePlot(clus_out, features = c("basal_sig1", "lumhr_sig1", "lumsec_sig1"), cols = c("grey89", "maroon"))
VlnPlot(clus_out, features = c("basal_sig1", "lumhr_sig1", "lumsec_sig1"), group.by = "integrated_snn_res.0.2", pt.size = 0)

VlnPlot(clus_out, features = c("KRT14","KRT17","KRT5","ACTG2","TUBB2B","COL17A1","LAMA3","SERPINB5"), group.by = "integrated_snn_res.0.2", pt.size = 0)
VlnPlot(clus_out, features = c("AREG","ANKRD30A","AGR2","AGR3","TMC5","DNAJC12","RASEF","MLPH","CA12","C3orf52","ACADSB","MAPK8"), group.by = "integrated_snn_res.0.2", pt.size = 0)
VlnPlot(clus_out, features = c("SLPI","LTF","KRT15","MMP7","CCL28","PIGR","RARRES1","ALDH1A3","LCN2","SNORC","S100A1","GABRP"), group.by = "integrated_snn_res.0.2", pt.size = 0, assay = "RNA")

DimPlot(clus_out, group.by = "integrated_snn_res.0.2", label = T, cols = colors_darj)
DimPlot(clus_out, group.by = "seclum_states", label = T,cols = colors_darj)
FeaturePlot(clus_out, features = c("rna_PCNA", "rna_KRT14", "rna_LALBA", "rna_CXCL8"))


# assign clusters ------
clus_out = tum
clus_out$cellstate_tk1 = "NULL"
clus_out$cellstate_tk1[clus_out$integrated_snn_res.0.2 %in% c("0")] = "ls_KRT23"
clus_out$cellstate_tk1[clus_out$integrated_snn_res.0.2 %in% c("1")] = "ls_SPINK5"
clus_out$cellstate_tk1[clus_out$integrated_snn_res.0.2 %in% c("2")] = "ls_HMOX1"
clus_out$cellstate_tk1[clus_out$integrated_snn_res.0.2 %in% c("3")] = "ls_MT1X"
clus_out$cellstate_tk1[clus_out$integrated_snn_res.0.2 %in% c("4")] = "ls_KIT"
clus_out$cellstate_tk1[clus_out$integrated_snn_res.0.2 %in% c("5")] = "ls_lactation"
clus_out$cellstate_tk1[clus_out$integrated_snn_res.0.2 %in% c("6")] = "ls_cd74"
clus_out$cellstate_tk1[clus_out$integrated_snn_res.0.2 %in% c("7")] = "ls_prol"

lvl = names(table(clus_out$cellstate_tk1))
clus_out$cellstate_tk1 = factor(clus_out$cellstate_tk1, levels = lvl)

clus_out$cellstate_tk2 = "NULL"
clus_out$cellstate_tk2[clus_out$integrated_snn_res.0.2 %in% c("0")] = "ls_KRT23"
clus_out$cellstate_tk2[clus_out$integrated_snn_res.0.2 %in% c("1")] = "ls_SPINK5"
clus_out$cellstate_tk2[clus_out$integrated_snn_res.0.2 %in% c("2")] = "ls_HMOX1"
clus_out$cellstate_tk2[clus_out$integrated_snn_res.0.2 %in% c("3", "4")] = "ls_KIT"
clus_out$cellstate_tk2[clus_out$integrated_snn_res.0.2 %in% c("5")] = "ls_lactation"
clus_out$cellstate_tk2[clus_out$integrated_snn_res.0.2 %in% c("6")] = "ls_cd74"
clus_out$cellstate_tk2[clus_out$integrated_snn_res.0.2 %in% c("7")] = "ls_prol"

lvl = names(table(clus_out$cellstate_tk2))
clus_out$cellstate_tk2 = factor(clus_out$cellstate_tk2, levels = lvl)

write_rds(clus_out, "clus_out.rds")
clus_out = read_rds("clus_out.rds")
Idents(clus_out) <- "cellstate_tk2"

clus_out_seu = FindAllMarkers(object = clus_out, assay = "RNA", logfc.threshold = 0)
write.csv(clus_out_seu, "clus_out_seu.csv")


markers_v2 = clus_out_seu %>% dplyr::filter(gene %in% clus_out_seu$gene)

top_markers = markers_v2 %>%
  dplyr::group_by(cluster) %>%
  top_n(7, avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()

DefaultAssay(clus_out) <- "integrated"
p1= DoHeatmap(clus_out, features = c(top_markers), size = 2,raster = FALSE, label=TRUE,
              group.by = "cellstate", slot = "scale.data", disp.min = -1.5, disp.max = 1.5, group.colors = colors_dark) + 
  scale_color_manual(values = colors_dark) +
  scale_fill_gradient2(low="steelblue1", high = "tomato3", mid = "white", midpoint = 0) + 
  theme(axis.text.y.left = element_text(size=6))

tavg = AverageExpression(clus_out, features = unique(c("EPCAM", "KIT", "TOP2A", top_markers)), return.seurat = T)
DoHeatmap(tavg, features = unique(c("EPCAM", "KIT", "TOP2A", top_markers)), assay = "integrated", label = TRUE,
          slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = colors_dark, raster = F)  +
  scale_fill_gradient2(low = "#29A7B8", mid = "white", high = "#E85041", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5)) -> tcell_heatmap
png("heatmap_all_labelled_lumsec_cellstatetk2_v1.png", width=12, height=12, units = "in", res = 300)
tcell_heatmap
dev.off()


ggsave(plot = p1, filename = "heatmap_cellstate_Tcells_v2.png", device = "png", width = 20, height = 18)

png("heatmap_cellstate_Tcells_v1.png", width=18, height=12, units = "in", res = 300) #
print(p1)
dev.off()

png("umap_Tcells_v2.png", width=9, height=7, units = "in", res = 300) #
DimPlot(clus_out, group.by = "cellstate_tk2", cols = colors_dark, label = T)
dev.off()

png("dimplot_split_v1.png", width=18, height=24, units = "in", res = 300) #
DimPlot(clus_out, split.by = "sheet_patient_id", group.by = "cellstate_tk2", cols = colors_dark, label = F, ncol = 6, pt.size = 1)
dev.off()

png("dimplot_split_v2.png", width=12, height=12, units = "in", res = 300) #
DimPlot(clus_out, split.by = "sheet_exp_proc", group.by = "cellstate_tk2", cols = colors_dark, label = F, ncol = 2, pt.size = 1)
dev.off()

v1 = VlnPlot(tum, features = c("rna_CD74", "rna_CCL20", "rna_HMOX1", "rna_S100A8", 
                               "rna_KIT", "rna_PLCG2", "rna_PTN", "rna_SCGB2A2", 
                               "rna_KRT23", "rna_CXCL8", "rna_CXCL1", 
                               "rna_LALBA", "rna_CSN2", "rna_TOP2A", 
                               "rna_STMN1", "rna_CYP24A1", "rna_DUSP4"), group.by = "cellstate_tk2",
             pt.size = 0, cols = colors_dark, y.max = 3)

png("VlnPlot_all_labelled_rescellstate2_cells_v1.png", width=11, height=11, units = "in", res = 200)
v1
dev.off()

v2 = VlnPlot(tum, features = c("nCount_RNA", "nFeature_RNA"),pt.size = 0, 
             cols = colors_dark, group.by = "cellstate_tk2")

png("VlnPlot_count_gene_all_labelled_res_cellstate2_cells_v1.png", width=6, height=3, units = "in", res = 200)
v2
dev.off()


# barplots ------
# **Barplots ---------------------------------------------------------------
col_t = colors_dark
df_meta = data.frame(clus_out@meta.data, 
                     celltype = clus_out$cellstate_tk2, 
                     stringsAsFactors = F) %>% 
  as_tibble() %>% clean_names()

df_meta = df_meta %>% dplyr::mutate(sheet_bmi_final = case_when(sheet_bmi_final == "Normal" ~ "normal",
                                                                TRUE ~ as.character(sheet_bmi_final)))

df_meta1 = df_meta
# create a summary table
func_bar(var_df = df_meta1, var_x = "sheet_patient_id", var_fill = "cellstate_tk1", var_pos = "stack", var_cols = col_t, celltype = "cellstate_tk1")
func_bar(var_df = df_meta1, var_x = "sheet_patient_id", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "cellstate_tk2")
func_bar(var_df = df_meta1, var_x = "sheet_patient_id", var_fill = "sheet_parity", var_pos = "fill", var_cols = col_t, celltype = "sheet_parity")
func_bar(var_df = df_meta1, var_x = "sheet_patient_id", var_fill = "sheet_parity", var_pos = "fill", var_cols = col_t, celltype = "sheet_parity")
func_bar(var_df = df_meta1, var_x = "sheet_parity", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "sheet_parity")
func_bar(var_df = df_meta1, var_x = "sheet_parity", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "sheet_parity")

func_bar(var_df = df_meta1, var_x = "sheet_exp_proc", var_fill = "sheet_parity", var_pos = "fill", var_cols = colors_dark, celltype = "sheet_parity")

func_bar(var_df = df_meta1, var_x = "sheet_patient_id", var_fill = "sheet_exp_proc", var_pos = "fill", var_cols = colors_neon, celltype = "sheet_exp_proc")

func_bar(var_df = df_meta1, var_x = "final_sample_id", var_fill = "celltype", var_pos = "stack", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "final_sample_id", var_fill = "celltype", var_pos = "fill", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_bmi_final", var_fill = "celltype", var_pos = "stack", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_bmi_final", var_fill = "celltype", var_pos = "fill", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_brca_updated", var_fill = "celltype", var_pos = "stack", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_brca_updated", var_fill = "celltype", var_pos = "fill", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_exp_proc", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_exp_proc", var_fill = "celltype", var_pos = "fill", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_cancer_risk", var_fill = "celltype", var_pos = "stack", var_cols = col_t, celltype = "ultra_res2")
func_bar(var_df = df_meta1, var_x = "sheet_cancer_risk", var_fill = "celltype", var_pos = "fill", var_cols = col_t, celltype = "ultra_res2")
func_bar(var_df = df_meta1, var_x = "sheet_figure1_grp_updated", var_fill = "celltype", var_pos = "stack", var_cols = col_t, celltype = "sheet_figure1_grp_updated")
func_bar(var_df = df_meta1, var_x = "sheet_figure1_grp_updated", var_fill = "celltype", var_pos = "fill", var_cols = col_t, celltype = "sheet_figure1_grp_updated")
func_bar(var_df = df_meta1, var_x = "sheet_menopause", var_fill = "celltype", var_pos = "stack", var_cols = col_t, celltype = "sheet_menopause")
func_bar(var_df = df_meta1, var_x = "sheet_menopause", var_fill = "celltype", var_pos = "fill", var_cols = col_t, celltype = "sheet_menopause")
func_bar(var_df = df_meta1, var_x = "sheet_menopause", var_fill = "celltype", var_pos = "stack", var_cols = col_t, celltype = "sheet_menopause")
func_bar(var_df = df_meta1, var_x = "sheet_menopause", var_fill = "celltype", var_pos = "fill", var_cols = col_t, celltype = "sheet_menopause")
func_bar(var_df = df_meta1, var_x = "sheet_cancer_history", var_fill = "cellstate", var_pos = "stack", var_cols = col_t, celltype = "sheet_cancer_history")
func_bar(var_df = df_meta1, var_x = "sheet_cancer_history", var_fill = "cellstate", var_pos = "fill", var_cols = col_t, celltype = "sheet_cancer_history")


# check tks clusters ------
tk_obj = load("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_v10/hbca_seclum.RObj")
tk_obj = hbca.seclum

tk_obj_df = data.frame(cellname = tk_obj$cellname, tk_states = tk_obj$seclum_states)
clus_out_meta = clus_out@meta.data
clus_out_meta = left_join(clus_out_meta, tk_obj_df)

clus_out@meta.data = clus_out_meta
clus_out$tk_states[(is.na(clus_out$tk_states))] = "unkn"
rownames(clus_out@meta.data) = clus_out$cellname

png("umap_tk_Tcells_v2.png", width=9, height=7, units = "in", res = 300) #
DimPlot(clus_out, group.by = "tk_states", label = T, cols = colors_dark, pt.size = 0.0001)
dev.off()

png("umap_tk_split_by_pt_v2.png", width=18, height=18, units = "in", res = 300) #
DimPlot(clus_out, group.by = "tk_states", split.by = "sheet_patient_id", ncol = 6, label = T, cols = colors_dark, pt.size = 0.0001)
dev.off()

png("umap_tk_split_by_pt_v2.png", width=18, height=24, units = "in", res = 300) #
DimPlot(tk_obj, group.by = "seclum_states", split.by = "sheet_patient_id", ncol = 6, cols = colors_dark, pt.size = 1)
dev.off()

# confusion matrix ------
cc = table(clus_out$cellstate_tk2, clus_out$tk_states)
g1 = ggplot(as.data.frame(cc), aes(Var1,Var2, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  #scale_fill_viridis_c()+
  scale_fill_gradient2(low="white",high="#009194") +
  labs(x = "Tapsi's labels",y = "tk's labels")
png(paste0(save_path,"/tapsi_tk_confusionmat_v2.png"), width=8, height=6, units = "in", res = 300)
print(g1)
dev.off()

# pathways -------
# pathway ------
genes = read.csv("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_v2/clus_out_seu.csv")
genes = genes[,-1]
genes$cluster = as.character(genes$cluster)

gsea_resultsc2hall = fun_gsea(genes = genes, set = "H", sub_cat = NULL, num_avg_logFC = 0.25)
gsea_resultsc2C5 = fun_gsea(genes = genes, set = "C5", sub_cat = "BP", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "fibro")
gsea_resultsc2C3 = fun_gsea(genes = genes, set = "C3", sub_cat = "TFT:GTRD", num_avg_logFC = 0.25, celltype = "lumsec", doplot = FALSE, ret_collap = FALSE)
gsea_resultsc2C8 = fun_gsea(genes = genes, set = "C8", sub_cat = NULL, num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "fibro")
gsea_resultsc2C2_reac = fun_gsea(genes = genes, set = "C2", sub_cat = "CP:REACTOME", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "fibro")
gsea_resultsc2C2_kegg = fun_gsea(genes = genes, set = "C2", sub_cat = "CP:KEGG", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "lumsec")

clus = levels(clus_out$cellstate_tk2)

gs = gsea_resultsc2hall
gs = gsea_resultsc2C2;gs_kegg = gsea_resultsc2C2_kegg;gs_reac = gsea_resultsc2C2_reac
gs = gs_kegg
gs = gs_reac
gs = gsea_resultsc2C5
gs = gsea_resultsc2C3
gs = gsea_resultsc2C7
gs = gsea_resultsc2C8

gs1 = gs %>% dplyr::select(-leadingEdge) %>% dplyr::group_by(cluster) %>% 
  dplyr::filter(padj <= 0.05) %>%
  #top_n(30, wt = abs(NES)) %>% 
  pull(pathway)
gs2 = gs %>% dplyr::filter(pathway %in% gs1)
gs2$cluster = factor(gs2$cluster, levels = clus, ordered = TRUE)
gs2 = gs2 %>% dplyr::filter(grepl('HALLMARK|KEGG|REACTOME|BIOCARTA|GO', pathway))

#gs2 = gs2 %>% mutate(pathway = str_split(pathway, "_", n=2, simplify = TRUE)[,2]) # for C7
gs2 = gs2 %>%  arrange(cluster, -NES)

library(forcats);library(circlize);library(ComplexHeatmap)

p1 = gs2 %>% 
  dplyr::select(pathway, NES, cluster) %>%
  arrange(cluster, NES) %>%
  spread(key = cluster, value = NES, fill=0) 

p1.1 = p1
colnames(p1.1) = c("pathway", clus)

# rearrange the datafram
test = p1

# for basal
test1 = test %>% dplyr::arrange(-`ls_cd74`) %>% dplyr::filter(`ls_cd74`>0) %>% pull(pathway)
test2 = test %>% dplyr::filter(!pathway %in% test1) %>% arrange(-ls_HMOX1) %>% dplyr::filter(ls_HMOX1>0, `ls_cd74`<0, `ls_KIT`<0, `ls_KRT23`<0, `ls_lactation`<0,  `ls_prol`<0, `ls_SPINK5`<0) %>% dplyr::top_n(25, `ls_HMOX1`) %>% pull(pathway)
test3 = test %>% dplyr::filter(!pathway %in% c(test1, test2)) %>% arrange(-`ls_KIT`) %>% dplyr::filter(ls_HMOX1<0, `ls_cd74`<0, `ls_KIT`>0, `ls_KRT23`<0, `ls_lactation`<0,  `ls_prol`<0, `ls_SPINK5`<0) %>% dplyr::top_n(25, `ls_KIT`) %>% pull(pathway)
test4 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3)) %>% arrange(-`ls_KRT23`) %>% dplyr::filter(ls_HMOX1<0, `ls_cd74`<0, `ls_KIT`<0, `ls_KRT23`>0, `ls_lactation`<0,  `ls_prol`<0, `ls_SPINK5`<0) %>% dplyr::top_n(25, `ls_KRT23`) %>% pull(pathway)
test5 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3, test4)) %>% arrange(-`ls_lactation`) %>% dplyr::filter(ls_HMOX1<0, `ls_cd74`<0, `ls_KIT`<0, `ls_KRT23`<0, `ls_lactation`>0,  `ls_prol`<0, `ls_SPINK5`<0) %>% dplyr::top_n(25, `ls_lactation`) %>% pull(pathway)
test6 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3, test4, test5)) %>% arrange(-`ls_prol`) %>% dplyr::filter(ls_HMOX1<0, `ls_cd74`<0, `ls_KIT`<0, `ls_KRT23`<0, `ls_lactation`<0,  `ls_prol`>0, `ls_SPINK5`<0) %>% dplyr::top_n(25, `ls_prol`) %>% pull(pathway)
test7 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3, test4, test5, test6)) %>% arrange(-`ls_SPINK5`) %>% dplyr::filter(ls_HMOX1<0, ls_cd74<0, `ls_KIT`<0, `ls_KRT23`<0, `ls_lactation`<0,  `ls_prol`<0, `ls_SPINK5`>0) %>% dplyr::top_n(25, `ls_SPINK5`) %>% pull(pathway)
seq_path = c(test1, test2, test3,test4, test5, test6, test7) #


# make heatmap 
p1.1 = p1.1 %>%
  dplyr::slice(match(seq_path, pathway))
p1 = p1.1
p2 = data.matrix(p1)
rownames(p2) = p1$pathway
p2 = p2[,-1]
col = colorRamp2(c(-2, 0, 2),c("steelblue",  "white",  "#F80F61")) 
# c("steelblue","white", "tomato3"))
ht_hall = Heatmap(p2, 
                  show_row_names = T, "NES", 
                  col = col,
                  #km = 3,
                  #clustering_method_rows = "complete",
                  #clustering_distance_rows = "maximum",
                  row_names_gp = gpar(cex = 0.4),
                  cluster_rows = F,
                  row_dend_side = "right",
                  row_names_side = "left",
                  width = unit(5, "cm"), 
                  cluster_columns = F, row_names_max_width = unit(20, "cm"),
                  column_names_gp = gpar(cex = 0.6));ht_hall

png("htmaps_pathways_bp.png", res = 300, units = "in", width = 12, height = 22)
print(ht_hall)
dev.off()


