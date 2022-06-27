# geo ----
setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_myeloid_v2";dir.create(odir)
options(future.globals.maxSize = 90000 * 1024^2)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "03_subset_myeloid_v1"
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_myeloid_v2/"
col_m = c("#f9f871", "#00c9ae", "#ff8863", "#6e3452", "#00aafd", "#b178de", "#ff86da", "#7d9300","red2","#C70E7B", "#FC6882", "#A6E000", "#1BB6AF", "#6C6C9D", "#172869", "#D64358", "#EAFB88", "#3C8C4D", "#DFCEE0","orange", "red", "darkorchid1","darkorchid4", colors_dark)

# core ----
setwd("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_myeloid_v2"; dir.create(odir)
library(future)
plan("multiprocess", workers = 100)
options(future.globals.maxSize = Inf)
setwd(paste0("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "03_subset_myeloid_v1"
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_myeloid_v2/"
col_m = c("#f9f871", "#00c9ae", "#ff8863", "#6e3452", "#00aafd", "#b178de", "#ff86da", "#7d9300","red2","#C70E7B", "#FC6882", "#A6E000", "#1BB6AF", "#6C6C9D", "#172869", "#D64358", "#EAFB88", "#3C8C4D", "#DFCEE0","orange", "red", "darkorchid1","darkorchid4", colors_dark)


# read original object --------------
seu_object = readRDS(file = paste0(wd, data_folder,"/tum_sub1.rds"))

# subset ---------
Idents(seu_object) <- "cellstate"
tum_sub = subset(seu_object, idents = c("ITGAX_XIST", "macro_HSP"), invert = TRUE) #endo_vas
tum_sub1 = tum_sub
tum_sub@meta.data$cellstate = tum_sub$cellstate
write_rds(tum_sub, path = glue("{save_path}/mye_subset_v1.rds"))


# method 9 ----------
# Log T integration with no scaling
library(future)
plan("multiprocess", workers = 50)
#plan("sequential")

# ** Clus function --------
# run clustering function
seu_object = read_rds(glue("{save_path}/t_subset_v1.rds"))
DefaultAssay(seu_object) <- "RNA"
dims = 30; k.param = 20; cols = colors_dark; var_split = "sheet_exp_proc"; var_scale = "nCount_RNA"
fts1 = c("nFeature_RNA","nCount_RNA", "percent.mito","percent.rb", "percent.stress")
fts2 = c("rna_EPCAM","rna_PTPRC","rna_PECAM1","rna_DCN","rna_MKI67",  "rna_CD8A", "rna_CD4", "rna_CD3D")
fts3 = c("rna_CD14", "rna_CXCL13","rna_LY6G","rna_S100A8", "rna_CD68", "rna_CXCL2", "rna_ISG15", "rna_CD86")
fts4 = c("rna_NOS2", "rna_MRC1", "rna_ITGAM", "rna_ARG1", "rna_C1QA", "rna_APOE", "rna_ARG1", "rna_CD33", "rna_CD68", "rna_CD16", "rna_CD163", 
         "rna_ITGAX", "rna_FCGR3A", "rna_CD1C", "rna_IL3RA", "rna_ALOX5AP")

ref_set = NULL

clus_out = func_cluster_vst(tum_int = seu_object, var_scale = var_scale, var_split = var_split, dims = dims, k.param = k.param, 
                            markers1 = fts1, markers2 = fts2, markers3 = fts3, markers4 = fts4, save_path = save_path,
                            cols = colors_dark, ref_set = ref_set)

char_ident = "integrated_snn_res.0.4"

clus_out = read_rds("~/projects/bcap/analysis/20201212_final2/03_subset_myeloid_v2/dim_30k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")
png("umap_myecells_integrated_snn_res.0.4.png", width=9, height=7, units = "in", res = 300) #
DimPlot(clus_out, group.by = "integrated_snn_res.0.4", cols = col_m, label = T)
DimPlot(clus_out, group.by = "cellstate", cols = col_m, label = T)
dev.off()


# ** Markers function (9)--------
dims = 30; k.param = 20;
char_ident = "integrated_snn_res.0.4"  #integrated_snn_res.0.1
ht_fts1 = c("PTPRC", "CD68", "CD86", "CD14", "TPSAB1", "CSF3R")

#clus_out = read_rds(path = glue("{save_path}/dim_{dims}k_{k.param}_tum.integrated.ns_SCT_v4_intgeratedassay_{char_ident}_vst_method9.rds")) 
mark_out = func_markers_vst(clus_out = clus_out, char_ident = char_ident, fts1 = ht_fts1, dims = dims, k.param = k.param, cols = colors_dark, 
                            save_path = save_path)

mks0_4 = FindMarkers(tum1, ident.1 = "0", "4", group.by = "integrated_snn_res.0.4")
write.csv(mks0_4, "mks0_4.csv")
mks1_7 = FindMarkers(tum1, ident.1 = "1", "7", group.by = "integrated_snn_res.0.4")
write.csv(mks1_7, "mks1_7.csv")

VlnPlot(tum1, features = "rna_CD14", group.by = "integrated_snn_res.0.4", pt.size = 0)

#clus_out = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_myeloid_v1/dim_30k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")
#clus_out = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_myeloid_v1/dim_30k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")
#markers = read.csv("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_myeloid_v1/dims_30k_20_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.4_vst_method_9.csv")
#markers = read.csv("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_myeloid_v1/dims_30k_20_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.4_vst_method_9.csv")
tum = clus_out
markers_v2 = markers %>% dplyr::filter(gene %in% markers$gene)

#tum.integrated.ns = NormalizeData(tum.integrated.ns, assay = "RNA")
#levels(markers_v2$cluster) = c("M1", "M2", "M3")
#markers_v2 = markers_v2 %>% dplyr::mutate(cluster_1 = factor(x = cluster, levels = c("M1", "M2", "M3"), ordered = TRUE))
#markers_v2 = markers_v2 %>% dplyr::mutate(cluster = factor(x = cluster, levels = unique(tum.integrated.ns$clusters1), ordered = TRUE))
top_markers = markers_v2 %>%
  #mutate(cluster_1 = factor(cluster, levels = c("M1", "M2", "M3"))) %>%
  dplyr::group_by(cluster) %>%
  #arrange(avg_logFC) %>% 
  top_n(7, avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()

#Idents(tum.integrated.ns)
#tum.integrated.ns$clusters1 <- factor(tum.integrated.ns$clusters1, levels = levels(Idents(tum.integrated.ns)))
message("Making Heatmaps")
DefaultAssay(tum) <- "integrated"
p1= DoHeatmap(tum, features = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B", top_markers), size = 2,raster = FALSE,
              group.by = char_ident, slot = "scale.data", disp.min = -1.5, disp.max = 1.5, group.colors = colors_dark) + 
  scale_color_manual(values = colors_dark) +
  scale_fill_gradient2(low="steelblue1", high = "tomato3", mid = "white", midpoint = 0) + 
  theme(axis.text.y.left = element_text(size=6))

png("heatmap_all_labelled_Tcells_v1.png", width=26, height=26, units = "in", res = 300)
p1
dev.off()
ggsave(plot = p1, filename = "figure_heatmap_cell_types_avg_integratedtop7.png", device = "png", width = 14, height = 10)


# manual checks ------
#tum = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_lymphendo_v1/dim_12k_20_nCount_RNA_tum.integrated.ns2_method_9.rds")

# filter other cells --------
Idents(tum1) <- "integrated_snn_res.0.4"

# filter again -----
tum_sub = subset(tum1, idents = c("14"), invert  = TRUE)
write_rds(tum_sub, "tum_sub_v2.rds")





# DimPlot(clus_out, group.by = "integrated_snn_res.0.4", label = T, split.by = "sheet_exp_proc", pt.size = 0.7, cols = colors_dark)
# DimPlot(clus_out, group.by = "integrated_snn_res.0.4", label = T, pt.size = 0.3, cols = colors_dark)
# VlnPlot(clus_out, group.by = "integrated_snn_res.0.4", 
#         features = c("rna_MS4A1","rna_CD79A", "rna_PTPRC", "rna_CD3D", "rna_EPCAM", "rna_CD68", "rna_KRT8", 
#                      "rna_COL1A1", "rna_RGS5", "rna_KRT5", "rna_NKG7", "rna_DCN", "rna_MKI67"), pt.size = 0, cols = colors_dark)
# 
# VlnPlot(clus_out, group.by = "integrated_snn_res.0.4", 
#         features = c("rna_CD1A", "rna_CLEC9A", "rna_IRF8", "rna_CD14", "rna_TPSAB1", "rna_CSF3R"), pt.size = 0)
# VlnPlot(clus_out, group.by = "integrated_snn_res.0.4", 
#         features = c("rna_CSF3R", "rna_FCGR3A", "rna_CXCL8", "rna_PI3", "rna_TPSAB1", "rna_THBS1"), pt.size = 0)
# VlnPlot(clus_out, group.by = "integrated_snn_res.0.4", 
#         features = c("rna_VCAN", "rna_CD52", "rna_C1QC", "rna_CD14", "rna_FCN1", "rna_S100A9"), pt.size = 0)
# 
# 
# reduction = "umap"
# naive_markers = c("rna_TCF7", "rna_SELL", "rna_LEF1", "rna_CCR7")
# inhibitory_markers = c("rna_LAG3", "rna_TIGIT", "rna_PDCD1", "rna_HAVCR2", "rna_CTLA4")
# effector_markers = c("rna_IL2", "rna_GZMA", "rna_GZMB", "rna_GZMK", "rna_GNLY", "rna_PRF1", "rna_IFNG", "rna_NKG7", "rna_CX3CR1")
# co_stimulatory_markers = c("rna_CD24", "rna_TNFRSF14", "rna_ICOS", "rna_TNFRSF9")
# treg_markers = c("rna_IL2RA", "rna_FOXP3", "rna_IKZF2")
# tiss_res_mem_markers = c("rna_HAVCR2", "rna_GZMB", "rna_VCAM1", "rna_PRF1", "rna_KLRC1", "rna_CCL4", "rna_LAG3", "rna_CCL3")
# tiss_eff_mem_markers = c("rna_KLRG1", "rna_LYAR", "rna_GZMK", "rna_GZMM", "rna_TXNIP", "rna_CD8B", "rna_FCRL6")
# nkt = c("rna_KLRD1", "rna_PRF1", "rna_NKG7", "rna_KIR2DL3", "rna_NCR1", "rna_GZMH")
# t = c("rna_CD3D", "rna_TRAC", "rna_CD28", "rna_IL7R", "rna_CD8B", "rna_GZMK")
# tissue_resident = c("rna_IL2RB", "rna_PRDM1", "rna_ITGAE", "rna_ITGAM", "rna_ITGAX", "rna_CXCR6","rna_CX3CR1", "rna_CD101", "rna_CD69", "rna_ITGA1", "rna_PRDM1")
# prol = c("rna_TOP2A", "rna_MKI67")
# 
# p1.1 = VlnPlot(tum, reduction = reduction, features = c("rna_CD3D", "rna_CD4", "rna_CD8A", "rna_CD8B", "rna_ITGAE", "rna_CD69","rna_PRDM1", "rna_GNLY", "rna_NKG7",
#                                                         "rna_CXCL13", "rna_FOXP3", "rna_CTLA4", "rna_IL7R", "rna_FOS", "rna_RGCC", "rna_GZMB", "rna_PRF1", "rna_STMN1", "rna_TUBB",
#                                                         "rna_KLRG1", "rna_GZMK", "rna_MKI67", "rna_TOP2A"))
# p2 = VlnPlot(tum, features = c(t, nkt), pt.size = 0)
# p3 = VlnPlot(tum, features = c(naive_markers, inhibitory_markers, effector_markers,co_stimulatory_markers), pt.size = 0)
# p4 = VlnPlot(tum, features = c(treg_markers, tiss_res_mem_markers, tiss_eff_mem_markers), pt.size = 0)
# p5 = VlnPlot(tum, features = c(tissue_resident), pt.size = 0)
# p6=VlnPlot(tum, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), group.by = "integrated_snn_res.0.4", pt.size = 0, cols = col_m)
# 
# VlnPlot(tum, features = c("rna_GZMK", "rna_GZMA", "rna_EOMES"), pt.size = 0)
# VlnPlot(tum, features = c("rna_GZMB", "rna_GNLY", "rna_NKG7", "rna_IFNG", "rna_PRF1", "rna_CD4", "rna_STAT4", "rna_TBX21"), pt.size = 0)
# VlnPlot(tum, features = c("rna_MX1", "rna_ISG15", "rna_IFIT3", "PDCD1"), pt.size = 0)
# 
# png(paste0(save_path,"/", reduction,dims,"k_", k.param,  "umap_all_labelled_Tcells.png"), width=16, height=16, units = "in", res = 300)
# print(p1.1)
# dev.off()
# 
# 
# png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_all_markers.png"), width=12, height=9, units = "in", res = 300)
# print(p2)
# dev.off()
# 
# png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_all_markers2.png"), width=12, height=14, units = "in", res = 300)
# print(p3)
# dev.off()
# 
# png(paste0(save_path,"/", reduction,dims,"k_", k.param,  "umap_all_markers3.png"), width=12, height=12, units = "in", res = 300)
# print(p4)
# dev.off()
# 
# png(paste0(save_path,"/", reduction,dims,"k_", k.param, "tissue_resident.png"), width=16, height=12, units = "in", res = 300)
# print(p5)
# dev.off()
# 
# png(paste0(save_path,"/", reduction,dims,"k_", k.param, "tissue_resident.png"), width=16, height=12, units = "in", res = 300)
# print(p5)
# dev.off()

# assign clusters ------
#tum_sub$integrated_snn_res.0.4[is.na(tum_sub$integrated_snn_res.0.4)] = "0"


tum_sub$cellstate = "NULL"
tum_sub$cellstate[tum_sub$integrated_snn_res.0.4 %in% c("0")] = "macro_ccl_m2"
tum_sub$cellstate[tum_sub$integrated_snn_res.0.4 %in% c("1")] = "macro_c3_apoe_m1"
tum_sub$cellstate[tum_sub$integrated_snn_res.0.4 %in% c("2")] = "cDC2_cd1c"
tum_sub$cellstate[tum_sub$integrated_snn_res.0.4 %in% c("3")] = "mono_cxcl5_cd14"
tum_sub$cellstate[tum_sub$integrated_snn_res.0.4 %in% c("4")] = "macro_macro_m2"
tum_sub$cellstate[tum_sub$integrated_snn_res.0.4 %in% c("5")] = "macro_apoc1_foam"
tum_sub$cellstate[tum_sub$integrated_snn_res.0.4 %in% c("6")] = "mono_cd52_cd16"
tum_sub$cellstate[tum_sub$integrated_snn_res.0.4 %in% c("7")] = "macro_c3_plcg2"
tum_sub$cellstate[tum_sub$integrated_snn_res.0.4 %in% c("8")] = "prol"
tum_sub$cellstate[tum_sub$integrated_snn_res.0.4 %in% c("9")] = "cDC1_clecl9a"
tum_sub$cellstate[tum_sub$integrated_snn_res.0.4 %in% c("10")] = "mDC"
tum_sub$cellstate[tum_sub$integrated_snn_res.0.4 %in% c("11")] = "macro_ifn"
tum_sub$cellstate[tum_sub$integrated_snn_res.0.4 %in% c("12")] = "mast"
tum_sub$cellstate[tum_sub$integrated_snn_res.0.4 %in% c("13")] = "pDC"

lvl = names(table(tum_sub$cellstate))
tum_sub$cellstate = factor(tum_sub$cellstate, levels = lvl)

DimPlot(tum_sub, group.by = "cellstate", label = T, pt.size = 0.3, cols = col_m)

# FINAL OBJECT ------
write_rds(tum_sub, "tum_sub2.rds")
tum_sub = read_rds("tum_sub2.rds")
Idents(tum_sub) <- "cellstate"

tum_sub_seu = FindAllMarkers(object = tum_sub, assay = "RNA", logfc.threshold = 0)
write.csv(tum_sub_seu, "tum_sub_cellstate2.csv")

tum_sub_seu = read.csv("tum_sub_cellstate2.csv")
markers_v2 = tum_sub_seu %>% dplyr::filter(gene %in% tum_sub_seu$gene)

top_markers = markers_v2 %>%
  dplyr::group_by(cluster) %>%
  top_n(7, avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% unique() %>% as.character()

DefaultAssay(tum_sub) <- "integrated"
Idents(tum_sub) <- "cellstate"

tum_sub1 = AverageExpression(tum_sub, features = top_markers, return.seurat = T)
levels(tum_sub1$orig.ident) <- levels(Idents(tum_sub))
p1= DoHeatmap(tum_sub1, features = c("CD68", "CD86", top_markers), size = 2,raster = FALSE, label=TRUE,
              group.by = "orig.ident", slot = "scale.data", disp.min = -1.5, disp.max = 1.5, group.colors = col_m) + 
  scale_color_manual(values = col_m) +
  scale_fill_gradient2(low="steelblue1", high = "tomato3", mid = "white", midpoint = 0) + 
  theme(axis.text.y.left = element_text(size=6))
p1= DoHeatmap(tum_sub1, features = top_markers, assay = "integrated", label = TRUE,
          slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = col_m, raster = F)  +
  scale_fill_gradient2(low = "#29A7B8", mid = "white", high = "#E85041", na.value = "00000000")
#ggsave(plot = p1, filename = "heatmap_cellstate_myecells_v2.png", device = "png", width = 20, height = 18)

png("heatmap_cellstate_myecells_v2.png", width=8, height=12, units = "in", res = 300) #
print(p1)
dev.off()

png("umap_myeloid_v2.png", width=9, height=7, units = "in", res = 300) #
DimPlot(tum_sub, group.by = "cellstate", cols = col_m, label = T)
dev.off()
png("umap_myeloid_v2.png", width=9, height=7, units = "in", res = 300) #
DimPlot(tum_sub, group.by = "cellstate", cols = col_m, label = F)
dev.off()

# barplots ------
# **Barplots ---------------------------------------------------------------
col_m = c("#f9f871", "#00c9ae", "#ff8863", "#6e3452", "#00aafd", "#b178de", "#ff86da", "#7d9300","red2","#C70E7B", "#FC6882", "#A6E000", "#1BB6AF", "#6C6C9D", "#172869", "#D64358", "#EAFB88", "#3C8C4D", "#DFCEE0","orange", "red", "darkorchid1","darkorchid4", colors_dark)
df_meta = data.frame(tum_sub@meta.data, 
                     celltype = tum_sub$cellstate, 
                     stringsAsFactors = F) %>% 
  as_tibble() %>% clean_names()

df_meta = df_meta %>% dplyr::mutate(sheet_bmi_final = case_when(sheet_bmi_final == "Normal" ~ "normal",
                                                                TRUE ~ as.character(sheet_bmi_final)))

df_meta1 = df_meta
# create a summary table
func_bar(var_df = df_meta1, var_x = "sheet_patient_id", var_fill = "cellstate", var_pos = "stack", var_cols = col_m, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_patient_id", var_fill = "cellstate", var_pos = "fill", var_cols = col_m, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "final_sample_id", var_fill = "cellstate", var_pos = "stack", var_cols = col_m, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "final_sample_id", var_fill = "cellstate", var_pos = "fill", var_cols = col_m, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_bmi_final", var_fill = "cellstate", var_pos = "stack", var_cols = col_m, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_bmi_final", var_fill = "cellstate", var_pos = "fill", var_cols = col_m, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_brca_updated", var_fill = "cellstate", var_pos = "stack", var_cols = col_m, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_brca_updated", var_fill = "cellstate", var_pos = "fill", var_cols = col_m, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_exp_proc", var_fill = "cellstate", var_pos = "stack", var_cols = col_m, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_exp_proc", var_fill = "cellstate", var_pos = "fill", var_cols = col_m, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_cancer_risk", var_fill = "cellstate", var_pos = "stack", var_cols = col_m, celltype = "ultra_res2")
func_bar(var_df = df_meta1, var_x = "sheet_cancer_risk", var_fill = "cellstate", var_pos = "fill", var_cols = col_m, celltype = "ultra_res2")
func_bar(var_df = df_meta1, var_x = "sheet_figure1_grp_updated", var_fill = "cellstate", var_pos = "stack", var_cols = col_m, celltype = "sheet_figure1_grp_updated")
func_bar(var_df = df_meta1, var_x = "sheet_figure1_grp_updated", var_fill = "cellstate", var_pos = "fill", var_cols = col_m, celltype = "sheet_figure1_grp_updated")
func_bar(var_df = df_meta1, var_x = "sheet_menopause", var_fill = "cellstate", var_pos = "stack", var_cols = col_m, celltype = "sheet_menopause")
func_bar(var_df = df_meta1, var_x = "sheet_menopause", var_fill = "cellstate", var_pos = "fill", var_cols = col_m, celltype = "sheet_menopause")
func_bar(var_df = df_meta1, var_x = "sheet_menopause", var_fill = "cellstate", var_pos = "stack", var_cols = col_m, celltype = "sheet_menopause")
func_bar(var_df = df_meta1, var_x = "sheet_menopause", var_fill = "cellstate", var_pos = "fill", var_cols = col_m, celltype = "sheet_menopause")
func_bar(var_df = df_meta1, var_x = "sheet_cancer_history", var_fill = "cellstate", var_pos = "stack", var_cols = col_m, celltype = "sheet_cancer_history")
func_bar(var_df = df_meta1, var_x = "sheet_cancer_history", var_fill = "cellstate", var_pos = "fill", var_cols = col_m, celltype = "sheet_cancer_history")

meta_df = read.csv("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/HCA_immune_metadata_cell_labels.txt",sep = "\t")
meta_df$cellname = rownames(meta_df)
meta_df_b = meta_df %>% dplyr::filter(final_group == "macro")

tum1 = clus_out
meta_df_orig = tum1@meta.data
head(meta_df_orig)

meta_df_mrg = left_join(meta_df_orig, meta_df)
tum1@meta.data = meta_df_mrg
rownames(tum1@meta.data) = tum1@meta.data$cellname
tum1@meta.data$second_pass_cluster_names = as.character(tum1@meta.data$second_pass_cluster_names)
DimPlot(tum1, group.by = "second_pass_cluster_names", cols = colors_dark, label = T)
DimPlot(tum1, group.by = "cellstate", cols = colors_dark, label = T)
DimPlot(tum1, group.by = "cellstate", cols = colors_dark, label = T, split.by = "sheet_exp_proc")
DimPlot(tum1, group.by = "integrated_snn_res.0.3", cols = colors_dark, label = T)
table(tum1@meta.data$second_pass_cluster_names)

tum2 = tum
meta_df_orig = tum2@meta.data
head(meta_df_orig)

meta_df_mrg = left_join(meta_df_orig, meta_df_b)
tum2@meta.data = meta_df_mrg
rownames(tum2@meta.data) = tum2@meta.data$cellname
tum2@meta.data$second_pass_cluster_names = as.character(tum2@meta.data$second_pass_cluster_names)
DimPlot(tum2, group.by = "second_pass_cluster_names", cols = colors_dark, label = T)
DimPlot(tum2, group.by = "cellstate", cols = colors_dark, label = T)
DimPlot(tum2, group.by = "integrated_snn_res.0.3", cols = colors_dark, label = T)
