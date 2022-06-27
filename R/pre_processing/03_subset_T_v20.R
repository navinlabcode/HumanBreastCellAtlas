# geo ----
setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_T_v2";dir.create(odir)
options(future.globals.maxSize = 90000 * 1024^2)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "03_subset_T_v1"
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v2/"
col_t = c("#C70E7B", "#FC6882", "#007BC3", "#54BCD1", "#EF7C12", "#F4B95A", "#009F3F", "#8FDA04", "#AF6125", "#F4E3C7", "#B25D91", "#EFC7E6",
          "#EF7C12", "#F4B95A", colors_dark)

# core ----
setwd("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_T_v2"; dir.create(odir)
library(future)
plan("multiprocess", workers = 75)
options(future.globals.maxSize = Inf)
setwd(paste0("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "03_subset_T_v1"
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v2/"
col_t = c("#C70E7B", "#FC6882", "#007BC3", "#54BCD1", "#EF7C12", "#F4B95A", "#009F3F", "#8FDA04", "#AF6125", "#F4E3C7", "#B25D91", "#EFC7E6",
          "#EF7C12", "#F4B95A", colors_dark)


# read original object --------------
tum = readRDS(file = paste0(wd, data_folder,"/clus_out.rds"))

# subset ---------
Idents(tum) <- "cellstate"
tum_sub = subset(tum, idents = c("CD8_KLRC1", "CD8_low_quality", "CD8_pip", "CD8_SAA1", "HSPs","MMP_T"), invert = TRUE) #endo_vas
tum_sub1 = tum_sub

write_rds(tum_sub, path = glue("{save_path}/t_subset_v1.rds"))

# method 9 ----------
# Log T integration with no scaling
library(future)
plan("multiprocess", workers = 20)
#plan("sequential")

# ** Clus function --------
# run clustering function
seu_object = read_rds(glue("{save_path}/t_subset_v1.rds"))
DefaultAssay(seu_object) <- "RNA"
dims = 18; k.param = 20; cols = colors_dark; var_split = "sheet_exp_proc"; var_scale = "nCount_RNA"
fts1 = c("nFeature_RNA","nCount_RNA", "percent.mito","percent.rb", "percent.stress")
fts2 = c("rna_EPCAM","rna_PTPRC","rna_PECAM1","rna_DCN","rna_KRT5",  "rna_CD86", "rna_COL1A1", "rna_PIP")
fts3 = c("rna_CD8A", "rna_CD8B", "rna_CD4", "rna_CCR7", "rna_IL7R", "rna_NKG7", "rna_GZMB", "rna_GZMK")
fts4 = c("rna_MS4A1","rna_CD40LG","rna_FOXP3","rna_IFNG","rna_LAG3","rna_TIGIT",  "rna_NCAM1", "rna_KLRB1", "rna_CD3D")

ref_set = NULL

clus_out = func_cluster_vst(tum_int = seu_object, var_scale = var_scale, var_split = var_split, dims = dims, k.param = k.param, 
                            markers1 = fts1, markers2 = fts2, markers3 = fts3, markers4 = fts4, save_path = save_path,
                            cols = colors_dark, ref_set = ref_set)

char_ident = "integrated_snn_res.0.4"
u1 = DimPlot(clus_out, group.by = "integrated_snn_res.0.6", cols = col_t, label = T)
ggsave(plot = u1, filename = "figure_uamp_integrated_snn_res.0.6.png", device = "png", width = 8, height = 6)


# ** Markers function (9)--------
dims = 18; k.param = 20;
char_ident = "integrated_snn_res.0.8"  #integrated_snn_res.0.1
ht_fts1 = c("CD3D", "CD8A", "CD8B")
#clus_out = read_rds(path = glue("{save_path}/dim_{dims}k_{k.param}_tum.integrated.ns_SCT_v4_intgeratedassay_{char_ident}_vst_method9.rds")) 
mark_out = func_markers_vst(clus_out = clus_out, char_ident = char_ident, fts1 = ht_fts1, dims = dims, k.param = k.param, cols = colors_dark, 
                            save_path = save_path)

clus_out = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v2/dim_18k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")
markers = read.csv("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v2/dims_18k_20_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.6_vst_method_9.csv")
# clus_out = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v2/dim_35k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")
# markers = read.csv("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v1/dims_30k_20_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.4_vst_method_9.csv")
tum = clus_out

markers_v2 = markers %>% dplyr::filter(gene %in% markers$gene)

top_markers = markers_v2 %>%
  dplyr::group_by(cluster) %>%
  top_n(10, avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()

#Idents(tum.integrated.ns)
#tum.integrated.ns$clusters1 <- factor(tum.integrated.ns$clusters1, levels = levels(Idents(tum.integrated.ns)))
message("Making Heatmaps")
Idents(tum) <- "integrated_snn_res.0.6"
DefaultAssay(tum) <- "integrated"
p1 = DoHeatmap(tum, features = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B", top_markers), size = 2,raster = FALSE,
              group.by = char_ident, slot = "scale.data", disp.min = -1.5, disp.max = 1.5, group.colors = colors_dark) + 
  scale_color_manual(values = colors_dark) +
  scale_fill_gradient2(low="steelblue1", high = "tomato3", mid = "white", midpoint = 0) + 
  theme(axis.text.y.left = element_text(size=6))

tum$cellstate = factor(tum$cellstate, levels = levels(markers_v2$cluster))
Idents(tum) <- "cellstate"
tavg = AverageExpression(tum, features = unique(c("CD3D", "CD3E", "CD4", "CD8A", "CD8B", top_markers)), return.seurat = T)
a1 = DoHeatmap(tavg, features = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B", top_markers), assay = "integrated", label = TRUE,
          slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = col_t, raster = F)  +
  scale_fill_gradient2(low = "#29A7B8", mid = "white", high = "#E85041", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5)) -> tcell_heatmap

png("heatmap_all_labelled_integrated_snn_res.0.6_Tcells_v1.png", width=26, height=26, units = "in", res = 100)
p1
dev.off()
ggsave(plot = a1, filename = "figure_heatmap_cell_states_avg_integratedtop7.png", device = "png", width = 10, height = 14)
ggsave(plot = p1, filename = "figure_heatmap_integrated_snn_res.0.6g_integratedtop7.png", device = "png", width = 10, height = 14)

u1 = DimPlot(tum, group.by = "cellstate", cols = col_t, label = T)
ggsave(plot = u1, filename = "figure_uamp_cell_states.png", device = "png", width = 8, height = 6)

# test markers ----------------------------------------
Idents(tum) <- "integrated_snn_res.0.8"
naive_markers = c("rna_TCF7", "rna_SELL", "rna_LEF1", "rna_CCR7")
inhibitory_markers = c("rna_LAG3", "rna_TIGIT", "rna_PDCD1", "rna_HAVCR2", "rna_CTLA4")
effector_markers = c("rna_IL2", "rna_GZMA", "rna_GZMB", "rna_GZMK", "rna_GNLY", "rna_PRF1", "rna_IFNG", "rna_NKG7", "rna_CX3CR1")
co_stimulatory_markers = c("rna_CD24", "rna_TNFRSF14", "rna_ICOS", "rna_TNFRSF9")
treg_markers = c("rna_IL2RA", "rna_FOXP3", "rna_IKZF2")
tiss_res_mem_markers = c("rna_HAVCR2", "rna_GZMB", "rna_VCAM1", "rna_PRF1", "rna_KLRC1", "rna_CCL4", "rna_LAG3", "rna_CCL3")
tiss_eff_mem_markers = c("rna_KLRG1", "rna_LYAR", "rna_GZMK", "rna_GZMM", "rna_TXNIP", "rna_CD8B", "rna_FCRL6")
nkt = c("rna_KLRD1", "rna_PRF1", "rna_NKG7", "rna_KIR2DL3", "rna_NCR1", "rna_GZMH")
t = c("rna_CD3D", "rna_TRAC", "rna_CD28", "rna_IL7R", "rna_CD8B", "rna_GZMK")
tissue_resident = c("rna_IL2RB", "rna_PRDM1", "rna_ITGAE", "rna_ITGAM", "rna_ITGAX", "rna_CXCR6","rna_CX3CR1", "rna_CD101", "rna_CD69", "rna_ITGA1", "rna_PRDM1")
prol = c("rna_TOP2A", "rna_MKI67")

p1.1 = VlnPlot(tum, features = c("rna_CD3D", "rna_CD4", "rna_CD8A", "rna_CD8B", "rna_ITGAE", "rna_CD69","rna_PRDM1", "rna_GNLY", "rna_NKG7",
                                 "rna_CXCL13", "rna_FOXP3", "rna_CTLA4", "rna_IL7R", "rna_FOS", "rna_RGCC", "rna_GZMB", "rna_PRF1", "rna_STMN1", "rna_TUBB",
                                 "rna_KLRG1", "rna_GZMK", "rna_MKI67", "rna_TOP2A"),pt.size = 0)
p2 = VlnPlot(tum, features = c(t, nkt), pt.size = 0)
p3 = VlnPlot(tum, features = c(naive_markers, inhibitory_markers, effector_markers,co_stimulatory_markers), pt.size = 0)
p4 = VlnPlot(tum, features = c(treg_markers, tiss_res_mem_markers, tiss_eff_mem_markers), pt.size = 0)
p5 = VlnPlot(tum, features = c(tissue_resident), pt.size = 0)
p6 = VlnPlot(tum, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), group.by = "integrated_snn_res.0.8", pt.size = 0, cols = colors_dark)

p7 = VlnPlot(tum, features = c("rna_KRT14", "rna_GZMA", "rna_EOMES"), pt.size = 0, group.by = "integrated_snn_res.0.8")
p8 = VlnPlot(tum, features = c("rna_GZMB", "rna_GNLY", "rna_NKG7", "rna_IFNG", "rna_PRF1", "rna_CD4", "rna_STAT4", "rna_TBX21"), pt.size = 0)
p9 = VlnPlot(tum, features = c("rna_MX1", "rna_ISG15", "rna_IFIT3", "PDCD1"), pt.size = 0)

FeaturePlot(tum, features = c("rna_CD4", "rna_CD8A"))
FeaturePlot(tum, features = c("rna_IL7R", "rna_CD8A"))
FeaturePlot(tum, features = c("rna_CCR7", "rna_CD4"))
FeaturePlot(tum, features = c("rna_GNLY", "rna_GZMB", "rna_NKG7","rna_PRF1"))
FeaturePlot(tum, features = c("rna_CCL5", "rna_GZMK", "rna_GZMA", "rna_NKG7","rna_PRF1", "rna_CD4"))
FeaturePlot(tum, features = c("rna_GZMK", "rna_GZMA","rna_CD4"))
DimPlot(tum, group.by = "integrated_snn_res.1")
VlnPlot(tum, features = c("nCount_RNA", "nFeature_RNA"), pt.size = 0)
DimPlot(tum, group.by= "sheet_exp_proc")

vln = VlnPlot(tum, features = c("rna_KLF2", "rna_SELL", 
                          "rna_CCL20", "rna_MAF",
                          "rna_CD40LG", "rna_LMNA","rna_GZMK", "rna_GZMA",
                          "rna_MMP3", "rna_SOD2",
                          "rna_ZNF683", "rna_PLCG2",
                          "rna_ITGA1", "rna_PRF1",
                          "rna_CD74", "rna_IFNG",
                          "rna_GPCPD1", "rna_CCDC57", "rna_LITAF",
                          "rna_CD52", "rna_NKG7","rna_GNLY", "rna_FCGR3A", "rna_NCAM1"), pt.size = 0, group.by = "cellstate", cols = col_t)

ggsave(plot = vln, filename = "violin_markrer_cellstate", device = "png", width = 18, height = 16)
p_all = cowplot::plot_grid(p1.1, p2, p3, p4, p6, p7, p8, p9)
ggsave(plot = p_all, filename = "violin_integrated_snn_res.0.8.png", device = "png", width = 32, height = 36)

VlnPlot(tum, features = c("rna_TRAC", "rna_CXCL13", "rna_IFNG", "rna_CCL3", "rna_CCL5"), pt.size = 0)

VlnPlot(tum, "FCGR3A", pt.size = 0)

tum_cd4 = subset(tum, idents = c("1",  "3", "4", "6", "9", "13", "14") )
mks146 = FindAllMarkers(tum_cd4, logfc.threshold = 0);write.csv(mks146, "mks_cd4.csv")

tum_cd8 = subset(tum, idents = c("0",  "2", "8", "12", "10") )
mks_cd8 = FindAllMarkers(tum_cd8, logfc.threshold = 0);write.csv(mks_cd8, "mks_cd8.csv")

mks9_12 = FindMarkers(tum, ident.1 = "9", ident.2 = "12", logfc.threshold = 0);write.csv(mks9_12, "mks9_12.csv")

# ----------------------------------------

# assign clusters ------
clus_out$integrated_snn_res.0.4[is.na(clus_out$integrated_snn_res.0.4)] = "0"


clus_out$cellstate = "NULL"
clus_out$cellstate[clus_out$integrated_snn_res.0.8 %in% c("0")] = "CD8_prf1_cyto"
clus_out$cellstate[clus_out$integrated_snn_res.0.8 %in% c("1")] = "CD4_ccl20"
clus_out$cellstate[clus_out$integrated_snn_res.0.8 %in% c("2")] = "CD8_gzmk/a"
clus_out$cellstate[clus_out$integrated_snn_res.0.8 %in% c("3")] = "CD4_maf"
clus_out$cellstate[clus_out$integrated_snn_res.0.8 %in% c("4")] = "CD4_cd40lg"
clus_out$cellstate[clus_out$integrated_snn_res.0.8 %in% c("5")] = "NK_cd56low"
clus_out$cellstate[clus_out$integrated_snn_res.0.8 %in% c("6")] = "CD4_naive"
clus_out$cellstate[clus_out$integrated_snn_res.0.8 %in% c("7")] = "NKT_cd16low"
clus_out$cellstate[clus_out$integrated_snn_res.0.8 %in% c("8")] = "CD4_CD8_znf683"
clus_out$cellstate[clus_out$integrated_snn_res.0.8 %in% c("9")] = "CD4_gzmk/a"
clus_out$cellstate[clus_out$integrated_snn_res.0.8 %in% c("10")] = "CD8_ifng"
clus_out$cellstate[clus_out$integrated_snn_res.0.8 %in% c("11")] = "NK_cd56high"
clus_out$cellstate[clus_out$integrated_snn_res.0.8 %in% c("12")] = "CD4_CD8_gzmk/a"
clus_out$cellstate[clus_out$integrated_snn_res.0.8 %in% c("13")] = "CD4_treg"
clus_out$cellstate[clus_out$integrated_snn_res.0.8 %in% c("14")] = "CD4_mmp3"
clus_out$cellstate[clus_out$integrated_snn_res.0.8 %in% c("15")] = "NKT_cd16high"
clus_out$cellstate[clus_out$integrated_snn_res.0.8 %in% c("16")] = "CD8_gd"
clus_out$cellstate[clus_out$integrated_snn_res.0.8 %in% c("17")] = "prol"

clus_out$cellstate = factor(clus_out$cellstate, c("CD4_naive", "CD4_treg", "CD4_ccl20", "CD4_maf", "CD4_cd40lg", "CD4_gzmk/a", "CD4_mmp3", "CD4_CD8_gzmk/a",
                                                  "CD4_CD8_znf683", "CD8_prf1_cyto", "CD8_gzmk/a", "CD8_ifng", "CD8_gd", "NKT_cd16low", "NKT_cd16high",
                                                  "NK_cd56low", "NK_cd56high", "prol"))
lvl = names(table(clus_out$cellstate))
#clus_out$cellstate = factor(clus_out$cellstate, levels = lvl)

write_rds(clus_out, "clus_out_updated.rds")
clus_out = read_rds("clus_out_updated.rds")
Idents(clus_out) <- "cellstate"

clus_out_seu = FindAllMarkers(object = clus_out, assay = "RNA", logfc.threshold = 0)
write.csv(clus_out_seu, "clus_out_updated.csv")

clus_out_seu = read.csv("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v2/dims_18k_20_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.8_vst_method_9_copy.csv")
markers_v2 = clus_out_seu %>% dplyr::filter(gene %in% clus_out_seu$gene)
markers_v2$cluster = factor(markers_v2$cluster, levels = lvl)

top_markers = markers_v2 %>% arrange(cluster) %>% 
  dplyr::group_by(cluster) %>%
  top_n(7, avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()

DefaultAssay(clus_out) <- "integrated"
tum_sub1 = AverageExpression(clus_out, features = unique(top_markers), return.seurat = T)
levels(tum_sub1$orig.ident) <- levels(Idents(clus_out))

p1 = DoHeatmap(tum_sub1, features = top_markers, assay = "integrated", label = TRUE,
              slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = col_t, raster = F)  +
  scale_fill_gradient2(low = "#29A7B8", mid = "white", high = "#E85041", na.value = "00000000")

p1 = DoHeatmap(clus_out, features = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B", top_markers), size = 2,raster = FALSE, label=TRUE,
              group.by = "cellstate", slot = "scale.data", disp.min = -1.5, disp.max = 1.5, group.colors = colors_dark) + 
  scale_color_manual(values = colors_dark) +
  scale_fill_gradient2(low="steelblue1", high = "tomato3", mid = "white", midpoint = 0) + 
  theme(axis.text.y.left = element_text(size=6))

ggsave(plot = p1, filename = "heatmap_cellstate_Tcells_v2.png", device = "png", width = 20, height = 18)

png("heatmap_cellstate_Tcells_v1.png", width=18, height=12, units = "in", res = 300) #
print(p1)
dev.off()

png("umap_Tcells_v1.png", width=9, height=7, units = "in", res = 300) #
DimPlot(clus_out, group.by = "cellstate", cols = col_t, label = T)
dev.off()

# barplots ------
# **Barplots ---------------------------------------------------------------
df_meta = data.frame(clus_out@meta.data, 
                     celltype = clus_out$cellstate, 
                     stringsAsFactors = F) %>% 
  as_tibble() %>% clean_names()

df_meta = df_meta %>% dplyr::mutate(sheet_bmi_final = case_when(sheet_bmi_final == "Normal" ~ "normal",
                                                                TRUE ~ as.character(sheet_bmi_final)))

df_meta1 = df_meta
# create a summary table
func_bar(var_df = df_meta1, var_x = "sheet_patient_id", var_fill = "cellstate", var_pos = "stack", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_patient_id", var_fill = "cellstate", var_pos = "fill", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "final_sample_id", var_fill = "cellstate", var_pos = "stack", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "final_sample_id", var_fill = "cellstate", var_pos = "fill", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_bmi_final", var_fill = "cellstate", var_pos = "stack", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_bmi_final", var_fill = "cellstate", var_pos = "fill", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_brca_updated", var_fill = "cellstate", var_pos = "stack", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_brca_updated", var_fill = "cellstate", var_pos = "fill", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_exp_proc", var_fill = "cellstate", var_pos = "stack", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_exp_proc", var_fill = "cellstate", var_pos = "fill", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_cancer_risk", var_fill = "cellstate", var_pos = "stack", var_cols = col_t, celltype = "ultra_res2")
func_bar(var_df = df_meta1, var_x = "sheet_cancer_risk", var_fill = "cellstate", var_pos = "fill", var_cols = col_t, celltype = "ultra_res2")
func_bar(var_df = df_meta1, var_x = "sheet_figure1_grp_updated", var_fill = "cellstate", var_pos = "stack", var_cols = col_t, celltype = "sheet_figure1_grp_updated")
func_bar(var_df = df_meta1, var_x = "sheet_figure1_grp_updated", var_fill = "cellstate", var_pos = "fill", var_cols = col_t, celltype = "sheet_figure1_grp_updated")
func_bar(var_df = df_meta1, var_x = "sheet_menopause", var_fill = "cellstate", var_pos = "stack", var_cols = col_t, celltype = "sheet_menopause")
func_bar(var_df = df_meta1, var_x = "sheet_menopause", var_fill = "cellstate", var_pos = "fill", var_cols = col_t, celltype = "sheet_menopause")
func_bar(var_df = df_meta1, var_x = "sheet_menopause", var_fill = "cellstate", var_pos = "stack", var_cols = col_t, celltype = "sheet_menopause")
func_bar(var_df = df_meta1, var_x = "sheet_menopause", var_fill = "cellstate", var_pos = "fill", var_cols = col_t, celltype = "sheet_menopause")
func_bar(var_df = df_meta1, var_x = "sheet_cancer_history", var_fill = "cellstate", var_pos = "stack", var_cols = col_t, celltype = "sheet_cancer_history")
func_bar(var_df = df_meta1, var_x = "sheet_cancer_history", var_fill = "cellstate", var_pos = "fill", var_cols = col_t, celltype = "sheet_cancer_history")

# compare with kerrigan --------
#clus_out = read_rds("clus_out.rds")
#tum = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v2/dim_12k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")
meta_df = read.csv("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/HCA_immune_metadata_cell_labels.txt",sep = "\t")
meta_df$cellname = rownames(meta_df)
meta_df_b = meta_df %>% dplyr::filter(final_group == "t")

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


# create a confusion matrix

cc = table(tum2$cellstate, tum2$second_pass_cluster_names)
g1 = ggplot(as.data.frame(cc), aes(Var1,Var2, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  #scale_fill_viridis_c()+
  scale_fill_gradient2(low="white",high="#009194") +
  labs(x = "Tapsi's labels",y = "Kerri's labels")
png(paste0(save_path,"/tapsi_kerri_confusionmat_v2.png"), width=8, height=6, units = "in", res = 300)
print(g1)
dev.off()


