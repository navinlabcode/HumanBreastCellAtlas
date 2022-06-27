# geo ----
setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_T_v1";dir.create(odir)
options(future.globals.maxSize = 90000 * 1024^2)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "figure_1_v1"
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v1/"

# core ----
setwd("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_T_v1"; dir.create(odir)
library(future)
plan("multiprocess", workers = 75)
options(future.globals.maxSize = Inf)
setwd(paste0("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "figure_1_v1"
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v1/"


# read original object --------------
tum = readRDS(file = paste0(wd, data_folder,"/final_use.rds"))

# subset ---------
Idents(tum) <- "final_group1"
tum_sub = subset(tum, idents = "T-cells") #endo_vas
tum_sub_cells = WhichCells(object = tum, idents = "T-cells") #subset epi_lumhr cells #83891 cells
tum_sub1 = tum_sub

# Idents(tum) <- "final_group1"
# tum_sub = subset(tum, idents = "Myeloid") #endo_vas
# tum_sub_cells = WhichCells(object = tum, idents = "Myeloid") #subset epi_lumhr cells #83891 cells
# write_rds(tum_sub, path = glue("{save_path}/myeloid_subset_v1.rds"))
# 
# Idents(tum) <- "final_group1"
# tum_sub = subset(tum, idents = "B-cells") #endo_vas
# tum_sub_cells = WhichCells(object = tum, idents = "B-cells") #subset epi_lumhr cells #83891 cells
# write_rds(tum_sub, path = glue("{save_path}/b_subset_v1.rds"))

# save objects
percent.stress <- Matrix::colSums(GetAssayData(object = tum_sub)[stress_sig, ])/Matrix::colSums(GetAssayData(object = tum_sub))
atum_sub <- AddMetaData(tum_sub, metadata = percent.stress, col.name = "percent.stress")
write_rds(tum_sub, path = glue("{save_path}/t_subset_v1.rds"))
write.csv(tum_sub_cells, file = glue("{save_path}/t_cells_v1.csv"))

png(paste0(save_path, "vlnplot_percent.stress_sample_id.png"), width=14, height=5, units = "in", res = 300)
VlnPlot(tum_sub, features = "percent.stress", group.by = "final_sample_id", pt.size = 0 )
dev.off()




# method 9 ----------
# Log T integration with no scaling
library(future)
plan("multiprocess", workers = 20)
#plan("sequential")

# ** Clus function --------
# run clustering function
#seu_object = read_rds(glue("{save_path}/t_subset_v1.rds"))
seu_object = read_rds(glue("{save_path}/tum_sub_v1.rds"))
DefaultAssay(seu_object) <- "RNA"
dims = 35; k.param = 20; cols = colors_dark; var_split = "sheet_exp_proc"; var_scale = "nCount_RNA"
fts1 = c("nFeature_RNA","nCount_RNA", "percent.mito","percent.rb", "percent.stress")
fts2 = c("rna_EPCAM","rna_PTPRC","rna_PECAM1","rna_DCN","rna_KRT5",  "rna_CD86", "rna_COL1A1", "rna_PIP")
fts3 = c("rna_CD8A", "rna_CD8B", "rna_CD4", "rna_CCR7", "rna_IL7R", "rna_NKG7", "rna_GZMB", "rna_GZMK")
fts4 = c("rna_MS4A1","rna_CD40LG","rna_FOXP3","rna_IFNG","rna_LAG3","rna_TIGIT",  "rna_NCAM1", "rna_KLRB1", "rna_CD3D")

ref_set = NULL

clus_out = func_cluster_vst(tum_int = seu_object, var_scale = var_scale, var_split = var_split, dims = dims, k.param = k.param, 
                            markers1 = fts1, markers2 = fts2, markers3 = fts3, markers4 = fts4, save_path = save_path,
                            cols = colors_dark, ref_set = ref_set)

char_ident = "integrated_snn_res.0.4"


# ** Markers function (9)--------
dims = 30; k.param = 20;
char_ident = "integrated_snn_res.0.4"  #integrated_snn_res.0.1
ht_fts1 = c("CD3D", "CD8A", "CD8B")
#clus_out = read_rds(path = glue("{save_path}/dim_{dims}k_{k.param}_tum.integrated.ns_SCT_v4_intgeratedassay_{char_ident}_vst_method9.rds")) 
mark_out = func_markers_vst(clus_out = clus_out, char_ident = char_ident, fts1 = ht_fts1, dims = dims, k.param = k.param, cols = colors_dark, 
                            save_path = save_path)

clus_out = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v1/dim_30k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")
markers = read.csv("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v1/dims_30k_20_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.4_vst_method_9.csv")
clus_out = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v1/dim_30k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")
markers = read.csv("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v1/dims_30k_20_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.4_vst_method_9.csv")
tum = clus_out

markers_v2 = markers %>% dplyr::filter(gene %in% markers$gene)

top_markers = markers_v2 %>%
  dplyr::group_by(cluster) %>%
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
Idents(clus_out) <- "integrated_snn_res.0.4"

# filter again -----
tum_sub = subset(tum, idents = c("13", "14", "16"), invert  = TRUE)
write_rds(tum_sub, "tum_sub_v1.rds")


DimPlot(clus_out, group.by = "integrated_snn_res.0.4", label = T, split.by = "sheet_exp_proc", pt.size = 0.7, cols = colors_dark)
DimPlot(clus_out, group.by = "integrated_snn_res.0.6", label = T, pt.size = 0.3, cols = colors_dark)
VlnPlot(clus_out, group.by = "integrated_snn_res.0.4", 
        features = c("rna_MS4A1","rna_CD79A", "rna_PTPRC", "rna_CD3D", "rna_EPCAM", "rna_CD68", "rna_KRT8", 
                     "rna_COL1A1", "rna_RGS5", "rna_KRT5", "rna_NKG7", "rna_DCN", "rna_MKI67"), pt.size = 0, cols = colors_dark)

VlnPlot(clus_out, group.by = "integrated_snn_res.0.4", 
        features = c("rna_CD3D","rna_CD4", "rna_CD8A", "rna_CD8B", "rna_NKG7", "rna_GZMB", "rna_GZMA", 
                     "rna_PRF1", "rna_CCR7", "rna_FOXP3", "rna_ITGAX", "rna_KLRD1", "rna_MKI67"), pt.size = 0, cols = colors_dark)

# cd4 makrers ---
cd4_seu = FindAllMarkers(object = subset(clus_out, idents = c("0", "5", "6", "13")), assay = "RNA", logfc.threshold = 0)
write.csv(cd4_seu, "cd4_seu.csv")
markers_v2 = cd4_seu %>% dplyr::filter(gene %in% cd4_seu$gene)

top_markers = markers_v2 %>%
  dplyr::group_by(cluster) %>%
  top_n(15, avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()
cd4 = subset(clus_out, idents = c("0", "5", "6", "13"))
DefaultAssay(cd4) <- "integrated"
p1 = DoHeatmap(cd4, features = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B", top_markers), size = 2,raster = FALSE,
              group.by = "integrated_snn_res.0.4", slot = "scale.data", disp.min = -1.5, disp.max = 1.5, group.colors = colors_dark) + 
  scale_color_manual(values = colors_dark) +
  scale_fill_gradient2(low="steelblue1", high = "tomato3", mid = "white", midpoint = 0) + 
  theme(axis.text.y.left = element_text(size=6))

png("heatmap_CD4_Tcells_v1.png", width=10, height=6, units = "in", res = 300)
p1
dev.off()

# ch8 markers ---
clus_out1 = clus_out
clus_out1$integrated_snn_res.0.4 = as.character(clus_out1$integrated_snn_res.0.4)
clus_out1$integrated_snn_res.0.4[clus_out1$integrated_snn_res.0.8 %in% c("6")] = "CD8_pip"
clus_out$integrated_snn_res.0.4 = factor(clus_out$integrated_snn_res.0.4, 
                                              levels = c(0, seq(1:16), "CD8_pip"))

Idents(clus_out) <- "integrated_snn_res.0.4"

cd8_seu = FindAllMarkers(object = subset(clus_out, idents = c("CD8_pip", "1", "2", "9", "10", "11", "15")), assay = "RNA", logfc.threshold = 0)
write.csv(cd8_seu, "cd8_seu.csv")

markers_v2 = cd8_seu %>% dplyr::filter(gene %in% cd8_seu$gene)

top_markers = markers_v2 %>%
  dplyr::group_by(cluster) %>%
  top_n(10, avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()

cd8 = subset(clus_out, idents = c("CD8_pip","1", "2", "9", "10", "11", "15"))
DefaultAssay(cd8) <- "integrated"
p1= DoHeatmap(cd8, features = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B", top_markers), size = 2,raster = FALSE,
              group.by = "integrated_snn_res.0.4", slot = "scale.data", disp.min = -1.5, disp.max = 1.5, group.colors = colors_dark) + 
  scale_color_manual(values = colors_dark) +
  scale_fill_gradient2(low="steelblue1", high = "tomato3", mid = "white", midpoint = 0) + 
  theme(axis.text.y.left = element_text(size=6))

png("heatmap_CD8_Tcells_v1.png", width=12, height=7.5, units = "in", res = 300)
p1
dev.off()

reduction = "umap"
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
p6 = VlnPlot(tum, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), group.by = "integrated_snn_res.0.4", pt.size = 0, cols = colors_dark)

VlnPlot(tum, features = c("rna_KRT14", "rna_GZMA", "rna_EOMES"), pt.size = 0, group.by = "integrated_snn_res.0.4")
VlnPlot(tum, features = c("rna_GZMB", "rna_GNLY", "rna_NKG7", "rna_IFNG", "rna_PRF1", "rna_CD4", "rna_STAT4", "rna_TBX21"), pt.size = 0)
VlnPlot(tum, features = c("rna_MX1", "rna_ISG15", "rna_IFIT3", "PDCD1"), pt.size = 0)
VlnPlot(tum, features = c("rna_CD4", "rna_TRAC", "rna_TRBC1", "rna_CD3D", "rna_CD3G", "rna_CD3E", "rna_KLRG1"), pt.size = 0)
FeaturePlot(tum, features = c("rna_KRT14"))


png(paste0(save_path,"/", reduction,dims,"k_", k.param,  "umap_all_labelled_Tcells.png"), width=16, height=16, units = "in", res = 300)
print(p1.1)
dev.off()


png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_all_markers.png"), width=12, height=9, units = "in", res = 300)
print(p2)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_all_markers2.png"), width=12, height=14, units = "in", res = 300)
print(p3)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param,  "umap_all_markers3.png"), width=12, height=12, units = "in", res = 300)
print(p4)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "tissue_resident.png"), width=16, height=12, units = "in", res = 300)
print(p5)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "QC.png"), width=16, height=12, units = "in", res = 300)
print(p6)
dev.off()

# assign clusters ------
clus_out$integrated_snn_res.0.4[is.na(clus_out$integrated_snn_res.0.4)] = "0"


clus_out$cellstate = "NULL"
clus_out$cellstate[clus_out$integrated_snn_res.0.4 %in% c("0")] = "CD4_il7r_ccl20"
clus_out$cellstate[clus_out$integrated_snn_res.0.4 %in% c("5")] = "CD4_naive"
clus_out$cellstate[clus_out$integrated_snn_res.0.4 %in% c("6")] = "CD4_gzma_k_ccl5"
clus_out$cellstate[clus_out$integrated_snn_res.0.4 %in% c("13")] = "CD4_treg"
clus_out$cellstate[clus_out$integrated_snn_res.0.4 %in% c("1")] = "CD8_low_quality"
clus_out$cellstate[clus_out$integrated_snn_res.0.4 %in% c("2")] = "CD8_exhaus_layn"
clus_out$cellstate[clus_out$integrated_snn_res.0.4 %in% c("9")] = "CD8_plgc2_znf683"
clus_out$cellstate[clus_out$integrated_snn_res.0.4 %in% c("10")] = "CD8_SAA1"
clus_out$cellstate[clus_out$integrated_snn_res.0.4 %in% c("11")] = "CD8_KLRC1"
clus_out$cellstate[clus_out$integrated_snn_res.0.4 %in% c("15")] = "CD8_KLF2"
clus_out$cellstate[clus_out$integrated_snn_res.0.4 %in% c("16")] = "prol"
clus_out$cellstate[clus_out$integrated_snn_res.0.4 %in% c("CD8_pip")] = "CD8_pip"
clus_out$cellstate[clus_out$integrated_snn_res.0.4 %in% c("3")] = "NKT_gzmk"
clus_out$cellstate[clus_out$integrated_snn_res.0.4 %in% c("8")] = "NKT_gzmb"
clus_out$cellstate[clus_out$integrated_snn_res.0.4 %in% c("4")] = "NK_FCERG_CD56low"
clus_out$cellstate[clus_out$integrated_snn_res.0.4 %in% c("12")] = "NK_XCL1_CD56high"
clus_out$cellstate[clus_out$integrated_snn_res.0.4 %in% c("7")] = "HSPs"
clus_out$cellstate[clus_out$integrated_snn_res.0.4 %in% c("14")] = "MMP_T"

lvl = names(table(clus_out$cellstate))
clus_out$cellstate = factor(clus_out$cellstate, levels = lvl)

write_rds(clus_out, "clus_out.rds")
clus_out = read_rds("clus_out.rds")
Idents(clus_out) <- "cellstate"

clus_out_seu = FindAllMarkers(object = clus_out, assay = "RNA", logfc.threshold = 0)
write.csv(clus_out_seu, "clus_out_seu.csv")

markers_v2 = clus_out_seu %>% dplyr::filter(gene %in% clus_out_seu$gene)

top_markers = markers_v2 %>%
  dplyr::group_by(cluster) %>%
  top_n(7, avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()

DefaultAssay(clus_out) <- "integrated"
p1= DoHeatmap(clus_out, features = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B", top_markers), size = 2,raster = FALSE, label=TRUE,
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
col_t = c("#C70E7B", "#FC6882", "#007BC3", "#54BCD1", "#EF7C12", "#F4B95A", "#009F3F", "#8FDA04", "#AF6125", "#F4E3C7", "#B25D91", "#EFC7E6",
          "#EF7C12", "#F4B95A", colors_dark)
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



