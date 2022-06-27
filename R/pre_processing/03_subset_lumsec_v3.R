# geo ----
setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_lumsec_v3";dir.create(odir)
options(future.globals.maxSize = 90000 * 1024^2)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "03_subset_lumsec_v2"
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_v3/"

# core ----
setwd("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_lumsec_v3"; dir.create(odir)
library(future)
plan("multiprocess", workers = 75)
options(future.globals.maxSize = Inf)
setwd(paste0("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "03_subset_lumsec_v2"
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_v3/"
# stress_sig = read.delim("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_fibroblasts_final1/stromal_stress_signature_032720.txt")
# stress_sig = as.character(stress_sig$stromal.stress.sig)


# read original object --------------
tum = read_rds(paste0(wd, data_folder,"/clus_out.rds"))

png(paste0(save_path, "vlnplot_nCount_RNA_sample_id.png"), width=4, height=3, units = "in", res = 300)
VlnPlot(tum, features = "nCount_RNA", group.by = "integrated_snn_res.0.4", pt.size = 0 )
dev.off()

# subset and filter ----------
Idents(tum) <- "cellstate_tk2"
tum_sub = subset(tum, idents = c("ls_SPINK5"), invert = T) #endo_vas
write_rds(tum_sub, "epi_lumsec_sub_v1.rds")

# method 9 ----------
# Log T integration with no scaling
library(future)
plan("multiprocess", workers = 85)
#plan("sequential")

# ** Clus function --------
# run clustering function
#seu_object = read_rds(glue("{save_path}/t_subset_v1.rds"))
seu_object = read_rds(glue("{save_path}/epi_lumsec_sub_v1.rds"))
DefaultAssay(seu_object) <- "RNA"
dims = 35; k.param = 20; cols = colors_dark; var_split = "sheet_exp_proc"; var_scale = "nCount_RNA"
fts1 = c("nFeature_RNA","nCount_RNA", "percent.mito","percent.rb", "percent.stress")
fts2 = c("rna_EPCAM","rna_PTPRC","rna_PECAM1","rna_DCN","rna_KRT5",  "rna_CD86", "rna_COL1A1", "rna_PIP")
fts3 = c("rna_S100A8", "rna_SERPINB4", "rna_FBLN5", "rna_CCL2", "rna_LTF","rna_SAA1", "rna_TNFAIP2", "rna_KYNU", "rna_NDRG1")
fts4 = c("rna_KRT15", "rna_KRT14", "rna_SCGB1D2", "rna_PTN", "rna_CITED4","rna_MYC", "rna_LEF3", "rna_CALML5", "rna_RGS2",  "rna_PCNA", "rna_MKI67", "rna_LTF")


ref_set = NULL

clus_out = func_cluster_vst(tum_int = seu_object, var_scale = var_scale, var_split = var_split, dims = dims, k.param = k.param, 
                            markers1 = fts1, markers2 = fts2, markers3 = fts3, markers4 = fts4, save_path = save_path,
                            cols = colors_dark, ref_set = ref_set)

char_ident = "integrated_snn_res.0.2"

p1 = DimPlot(clus_out, split.by = "sheet_patient_id", group.by = "cellstate_tk1", cols = colors_darj, label = F, shuffle = T, pt.size = 0.6, ncol = 10)
png("vln_allpts_v1.png", width=24, height=20, units = "in", res = 300)
p1
dev.off()

p1 = DimPlot(clus_out, split.by = "sheet_patient_id", group.by = "integrated_snn_res.0.2", 
             cols = colors_dark, label = F, pt.size = 0.6, ncol = 10)
png("dim_allpts_0.2v1.png", width=24, height=20, units = "in", res = 300)
p1
dev.off()

p2 = DimPlot(clus_out, group.by = "integrated_snn_res.0.2", cols = colors_dark, label = T)
png("umap0.2_v1.png", width=4, height=3, units = "in", res = 300)
p2
dev.off()
FeaturePlot(clus_out, features = c("rna_LALBA"))
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

clus_out = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_v3/dim_35k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")
markers = read.csv("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_v3/dims_35k_20_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.2_vst_method_9.csv")
#clus_out = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_v3/dim_35k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")
# markers = read.csv("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v1/dims_30k_20_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.4_vst_method_9.csv")
tum = clus_out

markers_v2 = markers %>% dplyr::filter(gene %in% markers$gene)

top_markers = markers_v2 %>%
  dplyr::group_by(cluster) %>%
  slice_max(avg_logFC, n=12) %>%
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
  theme(axis.text.y.left = element_text(size=7))

png("heatmap_all_labelled_res2_cells_v1.png", width=12, height=8, units = "in", res = 200)
p1
dev.off()

tavg = AverageExpression(tum, features = unique(c("EPCAM", "KIT", "TOP2A", top_markers)), return.seurat = T)
tavg = ScaleData(tavg, assay = "RNA")
DoHeatmap(tavg, features = unique(c("EPCAM", "KIT", "TOP2A", top_markers)), assay = "RNA", label = TRUE,
          slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = colors_dark, raster = F)  +
  scale_fill_gradient2(low = "#29A7B8", mid = "white", high = "#E85041", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5)) -> tcell_heatmap
png("heatmap_all_labelled_lumsec_0.2_v1.png", width=6, height=10, units = "in", res = 300)
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

write_rds(clus_out, "lumsec_clus_out.rds")
clus_out = read_rds("clus_out.rds")
Idents(clus_out) <- "cellstate_tk2"

clus_out_seu = FindAllMarkers(object = clus_out, assay = "RNA", logfc.threshold = 0)
write.csv(clus_out_seu, "clus_out_seu.csv")

DimPlot(clus_out, group.by = "cellstate_tk1", label = T,cols = colors_darj)

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


# check kevins clusters ------
kevin_obj = load("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_v10/hbca_seclum.RObj")
kevin_obj = hbca.seclum

kevin_obj_df = data.frame(cellname = kevin_obj$cellname, kevin_states = kevin_obj$seclum_states)
clus_out_meta = clus_out@meta.data
clus_out_meta = left_join(clus_out_meta, kevin_obj_df)

clus_out@meta.data = clus_out_meta
clus_out$kevin_states[(is.na(clus_out$kevin_states))] = "unkn"
rownames(clus_out@meta.data) = clus_out$cellname

png("umap_kevin_Tcells_v2.png", width=9, height=7, units = "in", res = 300) #
DimPlot(clus_out, group.by = "kevin_states", label = T, cols = colors_dark, pt.size = 0.0001)
dev.off()

png("umap_kevin_split_by_pt_v2.png", width=18, height=18, units = "in", res = 300) #
DimPlot(clus_out, group.by = "kevin_states", split.by = "sheet_patient_id", ncol = 6, label = T, cols = colors_dark, pt.size = 0.0001)
dev.off()

png("umap_kevin_split_by_pt_v2.png", width=18, height=24, units = "in", res = 300) #
DimPlot(kevin_obj, group.by = "seclum_states", split.by = "sheet_patient_id", ncol = 6, cols = colors_dark, pt.size = 1)
dev.off()

# confusion matrix ------
cc = table(clus_out$cellstate_tk2, clus_out$kevin_states)
g1 = ggplot(as.data.frame(cc), aes(Var1,Var2, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  #scale_fill_viridis_c()+
  scale_fill_gradient2(low="white",high="#009194") +
  labs(x = "Tapsi's labels",y = "Kevin's labels")
png(paste0(save_path,"/tapsi_kevin_confusionmat_v2.png"), width=8, height=6, units = "in", res = 300)
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


