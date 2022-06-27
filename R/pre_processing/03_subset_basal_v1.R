# geo ----
setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_basal_v1";dir.create(odir)
options(future.globals.maxSize = 90000 * 1024^2)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "figure_1_v1"
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_basal_v1/"

# core ----
setwd("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_basal_v1"; dir.create(odir)
library(future)
plan("multiprocess", workers = 50)
options(future.globals.maxSize = Inf)
setwd(paste0("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "figure_1_v1"
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_basal_v1/"
# stress_sig = read.delim("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_fibroblasts_final1/stromal_stress_signature_032720.txt")
# stress_sig = as.character(stress_sig$stromal.stress.sig)


# read original object --------------
tum = read_rds(paste0(wd, data_folder,"/final_use.rds"))
tum_sub = tum

# save objects
Idents(tum_sub) <- "final_group1"
tum_sub = subset(tum_sub, subset = final_group1 %in% c("Basal"))
write_rds(tum_sub, path = glue("{save_path}/basal_subset_v1.rds"))

png(paste0(save_path, "vlnplot_nCount_RNA_sample_id.png"), width=24, height=5, units = "in", res = 300)
VlnPlot(tum_sub, features = "nCount_RNA", group.by = "final_sample_id", pt.size = 0 )
dev.off()


# method 9 ----------
# Log T integration with no scaling
library(future)
plan("multiprocess", workers = 20)
#plan("sequential")

# ** Clus function --------
# run clustering function
#seu_object = read_rds(glue("{save_path}/t_subset_v1.rds"))
seu_object = read_rds(glue("{save_path}/basal_subset_v1.rds"))
DefaultAssay(seu_object) <- "RNA"
dims = 35; k.param = 20; cols = colors_dark; var_split = "sheet_exp_proc"; var_scale = "nCount_RNA"
fts1 = c("nFeature_RNA","nCount_RNA", "percent.mito","percent.rb", "percent.stress")
fts2 = c("rna_EPCAM","rna_PTPRC","rna_PECAM1","rna_DCN","rna_KRT5",  "rna_CD86", "rna_COL1A1", "rna_PIP")
fts3 = c("rna_KRT9", "rna_KRT8", "rna_KRT17", "rna_SLPI", "rna_ANKRD30A", "rna_TFF1", "rna_AGR3", "rna_EFHD1")
fts4 = c("rna_ESR1","rna_PGR","rna_AR","rna_BRCA1","rna_KRT5","rna_KRT18",  "rna_PCNA", "rna_MKI67", "rna_LTF")


ref_set = NULL

clus_out = func_cluster_vst(tum_int = seu_object, var_scale = var_scale, var_split = var_split, dims = dims, k.param = k.param, 
                            markers1 = fts1, markers2 = fts2, markers3 = fts3, markers4 = fts4, save_path = save_path,
                            cols = colors_dark, ref_set = ref_set)

char_ident = "integrated_snn_res.0.4"

p1 = DimPlot(clus_out, split.by = "sheet_patient_id", group.by = "integrated_snn_res.0.4", cols = colors_darj, label = F, shuffle = T, pt.size = 0.6, ncol = 10)
png("vln_allpts_v1.png", width=24, height=20, units = "in", res = 300)
p1
dev.off()


# ** Markers function (9)--------
dims = 35; k.param = 20;
char_ident = "integrated_snn_res.0.2"  #integrated_snn_res.0.1
ht_fts1 = c("EPCAM", "KRT17", "AR")
#clus_out = read_rds(path = glue("{save_path}/dim_{dims}k_{k.param}_tum.integrated.ns_SCT_v4_intgeratedassay_{char_ident}_vst_method9.rds")) 
mark_out = func_markers_vst(clus_out = clus_out, char_ident = char_ident, fts1 = ht_fts1, dims = dims, k.param = k.param, cols = colors_dark, 
                            save_path = save_path)

clus_out = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_basal_v1/dim_35k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")
markers = read.csv("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_basal_v1/dims_35k_20_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.2_vst_method_9.csv")
#clus_out = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_basal_v1/dim_35k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")
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
p1 = DoHeatmap(tum, features = c("EPCAM", "TAGLN", "MYLK", "ACTA2", "KRT14", top_markers), size = 2,raster = FALSE,
               group.by = char_ident, slot = "scale.data", disp.min = -1.5, disp.max = 1.5, group.colors = colors_darj) + 
  scale_color_manual(values = colors_darj) +
  scale_fill_gradient2(low="steelblue1", high = "tomato3", mid = "white", midpoint = 0) + 
  theme(axis.text.y.left = element_text(size=6))

png("heatmap_all_labelled_res4_cells_v1.png", width = 18, height = 12, units = "in", res = 200)
p1
dev.off()

tavg = AverageExpression(tum, features = unique(c("EPCAM", "TAGLN", "MYLK", "ACTA2", "KRT14", top_markers)), return.seurat = T)
DoHeatmap(tavg, features = unique(c("EPCAM", "TAGLN", "MYLK", "ACTA2", "KRT14", top_markers)), assay = "integrated", label = TRUE,
          slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = colors_dark, raster = F)  +
  scale_fill_gradient2(low = "#29A7B8", mid = "white", high = "#E85041", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=11)) -> tcell_heatmap
png("heatmap_all_labelled_basal_cellstatetk2_v1.png", width=12, height=22, units = "in", res = 300)
tcell_heatmap
dev.off()


v1 = VlnPlot(tum, features = c("rna_ESR1", "rna_PGR", "rna_AR", "rna_BRCA1", "rna_SLPI", "rna_KRT17", "rna_PIP", "rna_HES1", "rna_APOD", "rna_ACTB", 
                               "rna_EPCAM","rna_KRT9", "rna_KRT8", "rna_KRT17", "rna_SLPI", "rna_ANKRD30A", "rna_TFF1", "rna_AGR3", "rna_EFHD1"), pt.size = 0, cols = colors_darj)

f1 = FeaturePlot(tum, features = c("rna_ESR1", "rna_PGR", "rna_AR", "rna_BRCA1", "rna_SLPI", "rna_KRT17", "rna_PIP", "rna_HES1", "rna_APOD", "rna_ACTB", 
                                   "rna_EPCAM","rna_KRT9", "rna_KRT8", "rna_KRT17", "rna_SLPI", "rna_ANKRD30A", "rna_TFF1", "rna_AGR3", "rna_EFHD1"), pt.size = 0.001)

v1=VlnPlot(clus_out, features = c("KRT14","KRT17","KRT5","ACTG2","TUBB2B","COL17A1","LAMA3","SERPINB5"), group.by = "integrated_snn_res.0.2", pt.size = 0)
v2=VlnPlot(clus_out, features = c("AREG","ANKRD30A","AGR2","AGR3","TMC5","DNAJC12","RASEF","MLPH","CA12","C3orf52","ACADSB","MAPK8"), group.by = "integrated_snn_res.0.2", pt.size = 0)
v3=VlnPlot(clus_out, features = c("SLPI","LTF","KRT15","MMP7","CCL28","PIGR","RARRES1","ALDH1A3","LCN2","SNORC","S100A1","GABRP"), group.by = "integrated_snn_res.0.2", pt.size = 0, assay = "RNA")


png("VlnPlot_all_labelled_res2_cells_v1.png", width=11, height=6, units = "in", res = 200)
v1
dev.off()
png("VlnPlot_all_labelled_res2_cells_v2.png", width=11, height=6, units = "in", res = 200)
v2
dev.off()
png("VlnPlot_all_labelled_res2_cells_v3.png", width=11, height=6, units = "in", res = 200)
v3
dev.off()
png("VlnPlot_all_labelled_res2_qced.png", width=11, height=3, units = "in", res = 200)
VlnPlot(tum, c("nCount_RNA", "nFeature_RNA"), pt.size = 0, cols = colors_dark)
dev.off()

VlnPlot(tum, c("rna_KRT17", "rna_KRT5","rna_KRT14","rna_KRT19", "rna_DCN"), pt.size = 0, cols = colors_dark)

png("VlnPlot_all_labelled_res4_cells_v1.png", width=11, height=11, units = "in", res = 200)
v1
dev.off()

png("FeaturePlot_all_labelled_res4_cells_v1.png", width=11, height=11, units = "in", res = 200)
f1
dev.off()

# manual checks ------
#tum = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_lymphendo_v1/dim_12k_20_nCount_RNA_tum.integrated.ns2_method_9.rds")

# filter other cells --------
Idents(clus_out) <- "integrated_snn_res.0.2"

# filter again -----
tum_sub = subset(tum, idents = c("13", "3", "11", "10", "12", "5"), invert  = TRUE)
write_rds(tum_sub, "tum_sub_v1.rds")


DimPlot(clus_out, group.by = "integrated_snn_res.0.4", label = T, split.by = "sheet_exp_proc", pt.size = 0.7, cols = colors_dark)
DimPlot(clus_out, group.by = "integrated_snn_res.0.6", label = T, pt.size = 0.3, cols = colors_dark)
VlnPlot(clus_out, group.by = "integrated_snn_res.0.4", 
        features = c("rna_MS4A1","rna_CD79A", "rna_PTPRC", "rna_CD3D", "rna_EPCAM", "rna_CD68", "rna_KRT8", 
                     "rna_COL1A1", "rna_RGS5", "rna_KRT5", "rna_NKG7", "rna_DCN", "rna_MKI67"), pt.size = 0, cols = colors_dark)

VlnPlot(clus_out, group.by = "integrated_snn_res.0.4", 
        features = c("rna_CD3D","rna_CD4", "rna_CD8A", "rna_CD8B", "rna_NKG7", "rna_GZMB", "rna_GZMA", 
                     "rna_PRF1", "rna_CCR7", "rna_FOXP3", "rna_ITGAX", "rna_KLRD1", "rna_MKI67"), pt.size = 0, cols = colors_dark)

