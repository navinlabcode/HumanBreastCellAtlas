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
tum$sheet_patient_id[tum$sheet_ts_sampleid %in% c("hbca25_c")] = "pt37"

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
tum_sub <- AddMetaData(tum_sub, metadata = percent.stress, col.name = "percent.stress")
write_rds(tum_sub, path = glue("{save_path}/t_subset_v1.rds"))
write.csv(tum_sub_cells, file = glue("{save_path}/t_cells_v1.csv"))

png(paste0(save_path, "vlnplot_percent.stress_sample_id.png"), width=14, height=5, units = "in", res = 300)
VlnPlot(tum_sub, features = "percent.stress", group.by = "final_sample_id", pt.size = 0 )
dev.off()




# method 9 ----------
# Log T integration with no scaling
library(future)
plan("multiprocess", workers = 70)
#plan("sequential")

# ** Clus function --------
# run clustering function
seu_object = read_rds(glue("{save_path}/t_subset_v1.rds"))
#seu_object = read_rds(glue("{save_path}/tum_sub_v1.rds"))
DefaultAssay(seu_object) <- "RNA"
dims = 20; k.param = 20; cols = colors_dark; var_split = "sheet_exp_proc"; var_scale = "nCount_RNA"
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
dims = 20; k.param = 20;
char_ident = "integrated_snn_res.0.4"  #integrated_snn_res.0.1
ht_fts1 = c("CD3D", "CD8A", "CD8B")
#clus_out = read_rds(path = glue("{save_path}/dim_{dims}k_{k.param}_tum.integrated.ns_SCT_v4_intgeratedassay_{char_ident}_vst_method9.rds")) 
mark_out = func_markers_vst(clus_out = clus_out, char_ident = char_ident, fts1 = ht_fts1, dims = dims, k.param = k.param, cols = colors_dark, 
                            save_path = save_path)

clus_out = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v1/dim_35k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")
markers = read.csv("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v1/dims_35k_20_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.4_vst_method_9.csv")
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
Idents(tum) <- "integrated_snn_res.0.4"
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

tavg = AverageExpression(tum, features = unique(c("CD3D", "CD3G", "CD8A", "CD8B", top_markers)), return.seurat = T)
DoHeatmap(tavg, features = unique(c("CD3D", "CD3G", "CD8A", "CD8B", top_markers)), assay = "integrated", label = TRUE,
          slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = colors_dark, raster = F)  +
  scale_fill_gradient2(low = "#29A7B8", mid = "white", high = "#E85041", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5)) -> tcell_heatmap
png("heatmap_all_labelled_Tcells_0.4_v1.png", width=12, height=18, units = "in", res = 300)
tcell_heatmap
dev.off()

# filter other cells --------
Idents(clus_out) <- "integrated_snn_res.0.4"

# filter again -----
tum_sub = subset(tum, idents = c("12", "14", "17", "18", "2", "15"), invert  = TRUE)
write_rds(tum_sub, "tum_sub_v1.rds")


