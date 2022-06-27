# geo ----
setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_pericytes_v2";dir.create(odir)
options(future.globals.maxSize = 50000 * 1024^2)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "03_subset_pericytes_v1"
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_pericytes_v2/"
# stress_sig = read.delim("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_fibroblasts_final1/stromal_stress_signature_032720.txt")
# stress_sig = as.character(stress_sig$stromal.stress.sig)

# core ----
setwd("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_pericytes_v2"; dir.create(odir)
library(future)
plan("multiprocess", workers = 100)
options(future.globals.maxSize = Inf)
setwd(paste0("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "03_subset_pericytes_v1"
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_pericytes_v2/"
# stress_sig = read.delim("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_fibroblasts_final1/stromal_stress_signature_032720.txt")
# stress_sig = as.character(stress_sig$stromal.stress.sig)


# read original object --------------
#sftp://10.132.80.157/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_pericytes_v1/pericytes_final_use.rds
tum = readRDS(file = paste0(wd, data_folder,"/pericytes_final_use.rds"))

# subset ---------
Idents(tum) <- "integrated_snn_res.0.085"
tum_sub = subset(tum, idents = "3", invert = T) #subset object

# save objects
# percent.stress <- Matrix::colSums(GetAssayData(object = tum_sub)[stress_sig, ])/Matrix::colSums(GetAssayData(object = tum_sub))
# tum_sub <- AddMetaData(tum_sub, metadata = percent.stress, col.name = "percent.stress")
write_rds(tum_sub, path = glue("{save_path}/pericyte_subset_v1.rds"))

png(paste0(save_path, "vlnplot_percent.stress_sample_id.png"), width=24, height=5, units = "in", res = 300)
VlnPlot(tum_sub, features = "nCount_RNA", group.by = "final_sample_id", pt.size = 0 )
dev.off()

# method 9 ----------
# Log T integration with no scaling
library(future)
plan("multiprocess", workers = 45)
#plan("sequential")


# ** ONly clustering ---------
seu_object = read_rds(glue("{save_path}/pericyte_subset_v1.rds"))
DefaultAssay(seu_object) <- "RNA"
dims = 22; k.param = 20;reduction = "umap"
seu_object = RunUMAP(seu_object, reduction = "pca", dims = 1:dims)
seu_object = FindNeighbors(object = seu_object, dims = 1:dims, reduction = "pca", k.param = k.param)
seu_object = FindClusters(object = seu_object, resolution = c(1, 0.8, 0.6,0.4,0.3,0.2, 0.15,0.17, 0.1, 0.085, 0.090, 0.075, 0.05, 0.023), verbose = TRUE)
p1 = DimPlot(object = seu_object, reduction = reduction, label = TRUE, pt.size = 0.001, cols = cols, group.by = "integrated_snn_res.0.1")+  NoAxes()
p2 = DimPlot(object = seu_object, reduction = reduction, label = FALSE, pt.size = 0.001, cols = cols, group.by = "integrated_snn_res.0.2") +  NoAxes()
p3 = DimPlot(object = seu_object, reduction = reduction, group.by = "procedure",  pt.size = 0.00001, cols = c("tan1","steelblue1")) + NoAxes()
p4 = DimPlot(object = seu_object, reduction = reduction, group.by = "batch",  pt.size = 0.0001, cols = c("darkorange","darkorchid")) +  NoAxes()
p5 = DimPlot(object = seu_object, reduction = reduction, group.by = "sheet_patient_id", pt.size = 0.00013, cols = colors_dark) +  NoAxes()
p6 = DimPlot(object = seu_object, reduction = reduction, group.by = "exp_proc", pt.size = 0.00003, cols = colors_dark) +  NoAxes()


# ** Clus function --------
# run clustering function
seu_object = read_rds(glue("{save_path}/pericyte_subset_v1.rds"))
DefaultAssay(seu_object) <- "RNA"

dims = 15; k.param = 20; cols = colors_dark; var_split = "sheet_exp_proc"; var_scale = "nCount_RNA"
fts1 = c("nFeature_RNA","nCount_RNA", "percent.mito","percent.rb", "percent.stress")
fts2 = c("rna_EPCAM","rna_PTPRC","rna_PECAM1","rna_DCN","rna_KRT5",  "rna_LTF", "rna_KRT17", "rna_HPGD")
fts3 = c("rna_COL18A1", "rna_COL1A1", "rna_IFI27", "rna_SPARCL1", "rna_ZEB2", "rna_CPE", "rna_NR4A1", "rna_C11orf96")
fts4 = c("rna_MCAM","rna_ADAMTS1","rna_ADAMTS9","rna_NR4A1","rna_ZFHX3","rna_PDK4",  "rna_EPAS1", "rna_TAGLN", "rna_PDGFRB")

#ref_set = c("pt01", "pt06","pt09", "pt20", "pt21","pt25", "pt27", "pt28", "pt31")
ref_set = NULL

seu_object$regroup = seu_object$sheet_patient_id
seu_object$regroup[seu_object$sheet_patient_id == "pt08"] = "pt_cmb"
seu_object$regroup[seu_object$sheet_patient_id == "pt22"] = "pt_cmb"
seu_object$regroup[seu_object$sheet_patient_id == "pt27"] = "pt_cmb"
seu_object$regroup[seu_object$sheet_patient_id == "pt05"] = "pt_cmb"
seu_object$regroup[seu_object$sheet_patient_id == "pt38"] = "pt_cmb"
seu_object$regroup[seu_object$sheet_patient_id == "pt48"] = "pt_cmb"
seu_object$regroup[seu_object$sheet_patient_id == "pt40"] = "pt_cmb"
seu_object$regroup[seu_object$sheet_patient_id == "pt15"] = "pt_cmb"
seu_object$regroup[seu_object$sheet_patient_id == "pt19"] = "pt_cmb"
seu_object$regroup[seu_object$sheet_patient_id == "pt55"] = "pt_cmb"
seu_object$regroup[seu_object$sheet_patient_id == "pt05"] = "pt_cmb"


clus_out = func_cluster_vst(tum_int = seu_object, var_scale = var_scale, var_split = var_split, dims = dims, k.param = k.param, 
                            markers1 = fts1, markers2 = fts2, markers3 = fts3, markers4 = fts4, save_path = save_path,
                            cols = colors_dark, ref_set = ref_set)

clus_out = func_cluster_scale(tum_int = seu_object, var_scale = var_scale, var_split = var_split, dims = dims, k.param = k.param, 
                              markers1 = fts1, markers2 = fts2, markers3 = fts3, markers4 = fts4, save_path = save_path,
                              cols = colors_dark)

char_ident = "integrated_snn_res.0.4"

# method 10 ----------
# SCTransform, no integration
library(future)
plan("multiprocess", workers = 20)
#plan("sequential")

# ** Clus function --------
# run clustering function
seu_object = read_rds(glue("{save_path}/epi_basal_subset_v1.rds"))

dims = 18; k.param = 31; cols = colors_dark; var_scale = "nCount_RNA"
fts1 = c("nFeature_RNA","nCount_RNA", "percent.mito","percent.rb", "percent.stress")
fts2 = c("rna_EPCAM","rna_PTPRC","rna_PECAM1","rna_DCN","rna_KRT5",  "rna_LTF", "rna_KRT17", "rna_HPGD")
fts3 = c("rna_KRT9", "rna_KRT8", "rna_KRT17", "rna_SLPI", "rna_ANKRD30A", "rna_TFF1", "rna_AGR3", "rna_EFHD1")
fts4 = c("rna_ESR1","rna_PGR","rna_AR","rna_BRCA1","rna_KRT5","rna_KRT18",  "rna_EPCAM", "rna_SLPI", "rna_LTF")

clus_out = func_cluster_sct(tum_int = seu_object, var_scale = var_scale, dims = dims, k.param = k.param, 
                            markers1 = fts1, markers2 = fts2, 
                            markers3 = fts3, markers4 = fts4,cols = colors_dark,
                            save_path = save_path)


# CHECK CLUSTERING AND SAVE
write_rds(tum.integrated.ns, glue("{save_path}/dim_{dims}k_{k.param}_split_by_{var_split}_scaled_by_{var_scale}_tum.integrated.ns3_method_13.rds")) 

# read from geo -------
a = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/04_subset_epi_lumhrblasts_v1/dim_10k_30_nCount_RNA_tum.integrated.ns2_method_9.rds")

# read from core ---------
a = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/04_subset_epi_lumhrblasts_v1/dim_10k_30_nCount_RNA_tum.integrated.ns2_method_9.rds")

# ** Markers function (9)--------
dims = 15; k.param = 20;
char_ident = "integrated_snn_res.0.085"  #integrated_snn_res.0.1
ht_fts1 = c("RGS5")
#clus_out = read_rds(path = glue("{save_path}/dim_{dims}k_{k.param}_tum.integrated.ns_SCT_v4_intgeratedassay_{char_ident}_vst_method9.rds")) 
mark_out = func_markers_vst(clus_out = clus_out, char_ident = char_ident, fts1 = ht_fts1, dims = dims, k.param = k.param, cols = colors_dark, 
                            save_path = save_path)

# ** Markers function (13)--------
dims = 10; k.param = 30; cols = colors_dark; var_split = "patient_id"; var_scale = "percent.stress"
char_ident = "integrated_snn_res.0.1"  #integrated_snn_res.0.1
ht_fts1 = c("DCN", "COL1A1", "VIM")
mark_out = func_markers_scale(clus_out = d, var_split = var_split, var_scale = var_scale, char_ident = char_ident, dims = dims, 
                              fts1 = ht_fts1, k.param = k.param, cols = colors_dark, 
                              save_path = save_path)

# ** Markers function (10)--------
dims = 20; k.param = 31;  cols = colors_dark; var_scale = "percent.stress"
char_ident = "SCT_snn_res.0.1"  #integrated_snn_res.0.1
ht_fts1 = c("EPCAM")
mark_out = func_markers_sct(clus_out = tum.integrated.ns, char_ident = char_ident, var_scale = var_scale, dims = dims, fts1 = ht_fts1, k.param = k.param, cols = colors_dark, 
                            save_path = save_path)

# manual checks ------
tum = read_rds("~/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_pericytes_v1/dim_15k_20_nCount_RNA_tum.integrated.ns2_method_9.rds")
tum = FindClusters(object = tum.integrated.ns, resolution = c(1, 0.8, 0.6,0.4,0.3,0.2, 0.15,0.17, 0.1, 0.085, 0.090, 0.075, 0.05, 0.023), verbose = TRUE)

FeaturePlot(tum, features = c("rna_RGS5", "rna_TAGLN", "rna_ACTA2", "rna_MYLK", "rna_CSPG4", "rna_PDGFRB"))
FeaturePlot(tum, features = c("rna_DES", "rna_MCAM", "rna_EPHX3", "rna_KCNJ8", "rna_IFITM3", "rna_HIGD1B"))

DimPlot(tum, group.by = "integrated_snn_res.0.085")

# system("pdftoppm input.pdf outputname -png")
# save plots
png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale, "umap_all_labelled_method_9a.png"), width=6, height=5, units = "in", res = 300)
print(p1)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale, "umap_all_method_9a.png"), width=6, height=5, units = "in", res = 300)
print(p2)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_procedure_method_9a.png"), width=6, height=5, units = "in", res = 300)
print(p3)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_source_method_9a.png"), width=6, height=5, units = "in", res = 300)
print(p4)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_patient_method_9a.png"), width=6, height=5, units = "in", res = 300)
print(p5)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_sample_id_method_9a.png"), width=6, height=5, units = "in", res = 300)
print(p6)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_exp_proc_method_9a.png"), width=6, height=5, units = "in", res = 300)
print(p7)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_markers1_method_9.png"), width=16, height=15, units = "in", res = 300)
print(FeaturePlot(tum.integrated.ns, features = markers1, reduction = reduction, pt.size = 0.01) + NoAxes())
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_markers2_method_9.png"), width=16, height=15, units = "in", res = 300)
print(FeaturePlot(tum.integrated.ns, features = markers2, reduction = reduction, pt.size = 0.01) + NoAxes())
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_markers3_method_9.png"), width=16, height=15, units = "in", res = 300)
print(FeaturePlot(tum.integrated.ns, features = markers3, reduction = reduction, pt.size = 0.01) + NoAxes())
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_markers4_method_9.png"), width=16, height=15, units = "in", res = 300)
print(FeaturePlot(tum.integrated.ns, features = markers4, reduction = reduction, pt.size = 0.01) + NoAxes())
dev.off()


tum = read_rds("~/projects/bcap/analysis/20201212_final2/figure_1_v1/final_use.rds")


# heatmap -------

tum_sub_seu = read.csv("dims_15k_20_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.085_vst_method_9.csv")
markers_v2 = tum_sub_seu %>% dplyr::filter(gene %in% tum_sub_seu$gene)

top_markers = markers_v2 %>%
  dplyr::group_by(cluster) %>%
  arrange(cluster, avg_logFC) %>% 
  top_n(25, avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% unique() %>% as.character()

DefaultAssay(tum) <- "integrated"
Idents(tum) <- "integrated_snn_res.0.085"

tum_sub1 = AverageExpression(tum, features = top_markers, return.seurat = T)
#levels(tum_sub1$orig.ident) <- levels(Idents(tum_sub))
p1 = DoHeatmap(tum_sub1, features = top_markers, assay = "integrated", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = colors_dark, raster = T)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))

png("heatmap_0.085_FIBROaVG_v1.png", width=4, height=12, units = "in", res = 300) #
print(p1)
dev.off()

p1 = DoHeatmap(tum, features = top_markers, cells = sampled.cells, assay = "integrated", label = FALSE,
               slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = colors_dark, raster = F)  +
  scale_fill_gradient2(low = "#29A7B8", mid = "white", high = "#E85041", na.value = "00000000")
#ggsave(plot = p1, filename = "heatmap_cellstate_myecells_v2.png", device = "png", width = 20, height = 18)

png("heatmap_0.085_FIBRO_100CELLS_v1.png", width=8, height=12, units = "in", res = 300) #
print(p1)
dev.off()

tum = ScaleData(tum, assay = "RNA")
p1 = DoHeatmap(tum, features = top_markers, cells = sampled.cells, assay = "RNA", label = FALSE,
               slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = colors_dark, raster = F)  +
  scale_fill_gradient2(low = "#29A7B8", mid = "white", high = "#E85041", na.value = "00000000")
#ggsave(plot = p1, filename = "heatmap_cellstate_myecells_v2.png", device = "png", width = 20, height = 18)

png("heatmap_0.085_FIBRO_100CELLS_rna_v1.png", width=8, height=12, units = "in", res = 300) #
print(p1)
dev.off()

# ^^ Final object -------
write_rds(tum, "pericytes_final_use.rds")

#clus_out = read_rds("pericytes_final_use.rds")



























































































































