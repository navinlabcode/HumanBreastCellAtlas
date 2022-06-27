setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20200210_final1/")
odir = "02_markers"; dir.create(odir)
#plan("multiprocess", workers = 100)
options(future.globals.maxSize = 50000 * 1024^2)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20200210_final1/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20200210_final1/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20200210_final1/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20200210_final1/00_functions_final1.R")

# core ---
setwd("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20200210_final1/")
odir = "02_markers"; dir.create(odir)
#plan("multiprocess", workers = 100)
options(future.globals.maxSize = 50000 * 1024^2)
library(future)
plan("multiprocess", workers = 10)
plan(strategy = "multicore", workers = 20)
wd = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20200210_final1/"
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20200210_final1/00_load_packages.R")
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20200210_final1/00_functions_final1.R")


# set variables ----------------------------------------------------------------------------------------------------------------
dims = 20; k.param = 31; cols = colors_dark; reduction = "umap"
fts1 = c("EPCAM", "PTPRC", "PECAM1")

# find markers ----------------------------------------------------------------------------------------------------------------
# ** Markers function --------
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20200210_final1/02_markers/"
#save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20200210_final1/02_markers/"
char_ident = "integrated_snn_res.0.6"  #char_ident = "clusters1"
clus_out = read_rds(path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20200210_final1/02_cluster_integration_try1/umap20k_31tum.integrated.ns1_method09.rds") 
mark_out = func_markers_vst(clus_out = clus_out, char_ident = char_ident, fts1 = fts1, dims = dims, k.param = k.param, cols = colors_dark, 
                            save_path = save_path)

clus_out = ScaleData(clus_out, assay = "RNA")

# rename clusters ----------------------------------------------------------------------------------------------------------------

clus_out$clusters0 = as.character(clus_out$integrated_snn_res.0.6)
clus_out$clusters1 = as.character(clus_out$integrated_snn_res.0.6)
clus_out$clusters2 = as.character(clus_out$integrated_snn_res.0.6)
clus_out$clusters3 = as.character(clus_out$integrated_snn_res.0.6)

clus_out$clusters1[clus_out$clusters0 %in% c(0)] = "fibro"
clus_out$clusters1[clus_out$clusters0 %in% c(1)] = "fibro"
clus_out$clusters1[clus_out$clusters0 %in% c(2)] = "epi_basal"
clus_out$clusters1[clus_out$clusters0 %in% c(3)] = "epi_lumhr"
clus_out$clusters1[clus_out$clusters0 %in% c(4)] = "epi_lumsec"
clus_out$clusters1[clus_out$clusters0 %in% c(5)] = "endo_vas"
clus_out$clusters1[clus_out$clusters0 %in% c(6)] = "t_cd8"
clus_out$clusters1[clus_out$clusters0 %in% c(7)] = "endo_vas"
clus_out$clusters1[clus_out$clusters0 %in% c(8)] = "t_cd4"
clus_out$clusters1[clus_out$clusters0 %in% c(9)] = "pericyte"
clus_out$clusters1[clus_out$clusters0 %in% c(10)] = "macro"
clus_out$clusters1[clus_out$clusters0 %in% c(11)] = "pericyte"
clus_out$clusters1[clus_out$clusters0 %in% c(12)] = "endo_lymph"
clus_out$clusters1[clus_out$clusters0 %in% c(13)] = "macro"
clus_out$clusters1[clus_out$clusters0 %in% c(14)] = "doublet"
clus_out$clusters1[clus_out$clusters0 %in% c(15)] = "nk"
clus_out$clusters1[clus_out$clusters0 %in% c(16)] = "fibro"
clus_out$clusters1[clus_out$clusters0 %in% c(17)] = "b"
clus_out$clusters1[clus_out$clusters0 %in% c(18)] = "endo_vas"
clus_out$clusters1[clus_out$clusters0 %in% c(19)] = "macro"
clus_out$clusters1[clus_out$clusters0 %in% c(20)] = "fibro"
clus_out$clusters1[clus_out$clusters0 %in% c(21)] = "doublet"
clus_out$clusters1[clus_out$clusters0 %in% c(22)] = "endo_vas"
clus_out$clusters1[clus_out$clusters0 %in% c(23)] = "macro"
clus_out$clusters1[clus_out$clusters0 %in% c(24)] = "doublet"
clus_out$clusters1[clus_out$clusters0 %in% c(25)] = "apocrine"

clus_out$clusters2[clus_out$clusters0 %in% c(0)] = "fibro-1"
clus_out$clusters2[clus_out$clusters0 %in% c(1)] = "fibro-2"
clus_out$clusters2[clus_out$clusters0 %in% c(2)] = "epi_basal"
clus_out$clusters2[clus_out$clusters0 %in% c(3)] = "epi_lumhr"
clus_out$clusters2[clus_out$clusters0 %in% c(4)] = "epi_lumsec"
clus_out$clusters2[clus_out$clusters0 %in% c(5)] = "endo_vas-1"
clus_out$clusters2[clus_out$clusters0 %in% c(6)] = "t_cd8"
clus_out$clusters2[clus_out$clusters0 %in% c(7)] = "endo_vas-2"
clus_out$clusters2[clus_out$clusters0 %in% c(8)] = "t_cd4"
clus_out$clusters2[clus_out$clusters0 %in% c(9)] = "pericyte-1"
clus_out$clusters2[clus_out$clusters0 %in% c(10)] = "macro-1"
clus_out$clusters2[clus_out$clusters0 %in% c(11)] = "pericyte-2"
clus_out$clusters2[clus_out$clusters0 %in% c(12)] = "endo_lymph"
clus_out$clusters2[clus_out$clusters0 %in% c(13)] = "macro-2"
clus_out$clusters2[clus_out$clusters0 %in% c(14)] = "doublet-fib-endo"
clus_out$clusters2[clus_out$clusters0 %in% c(15)] = "nk"
clus_out$clusters2[clus_out$clusters0 %in% c(16)] = "fibro-3"
clus_out$clusters2[clus_out$clusters0 %in% c(17)] = "b"
clus_out$clusters2[clus_out$clusters0 %in% c(18)] = "endo_vas-3"
clus_out$clusters2[clus_out$clusters0 %in% c(19)] = "macro-3"
clus_out$clusters2[clus_out$clusters0 %in% c(20)] = "fibro-4"
clus_out$clusters2[clus_out$clusters0 %in% c(21)] = "doublet-fib-t"
clus_out$clusters2[clus_out$clusters0 %in% c(22)] = "endo_vas-4"
clus_out$clusters2[clus_out$clusters0 %in% c(23)] = "macro-4"
clus_out$clusters2[clus_out$clusters0 %in% c(24)] = "doublet-epi-endo"
clus_out$clusters2[clus_out$clusters0 %in% c(25)] = "apocrine"

clus_out$clusters3[clus_out$clusters0 %in% c(0)] = "stromal"
clus_out$clusters3[clus_out$clusters0 %in% c(1)] = "stromal"
clus_out$clusters3[clus_out$clusters0 %in% c(2)] = "epithelial"
clus_out$clusters3[clus_out$clusters0 %in% c(3)] = "epithelial"
clus_out$clusters3[clus_out$clusters0 %in% c(4)] = "epithelial"
clus_out$clusters3[clus_out$clusters0 %in% c(5)] = "stromal"
clus_out$clusters3[clus_out$clusters0 %in% c(6)] = "immune"
clus_out$clusters3[clus_out$clusters0 %in% c(7)] = "stromal"
clus_out$clusters3[clus_out$clusters0 %in% c(8)] = "immune"
clus_out$clusters3[clus_out$clusters0 %in% c(9)] = "stromal"
clus_out$clusters3[clus_out$clusters0 %in% c(10)] = "immune"
clus_out$clusters3[clus_out$clusters0 %in% c(11)] = "stromal"
clus_out$clusters3[clus_out$clusters0 %in% c(12)] = "stromal"
clus_out$clusters3[clus_out$clusters0 %in% c(13)] = "immune"
clus_out$clusters3[clus_out$clusters0 %in% c(14)] = "doublet"
clus_out$clusters3[clus_out$clusters0 %in% c(15)] = "immune"
clus_out$clusters3[clus_out$clusters0 %in% c(16)] = "stromal"
clus_out$clusters3[clus_out$clusters0 %in% c(17)] = "immune"
clus_out$clusters3[clus_out$clusters0 %in% c(18)] = "stromal"
clus_out$clusters3[clus_out$clusters0 %in% c(19)] = "immune"
clus_out$clusters3[clus_out$clusters0 %in% c(20)] = "stromal"
clus_out$clusters3[clus_out$clusters0 %in% c(21)] = "doublet"
clus_out$clusters3[clus_out$clusters0 %in% c(22)] = "stromal"
clus_out$clusters3[clus_out$clusters0 %in% c(23)] = "immune"
clus_out$clusters3[clus_out$clusters0 %in% c(24)] = "doublet"
clus_out$clusters3[clus_out$clusters0 %in% c(25)] = "stromal"

# save umaps ---
p1 = DimPlot(object = clus_out, reduction = reduction, label = FALSE, pt.size = 0.001, cols = cols, group.by = "clusters1")+  NoAxes()
p2 = DimPlot(object = clus_out, reduction = reduction, label = FALSE, pt.size = 0.001, cols = cols, group.by = "clusters2") +  NoAxes()
p3 = DimPlot(object = clus_out, reduction = reduction, group.by = "procedure",  pt.size = 0.0001, cols = c("tan1","steelblue1")) + NoAxes()
p4 = DimPlot(object = clus_out, reduction = reduction, group.by = "source",  pt.size = 0.0001, cols = c("darkorange","darkorchid")) +  NoAxes()
p5 = DimPlot(object = clus_out, reduction = reduction, group.by = "patient_id", pt.size = 0.0001, cols = colors_dark) +  NoAxes()
p6 = DimPlot(object = clus_out, reduction = reduction, group.by = "sample_id", pt.size = 0.0001, cols = colors_dark) +  NoAxes()
p7 = DimPlot(object = clus_out, reduction = reduction, group.by = "exp_proc", pt.size = 0.0001, cols = colors_dark) +  NoAxes()

# system("pdftoppm input.pdf outputname -png")
# save plots
png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_clusters1_method_9.png"), width=6, height=5, units = "in", res = 300)
print(p1)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_clusters2_method_9.png"), width=8, height=5, units = "in", res = 300)
print(p2)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_procedure_method_9.png"), width=5, height=5, units = "in", res = 300)
print(p3)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_source_method_9.png"), width=6, height=5, units = "in", res = 300)
print(p4)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_patient_method_9.png"), width=8, height=5, units = "in", res = 300)
print(p5)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_sample_id_method_9.png"), width=8, height=5, units = "in", res = 300)
print(p6)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_exp_proc_method_9.png"), width=8, height=5, units = "in", res = 300)
print(p7)
dev.off()

write_rds(x = clus_out, path = paste0(save_path, "clus_out.rds"))
clus_out = read_rds(paste0(save_path, "clus_out.rds"))

# barplots ----------------------------------------------------------------------------------------------------------------

# **Barplots ---------------------------------------------------------------
df_meta = data.frame(clus_out@meta.data, 
                     celltype = clus_out$clusters1, 
                     stringsAsFactors = F) %>% 
  tbl_df() %>% clean_names()

#df_meta = df_meta %>% mutate(orig_pt_ident = factor(orig_pt_ident, levels = c("HR+ Luminal Cells", "Secretory Luminal Cells", "Endothelial",  "Fibroblasts",  "Pericytes", "T cells", "Macrophages",
#                                                              "TumA", "TumB", "TumC", "TumD")))
df_meta1 = df_meta %>% mutate(patient_id =  factor(patient_id, levels = unique(patient_id[order(source, procedure)]), ordered = TRUE))
# fn = factor(f, levels=unique(f[order(a,b,f)]), ordered=TRUE)

# create a summary table
func_bar(var_df = df_meta1, var_x = "patient_id", var_fill = "source", var_pos = "fill", var_cols = colors_dark, celltype = "source")
func_bar(var_df = df_meta1, var_x = "patient_id", var_fill = "exp_proc", var_pos = "fill", var_cols = colors_dark[10:15], celltype = "exp_proc")
func_bar(var_df = df_meta1, var_x = "patient_id", var_fill = "procedure", var_pos = "fill", var_cols = colors_dark[4:7], celltype = "procedure")

func_bar(var_df = df_meta1, var_x = "sample_id", var_fill = "source", var_pos = "fill", var_cols = colors_dark, celltype = "source")
func_bar(var_df = df_meta1, var_x = "sample_id", var_fill = "exp_proc", var_pos = "fill", var_cols = colors_dark[10:15], celltype = "exp_proc")
func_bar(var_df = df_meta1, var_x = "sample_id", var_fill = "procedure", var_pos = "fill", var_cols = colors_dark[4:7], celltype = "procedure")

func_bar(var_df = df_meta1, var_x = "patient_id", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters1")
func_bar(var_df = df_meta1, var_x = "patient_id", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters1")
func_bar(var_df = df_meta1, var_x = "sample_id", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters1")
func_bar(var_df = df_meta1, var_x = "sample_id", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters1")
func_bar(var_df = df_meta1, var_x = "procedure", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters1")
func_bar(var_df = df_meta1, var_x = "procedure", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters1")
func_bar(var_df = df_meta1, var_x = "source", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters1")
func_bar(var_df = df_meta1, var_x = "source", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters1")
func_bar(var_df = df_meta1, var_x = "exp_proc", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters1")
func_bar(var_df = df_meta1, var_x = "exp_proc", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters1")

df_meta = data.frame(clus_out@meta.data, 
                     celltype = clus_out$clusters2, 
                     stringsAsFactors = F) %>% 
  tbl_df() %>% clean_names()

#df_meta = df_meta %>% mutate(orig_pt_ident = factor(orig_pt_ident, levels = c("HR+ Luminal Cells", "Secretory Luminal Cells", "Endothelial",  "Fibroblasts",  "Pericytes", "T cells", "Macrophages",
#                                                              "TumA", "TumB", "TumC", "TumD")))
df_meta1 = df_meta %>% mutate(patient_id =  factor(patient_id, levels = unique(patient_id[order(source, procedure)]), ordered = TRUE))

func_bar(var_df = df_meta1, var_x = "patient_id", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters2")
func_bar(var_df = df_meta1, var_x = "patient_id", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters2")
func_bar(var_df = df_meta1, var_x = "sample_id", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters2")
func_bar(var_df = df_meta1, var_x = "sample_id", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters2")
func_bar(var_df = df_meta1, var_x = "procedure", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters2")
func_bar(var_df = df_meta1, var_x = "procedure", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters2")
func_bar(var_df = df_meta1, var_x = "source", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters2")
func_bar(var_df = df_meta1, var_x = "source", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters2")
func_bar(var_df = df_meta1, var_x = "exp_proc", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters2")
func_bar(var_df = df_meta1, var_x = "exp_proc", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters2")

df_meta = data.frame(clus_out@meta.data, 
                     celltype = clus_out$clusters0, 
                     stringsAsFactors = F) %>% 
  tbl_df() %>% clean_names()

#df_meta = df_meta %>% mutate(orig_pt_ident = factor(orig_pt_ident, levels = c("HR+ Luminal Cells", "Secretory Luminal Cells", "Endothelial",  "Fibroblasts",  "Pericytes", "T cells", "Macrophages",
#                                                              "TumA", "TumB", "TumC", "TumD")))
df_meta1 = df_meta %>% mutate(patient_id =  factor(patient_id, levels = unique(patient_id[order(source, procedure)]), ordered = TRUE))

func_bar(var_df = df_meta1, var_x = "patient_id", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters0")
func_bar(var_df = df_meta1, var_x = "patient_id", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters0")
func_bar(var_df = df_meta1, var_x = "sample_id", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters0")
func_bar(var_df = df_meta1, var_x = "sample_id", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters0")
func_bar(var_df = df_meta1, var_x = "procedure", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters0")
func_bar(var_df = df_meta1, var_x = "procedure", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters0")
func_bar(var_df = df_meta1, var_x = "source", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters0")
func_bar(var_df = df_meta1, var_x = "source", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters0")
func_bar(var_df = df_meta1, var_x = "exp_proc", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters0")
func_bar(var_df = df_meta1, var_x = "exp_proc", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters0")

df_meta = data.frame(clus_out@meta.data, 
                     celltype = clus_out$clusters3, 
                     stringsAsFactors = F) %>% 
  tbl_df() %>% clean_names()

#df_meta = df_meta %>% mutate(orig_pt_ident = factor(orig_pt_ident, levels = c("HR+ Luminal Cells", "Secretory Luminal Cells", "Endothelial",  "Fibroblasts",  "Pericytes", "T cells", "Macrophages",
#                                                              "TumA", "TumB", "TumC", "TumD")))
df_meta1 = df_meta %>% mutate(patient_id =  factor(patient_id, levels = unique(patient_id[order(source, procedure)]), ordered = TRUE))

func_bar(var_df = df_meta1, var_x = "exp_proc", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters3")
func_bar(var_df = df_meta1, var_x = "exp_proc", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters3")

func_bar(var_df = df_meta1, var_x = "patient_id", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters3")
func_bar(var_df = df_meta1, var_x = "patient_id", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters3")
func_bar(var_df = df_meta1, var_x = "sample_id", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters3")
func_bar(var_df = df_meta1, var_x = "sample_id", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters3")
func_bar(var_df = df_meta1, var_x = "procedure", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters3")
func_bar(var_df = df_meta1, var_x = "procedure", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters3")
func_bar(var_df = df_meta1, var_x = "source", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters3")
func_bar(var_df = df_meta1, var_x = "source", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters3")
func_bar(var_df = df_meta1, var_x = "exp_proc", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters3")
func_bar(var_df = df_meta1, var_x = "exp_proc", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters3")

# gene plots ------
markers = read.csv("/volumes/lab/users/tapsi/projects/bcap/analysis/20200210_final1/02_markers/dims_20k_31_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.6_vst.csv")
genes = markers$gene
genes  = as.character(genes)
gg_feature_plot(clus_out, genes = unique(genes)[1], reduction = "umap", save_folder = save_path)
