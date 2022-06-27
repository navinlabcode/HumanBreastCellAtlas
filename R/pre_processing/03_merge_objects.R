# ---------
# setwd -------

# geo ---------
setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_merge_objects"; dir.create(odir)
library(future)
plan("multiprocess", workers = 30)
options(future.globals.maxSize = Inf)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_merge_objects/"

#cl = c("#c385ec", "#d98f00", "#00bfc1", "#f47c55", "#ff8ca5", "#405f8f")

# core ---------
setwd("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_merge_objects"; dir.create(odir)
library(future)
plan("multiprocess", workers = 30)
#plan("multisession", workers = 20)
options(future.globals.maxSize = Inf)
setwd(paste0("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_merge_objects"

# METHOD 9 ------------- 
###  logT and integrated by variable (var_split = "patient"), no scale by any vars
# ** NORM 
library(future)
plan("multiprocess", workers = 20)
#plan("sequential")

# Read objects --------
# geo
seu_object = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_remove_doublets_recluster/tum.integrated.ns_subsetted_1.rds")
seu_object_old = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20200210_final1/figure_1_v2/tum.rds")

# core
seu_object = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_remove_doublets_recluster/tum.integrated.ns_subsetted_1.rds")
seu_object_old = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20200210_final1/figure_1_v2/tum.rds")

# define cluster group 
seu_object@meta.data$cluster_updated = seu_object@meta.data$clusters1 
DefaultAssay(seu_object) <- "RNA"
seu_object@meta.data$batch = "b2" 
seu_object1 = DietSeurat(seu_object, counts = TRUE, scale.data = FALSE, dimreducs = TRUE, assays = "RNA")

seu_object_old@meta.data$cluster_updated = seu_object_old@meta.data$clusters1
seu_object_old@meta.data$batch = "b1" 
DefaultAssay(seu_object_old) <- "RNA"
Idents(seu_object_old) <- "tissue_location.x"
seu_object_old = subset(seu_object_old, idents = c("ipsi"),invert = TRUE) # remove ipsi samples
seu_object_old1 = DietSeurat(seu_object_old, counts = TRUE, scale.data = FALSE, dimreducs = TRUE, assays = "RNA")


# merge objects ---------------------------

# with integrated
tum.merged0 = merge(seu_object_old, seu_object, add.cell.ids = c("r1", "r2"), project = "hbca_cell")
write_rds(tum.merged0, glue("{save_path}tum.merged0.rds"))

# without integrated
tum.merged = merge(seu_object_old1, seu_object1, add.cell.ids = c("r1", "r2"), project = "hbca_cell")
write_rds(tum.merged, glue("{save_path}tum.merged.rds"))


# read objects ------------
tum.merged = read_rds(glue("{save_path}/tum.merged.rds"))
tum.merged = read_rds(glue("{save_path}/tum.merged0.rds"))

#dims = 18; k.param = 30; cols = colors_dark;var_split = "batch"
dims = 20; k.param = 30; cols = colors_dark;var_split = "batch"
fts1 = c("nFeature_RNA","nCount_RNA", "percent.mito","percent.rb", "percent.sg1", "percent.sg2", "percent.sg3")
fts2 = c("rna_EPCAM","rna_PTPRC","rna_PECAM1","rna_DCN","rna_KRT5",  "rna_LTF", "rna_KRT17", "rna_HPGD")
fts3 = c("rna_PDE3B","rna_ACACB","rna_ACSM1","rna_ACTA2","rna_CD4",  "rna_NKG7", "rna_CD19", "rna_MS4A1", "rna_RGS5")
fts4 = c("rna_VWF","rna_PRSS1","rna_CCL21","rna_FAP","COL1A1",  "rna_COL6A1", "rna_ESR2", "rna_KRT8", "rna_SERHL2")

clus_out = func_cluster_vst(tum_int = tum.merged, var_split = "batch", var_scale = "NULL", dims = dims, k.param = k.param, cols = colors_dark,  
                            markers1 = fts1, markers2 = fts2, markers3 = fts3, markers4 = fts4, save_path = save_path)
clus_out = func_cluster_rpca(tum_int = tum.merged, var_split = "batch", var_scale = "NULL", dims = dims, k.param = k.param, cols = colors_dark,  
                             markers1 = fts1, markers2 = fts2, markers3 = fts3, markers4 = fts4, save_path = save_path)
char_ident = "integrated_snn_res.0.2"

tum.integrated.ns = read_rds(path = glue("{save_path}/dim_20k_30_NULL_tum.integrated.ns2_method_9rcp.rds"))
Idents(tum.integrated.ns) <- "sample_id"
tum.integrated.ns = subset(tum.integrated.ns, idents = c("hbca_c02", "hbca_c11", "hbca_c04", "hbca_c08", "hbca_c13", "hbca_c06", "hbca_c15"), invert = T)
tum.integrated.ns = subset(tum.integrated.ns, idents = c("hbca_c15"), invert = T)
write_rds(tum.integrated.ns, "dim20_30_0.1_1.75_final1.rds")

# FINAL OBJECT MERGED ------
# tum.integrated.ns = read_rds("dim_20k_30_NULL_tum.integrated.ns2_method_9rcp_filtered.rds")
# tum.integrated.ns = read_rds("dim_20k_30_NULL_tum.integrated.ns2_method_9rcp.rds")
tum.integrated.ns = read_rds("dim_20k_30_NULL_final1.rds")

# ** Markers function --------
#save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201120_final_nuc2/02_cluster_integration_core/"
char_ident = "integrated_snn_res.0.6"  #integrated_snn_res.0.1
char_ident = "cluster_updated"  #integrated_snn_res.0.1
clus_out = read_rds("dim_20k_30_NULL_tum.integrated.ns2_method_9rcp_filtered.rds")
fts2 = c("rna_EPCAM","rna_PTPRC","rna_PECAM1","rna_DCN","rna_KRT5",  "rna_LTF", "rna_KRT17", "rna_CD68")
#clus_out = read_rds(path = glue("{save_path}/dim_{dims}k_{k.param}_{var_split}_tum.integrated.ns2_method_9.rds"))
mark_out = func_markers_vst(clus_out = clus_out, char_ident = char_ident, fts1 = fts2, dims = dims, k.param = k.param, cols = colors_dark, 
                            save_path = save_path)
char_ident = "integrated_snn_res.0.6"  #integrated_snn_res.0.1
mark_out = func_markers_vst(clus_out = clus_out, char_ident = char_ident, fts1 = fts2, dims = dims, k.param = k.param, cols = colors_dark, 
                            save_path = save_path)

tum.merged$clusters_final0 = tum.merged$cluster_updated
tum.merged$clusters_final0[tum.merged$cluster_updated %in% c("pericytes")] = "pericyte"
tum.merged$clusters_final0[tum.merged$cluster_updated %in% c("apocrine")] = "epi_lumhr"
tum.merged$clusters_final0[tum.merged$cluster_updated %in% c("t_cd4")] = "t"
tum.merged$clusters_final0[tum.merged$cluster_updated %in% c("t_cd8")] = "t"

p0 = DimPlot(object = tum.integrated.ns, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = colors_dark, group.by = "cluster_updated", shuffle = T)+  NoAxes()
png(paste0(save_path,"uamp_integrated_dims20_integrated_snn_res.0.6_dist0.3.png"), width=8, height=5, units = "in", res = 300)
print(p0)
dev.off()
p0 = DimPlot(object = tum.integrated.ns, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = colors_dark, group.by = "integrated_snn_res.0.8", shuffle = T)+  NoAxes()
png(paste0(save_path,"uamp_integrated_dims20_integrated_snn_res.0.8.png"), width=8, height=5, units = "in", res = 300)
print(p0)
dev.off()
p0 = DimPlot(object = tum.integrated.ns, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = colors_dark, group.by = "integrated_snn_res.0.4", shuffle = T)+  NoAxes()
png(paste0(save_path,"uamp_integrated_dims20_integrated_snn_res.0.6.png"), width=8, height=5, units = "in", res = 300)
print(p0)
dev.off()

p0 = DimPlot(object = tum.integrated.ns, reduction = "umap", label = TRUE, pt.size = 0.000001, group.by = "integrated_snn_res.4", shuffle = T)+  NoAxes()
png(paste0(save_path,"uamp_integrated_dims20_integrated_snn_res.0.6.png"), width=8, height=5, units = "in", res = 300)
print(p0)
dev.off()

#tum.integrated.ns1 = subset(x = tum.integrated.ns, subset = nCount_RNA > 600)
p1 = VlnPlot(object = tum.integrated.ns, pt.size = 0, cols = colors_dark, group.by = "integrated_snn_res.0.6", features = c("nCount_RNA", "nFeature_RNA"))
message("aggregating data ...2/6")
tum_count_agg = aggregate(nCount_RNA~integrated_snn_res.0.6, tum.integrated.ns@meta.data, summary)
tum_RNA_agg = aggregate(nFeature_RNA~integrated_snn_res.0.6, tum.integrated.ns@meta.data, summary)
tum_percent.rb_agg = aggregate(percent.rb~integrated_snn_res.0.6, tum.integrated.ns@meta.data, summary)
tum_percent.mito_agg = aggregate(percent.mito~integrated_snn_res.0.6, tum.integrated.ns@meta.data, summary)
tum_percent.sg1_agg = aggregate(percent.sg1~integrated_snn_res.0.6, tum.integrated.ns@meta.data, summary)
tum_percent.sg2_agg = aggregate(percent.sg2~integrated_snn_res.0.6, tum.integrated.ns@meta.data, summary)
tum_percent.sg3_agg = aggregate(percent.sg3~integrated_snn_res.0.6, tum.integrated.ns@meta.data, summary)
tum_summary = cbind(tum_count_agg, tum_RNA_agg, tum_percent.rb_agg, tum_percent.mito_agg, tum_percent.sg1_agg, tum_percent.sg2_agg, tum_percent.sg3_agg)

# save the aggreate summary of each filter
message("saving csv ...3/6")
write.csv(tum_summary, glue("{save_path}/tum_summary_integrated_snn_res.0.6.csv"))

# clean up small clusters ----
VlnPlot(object = tum.integrated.ns, pt.size = 0, cols = colors_dark, group.by = "integrated_snn_res.0.6", features = c("rna_KRT14", "rna_MUCL1"))

clus_out1 = RunUMAP(clus_out, dims = 1:22)
p9 = DimPlot(object = clus_out1, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = colors_dark, group.by = "cluster_updated", shuffle = T)+  NoAxes()
png(paste0(save_path,"uamp_integrated_dims22.png"), width=8, height=5, units = "in", res = 300)
print(p9)
dev.off()


p9 = DimPlot(object = tum.merged, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = colors_dark, group.by = "integrated_snn_res.0.6", shuffle = T)+  NoAxes()
png(paste0(save_path,"/", reduction,dims,"k_", k.param, "INTumap_merged1_RNA_integrated_snn_res.0.4.png"), width=8, height=5, units = "in", res = 300)
print(p9)
dev.off()

p0 = DimPlot(object = tum.merged, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = colors_dark, group.by = "integrated_snn_res.0.2", shuffle = T)+  NoAxes()
p000 = DimPlot(object = tum.merged, reduction = "umap", label = F, pt.size = 0.000001, cols = colors_dark, group.by = "clusters_final0", shuffle = T)+  NoAxes()
p00 = DimPlot(object = tum.merged, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = c("grey", "maroon"), order = "b1", group.by = "batch", shuffle = TRUE)+  NoAxes()
p01 = DimPlot(object = tum.merged, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = c("grey", "maroon"), order = "b2", group.by = "batch", shuffle = TRUE)+  NoAxes()
p00g = p00 + p01
p1 = DimPlot(object = tum.merged, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = colors_dark, group.by = "clusters_final0")+  NoAxes()
p2 = DimPlot(object = tum.merged, reduction = "umap", label = TRUE, pt.size = 0.001, cols = colors_dark, group.by = "clusters_final0", split.by = "batch")+  NoAxes()
p3 = DimPlot(object = tum.merged, reduction = "umap", label = F, pt.size = 0.001, cols = colors_dark, group.by = "clusters_final0", split.by = "exp_proc")+  NoAxes()
p4 = DimPlot(object = tum.merged, reduction = "umap", label = F, pt.size = 0.001, cols = colors_dark, group.by = "clusters_final0", split.by = "procedure")+  NoAxes()
p5 = DimPlot(object = tum.merged, reduction = "umap", label = F, pt.size = 0.001, cols = colors_dark, group.by = "clusters_final0", split.by = "source")+  NoAxes()
p6 = DimPlot(object = tum.merged, reduction = "umap", label = F, pt.size = 0.001, cols = colors_dark, group.by = "exp_proc", split.by = "tissue_location.x")+  NoAxes()
p7 = DimPlot(object = tum.merged, reduction = "umap", label = F, pt.size = 0.001, cols = colors_dark, group.by = "exp_proc", split.by = "tissue_location.y")+  NoAxes()

#tum.merged_umap = read_rds("tum.merged_umap.rds")
png(paste0(save_path,"/", reduction,dims,"k_", k.param, "INTumap_merged1_RNA_snn_res.0.6.png"), width=8, height=5, units = "in", res = 300)
print(p0)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "INTumap_merged1_clusters_final0.png"), width=10, height=5, units = "in", res = 300)
print(p000)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "INTumap_merged1_batch.png"), width=13, height=5, units = "in", res = 300)
print(p00g)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "INTumap_all_cluster_updated.png"), width=8, height=5, units = "in", res = 300)
print(p1)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "INTumap_all_cluster_updatedsplit_batch.png"), width=14, height=6.5, units = "in", res = 300)
print(p2)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "INTumap_all_cluster_updatedsplit_exp_proc.png"), width=24, height=5.5, units = "in", res = 300)
print(p3)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "INTumap_all_cluster_updatedsplit_procedure.png"), width=13, height=5.5, units = "in", res = 300)
print(p4)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "INTumap_all_cluster_updatedsplit_source.png"), width=14, height=5.5, units = "in", res = 300)
print(p5)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "INTumap_all_cluster_updated_split_exp_proc_tissue_location.png"), width=18, height=5.5, units = "in", res = 300)
print(p6)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "INTumap_all_cluster_updated_split_exp_proc_tissue_location_y.png"), width=18, height=5.5, units = "in", res = 300)
print(p7)
dev.off()











# MERGED MATRICES -----

# clustering on merged -----

tum.merged = read_rds(glue("{save_path}/tum.merged0.rds"))
tum.merged = tum.merged0
dims = 22; k.param = 30; cols = colors_dark;var_split = "batch";reduction = "umap"
tum.merged = NormalizeData(tum.merged, normalization.method = "LogNormalize", scale.factor = 10000)
tum.merged = FindVariableFeatures(tum.merged, selection.method = "vst", nfeatures = 2000)
tum.merged = ScaleData(tum.merged, assay = "RNA", vars.to.regress = NULL) # no scaling # no scaling
tum.merged = RunPCA(tum.merged, npcs = 75)
tum.merged = RunUMAP(object = tum.merged, reduction = "pca", dims = 1:dims)
write_rds(tum.merged, "tum.merged_umap_final_log.rds")

tum.merged = FindNeighbors(object = tum.merged, dims = 1:dims, reduction = "pca", k.param = k.param)
tum.merged = FindClusters(object = tum.merged, resolution = c(1, 0.8, 0.6,0.4,0.3,0.2, 0.15,0.17, 0.1, 0.085, 0.090, 0.075, 0.05, 0.023), verbose = TRUE)
write_rds(tum.merged, "tum.merged_umap_final_log_clus.rds")

tum.merged = read_rds(glue("{save_path}/tum.merged_umap_final_log_clus.rds"))
Idents(tum.merged) <- "cluster_updated"
tum.merged$clusters_final0 = tum.merged$cluster_updated
tum.merged$clusters_final0[tum.merged$cluster_updated %in% c("pericytes")] = "pericyte"
tum.merged$clusters_final0[tum.merged$cluster_updated %in% c("apocrine")] = "epi_lumhr"
tum.merged$clusters_final0[tum.merged$cluster_updated %in% c("t_cd4")] = "t"
tum.merged$clusters_final0[tum.merged$cluster_updated %in% c("t_cd8")] = "t"

p0 = DimPlot(object = tum.merged, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = colors_dark, group.by = "RNA_snn_res.0.6", shuffle = T)+  NoAxes()
p00 = DimPlot(object = tum.merged, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = c("grey", "maroon"), order = "b1", group.by = "batch", shuffle = TRUE)+  NoAxes()
p01 = DimPlot(object = tum.merged, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = c("grey", "maroon"), order = "b2", group.by = "batch", shuffle = TRUE)+  NoAxes()
p00g = p00 + p01
p1 = DimPlot(object = tum.merged, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = colors_dark, group.by = "clusters_final0")+  NoAxes()
p2 = DimPlot(object = tum.merged, reduction = "umap", label = TRUE, pt.size = 0.001, cols = colors_dark, group.by = "clusters_final0", split.by = "batch")+  NoAxes()
p3 = DimPlot(object = tum.merged, reduction = "umap", label = F, pt.size = 0.001, cols = colors_dark, group.by = "clusters_final0", split.by = "exp_proc")+  NoAxes()
p4 = DimPlot(object = tum.merged, reduction = "umap", label = F, pt.size = 0.001, cols = colors_dark, group.by = "clusters_final0", split.by = "procedure")+  NoAxes()
p5 = DimPlot(object = tum.merged, reduction = "umap", label = F, pt.size = 0.001, cols = colors_dark, group.by = "clusters_final0", split.by = "source")+  NoAxes()
p6 = DimPlot(object = tum.merged, reduction = "umap", label = F, pt.size = 0.001, cols = colors_dark, group.by = "exp_proc", split.by = "tissue_location.x")+  NoAxes()
p7 = DimPlot(object = tum.merged, reduction = "umap", label = F, pt.size = 0.001, cols = colors_dark, group.by = "exp_proc", split.by = "tissue_location.y")+  NoAxes()

#tum.merged_umap = read_rds("tum.merged_umap.rds")
png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_merged1_RNA_snn_res.0.6.png"), width=8, height=5, units = "in", res = 300)
print(p0)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_merged1_batch.png"), width=13, height=5, units = "in", res = 300)
print(p00g)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_all_cluster_updated.png"), width=6, height=5, units = "in", res = 300)
print(p1)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_all_cluster_updatedsplit_batch.png"), width=14, height=6.5, units = "in", res = 300)
print(p2)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_all_cluster_updatedsplit_exp_proc.png"), width=24, height=5.5, units = "in", res = 300)
print(p3)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_all_cluster_updatedsplit_procedure.png"), width=13, height=5.5, units = "in", res = 300)
print(p4)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_all_cluster_updatedsplit_source.png"), width=14, height=5.5, units = "in", res = 300)
print(p5)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_all_cluster_updated_split_exp_proc_tissue_location.png"), width=18, height=5.5, units = "in", res = 300)
print(p6)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_all_cluster_updated_split_exp_proc_tissue_location_y.png"), width=18, height=5.5, units = "in", res = 300)
print(p7)
dev.off()


# clustering on integrated merged -----

tum.merged = read_rds(glue("{save_path}/tum.merged0.rds"))
dims = 22; k.param = 30; cols = colors_dark;var_split = "batch";reduction = "umap"
tum.merged = NormalizeData(tum.merged, normalization.method = "LogNormalize", scale.factor = 10000)
tum.merged = FindVariableFeatures(tum.merged, selection.method = "vst", nfeatures = 2000)
tum.merged = ScaleData(tum.merged, assay = "RNA", vars.to.regress = NULL) # no scaling # no scaling
tum.merged = RunPCA(tum.merged, npcs = 75)
tum.merged = RunUMAP(object = tum.merged, reduction = "pca", dims = 1:dims)
write_rds(tum.merged, "tum.int_merged_umap_final_log22_30.rds")

tum.merged = FindNeighbors(object = tum.merged, dims = 1:dims, reduction = "pca", k.param = k.param)
tum.merged = FindClusters(object = tum.merged, resolution = c(1, 0.8, 0.6,0.4,0.3,0.2, 0.15,0.17, 0.1, 0.085, 0.090, 0.075, 0.05, 0.023), verbose = TRUE)
write_rds(tum.merged, "tum.int_merged_umap_final_log_clus22_30.rds")

tum.merged = read_rds(glue("{save_path}/tum.int_merged_umap_final_log_clus22_30.rds"))
Idents(tum.merged) <- "cluster_updated"
tum.merged$clusters_final0 = tum.merged$cluster_updated
tum.merged$clusters_final0[tum.merged$cluster_updated %in% c("pericytes")] = "pericyte"
tum.merged$clusters_final0[tum.merged$cluster_updated %in% c("apocrine")] = "epi_lumhr"
tum.merged$clusters_final0[tum.merged$cluster_updated %in% c("t_cd4")] = "t"
tum.merged$clusters_final0[tum.merged$cluster_updated %in% c("t_cd8")] = "t"

p0 = DimPlot(object = tum.merged, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = colors_dark, group.by = "integrated_snn_res.0.2", shuffle = T)+  NoAxes()
p00 = DimPlot(object = tum.merged, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = c("grey", "maroon"), order = "b1", group.by = "batch", shuffle = TRUE)+  NoAxes()
p01 = DimPlot(object = tum.merged, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = c("grey", "maroon"), order = "b2", group.by = "batch", shuffle = TRUE)+  NoAxes()
p00g = p00 + p01
p1 = DimPlot(object = tum.merged, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = colors_dark, group.by = "cluster_updated")+  NoAxes()
p2 = DimPlot(object = tum.merged, reduction = "umap", label = TRUE, pt.size = 0.001, cols = colors_dark, group.by = "cluster_updated", split.by = "batch")+  NoAxes()
p3 = DimPlot(object = tum.merged, reduction = "umap", label = F, pt.size = 0.001, cols = colors_dark, group.by = "cluster_updated", split.by = "exp_proc")+  NoAxes()
p4 = DimPlot(object = tum.merged, reduction = "umap", label = F, pt.size = 0.001, cols = colors_dark, group.by = "cluster_updated", split.by = "procedure")+  NoAxes()
p5 = DimPlot(object = tum.merged, reduction = "umap", label = F, pt.size = 0.001, cols = colors_dark, group.by = "cluster_updated", split.by = "source")+  NoAxes()
p6 = DimPlot(object = tum.merged, reduction = "umap", label = F, pt.size = 0.001, cols = colors_dark, group.by = "exp_proc", split.by = "tissue_location.x")+  NoAxes()
p7 = DimPlot(object = tum.merged, reduction = "umap", label = F, pt.size = 0.001, cols = colors_dark, group.by = "exp_proc", split.by = "tissue_location.y")+  NoAxes()

#tum.merged_umap = read_rds("tum.merged_umap.rds")
png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_int_merged1_integrated_snn_res.0.2.png"), width=8, height=5, units = "in", res = 300)
print(p0)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_int_merged1_batch.png"), width=13, height=5, units = "in", res = 300)
print(p00g)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_int_all_cluster_updated.png"), width=6, height=5, units = "in", res = 300)
print(p1)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_int_all_cluster_updatedsplit_batch.png"), width=14, height=6.5, units = "in", res = 300)
print(p2)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_int_all_cluster_updatedsplit_exp_proc.png"), width=24, height=5.5, units = "in", res = 300)
print(p3)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_int_all_cluster_updatedsplit_procedure.png"), width=13, height=5.5, units = "in", res = 300)
print(p4)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_int_all_cluster_updatedsplit_source.png"), width=14, height=5.5, units = "in", res = 300)
print(p5)
dev.off()
  
png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_int_all_cluster_updated_split_exp_proc_tissue_location.png"), width=18, height=5.5, units = "in", res = 300)
print(p6)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_int_all_cluster_updated_split_exp_proc_tissue_location_y.png"), width=18, height=5.5, units = "in", res = 300)
print(p7)
dev.off()


# **Barplots ---------------------------------------------------------------
library(gridExtra)
clus_out = tum_obj
df_meta = data.frame(clus_out@meta.data, 
                     celltype = clus_out$clusters_final, 
                     stringsAsFactors = F) %>% 
  tibble() %>% clean_names()

df_meta1 = df_meta %>% mutate(patient_id =  factor(patient_id, levels = unique(patient_id[order(batch)]), ordered = TRUE))
# fn = factor(f, levels=unique(f[order(a,b,f)]), ordered=TRUE)

# create a summary table
p1=func_bar(var_df = df_meta1, var_x = "patient_id", var_fill = "batch", var_pos = "fill", var_cols = colors_dark, celltype = "source")
p2=func_bar(var_df = df_meta1, var_x = "patient_id", var_fill = "procedure", var_pos = "fill", var_cols = colors_dark[4:7], celltype = "procedure")

p3=func_bar(var_df = df_meta1, var_x = "sample_id", var_fill = "batch", var_pos = "fill", var_cols = colors_dark, celltype = "source")
p4=func_bar(var_df = df_meta1, var_x = "sample_id", var_fill = "procedure", var_pos = "fill", var_cols = colors_dark[4:7], celltype = "procedure")

p5=func_bar(var_df = df_meta1, var_x = "patient_id", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters_final")
p6=func_bar(var_df = df_meta1, var_x = "patient_id", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters_final")

p7=func_bar(var_df = df_meta1, var_x = "sample_id", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters_final")
p8=func_bar(var_df = df_meta1, var_x = "sample_id", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters_final")

p9=func_bar(var_df = df_meta1, var_x = "celltype", var_fill = "sample_id", var_pos = "stack", var_cols = colors_dark, celltype = "clusters_final")
p10=func_bar(var_df = df_meta1, var_x = "celltype", var_fill = "sample_id", var_pos = "fill", var_cols = colors_dark, celltype = "clusters_final")
p11=func_bar(var_df = df_meta1, var_x = "procedure", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters_final")
p12=func_bar(var_df = df_meta1, var_x = "procedure", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters_final")
p13=func_bar(var_df = df_meta1, var_x = "batch", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters_final")
p14=func_bar(var_df = df_meta1, var_x = "batch", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters_final")

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "barplots_all.png"), width=32, height=40, units = "in", res = 300)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, nrow = 4)
dev.off()

# difference between left and right ------
df_meta2 = df_meta1 %>% dplyr::filter(sample_id %in%  c("hbca_n17","hbca_n18","hbca_n19","hbca_n20","hbca_n21","hbca_n22","hbca_n23","hbca_n24"))
df_meta2 = df_meta2 %>% dplyr::mutate(breast_side = if_else(sample_id %in% c("hbca_n17", "hbca_n19", "hbca_n21", "hbca_n23"), "left", "right"))

p1=func_bar(var_df = df_meta2, var_x = "patient_id", var_fill = "breast_side", var_pos = "fill", var_cols = colors_dark, celltype = "side_source")
p2=func_bar(var_df = df_meta2, var_x = "patient_id", var_fill = "procedure", var_pos = "fill", var_cols = colors_dark[4:7], celltype = "side_procedure")

p3=func_bar(var_df = df_meta2, var_x = "sample_id", var_fill = "breast_side", var_pos = "fill", var_cols = colors_dark, celltype = "side_source")
p4=func_bar(var_df = df_meta2, var_x = "sample_id", var_fill = "procedure", var_pos = "fill", var_cols = colors_dark[4:7], celltype = "side_procedure")

p5=func_bar(var_df = df_meta2, var_x = "patient_id", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "side_clusters_final")
p6=func_bar(var_df = df_meta2, var_x = "patient_id", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "side_clusters_final")

p7=func_bar(var_df = df_meta2, var_x = "sample_id", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "side_clusters_final")
p8=func_bar(var_df = df_meta2, var_x = "sample_id", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "side_clusters_final")

p9=func_bar(var_df = df_meta2, var_x = "celltype", var_fill = "sample_id", var_pos = "stack", var_cols = colors_dark, celltype = "side_clusters_final")
p10=func_bar(var_df = df_meta2, var_x = "celltype", var_fill = "sample_id", var_pos = "fill", var_cols = colors_dark, celltype = "side_clusters_final")
p11=func_bar(var_df = df_meta2, var_x = "celltype", var_fill = "breast_side", var_pos = "stack", var_cols = colors_dark, celltype = "side_clusters_final")
p12=func_bar(var_df = df_meta2, var_x = "celltype", var_fill = "breast_side", var_pos = "fill", var_cols = colors_dark, celltype = "side_clusters_final")
p13=func_bar(var_df = df_meta2, var_x = "breast_side", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "side_clusters_final")
p14=func_bar(var_df = df_meta2, var_x = "breast_side", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "side_clusters_final")

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "barplots_all_breast_side_compare.png"), width=28, height=40, units = "in", res = 300)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, nrow = 4)
dev.off()
