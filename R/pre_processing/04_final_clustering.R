
# ---------
# setwd -------

# geo ---------
setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/")
odir = "04_final_clustering"; dir.create(odir)
library(future)
plan("multiprocess", workers = 30)
options(future.globals.maxSize = Inf)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/04_final_clustering/"

#cl = c("#c385ec", "#d98f00", "#00bfc1", "#f47c55", "#ff8ca5", "#405f8f")

# core ---------
setwd("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/")
odir = "04_final_clustering"; dir.create(odir)
library(future)
plan("multiprocess", workers = 50)
#plan("multisession", workers = 20)
options(future.globals.maxSize = Inf)
setwd(paste0("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/04_final_clustering/"

# read object -----
tum = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_merge_objects/dim_20k_30final_tum.integrated.ns_RNA_v4_intgeratedassay_integrated_snn_res.0.6_vst_method_9.rds")
tum = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_merge_objects/dim_20k_30final_tum.integrated.ns_RNA_v4_intgeratedassay_integrated_snn_res.0.6_vst_method_9.rds")

# make plots -----

p0 = DimPlot(object = tum, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = colors_dark, group.by = "integrated_snn_res.0.6", shuffle = T)+  NoAxes()
p00 = DimPlot(object = tum, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = colors_dark, group.by = "cluster_updated", shuffle = T)+  NoAxes()
p1 = VlnPlot(tum, group.by = "integrated_snn_res.0.6", features = c("nCount_RNA", "nFeature_RNA"), pt.size = 0, cols = colors_dark)
p = cowplot::plot_grid(p0, p00, p1, nrow = 2)
png(paste0(save_path,"uamp_integrated_dims20_integrated_snn_res.0.6_dist0.1_QC.png"), width=18, height=15, units = "in", res = 300)
print(p)
dev.off()

# remove clusters -----
Idents(tum) <- "integrated_snn_res.0.6"

tum1 = subset(tum, idents = c("24", "25", "26", "27", "30", "31"), invert = TRUE)
write_rds(tum1, "tum_final.rds")

p00 = DimPlot(object = tum1, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = colors_dark, group.by = "cluster_updated", shuffle = T)+  NoAxes()
png(paste0(save_path,"uamp_integrated_dims20_integrated_snn_res.0.6_dist0.1_filtered.png"), width=8, height=5, units = "in", res = 300)
print(p00)
dev.off()


# recluster -------
tum = read_rds("tum_final.rds")
tum.integrated.ns = tum
dims=22; k.param=30; var_scale = "NULL"
DefaultAssay(tum.integrated.ns) <- "integrated"
tum.integrated.ns = ScaleData(tum.integrated.ns, assay = "integrated", vars.to.regress = NULL, features = rownames(tum.integrated.ns@assays$integrated@data)) # no scaling # no scaling
tum.integrated.ns = RunPCA(tum.integrated.ns, npcs = 100)
write_rds(tum.integrated.ns, glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_tum.integrated.ns11_method_9rcp.rds")) #tum_int = read_rds(glue("{odir}/dim_{dims}k_{k.param}_tum.integrated.ns11_method_9.rds"))
e1 = ElbowPlot(tum.integrated.ns, ndims = 30)
png(paste0(save_path,"elbow", dims,"k_", k.param, "elbow100_method_9rcp.png"), width=8, height=5, units = "in", res = 300)
print(e1)
dev.off()
tum.integrated.ns = RunUMAP(object = tum.integrated.ns, reduction = "pca", dims = 1:dims, min.dist = 0.1, spread = 1.5) #, umap.method = 'umap-learn', metric = 'manhattan'
tum.integrated.ns = FindNeighbors(object = tum.integrated.ns, dims = 1:dims, reduction = "pca", k.param = k.param)
tum.integrated.ns = FindClusters(object = tum.integrated.ns, resolution = c(1, 0.8, 0.6,0.4,0.3,0.2, 0.1), verbose = TRUE)
write_rds(tum.integrated.ns, glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_tum.integrated.ns2_method_9rcp_0.1_1.5.rds")) #tum_int = read_rds(glue("{odir}/dim_{dims}k_{k.param}_tum.integrated.ns_method_9.rds"))

#tum.integrated.ns = read_rds(glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_tum.integrated.ns2_method_9rcp_0.1_1.5.rds"))

p00 = DimPlot(object = tum.integrated.ns, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = colors_dark, group.by = "cluster_updated", shuffle = T)+  NoAxes()
p01 = DimPlot(object = tum.integrated.ns, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = colors_dark, group.by = "integrated_snn_res.0.6", shuffle = T)+  NoAxes()
png(paste0(save_path, "dim_",dims, "k_",k.param,"_",var_scale, "_uamp_integrated_snn_res.0.6_dist0.1_filtered2.png"), width=14, height=5, units = "in", res = 300)
print(plot_grid(p00, p01))
dev.off()


source("densmap.R")
out <- densMAP(as.matrix(tum.integrated.ns@assays$integrated@counts))

# meta data correction -----
meta_uniq = read.csv("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_merge_objects/meta_uniq.csv")
tum = tum.integrated.ns
meta_df = tum@meta.data
meta_df$cellname = rownames(meta_df)
write.csv(meta_df, "meta_df.csv");write_rds(meta_df, "meta_df.rds")
dim(meta_df)
colnames(meta_df)

meta_df1 = left_join(meta_df, meta_uniq, by = c("sample_id", "batch"))
dim(meta_df);dim(meta_df1)
rownames(meta_df1) = meta_df$cellname

# add back meta data
tum@meta.data = meta_df1
dim(tum@meta.data)
colnames(tum@meta.data)
tum[['patient_id.x']] <- NULL
tum[['patient_id.y']] <- NULL
tum[['ts_sampleid.x']] <- NULL
tum[['ts_ptid.x']] <- NULL
tum[['tissue_location.x']] <- NULL
tum[['brca_updated.x']] <- NULL
tum[['cancer_risk.x']] <- NULL
tum[['cancer_history.x']] <- NULL
tum[['age_updated.x']] <- NULL
tum[['gravida.1.x']] <- NULL
tum[['parity.1.x']] <- NULL
tum[['menopause.x']] <- NULL

tum[['fibrosis.x']] <- NULL
tum[['hyperplasia.x']] <- NULL
tum[['metaplasia.x']] <- NULL
tum[['density.1.x']] <- NULL
tum[['figure1_grp.x']] <- NULL

tum[['ts_sampleid.y']] <- NULL
tum[['ts_ptid.y']] <- NULL
tum[['tissue_location.y']] <- NULL
tum[['brca_updated.y']] <- NULL
tum[['cancer_risk.y']] <- NULL
tum[['cancer_history.y']] <- NULL
tum[['age_updated.y']] <- NULL
tum[['gravida.1.y']] <- NULL
tum[['parity.1.y']] <- NULL
tum[['menopause.y']] <- NULL

tum[['fibrosis.y']] <- NULL
tum[['hyperplasia.y']] <- NULL
tum[['metaplasia.y']] <- NULL
tum[['density.1.y']] <- NULL
tum[['figure1_grp.y']] <- NULL
colnames(tum@meta.data)

tracker = read_sheet("~/projects/bcap/hbca_sample_tracker.xlsx", sheet = 1)
tum@meta.data = left_join(tum@meta.data, tracker, by = c("final_sample_id" = "sheet_samp_id"))
View(tum@meta.data)

write_rds(tum, "final1.rds")
tum = read_rds("final1.rds")

# rename clusters ---
tum.integrated.ns = tum
tum.integrated.ns$clusters3 = as.character(tum.integrated.ns$integrated_snn_res.0.6)
tum.integrated.ns$clusters3[tum.integrated.ns$integrated_snn_res.0.6 %in% c(0, 3, 8, 10, 11, 13, 20)] = "fibro"
tum.integrated.ns$clusters3[tum.integrated.ns$integrated_snn_res.0.6 %in% c(1)] = "t"
tum.integrated.ns$clusters3[tum.integrated.ns$integrated_snn_res.0.6 %in% c(2, 25,23)] = "epi_lumhr" 
tum.integrated.ns$clusters3[tum.integrated.ns$integrated_snn_res.0.6 %in% c(7,9)] = "epi_basal"
tum.integrated.ns$clusters3[tum.integrated.ns$integrated_snn_res.0.6 %in% c(5,15)] = "epi_lumsec"
tum.integrated.ns$clusters3[tum.integrated.ns$integrated_snn_res.0.6 %in% c(4,12,18)] = "endo_vas"
tum.integrated.ns$clusters3[tum.integrated.ns$integrated_snn_res.0.6 %in% c(17)] = "endo_lymph"
tum.integrated.ns$clusters3[tum.integrated.ns$integrated_snn_res.0.6 %in% c(6)] = "pericytes"
tum.integrated.ns$clusters3[tum.integrated.ns$integrated_snn_res.0.6 %in% c(19,24)] = "b"
tum.integrated.ns$clusters3[tum.integrated.ns$integrated_snn_res.0.6 %in% c(14,16,21)] = "macro"
tum.integrated.ns$clusters3[tum.integrated.ns$integrated_snn_res.0.6 %in% c(22)] = "nk_t"
tum.integrated.ns$clusters3[tum.integrated.ns$integrated_snn_res.0.6 %in% c(26,27,28)] = "remove"
Idents(tum.integrated.ns) <- "clusters3"
tum1 = subset(tum.integrated.ns, idents = c("remove"), invert = TRUE)
#Idents(tum.integrated.ns) <- "integrated_snn_res.0.6"
#tum1 = subset(tum1, idents = c("26", "27", "28"), invert = TRUE)
write_rds(tum1, "final2.rds")


# recluster -------
tum.integrated.ns = read_rds("final2.rds")
dims=15; k.param=30; var_scale = "NULL"
DefaultAssay(tum.integrated.ns) <- "integrated"
tum.integrated.ns = ScaleData(tum.integrated.ns, assay = "integrated", vars.to.regress = "nCount_RNA", features = rownames(tum.integrated.ns@assays$integrated@data)) # no scaling # no scaling
tum.integrated.ns = RunPCA(tum.integrated.ns, npcs = 100)
tum.integrated.ns = RunUMAP(object = tum.integrated.ns, reduction = "pca", dims = 1:dims, min.dist = 0.1, spread = 1.5) #, umap.method = 'umap-learn', metric = 'manhattan'
tum.integrated.ns = FindNeighbors(object = tum.integrated.ns, dims = 1:dims, reduction = "pca", k.param = k.param)
tum.integrated.ns = FindClusters(object = tum.integrated.ns, resolution = c(1, 0.8, 0.6,0.4,0.3,0.2, 0.1), verbose = TRUE)
write_rds(tum.integrated.ns, glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_tum.integrated.ns2_method_9rcp_0.1_1.5final.rds")) #tum_int = read_rds(glue("{odir}/dim_{dims}k_{k.param}_tum.integrated.ns_method_9.rds"))

#tum.integrated.ns = read_rds(glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_tum.integrated.ns2_method_9rcp_0.1_1.5.rds"))

library(cowplot)
p00 = DimPlot(object = tum.integrated.ns, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = colors_dark, group.by = "clusters3", shuffle = T)+  NoAxes()
p01 = DimPlot(object = tum.integrated.ns, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = colors_dark, group.by = "integrated_snn_res.0.6", shuffle = T)+  NoAxes()
png(glue("{save_path}/uamp_integrated_dim_{dims}k_{k.param}_{var_scale}_integrated_snn_res.0.6_dist0.1_final.png"), width=14, height=5, units = "in", res = 300)
print(plot_grid(p00, p01))
dev.off()

write_rds(tum.integrated.ns, glue("{save_path}/final3.rds")) #tum_int = read_rds(glue("{odir}/dim_{dims}k_{k.param}_tum.integrated.ns_method_9.rds"))




# run markers -----
tum = read_rds("final2.rds")
fts1 = c("EPCAM", "PTPRC", "PECAM1")
char_ident = "integrated_snn_res.0.6"  #char_ident = "clusters1"
mark_out = func_markers_vst(clus_out = tum.integrated.ns, char_ident = char_ident, fts1 = fts1, dims = dims, k.param = k.param, cols = colors_dark, 
                            save_path = save_path)
char_ident = "clusters3"  #char_ident = "clusters1"
mark_out = func_markers_vst(clus_out = tum.integrated.ns, char_ident = char_ident, fts1 = fts1, dims = dims, k.param = k.param, cols = colors_dark, 
                            save_path = save_path)



