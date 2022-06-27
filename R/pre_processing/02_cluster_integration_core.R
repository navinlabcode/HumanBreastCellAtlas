# ---------
# setwd -------

# geo ---------
setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/")
odir = "02_cluster_integration_core"; dir.create(odir)
library(future)
plan("multiprocess", workers = 3)
options(future.globals.maxSize = Inf)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/02_cluster_integration_core/"

# core ---------
setwd("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/")
odir = "02_cluster_integration_core"; dir.create(odir)
library(future)
plan("multiprocess", workers = 20)
options(future.globals.maxSize = Inf)
setwd(paste0("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/02_cluster_integration_core/"

# METHOD 9 ------------- 
###  logT and integrated by variable (var_split = "patient"), no scale by any vars
# ** NORM 
library(future)
plan("multiprocess", workers = 20)
#plan("sequential")

# ** Clus function --------
seu_object = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/01_read_files/tum_v1.rds")
seu_object = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/01_read_files/tum_v1.rds")

# run clustering function
dims = 20; k.param = 30; cols = colors_dark
fts1 = c("nFeature_RNA","nCount_RNA", "percent.mito","percent.rb", "percent.sg1", "percent.sg2", "percent.sg3")
fts2 = c("rna_EPCAM","rna_PTPRC","rna_PECAM1","rna_DCN","rna_KRT5",  "rna_LTF", "rna_KRT17", "rna_HPGD", "rna_TPM2")
fts3 = c("rna_CD69","rna_CD3D","rna_CD68","rna_CD8A","rna_CD4",  "rna_NKG7", "rna_CD19", "rna_MS4A1", "rna_RGS5")
fts4 = c("rna_VWF","rna_PRSS1","rna_CCL21","rna_FAP","rna_COL1A1",  "rna_COL6A1", "rna_ESR2", "rna_KRT8", "rna_SERHL2")
clus_out = func_cluster_vst(tum_int = seu_object, var_scale = "NULL", var_split = "patient_id", dims = dims, k.param = k.param, cols = colors_dark,  
                            markers1 = fts1, markers2 = fts2, markers3 = fts3, markers4 = fts4, save_path = save_path)
clus_out = func_cluster_rpca(tum_int = seu_object, var_scale = "NULL", var_split = "exp_proc", dims = dims, k.param = k.param, cols = colors_dark,  
                            markers1 = fts1, markers2 = fts2, markers3 = fts3, markers4 = fts4, save_path = save_path)
char_ident = "integrated_snn_res.0.2"

# CHECK CLUSTERING AND SAVE
#tum.integrated.ns = read_rds(glue("{save_path}/dim_{dims}k_{k.param}_tum.integrated.ns2_method_9.rds"))
write_rds(clus_out, glue("{save_path}/dim_{dims}k_{k.param}_tum.integrated.ns2_method_9rcp.rds"))


# ** Markers function --------
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20200210_final1/02_cluster_integration/"
char_ident = "integrated_snn_res.0.2"  #integrated_snn_res.0.1
clus_out = read_rds(path = glue("{save_path}/dim_{dims}k_{k.param}_tum.integrated.ns2_method_9rcp.rds")) 
mark_out = func_markers_vst(clus_out = clus_out, char_ident = char_ident, fts1 = fts1, dims = dims, k.param = k.param, cols = colors_dark, 
                            save_path = save_path)

# METHOD 13 ------------- 
###  SCTransform, scale by a variable, integration

library(future)
plan("multiprocess", workers = 10)
#plan("sequential")

# ** Clus function --------
seu_object = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20200210_final1/01_read_files_core/tum_filt_v1.rds")

# run clustering function
dims = 20; k.param = 30; cols = colors_dark; var_split = "patient_id"; var_scale = "nCount_RNA"
fts1 = c("nFeature_RNA","nCount_RNA", "percent.mito","percent.rb", "percent.sg1", "percent.sg2", "percent.sg3")
fts2 = c("rna_EPCAM","rna_PTPRC","rna_PECAM1","rna_DCN","rna_KRT5",  "rna_LTF", "rna_KRT17", "rna_HPGD")
fts3 = c("rna_CD69","rna_CD3D","rna_CD68","rna_CD8A","rna_CD4",  "rna_NKG7", "rna_CD19", "rna_MS4A1")
fts4 = c("rna_VWF","rna_PRSS1","rna_CCL21","rna_FAP","rna_COL1A1",  "rna_COL6A1", "rna_ESR2", "rna_KRT8")
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20200210_final1/02_cluster_integration_core/"

clus_out = func_cluster_scale(tum_int = seu_object, var_split = var_split, dims = dims, k.param = k.param, 
                              markers1 = fts1, markers2 = fts2, markers3 = fts3, markers4 = fts4, save_path = save_path,
                              cols = colors_dark,
                              var_scale = var_scale)

char_ident = "integrated_snn_res.0.4"


# CHECK CLUSTERING AND SAVE
#tum.integrated.ns = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/ecis/analysis/20200120_ecis_tillsamp20/03_read_files_v2/all/dim_30k_20_tum.integrated.ns_SEL_v4vst_method_9.rds")
write_rds(tum.integrated.ns, glue("{save_path}/dim_{dims}k_{k.param}_split_by_{var_split}_scaled_by_{var_scale}_tum.integrated.ns3_method_13.rds")) 

# ** Markers function --------
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20200210_final1/02_cluster_integration/"
char_ident = "integrated_snn_res.0.4"  #integrated_snn_res.0.1
clus_out = read_rds(path = glue("{save_path}/dim_{dims}k_{k.param}_tum.integrated.ns_SCT_v4_intgeratedassay_{char_ident}_vst_method9.rds")) 
mark_out = func_markers_vst(clus_out = clus_out, char_ident = char_ident, fts1 = fts1, dims = dims, k.param = k.param, cols = colors_dark, 
                            save_path = save_path)

# extra
tum = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20200210_final1/02_cluster_integration/dim_30k_200_tum.integrated.ns3_method_9.rds")
length(tum)
DefaultAssay(tum) = "RNA"
FeaturePlot(tum, features = c("KRT5", "KRT14", "SLPI", "ANKRD30A"), reduction = "umap", pt.size = 0.001, max.cutoff = "q90")
DimPlot(tum_basal, group.by = "integrated_snn_res.0.6", reduction = "umap", pt.size = 0.01, cols = colors_dark)
DimPlot(tum, group.by = "source", reduction = "umap", pt.size = 0.01, cols = colors_dark)

tum_basal = subset(tum, subset = integrated_snn_res.0.6 == "5")
table(tum_basal$source)
table(tum_basal$exp_proc, tum_basal$source)
table(tum_basal$exp_proc, tum_basal$sample_id )


tum1 = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20200210_final1/02_cluster_integrationv2/dim_30k_200_tum.integrated.ns2_method_9.rds")
length(tum1$orig.ident)
DefaultAssay(tum1) = "RNA"
FeaturePlot(tum1, features = c("KRT5", "KRT14", "SLPI", "ANKRD30A"), reduction = "umap", pt.size = 0.001, max.cutoff = "q90")
DimPlot(tum_basal1, group.by = "integrated_snn_res.0.6", reduction = "umap", pt.size = 0.01, cols = colors_dark)
DimPlot(tum1, group.by = "source", reduction = "umap", pt.size = 0.01, cols = colors_dark)

tum_basal1 = subset(tum, subset = integrated_snn_res.0.6 == "1")
table(tum_basal1$source)
table(tum_basal1$exp_proc, tum_basal1$source)
table(tum_basal1$exp_proc, tum_basal1$sample_id )


# make umaps -------

devtools::install_github("JinmiaoChenLab/Rphenograph")
tum.integrated.ns = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20200210_final1/02_cluster_integrationv2/dim_30k_125_tum.integrated.ns2_method_9.rds")
tum.integrated.ns = tum
dims = 30;reduction = "umap"
tum.integrated.ns = RunUMAP(object = tum.integrated.ns, reduction = "pca", dims = 1:20, umap.method = "umap-learn", metric = "correlation")
DimPlot(object = tum.integrated.ns, reduction = "umap", label = TRUE, pt.size = 0.001, cols = colors_dark, group.by = "integrated_snn_res.0.6")
Idents(tum.integrated.ns) = "integrated_snn_res.0.6"
VlnPlot(tum.integrated.ns, features = c("percent.sg1", "percent.sg2", "percent.sg3", "percent.rb"), ncol = 2, pt.size = 0, cols = colors_dark)


