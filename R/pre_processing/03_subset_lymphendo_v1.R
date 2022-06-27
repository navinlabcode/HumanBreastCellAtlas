# geo ----
setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_lymphendo_v1";dir.create(odir)
options(future.globals.maxSize = 90000 * 1024^2)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "figure_1_v1"
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_lymphendo_v1/"
stress_sig = read.delim("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_fibroblasts_final1/stromal_stress_signature_032720.txt")
stress_sig = as.character(stress_sig$stromal.stress.sig)

# core ----
setwd("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_lymphendo_v1"; dir.create(odir)
library(future)
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 50000 * 1024^2)
setwd(paste0("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "figure_1_v1"
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lymphendo_v1/"
stress_sig = read.delim("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_fibroblasts_final1/stromal_stress_signature_032720.txt")
stress_sig = as.character(stress_sig$stromal.stress.sig)


# read original object --------------
tum = readRDS(file = paste0(wd, data_folder,"/final_use.rds"))

# subset ---------
Idents(tum) <- "final_group1"
tum_sub = subset(tum, idents = "Lymphatic") #endo_vas
tum_sub_cells = WhichCells(object = tum, idents = "Lymphatic") #subset epi_lumhr cells #83891 cells

#fibro = subset(tum, idents = "Pericytes") #endo_vas #Pericytes

# save objects
percent.stress <- Matrix::colSums(GetAssayData(object = tum_sub)[stress_sig, ])/Matrix::colSums(GetAssayData(object = tum_sub))
tum_sub <- AddMetaData(tum_sub, metadata = percent.stress, col.name = "percent.stress")
write_rds(tum_sub, path = glue("{save_path}/endo_lymph_subset_v1.rds"))
write.csv(tum_sub_cells, file = glue("{save_path}/endo_lymph_cells_v1.csv"))

png(paste0(save_path, "vlnplot_percent.stress_sample_id.png"), width=14, height=5, units = "in", res = 300)
VlnPlot(tum_sub, features = "percent.stress", group.by = "final_sample_id", pt.size = 0 )
dev.off()


# method 9 ----------
# Log T integration with no scaling
library(future)
plan("multiprocess", workers = 25)
#plan("sequential")

# ** Clus function --------
# run clustering function
seu_object = read_rds(glue("{save_path}/endo_lymph_subset_v1.rds"))
DefaultAssay(seu_object) <- "RNA"
dims = 12; k.param = 20; cols = colors_dark; var_split = "batch"; var_scale = "nCount_RNA"
fts1 = c("nFeature_RNA","nCount_RNA", "percent.mito","percent.rb", "percent.stress")
fts2 = c("rna_EPCAM","rna_PTPRC","rna_PECAM1","rna_DCN","rna_KRT5",  "rna_NCAM1", "rna_VWF", "rna_CLU")
fts3 = c("rna_CDH5", "rna_GJA5", "rna_BMX", "rna_ACKR1", "rna_CA4", "rna_PROX1", "rna_PDPN", "rna_CCL21")
fts4 = c("rna_MCAM","rna_TFF3","rna_MMRN1","rna_LYVE1","rna_CA4","rna_BTNL9",  "rna_FABP4", "rna_TIMP4", "rna_LPL")

#ref_set = c("pt01", "pt06","pt09", "pt20", "pt21","pt25", "pt27", "pt28", "pt31")
ref_set = NULL

# seu_object$regroup = seu_object$patient_id
# seu_object$regroup[seu_object$patient_id == "pt14"] = "pt_cmb"
# seu_object$regroup[seu_object$patient_id == "pt15"] = "pt_cmb"
# seu_object$regroup[seu_object$patient_id == "pt16"] = "pt_cmb"
# seu_object$regroup[seu_object$patient_id == "pt18"] = "pt_cmb"
# seu_object$regroup[seu_object$patient_id == "pt19"] = "pt_cmb"
# seu_object$regroup[seu_object$patient_id == "pt28"] = "pt_cmb"
# seu_object$regroup[seu_object$patient_id == "pt33"] = "pt_cmb"


clus_out = func_cluster_vst(tum_int = seu_object, var_scale = var_scale, var_split = var_split, dims = dims, k.param = k.param, 
                            markers1 = fts1, markers2 = fts2, markers3 = fts3, markers4 = fts4, save_path = save_path,
                            cols = colors_dark, ref_set = ref_set)

char_ident = "integrated_snn_res.0.1"


# ** Markers function (9)--------
dims = 12; k.param = 20;
char_ident = "integrated_snn_res.0.2"  #integrated_snn_res.0.1
ht_fts1 = c("PECAM1")
#clus_out = read_rds(path = glue("{save_path}/dim_{dims}k_{k.param}_tum.integrated.ns_SCT_v4_intgeratedassay_{char_ident}_vst_method9.rds")) 
mark_out = func_markers_vst(clus_out = clus_out, char_ident = char_ident, fts1 = ht_fts1, dims = dims, k.param = k.param, cols = colors_dark, 
                            save_path = save_path)

clus_out = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_lymphendo_v1/dim_12k_20_tum.integrated.ns_RNA_v4_intgeratedassay_integrated_snn_res.0.2_vst_method_9.rds")
markers = read.csv("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_lymphendo_v1/dims_12k_20_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.2_vst_method_9.csv")
# manual checks ------
tum = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_lymphendo_v1/dim_12k_20_nCount_RNA_tum.integrated.ns2_method_9.rds")

# filter other cells --------
Idents(tum) <- "integrated_snn_res.0.2"
tum_sub = subset
DimPlot(tum, group.by = "integrated_snn_res.0.2", label = T)
reduction = "umap"
p1=FeaturePlot(tum, reduction = reduction, features = c("rna_KRT14", "rna_TAGLN", "rna_ACTA2", "rna_MYLK", "rna_KRT5", "rna_KIT"))
p2=FeaturePlot(tum, reduction = reduction, features = c("rna_ACTA2", "rna_PDGFRB", "rna_FAP", "rna_TGFB1", "rna_PDGFA","rna_PDGFB", "rna_FGF2", 
                                                                   "rna_S100A8", "rna_VIM", "rna_DDR2", "rna_HGF", "rna_TIMP1", "rna_ADAM10"))
p3=FeaturePlot(tum, reduction = reduction, features = c("rna_ACTA2", "rna_MYH11", "rna_TAGLN", "rna_HHIP", "rna_ASPN",  "rna_JUNB", "rna_SCX"))
p4=FeaturePlot(tum, reduction = reduction, features = c("rna_COL1A1", "rna_CA4", "rna_CD36", "rna_ACKR1", "rna_CREM"))

# add module score
p5=VlnPlot(tum, features = c("nCount_RNA", "nFeature_RNA"), group.by = "integrated_snn_res.0.1", pt.size = 0)

# system("pdftoppm input.pdf outputname -png")
# save plots
png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale, "umap_mk1_9a.png"), width=12, height=12, units = "in", res = 300)
print(p1)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale, "umap_mk2_9a.png"), width=12, height=12, units = "in", res = 300)
print(p2)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_mk3_9a.png"), width=12, height=12, units = "in", res = 300)
print(p3)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_mk4_9a.png"), width=12, height=12, units = "in", res = 300)
print(p4)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_mk5violin_9a.png"), width=12, height=5, units = "in", res = 300)
print(p5)
dev.off()

