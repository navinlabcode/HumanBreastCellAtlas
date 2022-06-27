# geo ----
setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_B_v1";dir.create(odir)
options(future.globals.maxSize = 90000 * 1024^2)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "figure_1_v1"
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_B_v1/"
# stress_sig = read.delim("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_fibroblasts_final1/stromal_stress_signature_032720.txt")
# stress_sig = as.character(stress_sig$stromal.stress.sig)

# core ----
setwd("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_B_v1"; dir.create(odir)
library(future)
plan("multiprocess", workers = 50)
options(future.globals.maxSize = Inf)
setwd(paste0("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "figure_1_v1"
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_B_v1/"
# stress_sig = read.delim("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis//03_subset_fibroblasts_final1/stromal_stress_signature_032720.txt")
# stress_sig = as.character(stress_sig$stromal.stress.sig)


# read original object --------------
tum = readRDS(file = paste0(wd, data_folder,"/final_use.rds"))

# subset ---------
Idents(tum) <- "final_group1"
tum_sub = subset(tum, idents = "B-cells") #endo_vas
tum_sub_cells = WhichCells(object = tum, idents = "B-cells") #subset epi_lumhr cells #83891 cells
write_rds(tum_sub, path = glue("{save_path}/b_subset_v1.rds"))

# save objects
percent.stress <- Matrix::colSums(GetAssayData(object = tum_sub)[stress_sig, ])/Matrix::colSums(GetAssayData(object = tum_sub))
tum_sub <- AddMetaData(tum_sub, metadata = percent.stress, col.name = "percent.stress")
write_rds(tum_sub, path = glue("{save_path}/b_subset_v1.rds"))
write.csv(tum_sub_cells, file = glue("{save_path}/b_cells_v1.csv"))

png(paste0(save_path, "vlnplot_percent.stress_sample_id.png"), width=14, height=5, units = "in", res = 300)
VlnPlot(tum_sub, features = "percent.stress", group.by = "final_sample_id", pt.size = 0 )
dev.off()




# method 9 ----------
# Log T integration with no scaling
library(future)
plan("multiprocess", workers = 50)
#plan("sequential")

# ** Clus function --------
# run clustering function
seu_object = read_rds(glue("{save_path}/b_subset_v1.rds"))
DefaultAssay(seu_object) <- "RNA"
dims = 12; k.param = 20; cols = colors_dark; var_split = "sheet_exp_proc"; var_scale = "nCount_RNA"
fts1 = c("nFeature_RNA","nCount_RNA", "percent.mito","percent.rb", "percent.stress")
fts2 = c("rna_EPCAM","rna_PTPRC","rna_PECAM1","rna_DCN","rna_KRT5",  "rna_CD86", "rna_COL1A1", "rna_CLU")
fts3 = c("rna_PTPRC","rna_CD79A","rna_CD19","rna_IGKC","rna_IGLC3","rna_IGHG3", "rna_CD79B", "rna_IGHM")
fts4 = c("rna_IRF8","rna_MS4A1","rna_IGHD","rna_CD27","rna_MS4A1", "rna_IGHG1","rna_JCHAIN", "rna_IGHA1","rna_IGHM", "rna_IGLC2", "rna_IGLC3")

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

char_ident = "integrated_snn_res.0.4"


# ** Markers function (9)--------
dims = 12; k.param = 20;
char_ident = "integrated_snn_res.0.4"  #integrated_snn_res.0.1
ht_fts1 = c("MS4A1", "JCHAIN", "IGHA1")
#clus_out = read_rds(path = glue("{save_path}/dim_{dims}k_{k.param}_tum.integrated.ns_SCT_v4_intgeratedassay_{char_ident}_vst_method9.rds")) 
mark_out = func_markers_vst(clus_out = clus_out, char_ident = char_ident, fts1 = ht_fts1, dims = dims, k.param = k.param, cols = colors_dark, 
                            save_path = save_path)

clus_out = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_B_v1/dim_12k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")
#markers = read.csv("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_lymphendo_v1/dims_12k_20_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.2_vst_method_9.csv")



# manual checks ------
#tum = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_lymphendo_v1/dim_12k_20_nCount_RNA_tum.integrated.ns2_method_9.rds")

Idents(tum) <- "integrated_snn_res.0.2"
DimPlot(clus_out, group.by = "integrated_snn_res.0.4", label = T, split.by = "sheet_exp_proc")
reduction = "umap"

VlnPlot(tum, features = c("rna_MS4A1","rna_CD79A", "rna_PTPRC", "rna_CD3D", "rna_EPCAM", "rna_CD68", "rna_KRT8", "rna_COL1A1", "rna_RGS5", "rna_KRT5", "rna_NKG7"), 
        pt.size = 0, group.by = "integrated_snn_res.0.4")


b = c("rna_MS4A1", "rna_BANK1", "rna_PAX5", "rna_CD24", "rna_CD19", "rna_CR2", "rna_TNFRSF13C")
plasma = c("rna_IGLC2", "rna_IGHA1", "rna_IGHG3", "rna_MZB1", "rna_SDC1", "rna_TNFRSF17", "rna_LEU1", "rna_PRDM1", "rna_CD1D")
prol = c("rna_TOP2A", "rna_MKI67")
naive = c("rna_IGHD", "rna_MS4A1", "rna_NT5E", "rna_RGS13", "rna_CD24", "rna_NEIL1") #NT5E is cd73
memory = c("rna_CD27","rna_IGHD")

p1=FeaturePlot(tum, reduction = reduction, features = c(b, plasma, prol, naive, memory), order = TRUE, pt.size = 0.001)
p2=FeaturePlot(tum, reduction = reduction, features = c("rna_ACTA2", "rna_PDGFRB", "rna_FAP", "rna_TGFB1", "rna_PDGFA","rna_PDGFB", "rna_FGF2", 
                                                        "rna_S100A8", "rna_VIM", "rna_DDR2", "rna_HGF", "rna_TIMP1", "rna_ADAM10"))
p3=FeaturePlot(tum, reduction = reduction, features = c("rna_ACTA2", "rna_MYH11", "rna_TAGLN", "rna_HHIP", "rna_ASPN",  "rna_JUNB", "rna_SCX"))
p4=FeaturePlot(tum, reduction = reduction, features = c("rna_COL1A1", "rna_CA4", "rna_CD36", "rna_ACKR1", "rna_CREM"))

# add module score
p5=VlnPlot(tum, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), group.by = "integrated_snn_res.0.4", pt.size = 0)

# system("pdftoppm input.pdf outputname -png")
# save plots
png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale, "umap_mk1_9a.png"), width=12, height=16, units = "in", res = 300)
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


# REMOVE CLUSTERS --------
# _ ------------------------- --------------------------------
tum = clus_out
Idents(tum) <- "integrated_snn_res.0.4"
tum_sub = subset(tum, idents = c("9", "10", "12", "7", "8", "1"), invert = TRUE)
DimPlot(tum_sub, group.by = "integrated_snn_res.0.4", label = T)

write_rds(tum_sub, "tum_sub.rds")

# redo the clustering ---



# method 9 ----------
# Log T integration with no scaling
library(future)
plan("multiprocess", workers = 50)
#plan("sequential")

# ** Clus function --------
# run clustering function
seu_object = read_rds(glue("{save_path}/tum_sub.rds"))
DefaultAssay(seu_object) <- "RNA"
dims = 20; k.param = 20; cols = colors_dark; var_split = "sheet_exp_proc"; var_scale = "nCount_RNA"
fts1 = c("nFeature_RNA","nCount_RNA", "percent.mito","percent.rb", "percent.stress")
fts2 = c("rna_EPCAM","rna_PTPRC","rna_PECAM1","rna_DCN","rna_KRT5",  "rna_CD86", "rna_COL1A1", "rna_CLU")
fts3 = c("rna_PTPRC","rna_CD79A","rna_CD19","rna_IGKC","rna_IGLC3","rna_IGHG3", "rna_CD79B", "rna_IGHM")
fts4 = c("rna_IRF8","rna_MS4A1","rna_IGHD","rna_CD27","rna_MS4A1", "rna_IGHG1","rna_JCHAIN", "rna_IGHA1","rna_IGHM", "rna_IGLC2", "rna_IGLC3")

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

char_ident = "integrated_snn_res.0.4"


# ** Markers function (9)--------
dims = 20; k.param = 20;
char_ident = "integrated_snn_res.0.3"  #integrated_snn_res.0.1
ht_fts1 = c("MS4A1", "JCHAIN", "IGHA1")
#clus_out = read_rds(path = glue("{save_path}/dim_{dims}k_{k.param}_tum.integrated.ns_SCT_v4_intgeratedassay_{char_ident}_vst_method9.rds")) 
mark_out = func_markers_vst(clus_out = clus_out, char_ident = char_ident, fts1 = ht_fts1, dims = dims, k.param = k.param, cols = colors_dark, 
                            save_path = save_path)

clus_out = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_B_v1/dim_20k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")
clus_out = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_B_v1/dim_20k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")
#markers = read.csv("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_lymphendo_v1/dims_12k_20_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.2_vst_method_9.csv")



# manual checks ------
#tum = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_lymphendo_v1/dim_12k_20_nCount_RNA_tum.integrated.ns2_method_9.rds")
tum = clus_out
Idents(tum) <- "integrated_snn_res.0.3"
DimPlot(clus_out, group.by = "integrated_snn_res.0.3", label = T)
reduction = "umap"

VlnPlot(clus_out, features = c("rna_MS4A1","rna_CD79A", "rna_PTPRC", "rna_CD27", "rna_EPCAM", "rna_CD68", "rna_KRT8", "rna_COL1A1", "rna_RGS5", "rna_KRT5", "rna_NKG7"), 
        pt.size = 0, group.by = "integrated_snn_res.0.3")

VlnPlot(clus_out, features = c("rna_IGK","rna_IGM", "rna_PTPRC", "rna_CD3D", "rna_EPCAM", "rna_CD68", "rna_KRT8", "rna_COL1A1", "rna_RGS5", "rna_KRT5", "rna_NKG7"), 
        pt.size = 0, group.by = "integrated_snn_res.0.4")


b = c("rna_MS4A1", "rna_BANK1", "rna_PAX5", "rna_CD24", "rna_CD19", "rna_CR2", "rna_TNFRSF13C")
plasma = c("rna_IGLC2", "rna_IGHA1", "rna_IGHG3", "rna_MZB1", "rna_SDC1", "rna_TNFRSF17", "rna_LEU1", "rna_PRDM1", "rna_CD1D")
prol = c("rna_TOP2A", "rna_MKI67")
naive = c("rna_IGHD", "rna_MS4A1", "rna_NT5E", "rna_RGS13", "rna_CD24", "rna_NEIL1") #NT5E is cd73
memory = c("rna_CD27","rna_IGHD")

p1=FeaturePlot(tum, reduction = reduction, features = c(b, plasma, prol, naive, memory), order = TRUE, pt.size = 0.001)
p2=FeaturePlot(tum, reduction = reduction, features = c("rna_ACTA2", "rna_PDGFRB", "rna_FAP", "rna_TGFB1", "rna_PDGFA","rna_PDGFB", "rna_FGF2", 
                                                        "rna_S100A8", "rna_VIM", "rna_DDR2", "rna_HGF", "rna_TIMP1", "rna_ADAM10"))
p3=FeaturePlot(tum, reduction = reduction, features = c("rna_ACTA2", "rna_MYH11", "rna_TAGLN", "rna_HHIP", "rna_ASPN",  "rna_JUNB", "rna_SCX"))
p4=FeaturePlot(tum, reduction = reduction, features = c("rna_COL1A1", "rna_CA4", "rna_CD36", "rna_ACKR1", "rna_CREM"))

# add module score
p5=VlnPlot(tum, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), group.by = "integrated_snn_res.0.4", pt.size = 0)

# system("pdftoppm input.pdf outputname -png")
# save plots
png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale, "umap_mk1_9a.png"), width=12, height=16, units = "in", res = 300)
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

# Rename clusters -------
tum$cellstate = "NULL"
tum$cellstate[tum$integrated_snn_res.0.3 == "1"] = "b_naive_activated_mhc2"
tum$cellstate[tum$integrated_snn_res.0.3 == "3"] = "b_naive_cd27"
tum$cellstate[tum$integrated_snn_res.0.3 == "5"] = "b_naive_tcl1a"
tum$cellstate[tum$integrated_snn_res.0.3 == "6"] = "b_naive_myc"

tum$cellstate[tum$integrated_snn_res.0.3 == "0"] = "b_plasma_igl"
tum$cellstate[tum$integrated_snn_res.0.3 == "2"] = "b_plasma_hsp"
tum$cellstate[tum$integrated_snn_res.0.3 == "4"] = "b_plasma_igh_igk"

tum$cellmajor = "NULL"
tum$cellmajor[tum$integrated_snn_res.0.3 %in% c("1", "3", "5", "6")] = "b_naive"
tum$cellmajor[tum$integrated_snn_res.0.3 %in% c("0", "2", "4")] = "b_plasma"


p6 = VlnPlot(tum, features = c("rna_IGLL5", "rna_IGLC2", 
                          "rna_HLA-DRA", "rna_HLA-DRB1",
                          "rna_HSPA6", "rna_HSPA1A",
                          "rna_CD27",
                          "rna_IGHG4", "rna_IGHG2",
                          "rna_TCL1A", "rna_SELL",
                          "rna_MYC", "rna_BCL2A1", "rna_EGR3"), group.by = "cellstate",
             pt.size = 0,cols = c("#5a189a", "#f4a261","#2a9d8f", "#f15bb5","#fee440", "#00bbf9","#00f5d4"))

png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"vln_all.png"), width=9, height=14, units = "in", res = 300)
print(p6)
dev.off()

p7 = DimPlot(tum, group.by = "cellstate", label = T, cols = c("#5a189a", "#f4a261","#2a9d8f", "#f15bb5","#fee440", "#00bbf9","#00f5d4"))
png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_final_9a.png"), width=6, height=3.5, units = "in", res = 300)
print(p7)
dev.off()

tum$cellstate = factor(tum$cellstate, levels = c("b_plasma_igl", "b_plasma_hsp", "b_plasma_igh_igk", 
                                                 "b_naive_activated_mhc2", "b_naive_cd27", "b_naive_tcl1a", "b_naive_myc"))
write_rds(tum, "b_cell_final.rds")
Idents(b_cell_final) <- "cellstate"

# ** Markers function (9)--------
dims = 20; k.param = 20;
char_ident = "cellstate"  #integrated_snn_res.0.1
mark_out = func_markers_vst(clus_out = clus_out, char_ident = char_ident, fts1 = ht_fts1, dims = dims, k.param = k.param, cols = colors_dark, 
                            save_path = save_path)



