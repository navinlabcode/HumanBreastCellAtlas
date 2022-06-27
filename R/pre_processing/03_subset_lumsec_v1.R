# geo ----
setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_lumsec_v1";dir.create(odir)
options(future.globals.maxSize = 90000 * 1024^2)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "figure_1_v1"
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_v1/"

# core ----
setwd("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_lumsec_v1"; dir.create(odir)
library(future)
plan("multiprocess", workers = 75)
options(future.globals.maxSize = Inf)
setwd(paste0("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "figure_1_v1"
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_v1/"


# read original object --------------
tum = read_rds(path = paste0(wd, data_folder,"/final_use.rds"))
tum$sheet_patient_id[tum$sheet_ts_sampleid %in% c("hbca25_c")] = "pt37"
tum_sub = tum

# save objects
Idents(tum_sub) <- "final_group1"
tum_sub = subset(tum_sub, subset = final_group1 %in% c("LumSec"))
write_rds(tum_sub, path = glue("{save_path}/lumsec_subset_v1.rds"))

png(paste0(save_path, "vlnplot_nCount_RNA_sample_id.png"), width=24, height=5, units = "in", res = 300)
VlnPlot(tum_sub, features = "nCount_RNA", group.by = "final_sample_id", pt.size = 0 )
dev.off()


# method 9 ----------
# Log T integration with no scaling
library(future)
plan("multiprocess", workers = 75)
#plan("sequential")

# ** Clus function --------
# run clustering function
seu_object = read_rds(glue("{save_path}/lumsec_subset_v1.rds"))
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

p1 = DimPlot(clus_out, split.by = "sheet_patient_id", group.by = "integrated_snn_res.0.4", 
             cols = colors_dark, label = F, shuffle = T, pt.size = 0.6, ncol = 10)
png("vln_allpts_v1.png", width=24, height=20, units = "in", res = 300)
p1
dev.off()

DimPlot(clus_out, group.by = "integrated_snn_res.0.4", cols = colors_dark, label = T)
VlnPlot(clus_out, group.by = "integrated_snn_res.0.4", features = c("nCount_RNA", "percent.mito"), pt.size = 0)
VlnPlot(clus_out, group.by = "integrated_snn_res.0.4", features = c("EPCAM", "SLPI", "LTF", "KRT5", 
"KRT6B", "KRT14", "KRT15"), pt.size = 0, cols = colors_dark, assay = "RNA")

# add module score ------
basal_sig = list(c("KRT14","KRT17","KRT5","ACTG2","TUBB2B","COL17A1","LAMA3","SERPINB5"))
lumhr_sig = list(c("AREG","ANKRD30A","AGR2","AGR3","TMC5","DNAJC12","RASEF","MLPH","CA12","C3orf52","ACADSB","MAPK8"))
lumsec_sig = list(c("SLPI","LTF","KRT15","MMP7","CCL28","PIGR","RARRES1","ALDH1A3","LCN2","SNORC","S100A1","GABRP"))
clus_out = AddModuleScore(object = clus_out, features = basal_sig, assay = "RNA", name = "basal_sig")
clus_out = AddModuleScore(object = clus_out, features = lumhr_sig, assay = "RNA", name = "lumhr_sig")
clus_out = AddModuleScore(object = clus_out, features = lumsec_sig, assay = "RNA", name = "lumsec_sig")
FeaturePlot(clus_out, features = c("basal_sig1", "lumhr_sig1", "lumsec_sig1")) + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
VlnPlot(clus_out, features = c("basal_sig1", "lumhr_sig1", "lumsec_sig1"), group.by = "integrated_snn_res.0.4", 
        pt.size = 0, cols = colors_dark)

VlnPlot(clus_out, features = c("KRT14","KRT17","KRT5","ACTG2","TUBB2B","COL17A1","LAMA3","SERPINB5"), group.by = "integrated_snn_res.0.2", pt.size = 0)
VlnPlot(clus_out, features = c("AREG","ANKRD30A","AGR2","AGR3","TMC5","DNAJC12","RASEF","MLPH","CA12","C3orf52","ACADSB","MAPK8"), group.by = "integrated_snn_res.0.2", pt.size = 0)
VlnPlot(clus_out, features = c("SLPI","LTF","KRT15","MMP7","CCL28","PIGR","RARRES1","ALDH1A3","LCN2","SNORC","S100A1","GABRP"), group.by = "integrated_snn_res.0.2", pt.size = 0)

df_meta = data.frame(clus_out@meta.data, 
                     celltype = clus_out$integrated_snn_res.0.4, 
                     stringsAsFactors = F) %>% 
  as_tibble() %>% clean_names()
func_bar(var_df = df_meta, var_x = "sheet_patient_id", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "celltype")
func_bar(var_df = df_meta, var_x = "sheet_exp_proc", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "celltype")


# ** Markers function (9)--------
dims = 35; k.param = 20;
char_ident = "integrated_snn_res.0.4"  #integrated_snn_res.0.1
ht_fts1 = c("EPCAM", "KRT17", "AR")
#clus_out = read_rds(file = glue("{save_path}/dim_35k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")) 
mark_out = func_markers_vst(clus_out = clus_out, char_ident = char_ident, fts1 = ht_fts1, dims = dims, k.param = k.param, cols = colors_dark, 
                            save_path = save_path)

clus_out = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_v1/dim_35k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")
markers = read.csv("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_v1/dims_35k_20_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.4_vst_method_9.csv")
# clus_out = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v1/dim_30k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")
# markers = read.csv("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v1/dims_30k_20_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.4_vst_method_9.csv")
tum = clus_out

markers_v2 = markers %>% dplyr::filter(gene %in% markers$gene)

top_markers = markers_v2 %>%
  dplyr::group_by(cluster) %>%
  top_n(10, avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character() %>% unique()

#Idents(tum.integrated.ns)
#tum.integrated.ns$clusters1 <- factor(tum.integrated.ns$clusters1, levels = levels(Idents(tum.integrated.ns)))
message("Making Heatmaps")
DefaultAssay(tum) <- "integrated"
Idents(tum) <- "integrated_snn_res.0.4"
p1 = DoHeatmap(tum, features = c("EPCAM", top_markers), assay = "integrated", label = TRUE, size = 6,raster = TRUE,
              group.by = char_ident, slot = "scale.data", disp.min = -1.5, disp.max = 1.5, group.colors = colors_dark) + 
  scale_color_manual(values = colors_dark) +
  scale_fill_gradient2(low="steelblue1", high = "tomato3", mid = "white", midpoint = 0) + 
  theme(axis.text.y.left = element_text(size=5))

png("heatmap_all_labelled_res4_cells_v1.png", width=12, height=10, units = "in", res = 300)
p1
dev.off()

seu_avg = AverageExpression(tum, features = top_markers, return.seurat = T)
DoHeatmap(seu_avg, features = c(top_markers), assay = "integrated", label = TRUE,
          slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = colors_dark, raster = F)  +
  scale_fill_gradient2(low = "#29A7B8", mid = "white", high = "#E85041", na.value = "00000000") +
  theme(axis.text.y = element_text(size=7.5)) -> avg_heatmap

png("heatmap_Avg_all_labelled_res4_cells_v1.png", width=6, height=14, units = "in", res = 300)
avg_heatmap
dev.off()

v1 = VlnPlot(tum, features = c("rna_DUSP4", "rna_SPINK5", "rna_KRT14", "rna_MGP", "rna_SCGB2A2", "rna_CLDN1", "rna_CXCL1", "rna_PI3", "rna_KRT15", "rna_MT1X", "rna_HMOX1", 
        "rna_PTN", "rna_CA3", "rna_SLPI", "rna_CCL20", "rna_IFI27", "rna_FABP5", "rna_LALBA", "rna_PLCG2", "rna_KIT", "rna_HES1", "rna_XIST","rna_VEGFA"),pt.size = 0, cols = colors_darj)

png("VlnPlot_all_labelled_res2_cells_v1.png", width=11, height=11, units = "in", res = 200)
v1
dev.off()

v2 = VlnPlot(tum, features = c("nCount_RNA", "nFeature_RNA"),pt.size = 0, cols = colors_darj)

png("VlnPlot_count_gene_all_labelled_res2_cells_v1.png", width=6, height=3, units = "in", res = 200)
v2
dev.off()




