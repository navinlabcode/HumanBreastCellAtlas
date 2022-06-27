# geo ----
setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_lumsec_20samp_v1";dir.create(odir)
options(future.globals.maxSize = 90000 * 1024^2)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "03_subset_lumsec_v2"
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_20samp_v1/"

# core ----
setwd("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_lumsec_20samp_v1"; dir.create(odir)
library(future)
plan("multiprocess", workers = 50)
options(future.globals.maxSize = Inf)
setwd(paste0("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "03_subset_lumsec_v2"
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_20samp_v1/"
# stress_sig = read.delim("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_fibroblasts_final1/stromal_stress_signature_032720.txt")
# stress_sig = as.character(stress_sig$stromal.stress.sig)


# read original object --------------
tum = read_rds(paste0(wd, "03_subset_lumsec_v2/epi_lumsec_sub_v1.rds"))

png(paste0(save_path, "vlnplot_nCount_RNA_sample_id.png"), width=14, height=3, units = "in", res = 300)
VlnPlot(tum, features = "nCount_RNA", group.by = "final_sample_id", pt.size = 0 )
dev.off()

# subset and filter ----------
Idents(tum) <- "batch"
tum_sub = subset(tum, idents = c("b2"), invert = F) #endo_vas
Idents(tum_sub) <- "sheet_exp_proc"
tum_sub = subset(tum_sub, idents = c("3hr_digestion"), invert = F) #endo_vas
Idents(tum_sub) <- "sheet_figure1_grp_updated"
tum_sub = subset(tum_sub, idents = c("Reduction Mammoplasty"), invert = F) #endo_vas

write_rds(tum_sub, "epi_lumsec_sub_3hr_b2_reduction_v1.rds")
tum_sub = read_rds("epi_lumsec_sub_3hr_b2_reduction_v1.rds")

# recluster ----
seu_obj = tum_sub
DefaultAssay(seu_obj) <- "RNA"
seu_obj = NormalizeData(seu_obj) %>% FindVariableFeatures(nfeatures = 5000) 
all.genes <- rownames(seu_obj)
seu_obj <- ScaleData(seu_obj, features = all.genes)
seu_obj <- RunPCA(seu_obj, features = VariableFeatures(object = seu_obj))
seu_obj <- FindNeighbors(seu_obj, dims = 1:35)
seu_obj <- FindClusters(seu_obj, resolution = c(0.5, 0.1, 0.05, 0.075, 0.2, 0.15, 0.3))
seu_obj <- RunUMAP(seu_obj, dims = 1:35)
DimPlot(seu_obj, reduction = "umap", group.by = "RNA_snn_res.0.2")
DimPlot(seu_obj, reduction = "umap", group.by = "sheet_exp_proc")
DimPlot(seu_obj, reduction = "umap", group.by = "RNA_snn_res.0.2", split.by = "sheet_patient_id")
seu_obj$RNA_snn_res.0.2 = factor(seu_obj$RNA_snn_res.0.2, levels = c("0", "1", "2", "3", "4", "5", "6"))
Idents(seu_obj) <- "RNA_snn_res.0.2"

FeaturePlot(seu_obj, "rna_LALBA")

# remove and re run
seu_obj1 = subset(seu_obj, idents = c("0"), invert = T)
DimPlot(seu_obj1, reduction = "umap", group.by = "RNA_snn_res.0.2")
seu_obj = seu_obj1
# go back to line 47 recluster


mks = FindAllMarkers(seu_obj,assay = "RNA")
write.csv(mks, "mks.csv")
mks %>% 
  #arrange(cluster) %>%
  filter(!str_detect(gene, "^RP")) %>%
  filter(!str_detect(gene, "^MT-")) %>%
  filter(!str_detect(gene, "^hbb")) %>%
  group_by(cluster) %>%  
  slice_max(n = 20, order_by = avg_logFC) -> top10

DoHeatmap(seu_obj, features = top10$gene, group.colors = colors_dark) + NoLegend()
tavg = AverageExpression(seu_obj, features = unique(c("EPCAM", "KIT", "TOP2A", top10$gene)), return.seurat = T)
DoHeatmap(tavg, features = unique(c("EPCAM", "KIT", "TOP2A", top10$gene)), assay = "RNA", label = TRUE,
          slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = colors_dark, raster = F)  +
  scale_fill_gradient2(low = "#29A7B8", mid = "white", high = "#E85041", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5)) -> tcell_heatmap
png("heatmap_all_labelled_lumsec_0.2_v3.png", width=5, height=12, units = "in", res = 300)
tcell_heatmap
dev.off()

# check my clusters ------
tk_obj = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_v2/clus_out.rds")
kevin_obj = load("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_v10/hbca_seclum.RObj")
kevin_obj = hbca.seclum

tk_obj_df = data.frame(cellname = tk_obj$cellname, tk_states = tk_obj$cellstate_tk2)
clus_out_meta = seu_obj@meta.data
clus_out_meta = left_join(clus_out_meta, tk_obj_df)

seu_obj@meta.data = clus_out_meta
seu_obj$RNA_snn_res.0.2 = as.character(seu_obj$RNA_snn_res.0.2)
seu_obj$RNA_snn_res.0.2[(is.na(seu_obj$RNA_snn_res.0.2))] = "unkn"
rownames(seu_obj@meta.data) = seu_obj$cellname

kevin_obj_df = data.frame(cellname = kevin_obj$cellname, kevin_states = kevin_obj$seclum_states)
clus_out_meta = seu_obj@meta.data
clus_out_meta = left_join(clus_out_meta, kevin_obj_df)

seu_obj@meta.data = clus_out_meta
seu_obj$kevin_states = as.character(seu_obj$kevin_states)
seu_obj$kevin_states[(is.na(seu_obj$kevin_states))] = "unkn"
rownames(seu_obj@meta.data) = seu_obj$cellname


p00 = DimPlot(tk_obj, reduction  = "umap", group.by = "cellstate_tk2", ncol = 1, label = T, cols = colors_dark, pt.size = 0.0001)
p0 = DimPlot(seu_obj, reduction  = "umap", group.by = "kevin_states", ncol = 1, label = T, cols = colors_dark, pt.size = 0.0001)
p1 = DimPlot(seu_obj, reduction  = "umap", group.by = "tk_states", ncol = 1, label = T, cols = colors_dark, pt.size = 0.0001)
p2 = DimPlot(seu_obj, reduction  = "umap", group.by = "RNA_snn_res.0.2", ncol = 1, label = T, cols = colors_dark, pt.size = 0.0001)
p3 = DimPlot(seu_obj, reduction  = "umap", group.by = "tk_states", split.by = "final_sample_id", ncol = 6, label = F, cols = colors_dark, pt.size = 1)
png("umap_tk_v3.png", width=20, height=4, units = "in", res = 300) #
cowplot::plot_grid(p00, p0, p1, p2, ncol = 4)
dev.off()

png("umap_tk_split_by_pt_v3.png", width=11, height=11, units = "in", res = 300) #
p3
dev.off()

png("umap_rna0.2_split_by_pt_v3.png", width=11, height=11, units = "in", res = 300) #
DimPlot(seu_obj, reduction  = "umap", group.by = "RNA_snn_res.0.2", split.by = "final_sample_id", ncol = 6, label = F, cols = colors_dark, pt.size = 1)
dev.off()

# END HERE -----------



# assign clusters ------
clus_out = tum
clus_out$cellstate_tk1 = "NULL"
clus_out$cellstate_tk1[clus_out$integrated_snn_res.0.2 %in% c("0")] = "ls_KRT23"
clus_out$cellstate_tk1[clus_out$integrated_snn_res.0.2 %in% c("1")] = "ls_SPINK5"
clus_out$cellstate_tk1[clus_out$integrated_snn_res.0.2 %in% c("2")] = "ls_HMOX1"
clus_out$cellstate_tk1[clus_out$integrated_snn_res.0.2 %in% c("3")] = "ls_MT1X"
clus_out$cellstate_tk1[clus_out$integrated_snn_res.0.2 %in% c("4")] = "ls_KIT"
clus_out$cellstate_tk1[clus_out$integrated_snn_res.0.2 %in% c("5")] = "ls_lactation"
clus_out$cellstate_tk1[clus_out$integrated_snn_res.0.2 %in% c("6")] = "ls_cd74"
clus_out$cellstate_tk1[clus_out$integrated_snn_res.0.2 %in% c("7")] = "ls_prol"

lvl = names(table(clus_out$cellstate_tk1))
clus_out$cellstate_tk1 = factor(clus_out$cellstate_tk1, levels = lvl)

clus_out$cellstate_tk2 = "NULL"
clus_out$cellstate_tk2[clus_out$integrated_snn_res.0.2 %in% c("0")] = "ls_KRT23"
clus_out$cellstate_tk2[clus_out$integrated_snn_res.0.2 %in% c("1")] = "ls_SPINK5"
clus_out$cellstate_tk2[clus_out$integrated_snn_res.0.2 %in% c("2")] = "ls_HMOX1"
clus_out$cellstate_tk2[clus_out$integrated_snn_res.0.2 %in% c("3", "4")] = "ls_KIT"
clus_out$cellstate_tk2[clus_out$integrated_snn_res.0.2 %in% c("5")] = "ls_lactation"
clus_out$cellstate_tk2[clus_out$integrated_snn_res.0.2 %in% c("6")] = "ls_cd74"
clus_out$cellstate_tk2[clus_out$integrated_snn_res.0.2 %in% c("7")] = "ls_prol"

lvl = names(table(clus_out$cellstate_tk2))
clus_out$cellstate_tk2 = factor(clus_out$cellstate_tk2, levels = lvl)

write_rds(clus_out, "clus_out.rds")
clus_out = read_rds("clus_out.rds")
Idents(clus_out) <- "cellstate_tk2"

clus_out_seu = FindAllMarkers(object = clus_out, assay = "RNA", logfc.threshold = 0)
write.csv(clus_out_seu, "clus_out_seu.csv")


markers_v2 = clus_out_seu %>% dplyr::filter(gene %in% clus_out_seu$gene)

top_markers = markers_v2 %>%
  dplyr::group_by(cluster) %>%
  top_n(7, avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()

DefaultAssay(clus_out) <- "integrated"
p1= DoHeatmap(clus_out, features = c(top_markers), size = 2,raster = FALSE, label=TRUE,
              group.by = "cellstate", slot = "scale.data", disp.min = -1.5, disp.max = 1.5, group.colors = colors_dark) + 
  scale_color_manual(values = colors_dark) +
  scale_fill_gradient2(low="steelblue1", high = "tomato3", mid = "white", midpoint = 0) + 
  theme(axis.text.y.left = element_text(size=6))

tavg = AverageExpression(clus_out, features = unique(c("EPCAM", "KIT", "TOP2A", top_markers)), return.seurat = T)
DoHeatmap(tavg, features = unique(c("EPCAM", "KIT", "TOP2A", top_markers)), assay = "integrated", label = TRUE,
          slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = colors_dark, raster = F)  +
  scale_fill_gradient2(low = "#29A7B8", mid = "white", high = "#E85041", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5)) -> tcell_heatmap
png("heatmap_all_labelled_lumsec_cellstatetk2_v1.png", width=12, height=12, units = "in", res = 300)
tcell_heatmap
dev.off()


ggsave(plot = p1, filename = "heatmap_cellstate_Tcells_v2.png", device = "png", width = 20, height = 18)

png("heatmap_cellstate_Tcells_v1.png", width=18, height=12, units = "in", res = 300) #
print(p1)
dev.off()

png("umap_Tcells_v2.png", width=9, height=7, units = "in", res = 300) #
DimPlot(clus_out, group.by = "cellstate_tk2", cols = colors_dark, label = T)
dev.off()

png("dimplot_split_v1.png", width=18, height=24, units = "in", res = 300) #
DimPlot(clus_out, split.by = "sheet_patient_id", group.by = "cellstate_tk2", cols = colors_dark, label = F, ncol = 6, pt.size = 1)
dev.off()

png("dimplot_split_v2.png", width=12, height=12, units = "in", res = 300) #
DimPlot(clus_out, split.by = "sheet_exp_proc", group.by = "cellstate_tk2", cols = colors_dark, label = F, ncol = 2, pt.size = 1)
dev.off()

v1 = VlnPlot(tum, features = c("rna_CD74", "rna_CCL20", "rna_HMOX1", "rna_S100A8", 
                               "rna_KIT", "rna_PLCG2", "rna_PTN", "rna_SCGB2A2", 
                               "rna_KRT23", "rna_CXCL8", "rna_CXCL1", 
                               "rna_LALBA", "rna_CSN2", "rna_TOP2A", 
                               "rna_STMN1", "rna_CYP24A1", "rna_DUSP4"), group.by = "cellstate_tk2",
             pt.size = 0, cols = colors_dark, y.max = 3)

png("VlnPlot_all_labelled_rescellstate2_cells_v1.png", width=11, height=11, units = "in", res = 200)
v1
dev.off()

v2 = VlnPlot(tum, features = c("nCount_RNA", "nFeature_RNA"),pt.size = 0, 
             cols = colors_dark, group.by = "cellstate_tk2")

png("VlnPlot_count_gene_all_labelled_res_cellstate2_cells_v1.png", width=6, height=3, units = "in", res = 200)
v2
dev.off()


# barplots ------
# **Barplots ---------------------------------------------------------------
col_t = colors_dark
df_meta = data.frame(clus_out@meta.data, 
                     celltype = clus_out$cellstate_tk2, 
                     stringsAsFactors = F) %>% 
  as_tibble() %>% clean_names()

df_meta = df_meta %>% dplyr::mutate(sheet_bmi_final = case_when(sheet_bmi_final == "Normal" ~ "normal",
                                                                TRUE ~ as.character(sheet_bmi_final)))

df_meta1 = df_meta
# create a summary table
func_bar(var_df = df_meta1, var_x = "sheet_patient_id", var_fill = "cellstate_tk1", var_pos = "stack", var_cols = col_t, celltype = "cellstate_tk1")
func_bar(var_df = df_meta1, var_x = "sheet_patient_id", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "cellstate_tk2")
func_bar(var_df = df_meta1, var_x = "sheet_patient_id", var_fill = "sheet_parity", var_pos = "fill", var_cols = col_t, celltype = "sheet_parity")
func_bar(var_df = df_meta1, var_x = "sheet_patient_id", var_fill = "sheet_parity", var_pos = "fill", var_cols = col_t, celltype = "sheet_parity")
func_bar(var_df = df_meta1, var_x = "sheet_parity", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "sheet_parity")
func_bar(var_df = df_meta1, var_x = "sheet_parity", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "sheet_parity")

func_bar(var_df = df_meta1, var_x = "sheet_exp_proc", var_fill = "sheet_parity", var_pos = "fill", var_cols = colors_dark, celltype = "sheet_parity")

func_bar(var_df = df_meta1, var_x = "sheet_patient_id", var_fill = "sheet_exp_proc", var_pos = "fill", var_cols = colors_neon, celltype = "sheet_exp_proc")

func_bar(var_df = df_meta1, var_x = "final_sample_id", var_fill = "celltype", var_pos = "stack", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "final_sample_id", var_fill = "celltype", var_pos = "fill", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_bmi_final", var_fill = "celltype", var_pos = "stack", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_bmi_final", var_fill = "celltype", var_pos = "fill", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_brca_updated", var_fill = "celltype", var_pos = "stack", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_brca_updated", var_fill = "celltype", var_pos = "fill", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_exp_proc", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_exp_proc", var_fill = "celltype", var_pos = "fill", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_cancer_risk", var_fill = "celltype", var_pos = "stack", var_cols = col_t, celltype = "ultra_res2")
func_bar(var_df = df_meta1, var_x = "sheet_cancer_risk", var_fill = "celltype", var_pos = "fill", var_cols = col_t, celltype = "ultra_res2")
func_bar(var_df = df_meta1, var_x = "sheet_figure1_grp_updated", var_fill = "celltype", var_pos = "stack", var_cols = col_t, celltype = "sheet_figure1_grp_updated")
func_bar(var_df = df_meta1, var_x = "sheet_figure1_grp_updated", var_fill = "celltype", var_pos = "fill", var_cols = col_t, celltype = "sheet_figure1_grp_updated")
func_bar(var_df = df_meta1, var_x = "sheet_menopause", var_fill = "celltype", var_pos = "stack", var_cols = col_t, celltype = "sheet_menopause")
func_bar(var_df = df_meta1, var_x = "sheet_menopause", var_fill = "celltype", var_pos = "fill", var_cols = col_t, celltype = "sheet_menopause")
func_bar(var_df = df_meta1, var_x = "sheet_menopause", var_fill = "celltype", var_pos = "stack", var_cols = col_t, celltype = "sheet_menopause")
func_bar(var_df = df_meta1, var_x = "sheet_menopause", var_fill = "celltype", var_pos = "fill", var_cols = col_t, celltype = "sheet_menopause")
func_bar(var_df = df_meta1, var_x = "sheet_cancer_history", var_fill = "cellstate", var_pos = "stack", var_cols = col_t, celltype = "sheet_cancer_history")
func_bar(var_df = df_meta1, var_x = "sheet_cancer_history", var_fill = "cellstate", var_pos = "fill", var_cols = col_t, celltype = "sheet_cancer_history")


# check tks clusters ------
tk_obj = load("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_v10/hbca_seclum.RObj")
tk_obj = hbca.seclum

tk_obj_df = data.frame(cellname = tk_obj$cellname, tk_states = tk_obj$seclum_states)
clus_out_meta = clus_out@meta.data
clus_out_meta = left_join(clus_out_meta, tk_obj_df)

clus_out@meta.data = clus_out_meta
clus_out$tk_states[(is.na(clus_out$tk_states))] = "unkn"
rownames(clus_out@meta.data) = clus_out$cellname

png("umap_tk_Tcells_v2.png", width=9, height=7, units = "in", res = 300) #
DimPlot(clus_out, group.by = "tk_states", label = T, cols = colors_dark, pt.size = 0.0001)
dev.off()

png("umap_tk_split_by_pt_v2.png", width=18, height=18, units = "in", res = 300) #
DimPlot(clus_out, group.by = "tk_states", split.by = "sheet_patient_id", ncol = 6, label = T, cols = colors_dark, pt.size = 0.0001)
dev.off()

png("umap_tk_split_by_pt_v2.png", width=18, height=24, units = "in", res = 300) #
DimPlot(tk_obj, group.by = "seclum_states", split.by = "sheet_patient_id", ncol = 6, cols = colors_dark, pt.size = 1)
dev.off()

# confusion matrix ------
cc = table(clus_out$cellstate_tk2, clus_out$tk_states)
g1 = ggplot(as.data.frame(cc), aes(Var1,Var2, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  #scale_fill_viridis_c()+
  scale_fill_gradient2(low="white",high="#009194") +
  labs(x = "Tapsi's labels",y = "tk's labels")
png(paste0(save_path,"/tapsi_tk_confusionmat_v2.png"), width=8, height=6, units = "in", res = 300)
print(g1)
dev.off()

# pathways -------
# pathway ------
genes = read.csv("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_v2/clus_out_seu.csv")
genes = genes[,-1]
genes$cluster = as.character(genes$cluster)

gsea_resultsc2hall = fun_gsea(genes = genes, set = "H", sub_cat = NULL, num_avg_logFC = 0.25)
gsea_resultsc2C5 = fun_gsea(genes = genes, set = "C5", sub_cat = "BP", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "fibro")
gsea_resultsc2C3 = fun_gsea(genes = genes, set = "C3", sub_cat = "TFT:GTRD", num_avg_logFC = 0.25, celltype = "lumsec", doplot = FALSE, ret_collap = FALSE)
gsea_resultsc2C8 = fun_gsea(genes = genes, set = "C8", sub_cat = NULL, num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "fibro")
gsea_resultsc2C2_reac = fun_gsea(genes = genes, set = "C2", sub_cat = "CP:REACTOME", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "fibro")
gsea_resultsc2C2_kegg = fun_gsea(genes = genes, set = "C2", sub_cat = "CP:KEGG", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "lumsec")

clus = levels(clus_out$cellstate_tk2)

gs = gsea_resultsc2hall
gs = gsea_resultsc2C2;gs_kegg = gsea_resultsc2C2_kegg;gs_reac = gsea_resultsc2C2_reac
gs = gs_kegg
gs = gs_reac
gs = gsea_resultsc2C5
gs = gsea_resultsc2C3
gs = gsea_resultsc2C7
gs = gsea_resultsc2C8

gs1 = gs %>% dplyr::select(-leadingEdge) %>% dplyr::group_by(cluster) %>% 
  dplyr::filter(padj <= 0.05) %>%
  #top_n(30, wt = abs(NES)) %>% 
  pull(pathway)
gs2 = gs %>% dplyr::filter(pathway %in% gs1)
gs2$cluster = factor(gs2$cluster, levels = clus, ordered = TRUE)
gs2 = gs2 %>% dplyr::filter(grepl('HALLMARK|KEGG|REACTOME|BIOCARTA|GO', pathway))

#gs2 = gs2 %>% mutate(pathway = str_split(pathway, "_", n=2, simplify = TRUE)[,2]) # for C7
gs2 = gs2 %>%  arrange(cluster, -NES)

library(forcats);library(circlize);library(ComplexHeatmap)

p1 = gs2 %>% 
  dplyr::select(pathway, NES, cluster) %>%
  arrange(cluster, NES) %>%
  spread(key = cluster, value = NES, fill=0) 

p1.1 = p1
colnames(p1.1) = c("pathway", clus)

# rearrange the datafram
test = p1

# for basal
test1 = test %>% dplyr::arrange(-`ls_cd74`) %>% dplyr::filter(`ls_cd74`>0) %>% pull(pathway)
test2 = test %>% dplyr::filter(!pathway %in% test1) %>% arrange(-ls_HMOX1) %>% dplyr::filter(ls_HMOX1>0, `ls_cd74`<0, `ls_KIT`<0, `ls_KRT23`<0, `ls_lactation`<0,  `ls_prol`<0, `ls_SPINK5`<0) %>% dplyr::top_n(25, `ls_HMOX1`) %>% pull(pathway)
test3 = test %>% dplyr::filter(!pathway %in% c(test1, test2)) %>% arrange(-`ls_KIT`) %>% dplyr::filter(ls_HMOX1<0, `ls_cd74`<0, `ls_KIT`>0, `ls_KRT23`<0, `ls_lactation`<0,  `ls_prol`<0, `ls_SPINK5`<0) %>% dplyr::top_n(25, `ls_KIT`) %>% pull(pathway)
test4 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3)) %>% arrange(-`ls_KRT23`) %>% dplyr::filter(ls_HMOX1<0, `ls_cd74`<0, `ls_KIT`<0, `ls_KRT23`>0, `ls_lactation`<0,  `ls_prol`<0, `ls_SPINK5`<0) %>% dplyr::top_n(25, `ls_KRT23`) %>% pull(pathway)
test5 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3, test4)) %>% arrange(-`ls_lactation`) %>% dplyr::filter(ls_HMOX1<0, `ls_cd74`<0, `ls_KIT`<0, `ls_KRT23`<0, `ls_lactation`>0,  `ls_prol`<0, `ls_SPINK5`<0) %>% dplyr::top_n(25, `ls_lactation`) %>% pull(pathway)
test6 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3, test4, test5)) %>% arrange(-`ls_prol`) %>% dplyr::filter(ls_HMOX1<0, `ls_cd74`<0, `ls_KIT`<0, `ls_KRT23`<0, `ls_lactation`<0,  `ls_prol`>0, `ls_SPINK5`<0) %>% dplyr::top_n(25, `ls_prol`) %>% pull(pathway)
test7 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3, test4, test5, test6)) %>% arrange(-`ls_SPINK5`) %>% dplyr::filter(ls_HMOX1<0, ls_cd74<0, `ls_KIT`<0, `ls_KRT23`<0, `ls_lactation`<0,  `ls_prol`<0, `ls_SPINK5`>0) %>% dplyr::top_n(25, `ls_SPINK5`) %>% pull(pathway)
seq_path = c(test1, test2, test3,test4, test5, test6, test7) #


# make heatmap 
p1.1 = p1.1 %>%
  dplyr::slice(match(seq_path, pathway))
p1 = p1.1
p2 = data.matrix(p1)
rownames(p2) = p1$pathway
p2 = p2[,-1]
col = colorRamp2(c(-2, 0, 2),c("steelblue",  "white",  "#F80F61")) 
# c("steelblue","white", "tomato3"))
ht_hall = Heatmap(p2, 
                  show_row_names = T, "NES", 
                  col = col,
                  #km = 3,
                  #clustering_method_rows = "complete",
                  #clustering_distance_rows = "maximum",
                  row_names_gp = gpar(cex = 0.4),
                  cluster_rows = F,
                  row_dend_side = "right",
                  row_names_side = "left",
                  width = unit(5, "cm"), 
                  cluster_columns = F, row_names_max_width = unit(20, "cm"),
                  column_names_gp = gpar(cex = 0.6));ht_hall

png("htmaps_pathways_bp.png", res = 300, units = "in", width = 12, height = 22)
print(ht_hall)
dev.off()


