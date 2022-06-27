# geo ----
setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_T_v2";dir.create(odir)
options(future.globals.maxSize = 90000 * 1024^2)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "03_subset_T_v1"
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v2/"

# core ----
setwd("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_T_v2"; dir.create(odir)
library(future)
plan("multiprocess", workers = 50)
options(future.globals.maxSize = Inf)
setwd(paste0("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "03_subset_T_v1"
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v2/"


# read original object --------------
tum = readRDS(file = paste0(wd, data_folder,"/tum_sub_v1.rds"))

# subset ---------
tcells = as.numeric(tum@assays$RNA@counts["CD3E",] > 0) + as.numeric(tum@assays$RNA@counts["CD3G",] > 0) + 
  as.numeric(tum@assays$RNA@counts["CD3D",] > 0) + as.numeric(tum@assays$RNA@counts["CD8A",] > 0) + 
  as.numeric(tum@assays$RNA@counts["CD8B",] > 0) + as.numeric(tum@assays$RNA@counts["CD4",] > 0) +
  as.numeric(tum@assays$RNA@counts["FCGR3A",] > 0) + as.numeric(tum@assays$RNA@counts["MKI67",] > 0) + 
  as.numeric(tum@assays$RNA@counts["FCER1G",] > 0) + 
  as.numeric(tum@assays$RNA@counts["NKG7",] > 0) + as.numeric(tum@assays$RNA@counts["IL7R",] > 0)
tcells2 = row.names(tum@meta.data[tcells > 1,])    

tum2 = SubsetData(tum, cells = tcells2)
p1 = DimPlot(tum2, group.by = "integrated_snn_res.0.4")
p2 = DimPlot(tum, group.by = "integrated_snn_res.0.4")
p0 = cowplot::plot_grid(p1, p2) 
png(paste0(save_path, "filtered_cd4cd8_featureplot.png"), width=16, height=8, units = "in", res = 300)
p0
dev.off()

# assign CD4 cd8
cd4 = as.numeric(tum2@assays$RNA@counts["CD3E",] > 0) + as.numeric(tum2@assays$RNA@counts["CD3G",] > 0) + 
  as.numeric(tum2@assays$RNA@counts["CD3D",] > 0) + as.numeric(tum2@assays$RNA@counts["CD4",] > 0) + 
  as.numeric(tum2@assays$RNA@counts["MKI67",] > 0)
cd4 = row.names(tum2@meta.data[cd4 > 1,]) 
cd8 = as.numeric(tum2@assays$RNA@counts["CD3E",] > 0) + as.numeric(tum2@assays$RNA@counts["CD3G",] > 0) + 
  as.numeric(tum2@assays$RNA@counts["CD3D",] > 0) + as.numeric(tum2@assays$RNA@counts["CD8A",] > 0) + 
  as.numeric(tum2@assays$RNA@counts["MKI67",] > 0) + as.numeric(tum2@assays$RNA@counts["CD8B",] > 0)
cd8 = row.names(tum2@meta.data[cd8 > 1,])

CD8_count_cells = colnames(tum2@assays$RNA@counts[,(tum2@assays$RNA@counts["CD8A",] >0 | tum2@assays$RNA@counts["CD8B",] >0) & !(tum2@assays$RNA@counts["CD4",] >0)])
CD4_count_cells = colnames(tum2@assays$RNA@counts[,!(tum2@assays$RNA@counts["CD8A",] >0 | tum2@assays$RNA@counts["CD8B",] >0) & (tum2@assays$RNA@counts["CD4",] >0)])


tum2@meta.data$CD4CD8_count = "unknown"
tum2@meta.data[CD8_count_cells,]$CD4CD8_count = "CD8"
tum2@meta.data[CD4_count_cells,]$CD4CD8_count = "CD4"

tum2 = SetIdent(tum2, value="CD4CD8_count")
tum2.CD4CD8_markers = FindAllMarkers(tum2, nthreads=20)
write.csv(tum2.CD4CD8_markers, "tum2.CD4CD8_markers.csv")

tum2.CD4CD8_markers %>% 
  filter(cluster == "CD8") %>% 
  filter(avg_logFC >0) %>% 
  pull(gene) ->
  CD8_genes

tum2.CD4CD8_markers %>% 
  filter(cluster == "CD4") %>% 
  filter(avg_logFC >0) %>% 
  pull(gene) ->
  CD4_genes

tum2.CD4CD8_markers %>% 
  filter(cluster == "unknown") %>% 
  filter(avg_logFC >0) %>% 
  pull(gene) ->
  unknown_genes

#score for CD4/CD8
tum2 = AddModuleScore(tum2, features = list(CD4_genes), name = "CD4_score", assay = "integrated", nbin=20)
tum2 = AddModuleScore(tum2, features = list(CD8_genes), name = "CD8_score", assay = "integrated", nbin=20)

#add final CD4CD8 labels
tum2@meta.data$CD4CD8 = "unknown"
tum2@meta.data[(tum2@meta.data$CD8_score1 > 0) & (tum2@meta.data$CD8_score1 > tum2@meta.data$CD4_score1),]$CD4CD8 = "CD8"
tum2@meta.data[(tum2@meta.data$CD4_score1 > 0) & (tum2@meta.data$CD4_score1 > tum2@meta.data$CD8_score1),]$CD4CD8 = "CD4"

#switch cells with high scores but clear opposite expression back to unknown
tum2@meta.data[tum2@meta.data$CD4CD8 =="CD4" & (row.names(tum2@meta.data) %in% CD8_count_cells),]$CD4CD8 = "unknown"
tum2@meta.data[tum2@meta.data$CD4CD8 =="CD8" & (row.names(tum2@meta.data) %in% CD4_count_cells),]$CD4CD8 = "unknown"

p1 = DimPlot(tum2, group.by = "CD4CD8")
p2 = FeaturePlot(tum2, c("rna_CD4", "rna_CD8A", "rna_CD8B", "rna_GNLY"))
p0 = cowplot::plot_grid(p1, p2)
png(paste0(save_path, "cd4cd8_featureplot.png"), width=16, height=8, units = "in", res = 300)
p0
dev.off()

# save objects
write_rds(tum2, path = glue("{save_path}/t_subset_v1.rds"))

png(paste0(save_path, "vlnplot_counts_sheet_patient_id.png"), width=14, height=5, units = "in", res = 300)
VlnPlot(tum2, features = "nCount_RNA", group.by = "sheet_patient_id", pt.size = 0 )
dev.off()


# method 9 ----------
# Log T integration with no scaling
library(future)
plan("multiprocess", workers = 70)
#plan("sequential")

# ** Clus function --------
# run clustering function
seu_object = read_rds(glue("{save_path}/t_subset_v1.rds"))
DefaultAssay(seu_object) <- "RNA"
dims = 35; k.param = 20; cols = colors_dark; var_split = "sheet_exp_proc"; var_scale = "nCount_RNA"
fts1 = c("nFeature_RNA","nCount_RNA", "percent.mito","percent.rb", "percent.stress")
fts2 = c("rna_EPCAM","rna_PTPRC","rna_PECAM1","rna_DCN","rna_KRT5",  "rna_CD86", "rna_COL1A1", "rna_PIP")
fts3 = c("rna_CD8A", "rna_CD8B", "rna_CD4", "rna_CCR7", "rna_IL7R", "rna_NKG7", "rna_GZMB", "rna_GZMK")
fts4 = c("rna_MS4A1","rna_CD40LG","rna_FOXP3","rna_IFNG","rna_LAG3","rna_TIGIT",  "rna_NCAM1", "rna_KLRB1", "rna_CD3D")

ref_set = NULL

clus_out = func_cluster_vst(tum_int = seu_object, var_scale = var_scale, var_split = var_split, dims = dims, k.param = k.param, 
                            markers1 = fts1, markers2 = fts2, markers3 = fts3, markers4 = fts4, save_path = save_path,
                            cols = colors_dark, ref_set = ref_set)

char_ident = "integrated_snn_res.0.6"


# ** Markers function (9)--------
dims = 35; k.param = 20;
char_ident = "integrated_snn_res.0.6"  #integrated_snn_res.0.1
ht_fts1 = c("CD3D", "CD8A", "CD8B")
#clus_out = read_rds(path = glue("{save_path}/dim_{dims}k_{k.param}_tum.integrated.ns_SCT_v4_intgeratedassay_{char_ident}_vst_method9.rds")) 
mark_out = func_markers_vst(clus_out = clus_out, char_ident = char_ident, fts1 = ht_fts1, dims = dims, k.param = k.param, cols = colors_dark, 
                            save_path = save_path)

clus_out = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v2/dim_35k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")
markers = read.csv("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v2/dims_35k_20_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.6_vst_method_9.csv")
# clus_out = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v2/dim_30k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")
# markers = read.csv("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_T_v2/dims_30k_20_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.4_vst_method_9.csv")
tum = clus_out

markers_v2 = markers %>% dplyr::filter(gene %in% markers$gene)

top_markers = markers_v2 %>%
  dplyr::group_by(cluster) %>%
  slice_max(avg_logFC, n=7) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()

FeaturePlot(tum, "rna_CD27", order = T)
#IL17F

VlnPlot(tum, "rna_CD69", pt.size = 0, group.by = "integrated_snn_res.0.6")

#Idents(tum.integrated.ns)
#tum.integrated.ns$clusters1 <- factor(tum.integrated.ns$clusters1, levels = levels(Idents(tum.integrated.ns)))
message("Making Heatmaps")
Idents(tum) <- "integrated_snn_res.0.6"
DefaultAssay(tum) <- "integrated"
p1 = DoHeatmap(tum, features = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B", top_markers), size = 2,raster = FALSE,
              group.by = char_ident, slot = "scale.data", disp.min = -1.5, disp.max = 1.5, group.colors = colors_dark) + 
  scale_color_manual(values = colors_dark) +
  scale_fill_gradient2(low="steelblue1", high = "tomato3", mid = "white", midpoint = 0) + 
  theme(axis.text.y.left = element_text(size=6))

png("heatmap_allcells_labelled_Tcells_0.6_v1.png", width=26, height=26, units = "in", res = 300)
p1
dev.off()
ggsave(plot = p1, filename = "figure_heatmap_cell_types_avg_integratedtop7.png", device = "png", width = 14, height = 10)
  
tavg = AverageExpression(tum, features = unique(c("CD3D", "CD3G", "CD4", "CD8A", "CD8B", top_markers)), return.seurat = T)
DoHeatmap(tavg, features = unique(c("CD3D", "CD3G", "CD4","CD8A", "CD8B", top_markers)), assay = "integrated", label = TRUE,
          slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = colors_dark, raster = F)  +
  scale_fill_gradient2(low = "#29A7B8", mid = "white", high = "#E85041", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5)) -> tcell_heatmap
png("heatmap_all_labelled_Tcells_0.6_v1.png", width=16, height=24, units = "in", res = 300)
tcell_heatmap
dev.off()

VlnPlot(tum, "rna_TRGV10")
DimPlot(clus_out, group.by = "integrated_snn_res.0.6", cols = colors_dark)


# other plots ---
p0 = DimPlot(tum, group.by = "integrated_snn_res.0.6", cols = colors_dark, label = T)
p1 = DimPlot(tum, group.by = "sheet_patient_id", cols = colors_dark, label = T)
p2 = DimPlot(tum, group.by = "CD4CD8", cols = colors_darj, label = T)
p3 = FeaturePlot(tum, features = c("rna_CD4", "rna_CD8A"), blend = T)
p = cowplot::plot_grid(p0, p1, p2, p3)
png("DimPlots_all_labelled_Tcells_0.6_v1.png", width=12, height=12, units = "in", res = 300)
p
dev.off()

Idents(tum) <- "integrated_snn_res.0.6"

reduction = "umap"
naive_markers = c("rna_TCF7", "rna_SELL", "rna_LEF1", "rna_CCR7")
inhibitory_markers = c("rna_LAG3", "rna_TIGIT", "rna_PDCD1", "rna_HAVCR2", "rna_CTLA4")
effector_markers = c("rna_IL2", "rna_GZMA", "rna_GZMB", "rna_GZMK", "rna_GNLY", "rna_PRF1", "rna_IFNG", "rna_NKG7", "rna_CX3CR1")
co_stimulatory_markers = c("rna_CD24", "rna_TNFRSF14", "rna_ICOS", "rna_TNFRSF9")
treg_markers = c("rna_IL2RA", "rna_FOXP3", "rna_IKZF2")
tiss_res_mem_markers = c("rna_HAVCR2", "rna_GZMB", "rna_VCAM1", "rna_PRF1", "rna_KLRC1", "rna_CCL4", "rna_LAG3", "rna_CCL3")
tiss_eff_mem_markers = c("rna_KLRG1", "rna_LYAR", "rna_GZMK", "rna_GZMM", "rna_TXNIP", "rna_CD8B", "rna_FCRL6")
nkt = c("rna_KLRD1", "rna_PRF1", "rna_NKG7", "rna_KIR2DL3", "rna_NCR1", "rna_GZMH")
t = c("rna_CD3D", "rna_TRAC", "rna_CD28", "rna_IL7R", "rna_CD8B", "rna_GZMK")
tissue_resident = c("rna_IL2RB", "rna_PRDM1", "rna_ITGAE", "rna_ITGAM", "rna_ITGAX", "rna_CXCR6","rna_CX3CR1", "rna_CD101", "rna_CD69", "rna_ITGA1", "rna_PRDM1")
prol = c("rna_TOP2A", "rna_MKI67")

p1.1 = VlnPlot(tum, features = c("rna_CD3D", "rna_CD4", "rna_CD8A", "rna_CD8B", "rna_ITGAE", "rna_CD69","rna_PRDM1", "rna_GNLY", "rna_NKG7",
                                 "rna_CXCL13", "rna_FOXP3", "rna_CTLA4", "rna_IL7R", "rna_FOS", "rna_RGCC", "rna_GZMB", "rna_PRF1", "rna_STMN1", "rna_TUBB",
                                 "rna_KLRG1", "rna_GZMK", "rna_MKI67", "rna_TOP2A"),pt.size = 0)
p2 = VlnPlot(tum, features = c(t, nkt), pt.size = 0)
p3 = VlnPlot(tum, features = c(naive_markers, inhibitory_markers, effector_markers,co_stimulatory_markers), pt.size = 0)
p4 = VlnPlot(tum, features = c(treg_markers, tiss_res_mem_markers, tiss_eff_mem_markers), pt.size = 0)
p5 = VlnPlot(tum, features = c(tissue_resident), pt.size = 0)
p6 = VlnPlot(tum, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), group.by = "integrated_snn_res.0.6", pt.size = 0, cols = colors_dark)

VlnPlot(tum, features = c("rna_KRT14", "rna_GZMA", "rna_EOMES"), pt.size = 0, group.by = "integrated_snn_res.0.4")
VlnPlot(tum, features = c("rna_GZMB", "rna_GNLY", "rna_NKG7", "rna_IFNG", "rna_PRF1", "rna_CD4", "rna_STAT4", "rna_TBX21"), pt.size = 0)
VlnPlot(tum, features = c("rna_MX1", "rna_ISG15", "rna_IFIT3", "PDCD1"), pt.size = 0)
VlnPlot(tum, features = c("rna_CD4", "rna_TRAC", "rna_TRBC1", "rna_CD3D", "rna_CD3G", "rna_CD3E", "rna_KLRG1"), pt.size = 0)
FeaturePlot(tum, features = c("rna_KRT14"))

Idents(clus_out) <- "integrated_snn_res.0.6"
VlnPlot(clus_out, "rna_CCR7", pt.size = 0, cols = colors_dark)
FeaturePlot(clus_out, "rna_SELL", order = TRUE)  
  # & scale_colour_gradientn(colours = rev(brewer.pal(n = 2, name = "YlGn")))
  DimPlot(clus_out, group.by = "CD4CD8", pt.size = 0, cols = colors_dark)
DimPlot(clus_out, group.by = "integrated_snn_res.0.6", pt.size = 0, cols = colors_dark, label = T)
VlnPlot(clus_out, "nCount_RNA", pt.size = 0, cols = colors_dark)


library(IKAP);library(clustree)
Seurat_obj <- IKAP(clus_out, out.dir = "./IKAP")
clustree(clus_out, prefix = "testclus_", layout = "sugiyama")
clustree(clus_out, prefix = "testclus_", node_colour = "sc3_stability")
clustree(clus_out, prefix = "testclus_",
         node_colour = "CD4", node_colour_aggr = "median")
clus_out@meta.data$testclus_0.4 = clus_out@meta.data$integrated_snn_res.0.4
clus_out@meta.data$testclus_0.6 = clus_out@meta.data$integrated_snn_res.0.6
clus_out@meta.data$testclus_0.2 = clus_out@meta.data$integrated_snn_res.0.2


# single R
# Say this was part of your code for the SingleR run:
clus_out_test = clus_out
DefaultAssay(clus_out_test) <- "RNA"
results <- SingleR(test = as.SingleCellExperiment(clus_out))

# Copy over the labels and pruned.labels (Note: any other column of the results could be used as well)
seurat_obj$SingleR.pruned.calls <- results$pruned.labels
seurat_obj$SingleR.calls <- results$labels

png(paste0(save_path,"/", reduction,dims,"k_", k.param,  "umap_all_labelled_Tcells.png"), width=16, height=16, units = "in", res = 300)
print(p1.1)
dev.off()


png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_all_markers.png"), width=12, height=9, units = "in", res = 300)
print(p2)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_all_markers2.png"), width=12, height=14, units = "in", res = 300)
print(p3)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param,  "umap_all_markers3.png"), width=12, height=12, units = "in", res = 300)
print(p4)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "tissue_resident.png"), width=16, height=12, units = "in", res = 300)
print(p5)
dev.off()

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "QC.png"), width=16, height=12, units = "in", res = 300)
print(p6)
dev.off()

# assign clusters ------
Idents(clus_out) <- "integrated_snn_res.0.6"
mks7_13 = FindMarkers(clus_out, ident.1 = "13", ident.2 = "7")

clus_out$cellstate = "NULL"
clus_out$cellstate[clus_out$integrated_snn_res.0.6 %in% c("0")] = "CD4_Th_like"
clus_out$cellstate[clus_out$integrated_snn_res.0.6 %in% c("1")] = "CD4_Th"
clus_out$cellstate[clus_out$integrated_snn_res.0.6 %in% c("2")] = "CD8_TRM_like"
clus_out$cellstate[clus_out$integrated_snn_res.0.6 %in% c("3")] = "CD4_naive_like"
clus_out$cellstate[clus_out$integrated_snn_res.0.6 %in% c("4")] = "CD4_Teff_memory"
clus_out$cellstate[clus_out$integrated_snn_res.0.6 %in% c("5")] = "NK"
clus_out$cellstate[clus_out$integrated_snn_res.0.6 %in% c("6")] = "NKT"
clus_out$cellstate[clus_out$integrated_snn_res.0.6 %in% c("7")] = "CD8_TRM_like"
clus_out$cellstate[clus_out$integrated_snn_res.0.6 %in% c("8")] = "CD8_Teff_memory"
clus_out$cellstate[clus_out$integrated_snn_res.0.6 %in% c("9")] = "doublets"
clus_out$cellstate[clus_out$integrated_snn_res.0.6 %in% c("10")] = "stressed"
clus_out$cellstate[clus_out$integrated_snn_res.0.6 %in% c("11")] = "CD4-CD8_Tresident"
clus_out$cellstate[clus_out$integrated_snn_res.0.6 %in% c("12")] = "low_quality"
clus_out$cellstate[clus_out$integrated_snn_res.0.6 %in% c("13")] = "CD8_TRM_like"
clus_out$cellstate[clus_out$integrated_snn_res.0.6 %in% c("14")] = "CD8_Teff_memory"
clus_out$cellstate[clus_out$integrated_snn_res.0.6 %in% c("15")] = "NK/ILCs"
clus_out$cellstate[clus_out$integrated_snn_res.0.6 %in% c("16")] = "CD4_treg"
clus_out$cellstate[clus_out$integrated_snn_res.0.6 %in% c("17")] = "CD8_TRM_like"
clus_out$cellstate[clus_out$integrated_snn_res.0.6 %in% c("18")] = "doublets"
clus_out$cellstate[clus_out$integrated_snn_res.0.6 %in% c("19")] = "not_sure"
clus_out$cellstate[clus_out$integrated_snn_res.0.6 %in% c("20")] = "t_prol"


lvl = names(table(clus_out$cellstate))
clus_out$cellstate = factor(clus_out$cellstate, levels = lvl)

write_rds(clus_out, "clus_out.rds")
clus_out = read_rds("clus_out.rds")
Idents(clus_out) <- "cellstate"

clus_out_seu = FindAllMarkers(object = clus_out, assay = "RNA", logfc.threshold = 0)
write.csv(clus_out_seu, "clus_out_seu.csv")

markers_v2 = clus_out_seu %>% dplyr::filter(gene %in% clus_out_seu$gene)

top_markers = markers_v2 %>%
  dplyr::group_by(cluster) %>%
  top_n(7, avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()

DefaultAssay(clus_out) <- "integrated"
p1= DoHeatmap(clus_out, features = c("CD3D", "CD3E", "CD4", "CD8A", "CD8B", top_markers), size = 2,raster = FALSE, label=TRUE,
              group.by = "cellstate", slot = "scale.data", disp.min = -1.5, disp.max = 1.5, group.colors = colors_dark) + 
  scale_color_manual(values = colors_dark) +
  scale_fill_gradient2(low="steelblue1", high = "tomato3", mid = "white", midpoint = 0) + 
  theme(axis.text.y.left = element_text(size=6))

ggsave(plot = p1, filename = "heatmap_cellstate_Tcells_v2.png", device = "png", width = 20, height = 18)

png("heatmap_cellstate_Tcells_v1.png", width=18, height=12, units = "in", res = 300) #
print(p1)
dev.off()

png("umap_Tcells_v1.png", width=9, height=7, units = "in", res = 300) #
DimPlot(clus_out, group.by = "cellstate", cols = col_t, label = T)
dev.off()

# barplots ------
# **Barplots ---------------------------------------------------------------
col_t = c("#C70E7B", "#FC6882", "#007BC3", "#54BCD1", "#EF7C12", "#F4B95A", "#009F3F", "#8FDA04", "#AF6125", "#F4E3C7", "#B25D91", "#EFC7E6",
          "#EF7C12", "#F4B95A", colors_dark)
df_meta = data.frame(clus_out@meta.data, 
                     celltype = clus_out$integrated_snn_res.0.6, 
                     stringsAsFactors = F) %>% 
  as_tibble() %>% clean_names()

df_meta = df_meta %>% dplyr::mutate(sheet_bmi_final = case_when(sheet_bmi_final == "Normal" ~ "normal",
                                                                TRUE ~ as.character(sheet_bmi_final)))

df_meta1 = df_meta
# create a summary table
func_bar(var_df = df_meta1, var_x = "sheet_patient_id", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_patient_id", var_fill = "cellstate", var_pos = "fill", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "final_sample_id", var_fill = "cellstate", var_pos = "stack", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "final_sample_id", var_fill = "cellstate", var_pos = "fill", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_bmi_final", var_fill = "cellstate", var_pos = "stack", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_bmi_final", var_fill = "cellstate", var_pos = "fill", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_brca_updated", var_fill = "cellstate", var_pos = "stack", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_brca_updated", var_fill = "cellstate", var_pos = "fill", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_exp_proc", var_fill = "cellstate", var_pos = "stack", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_exp_proc", var_fill = "cellstate", var_pos = "fill", var_cols = col_t, celltype = "celltype2")
func_bar(var_df = df_meta1, var_x = "sheet_cancer_risk", var_fill = "cellstate", var_pos = "stack", var_cols = col_t, celltype = "ultra_res2")
func_bar(var_df = df_meta1, var_x = "sheet_cancer_risk", var_fill = "cellstate", var_pos = "fill", var_cols = col_t, celltype = "ultra_res2")
func_bar(var_df = df_meta1, var_x = "sheet_figure1_grp_updated", var_fill = "cellstate", var_pos = "stack", var_cols = col_t, celltype = "sheet_figure1_grp_updated")
func_bar(var_df = df_meta1, var_x = "sheet_figure1_grp_updated", var_fill = "cellstate", var_pos = "fill", var_cols = col_t, celltype = "sheet_figure1_grp_updated")
func_bar(var_df = df_meta1, var_x = "sheet_menopause", var_fill = "cellstate", var_pos = "stack", var_cols = col_t, celltype = "sheet_menopause")
func_bar(var_df = df_meta1, var_x = "sheet_menopause", var_fill = "cellstate", var_pos = "fill", var_cols = col_t, celltype = "sheet_menopause")
func_bar(var_df = df_meta1, var_x = "sheet_menopause", var_fill = "cellstate", var_pos = "stack", var_cols = col_t, celltype = "sheet_menopause")
func_bar(var_df = df_meta1, var_x = "sheet_menopause", var_fill = "cellstate", var_pos = "fill", var_cols = col_t, celltype = "sheet_menopause")
func_bar(var_df = df_meta1, var_x = "sheet_cancer_history", var_fill = "cellstate", var_pos = "stack", var_cols = col_t, celltype = "sheet_cancer_history")
func_bar(var_df = df_meta1, var_x = "sheet_cancer_history", var_fill = "cellstate", var_pos = "fill", var_cols = col_t, celltype = "sheet_cancer_history")

func_bar(var_df = df_meta1, var_x = "celltype", var_fill = "sheet_patient_id", var_pos = "fill", var_cols = colors_dark, celltype = "celltype2")


