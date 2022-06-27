# geo ----
setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_lumsec_labeltransfer_v1";dir.create(odir)
options(future.globals.maxSize = 90000 * 1024^2)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "03_subset_lumsec_v2"
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_labeltransfer_v1/"

# core ----
setwd("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_lumsec_labeltransfer_v1"; dir.create(odir)
library(future)
plan("multiprocess", workers = 75)
options(future.globals.maxSize = Inf)
setwd(paste0("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "03_subset_lumsec_20samp_v1"
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_labeltransfer_v1/"
# stress_sig = read.delim("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_fibroblasts_final1/stromal_stress_signature_032720.txt")
# stress_sig = as.character(stress_sig$stromal.stress.sig)


# reference
seu_b2_obj = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_20samp_v2/seu_obj.rds")
Idents(seu_b2_obj) <- "RNA_snn_res.0.3"
DimPlot(seu_b2_obj, reduction = "umap", group.by = "RNA_snn_res.0.3")
FeaturePlot(seu_b2_obj, c("rna_KRT23", "rna_CLDN3", "rna_HMOX1", "rna_KIT", "rna_MAFB", "rna_TOP2A", "rna_CXCL3", "rna_CD74"))
seu_b2_obj$cellstate_tk3 = "NULL"
seu_b2_obj$cellstate_tk3[seu_b2_obj$RNA_snn_res.0.3 %in% c("0")] = "ls_KRT23"
seu_b2_obj$cellstate_tk3[seu_b2_obj$RNA_snn_res.0.3 %in% c("1")] = "ls_CLDN3"
seu_b2_obj$cellstate_tk3[seu_b2_obj$RNA_snn_res.0.3 %in% c("2")] = "ls_KIT"
seu_b2_obj$cellstate_tk3[seu_b2_obj$RNA_snn_res.0.3 %in% c("3")] = "ls_prol"
seu_b2_obj$cellstate_tk3[seu_b2_obj$RNA_snn_res.0.3 %in% c("4")] = "ls_CXCL3"
seu_b2_obj$cellstate_tk3[seu_b2_obj$RNA_snn_res.0.3 %in% c("5")] = "ls_stress"
seu_b2_obj$cellstate_tk3[seu_b2_obj$RNA_snn_res.0.3 %in% c("6")] = "ls_cd74"

write_rds(seu_b2_obj, "seu_b2_obj.rds")
seu_b2_obj = read_rds("seu_b2_obj.rds")

# query
seu_comb_obj = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_v3/dim_35k_20_tum.integrated.ns_RNA_v4_intgeratedassay_integrated_snn_res.0.2_vst_method_9.rds")
FeaturePlot(seu_comb_obj, c("rna_KRT23", "rna_CLDN3", "rna_HMOX1", "rna_KIT", "rna_MAFB", "rna_TOP2A", "rna_CXCL3", "rna_CD74", "rna_LALBA"))
DimPlot(seu_comb_obj, group.by = "cellstate_tk2", cols = colors_dark)
DimPlot(seu_comb_obj, group.by = "cellstate_tk2", cols = colors_dark)

# transfer labels
comb.anchors <- FindTransferAnchors(reference = seu_b2_obj, query = seu_comb_obj, dims = 1:30);write_rds(comb.anchors, "comb.anchors.rds")
comb.prediction <- TransferData(anchorset = comb.anchors, refdata = seu_b2_obj$cellstate_tk3, dims = 1:20);write_rds(comb.prediction, "comb.prediction.rds")
seu_comb_obj <- AddMetaData(seu_comb_obj, metadata = comb.prediction);write_rds(seu_comb_obj, "seu_comb_obj_new.rds")

seu_comb_obj_new = read_rds("seu_comb_obj_new.rds")
DimPlot(seu_comb_obj_new, group.by = "predicted.id", label = T, label.size = 6, cols = colors_dark)

#markers ------

# make heatmap for combined object, original labels 
Idents(seu_b2_obj) <- "cellstate_tk3"
clus_out_seu = FindAllMarkers(object = seu_b2_obj, assay = "RNA", logfc.threshold = 0)
write.csv(clus_out_seu, "clus_out_seu.csv")
clus_out_seu = read.csv("clus_out_seu.csv")
markers_v2 = clus_out_seu %>% dplyr::filter(gene %in% clus_out_seu$gene)

top_markers = markers_v2 %>%
  dplyr::group_by(cluster) %>%
  top_n(12, avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()

DefaultAssay(seu_b2_obj) <- "integrated"
p1= DoHeatmap(seu_b2_obj, features = c(top_markers), size = 2,raster = FALSE, label=TRUE,
              group.by = "cellstate_tk3", slot = "scale.data", disp.min = -1.5, disp.max = 1.5, group.colors = colors_dark) + 
  scale_color_manual(values = colors_dark) +
  scale_fill_gradient2(low="steelblue1", high = "tomato3", mid = "white", midpoint = 0) + 
  theme(axis.text.y.left = element_text(size=6))

tavg = AverageExpression(seu_b2_obj, features = unique(c("EPCAM", "KIT", "TOP2A", top_markers)), return.seurat = T)
DoHeatmap(tavg, features = unique(c("EPCAM", "KIT", "TOP2A", top_markers)), assay = "integrated", label = TRUE,
          slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = colors_dark, raster = F)  +
  scale_fill_gradient2(low = "#29A7B8", mid = "white", high = "#E85041", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5)) -> tcell_heatmap
png("heatmap_all_labelled_lumsec_cellstatetk3_v1.png", width=12, height=12, units = "in", res = 300)
tcell_heatmap
dev.off()

