# geo ----------
setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/")
odir = "figure_4_v1";dir.create(odir)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")

p_load(future)
plan("multiprocess", workers = 40)
options(future.globals.maxSize = 90000 * 1024^2)

data_folder = "figure_4_v1"
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/figure_4_v1/"
setwd(save_path)
col_fibro = c("#5A719A", "#AB5A6B","#E7B356")

# core --------
setwd("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/")
odir = "figure_4_v1";dir.create(odir)
setwd(paste0("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")

p_load(future)
plan("multiprocess", workers = 100)
options(future.globals.maxSize = Inf)

data_folder = "figure_4_v1"
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_4_v1/"
setwd(save_path)
col_fibro = c("#5A719A", "#AB5A6B","#E7B356")

# ************************* -----
# 1. MAIN figure -----
# ** fibro cell figure -----

# **read data ----
fibro = load(file = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_4_v1/hbca_fibros.RObj")
fibro = hbca.fibros
fibro$sheet_patient_id[fibro$sheet_ts_sampleid %in% c("hbca25_c")] = "pt37"

# **annotate clusters ----
fibro$fibro_state = fibro$cell_states
DefaultAssay(fibro) <- 'RNA'
fibro = ScaleData(fibro)


# DefaultAssay(fibro) <- 'integrated'
# fibro = ScaleData(fibro)
# fibro = RunPCA(fibro)
# fibro1 = FindNeighbors(fibro, dims = 1:25, reduction = "pca", k.param = 20)
# fibro1 = FindClusters(fibro1, resolution = c(1, 0.8, 0.6,0.4,0.3,0.2, 0.15,0.17, 0.1, 0.085, 0.090, 0.075, 0.05, 0.023), verbose = TRUE)
# 
# DimPlot(fibro2, group.by = "integrated_snn_res.0.6", cols = colors_dark)
# FeaturePlot(fibro1, c("rna_MMP3", "rna_COL1A1", "rna_FAP", "rna_POSTN", "rna_CXCL14", "rna_APOC1", 
#                       "rna_FN1", "rna_ACTA2", "rna_PDGFRA", "rna_SPRY1", "rna_LOX", "rna_PPARG"), order = T)
# FeaturePlot(fibro1, c("rna_FABP4"), order = T)
# write_rds(fibro1, "fibro1.rds")
# 
# Idents(fibro1) <- "integrated_snn_res.0.8"
# mks = FindAllMarkers(object = fibro1, assay = "RNA", logfc.threshold = 0.1, only.pos = T)
# mks = read.csv("mks.csv")
# mks = mks[, -1]
# top_markers = mks %>% arrange(cluster) %>% 
#   dplyr::group_by(cluster) %>%
#   top_n(7, avg_logFC) %>%
#   dplyr::filter(avg_logFC > 0) %>%
#   pull(gene) %>% as.character()
# 
# d1 = DoHeatmap(fibro1, top_markers, group.colors = colors_dark)
# ggsave(plot = d1, filename = "d1.png", device = "png", dpi = 300, units = "in")
# 

fibro$final_group = "NULL"
fibro$final_group[fibro$fibro_state %in% c("Fibro_COL1A1")] = "Fibro_COL1A1"
fibro$final_group[fibro$fibro_state %in% c("Fibro_PCOLCE2")] = "Fibro_PCOLCE2"
fibro$final_group[fibro$fibro_state %in% c("Fibro_PDPN")] = "Fibro_PDPN"

fibro$final_group = factor(fibro$final_group, levels = c("Fibro_COL1A1", "Fibro_PCOLCE2", "Fibro_PDPN"))


# ^^ Final object -------
write_rds(fibro, "fibro_final_use.rds")

#clus_out = read_rds("fibro_final_use.rds")
fibro = clus_out


DimPlot(fibro, group.by = "final_group", cols = c("#5A719A", "#AB5A6B","#E7B356"), label = T)

# markers ----
Idents(fibro) <- "final_group"
mks = FindAllMarkers(object = fibro, assay = "RNA", logfc.threshold = 0.1, only.pos = T)
write.csv(mks, "fibro_mks.csv")
mks = read.csv("fibro_mks.csv")

top_markers = mks %>% arrange(cluster) %>% 
  dplyr::group_by(cluster) %>%
  slice_max(order_by = avg_logFC, n = 7) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()

DoHeatmap(fibro, top_markers, group.colors = col_fibro)

main_markers = c("COL1A1", "PCOLCE2","PDPN")

DefaultAssay(fibro) <- "integrated"
Idents(fibro) <- "final_group"

tum_sub1 = AverageExpression(fibro, features = unique(top_markers), return.seurat = T, assays = "RNA")
Idents(tum_sub1)

p1 = DoHeatmap(tum_sub1, features = top_markers, assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = col_fibro, raster = F)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
ggsave(plot = p1, filename = "heatmap_cellstate_fibro_allcells_v3.pdf", width = 4, height = 6)

sampled.cells <- WhichCells(fibro, downsample = 100)
sampled.cells2 = fibro[, sampled.cells]
table(sampled.cells2$fibro_state)

fibro = ScaleData(fibro, assay = "RNA")
fibro = ScaleData(fibro, assay = "integrated")
p2 = DoHeatmap(fibro, cells = sampled.cells, features = top_markers, assay = "RNA", label = FALSE,
               slot = "scale.data", disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = col_fibro, raster = F)  +
  scale_fill_gradient2(low = "#29A7B8", mid = "white", high = "#E85041", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
#ggsave(plot = p1, filename = "heatmap_cellstate_myecells_v2.png", device = "png", width = 20, height = 18)

ggsave(plot = p2, filename = "heatmap_cellstate_fibro_100cellsRNA_v1.pdf", width = 8, height = 6)


# ** heatmap chemokines ------
vec_chemokines = rownames(fibro@assays$RNA@scale.data)[str_detect(rownames(fibro@assays$RNA@scale.data), '^CXCL|^CCL')]
vec_chemokines = vec_chemokines[!str_detect(vec_chemokines, 'CCL15-CCL14')]
# vec_krts2 = vec_krts[!vec_krts %in% vec_krts1]
# final_vec_krt = c(vec_krts1,vec_krts2) %>% unique()  

chemo_obj = AverageExpression(fibro, features = vec_chemokines, return.seurat = T, assays = "RNA")

test = chemo_obj@assays$RNA@scale.data
test = cbind(rownames(test), test)
test1 = test[order(test[,2], decreasing = T),]
test2 = test1[order(test1[c(22:40),3], decreasing = T),]
test3 = test[order(test[,4], decreasing = T),]

# ht = Heatmap(chemo_obj@assays$RNA@scale.data)
d = dist(chemo_obj@assays$RNA@scale.data, method = "euclidean")
out = hclust(d, method = "complete")
# plot(out)
# names(sort(cutree(out, h=2.5)))

rownames(chemo_obj)[out$order]

# test0 = data.frame(test)
# test0$gene = rownames(test0)
# 
# test$Fibro_COL1A1 = as.numeric(test$Fibro_COL1A1)
# test$Fibro_PCOLCE2 = as.numeric(test$Fibro_PCOLCE2)
# test$Fibro_PDPN = as.numeric(test$Fibro_PDPN)
# test$V1 = as.character(test$V1)
# 
# test1 = test0 %>% dplyr::arrange(-Fibro_COL1A1) %>% dplyr::filter(Fibro_COL1A1>0, Fibro_PCOLCE2<0, Fibro_PDPN<0) %>% pull(gene)
# test2 = test0 %>% dplyr::filter(!gene %in% test1) %>% arrange(-Fibro_PCOLCE2) %>% dplyr::filter(Fibro_COL1A1<0, Fibro_PCOLCE2>0, Fibro_PDPN<0) %>% pull(gene)
# test3 = test0 %>% dplyr::filter(!gene %in% c(test1, test2)) %>% arrange(-Fibro_PDPN) %>% dplyr::filter(Fibro_COL1A1<0, Fibro_PCOLCE2<0, Fibro_PDPN>0) %>% pull(gene)
# seq_genes = c(test1, test2, test3)
# seq_genes2 = vec_chemokines[!vec_chemokines %in% seq_genes]
# final_genes = c(seq_genes, seq_genes2)

p1 = DoHeatmap(fibro, #cells = sampled.cells, 
               features = rownames(chemo_obj)[out$order], 
               assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = col_fibro, raster = T)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
p1
ggsave(plot = p1, filename = "heatmap_chemokines_celltype_v2.pdf", width = 10, height = 6)


fibro = AddModuleScore(fibro,features = list(vec_chemokines), assay = "RNA", name = "chemokines")
v1=VlnPlot(fibro, features = "chemokines1", pt.size = 0, group.by = "final_group", cols = col_fibro)
ggsave(plot = v1, filename = "vlnplot_chemokines_v1.pdf", width = 4, height = 3)

# ** heatmap MMPs ------
vec_chemokines = rownames(fibro@assays$RNA@counts)[str_detect(rownames(fibro@assays$RNA@counts), '^MMP')]
vec_chemokines = vec_chemokines[!str_detect(vec_chemokines, "MMP25-AS1")]
#vec_krts1 = c(rownames(test1)[1:21], rownames(test2))
# vec_krts2 = vec_krts[!vec_krts %in% vec_krts1]
# final_vec_krt = c(vec_krts1,vec_krts2) %>% unique()  

chemo_obj = AverageExpression(fibro, features = vec_chemokines, return.seurat = T, assays = "RNA")

test = chemo_obj@assays$RNA@scale.data
test = cbind(rownames(test), test)
test1 = test[order(test[,2], decreasing = T),]
test2 = test1[order(test1[c(22:40),3], decreasing = T),]
test3 = test[order(test[,4], decreasing = T),]

# ht = Heatmap(chemo_obj@assays$RNA@scale.data)
d = dist(chemo_obj@assays$RNA@scale.data, method = "euclidean")
out = hclust(d, method = "complete")
plot(out)
names(sort(cutree(out, h=2.5)))

rownames(chemo_obj)[out$order]

save.image("cbstillumap.RData")

# test0 = data.frame(test)
# test0$gene = rownames(test0)
# 
# test$Fibro_COL1A1 = as.numeric(test$Fibro_COL1A1)
# test$Fibro_PCOLCE2 = as.numeric(test$Fibro_PCOLCE2)
# test$Fibro_PDPN = as.numeric(test$Fibro_PDPN)
# test$V1 = as.character(test$V1)
# 
# test1 = test0 %>% dplyr::arrange(-Fibro_COL1A1) %>% dplyr::filter(Fibro_COL1A1>0, Fibro_PCOLCE2<0, Fibro_PDPN<0) %>% pull(gene)
# test2 = test0 %>% dplyr::filter(!gene %in% test1) %>% arrange(-Fibro_PCOLCE2) %>% dplyr::filter(Fibro_COL1A1<0, Fibro_PCOLCE2>0, Fibro_PDPN<0) %>% pull(gene)
# test3 = test0 %>% dplyr::filter(!gene %in% c(test1, test2)) %>% arrange(-Fibro_PDPN) %>% dplyr::filter(Fibro_COL1A1<0, Fibro_PCOLCE2<0, Fibro_PDPN>0) %>% pull(gene)
# seq_genes = c(test1, test2, test3)
# seq_genes2 = vec_chemokines[!vec_chemokines %in% seq_genes]
# final_genes = c(seq_genes, seq_genes2)

p1 = DoHeatmap(chemo_obj, #cells = sampled.cells, 
               features = rownames(chemo_obj)[out$order], 
               assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = col_fibro, raster = T)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
p1
ggsave(plot = p1, filename = "heatmap_mmps_celltype_v3.pdf", width = 4, height = 6)


sampled.cells <- sample(x = Cells(fibro), size = 100, replace = F)

p1 = DoHeatmap(fibro, cells = sampled.cells, 
               features = rownames(chemo_obj)[out$order], 
               assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = col_fibro, raster = T)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
p1
ggsave(plot = p1, filename = "heatmap_MMPs_100cells_celltype_v3.pdf", width = 6, height = 4)


fibro = AddModuleScore(fibro,features = list(vec_chemokines), assay = "RNA", name = "mmps")
v1=VlnPlot(fibro, features = "mmps1", pt.size = 0, group.by = "final_group", cols = col_fibro)
ggsave(plot = v1, filename = "vlnplot_mmps_v1.pdf", width = 4, height = 3)



# ** heatmap Ils ------
vec_chemokines = rownames(fibro@assays$RNA@scale.data)[str_detect(rownames(fibro@assays$RNA@scale.data), '^IL')]
#vec_chemokines = vec_chemokines[!str_detect(vec_chemokines, "MMP25-AS1")]
#vec_krts1 = c(rownames(test1)[1:21], rownames(test2))
# vec_krts2 = vec_krts[!vec_krts %in% vec_krts1]
# final_vec_krt = c(vec_krts1,vec_krts2) %>% unique()  

chemo_obj = AverageExpression(fibro, features = vec_chemokines, return.seurat = T, assays = "RNA")

# test = chemo_obj@assays$RNA@scale.data
# test = cbind(rownames(test), test)
# test1 = test[order(test[,2], decreasing = T),]
# test2 = test1[order(test1[c(22:40),3], decreasing = T),]
# test3 = test[order(test[,4], decreasing = T),]

# ht = Heatmap(chemo_obj@assays$RNA@scale.data)
d = dist(chemo_obj@assays$RNA@scale.data, method = "euclidean")
out = hclust(d, method = "complete")
# plot(out)
# names(sort(cutree(out, h=2.5)))
# 
# rownames(chemo_obj)[out$order]


# test0 = data.frame(test)
# test0$gene = rownames(test0)
# 
# test$Fibro_COL1A1 = as.numeric(test$Fibro_COL1A1)
# test$Fibro_PCOLCE2 = as.numeric(test$Fibro_PCOLCE2)
# test$Fibro_PDPN = as.numeric(test$Fibro_PDPN)
# test$V1 = as.character(test$V1)
# 
# test1 = test0 %>% dplyr::arrange(-Fibro_COL1A1) %>% dplyr::filter(Fibro_COL1A1>0, Fibro_PCOLCE2<0, Fibro_PDPN<0) %>% pull(gene)
# test2 = test0 %>% dplyr::filter(!gene %in% test1) %>% arrange(-Fibro_PCOLCE2) %>% dplyr::filter(Fibro_COL1A1<0, Fibro_PCOLCE2>0, Fibro_PDPN<0) %>% pull(gene)
# test3 = test0 %>% dplyr::filter(!gene %in% c(test1, test2)) %>% arrange(-Fibro_PDPN) %>% dplyr::filter(Fibro_COL1A1<0, Fibro_PCOLCE2<0, Fibro_PDPN>0) %>% pull(gene)
# seq_genes = c(test1, test2, test3)
# seq_genes2 = vec_chemokines[!vec_chemokines %in% seq_genes]
# final_genes = c(seq_genes, seq_genes2)

p1 = DoHeatmap(chemo_obj, #cells = sampled.cells, 
               features = rownames(chemo_obj)[out$order], 
               assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = col_fibro, raster = T)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
p1
ggsave(plot = p1, filename = "heatmap_ilss_celltype_v3.pdf", width = 4, height = 12)


sampled.cells <- sample(x = Cells(fibro), size = 100, replace = F)

p1 = DoHeatmap(fibro, cells = sampled.cells, 
               features = rownames(chemo_obj)[out$order], 
               assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = col_fibro, raster = T)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
p1
ggsave(plot = p1, filename = "heatmap_ilss_100cells_celltype_v3.pdf", width = 6, height = 12)


fibro = AddModuleScore(fibro,features = list(vec_chemokines), assay = "RNA", name = "ils")
v1=VlnPlot(fibro, features = "ils1", pt.size = 0, group.by = "final_group", cols = col_fibro)
ggsave(plot = v1, filename = "vlnplot_ils_v1.pdf", width = 4, height = 3)


# ** heatmap collagen ------
fibro = clus_out
fibro = ScaleData(fibro, assay = "RNA")
Idents(fibro) <- "fibro_state"
vec_chemokines = rownames(fibro@assays$RNA@scale.data)[str_detect(rownames(fibro@assays$RNA@scale.data), '^COL')]
#vec_chemokines = vec_chemokines[!str_detect(vec_chemokines, "MMP25-AS1")]
#vec_krts1 = c(rownames(test1)[1:21], rownames(test2))
# vec_krts2 = vec_krts[!vec_krts %in% vec_krts1]
# final_vec_krt = c(vec_krts1,vec_krts2) %>% unique()  

chemo_obj = AverageExpression(fibro, features = vec_chemokines, return.seurat = T, assays = "RNA")

test = chemo_obj@assays$RNA@scale.data
test = cbind(rownames(test), test)
test1 = test[order(test[,2], decreasing = T),]
test2 = test1[order(test1[c(22:40),3], decreasing = T),]
test3 = test[order(test[,4], decreasing = T),]

# ht = Heatmap(chemo_obj@assays$RNA@scale.data)
d = dist(chemo_obj@assays$RNA@scale.data, method = "euclidean")
out = hclust(d, method = "complete")
plot(out)
names(sort(cutree(out, h=2.5)))

rownames(chemo_obj)[out$order]

save.image("cbstillumap.RData")

# test0 = data.frame(test)
# test0$gene = rownames(test0)
# 
# test$Fibro_COL1A1 = as.numeric(test$Fibro_COL1A1)
# test$Fibro_PCOLCE2 = as.numeric(test$Fibro_PCOLCE2)
# test$Fibro_PDPN = as.numeric(test$Fibro_PDPN)
# test$V1 = as.character(test$V1)
# 
# test1 = test0 %>% dplyr::arrange(-Fibro_COL1A1) %>% dplyr::filter(Fibro_COL1A1>0, Fibro_PCOLCE2<0, Fibro_PDPN<0) %>% pull(gene)
# test2 = test0 %>% dplyr::filter(!gene %in% test1) %>% arrange(-Fibro_PCOLCE2) %>% dplyr::filter(Fibro_COL1A1<0, Fibro_PCOLCE2>0, Fibro_PDPN<0) %>% pull(gene)
# test3 = test0 %>% dplyr::filter(!gene %in% c(test1, test2)) %>% arrange(-Fibro_PDPN) %>% dplyr::filter(Fibro_COL1A1<0, Fibro_PCOLCE2<0, Fibro_PDPN>0) %>% pull(gene)
# seq_genes = c(test1, test2, test3)
# seq_genes2 = vec_chemokines[!vec_chemokines %in% seq_genes]
# final_genes = c(seq_genes, seq_genes2)

p1 = DoHeatmap(chemo_obj, #cells = sampled.cells, 
               features = rownames(chemo_obj)[out$order], 
               assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = col_fibro, raster = T)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
p1
ggsave(plot = p1, filename = "heatmap_collagens_celltype_v3.pdf", width = 4, height = 6)

sampled.cells <- sample(x = Cells(fibro), size = 100, replace = F)

p1 = DoHeatmap(fibro, cells = sampled.cells, 
               features = rownames(chemo_obj)[out$order], 
               assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = col_fibro, raster = T)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
p1
ggsave(plot = p1, filename = "heatmap_collagens_100cells_celltype_v3.pdf", width = 6, height = 4)

fibro = AddModuleScore(fibro,features = list(vec_chemokines), assay = "RNA", name = "collagens")
v1=VlnPlot(fibro, features = "collagens1", pt.size = 0, group.by = "final_group", cols = col_fibro)
ggsave(plot = v1, filename = "vlnplot_collagens_v1.pdf", width = 4, height = 3)
FeaturePlot(fibro, features = c("collagens1", "mmps1", "chemokines1"))
DimPlot(fibro, label = T)
# **Umap ------

umap_tx_fibro = fibro@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(tum_meta = fibro@meta.data)

umap_tx_fibro  =  umap_tx_fibro %>% mutate(major_group = factor(tum_meta.final_group, levels = c("Fibro_COL1A1", "Fibro_PCOLCE2", "Fibro_PDPN")))

fill_col = col_fibro
library(ggpubr);library(grid)
# grob <- grobTree(textGrob("113157 cells", x=0.7,  y=0.95, hjust=0,
#                           gp=gpar(col="black", fontsize=10, fontface="italic")))
umap_tx_fibro1 <- umap_tx_basal[sample(x = 1:nrow(x = umap_tx_basal)), ]
p1 = ggplot(umap_tx_fibro, aes(x=UMAP_1, y=UMAP_2)) + #annotation_custom(grob) + 
  geom_jitter(shape = 21, size = 1, aes(fill = major_group),  color = "black", alpha = 0.9, stroke = 0.04) + 
  scale_fill_manual(values = col_fibro) +
  xlab("UMAP-1") + ylab("UMAP-2") +
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(),
    legend.position = "none",
    axis.line.x = element_line(colour = "black"), #,arrow=arrow
    axis.line.y = element_line(colour = "black"),#,arrow=arrow
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),axis.text.y = element_blank());p1
p_fixed <- set_panel_size(p1,
                          width  = unit(1.5, "in"),
                          height = unit(1.5, "in"))

ggsave(plot = p_fixed, filename = "fibro_major_umap.png", device = "png", path = save_path, width = 3, height = 3, dpi = 300, units = "in")
ggsave(plot = p_fixed, filename = "fibro_major_umap.pdf", device = "pdf", path = save_path, width = 3, height = 3, dpi = 300, units = "in")

p03 = DimPlot(object = fibro, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = col_fibro, group.by = "final_group", shuffle = T)+  NoAxes()
pdf(glue("{save_path}/fibro_major_uamp_1.pdf"), width=6, height=5)
print(plot_grid(p03))
dev.off()

write_rds(fibro, "fibro_final.rds")
#basal = read_rds("basal_final.rds")

# **Barplots ---------------------------------------------------------------
fill_col = col_fibro
tum = fibro
df_meta = data.frame(tum@meta.data, 
                     celltype = tum$final_group, 
                     stringsAsFactors = F) %>% 
  as_tibble() %>% clean_names()

df_meta1 = df_meta
# create a summary table
func_bar(var_df = df_meta1, var_x = "sheet_patient_id", var_fill = "celltype", var_pos = "fill", var_cols = fill_col, celltype = "fibro_cells")
func_bar(var_df = df_meta1, var_x = "sheet_figure1_grp_updated", var_fill = "celltype", var_pos = "fill", var_cols = fill_col, celltype = "fibro_cells")

df_meta2 = df_meta1 %>% dplyr::filter(sheet_tissue_location %in% c("left", "right"), !sheet_patient_id == "pt45")
df_meta3 = df_meta2 %>% mutate(sheet_patient_id =  factor(sheet_patient_id, levels = unique(sheet_patient_id[order(sheet_patient_id, sheet_tissue_location)]), ordered = TRUE))
df_meta3 = df_meta3 %>% mutate(key = paste0(sheet_patient_id, "_", sheet_tissue_location))

p1 = ggplot(df_meta3, aes(sheet_tissue_location, fill = celltype)) +
  geom_bar(position = "fill", width = 1) + 
  xlab("") + ylab("") + scale_fill_manual(values = fill_col) +
  theme(axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
        panel.background = element_blank(),axis.title.x = element_blank(), 
        axis.text.x = element_blank(), axis.ticks = element_blank(), 
        legend.title = element_blank()) + facet_wrap(~sheet_patient_id, ncol = 22, strip.position = "bottom") + 
  theme(strip.background = element_blank(), panel.spacing.x=unit(0.2, "lines"), strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(2,0,0,0, "pt")))
ggsave(p1, file = paste0(save_path, "/figure4_b_tissue_location_patient.pdf"), width = 7, height = 2.5)

func_bar(var_df = df_meta3, var_x = "key", var_fill = "celltype", var_pos = "fill", var_cols = fill_col, celltype = "fibro_cells_leftright")
func_bar(var_df = df_meta3, var_x = "sheet_tissue_location", var_fill = "celltype", var_pos = "fill", var_cols = fill_col, celltype = "fibro_cells_leftright")

# df_meta2 = df_meta1 %>% dplyr::filter(sheet_bmi_final %in% c("normal","Normal","obese", "overweight"))
# df_meta2 = df_meta2 %>% dplyr::mutate(bmi = case_when(sheet_bmi_final %in% c("normal", "Normal") ~ "normal",
#                                                       sheet_bmi_final %in% c("obese") ~ "obese",
#                                                       sheet_bmi_final %in% c("overweight") ~ "overweight",
#                                                       sheet_bmi_final %in% c("unknown") ~ "unknown",
#                                                       is.na(sheet_bmi_final)  ~ "unknown"))
# df_meta2$bmi = factor(df_meta2$bmi, levels = c("normal", "overweight", "obese"), ordered = T) 
# df_meta3 = df_meta2 %>% mutate(sheet_patient_id =  factor(sheet_patient_id, levels = unique(sheet_patient_id[order(sheet_patient_id, bmi)]), ordered = TRUE))
# df_meta3 = df_meta3 %>% mutate(key = paste0(sheet_patient_id, "_", bmi))
# 
# func_bar(var_df = df_meta3, var_x = "bmi", var_fill = "celltype", var_pos = "fill", var_cols = col_grp, celltype = "b_cells_bmi_average")
# 
# df_meta2 = df_meta2 %>% dplyr::filter(!sheet_figure1_grp_updated %in% c("unknown"), !is.na(sheet_figure1_grp_updated))
# df_meta2$sheet_figure1_grp_updated = factor(df_meta2$sheet_figure1_grp_updated, levels = c("Reduction Mammoplasty", "Prophylatic Mastectomy", "Cancer Mastectomy"), ordered = T)
# df_meta3 = df_meta2 %>% mutate(sheet_patient_id =  factor(sheet_patient_id, levels = unique(sheet_patient_id[order(sheet_patient_id, sheet_figure1_grp_updated)]), ordered = TRUE))
# df_meta3 = df_meta3 %>% mutate(key = paste0(sheet_patient_id, "_", sheet_figure1_grp_updated))
# 
# func_bar(var_df = df_meta3, var_x = "sheet_figure1_grp_updated", var_fill = "celltype", var_pos = "fill", var_cols = col_grp, celltype = "b_cells_sheet_figure1_grp_updated_average")

# **Clinical ---------
df_meta2_sel = df_meta1 %>% dplyr::select(final_group1, final_group, final_sample_id, ethnicity, age, sheet_patient_id, sheet_source, sheet_exp_proc, sheet_age_updated, sheet_brca_updated, sheet_parity, sheet_gravida, sheet_menopause, sheet_density, sheet_bmi_final, sheet_cancer_risk,sheet_figure1_grp_updated)
colnames(df_meta2_sel) = c("cellgroup", "celltype", "sample_id", "ethnicity", "age_orig", "patient_id",  "source", "exp_proc", "age", "brca", "parity", "gravida", "menopause", "density", "bmi", "cancer_risk", "surgery")
# df_mrg = df_meta2_sel %>% mutate(pt_x = paste0(patient_id, "_", type))
write_rds(df_meta2_sel, "fibro_df_cell_meta2_sel.rds")
df_meta2_sel = read_rds("fibro_df_cell_meta2_sel.rds")
# annotate cell group, age
df_meta2_sel = df_meta2_sel %>% dplyr::mutate(age2 = case_when(age_orig >= 50 ~ "old",
                                                               age_orig < 50 ~ "young",
                                                               age_orig ==  "unknown" ~ "unknown",
                                                               is.na(age_orig)  ~ "unknown"))

df_meta2_sel = df_meta2_sel %>% dplyr::mutate(density2 = case_when(density %in% c("Density_1") ~ "low",
                                                                   density %in% c("Density_2", "Density_3") ~ "high",
                                                                   density ==  "unknown" ~ "unknown",
                                                                   is.na(density)  ~ "unknown"))

df_meta2_sel = df_meta2_sel %>% dplyr::mutate(bmi = case_when(bmi %in% c("normal", "Normal") ~ "normal",
                                                              bmi %in% c("obese") ~ "obese",
                                                              bmi %in% c("overweight") ~ "overweight",
                                                              bmi %in% c("unknown") ~ "unknown",
                                                              is.na(bmi)  ~ "unknown"))

df_meta2_sel = df_meta2_sel %>% dplyr::mutate(bmi2 = case_when(bmi %in% c("normal", "Normal") ~ "normal",
                                                               bmi %in% c("obese") ~ "overweight",
                                                               bmi %in% c("overweight") ~ "overweight",
                                                               bmi %in% c("unknown") ~ "unknown",
                                                               is.na(bmi)  ~ "unknown"))

df_meta2_sel = df_meta2_sel %>% dplyr::mutate(exp_proc_updated = case_when(exp_proc %in% c("2hr_digstion","3hr_digestion","4hr_digestion","6hr_digestion", "short_digestion") ~ "short",
                                                                           exp_proc %in% c("overnight_digestion") ~ "overnight",
                                                                           exp_proc %in% c("unknown") ~ "unknown",
                                                                           is.na(exp_proc)  ~ "unknown"))

df_meta2_sel = df_meta2_sel %>% dplyr::mutate(germline_cancer_risk = case_when(brca %in% c("brca_neg") ~ "normal",
                                                                               brca %in% c("brca_pos", "hereditary") ~ "high",
                                                                               brca %in% c("unknown") ~ "unknown",
                                                                               is.na(brca)  ~ "unknown"))

df_meta2_sel = df_meta2_sel %>% dplyr::mutate(ethnicity2 = case_when(ethnicity %in% c("black") ~ "african-american",
                                                                     ethnicity %in% c("hispanic") ~ "latino",
                                                                     ethnicity %in% c("white") ~ "caucasian",
                                                                     ethnicity %in% c("unknown") ~ "unknown",
                                                                     is.na(ethnicity)  ~ "unknown"))

df_meta2_sel[is.na(df_meta2_sel)] = "unknown"
table(is.na(df_meta2_sel))

df_mrg_nuc1 = df_meta2_sel %>% group_by(sample_id) %>% mutate(total_cells_per_samp = n()) %>% ungroup() %>% 
  group_by(sample_id, celltype) %>% mutate(cell_per_samp = n(), prop_cell_per_samp = cell_per_samp/total_cells_per_samp) %>% ungroup() %>% 
  group_by(sample_id, cellgroup) %>% mutate(cell_per_grp_samp = n(), prop_cell_per_grp_samp = cell_per_samp/cell_per_grp_samp) %>% ungroup()

df_mrg_nuc2 = df_meta2_sel %>% group_by(patient_id) %>% mutate(total_cells_per_pt = n()) %>% ungroup() %>% 
  group_by(patient_id, celltype) %>% mutate(cell_per_pt = n(), prop_cell_per_pt = cell_per_pt/total_cells_per_pt) %>% ungroup() %>% 
  group_by(patient_id, cellgroup) %>% mutate(cell_per_grp_pt = n(), prop_cell_per_grp_pt = cell_per_pt/cell_per_grp_pt) %>% ungroup()

df_mrg_nuc2$celltype = factor(df_mrg_nuc2$celltype, levels = c("Fibro_COL1A1", "Fibro_PCOLCE2", "Fibro_PDPN"))

write.csv(df_mrg_nuc1, "df_mrg_nuc1_samp_fibro.csv")
write.csv(df_mrg_nuc2, "df_mrg_nuc2_pt_fibro.csv")

save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_4_v1/fibro_all/"
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/figure_4_v1/fibro_all/"
col_pal = fill_col

do_clinical_comp(df_mrg_nuc2 = df_mrg_nuc2)

df_meta2_sel_short =  df_meta2_sel %>% dplyr::filter(exp_proc_updated == "short")
df_meta2_sel_over =  df_meta2_sel %>% dplyr::filter(exp_proc_updated == "overnight")

df_mrg_nuc2_short = df_meta2_sel_short %>% group_by(patient_id) %>% mutate(total_cells_per_pt = n()) %>% ungroup() %>% 
  group_by(patient_id, celltype) %>% mutate(cell_per_pt = n(), prop_cell_per_pt = cell_per_pt/total_cells_per_pt) %>% ungroup() %>% 
  group_by(patient_id, cellgroup) %>% mutate(cell_per_grp_pt = n(), prop_cell_per_grp_pt = cell_per_pt/cell_per_grp_pt) %>% ungroup()
df_mrg_nuc2_over = df_meta2_sel_over %>% group_by(patient_id) %>% mutate(total_cells_per_pt = n()) %>% ungroup() %>% 
  group_by(patient_id, celltype) %>% mutate(cell_per_pt = n(), prop_cell_per_pt = cell_per_pt/total_cells_per_pt) %>% ungroup() %>% 
  group_by(patient_id, cellgroup) %>% mutate(cell_per_grp_pt = n(), prop_cell_per_grp_pt = cell_per_pt/cell_per_grp_pt) %>% ungroup()

save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_4_v1/fibro_short/"
do_clinical_comp(df_mrg_nuc2 = df_mrg_nuc2_short)
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_4_v1/fibro_over/"
do_clinical_comp(df_mrg_nuc2 = df_mrg_nuc2_over)



# cell cycle analysis Cell---------
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
DefaultAssay(tum) <- "RNA"
clus_out <- CellCycleScoring(tum, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

RidgePlot(clus_out, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2, group.by = "final_group", cols = fill_col)
DimPlot(clus_out, group.by = "Phase", cols = colors_bottle2)
DimPlot(clus_out, group.by = "final_group", cols = fill_col)
FeaturePlot(clus_out, features = c("S.Score", "G2M.Score"))


v1 = VlnPlot(clus_out, features = c("S.Score", "G2M.Score"), group.by = "final_group", cols = fill_col, pt.size = 0)
ggsave(v1, file =  "cellcycle_cellstate.pdf", width = 4, height = 3)

# pathway ------
Idents(clus_out) <- "final_group"
all_fibro = FindAllMarkers(object = clus_out, assay = "RNA", logfc.threshold = 0, only.pos = FALSE)
write.csv(all_fibro, paste0(save_path, "fibro_markers.csv"))

genes = read.csv("fibro_markers.csv")
genes$cluster = as.character(genes$cluster)

gsea_resultsc2hall = fun_gsea(genes = genes, set = "H", sub_cat = NULL, num_avg_logFC = 0.25)
gsea_resultsc2C5 = fun_gsea(genes = genes, set = "C5", sub_cat = "BP", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "lumsec")
gsea_resultsc2C3 = fun_gsea(genes = genes, set = "C3", sub_cat = "TFT:GTRD", num_avg_logFC = 0.25, celltype = "lumsec", doplot = FALSE, ret_collap = FALSE)
gsea_resultsc2C8 = fun_gsea(genes = genes, set = "C8", sub_cat = NULL, num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "lumsec")
gsea_resultsc2C2_reac = fun_gsea(genes = genes, set = "C2", sub_cat = "CP:REACTOME", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "lumsec")
gsea_resultsc2C2_kegg = fun_gsea(genes = genes, set = "C2", sub_cat = "CP:KEGG", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "fibro")

clus = c("Fibro_COL1A1", "Fibro_PCOLCE2", "Fibro_PDPN")

gs = gsea_resultsc2hall
gs = gsea_resultsc2C2;gs_kegg = gsea_resultsc2C2_kegg;gs_reac = gsea_resultsc2C2_reac
gs = gs_kegg
gs = gs_reac
gs = gsea_resultsc2C5
gs = gsea_resultsc2C3
gs = gsea_resultsc2C7
gs = gsea_resultsc2C8

gs = rbind(gsea_resultsc2C5, gsea_resultsc2C2_kegg, gsea_resultsc2C2_reac)
#gs = gs %>% dplyr::filter(cluster %in% c("lumhr_egln3","lumhr_fasn","lumhr_pip"))

gs1 = gs %>% dplyr::select(-leadingEdge) %>% dplyr::group_by(cluster) %>% 
  dplyr::filter(padj <= 0.05) %>%
  top_n(20, wt = NES) %>% 
  pull(pathway)
gs2 = gs %>% dplyr::filter(pathway %in% gs1)
gs2$cluster = factor(gs2$cluster, levels = clus, ordered = TRUE)
gs2 = gs2 %>% dplyr::filter(grepl('HALLMARK|KEGG|REACTOME|BIOCARTA|GO', pathway))
#gs2 = gs2 %>% dplyr::filter(grepl('EPITHELIAL|ESTROGEN|MAMMARY|IMMUNE|PROGESTRONE|PROLIFERATION|STEM|HORMONE|LACTATION|MILK|SECRETION|SECRETORY', pathway))

#gs2 = gs2 %>% mutate(pathway = str_split(pathway, "_", n=2, simplify = TRUE)[,2]) # for C7
gs2 = gs2 %>%  arrange(cluster, -NES)

library(forcats);library(circlize);library(ComplexHeatmap)

p1 = gs2 %>% 
  dplyr::select(pathway, NES, cluster) %>%
  arrange(cluster, NES) %>%
  spread(key = cluster, value = NES, fill=0) 
#anova(p1[,-1])

p1.1 = p1
colnames(p1.1) = c("pathway", "Fibro_COL1A1", "Fibro_PCOLCE2","Fibro_PDPN" )

# rearrange the datafram
test = p1

# for basal
test1 = test %>% dplyr::arrange(-`Fibro_COL1A1`) %>% dplyr::filter(Fibro_COL1A1>0, Fibro_PCOLCE2<0, Fibro_PDPN<0) %>% pull(pathway)
test2 = test %>% dplyr::filter(!pathway %in% test1) %>% arrange(-Fibro_PCOLCE2) %>% dplyr::filter(Fibro_COL1A1<0, Fibro_PCOLCE2>0, Fibro_PDPN<0) %>% dplyr::top_n(25, Fibro_PCOLCE2) %>% pull(pathway)
test3 = test %>% dplyr::filter(!pathway %in% c(test1, test2)) %>% arrange(-Fibro_PDPN) %>% dplyr::filter(Fibro_COL1A1<0, Fibro_PCOLCE2<0, Fibro_PDPN>0) %>% dplyr::top_n(25, Fibro_PDPN) %>% pull(pathway)
seq_path = c(test1, test2, test3) #

#seq_path = test$pathway
# make heatmap 
p1.1 = p1.1 %>%
  dplyr::slice(match(seq_path, pathway))
p1 = p1.1
p2 = data.matrix(p1)
rownames(p2) = p1$pathway
p2 = p2[,-1]

# paths = ((abs(p2) > 1) * 1) %>% rowMax()
# p2.2 = p2[paths == 1, ]

col = colorRamp2(c(-2, 0, 2),c("steelblue",  "white",  "#F80F61")) 
# c("steelblue","white", "tomato3"))
ht_hall = Heatmap(p2, 
                  show_row_names = T, "NES", 
                  col = col,
                  #km = 3,
                  #clustering_method_rows = "complete",
                  #clustering_distance_rows = "maximum",
                  row_names_gp = gpar(cex = 0.6),
                  cluster_rows = T,
                  cluster_columns = F,
                  #column_dend_side = "top",
                  #show_column_dend = T,
                  row_dend_side = "right",
                  row_names_side = "left",
                  width = unit(5, "cm"), 
                  row_names_max_width = unit(20, "cm"),
                  column_names_gp = gpar(cex = 0.8));ht_hall


png("fibro_htmaps_pathways_top20padj0.05_v1.png", res = 300, units = "in", width = 10, height = 14)
print(ht_hall)
dev.off()

# cluster profiler -----

#ClusterProfiler
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationHub")
library("clusterProfiler")
library("org.Hs.eg.db")
library("AnnotationHub")

top100 <- mks %>% group_by(cluster) %>% dplyr::filter(avg_logFC>0.25) %>% top_n(n = 200, wt = avg_logFC)
top100pval <- subset(top100, rowSums(top100[6] < 0.05) > 0)

df <- top100pval[,7:8]
dfsample <- split(df$gene,df$cluster)
length(dfsample)

#The output of length(dfsample) returns how many clusters you have
#I'm sure there's a better way but you have to make a line like below for each cluster

dfsample$Fibro_COL1A1 = bitr(dfsample$Fibro_COL1A1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$Fibro_PCOLCE2 = bitr(dfsample$Fibro_PCOLCE2, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$Fibro_PDPN = bitr(dfsample$Fibro_PDPN, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")


#do the same here, a line like below for each cluster
genelist <- list("Fibro_COL1A1" = dfsample$Fibro_COL1A1$ENTREZID, 
                 "Fibro_PCOLCE2" = dfsample$Fibro_PCOLCE2$ENTREZID,
                 "Fibro_PDPN" = dfsample$Fibro_PDPN$ENTREZID)
GOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db")
pdf("fibro_GOclusterplot_v1_200.pdf", width = 10, height = 5)
dotplot(GOclusterplot, showCategory = 10)
dev.off()

KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG")
pdf("fibro_keggclusterplot_v1_200.pdf", width = 10, height = 5)
dotplot(KEGGclusterplot , showCategory = 10)
dev.off()

KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichDO")
dotplot(KEGGclusterplot , showCategory = 10)

# ** signatures -------

pre_adipo = read.table("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_4_v1/fibroblast_sigs/fibro_nuc_preadipo_sig.txt") %>% dplyr::filter(avg_logFC>0.75) %>% pull(gene) %>% as.character()
pre_adipo[19] = "MARCH1"

# plot signatures
pre_adipo_lst = list(pre_adipo)
clus_out <- AddModuleScore(clus_out, features = pre_adipo_lst, assay = "RNA", name = "pre_adipo_lst_1", nbin = 20)
v3 = VlnPlot(clus_out, features = c("pre_adipo_lst_11"), group.by = "fibro_state", cols = col_fibro, pt.size = 0)
ggsave(v3, file =  "pre_adipo_score.pdf", width = 4, height = 3)

ls.files = list.files(full.names = T, path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_4_v1/fibroblast_sigs/",
                      pattern = ".txt", )
for (i in 4:length(ls.files)) {
  
  
  file_input = read.table(file = ls.files[i+1]) %>% pull(V1)
  file_input_lst = list(file_input)
  bs.name = basename(ls.files[i+1])
  bs.name = str_remove(bs.name, pattern = ".txt")
  
  message("working on - ", bs.name)
  clus_out <- AddModuleScore(clus_out, features = file_input_lst, assay = "RNA", name = bs.name, nbin = 20)
  v3 = VlnPlot(clus_out, features = paste0(bs.name, "1"), group.by = "fibro_state", cols = col_fibro, pt.size = 0)
  ggsave(v3, file =  paste0(bs.name, ".pdf"), width = 4, height = 3)
  
}

g1 = FeaturePlot(fibro, features = c("rna_ACTA2", "rna_TAGLN", "rna_MYL9"), ncol = 3, combine = T)
g2 = g1 + theme_cowplot()
ggsave(plot = g2, filename = "featureplot_myofibro.png", width = 7, height = 2)

to_plots <- c("rna_COL1A1", "rna_WISP2")
DefaultAssay(clus_out) <- "RNA"
plots <- lapply(to_plots, function(to_plot) {               
  plt <- FeaturePlot(clus_out, to_plot, order = TRUE) + 
    ggtitle(to_plot) +
    theme(title = element_text(size = 7)) +
    scale_color_viridis_c(direction = -1) +
    theme_cowplot() & NoLegend() & NoAxes()
  
  return(plt)
})

names(plots) <- to_plots
#plots_flat <- do.call(, plots)

# Stitch all the plots into 1.
plt <- patchwork::wrap_plots(plots, ncol=1)

# Save
ggsave("feature_MMP3.pdf", plot = plt, width = 2, height = 4, limitsize = FALSE)  



# ************************* -----
# 1. SUPP figure -----
# ** Full heatmap -----------
# supp mean heatmap for fig 1
DefaultAssay(clus_out) <- "RNA"
clus_out = ScaleData(clus_out)
all_markers = read.csv("./fibro_mks.csv")
all_markers = all_markers[,-1]

all_markers_sel2 = all_markers
lvl = levels(clus_out$final_group)
all_markers_sel2$cluster = factor(all_markers_sel2$cluster, levels = lvl)


top_basal_markers = all_markers %>% 
  group_by(cluster) %>%
  filter(!str_detect(gene, "^RP|^HSP")) %>% 
  dplyr::slice_max(n = 7, order_by = avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()

genes_sel_df = all_markers_sel2 %>% arrange(cluster) %>% 
  group_by(cluster) %>%
  filter(!str_detect(gene, "^RP|^HSP")) %>% 
  dplyr::slice_max(n = 7, order_by = avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  dplyr::select(gene, cluster) %>% mutate(id = paste0(gene, "-", cluster))

main_markers = c("POSTN", "PCOLCE2", "CXCL1")

all_markers_sel_fig = data.frame(gene = c(top_basal_markers), cluster1 = c(rep("Fibro_COL1A1", 7),rep("Fibro_PCOLCE2", 7),rep("Fibro_PDPN", 7)))
DefaultAssay(clus_out) <- "integrated"
clus_out = ScaleData(clus_out)
Idents(clus_out) <- "final_group"
all_markers_sel_fig_genes = unique(c(main_markers, top_basal_markers))
col_grp = col_fibro
names(col_grp) = lvl
df_fibro1 = run_mean_htmaps_fig1_supp(df = all_markers_sel_fig, seu = clus_out, cell = "fibro_cells", genes_df = all_markers_sel_fig_genes, #sel_genes
                                      genes_df_cann = main_markers, genes_sel_df = genes_sel_df, 
                                      clus_levels = lvl)


v1 = VlnPlot(clus_out, features = c("FAP"), pt.size = 0, assay = "RNA", group.by = "final_group", cols = fill_col)
ggsave(v1, file =  "vlnplot_progenitor_genes.pdf", width = 8, height = 7)
v1 = FeaturePlot(clus_out, features = c("rna_FAP"))
ggsave(v1, file =  "vlnplot_progenitor_genes.pdf", width = 8, height = 7)

# ** > fibro ----
to_plots <- c("COL1A1", "MMP3", "LUM")
to_plots <- c("FBLN1", "SERPINF1", "COL1A2")
DefaultAssay(clus_out) <- "RNA"
plots <- lapply(to_plots, function(to_plot) {               
  plt <- FeaturePlot(clus_out, to_plot, order = T, raster = F) + 
    ggtitle(to_plot) +
    theme(title = element_text(size = 7)) +
    scale_color_viridis_c(direction = -1) +
    theme_cowplot() & NoLegend() & NoAxes()
  
  return(plt)
})

names(plots) <- to_plots
#plots_flat <- do.call(, plots)
plots[]
# Stitch all the plots into 1.
plt <- patchwork::wrap_plots(plots, ncol=1)

# Save
ggsave("feature_fibro_canonical.pdf", plot = plt, width = 1.5, height = 6, limitsize = FALSE)  
ggsave("feature_fibro_canonical21.pdf", plot = plt, width = 2, height = 6, limitsize = FALSE)  
ggsave("feature_fibro_canonical31.pdf", plot = plt, width = 2, height = 6, limitsize = FALSE)  
ggsave("feature_fibro_canonical41.pdf", plot = plt, width = 2, height = 6, limitsize = FALSE)  



# *********************** ----------
# 3. ADIPOCYTES ---------
col_adipo = "ivory"
# **Umap ------
adipo = load("hbca_adipo_nucs_final.RObj")
adipo = adipo.nucs.filtered.integrated
umap_tx_adipo = adipo@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(tum_meta = adipo@meta.data)

umap_tx_adipo  =  umap_tx_adipo %>% mutate(major_group = factor(tum_meta.cell_states, levels = c("Adipo")))
dim(umap_tx_adipo) #4197

fill_col = col_adipo
library(ggpubr);library(grid)
umap_tx_adipo1 <- umap_tx_adipo[sample(x = 1:nrow(x = umap_tx_adipo)), ]
p1 = ggplot(umap_tx_adipo, aes(x=UMAP_1, y=UMAP_2)) + #annotation_custom(grob) + 
  geom_jitter(shape = 21, size = 1, aes(fill = major_group),  color = "black", alpha = 0.9, stroke = 0.04) + 
  scale_fill_manual(values = "#fcf7e3") + #cornsilk ivory2 f9f7f2
  xlab("UMAP-1") + ylab("UMAP-2") +
  theme(
    panel.grid = element_blank(), 
    panel.background = element_blank(),
    legend.position = "none",
    axis.line.x = element_line(colour = "black"), #,arrow=arrow
    axis.line.y = element_line(colour = "black"),#,arrow=arrow
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),axis.text.y = element_blank());p1
p_fixed <- set_panel_size(p1,
                          width  = unit(1.5, "in"),
                          height = unit(1.5, "in"))

ggsave(plot = p_fixed, filename = "adipo_major_umap.png", device = "png", path = save_path, width = 3, height = 3, dpi = 300, units = "in")
ggsave(plot = p_fixed, filename = "adipo_major_umap.pdf", device = "pdf", path = save_path, width = 3, height = 3, dpi = 300, units = "in")

p03 = DimPlot(object = adipo, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = col_adipo, group.by = "cell_states", shuffle = T)+  NoAxes()
pdf(glue("{save_path}/adipo_major_uamp_1.pdf"), width=6, height=5)
print(plot_grid(p03))
dev.off()
