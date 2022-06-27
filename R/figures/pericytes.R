# geo ----------
setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/")
odir = "figure_5_v1";dir.create(odir)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")

p_load(future)
plan("multiprocess", workers = 40)
options(future.globals.maxSize = 90000 * 1024^2)

data_folder = "figure_5_v1"
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/figure_5_v1/"
setwd(save_path)
col_peri = c("#ECA141", "#D94087","#66528E")

# core --------
setwd("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/")
odir = "figure_5_v1";dir.create(odir)
setwd(paste0("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")

p_load(future)
plan("multiprocess", workers = 100)
options(future.globals.maxSize = Inf)

data_folder = "figure_5_v1"
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_5_v1/"
setwd(save_path)
col_peri = c("#ECA141", "#D94087","#66528E")

# ************************* -----
# 1. MAIN figure -----
# ** peri cell figure -----

# **read data ----
peri = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_pericytes_v2/dim_15k_20_tum.integrated.ns_RNA_v4_intgeratedassay_integrated_snn_res.0.085_vst_method_9.rds")
peri$sheet_patient_id[peri$sheet_ts_sampleid %in% c("hbca25_c")] = "pt37"

# **annotate clusters ----
peri$peri_state = peri$integrated_snn_res.0.085
DefaultAssay(peri) <- 'RNA'
#peri = ScaleData(peri)

peri$final_group = "NULL"
peri$final_group[peri$peri_state %in% c("0")] = "peri_imm"
peri$final_group[peri$peri_state %in% c("1")] = "peri_CREM"
peri$final_group[peri$peri_state %in% c("2")] = "peri_myo"

peri$final_group = factor(peri$final_group, levels = c("peri_imm", "peri_CREM", "peri_myo"))

FeaturePlot(peri, c('rna_COL6A3', "rna_PDPN", "rna_PCOLCE2"))
FeaturePlot(peri, c('rna_COL6A3', "rna_CREM","rna_CXCL3",  "rna_ACTA2"), order = T)
FeaturePlot(peri, c('rna_CMSS1', "rna_APOC1", "rna_ADIPOQ", "rna_FABP4"), order = T)

FeaturePlot(peri, c('rna_VEGFA', "rna_VEGFB", "rna_MYL9", "rna_PROCR"))
d1 = FeaturePlot(peri, c('rna_BMP4', "rna_MCAM", "rna_PECAM1",  "rna_PROCR"), order = T)


# ^^ Final object -------
write_rds(peri, "peri_final_use.rds")

#clus_out = read_rds("peri_final_use.rds")
peri = clus_out

DimPlot(peri, group.by = "final_group", cols = col_peri, label = T)

# markers ----
Idents(peri) <- "final_group"
mks = FindAllMarkers(object = peri, assay = "RNA", logfc.threshold = 0.1, only.pos = T)
write.csv(mks, "peri_mks.csv")
mks = read.csv("peri_mks.csv")
mks$cluster = factor(mks$cluster, levels = c("peri_imm", "peri_CREM", "peri_myo"))

top_markers = mks %>% arrange(cluster) %>% 
  dplyr::group_by(cluster) %>%
  slice_max(order_by = avg_logFC, n = 7) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()

DoHeatmap(peri, top_markers, group.colors = col_peri)

main_markers = c("COL1A1", "PCOLCE2","PDPN")

DefaultAssay(peri) <- "RNA"
Idents(peri) <- "final_group"

tum_sub1 = AverageExpression(peri, features = unique(top_markers), return.seurat = T, assays = "RNA")
Idents(tum_sub1)

p1 = DoHeatmap(tum_sub1, features = top_markers, assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = col_peri, raster = F)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
ggsave(plot = p1, filename = "heatmap_cellstate_peri_allcells_v3.pdf", width = 4, height = 6)

sampled.cells <- WhichCells(peri, downsample = 100)
sampled.cells2 = peri[, sampled.cells]
table(sampled.cells2$peri_state)

peri = ScaleData(peri, assay = "RNA")
peri = ScaleData(peri, assay = "integrated")
p2 = DoHeatmap(peri, cells = sampled.cells, features = top_markers, assay = "RNA", label = FALSE,
               slot = "scale.data", disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = col_peri, raster = F)  +
  scale_fill_gradient2(low = "#29A7B8", mid = "white", high = "#E85041", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
#ggsave(plot = p1, filename = "heatmap_cellstate_myecells_v2.png", device = "png", width = 20, height = 18)

ggsave(plot = p2, filename = "heatmap_cellstate_peri_100cellsRNA_v1.pdf", width = 8, height = 6)


# ** heatmap chemokines ------
vec_chemokines = rownames(peri@assays$RNA@scale.data)[str_detect(rownames(peri@assays$RNA@scale.data), '^CXCL|^CCL')]
vec_chemokines = vec_chemokines[!str_detect(vec_chemokines, 'CCL15-CCL14')]
# vec_krts2 = vec_krts[!vec_krts %in% vec_krts1]
# final_vec_krt = c(vec_krts1,vec_krts2) %>% unique()  

chemo_obj = AverageExpression(peri, features = vec_chemokines, return.seurat = T, assays = "RNA")

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
# test$peri_COL1A1 = as.numeric(test$peri_COL1A1)
# test$peri_PCOLCE2 = as.numeric(test$peri_PCOLCE2)
# test$peri_PDPN = as.numeric(test$peri_PDPN)
# test$V1 = as.character(test$V1)
# 
# test1 = test0 %>% dplyr::arrange(-peri_COL1A1) %>% dplyr::filter(peri_COL1A1>0, peri_PCOLCE2<0, peri_PDPN<0) %>% pull(gene)
# test2 = test0 %>% dplyr::filter(!gene %in% test1) %>% arrange(-peri_PCOLCE2) %>% dplyr::filter(peri_COL1A1<0, peri_PCOLCE2>0, peri_PDPN<0) %>% pull(gene)
# test3 = test0 %>% dplyr::filter(!gene %in% c(test1, test2)) %>% arrange(-peri_PDPN) %>% dplyr::filter(peri_COL1A1<0, peri_PCOLCE2<0, peri_PDPN>0) %>% pull(gene)
# seq_genes = c(test1, test2, test3)
# seq_genes2 = vec_chemokines[!vec_chemokines %in% seq_genes]
# final_genes = c(seq_genes, seq_genes2)

p1 = DoHeatmap(peri, #cells = sampled.cells, 
               features = rownames(chemo_obj)[out$order], 
               assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = col_peri, raster = T)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
p1
ggsave(plot = p1, filename = "heatmap_chemokines_celltype_v2.pdf", width = 10, height = 6)


peri = AddModuleScore(peri,features = list(vec_chemokines), assay = "RNA", name = "chemokines")
v1=VlnPlot(peri, features = "chemokines1", pt.size = 0, group.by = "final_group", cols = col_peri)
ggsave(plot = v1, filename = "vlnplot_chemokines_v1.pdf", width = 4, height = 3)


# ** heatmap collagen ------
peri = ScaleData(peri, assay = "RNA")
Idents(peri) <- "peri_state"
vec_chemokines = rownames(peri@assays$RNA@scale.data)[str_detect(rownames(peri@assays$RNA@scale.data), '^COL')]
#vec_chemokines = vec_chemokines[!str_detect(vec_chemokines, "MMP25-AS1")]
#vec_krts1 = c(rownames(test1)[1:21], rownames(test2))
# vec_krts2 = vec_krts[!vec_krts %in% vec_krts1]
# final_vec_krt = c(vec_krts1,vec_krts2) %>% unique()  

chemo_obj = AverageExpression(peri, features = vec_chemokines, return.seurat = T, assays = "RNA")

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
# test$peri_COL1A1 = as.numeric(test$peri_COL1A1)
# test$peri_PCOLCE2 = as.numeric(test$peri_PCOLCE2)
# test$peri_PDPN = as.numeric(test$peri_PDPN)
# test$V1 = as.character(test$V1)
# 
# test1 = test0 %>% dplyr::arrange(-peri_COL1A1) %>% dplyr::filter(peri_COL1A1>0, peri_PCOLCE2<0, peri_PDPN<0) %>% pull(gene)
# test2 = test0 %>% dplyr::filter(!gene %in% test1) %>% arrange(-peri_PCOLCE2) %>% dplyr::filter(peri_COL1A1<0, peri_PCOLCE2>0, peri_PDPN<0) %>% pull(gene)
# test3 = test0 %>% dplyr::filter(!gene %in% c(test1, test2)) %>% arrange(-peri_PDPN) %>% dplyr::filter(peri_COL1A1<0, peri_PCOLCE2<0, peri_PDPN>0) %>% pull(gene)
# seq_genes = c(test1, test2, test3)
# seq_genes2 = vec_chemokines[!vec_chemokines %in% seq_genes]
# final_genes = c(seq_genes, seq_genes2)

p1 = DoHeatmap(chemo_obj, #cells = sampled.cells, 
               features = rownames(chemo_obj)[out$order], 
               assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = col_peri, raster = T)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
p1
ggsave(plot = p1, filename = "heatmap_collagens_celltype_v3.pdf", width = 4, height = 6)

sampled.cells <- sample(x = Cells(peri), size = 100, replace = F)

p1 = DoHeatmap(peri, cells = sampled.cells, 
               features = rownames(chemo_obj)[out$order], 
               assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = col_peri, raster = T)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
p1
ggsave(plot = p1, filename = "heatmap_collagens_100cells_celltype_v3.pdf", width = 6, height = 4)

peri = AddModuleScore(peri,features = list(vec_chemokines), assay = "RNA", name = "collagens")
v1=VlnPlot(peri, features = "collagens1", pt.size = 0, group.by = "final_group", cols = col_peri)
ggsave(plot = v1, filename = "vlnplot_collagens_v1.pdf", width = 4, height = 3)

FeaturePlot(peri, features = c("collagens1",  "chemokines1"))
DimPlot(peri, label = T)
# **Umap ------

umap_tx_peri = peri@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(tum_meta = peri@meta.data)

umap_tx_peri  =  umap_tx_peri %>% mutate(major_group = factor(tum_meta.final_group, levels = c("peri_imm", "peri_CREM", "peri_myo")))

fill_col = col_peri
library(ggpubr);library(grid)
# grob <- grobTree(textGrob("27720 cells", x=0.7,  y=0.95, hjust=0,
#                           gp=gpar(col="black", fontsize=10, fontface="italic")))
#umap_tx_peri1 <- umap_tx_basal[sample(x = 1:nrow(x = umap_tx_basal)), ]
p1 = ggplot(umap_tx_peri, aes(x=UMAP_1, y=UMAP_2)) + #annotation_custom(grob) + 
  geom_jitter(shape = 21, size = 1, aes(fill = major_group),  color = "black", alpha = 0.9, stroke = 0.04) + 
  scale_fill_manual(values = col_peri) +
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

ggsave(plot = p_fixed, filename = "peri_major_umap.png", device = "png", path = save_path, width = 3, height = 3, dpi = 300, units = "in")
ggsave(plot = p_fixed, filename = "peri_major_umap.pdf", device = "pdf", path = save_path, width = 3, height = 3, dpi = 300, units = "in")

p03 = DimPlot(object = peri, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = col_peri, group.by = "final_group", shuffle = T)+  NoAxes()
pdf(glue("{save_path}/peri_major_uamp_1.pdf"), width=6, height=5)
print(plot_grid(p03))
dev.off()

write_rds(peri, "peri_final.rds")
#peri = read_rds("peri_final.rds")

peri = clus_out

# ** Violi n plots -----------
v1 = VlnPlot(peri, features = "LIF", assay = "RNA", group.by = "final_group", pt.size = 0, cols = col_peri)
ggsave(plot = v1, filename = "periviolin_LIF.pdf", device = "pdf", path = save_path, width = 5, height = 3, dpi = 300, units = "in")

v1 = FeaturePlot(peri, features = "LIF", order = T)
ggsave(plot = v1, filename = "perifeature_LIF.pdf", device = "pdf", path = save_path, width = 4, height = 4, dpi = 300, units = "in")

# **Barplots ---------------------------------------------------------------
fill_col = col_peri
tum = peri
df_meta = data.frame(tum@meta.data, 
                     celltype = tum$final_group, 
                     stringsAsFactors = F) %>% 
  as_tibble() %>% clean_names()

df_meta1 = df_meta
# create a summary table
func_bar(var_df = df_meta1, var_x = "sheet_patient_id", var_fill = "celltype", var_pos = "fill", var_cols = fill_col, celltype = "peri_cells")
func_bar(var_df = df_meta1, var_x = "sheet_figure1_grp_updated", var_fill = "celltype", var_pos = "fill", var_cols = fill_col, celltype = "peri_cells")

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

func_bar(var_df = df_meta3, var_x = "key", var_fill = "celltype", var_pos = "fill", var_cols = fill_col, celltype = "peri_cells_leftright")
func_bar(var_df = df_meta3, var_x = "sheet_tissue_location", var_fill = "celltype", var_pos = "fill", var_cols = fill_col, celltype = "peri_cells_leftright")

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
write_rds(df_meta2_sel, "peri_df_cell_meta2_sel.rds")
df_meta2_sel = read_rds("peri_df_cell_meta2_sel.rds")
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

df_mrg_nuc2$celltype = factor(df_mrg_nuc2$celltype, levels = c("peri_imm", "peri_CREM", "peri_myo"))

write.csv(df_mrg_nuc1, "df_mrg_nuc1_samp_peri.csv")
write.csv(df_mrg_nuc2, "df_mrg_nuc2_pt_peri.csv")

save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_5_v1/peri_all/"
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/figure_5_v1/peri_all/"
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

save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_5_v1/peri_short/"
do_clinical_comp(df_mrg_nuc2 = df_mrg_nuc2_short)
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_5_v1/peri_over/"
do_clinical_comp(df_mrg_nuc2 = df_mrg_nuc2_over)

# pathway ------
Idents(clus_out) <- "final_group"
all_peri = FindAllMarkers(object = clus_out, assay = "RNA", logfc.threshold = 0, only.pos = FALSE)
write.csv(all_peri, paste0(save_path, "peri_markers.csv"))

genes = read.csv("peri_markers.csv")
genes$cluster = as.character(genes$cluster)

gsea_resultsc2hall = fun_gsea(genes = genes, set = "H", sub_cat = NULL, num_avg_logFC = 0.25)
gsea_resultsc2C5 = fun_gsea(genes = genes, set = "C5", sub_cat = "BP", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "lumsec")
gsea_resultsc2C3 = fun_gsea(genes = genes, set = "C3", sub_cat = "TFT:GTRD", num_avg_logFC = 0.25, celltype = "lumsec", doplot = FALSE, ret_collap = FALSE)
gsea_resultsc2C8 = fun_gsea(genes = genes, set = "C8", sub_cat = NULL, num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "lumsec")
gsea_resultsc2C2_reac = fun_gsea(genes = genes, set = "C2", sub_cat = "CP:REACTOME", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "lumsec")
gsea_resultsc2C2_kegg = fun_gsea(genes = genes, set = "C2", sub_cat = "CP:KEGG", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "peri")

clus = c("peri_COL1A1", "peri_PCOLCE2", "peri_PDPN")

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
colnames(p1.1) = c("pathway", "peri_COL1A1", "peri_PCOLCE2","peri_PDPN" )

# rearrange the datafram
test = p1

# for basal
test1 = test %>% dplyr::arrange(-`peri_COL1A1`) %>% dplyr::filter(peri_COL1A1>0, peri_PCOLCE2<0, peri_PDPN<0) %>% pull(pathway)
test2 = test %>% dplyr::filter(!pathway %in% test1) %>% arrange(-peri_PCOLCE2) %>% dplyr::filter(peri_COL1A1<0, peri_PCOLCE2>0, peri_PDPN<0) %>% dplyr::top_n(25, peri_PCOLCE2) %>% pull(pathway)
test3 = test %>% dplyr::filter(!pathway %in% c(test1, test2)) %>% arrange(-peri_PDPN) %>% dplyr::filter(peri_COL1A1<0, peri_PCOLCE2<0, peri_PDPN>0) %>% dplyr::top_n(25, peri_PDPN) %>% pull(pathway)
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


png("peri_htmaps_pathways_top20padj0.05_v1.png", res = 300, units = "in", width = 10, height = 14)
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

genes = mks
top100 <- genes %>% group_by(cluster) %>% dplyr::filter(avg_logFC>0.25) %>% top_n(n = 200, wt = avg_logFC)
top100pval <- subset(top100, rowSums(top100[6] < 0.05) > 0)

df <- top100pval[,7:8]
dfsample <- split(df$gene,df$cluster)
length(dfsample)

#The output of length(dfsample) returns how many clusters you have
#I'm sure there's a better way but you have to make a line like below for each cluster

dfsample$peri_CREM = bitr(dfsample$peri_CREM, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$peri_imm = bitr(dfsample$peri_imm, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dfsample$peri_myo = bitr(dfsample$peri_myo, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")


#do the same here, a line like below for each cluster
genelist <- list("peri_imm" = dfsample$peri_imm$ENTREZID, 
                 "peri_CREM" = dfsample$peri_CREM$ENTREZID,
                 "peri_myo" = dfsample$peri_myo$ENTREZID)
GOclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichGO", OrgDb = "org.Hs.eg.db")
pdf("peri_GOclusterplot_v1.pdf", width = 10, height = 5)
dotplot(GOclusterplot, showCategory = 10)
dev.off()

KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichKEGG")
pdf("peri_keggclusterplot_v1.pdf", width = 10, height = 5)
dotplot(KEGGclusterplot , showCategory = 10)
dev.off()

KEGGclusterplot <- compareCluster(geneCluster = genelist, fun = "enrichDO")
dotplot(KEGGclusterplot , showCategory = 10)


# ************************* -----
# 1. SUPP figure -----
# ** Full heatmap -----------
# supp mean heatmap for fig 1
DefaultAssay(clus_out) <- "RNA"
clus_out = ScaleData(clus_out)
all_markers = read.csv("./peri_mks.csv")
all_markers = all_markers[,-1]
all_markers = all_markers %>% dplyr::filter(!gene == "PLCG2")
all_markers$cluster = factor(all_markers$cluster, levels = lvl)

all_markers_sel2 = all_markers
lvl = levels(clus_out$final_group)
all_markers_sel2$cluster = factor(all_markers_sel2$cluster, levels = lvl)
all_markers_sel2 = all_markers_sel2 %>% dplyr::filter(!gene == "PLCG2")

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

main_markers = c("SERPINE1", "CREM", "ACTA2")

all_markers_sel_fig = data.frame(gene = c(top_basal_markers), cluster1 = c(rep("peri_imm", 7),rep("peri_CREM", 7),rep("peri_myo", 7)))
DefaultAssay(clus_out) <- "integrated"
#clus_out = ScaleData(clus_out)
Idents(clus_out) <- "final_group"
all_markers_sel_fig_genes = unique(c(main_markers, top_basal_markers))
col_grp = col_peri
names(col_grp) = lvl
df_fibro1 = run_mean_htmaps_fig1_supp(df = all_markers_sel_fig, seu = clus_out, cell = "peri_cells", genes_df = all_markers_sel_fig_genes, #sel_genes
                                      genes_df_cann = main_markers, genes_sel_df = genes_sel_df, 
                                      clus_levels = lvl)


# v1 = VlnPlot(clus_out, features = c("FAP"), pt.size = 0, assay = "RNA", group.by = "final_group", cols = fill_col)
# ggsave(v1, file =  "vlnplot_progenitor_genes.pdf", width = 8, height = 7)
# v1 = FeaturePlot(clus_out, features = c("rna_FAP"))
# ggsave(v1, file =  "vlnplot_progenitor_genes.pdf", width = 8, height = 7)

# ** > lumhr ----
to_plots <- c("RGS5", "MCAM")
to_plots <- c("SERPINE1", "CREM", "EGR1")
DefaultAssay(clus_out) <- "RNA"
plots <- lapply(to_plots, function(to_plot) {               
  plt <- FeaturePlot(clus_out, to_plot, order = T) + 
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
ggsave("feature_peri_canonical.pdf", plot = plt, width = 2.5, height = 6, limitsize = FALSE)  
ggsave("feature_peri_states.pdf", plot = plt, width = 2.5, height = 8, limitsize = FALSE)  

