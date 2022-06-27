# geo ----------
setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/")
odir = "figure_6_v1";dir.create(odir)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")

p_load(future)
plan("multiprocess", workers = 40)
options(future.globals.maxSize = 90000 * 1024^2)

data_folder = "03_subset_myeloid_v2"
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/figure_6_v1/"
setwd(save_path)
col_m = c("#f9f871", "#00c9ae", "#ff8863", "#6e3452", "#00aafd", "#b178de", "#ff86da", "#7d9300","red2","#C70E7B", "#FC6882", "#A6E000", "#1BB6AF", "#6C6C9D", "#172869", "#D64358", "#EAFB88", "#3C8C4D", "#DFCEE0","orange", "red", "darkorchid1","darkorchid4", colors_dark)


# core --------
setwd("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/")
odir = "figure_6_v1";dir.create(odir)
setwd(paste0("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")

p_load(future)
plan("multiprocess", workers = 40)
options(future.globals.maxSize = Inf)

data_folder = "03_subset_myeloid_v2"
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_6_v1/"
setwd(save_path)
col_m = c("#f9f871", "#00c9ae", "#ff8863", "#6e3452", "#00aafd", "#b178de", "#ff86da", "#7d9300","red2","#C70E7B", "#FC6882", "#A6E000", "#1BB6AF", "#6C6C9D", "#172869", "#D64358", "#EAFB88", "#3C8C4D", "#DFCEE0","orange", "red", "darkorchid1","darkorchid4", colors_dark)

# ************************* -----
# ** Myeloid cell figure -----

# **read data ----
myeloid = readRDS(file = paste0(wd, data_folder,"/tum_sub2.rds"))
myeloid$sheet_patient_id[myeloid$sheet_ts_sampleid %in% c("hbca25_c")] = "pt37"

DimPlot(myeloid, group.by = "cellstate", cols = col_m)

# **Rename clusters -----
# tum.integrated.ns = myeloid
# tum.integrated.ns$final_group = "NULL"
# tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("cDC1_clecl9a")] = "cDC1"
# tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("cDC2")] = "cDC2"
# tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("mDC")] = "mDC" 
# tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("pDC")] = "pDC"
# tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("mast")] = "mast"
# 
# tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("mono_cd52_cd16")] = "mono_cd16"
# tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("mono_cxcl5_cd14")] = "mono_cd14"
# tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("macro_apoc1_foam")] = "macro_apoc1"
# tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("macro_c3_apoe_m1")] = "macro_c3_apoe_m1"
# tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("macro_c3_plcg2")] = "macro_c3_plcg2_m1"
# tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("macro_ccl_m2")] = "macro_ccl_m2"
# tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("macro_c3_plcg2")] = "macro_c3_plcg2_m1"
# tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("macro_macro_m2")] = "macro_MARCO_m2"
# tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("macro_ifn")] = "macro_ifn"
# tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("prol")] = "mye_prol"
# tum.integrated.ns$final_group1 = factor(tum.integrated.ns$final_group, levels = c("cDC1", "cDC2", "mDC",
#                                                                                   "pDC", "mast", "mono_cd16",
#                                                                                   "mono_cd14", "macro_apoc1",
#                                                                                   "macro_c3_apoe_m1",  "macro_c3_plcg2_m1","macro_ccl_m2",
#                                                                                   "macro_MARCO_m2", "macro_ifn",
#                                                                                   "mye_prol"))

tum.integrated.ns = myeloid
tum.integrated.ns$final_group = "NULL"
tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("cDC1_clecl9a")] = "cDC1"
tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("cDC2")] = "cDC2"
tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("mDC")] = "mDC" 
tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("pDC")] = "pDC"
tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("mast")] = "mast"

tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("mono_cd52_cd16")] = "mono_naive"
tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("mono_cxcl5_cd14")] = "mono_active"
tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("macro_apoc1_foam")] = "macro_lipo"
tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("macro_c3_apoe_m1")] = "macro_m1"
tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("macro_c3_plcg2")] = "macro_m1"
tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("macro_ccl_m2")] = "macro_m2"
tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("macro_c3_plcg2")] = "macro_m1"
tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("macro_macro_m2")] = "macro_m2"
tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("macro_ifn")] = "macro_ifn_m2"
tum.integrated.ns$final_group[tum.integrated.ns$cellstate %in% c("prol")] = "mye_prol"
tum.integrated.ns$final_group1 = factor(tum.integrated.ns$final_group, levels = c("cDC1", "cDC2", "mDC",
                                                                                  "pDC", "mast", "mono_naive",
                                                                                  "mono_active", "macro_lipo",
                                                                                  "macro_m1",  "macro_m2", "macro_ifn_m2",
                                                                                  "mye_prol"))

tum.integrated.ns$cell_group = "NULL"
tum.integrated.ns$cell_group[tum.integrated.ns$cellstate %in% c("cDC1_clecl9a")] = "cDC"
tum.integrated.ns$cell_group[tum.integrated.ns$cellstate %in% c("cDC2")] = "cDC"
tum.integrated.ns$cell_group[tum.integrated.ns$cellstate %in% c("mDC")] = "mDC" 
tum.integrated.ns$cell_group[tum.integrated.ns$cellstate %in% c("pDC")] = "pDC"
tum.integrated.ns$cell_group[tum.integrated.ns$cellstate %in% c("mast")] = "mast"

tum.integrated.ns$cell_group[tum.integrated.ns$cellstate %in% c("mono_cd52_cd16")] = "monocytes"
tum.integrated.ns$cell_group[tum.integrated.ns$cellstate %in% c("mono_cxcl5_cd14")] = "monocytes"
tum.integrated.ns$cell_group[tum.integrated.ns$cellstate %in% c("macro_apoc1_foam")] = "macrophages"
tum.integrated.ns$cell_group[tum.integrated.ns$cellstate %in% c("macro_c3_apoe_m1")] = "macrophages"
tum.integrated.ns$cell_group[tum.integrated.ns$cellstate %in% c("macro_c3_plcg2")] = "macrophages"
tum.integrated.ns$cell_group[tum.integrated.ns$cellstate %in% c("macro_ccl_m2")] = "macrophages"
tum.integrated.ns$cell_group[tum.integrated.ns$cellstate %in% c("macro_c3_plcg2")] = "macrophages"
tum.integrated.ns$cell_group[tum.integrated.ns$cellstate %in% c("macro_macro_m2")] = "macrophages"
tum.integrated.ns$cell_group[tum.integrated.ns$cellstate %in% c("macro_ifn")] = "macrophages"
tum.integrated.ns$cell_group[tum.integrated.ns$cellstate %in% c("prol")] = "mye_prol"
tum.integrated.ns$cell_group1 = factor(tum.integrated.ns$cell_group, levels = c("cDC", "mDC",
                                                                                  "pDC", "mast", "monocytes",
                                                                                  "macrophages", "mye_prol"))

write_rds(tum.integrated.ns, "mye_final_use.rds")
tum = read_rds("mye_final_use.rds")
tum = tum.integrated.ns

tum$sheet_patient_id[tum$sheet_ts_sampleid %in% c("hbca25_c")] = "pt37"

# **Umap ------
tum1 = RunUMAP(tum, dims = 1:12, cols = col_m)
tum2 = RunUMAP(tum, dims = 1:13) #18
DimPlot(tum, cols = col_m, group.by = "final_group1")

umap_tx = tum@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(tum_meta = tum@meta.data)

umap_tx  =  umap_tx %>% mutate(major_group = factor(tum_meta.final_group1, levels = c("cDC1", "cDC2", "mDC",
                                                                                      "pDC", "mast", "mono_naive",
                                                                                      "mono_active", "macro_lipo",
                                                                                      "macro_m1",  "macro_m2", "macro_ifn_m2",
                                                                                      "mye_prol")))

# plot ggplot
# figure 1b umap 
ggthemes::theme_excel()
# fill_col = col_m
# fill_col = c("#FFADD6", "#EA3FE5",  "#783C8C", "#19B6AF", "#0176BC", "#0044DC", "#9DF009", "#009E3F", "#FFDC15", "#F4B084", "#FF6300", "#DD1F00", "#641412", "black")
fill_col = c("blue", "#1f83b4" ,"cyan",  "#2ca030" , "green" , "yellow",  "#ffaa0e" , 
             "tomato1", "#d63a3a" , "#ba43b4" , "pink", "black") ##bcbd22
library(ggpubr);library(grid)
grob <- grobTree(textGrob("21041 cells", x=0.7,  y=0.95, hjust=0,
                          gp=gpar(col="black", fontsize=10, fontface="italic")))
umap_tx1 <- umap_tx[sample(x = 1:nrow(x = umap_tx)), ]
p1 = ggplot(umap_tx1, aes(x=UMAP_1, y=UMAP_2)) + #annotation_custom(grob) + 
  #geom_point(shape = 16, size = 0.005, aes(color = major_group)) + 
  #geom_point(shape = 21, size = 1, aes(fill = major_group), alpha = 0.5, stroke = 0.05, position = position_jitter(w = 0.1, h = 0.1)) + #, position = "dodge"
  geom_jitter(shape = 21, size = 0.5, aes(fill = major_group),  color = "black", alpha = 0.9, stroke = 0.035) + 
  #geom_point(shape = 19, size = 0.5, aes(color = major_group)) + 
  scale_fill_manual(values = fill_col) +
  #ggthemes::calc_pal()(12) +
  #rcartocolor::scale_fill_carto_d(palette = 2) +
  #scale_color_manual(values = fig1_colors) +
  #guides(fill = guide_legend(override.aes = list(color = fill_col)),
  #color = guide_legend(override.aes = list(shape = 21))) +
  xlab("UMAP-1") + ylab("UMAP-2") +
  theme(
          panel.grid = element_blank(), 
  #       panel.background = element_rect(fill = NULL, colour = NULL),
         panel.background = element_blank(),
          legend.position = "none",
        #axis.title = element_blank(),
        #axis.line = element_line(arrow=arrow)
        axis.line.x = element_line(colour = "black"), #,arrow=arrow
        axis.line.y = element_line(colour = "black"),#,arrow=arrow
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank());p1
p_fixed <- set_panel_size(p1,
                          width  = unit(1.5, "in"),
                          height = unit(1.5, "in"))

ggsave(plot = p_fixed, filename = "myeloidl_umap.png", device = "png", path = save_path, width = 3, height = 3, dpi = 300, units = "in")

p03 = DimPlot(object = tum, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = fill_col, group.by = "final_group1", shuffle = T)+  NoAxes()
png(glue("{save_path}/mye_uamp_1.png"), width=6, height=5, units = "in", res = 300)
print(plot_grid(p03))
dev.off()

# **markers ----
# FINAL OBJECT 
# tum = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/figure_6_v1/mye_final_use.rds")
# tum = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_6_v1/mye_final_use.rds")
Idents(tum) = "final_group1" #major_group

all_markers = FindAllMarkers(object = tum, assay = "RNA", logfc.threshold = 0.25, only.pos = TRUE)
write.csv(all_markers, paste0(save_path, "myeloid_fig_markers_final_group1.csv"))
# all_markers = read.csv("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/figure_6_v1/myeloid_fig_markers_final_group1.csv")
# all_markers = all_markers[,-1]

allcell_markers = FindAllMarkers(object = tum, assay = "RNA", logfc.threshold = 0, only.pos = FALSE)
write.csv(allcell_markers, paste0(save_path, "allcell_markers.csv"))
allcell_markers = read.csv("allcell_markers.csv")
allcell_markers = allcell_markers[,-1]


# **Heatmap ------
all_markers_sel2 = all_markers
lvl = levels(tum$final_group1)
all_markers_sel2$cluster = factor(all_markers_sel2$cluster, levels = lvl)


# ## run function to select markers 
# allcell_markers %>%
#   mutate(cluster1 = forcats::fct_relevel(cluster, levels = c("cDC1", "cDC2", "mDC",
#                                                              "pDC", "mast", "mono_naive",
#                                                              "mono_active", "macro_lipo",
#                                                              "macro_m1",  "macro_m2", "macro_ifn_m2",
#                                                              "mye_prol"))) %>% 
#   arrange(cluster1) %>% 
#   filter(!str_detect(gene, "^RP")) %>% 
#   filter(!str_detect(gene, "^MT-")) %>% 
#   filter(!str_detect(gene, "^HBB")) -> all_markers_sel
# 
df_thresh = data.frame(
  cluster = levels(umap_tx$major_group),
  thresh_pct1 = c(0, 0, 0,
                  0, 0, 0,
                  0, 0, 0,0, 0, 0, 0, 0), stringsAsFactors = F
) %>% data.frame(row.names = .$cluster)
# 
# p_load(magrittr)
# all_markers_sel %<>% as_tibble()
# 
# all_markers_sel2 = func_select_markers_thresh(all_markers_sel = all_markers_sel, pct_var = 0.25)
# 


all_markers_sel3 = func_select_markers_thresh(all_markers_sel = all_markers_sel2, pct_var = 0.25)

top_markers = allcell_markers %>% arrange(cluster) %>% 
  dplyr::group_by(cluster) %>%
  top_n(7, avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()

genes_sel_df = allcell_markers %>% arrange(cluster) %>% 
  group_by(cluster) %>% 
  top_n(7, avg_logFC) %>% 
  dplyr::filter(avg_logFC > 0) %>%
  dplyr::select(gene, cluster) %>% mutate(id = paste0(gene, "-", cluster))

DefaultAssay(tum) <- "integrated"
Idents(tum) <- "final_group1"
tum_sub1 = AverageExpression(tum, features = unique(top_markers), return.seurat = T)
levels(tum_sub1$orig.ident) <- levels(Idents(tum))

p1 = DoHeatmap(tum_sub1, features = top_markers, assay = "integrated", label = TRUE,
               slot = "scale.data", disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = fill_col, raster = T)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  rotate_x_text(angle = 90)  + 
  scale_y_discrete(limits=rev) +   
  scale_x_discrete(limits=rev) +
  coord_flip() +
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5));p1
ggsave(plot = p1, filename = "heatmap_cellstate_myeloidcells_v2_paper.pdf", width = 8, height = 2.5)


main_markers = c("PTPRC")
Groups_set <- c("cDC1", "cDC2", "mDC",
                "pDC", "mast", "mono_naive",
                "mono_active", "macro_lipo",
                "macro_m1",  "macro_m2", "macro_ifn_m2",
                "mye_prol")
col_grp <- fill_col[1:12]
names(col_grp) <- Groups_set

df_fibro1 = run_mean_htmaps_fig6_horizontal(df = allcell_markers, seu = tum, cell = "all_cells", genes_df = top_markers, #sel_genes #all_markers_sel2
                                            genes_df_cann = main_markers, genes_sel_df = genes_sel_df, 
                                            clus_levels = c("cDC1", "cDC2", "mDC",
                                                            "pDC", "mast", "mono_naive",
                                                            "mono_active", "macro_lipo",
                                                            "macro_m1",  "macro_m2", "macro_ifn_m2",
                                                            "mye_prol"))

# **Barplots ---------------------------------------------------------------
df_meta = data.frame(tum@meta.data, 
                     celltype = tum$final_group1, 
                     stringsAsFactors = F) %>% 
  as_tibble() %>% clean_names()

#df_meta = df_meta %>% mutate(orig_pt_ident = factor(orig_pt_ident, levels = c("HR+ Luminal Cells", "Secretory Luminal Cells", "Endothelial",  "Fibroblasts",  "Pericytes", "T cells", "Macrophages",
#                                                              "TumA", "TumB", "TumC", "TumD")))
#df_meta1 = df_meta %>% mutate(patient_id =  factor(patient_id, levels = unique(patient_id[order(source, procedure)]), ordered = TRUE))
# fn = factor(f, levels=unique(f[order(a,b,f)]), ordered=TRUE)
df_meta1 = df_meta
# create a summary table
func_bar(var_df = df_meta1, var_x = "sheet_patient_id", var_fill = "celltype", var_pos = "fill", var_cols = col_grp, celltype = "cells")
func_bar(var_df = df_meta1, var_x = "sheet_figure1_grp_updated", var_fill = "celltype", var_pos = "fill", var_cols = col_grp, celltype = "cells")

df_meta2 = df_meta1 %>% dplyr::filter(sheet_tissue_location %in% c("left", "right"), !sheet_patient_id == "pt45")
df_meta3 = df_meta2 %>% mutate(sheet_patient_id =  factor(sheet_patient_id, levels = unique(sheet_patient_id[order(sheet_patient_id, sheet_tissue_location)]), ordered = TRUE))
df_meta3 = df_meta3 %>% mutate(key = paste0(sheet_patient_id, "_", sheet_tissue_location))

p1 = ggplot(df_meta3, aes(sheet_tissue_location, fill = celltype)) +
  geom_bar(position = "fill", width = 1) + 
  xlab("") + ylab("") + scale_fill_manual(values = col_grp) +
  theme(axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
        panel.background = element_blank(),axis.title.x = element_blank(), 
        axis.text.x = element_blank(), axis.ticks = element_blank(), 
        legend.title = element_blank()) + facet_wrap(~sheet_patient_id, ncol = 22, strip.position = "bottom") + 
  theme(strip.background = element_blank(), panel.spacing.x=unit(0.2, "lines"), strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(2,0,0,0, "pt")))
ggsave(p1, file = paste0(save_path, "/figure6_myeloid_tissue_location_patient.pdf"), width = 7, height = 2.5)

func_bar(var_df = df_meta3, var_x = "key", var_fill = "celltype", var_pos = "fill", var_cols = col_grp, celltype = "cells_leftright")
func_bar(var_df = df_meta3, var_x = "sheet_tissue_location", var_fill = "celltype", var_pos = "fill", var_cols = col_grp, celltype = "cells_leftright")

df_meta2 = df_meta1 %>% dplyr::filter(sheet_bmi_final %in% c("normal","Normal","obese", "overweight"))
df_meta2 = df_meta2 %>% dplyr::mutate(bmi = case_when(sheet_bmi_final %in% c("normal", "Normal") ~ "normal",
                                                      sheet_bmi_final %in% c("obese") ~ "obese",
                                                      sheet_bmi_final %in% c("overweight") ~ "overweight",
                                                      sheet_bmi_final %in% c("unknown") ~ "unknown",
                                                      is.na(sheet_bmi_final)  ~ "unknown"))
df_meta2$bmi = factor(df_meta2$bmi, levels = c("normal", "overweight", "obese"), ordered = T) 
df_meta3 = df_meta2 %>% mutate(sheet_patient_id =  factor(sheet_patient_id, levels = unique(sheet_patient_id[order(sheet_patient_id, bmi)]), ordered = TRUE))
df_meta3 = df_meta3 %>% mutate(key = paste0(sheet_patient_id, "_", bmi))

func_bar(var_df = df_meta3, var_x = "bmi", var_fill = "celltype", var_pos = "fill", var_cols = col_grp, celltype = "cells_bmi_average")

df_meta2 = df_meta2 %>% dplyr::filter(!sheet_figure1_grp_updated %in% c("unknown"), !is.na(sheet_figure1_grp_updated))
df_meta2$sheet_figure1_grp_updated = factor(df_meta2$sheet_figure1_grp_updated, levels = c("Reduction Mammoplasty", "Prophylatic Mastectomy", "Cancer Mastectomy"), ordered = T)
df_meta3 = df_meta2 %>% mutate(sheet_patient_id =  factor(sheet_patient_id, levels = unique(sheet_patient_id[order(sheet_patient_id, sheet_figure1_grp_updated)]), ordered = TRUE))
df_meta3 = df_meta3 %>% mutate(key = paste0(sheet_patient_id, "_", sheet_figure1_grp_updated))

func_bar(var_df = df_meta3, var_x = "sheet_figure1_grp_updated", var_fill = "celltype", var_pos = "fill", var_cols = col_grp, celltype = "cells_sheet_figure1_grp_updated_average")


# **Clinical ---------
df_meta2_sel = df_meta1 %>% dplyr::select(cell_group1, final_group1, final_sample_id, ethnicity, age, sheet_patient_id, sheet_source, sheet_exp_proc, sheet_age_updated, sheet_brca_updated, sheet_parity, sheet_gravida, sheet_menopause, sheet_density, sheet_bmi_final, sheet_cancer_risk,sheet_figure1_grp_updated)
colnames(df_meta2_sel) = c("cellgroup", "celltype", "sample_id", "ethnicity", "age_orig", "patient_id",  "source", "exp_proc", "age", "brca", "parity", "gravida", "menopause", "density", "bmi", "cancer_risk", "surgery")
# df_mrg = df_meta2_sel %>% mutate(pt_x = paste0(patient_id, "_", type))
write_rds(df_meta2_sel, "myeloid_df_cell_meta2_sel.rds")

df_meta2_sel = read_rds("myeloid_df_cell_meta2_sel.rds")

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


df_mrg_nuc1 = df_meta2_sel %>% group_by(sample_id) %>% mutate(total_cells_per_samp = n()) %>% ungroup() %>% 
  group_by(sample_id, celltype) %>% mutate(cell_per_samp = n(), prop_cell_per_samp = cell_per_samp/total_cells_per_samp) %>% ungroup() %>% 
  group_by(sample_id, cellgroup) %>% mutate(cell_per_grp_samp = n(), prop_cell_per_grp_samp = cell_per_samp/cell_per_grp_samp) %>% ungroup()

df_mrg_nuc2 = df_meta2_sel %>% group_by(patient_id) %>% mutate(total_cells_per_pt = n()) %>% ungroup() %>% 
  group_by(patient_id, celltype) %>% mutate(cell_per_pt = n(), prop_cell_per_pt = cell_per_pt/total_cells_per_pt) %>% ungroup() %>% 
  group_by(patient_id, cellgroup) %>% mutate(cell_per_grp_pt = n(), prop_cell_per_grp_pt = cell_per_pt/cell_per_grp_pt) %>% ungroup()


write.csv(df_mrg_nuc1, "df_mrg_nuc1_samp_myeloid.csv")
write.csv(df_mrg_nuc2, "df_mrg_nuc2_pt_myeloid.csv")

save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_6_v1/myeloid/"
col_pal = col_grp
do_clinical_comp(df_mrg_nuc2 = df_mrg_nuc2)

df_meta2_sel_short =  df_meta2_sel %>% dplyr::filter(exp_proc_updated == "short")
df_meta2_sel_over =  df_meta2_sel %>% dplyr::filter(exp_proc_updated == "overnight")

df_mrg_nuc2_short = df_meta2_sel_short %>% group_by(patient_id) %>% mutate(total_cells_per_pt = n()) %>% ungroup() %>% 
  group_by(patient_id, celltype) %>% mutate(cell_per_pt = n(), prop_cell_per_pt = cell_per_pt/total_cells_per_pt) %>% ungroup() %>% 
  group_by(patient_id, cellgroup) %>% mutate(cell_per_grp_pt = n(), prop_cell_per_grp_pt = cell_per_pt/cell_per_grp_pt) %>% ungroup()
df_mrg_nuc2_over = df_meta2_sel_over %>% group_by(patient_id) %>% mutate(total_cells_per_pt = n()) %>% ungroup() %>% 
  group_by(patient_id, celltype) %>% mutate(cell_per_pt = n(), prop_cell_per_pt = cell_per_pt/total_cells_per_pt) %>% ungroup() %>% 
  group_by(patient_id, cellgroup) %>% mutate(cell_per_grp_pt = n(), prop_cell_per_grp_pt = cell_per_pt/cell_per_grp_pt) %>% ungroup()

do_clinical_comp(df_mrg_nuc2 = df_mrg_nuc2_short)
do_clinical_comp(df_mrg_nuc2 = df_mrg_nuc2_over)


# ************************* -----
# ** B cell figure -----

# read data ----
b = readRDS(file = paste0(wd, "03_subset_B_v2/dim_12k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")) 
b$sheet_patient_id[b$sheet_ts_sampleid %in% c("hbca25_c")] = "pt37"
DimPlot(b, group.by = "cellstate", cols = col_m)


FeaturePlot(b, c("rna_MS4A1", "rna_JCHAIN"))

FeaturePlot(b, c("rna_SELL", "rna_CD27", "rna_IGHM", "rna_IGHD"), order = T)

# plasma B
FeaturePlot(b, c("rna_IGHM", "rna_SDC1", "rna_CD38", "rna_SLAMF7"))
# plasma B cells
FeaturePlot(b, c("rna_IGHA1","rna_IGHA2", "rna_IGHG1", "rna_IGHG2"))

# plasma B cells
FeaturePlot(b, c("rna_KLF2","rna_HLA-DQA2", "rna_TCL1A", "rna_BCL7A"))
FeaturePlot(b, c("rna_MS4A1","rna_CD38", "rna_MME","rna_TCL1A", "rna_BCL7A"))

# activated naive
FeaturePlot(b, c("rna_HLA-DRB1","rna_HLA-DRA", "rna_HLA-DPA1", "rna_HLA-DQA1"))
# activated naive
FeaturePlot(b, c("rna_MYC","rna_EGR3", "rna_KRT17", "rna_BCL2A1"))
FeaturePlot(b, c("rna_CD70","rna_CD69", "rna_NR4A3", "rna_NR4A1"))


b$final_group = "NULL"
b$final_group[b$cellstate %in% c("b_naive_cd27")] = "bmem_switched"
b$final_group[b$cellstate %in% c("b_naive_activated_mhc2")] = "b_naive"
b$final_group[b$cellstate %in% c("b_naive_myc")] = "bmem_unswitched"
b$final_group[b$cellstate %in% c("b_naive_tcl1a")] = "b_naive"

b$final_group[b$cellstate %in% c("b_plasma_igh_igk")] = "plasma_IgG"
b$final_group[b$cellstate %in% c("b_plasma_igl")] = "plasma_IgA"

write_rds(b, "b_final_use.rds")
#b = read_rds("b_final_use.rds")


DimPlot(b, group.by = "final_group", cols = colors_dark, label = T)
Idents(b) <- "final_group"
mks = FindAllMarkers(b, assay = "RNA")


top_markers = mks %>% arrange(cluster) %>% 
  dplyr::group_by(cluster) %>%
  top_n(7, avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()

DoHeatmap(b, top_markers, group.colors = colors_dark)

# rna_CXCL12 in all other cells
FeaturePlot(b, c("rna_CXCL12", "rna_CXCR4", "rna_CCR6","rna_CXCR5"))
FeaturePlot(b, c("rna_CD24", "rna_CD27", "rna_PTPRC"), order = T)
FeaturePlot(b, c("rna_CD40", "rna_CD83", "rna_FCER2", "rna_CXCR4"), order = T)
FeaturePlot(b, c("rna_CD38", "rna_MME", "rna_CD24", "rna_IGHD"), order = T)
FeaturePlot(b, c("rna_CD86", "rna_FAS", "rna_CD27", "rna_IGHD"), order = T)


FeaturePlot(b, c("rna_CD19", "rna_MS4A1", "rna_CD27", "rna_IGHD"), order = T, pt.size = 0.01)
RidgePlot(b, group.by = "final_group", features = c("rna_CD19", "rna_MS4A1", "rna_CD27", "rna_IGHD"), 
          cols = c("#f9f871", "#00c9ae", "#ff8863",  "#00aafd", "#ff86da"), ncol = 4)
FeaturePlot(b, c("rna_CD19", "rna_MS4A1", "rna_CD27", "rna_IGHD"), order = T, pt.size = 0.01)
RidgePlot(b, group.by = "final_group", features = c("rna_IGHM"), 
          cols = c("#f9f871", "#00c9ae", "#ff8863",  "#00aafd", "#ff86da"), ncol = 4)

RidgePlot(b, group.by = "final_group", features = c("rna_MME", "rna_CD24", "rna_CD38", "rna_IGHD"), 
          cols = c("#f9f871", "#00c9ae", "#ff8863",  "#00aafd", "#ff86da"), ncol = 4)
RidgePlot(b, group.by = "final_group", features = c("rna_CD19", "rna_MS4A1", "rna_CD27", "rna_IGHD"), 
          cols = c("#f9f871", "#00c9ae", "#ff8863",  "#00aafd", "#ff86da"), ncol = 4)
RidgePlot(b, group.by = "final_group", features = c("rna_CD19", "rna_MS4A1", "rna_CD27", "rna_IGHD"), 
          cols = c("#f9f871", "#00c9ae", "#ff8863",  "#00aafd", "#ff86da"), ncol = 4)


table(b$final_group, b$sheet_brca_updated)
table(b$final_group, b$sheet_patient_id)

# from umap

# **Umap ------
tum = b
umap_tx = tum@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(tum_meta = tum@meta.data)

umap_tx  =  umap_tx %>% mutate(major_group = factor(tum_meta.final_group, levels = c("b_naive","bmem_switched",
                                                                                   "bmem_unswitched","plasma_IgA","plasma_IgG")))

#fill_col = c("blue", "#1f83b4" ,"cyan",  "#2ca030" , "#bcbd22" ,"green",  "#ffaa0e" , "yellow", "tomato1", "#d63a3a" , "#ba43b4" , "#6f63bb", "pink", "black")
fill_col = c("#f9f871", "#00c9ae", "#ff8863",  "#00aafd", "#ff86da")
library(ggpubr);library(grid)
grob <- grobTree(textGrob("6977 cells", x=0.7,  y=0.95, hjust=0,
                          gp=gpar(col="black", fontsize=10, fontface="italic")))
umap_tx1 <- umap_tx[sample(x = 1:nrow(x = umap_tx)), ]
p1 = ggplot(umap_tx1, aes(x=UMAP_1, y=UMAP_2)) + #annotation_custom(grob) + 
  #geom_point(shape = 16, size = 0.005, aes(color = major_group)) + 
  #geom_point(shape = 21, size = 1, aes(fill = major_group), alpha = 0.5, stroke = 0.05, position = position_jitter(w = 0.1, h = 0.1)) + #, position = "dodge"
  geom_jitter(shape = 21, size = 1, aes(fill = major_group),  color = "black", alpha = 0.9, stroke = 0.04) + 
  #geom_point(shape = 19, size = 0.5, aes(color = major_group)) + 
  scale_fill_manual(values = fill_col) +
  #ggthemes::calc_pal()(12) +
  #rcartocolor::scale_fill_carto_d(palette = 2) +
  #scale_color_manual(values = fig1_colors) +
  #guides(fill = guide_legend(override.aes = list(color = fill_col)),
  #color = guide_legend(override.aes = list(shape = 21))) +
  xlab("UMAP-1") + ylab("UMAP-2") +
  theme(
    panel.grid = element_blank(), 
    #       panel.background = element_rect(fill = NULL, colour = NULL),
    panel.background = element_blank(),
    legend.position = "none",
    #axis.title = element_blank(),
    #axis.line = element_line(arrow=arrow)
    axis.line.x = element_line(colour = "black"), #,arrow=arrow
    axis.line.y = element_line(colour = "black"),#,arrow=arrow
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),axis.text.y = element_blank());p1
p_fixed <- set_panel_size(p1,
                          width  = unit(1.5, "in"),
                          height = unit(1.5, "in"))

ggsave(plot = p_fixed, filename = "b_umap.png", device = "png", path = save_path, width = 3, height = 3, dpi = 300, units = "in")
ggsave(plot = p_fixed, filename = "b_umap.pdf", device = "pdf", path = save_path, width = 3, height = 3, dpi = 300, units = "in")

p03 = DimPlot(object = tum, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = fill_col, group.by = "final_group", shuffle = T)+  NoAxes()
pdf(glue("{save_path}/b_uamp_1.pdf"), width=6, height=5)
print(plot_grid(p03))
dev.off()

# **markers ----
# FINAL OBJECT 
# tum = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/figure_6_v1/mye_final_use.rds")
# tum = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_6_v1/mye_final_use.rds")
tum$final_group = factor(tum$final_group, levels = levels(umap_tx$major_group))
Idents(tum) = "final_group"

all_markers = FindAllMarkers(object = tum, assay = "RNA", logfc.threshold = 0.25, only.pos = TRUE)
write.csv(all_markers, paste0(save_path, "b_fig_markers_final_group.csv"))
# all_markers = read.csv("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/figure_6_v1/myeloid_fig_markers_final_group1.csv")
# all_markers = all_markers[,-1]

allcell_markers = FindAllMarkers(object = tum, assay = "RNA", logfc.threshold = 0, only.pos = FALSE)
write.csv(allcell_markers, paste0(save_path, "b_allcell_markers.csv"))
allcell_markers = read.csv("b_allcell_markers.csv")
allcell_markers = allcell_markers[,-1]

# **WRITE Object ------
write_rds(tum, "b_final.rds")
tum = read_rds("b_final.rds")

# **Heatmap ------
all_markers_sel2 = allcell_markers
lvl = levels(tum$final_group)
all_markers_sel2$cluster = factor(all_markers_sel2$cluster, levels = lvl)


top_markers = all_markers_sel2 %>% arrange(cluster) %>% 
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(n = 7, order_by = avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()

genes_sel_df = all_markers_sel2 %>% arrange(cluster) %>% 
  group_by(cluster) %>% 
  dplyr::slice_max(n = 7, order_by = avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  dplyr::select(gene, cluster) %>% mutate(id = paste0(gene, "-", cluster))

DefaultAssay(tum) <- "integrated"
Idents(tum) <- "final_group"
tum_sub1 = AverageExpression(tum, features = unique(top_markers), return.seurat = T, assays = "integrated")
levels(tum_sub1$orig.ident) <- levels(Idents(tum))

p1 = DoHeatmap(tum_sub1, features = c(top_markers), assay = "integrated", label = TRUE,
               slot = "scale.data", disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = fill_col, raster = T)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  rotate_x_text(angle = 90)  + 
  scale_y_discrete(limits=rev) +   
  scale_x_discrete(limits=rev) +
  coord_flip() +
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5));p1
ggsave(plot = p1, filename = "heatmap_cellstate_bcells_v2_paper.pdf", width = 8, height = 2.5)


main_markers = c("MS4A1", "JCHAIN")
Groups_set <- c("b_naive", "bmem_switched", "bmem_unswitched", "plasma_IgA", "plasma_IgG")
col_grp <- fill_col[1:5]
names(col_grp) <- Groups_set

df_fibro1 = run_mean_htmaps_fig6_horizontal(df = all_markers_sel2, seu = tum, cell = "b_cells", genes_df = top_markers, #sel_genes
                                            genes_df_cann = main_markers, genes_sel_df = genes_sel_df, 
                                            clus_levels = c("b_naive", "bmem_switched", "bmem_unswitched", "plasma_IgA", "plasma_IgG"))

save.image(file = "till_bcells.RData")


# **Barplots ---------------------------------------------------------------
df_meta = data.frame(tum@meta.data, 
                     celltype = tum$final_group, 
                     stringsAsFactors = F) %>% 
  as_tibble() %>% clean_names()

df_meta1 = df_meta
# create a summary table
func_bar(var_df = df_meta1, var_x = "sheet_patient_id", var_fill = "celltype", var_pos = "fill", var_cols = col_grp, celltype = "b_cells")
func_bar(var_df = df_meta1, var_x = "sheet_figure1_grp_updated", var_fill = "celltype", var_pos = "fill", var_cols = col_grp, celltype = "b_cells")

df_meta2 = df_meta1 %>% dplyr::filter(sheet_tissue_location %in% c("left", "right"), !sheet_patient_id == "pt45")
df_meta3 = df_meta2 %>% mutate(sheet_patient_id =  factor(sheet_patient_id, levels = unique(sheet_patient_id[order(sheet_patient_id, sheet_tissue_location)]), ordered = TRUE))
df_meta3 = df_meta3 %>% mutate(key = paste0(sheet_patient_id, "_", sheet_tissue_location))

p1 = ggplot(df_meta3, aes(sheet_tissue_location, fill = celltype)) +
  geom_bar(position = "fill", width = 1) + 
  xlab("") + ylab("") + scale_fill_manual(values = col_grp) +
  theme(axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
        panel.background = element_blank(),axis.title.x = element_blank(), 
        axis.text.x = element_blank(), axis.ticks = element_blank(), 
        legend.title = element_blank()) + facet_wrap(~sheet_patient_id, ncol = 22, strip.position = "bottom") + 
  theme(strip.background = element_blank(), panel.spacing.x=unit(0.2, "lines"), strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(2,0,0,0, "pt")))
ggsave(p1, file = paste0(save_path, "/figure6_b_tissue_location_patient.pdf"), width = 7, height = 2.5)

func_bar(var_df = df_meta3, var_x = "key", var_fill = "celltype", var_pos = "fill", var_cols = col_grp, celltype = "b_cells_leftright")
func_bar(var_df = df_meta3, var_x = "sheet_tissue_location", var_fill = "celltype", var_pos = "fill", var_cols = col_grp, celltype = "b_cells_leftright")

df_meta2 = df_meta1 %>% dplyr::filter(sheet_bmi_final %in% c("normal","Normal","obese", "overweight"))
df_meta2 = df_meta2 %>% dplyr::mutate(bmi = case_when(sheet_bmi_final %in% c("normal", "Normal") ~ "normal",
                                                      sheet_bmi_final %in% c("obese") ~ "obese",
                                                      sheet_bmi_final %in% c("overweight") ~ "overweight",
                                                      sheet_bmi_final %in% c("unknown") ~ "unknown",
                                                      is.na(sheet_bmi_final)  ~ "unknown"))
df_meta2$bmi = factor(df_meta2$bmi, levels = c("normal", "overweight", "obese"), ordered = T) 
df_meta3 = df_meta2 %>% mutate(sheet_patient_id =  factor(sheet_patient_id, levels = unique(sheet_patient_id[order(sheet_patient_id, bmi)]), ordered = TRUE))
df_meta3 = df_meta3 %>% mutate(key = paste0(sheet_patient_id, "_", bmi))

func_bar(var_df = df_meta3, var_x = "bmi", var_fill = "celltype", var_pos = "fill", var_cols = col_grp, celltype = "b_cells_bmi_average")

df_meta2 = df_meta2 %>% dplyr::filter(!sheet_figure1_grp_updated %in% c("unknown"), !is.na(sheet_figure1_grp_updated))
df_meta2$sheet_figure1_grp_updated = factor(df_meta2$sheet_figure1_grp_updated, levels = c("Reduction Mammoplasty", "Prophylatic Mastectomy", "Cancer Mastectomy"), ordered = T)
df_meta3 = df_meta2 %>% mutate(sheet_patient_id =  factor(sheet_patient_id, levels = unique(sheet_patient_id[order(sheet_patient_id, sheet_figure1_grp_updated)]), ordered = TRUE))
df_meta3 = df_meta3 %>% mutate(key = paste0(sheet_patient_id, "_", sheet_figure1_grp_updated))

func_bar(var_df = df_meta3, var_x = "sheet_figure1_grp_updated", var_fill = "celltype", var_pos = "fill", var_cols = col_grp, celltype = "b_cells_sheet_figure1_grp_updated_average")

# **Clinical ---------
df_meta2_sel = df_meta1 %>% dplyr::select(celltype, final_group, final_sample_id, ethnicity, age, sheet_patient_id, sheet_source, sheet_exp_proc, sheet_age_updated, sheet_brca_updated, sheet_parity, sheet_gravida, sheet_menopause, sheet_density, sheet_bmi_final, sheet_cancer_risk,sheet_figure1_grp_updated)
colnames(df_meta2_sel) = c("cellgroup", "celltype", "sample_id", "ethnicity", "age_orig", "patient_id",  "source", "exp_proc", "age", "brca", "parity", "gravida", "menopause", "density", "bmi", "cancer_risk", "surgery")
# df_mrg = df_meta2_sel %>% mutate(pt_x = paste0(patient_id, "_", type))
write_rds(df_meta2_sel, "b_df_cell_meta2_sel.rds")

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


df_mrg_nuc1 = df_meta2_sel %>% group_by(sample_id) %>% mutate(total_cells_per_samp = n()) %>% ungroup() %>% 
  group_by(sample_id, celltype) %>% mutate(cell_per_samp = n(), prop_cell_per_samp = cell_per_samp/total_cells_per_samp) %>% ungroup() %>% 
  group_by(sample_id, cellgroup) %>% mutate(cell_per_grp_samp = n(), prop_cell_per_grp_samp = cell_per_samp/cell_per_grp_samp) %>% ungroup()

df_mrg_nuc2 = df_meta2_sel %>% group_by(patient_id) %>% mutate(total_cells_per_pt = n()) %>% ungroup() %>% 
  group_by(patient_id, celltype) %>% mutate(cell_per_pt = n(), prop_cell_per_pt = cell_per_pt/total_cells_per_pt) %>% ungroup() %>% 
  group_by(patient_id, cellgroup) %>% mutate(cell_per_grp_pt = n(), prop_cell_per_grp_pt = cell_per_pt/cell_per_grp_pt) %>% ungroup()


write.csv(df_mrg_nuc1, "df_mrg_nuc1_samp_B.csv")
write.csv(df_mrg_nuc2, "df_mrg_nuc2_pt_B.csv")

save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_6_v1/B/"
col_pal = col_grp
do_clinical_comp(df_mrg_nuc2 = df_mrg_nuc2)

df_meta2_sel_short =  df_meta2_sel %>% dplyr::filter(exp_proc_updated == "short")
df_meta2_sel_over =  df_meta2_sel %>% dplyr::filter(exp_proc_updated == "overnight")

df_mrg_nuc2_short = df_meta2_sel_short %>% group_by(patient_id) %>% mutate(total_cells_per_pt = n()) %>% ungroup() %>% 
  group_by(patient_id, celltype) %>% mutate(cell_per_pt = n(), prop_cell_per_pt = cell_per_pt/total_cells_per_pt) %>% ungroup() %>% 
  group_by(patient_id, cellgroup) %>% mutate(cell_per_grp_pt = n(), prop_cell_per_grp_pt = cell_per_pt/cell_per_grp_pt) %>% ungroup()
df_mrg_nuc2_over = df_meta2_sel_over %>% group_by(patient_id) %>% mutate(total_cells_per_pt = n()) %>% ungroup() %>% 
  group_by(patient_id, celltype) %>% mutate(cell_per_pt = n(), prop_cell_per_pt = cell_per_pt/total_cells_per_pt) %>% ungroup() %>% 
  group_by(patient_id, cellgroup) %>% mutate(cell_per_grp_pt = n(), prop_cell_per_grp_pt = cell_per_pt/cell_per_grp_pt) %>% ungroup()

do_clinical_comp(df_mrg_nuc2 = df_mrg_nuc2_short)
do_clinical_comp(df_mrg_nuc2 = df_mrg_nuc2_over)


# ************************* -----
# ** T cell figure -----

# read data ----
t = readRDS(file = paste0(wd, "03_subset_T_v2/T_clus_out.rds")) 
DimPlot(t, group.by = "cellstate", cols = col_m)

Idents(t) <- "cellstate"
t_sub = subset(t, idents = c("doublets", "low_quality", "stressed"), invert = TRUE) #endo_vas
write_rds(t_sub, "t_sub.rds")

t_sub = read_rds("t_sub.rds")
t_sub = RunUMAP(t_sub, dims = 1:30)
DimPlot(t_sub, group.by = "cellstate", cols = col_m, label = T)
DimPlot(t_sub, group.by = "CD4CD8", cols = col_m, label = T)
DefaultAssay(t_sub) <- "RNA"
dims = 30; k.param = 20; cols = colors_dark; var_split = "sheet_exp_proc"; var_scale = "nCount_RNA"
fts1 = c("nFeature_RNA","nCount_RNA", "percent.mito","percent.rb", "percent.stress")
fts2 = c("rna_EPCAM","rna_PTPRC","rna_PECAM1","rna_DCN","rna_KRT5",  "rna_CD86", "rna_COL1A1", "rna_PIP")
fts3 = c("rna_CD8A", "rna_CD8B", "rna_CD4", "rna_CCR7", "rna_IL7R", "rna_NKG7", "rna_GZMB", "rna_GZMK")
fts4 = c("rna_MS4A1","rna_CD40LG","rna_FOXP3","rna_IFNG","rna_LAG3","rna_TIGIT",  "rna_NCAM1", "rna_KLRB1", "rna_CD3D")

ref_set = NULL
clus_out = func_cluster_vst(tum_int = t_sub, var_scale = var_scale, var_split = var_split, dims = dims, k.param = k.param, 
                            markers1 = fts1, markers2 = fts2, markers3 = fts3, markers4 = fts4, save_path = save_path,
                            cols = colors_dark, ref_set = ref_set)
clus_out = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_6_v1/dim_35k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")

p1 = DimPlot(clus_out, group.by = "integrated_snn_res.0.6", cols = col_m, label = T)
p2 = DimPlot(clus_out, group.by = "cellstate", cols = col_m, label = T)
p = cowplot::plot_grid(p1, p2)
char_ident = "integrated_snn_res.0.6"
dims = 35; k.param = 20;
char_ident = "integrated_snn_res.0.6"  #integrated_snn_res.0.1
ht_fts1 = c("CD3D", "CD8A", "CD8B")
#clus_out = read_rds(path = glue("{save_path}/dim_{dims}k_{k.param}_tum.integrated.ns_SCT_v4_intgeratedassay_{char_ident}_vst_method9.rds")) 
mark_out = func_markers_vst(clus_out = clus_out, char_ident = char_ident, fts1 = ht_fts1, dims = dims, k.param = k.param, cols = colors_dark, 
                            save_path = save_path)

mks = read.csv("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_6_v1/dims_35k_20_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.6_vst_method_9.csv")
top_markers = mks %>% arrange(cluster) %>% 
  dplyr::group_by(cluster) %>%
  top_n(7, avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()

DoHeatmap(b, top_markers, group.colors = colors_dark)
Idents(clus_out) <- "integrated_snn_res.0.6"
mks2_4 = FindMarkers(clus_out, ident.1 = "2", ident.2 = "4")
mks2_4$gene = rownames(mks2_4)
View(mks2_4)
top_markers_mks2_4 = mks2_4 %>% 
  top_n(20, avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()
bot_markers_mks2_4 = mks2_4 %>% 
  top_n(20, -avg_logFC) %>%
  #dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()

clus_out_2_4 = subset(clus_out, idents = c("2", "4"))
DoHeatmap(clus_out_2_4, c(top_markers_mks2_4,bot_markers_mks2_4),  group.colors = colors_dark)


clus_out_7_16_17 = subset(clus_out, idents = c("7", "16", "17"))
mks7_16_17 = FindAllMarkers(clus_out_7_16_17)
mks7_16_17$gene = rownames(mks7_16_17)
top_markers_mks7_16_17 = mks7_16_17 %>% arrange(cluster) %>% 
  dplyr::group_by(cluster) %>%
  top_n(15, avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()
DoHeatmap(clus_out_7_16_17, top_markers_mks7_16_17,  group.colors = col_m)



clus_out_9_11_13_20 = subset(clus_out, idents = c("9", "11", "13", "20"))
mks9_11_13_20 = FindAllMarkers(clus_out_9_11_13_20)
mks9_11_13_20$gene = rownames(mks9_11_13_20)
top_markers_mks9_11_13_20 = mks9_11_13_20 %>% arrange(cluster) %>% 
  dplyr::group_by(cluster) %>%
  top_n(15, avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()
DoHeatmap(clus_out_9_11_13_20, top_markers_mks9_11_13_20,  group.colors = col_m)


FeaturePlot(clus_out, c("rna_IL7R"))
DimPlot(clus_out, group.by = "integrated_snn_res.0.4", cols = colors_dark)
DimPlot(clus_out, group.by = "CD4CD8")
DimPlot(clus_out, group.by = "sheet_patient_id", cols = colors_dark)
VlnPlot(clus_out, features = "rna_PLCG2", group.by = "sheet_patient_id", cols = colors_dark, pt.size = 0) +
  facet_wrap(~clus_out$CD4CD8, scales = "free")

FeaturePlot(b, c("rna_CD70","rna_CD69", "rna_NR4A3", "rna_NR4A1"))
VlnPlot(clus_out, features = "nCount_RNA", pt.size = 0, cols = col_m)

b$final_group = "NULL"
b$final_group[b$cellstate %in% c("b_naive_cd27")] = "bmem_switched"
b$final_group[b$cellstate %in% c("b_naive_activated_mhc2")] = "b_naive"
b$final_group[b$cellstate %in% c("b_naive_myc")] = "bmem_unswitched"
b$final_group[b$cellstate %in% c("b_naive_tcl1a")] = "b_naive"

b$final_group[b$cellstate %in% c("b_plasma_igh_igk")] = "plasma_IgG"
b$final_group[b$cellstate %in% c("b_plasma_igl")] = "plasma_IgA"

write_rds(b, "b_final_use.rds")
#tum = read_rds("mye_final_use.rds")


DimPlot(b, group.by = "final_group", cols = colors_dark, label = T)
Idents(b) <- "final_group"
mks = FindAllMarkers(b, assay = "RNA")


top_markers = mks %>% arrange(cluster) %>% 
  dplyr::group_by(cluster) %>%
  top_n(7, avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()

Idents(clus_out) <- "integrated_snn_res.0.6"
tum_sub1 = AverageExpression(clus_out, features = unique(top_markers), return.seurat = T, assays = "RNA")
p  = DoHeatmap(tum_sub1, top_markers, group.colors = col_m, draw.lines = F)
png(glue("{save_path}/t_htmap1.png"), width=6, height=15, units = "in", res = 200)
print(p)
dev.off()

FeaturePlot(clus_out, "rna_PLAC8")
FeaturePlot(clus_out, "rna_TGFBR2", order = T)

tum.integrated.ns = clus_out
tum.integrated.ns$final_group = "NULL"
tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c("0")] = "CD4_Th"
tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c("1")] = "CD4_Th17"
tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c("2")] = "CD8_Trm" 
tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c("3")] = "NK"
tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c("4")] = "CD8_Trm"

tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.4 %in% c("7")] = "NKT" #6
tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.4 %in% c("5")] = "CD8_Tem" #7
tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c("17")] = "CD8_Tem" #7

tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c("5", "15")] = "CD4_naive"
tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c("10")] = "NK/ILCs"
tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c("11", "13", "20", "9")] = "CD4_Tem"
tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c("12")] = "CD4_Treg"
tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c("14")] = "CD8_Trm"
tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c("18")] = "GD?"
tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c("19")] = "t_prol"

tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c("8") & tum.integrated.ns$CD4CD8 %in% c("CD8")] = "CD8_znf683"
tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c("8") & tum.integrated.ns$CD4CD8 %in% c("CD4")] = "CD4_PLCG2"

tum.integrated.ns$final_group[tum.integrated.ns$final_group %in% c("NULL")] = tum.integrated.ns$cellstate
tum.integrated.ns$final_group[tum.integrated.ns$final_group %in% c("CD4_naive_like")] = "CD4_naive"
tum.integrated.ns$final_group[tum.integrated.ns$final_group %in% c("CD4_Teff_memory")] = "CD4_Tem"
tum.integrated.ns$final_group[tum.integrated.ns$final_group %in% c("CD4_Th")] = "CD4_Th1"
tum.integrated.ns$final_group[tum.integrated.ns$final_group %in% c("CD4_Th_like")] = "CD4_Th1like"
tum.integrated.ns$final_group[tum.integrated.ns$final_group %in% c("CD4-CD8_Tresident")] = "CD4_PLCG2"
tum.integrated.ns$final_group[tum.integrated.ns$final_group %in% c("CD8_TRM_like")] = "CD8_Trm"
tum.integrated.ns$final_group[tum.integrated.ns$final_group %in% c("CD8_Teff_memory")] = "CD8_Tem"


tum.integrated.ns$final_group1 = factor(tum.integrated.ns$final_group, levels = c("CD4_Th1", "CD4_Th1like", "CD4_naive",
                                                                                  "CD4_Tem", "CD4_Treg", "CD4_PLCG2",
                                                                                  "CD8_znf683", "CD8_Trm",
                                                                                  "CD8_Tem",  "NKT", "NK/ILCs",
                                                                                  "NK", "GD?", "t_prol"))


# DimPlot(tum.integrated.ns, group.by = "final_group", cols = colors_dark, label = T)
# table(tum.integrated.ns$final_group1)
# table(is.na(tum.integrated.ns$final_group))

write_rds(tum.integrated.ns, "t_final_v1.rds")
tum.integrated.ns = read_rds("t_final_v1.rds")

# from umap

# **Umap ------
tum = tum.integrated.ns
umap_tx = tum@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(tum_meta = tum@meta.data)

umap_tx  =  umap_tx %>% mutate(major_group = factor(tum_meta.final_group1, levels = c("CD4_Th", "CD4_Thlike", "CD4_naive",
                                                                                      "CD4_Tem", "CD4_Treg", "CD4_PLCG2",
                                                                                      "CD8_znf683", "CD8_Trm",
                                                                                      "CD8_Tem",  "NKT", "NK/ILCs",
                                                                                      "NK", "GD?", "t_prol")))

fill_col = c("#C70E7B", "#FC6882", "#007BC3", "#54BCD1", "#EF7C12", "#F4B95A", "#009F3F", "#8FDA04", "#AF6125", "#F4E3C7", "#B25D91", "#EFC7E6",
                     "red", "black", colors_dark)


library(ggpubr);library(grid)
grob <- grobTree(textGrob("6977 cells", x=0.7,  y=0.95, hjust=0,
                          gp=gpar(col="black", fontsize=10, fontface="italic")))
umap_tx1 <- umap_tx[sample(x = 1:nrow(x = umap_tx)), ]
p1 = ggplot(umap_tx1, aes(x=UMAP_1, y=UMAP_2)) + #annotation_custom(grob) + 
  #geom_point(shape = 16, size = 0.005, aes(color = major_group)) + 
  geom_jitter(shape = 21, size = 0.5, aes(fill = major_group),  color = "black", alpha = 0.9, stroke = 0.025) + 
  #geom_jitter(shape = 21, size = 1, aes(fill = major_group),  color = "black", alpha = 0.9, stroke = 0.04) + 
  #geom_point(shape = 19, size = 0.5, aes(color = major_group)) + 
  scale_fill_manual(values = fill_col) +
  #ggthemes::calc_pal()(12) +
  #rcartocolor::scale_fill_carto_d(palette = 2) +
  #scale_color_manual(values = fig1_colors) +
  #guides(fill = guide_legend(override.aes = list(color = fill_col)),
  #color = guide_legend(override.aes = list(shape = 21))) +
  xlab("UMAP-1") + ylab("UMAP-2") +
  theme(
    panel.grid = element_blank(), 
    #       panel.background = element_rect(fill = NULL, colour = NULL),
    panel.background = element_blank(),
    legend.position = "none",
    #axis.title = element_blank(),
    #axis.line = element_line(arrow=arrow)
    axis.line.x = element_line(colour = "black"), #,arrow=arrow
    axis.line.y = element_line(colour = "black"),#,arrow=arrow
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),axis.text.y = element_blank());p1
p_fixed <- set_panel_size(p1,
                          width  = unit(1.5, "in"),
                          height = unit(1.5, "in"))

ggsave(plot = p_fixed, filename = "t_umap.png", device = "png", path = save_path, width = 3, height = 3, dpi = 300, units = "in")
ggsave(plot = p_fixed, filename = "t_umap.pdf", device = "pdf", path = save_path, width = 3, height = 3, dpi = 300, units = "in")

p03 = DimPlot(object = tum, reduction = "umap", label = F, pt.size = 0.000001, cols = fill_col, group.by = "final_group1", shuffle = T)+  NoAxes()
pdf(glue("{save_path}/t_uamp_1.pdf"), width=6, height=5)
print(plot_grid(p03))
dev.off()

# **markers ----
# FINAL OBJECT 
# tum = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/figure_6_v1/mye_final_use.rds")
# tum = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_6_v1/mye_final_use.rds")
# tum$final_group = 
#   tum$final_group = factor(tum$final_group, levels = levels(umap_tx$major_group))
Idents(tum) = "final_group1"

all_markers = FindAllMarkers(object = tum, assay = "RNA", logfc.threshold = 0.25, only.pos = TRUE)
write.csv(all_markers, paste0(save_path, "t_fig_markers_final_group1.csv"))
# all_markers = read.csv("t_fig_markers_final_group1.csv")
# all_markers = all_markers[,-1]

allcell_markers = FindAllMarkers(object = tum, assay = "RNA", logfc.threshold = 0, only.pos = FALSE)
write.csv(allcell_markers, paste0(save_path, "t_allcell_markers.csv"))
allcell_markers = read.csv("t_allcell_markers.csv")
allcell_markers = allcell_markers[,-1]


# **Heatmap ------
all_markers_sel2 = all_markers
lvl = levels(umap_tx$major_group)
all_markers_sel2$cluster = factor(all_markers_sel2$cluster, levels = lvl)


top_markers = all_markers_sel2 %>% arrange(cluster) %>% 
  dplyr::group_by(cluster) %>%
  filter(!str_detect(gene, "^RP|^HSP")) %>% 
  dplyr::filter(avg_logFC > 0) %>%
  dplyr::slice_max(n = 7, order_by = avg_logFC) %>%
  pull(gene) %>% as.character()

top_markers = top_markers[!top_markers %in% c("IGKC", "SYTL3")]

genes_sel_df = all_markers_sel2 %>% arrange(cluster) %>% 
  group_by(cluster) %>% 
  filter(!str_detect(gene, "^RP|^HSP")) %>% 
  dplyr::filter(avg_logFC > 0) %>%
  dplyr::slice_max(n = 7, order_by = avg_logFC) %>%
  dplyr::select(gene, cluster) %>% mutate(id = paste0(gene, "-", cluster))

genes_sel_df = genes_sel_df %>% dplyr::filter(!gene %in% c("IGKC", "SYTL3"))

DefaultAssay(tum) <- "RNA"
Idents(tum) <- "final_group1"
tum_sub1 = AverageExpression(object = tum, features = unique(top_markers), return.seurat = T, assays = "integrated")
levels(tum_sub1$orig.ident) <- levels(Idents(tum))

library(ggpubr)
p_load(ggthemr)
p1 = DoHeatmap(tum_sub1, features = top_markers, assay = "integrated", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = fill_col, raster = T)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5)) +
  rotate_x_text(angle = 90) + coord_flip() + 
  scale_y_discrete(limits=rev) +   scale_x_discrete(limits=rev) 
ggsave(plot = p1, filename = "heatmap_cellstate_tcells_v2.pdf", width = 10, height = 2.5)


main_markers = c("CD4", "CD8B")
Groups_set <- levels(umap_tx$major_group)
col_grp <- fill_col[1:14]
names(col_grp) <- Groups_set

df_fibro1 = run_mean_htmaps_fig6_horizontal(df = all_markers_sel3, seu = tum, cell = "t_cells", genes_df = top_markers, #sel_genes
                                            genes_df_cann = main_markers, genes_sel_df = genes_sel_df, 
                                            clus_levels = levels(umap_tx$major_group))



# **Barplots ---------------------------------------------------------------
tum = read_rds("t_final_v1.rds")
df_meta = data.frame(tum@meta.data, 
                     celltype = tum$final_group1, 
                     stringsAsFactors = F) %>% 
  as_tibble() %>% clean_names()

df_meta1 = df_meta
# create a summary table
func_bar(var_df = df_meta1, var_x = "sheet_patient_id", var_fill = "celltype", var_pos = "fill", var_cols = col_grp, celltype = "t_cells")
func_bar(var_df = df_meta1, var_x = "sheet_figure1_grp_updated", var_fill = "celltype", var_pos = "fill", var_cols = col_grp, celltype = "t_cells")

df_meta2 = df_meta1 %>% dplyr::filter(sheet_tissue_location %in% c("left", "right"), !sheet_patient_id == "pt45")
df_meta3 = df_meta2 %>% mutate(sheet_patient_id =  factor(sheet_patient_id, levels = unique(sheet_patient_id[order(sheet_patient_id, sheet_tissue_location)]), ordered = TRUE))
df_meta3 = df_meta3 %>% mutate(key = paste0(sheet_patient_id, "_", sheet_tissue_location))

p1 = ggplot(df_meta3, aes(sheet_tissue_location, fill = celltype)) +
  geom_bar(position = "fill", width = 1) + 
  xlab("") + ylab("") + scale_fill_manual(values = col_grp) +
  theme(axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
        panel.background = element_blank(),axis.title.x = element_blank(), 
        axis.text.x = element_blank(), axis.ticks = element_blank(), 
        legend.title = element_blank()) + facet_wrap(~sheet_patient_id, ncol = 22, strip.position = "bottom") + 
  theme(strip.background = element_blank(), panel.spacing.x=unit(0.2, "lines"), strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(2,0,0,0, "pt")))
ggsave(p1, file = paste0(save_path, "/figure6_b_tissue_location_patient_t.pdf"), width = 7, height = 2.5)

func_bar(var_df = df_meta3, var_x = "key", var_fill = "celltype", var_pos = "fill", var_cols = col_grp, celltype = "t_cells_leftright")
func_bar(var_df = df_meta3, var_x = "sheet_tissue_location", var_fill = "celltype", var_pos = "fill", var_cols = col_grp, celltype = "t_cells_leftright")

df_meta2 = df_meta1 %>% dplyr::filter(sheet_bmi_final %in% c("normal","Normal","obese", "overweight"))
df_meta2 = df_meta2 %>% dplyr::mutate(bmi = case_when(sheet_bmi_final %in% c("normal", "Normal") ~ "normal",
                                                      sheet_bmi_final %in% c("obese") ~ "obese",
                                                      sheet_bmi_final %in% c("overweight") ~ "overweight",
                                                      sheet_bmi_final %in% c("unknown") ~ "unknown",
                                                      is.na(sheet_bmi_final)  ~ "unknown"))
df_meta2$bmi = factor(df_meta2$bmi, levels = c("normal", "overweight", "obese"), ordered = T) 
df_meta3 = df_meta2 %>% mutate(sheet_patient_id =  factor(sheet_patient_id, levels = unique(sheet_patient_id[order(sheet_patient_id, bmi)]), ordered = TRUE))
df_meta3 = df_meta3 %>% mutate(key = paste0(sheet_patient_id, "_", bmi))

func_bar(var_df = df_meta3, var_x = "bmi", var_fill = "celltype", var_pos = "fill", var_cols = col_grp, celltype = "t_cells_bmi_average")

df_meta2 = df_meta2 %>% dplyr::filter(!sheet_figure1_grp_updated %in% c("unknown"), !is.na(sheet_figure1_grp_updated))
df_meta2$sheet_figure1_grp_updated = factor(df_meta2$sheet_figure1_grp_updated, levels = c("Reduction Mammoplasty", "Prophylatic Mastectomy", "Cancer Mastectomy"), ordered = T)
df_meta3 = df_meta2 %>% mutate(sheet_patient_id =  factor(sheet_patient_id, levels = unique(sheet_patient_id[order(sheet_patient_id, sheet_figure1_grp_updated)]), ordered = TRUE))
df_meta3 = df_meta3 %>% mutate(key = paste0(sheet_patient_id, "_", sheet_figure1_grp_updated))

func_bar(var_df = df_meta3, var_x = "sheet_figure1_grp_updated", var_fill = "celltype", var_pos = "fill", var_cols = col_grp, celltype = "t_cells_sheet_figure1_grp_updated_average")

# **Clinical ---------
df_meta2_sel = df_meta1 %>% dplyr::select(celltype, final_group1, final_sample_id, ethnicity, age, sheet_patient_id, sheet_source, sheet_exp_proc, sheet_age_updated, sheet_brca_updated, sheet_parity, sheet_gravida, sheet_menopause, sheet_density, sheet_bmi_final, sheet_cancer_risk,sheet_figure1_grp_updated)
colnames(df_meta2_sel) = c("cellgroup", "celltype", "sample_id", "ethnicity", "age_orig", "patient_id",  "source", "exp_proc", "age", "brca", "parity", "gravida", "menopause", "density", "bmi", "cancer_risk", "surgery")
# df_mrg = df_meta2_sel %>% mutate(pt_x = paste0(patient_id, "_", type))
write_rds(df_meta2_sel, "t_df_cell_meta2_sel.rds")

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


df_mrg_nuc1 = df_meta2_sel %>% group_by(sample_id) %>% mutate(total_cells_per_samp = n()) %>% ungroup() %>% 
  group_by(sample_id, celltype) %>% mutate(cell_per_samp = n(), prop_cell_per_samp = cell_per_samp/total_cells_per_samp) %>% ungroup() %>% 
  group_by(sample_id, cellgroup) %>% mutate(cell_per_grp_samp = n(), prop_cell_per_grp_samp = cell_per_samp/cell_per_grp_samp) %>% ungroup()

df_mrg_nuc2 = df_meta2_sel %>% group_by(patient_id) %>% mutate(total_cells_per_pt = n()) %>% ungroup() %>% 
  group_by(patient_id, celltype) %>% mutate(cell_per_pt = n(), prop_cell_per_pt = cell_per_pt/total_cells_per_pt) %>% ungroup() %>% 
  group_by(patient_id, cellgroup) %>% mutate(cell_per_grp_pt = n(), prop_cell_per_grp_pt = cell_per_pt/cell_per_grp_pt) %>% ungroup()


write.csv(df_mrg_nuc1, "df_mrg_nuc1_samp_T.csv")
write.csv(df_mrg_nuc2, "df_mrg_nuc2_pt_T.csv")

save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_6_v1/B/"
col_pal = col_grp
do_clinical_comp(df_mrg_nuc2 = df_mrg_nuc2)

# dotplots exhaustion ------------
exh_genes = c("CTLA4", "TIGIT", "LAG3",  "CD27", "ICOS", "TNFRSF4", "CD40LG", "TNFRSF18", "CD226", "BTLA", "VSIR", #"CD96", 
"TIGIT", "LAG3", "PDCD1", "CD28", "CTLA4", "HAVCR2", "TNFRSF14")
pal <- wes_palette("Zissou1", 5, type = "continuous")
DefaultAssay(tum.integrated.ns) = "RNA"
p12 = DotPlot(tum.integrated.ns, features = unique(as.character(exh_genes)), group.by = "final_group1", assay = "RNA") +
  scale_color_gradientn(colours = pal) + coord_flip() 
p13 = p12 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
                  axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size = 6)) 
ggsave(plot = p13, filename = "figure_dotplot_exhaustion.pdf", width = 5, height = 5)

# dotplots checkpoint ------------
chk_genes = c("CD27", "CD40LG", "TNFRSF4", "ICOS", "TNFRSF18",
               "CD28",
               "PDCD1", "CTLA4", "HAVCR2", "BTLA", "LAG3", "CD226", "VSIR")
pal <- wes_palette("Zissou1", 5, type = "continuous")
DefaultAssay(tum.integrated.ns) = "RNA"
p12 = DotPlot(tum.integrated.ns, features = unique(as.character(chk_genes)), group.by = "final_group1", assay = "RNA") +
  #scale_color_gradientn(colours = pal) + 
  scale_color_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") +
  rotate_x_text(angle = 90)  + 
  scale_y_discrete(limits=rev) +   scale_x_discrete(limits=rev) 
  #+ coord_flip() 
p13 = p12 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
                  axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size = 6)) 
ggsave(plot = p13, filename = "figure_dotplot_checkpoint.pdf", width = 6.5, height = 2.3)

# dotplots CANNONICAL ------------
can_genes = c("CD8A","SELL","CD3D","CD4","IL7R","FOXP3","PLCG2","CD2","CD52","ITGA1","GNLY","TRDC","GZMK","MKI67")
DefaultAssay(tum.integrated.ns) = "RNA"
p12 = DotPlot(tum.integrated.ns, features = unique(as.character(can_genes)), group.by = "final_group1", assay = "RNA",dot.scale = 2.5) +
  scale_color_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") +
  rotate_x_text(angle = 90)  + 
  scale_y_discrete(limits=rev) +   scale_x_discrete(limits=rev) 
p13 = p12 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
                  axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size = 6)) 
ggsave(plot = p13, filename = "figure_dotplot_cann.pdf", width = 6.5, height = 2)

# combine bar plots ----------------
df_meta_b1 = df_meta_b %>% dplyr::select(patient_id, celltype)
df_meta_myeloid1 = df_meta_myeloid %>% dplyr::select(patient_id, celltype)
df_meta_t1 = df_meta_t %>% dplyr::select(patient_id, celltype)
df_meta_all = rbind(df_meta_b1, df_meta_myeloid1, df_meta_t1)
df_meta_all = df_meta_all %>% dplyr::mutate(immunegroup = case_when(celltype %in% c("b_naive","bmem_switched", "bmem_unswitched","plasma_IgA","plasma_IgG") ~ "B",
                                                                    celltype %in% c("cDC1","cDC2","mDC","pDC","mast","mono_naive","mono_active","macro_lipo",
                                                                                    "macro_m1","macro_m2","macro_ifn_m2","mye_prol") ~ "Myeloid",
                                            celltype %in% c("CD4_Th1","CD4_Th1like", "CD4_naive","CD4_Tem","CD4_Treg","CD4_PLCG2","CD8_znf683","CD8_Trm","CD8_Tem","NKT",
                                            "NK", "NK/ILCs", "GD?","t_prol") ~ "NK/T"))

func_bar(var_df = df_meta_all, var_x = "patient_id", var_fill = "immunegroup", var_pos = "fill", var_cols = c("#00B0F0FF", "#F8A19FFF", "#B10DA1FF"), celltype = "immunegroup")

df_imm_non_immune = data.frame(immune = 63815, non_immune = 535941, type = "cell")
df_imm_non_immune$type = as.character(df_imm_non_immune$type)
df_imm_non_immune[2, ] = c(17289, 102735, 'nuclei')

ggplot(data_long, aes(fill=variable, y=value, x=type)) + 
  geom_bar(position="dodge", stat="identity")

data_long <- reshape2::melt(df_imm_non_immune, id.vars = "type")
data_long$value = as.numeric(data_long$value)
ggplot(data_long,aes(fill=variable, y=value, x=type)) +
  geom_bar(position = "fill", stat = "identity") + 
  xlab("") + ylab("") + scale_fill_manual(values = c("#00B0F0FF", "#F8A19FFF", "#B10DA1FF")) +
  theme(axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
        panel.background = element_blank(),axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(2,0,0,0, "pt")), axis.ticks = element_blank(), 
        legend.title = element_blank())

