# geo ----------
setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/")
odir = "figure_1_v1";dir.create(odir)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")

p_load(future)
plan("multiprocess", workers = 40)
options(future.globals.maxSize = 90000 * 1024^2)

data_folder = "04_final_clustering"
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/figure_1_v1/"
setwd(save_path)

# core --------
setwd("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/")
odir = "figure_1_v1";dir.create(odir)
setwd(paste0("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")

p_load(future)
plan("multiprocess", workers = 40)
options(future.globals.maxSize = Inf)

data_folder = "04_final_clustering"
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_1_v1/"
setwd(save_path)

# ************************* -----
# **cell data figure -----

# read data ----
# 15, 30
tum = readRDS(file = paste0(wd, data_folder,"/final3.rds")) #dim_15k_30_NULL_tum.integrated.ns2_method_9rcp_0.1_1.5final.rds
# tum1 = tum
tum1 = RunUMAP(tum, spread = 1.85, dims = 1:15, n.neighbors = 30, min.dist = 0.1)
# DimPlot(tum1)
write.rds(tum1, "tum1_85.rds")

library(cowplot)
p00 = DimPlot(object = tum, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = colors_dark, group.by = "integrated_snn_res.1", shuffle = T)+  NoAxes()
p01 = DimPlot(object = tum, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = colors_dark, group.by = "integrated_snn_res.0.6", shuffle = T)+  NoAxes()
p02 = DimPlot(object = tum, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = colors_dark, group.by = "major_group", shuffle = T)+  NoAxes()
p03 = DimPlot(object = tum, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = colors_dark, group.by = "integrated_snn_res.0.4", shuffle = T)+  NoAxes()
png(glue("{save_path}/uamp_1.png"), width=16, height=15, units = "in", res = 300)
print(plot_grid(p00, p01, p02, p03))
dev.off()

# rename clusters -----
tum.integrated.ns = tum
tum.integrated.ns$final_group = "NULL"
tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c(0, 2, 8, 11, 16, 12, 20, 22)] = "fibro"
tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c(1)] = "t"
tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c(4, 19)] = "epi_lumhr" 
tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c(3)] = "epi_basal"
tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c(5,18)] = "epi_lumsec"
tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c(6,10,17)] = "endo_vas"
tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c(15)] = "endo_lymph"
tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c(7)] = "pericytes"
tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c(13)] = "b"
tum.integrated.ns$final_group[tum.integrated.ns$integrated_snn_res.0.6 %in% c(14,9,21)] = "macro"
tum.integrated.ns$final_group1 = factor(tum.integrated.ns$final_group, levels = c("epi_basal", "epi_lumhr", "epi_lumsec",
                                                                      "fibro", "endo_lymph", "endo_vas",
                                                                      "pericytes", "macro",
                                                                      "t",  "b"), labels = c("Basal", "LumHR", "LumSec",
                                                                                                     "Fibroblasts", "Lymphatic", "Vascular",
                                                                                                     "Pericytes", "Myeloid",
                                                                                                     "T-cells",  "B-cells"))

write_rds(tum.integrated.ns, "final_use.rds")
tum = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/figure_1_v1/final_use.rds")



library(cowplot)
p00 = DimPlot(object = tum, reduction = "umap", label = FALSE, pt.size = 0.01, cols = colors_dark, group.by = "final_group1", shuffle = FALSE)
png(glue("{save_path}/uamp_3.png"), width=16, height=14, units = "in", res = 300)
print(p00)
dev.off()

source("densmap.R")
dm <- reticulate::import('densmap')
#data <- dm$load_trial$load_trial()
#densMAP <- dm$densMAP("final_dens"=FALSE)
#mdm.full <- densMAP$fit_transform(data)
#mdm.full
out <- dm$densMAP(as.matrix(t(tum@assays$integrated@data)), min_dist = 0.1, spread = 1, n_epochs=500, dens_frac=0.3, dens_lambda=0.5)
out1  = out
umap_tx1 = out[[1]] %>% 
  as.data.frame()
out = read_rds("out.rds")

# convert to table -----
# from densmap
library(huxtable)
umap_tx = out[[1]] %>% 
  as.data.frame() 
colnames(umap_tx)  = c("UMAP_1", "UMAP_2")
umap_tx = umap_tx %>% cbind(tum_meta = tum@meta.data)


# from umap
umap_tx = tum@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(tum_meta = tum@meta.data)

# set colors ----
# paletteer::paletteer_d("pals::alphabet2", 26)
# paletteer::paletteer_d("Redmonder::qMSOStd", 10)
# cc = paletteer::paletteer_d("Redmonder::qMSOStd", 10)[1:10]

fig1_colors = c("#C00000FF", "#FF0000FF", "#FFC000FF", "#FFFF00FF", "#92D050FF", "#00B050FF", "saddlebrown","#F8A19FFF", "#B10DA1FF","#00B0F0FF", "#0070C0FF", "#002060FF", "#7030A0FF")
fig1_colors2 = c("#971F2A","#E2734E","#F1AF71","#6E3380","#668C36","#C0E093","#3866A7","#B3337A","#E95FA5", "#E8B8DA", "black")
#co = c("#f9f871", "#00c9ae", "#ff8863", "#6e3452", "#00aafd", "#b178de", "#ff86da", "#7d9300","red2");show_col(co)

umap_tx  =  umap_tx %>% mutate(major_group = factor(tum_meta.final_group1, levels = c("Basal", "LumHR", "LumSec",
                                                                                      "Fibroblasts", "Lymphatic", "Vascular",
                                                                                      "Pericytes", "Myeloid",
                                                                                      "T-cells",  "B-cells"), labels = c("Basal", "LumHR", "LumSec",
                                                                                                                         "Fibroblasts", "Lymphatic", "Vascular",
                                                                                                                         "Pericytes", "Myeloid",
                                                                                                                         "T-cells",  "B-cells")))

# plot ggplot
# figure 1b umap -----
fill_col = fig1_colors
library(ggpubr);library(grid)
grob <- grobTree(textGrob("535941 cells", x=0.7,  y=0.95, hjust=0,
                          gp=gpar(col="black", fontsize=10, fontface="italic")))
p1 = ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) + #annotation_custom(grob) + 
  #geom_point(shape = 16, size = 0.005, aes(color = major_group)) + 
  #geom_point(shape = 21, size = 1, aes(fill = major_group), alpha = 0.5, stroke = 0.05, position = position_jitter(w = 0.1, h = 0.1)) + #, position = "dodge"
  geom_jitter(shape = 21, size = 2, aes(fill = major_group),  color = "black", alpha = 0.5, stroke = 0.2, position = "dodge") + 
  scale_fill_manual(values = fig1_colors) + 
  #scale_color_manual(values = fig1_colors) +
  #guides(fill = guide_legend(override.aes = list(color = fill_col)),
         #color = guide_legend(override.aes = list(shape = 21))) +
  xlab("UMAP-1") + ylab("UMAP-2") +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.border = element_blank(),
        legend.position = "none",
        #axis.title = element_blank(),
        #axis.line = element_line(arrow=arrow)
        axis.line.x = element_line(colour = "black"), #,arrow=arrow
        axis.line.y = element_line(colour = "black"),#,arrow=arrow
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank());p1

ggsave(plot = p1, filename = "figure1_cell_umap.png", device = "png", path = save_path, width = 3.5, height = 3, dpi = 300, units = "in")
ggsave(plot = p1, filename = "figure1_cell_umap_large1.png", device = "png", path = save_path, width = 13, height = 13, dpi = 300, units = "in")
ggsave(plot = p1, filename = "figure1_cell_umap_largedots.png", device = "png", path = save_path, width = 3.5, height = 3, dpi = 300, units = "in")
ggsave(plot = p1, filename = "figure1_cell_umap_large_largedots3.png", device = "png", path = save_path, width = 18, height = 18, dpi = 300, units = "in")
#ggsave(plot = p1, filename = "figure1_umap.pdf", device = "pdf", path = save_path, width = 3, height = 3)
#cairo_ps(filename = "figure1_umap.eps", width = 3.5, height = 3, fallback_resolution = 300); p1; dev.off()

# figure s1c -----
p2 = umap_tx %>%  
  mutate(sampleid2 = fct_reorder(tum_meta.sample_id, as.numeric(factor(tum_meta.figure1_grp_updated)))) %>% 
  ggplot(aes(x = sampleid2, fill = major_group)) +
  geom_bar(position = "fill", width = 1) + 
  #coord_polar() +
  ylab("Cellular Proportions") + xlab("Samples") +  
  scale_fill_manual(values = fig1_colors)+
  theme(panel.background = element_blank()) +
  theme(legend.title = element_blank(), 
        legend.key.size = unit(0.5,"line"), 
        plot.background = element_blank() ,
        #axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank());p2
ggsave(plot = p2, filename = "figure_s1c_cell_barplot_groups.png", device = "png", path = save_path, width = 3.5, height = 3, dpi = 300, units = "in")
ggsave(plot = p2, filename = "figure_barplot_fill.pdf", device = "pdf", path = save_path, width = 5, height = 3)
cairo_ps(filename = "figure_barplot_fill.eps", width = 5, height = 3,fallback_resolution = 300); p2; dev.off()

p2_header = umap_tx %>%  
  mutate(sampleid2 = fct_reorder(tum_meta.sample_id, as.numeric(factor(tum_meta.figure1_grp_updated)))) %>% 
  ggplot(aes(x = sampleid2, fill = tum_meta.figure1_grp_updated)) +
  geom_bar(position = "fill", width = 1) + 
  #coord_polar() +
  ylab("Cellular Proportions") + xlab("Samples") +  
  scale_fill_manual(values = beyonce_palette(78)[2:4])+
  theme(panel.background = element_blank()) +
  theme(legend.title = element_blank(),
        legend.direction = "vertical",
        legend.key.size = unit(0.5,"line"), 
        plot.background = element_blank() ,
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank());p2_header
ggsave(plot = p2_header, filename = "figure_s1c_cell_barplot_groups_header.png", device = "png", path = save_path, width = 3.5, height = 2, dpi = 300, units = "in")

# markers ----
# FINAL OBJECT -----
tum = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/figure_1_v1/final_use.rds")
tum = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_1_v1/final_use.rds")
Idents(tum) = "final_group1" #major_group
table(tum$final_group1)

all_markers = FindAllMarkers(object = tum, assay = "RNA", logfc.threshold = 0.25, only.pos = TRUE)
write.csv(all_markers, paste0(save_path, "fig_markers_final_group1.csv"))
all_markers = read.csv("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/figure_1_v1/fig_markers_final_group1.csv")
all_markers = read.csv("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_1_v1/fig_markers_final_group1.csv")
all_markers = all_markers[,-1]

allcell_markers = FindAllMarkers(object = tum, assay = "RNA", logfc.threshold = 0, only.pos = FALSE)
write.csv(allcell_markers, paste0(save_path, "allcell_markers.csv"))
allcell_markers = read.csv("allcell_markers.csv")
allcell_markers = allcell_markers[,-1]

# markers for heatmap

## run function to select markers ----
allcell_markers %>%
  mutate(cluster1 = forcats::fct_relevel(cluster, levels = c("Basal", "LumHR", "LumSec",
                                                             "Fibroblasts", "Lymphatic", "Vascular",
                                                             "Pericytes", "Myeloid",
                                                             "T-cells", "B-cells"))) %>% 
  arrange(cluster1) %>% 
  filter(!str_detect(gene, "^RP")) %>% 
  filter(!str_detect(gene, "^MT-")) %>% 
  filter(!str_detect(gene, "^HBB")) -> all_markers_sel

df_thresh = data.frame(
  cluster = c("Basal", "LumHR", "LumSec", 
              "Fibroblasts", "Lymphatic", "Vascular","Pericytes",
              "Myeloid", "T-cells", "B-cells"),
  thresh_pct1 = c(0.40, 0.40, 0.40, 
                  0.20, 0.20, 0.20, 0.20,
                  0.30, 0.30, 0.30), stringsAsFactors = F
) %>% data.frame(row.names = .$cluster)

p_load(magrittr)
all_markers_sel %<>% as_tibble()

all_markers_sel2 = func_select_markers_thresh(all_markers_sel = all_markers_sel, pct_var = 0.25)

view(all_markers_sel2)

write.csv(all_markers_sel2, "all_markers_sel2_cells.csv")
all_markers_sel2 %>%
  #filter(pct.2 < 0.25) %>% 
  # filter(pct.1 > 0.50) %>% 
   group_by(cluster1) %>%
  #arrange(avg_logFC)
  #top_n(50, avg_logFC) %>%
  top_n(10, avg_logFC) %>% 
  pull(gene) %>% as.character() ->
  #tum_markers_top
   tum_markers_top_rna

all_markers_sel2 %>%
  # filter(pct.2 < 0.30) %>% 
  # filter(pct.1 > 0.75) %>% 
  group_by(cluster1) %>% 
  #top_n(50, avg_logFC) %>%
  top_n(10, avg_logFC) %>% 
  dplyr::select(gene, cluster1) %>% mutate(id = paste0(gene, "-", cluster1)) -> genes_sel_df
#main_markers = c("EPCAM", "DCN", "PECAM1", "RGS5", "PTPRC", "CD68", "CD3D", "CD79A")

data.frame(new_cluster_factor = tum@meta.data$final_group1, cell = tum$cellname) %>%
  #filter(!str_detect(cell, "untreated")) %>% 
  group_by(new_cluster_factor) %>% 
  sample_n(500) %>% 
  mutate(cell = as.character(cell)) %>% 
  pull(cell) ->
  tum_cells

# figure 1d markers ----
#main_markers = c("EPCAM", "DCN", "PECAM1", "RGS5", "PTPRC", "CD68", "CD3D", "NKG7", "CD79A")
DefaultAssay(tum) = "RNA"
p10 = DoHeatmap(tum, group.by = "final_group1", features = as.character(tum_markers_top_rna), cells = tum_cells, label=TRUE, raster = TRUE, 
                group.colors = fig1_colors, disp.min = -1.0, disp.max = 1, size = 2) + 
  theme(text = element_text(size = 8)) + 
  #scale_fill_gradient2(low="navyblue", mid="black",high="cornsilk", midpoint = -0.2)
  scale_fill_gradient2(low="steelblue4", mid="white",high="tomato4", midpoint = 0)
  #scale_fill_gradient2(low="dodgerblue", mid="white",high="maroon", midpoint = 0);p10
p10
ggsave(plot = p10, filename = paste0(save_path,"/figure_heatmap_cell_types_redblue100L_cells_integratedpct0.40_0.20_0.25_pct25_top15.png"), device = "png", width = 6, height = 12)


# dotplot ---
pal <- wes_palette("Zissou1", 5, type = "continuous")
p12 = DotPlot(tum, features = unique(as.character(tum_markers_top_rna)), dot.scale = 1.5, group.by = "major_group_factor") +
  scale_color_gradientn(colours = pal) + coord_flip() 
p13 = p12 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
            axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size = 6)) 
ggsave(plot = p13, filename = "figure_dotplot_cell_types.png", device = "png", width = 5, height = 6.5)
ggsave(plot = p13, filename = "figure_dotplot_cell_types.pdf", device = "pdf", width = 5, height = 6.5)

# big heatmap -----
Groups_set <- c("Basal", "LumHR", "LumSec",
  "Fibroblasts", "Lymphatic", "Vascular",
  "Pericytes", "Myeloid",
  "T-cells",  "B-cells")
col_grp <- fig1_colors[1:10]
names(col_grp) <- Groups_set

Protocol = c("short_digestion", "6hr_digestion", "overnight_digestion")
col_Protocol = c("bisque", "darkorange1", "deeppink2")
names(col_Protocol) <- Protocol

meno = c("pre", "post", "unknown")
col_meno = c("#EA879C", "#8E2322", "gray90")
names(col_meno) <- meno

paridy = c("0", "1", "unknown")
col_paridy = c("#CFA2D2", "#5D36B1", "gray90")
names(col_paridy) <- paridy

procedure = c("Reduction Mammoplasty", "Prophylatic Mastectomy", "Cancer Mastectomy")
col_procedure = c("#F5EC54", "#84C876", "#539E59")
names(col_procedure) <- procedure

df_fibro = run_htmaps(df = all_markers, seu = tum, cell = "all_cells", genes_df = tum_markers_top_rna, 
                      clus_levels = c("Basal", "LumHR", "LumSec",
                                      "Fibroblasts", "Lymphatic", "Vascular",
                                      "Pericytes", "Myeloid",
                                      "T-cells", "B-cells"))

df_meta = data.frame(seu@meta.data, 
                      celltype = seu$major_group_factor, 
                      stringsAsFactors = F) %>% tbl_df() %>% clean_names()
df_meta_meno = df_meta %>% dplyr::select(celltype, menopause_y)
df_meta_meno1 = table(df_meta_meno$menopause_y, df_meta_meno$celltype)
df_meta_meno2 = prop.table(df_meta_meno1, margin = 2)

df_meta_paridy = df_meta %>% dplyr::select(celltype, parity_1_y)
df_meta_paridy1 = table(df_meta_paridy$parity_1_y, df_meta_meno$celltype)
df_meta_paridy2 = prop.table(df_meta_paridy1, margin = 2)

df_meta_grp = df_meta %>% dplyr::select(celltype, figure1_grp_updated)
df_meta_grp1 = table(df_meta_grp$figure1_grp_updated, df_meta_grp$celltype)
df_meta_grp2 = prop.table(df_meta_grp1, margin = 2)

df_meta_dig = df_meta %>% dplyr::select(celltype, exp_proc)
df_meta_dig1 = table(df_meta_dig$exp_proc, df_meta_grp$celltype)
df_meta_dig2 = prop.table(df_meta_dig1, margin = 2)

df_meta_grp = df_meta %>% dplyr::select(celltype, figure1_grp_updated)
df_meta_grp1 = table(df_meta_grp$figure1_grp_updated, df_meta_grp$celltype)
df_meta_grp2 = prop.table(df_meta_grp1, margin = 2)


all_markers %>%
  mutate(cluster1 = forcats::fct_relevel(cluster, levels = c("Basal", "LumHR", "LumSec",
                                                             "Fibroblasts", "Lymphatic", "Vascular",
                                                             "Pericytes", "Myeloid",
                                                             "T-cells", "B-cells"))) %>% 
  arrange(cluster1) %>% 
  filter(!str_detect(gene, "^RP")) %>% 
  filter(!str_detect(gene, "^MT-")) %>% 
  filter(!str_detect(gene, "^HBB")) -> all_markers_sel_fig


df_fibro1 = run_mean_htmaps_fig1(df = all_markers_sel_fig, seu = tum, cell = "all_cells", genes_df = tum_markers_top_rna, #sel_genes
                                genes_df_cann = main_markers, genes_sel_df = genes_sel_df, 
                                clus_levels = c("Basal", "LumHR", "LumSec",
                                                "Fibroblasts", "Lymphatic", "Vascular",
                                                "Pericytes", "Myeloid",
                                                "T-cells",  "B-cells"))

df_fibro1 = run_mean_htmaps_fig1_horizontal(df = all_markers_sel_fig, seu = tum, cell = "all_cells", genes_df = tum_markers_top_rna, #sel_genes
                                 genes_df_cann = main_markers, genes_sel_df = genes_sel_df, 
                                 clus_levels = c("Basal", "LumHR", "LumSec",
                                                 "Fibroblasts", "Lymphatic", "Vascular",
                                                 "Pericytes", "Myeloid",
                                                 "T-cells",  "B-cells"))


# supp mean heatmap for fig 1
df_fibro1 = run_mean_htmaps_fig1_supp(df = all_markers_sel_fig, seu = tum, cell = "all_cells", genes_df = tum_markers_top_rna, #sel_genes
                                      genes_df_cann = main_markers, genes_sel_df = genes_sel_df, 
                                      clus_levels = c("Basal", "LumHR", "LumSec",
                                                      "Fibroblasts", "Lymphatic", "Vascular",
                                                      "Pericytes", "Myeloid",
                                                      "T-cells",  "B-cells"))
# supp horizontal mean heatmap for fig 1
df_fibro1 = run_mean_htmaps_fig1_horizontal_supp(df = all_markers_sel_fig, seu = tum, cell = "all_cells", genes_df = tum_markers_top_rna, #sel_genes
                                                 genes_df_cann = main_markers, genes_sel_df = genes_sel_df, 
                                                 clus_levels = c("Basal", "LumHR", "LumSec",
                                                                 "Fibroblasts", "Lymphatic", "Vascular",
                                                                 "Pericytes", "Myeloid",
                                                                 "T-cells",  "B-cells"))

# pathway ------
genes = read.csv("/volumes/lab/users/tapsi/projects/bcap/analysis/20200210_final1/figure_1_v2/allcell_markers.csv")
genes = genes[,-1]
genes$cluster = as.character(genes$cluster)

gsea_resultsc2hall = fun_gsea(genes = genes, set = "H", sub_cat = NULL, num_avg_logFC = 0.25)
gsea_resultsc2C5 = fun_gsea(genes = genes, set = "C5", sub_cat = "BP", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "fibro")
gsea_resultsc2C3 = fun_gsea(genes = genes, set = "C3", sub_cat = "TFT:GTRD", num_avg_logFC = 0.25, celltype = "lumsec", doplot = FALSE, ret_collap = FALSE)
gsea_resultsc2C8 = fun_gsea(genes = genes, set = "C8", sub_cat = NULL, num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "fibro")
gsea_resultsc2C2_reac = fun_gsea(genes = genes, set = "C2", sub_cat = "CP:REACTOME", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "fibro")
gsea_resultsc2C2_kegg = fun_gsea(genes = genes, set = "C2", sub_cat = "CP:KEGG", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "fibro")

clus = c("fibro-1", "fibro-2", "fibro-3")

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
colnames(p1.1) = c("pathway", "fibro-1", "fibro-2", "fibro-3")

# rearrange the datafram
test = p1

# for basal
test1 = test %>% dplyr::arrange(-`fibro-1`) %>% dplyr::filter(`fibro-1`>0) %>% pull(pathway)
test2 = test %>% dplyr::filter(!pathway %in% test1) %>% arrange(-`fibro-2`) %>% dplyr::filter(`fibro-1`<0, `fibro-2`>0, `fibro-3`<0) %>% dplyr::top_n(25, `fibro-2`) %>% pull(pathway)
test3 = test %>% dplyr::filter(!pathway %in% c(test1, test2)) %>% arrange(-`fibro-3`) %>% dplyr::filter(`fibro-1`<0, `fibro-2`<0, `fibro-3`>0) %>% dplyr::top_n(25, `fibro-3`) %>% pull(pathway)
seq_path = c(test1, test2, test3) #, test4, test5


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

png("/htmaps_pathways_kegg_reac_biocarta.png", res = 300, units = "in", width = 12, height = 12)
print(ha)
dev.off()

png("fibro_htmaps_pathways_GO_unique_ones_NES.png", res = 300, units = "in", width = 10, height = 24)
print(ht_hall)
dev.off()

png("htmaps_pathways_hallmark_NES.png", res = 300, units = "in", width = 6, height = 6)
print(ht_hall)
dev.off()

png("lum_hr_htmaps_pathways_kegg_NES.png", res = 300, units = "in", width = 6, height = 4)
print(ht_hall)
dev.off()

png("lum_hr_htmaps_pathways_reactome_NES_unique_top25.png", res = 300, units = "in", width = 10, height = 10)
print(ht_hall)
dev.off()

png("htmaps_pathways_immune_NES.png", res = 300, units = "in", width = 10, height = 20)
print(ht_hall)
dev.off()

png("lum_sec_htmaps_pathways_c3_NES.png", res = 300, units = "in", width = 10, height = 6)
print(ht_hall)
dev.off()

png("htmaps_pathways_atlas_NES.png", res = 300, units = "in", width = 10, height = 20)
print(ht_hall)
dev.off()


# **Barplots ---------------------------------------------------------------
df_meta = data.frame(tum@meta.data, 
                     celltype = tum$final_group1, 
                     stringsAsFactors = F) %>% 
  as_tibble() %>% clean_names()
write_rds(df_meta, "df_meta_cell.rds")
#df_meta = df_meta %>% mutate(orig_pt_ident = factor(orig_pt_ident, levels = c("HR+ Luminal Cells", "Secretory Luminal Cells", "Endothelial",  "Fibroblasts",  "Pericytes", "T cells", "Macrophages",
#                                                              "TumA", "TumB", "TumC", "TumD")))
#df_meta1 = df_meta %>% mutate(patient_id =  factor(patient_id, levels = unique(patient_id[order(source, procedure)]), ordered = TRUE))
# fn = factor(f, levels=unique(f[order(a,b,f)]), ordered=TRUE)
df_meta1 = df_meta
# create a summary table
func_bar(var_df = df_meta1, var_x = "sheet_patient_id", var_fill = "celltype", var_pos = "fill", var_cols = fig1_colors, celltype = "cells")
func_bar(var_df = df_meta1, var_x = "sheet_figure1_grp_updated", var_fill = "celltype", var_pos = "fill", var_cols = fig1_colors, celltype = "cells")

df_meta2 = df_meta1 %>% dplyr::filter(sheet_tissue_location %in% c("left", "right"), !sheet_patient_id == "pt45")
df_meta3 = df_meta2 %>% mutate(sheet_patient_id =  factor(sheet_patient_id, levels = unique(sheet_patient_id[order(sheet_patient_id, sheet_tissue_location)]), ordered = TRUE))
df_meta3 = df_meta3 %>% mutate(key = paste0(sheet_patient_id, "_", sheet_tissue_location))

p1 = ggplot(df_meta3, aes(sheet_tissue_location, fill = celltype)) +
  geom_bar(position = "fill", width = 1) + 
  xlab("") + ylab("") + scale_fill_manual(values = var_cols) +
  theme(axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
        panel.background = element_blank(),axis.title.x = element_blank(), 
        axis.text.x = element_blank(), axis.ticks = element_blank(), 
        legend.title = element_blank()) + facet_wrap(~sheet_patient_id, ncol = 22, strip.position = "bottom") + 
  theme(strip.background = element_blank(), panel.spacing.x=unit(0.2, "lines"), strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(2,0,0,0, "pt")))
ggsave(p1, file = paste0(save_path, "/figure1_tissue_location_patient.pdf"), width = 7, height = 2.5)

func_bar(var_df = df_meta3, var_x = "key", var_fill = "celltype", var_pos = "fill", var_cols = fig1_colors, celltype = "cells_leftright")
func_bar(var_df = df_meta3, var_x = "sheet_tissue_location", var_fill = "celltype", var_pos = "fill", var_cols = fig1_colors, celltype = "cells_leftright")

df_meta2 = df_meta1 %>% dplyr::filter(sheet_bmi_final %in% c("normal","Normal","obese", "overweight"))
df_meta2 = df_meta2 %>% dplyr::mutate(bmi = case_when(sheet_bmi_final %in% c("normal", "Normal") ~ "normal",
                                                      sheet_bmi_final %in% c("obese") ~ "obese",
                                                      sheet_bmi_final %in% c("overweight") ~ "overweight",
                                                      sheet_bmi_final %in% c("unknown") ~ "unknown",
                                                              is.na(sheet_bmi_final)  ~ "unknown"))
df_meta2$bmi = factor(df_meta2$bmi, levels = c("normal", "overweight", "obese"), ordered = T) 
df_meta3 = df_meta2 %>% mutate(sheet_patient_id =  factor(sheet_patient_id, levels = unique(sheet_patient_id[order(sheet_patient_id, bmi)]), ordered = TRUE))
df_meta3 = df_meta3 %>% mutate(key = paste0(sheet_patient_id, "_", bmi))

func_bar(var_df = df_meta3, var_x = "bmi", var_fill = "celltype", var_pos = "fill", var_cols = fig1_colors, celltype = "cells_bmi_average")

df_meta2 = df_meta2 %>% dplyr::filter(!sheet_figure1_grp_updated %in% c("unknown"), !is.na(sheet_figure1_grp_updated))
df_meta2$sheet_figure1_grp_updated = factor(df_meta2$sheet_figure1_grp_updated, levels = c("Reduction Mammoplasty", "Prophylatic Mastectomy", "Cancer Mastectomy"), ordered = T)
df_meta3 = df_meta2 %>% mutate(sheet_patient_id =  factor(sheet_patient_id, levels = unique(sheet_patient_id[order(sheet_patient_id, sheet_figure1_grp_updated)]), ordered = TRUE))
df_meta3 = df_meta3 %>% mutate(key = paste0(sheet_patient_id, "_", sheet_figure1_grp_updated))

func_bar(var_df = df_meta3, var_x = "sheet_figure1_grp_updated", var_fill = "celltype", var_pos = "fill", var_cols = fig1_colors, celltype = "cells_sheet_figure1_grp_updated_average")

# ************************* -----
# **nuclei data figure -----

# ************************* -----

# read data ----
# 15, 30
load("hbca_nuclei_final.RObj") #dim_15k_30_NULL_tum.integrated.ns2_method_9rcp_0.1_1.5final.rds
tum = hbca.nuclei.final
tum1 = RunUMAP(tum, spread = 1.85, dims = 1:15, n.neighbors = 30, min.dist = 0.1)
# DimPlot(tum1)
write.rds(tum1, "tum1_85.rds")

library(cowplot)
p00 = DimPlot(object = tum, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = colors_dark, group.by = "patient_id", shuffle = T)+  NoAxes()
p02 = DimPlot(object = tum, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = colors_dark, group.by = "cell_type_final", shuffle = T)+  NoAxes()
png(glue("{save_path}/nuclei_umap_1.png"), width=16, height=15, units = "in", res = 300)
print(plot_grid(p00, p02))
dev.off()

# rename clusters -----
tum.integrated.ns = tum
tum.integrated.ns$cell_type_final2 = Idents(tum.integrated.ns)
tum.integrated.ns$final_group = "NULL"
tum.integrated.ns$final_group[tum.integrated.ns$cell_type_final2 %in% c("adipo")] = "Adipocytes"
tum.integrated.ns$final_group[tum.integrated.ns$cell_type_final2 %in% c("endo_lymph")] = "Lymphatic"
tum.integrated.ns$final_group[tum.integrated.ns$cell_type_final2 %in% c("endo_vasc")] = "Vascular" 
tum.integrated.ns$final_group[tum.integrated.ns$cell_type_final2 %in% c("epi_basal")] = "Basal"
tum.integrated.ns$final_group[tum.integrated.ns$cell_type_final2 %in% c("epi_lumhr")] = "LumHR"
tum.integrated.ns$final_group[tum.integrated.ns$cell_type_final2 %in% c("epi_lumsec")] = "LumSec"
tum.integrated.ns$final_group[tum.integrated.ns$cell_type_final2 %in% c("fibro")] = "Fibroblasts"
tum.integrated.ns$final_group[tum.integrated.ns$cell_type_final2 %in% c("pericytes")] = "Pericytes"
tum.integrated.ns$final_group[tum.integrated.ns$cell_type_final2 %in% c("immune_mast")] = "Mast"
tum.integrated.ns$final_group[tum.integrated.ns$cell_type_final2 %in% c("immune_myeloid")] = "Myeloid"
tum.integrated.ns$final_group[tum.integrated.ns$cell_type_final2 %in% c("immune_t")] = "T-cells"
tum.integrated.ns$final_group1 = factor(tum.integrated.ns$final_group, levels = c("Basal", "LumHR", "LumSec",
                                        "Fibroblasts", "Lymphatic", "Vascular",
                                        "Pericytes", "Myeloid",
                                        "T-cells", "Mast", "Adipocytes"), labels = c("Basal", "LumHR", "LumSec","Fibroblasts", "Lymphatic", "Vascular","Pericytes", "Myeloid", "T-cells",  "Mast", "Adipocytes"))

write_rds(tum.integrated.ns, "final_nuclei_use.rds")

tum$patient_id[tum$sample_id %in% c("hbca_n04")] = "pt37"
tum = read_rds("final_nuclei_use.rds")
write_rds(tum, "final_nuclei_use.rds")


colors_dark_nuclei = colors_dark

library(cowplot)
p00 = DimPlot(object = tum, reduction = "umap", label = FALSE, pt.size = 0.01, cols = colors_dark, group.by = "final_group1", shuffle = FALSE)
png(glue("{save_path}/nuclei_uamp_3.png"), width=16, height=14, units = "in", res = 300)
print(p00)
dev.off()

source("densmap.R")
dm <- reticulate::import('densmap')
#data <- dm$load_trial$load_trial()
#densMAP <- dm$densMAP("final_dens"=FALSE)
#mdm.full <- densMAP$fit_transform(data)
#mdm.full
out <- dm$densMAP(as.matrix(t(tum@assays$integrated@data)), min_dist = 0.1, spread = 1, n_epochs=500, dens_frac=0.3, dens_lambda=0.5)
out1  = out
umap_tx1 = out[[1]] %>% 
  as.data.frame()
out = read_rds("out.rds")

# convert to table -----
# from densmap
library(huxtable)
umap_tx = out[[1]] %>% 
  as.data.frame() 
colnames(umap_tx)  = c("UMAP_1", "UMAP_2")
umap_tx = umap_tx %>% cbind(tum_meta = tum@meta.data)


# from umap
umap_tx = tum@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(tum_meta = tum@meta.data)

# set colors ----
# paletteer::paletteer_d("pals::alphabet2", 26)
# paletteer::paletteer_d("Redmonder::qMSOStd", 10)
# cc = paletteer::paletteer_d("Redmonder::qMSOStd", 10)[1:10]

fig1_colors = c("#C00000FF", "#FF0000FF", "#FFC000FF", "#FFFF00FF", "#92D050FF", "#00B050FF", "saddlebrown","#F8A19FFF", "#B10DA1FF", "blue", "black")
#co = c("#f9f871", "#00c9ae", "#ff8863", "#6e3452", "#00aafd", "#b178de", "#ff86da", "#7d9300","red2");show_col(co)

umap_tx  =  umap_tx %>% mutate(major_group = factor(tum_meta.final_group1, levels = c("Basal", "LumHR", "LumSec","Fibroblasts", "Lymphatic", "Vascular","Pericytes", "Myeloid", "T-cells",  "Mast", "Adipocytes"), 
                                                    labels = c("Basal", "LumHR", "LumSec","Fibroblasts", "Lymphatic", "Vascular","Pericytes", "Myeloid", "T-cells",  "Mast", "Adipocytes")))

# plot ggplot
# figure 1b umap -----
fill_col = fig1_colors
library(ggpubr);library(grid)
grob <- grobTree(textGrob("120024 cells", x=0.7,  y=0.95, hjust=0,
                          gp=gpar(col="black", fontsize=10, fontface="italic")))
p1 = ggplot(umap_tx, aes(x=UMAP_1, y=UMAP_2)) + #annotation_custom(grob) + 
  #geom_point(shape = 16, size = 0.005, aes(color = major_group)) + 
  #geom_point(shape = 21, size = 1, aes(fill = major_group), alpha = 0.5, stroke = 0.05, position = position_jitter(w = 0.1, h = 0.1)) + #, position = "dodge"
  geom_jitter(shape = 21, size = 2, aes(fill = major_group),  color = "black", alpha = 0.5, stroke = 0.2, position = "dodge") + 
  scale_fill_manual(values = fig1_colors) + 
  #scale_color_manual(values = fig1_colors) +
  #guides(fill = guide_legend(override.aes = list(color = fill_col)),
  #color = guide_legend(override.aes = list(shape = 21))) +
  xlab("UMAP-1") + ylab("UMAP-2") +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.border = element_blank(),
        legend.position = "none",
        #axis.title = element_blank(),
        #axis.line = element_line(arrow=arrow)
        axis.line.x = element_line(colour = "black"), #,arrow=arrow
        axis.line.y = element_line(colour = "black"),#,arrow=arrow
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank());p1

ggsave(plot = p1, filename = "figure1_nuc_umap.png", device = "png", path = save_path, width = 3.5, height = 3, dpi = 300, units = "in")
ggsave(plot = p1, filename = "figure1_nuc_umap_large.png", device = "png", path = save_path, width = 13, height = 13, dpi = 300, units = "in")
ggsave(plot = p1, filename = "figure1_nuc_umap_largedots.png", device = "png", path = save_path, width = 3.5, height = 3, dpi = 300, units = "in")
ggsave(plot = p1, filename = "figure1_nuc_umap_large_largedots3.png", device = "png", path = save_path, width = 18, height = 18, dpi = 300, units = "in")
#ggsave(plot = p1, filename = "figure1_umap.pdf", device = "pdf", path = save_path, width = 3, height = 3)
#cairo_ps(filename = "figure1_umap.eps", width = 3.5, height = 3, fallback_resolution = 300); p1; dev.off()

# figure s1c -----
p2 = umap_tx %>%  
  mutate(sampleid2 = fct_reorder(tum_meta.sample_id, as.numeric(factor(tum_meta.figure1_grp_updated)))) %>% 
  ggplot(aes(x = sampleid2, fill = major_group)) +
  geom_bar(position = "fill", width = 1) + 
  #coord_polar() +
  ylab("Cellular Proportions") + xlab("Samples") +  
  scale_fill_manual(values = fig1_colors)+
  theme(panel.background = element_blank()) +
  theme(legend.title = element_blank(), 
        legend.key.size = unit(0.5,"line"), 
        plot.background = element_blank() ,
        #axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank());p2
ggsave(plot = p2, filename = "figure_s1c_cell_barplot_groups.png", device = "png", path = save_path, width = 3.5, height = 3, dpi = 300, units = "in")
ggsave(plot = p2, filename = "figure_barplot_fill.pdf", device = "pdf", path = save_path, width = 5, height = 3)
cairo_ps(filename = "figure_barplot_fill.eps", width = 5, height = 3,fallback_resolution = 300); p2; dev.off()

p2_header = umap_tx %>%  
  mutate(sampleid2 = fct_reorder(tum_meta.sample_id, as.numeric(factor(tum_meta.figure1_grp_updated)))) %>% 
  ggplot(aes(x = sampleid2, fill = tum_meta.figure1_grp_updated)) +
  geom_bar(position = "fill", width = 1) + 
  #coord_polar() +
  ylab("Cellular Proportions") + xlab("Samples") +  
  scale_fill_manual(values = beyonce_palette(78)[2:4])+
  theme(panel.background = element_blank()) +
  theme(legend.title = element_blank(),
        legend.direction = "vertical",
        legend.key.size = unit(0.5,"line"), 
        plot.background = element_blank() ,
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank());p2_header
ggsave(plot = p2_header, filename = "figure_s1c_cell_barplot_groups_header.png", device = "png", path = save_path, width = 3.5, height = 2, dpi = 300, units = "in")

# markers ----
# FINAL OBJECT -----
# tum = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/final_use.rds")
# tum = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_1_v1/final_use.rds")
# Idents(tum) = "final_group1" #major_group

all_nucposmarkers = FindAllMarkers(object = tum, assay = "RNA", logfc.threshold = 0.25, only.pos = TRUE)
write.csv(all_nucposmarkers, paste0(save_path, "all_nucposmarkers.csv"))
allnuc_markers = FindAllMarkers(object = tum, assay = "RNA", logfc.threshold = 0, only.pos = FALSE)
write.csv(allnuc_markers, paste0(save_path, "allnuc_markers.csv"))

all_nucposmarkers = read.csv("all_nucposmarkers.csv")
all_nucposmarkers = all_nucposmarkers[,-1]

allnuc_markers = read.csv("allnuc_markers.csv")
allnuc_markers = allnuc_markers[,-1]

## run function to select markers ----
allnuc_markers %>%
  # mutate(cluster1 = factor(cluster, levels = c("adipo","endo_lymph","endo_vasc","epi_basal","epi_lumhr","epi_lumsec","fibro",
  #                                              "immune_mast", "immune_myeloid","immune_t","pericytes"), 
  #                          labels = c("Adipocytes", "Lymphatic", "Vascular", "Basal", "LumHR", "LumSec","Fibroblasts", "Mast" , "Myeloid", "T-cells","Pericytes"))) %>% 
  # 
  mutate(cluster1 = forcats::fct_relevel(cluster, levels = c("Basal", "LumHR", "LumSec","Fibroblasts", "Lymphatic", "Vascular","Pericytes", "Myeloid", "T-cells",  "Mast", "Adipocytes"))) %>% 
  arrange(cluster1) %>% 
  filter(!str_detect(gene, "^RP")) %>% 
  filter(!str_detect(gene, "^MT-")) %>% 
  filter(!str_detect(gene, "^HBB")) -> all_markers_sel

df_thresh = data.frame(
  cluster = c("Basal", "LumHR", "LumSec",
              "Fibroblasts", "Lymphatic", "Vascular","Pericytes", 
              "Myeloid", "T-cells",  "Mast", "Adipocytes"),
  thresh_pct1 = c(0.40, 0.40, 0.40, 
                  0.20, 0.20, 0.20, 0.15,
                  0.30, 0.30, 0.30, 0.20), stringsAsFactors = F
) %>% data.frame(row.names = .$cluster)

p_load(magrittr)
all_markers_sel %<>% as_tibble()

all_markers_sel2 = func_select_markers_thresh(all_markers_sel = all_markers_sel, pct_var = 0.25)
write.csv(all_markers_sel2, "all_markers_sel2_nuclei.csv")

all_markers_sel2 %>%
  group_by(cluster1) %>% 
  #top_n(50, avg_logFC) %>%
  top_n(10, avg_logFC) %>% 
  pull(gene) %>% as.character() ->
  #tum_markers_top
  tum_markers_top_rna

all_markers_sel2 %>%
  group_by(cluster1) %>% 
  #top_n(50, avg_logFC) %>%
  top_n(10, avg_logFC) %>% 
  dplyr::select(gene, cluster1) %>% mutate(id = paste0(gene, "-", cluster1)) -> genes_sel_df

tum$cellname = rownames(tum@meta.data)
data.frame(new_cluster_factor = tum@meta.data$final_group1, cell = tum$cellname) %>%
  #filter(!str_detect(cell, "untreated")) %>% 
  group_by(new_cluster_factor) %>% 
  sample_n(500) %>% 
  mutate(cell = as.character(cell)) %>% 
  pull(cell) ->
  tum_cells

# figure 1d markers ----
main_markers = c("EPCAM", "DCN", "PECAM1", "RGS5", "PTPRC", "CD68", "CD3D", "NKG7", "CD79A")

DefaultAssay(tum) = "RNA"
tum = ScaleData(tum, assay = "RNA")

p10 = DoHeatmap(tum, group.by = "final_group1", features = as.character(tum_markers_top_rna), cells = tum_cells, label=TRUE, raster = TRUE, 
                group.colors = fig1_colors, disp.min = -1.0, disp.max = 1, size = 2) + 
  theme(text = element_text(size = 8)) + 
  #scale_fill_gradient2(low="navyblue", mid="black",high="cornsilk", midpoint = -0.2)
  scale_fill_gradient2(low="steelblue4", mid="white",high="tomato4", midpoint = 0)
#scale_fill_gradient2(low="dodgerblue", mid="white",high="maroon", midpoint = 0);p10
ggsave(plot = p10, filename = paste0(save_path,"/figure_heatmap_nuc_types_redblue100L_nuc_rnapct_0.40_0.20_0.25_pct25_top15.png"), device = "png", width = 6, height = 12)
# ggsave(plot = p10, filename = paste0(save_path,"/figure_heatmap_cell_types_redblue100L_cells_integrated.pdf"), device = "pdf", width = 4, height = 6)
# ggsave(plot = p10, filename = "/figure_heatmap_cell_types_redblue100cells_rna.pdf", device = "pdf", path = save_path, width = 8, height = 6)


# dotplot ---
pal <- wes_palette("Zissou1", 5, type = "continuous")
p12 = DotPlot(tum, features = unique(as.character(tum_markers_top_rna)), dot.scale = 1.5, group.by = "major_group_factor") +
  scale_color_gradientn(colours = pal) + coord_flip() 
p13 = p12 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
                  axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size = 6)) 
ggsave(plot = p13, filename = "figure_dotplot_cell_types.png", device = "png", width = 5, height = 6.5)
ggsave(plot = p13, filename = "figure_dotplot_cell_types.pdf", device = "pdf", width = 5, height = 6.5)

# big heatmap -----
#Groups_set <- c("epi_myo", "epi_lumHR", "epi_lumSec", "fibro", "endo_lymph", "endo_vasc","pericytes","immune_myeloid","immune_t", "immune_nk", "immune_b")
Groups_set <- c("Basal", "LumHR", "LumSec",
                "Fibroblasts", "Lymphatic", "Vascular",
                "Pericytes", "Myeloid",
                "T-cells",  "Mast", "Adipocytes")
col_grp <- fig1_colors[1:11]
names(col_grp) <- Groups_set

Protocol = c("short_digestion", "6hr_digestion", "overnight_digestion")
col_Protocol = c("bisque", "darkorange1", "deeppink2")
names(col_Protocol) <- Protocol

meno = c("pre", "post", "unknown")
col_meno = c("#EA879C", "#8E2322", "gray90")
names(col_meno) <- meno

paridy = c("0", "1", "unknown")
col_paridy = c("#CFA2D2", "#5D36B1", "gray90")
names(col_paridy) <- paridy

procedure = c("Reduction Mammoplasty", "Prophylatic Mastectomy", "Cancer Mastectomy")
col_procedure = c("#F5EC54", "#84C876", "#539E59")
names(col_procedure) <- procedure

df_fibro = run_htmaps(df = all_markers, seu = tum, cell = "all_cells", genes_df = tum_markers_top_rna, 
                      clus_levels = c("Basal", "LumHR", "LumSec",
                                      "Fibroblasts", "Lymphatic", "Vascular",
                                      "Pericytes", "Myeloid",
                                      "T-cells", "B-cells"))

df_meta = data.frame(seu@meta.data, 
                     celltype = seu$major_group_factor, 
                     stringsAsFactors = F) %>% tbl_df() %>% clean_names()
df_meta_meno = df_meta %>% dplyr::select(celltype, menopause_y)
df_meta_meno1 = table(df_meta_meno$menopause_y, df_meta_meno$celltype)
df_meta_meno2 = prop.table(df_meta_meno1, margin = 2)

df_meta_paridy = df_meta %>% dplyr::select(celltype, parity_1_y)
df_meta_paridy1 = table(df_meta_paridy$parity_1_y, df_meta_meno$celltype)
df_meta_paridy2 = prop.table(df_meta_paridy1, margin = 2)

df_meta_grp = df_meta %>% dplyr::select(celltype, figure1_grp_updated)
df_meta_grp1 = table(df_meta_grp$figure1_grp_updated, df_meta_grp$celltype)
df_meta_grp2 = prop.table(df_meta_grp1, margin = 2)

df_meta_dig = df_meta %>% dplyr::select(celltype, exp_proc)
df_meta_dig1 = table(df_meta_dig$exp_proc, df_meta_grp$celltype)
df_meta_dig2 = prop.table(df_meta_dig1, margin = 2)

df_meta_grp = df_meta %>% dplyr::select(celltype, figure1_grp_updated)
df_meta_grp1 = table(df_meta_grp$figure1_grp_updated, df_meta_grp$celltype)
df_meta_grp2 = prop.table(df_meta_grp1, margin = 2)


all_nucposmarkers %>%
  mutate(cluster1 = forcats::fct_relevel(cluster, levels = c("Basal", "LumHR", "LumSec",
                                                             "Fibroblasts", "Lymphatic", "Vascular",
                                                             "Pericytes", "Myeloid",
                                                             "T-cells", "Mast", "Adipocytes"))) %>% 
  arrange(cluster1) %>% 
  filter(!str_detect(gene, "^RP")) %>% 
  filter(!str_detect(gene, "^MT-")) %>% 
  filter(!str_detect(gene, "^HBB")) -> all_markers_sel_fig


# normal mean heatmap for fig 1
df_fibro1 = run_mean_htmaps_fig1(df = all_markers_sel_fig, seu = tum, cell = "all_cells", genes_df = tum_markers_top_rna, #sel_genes
                                 genes_df_cann = main_markers, genes_sel_df = genes_sel_df, 
                                 clus_levels = c("Basal", "LumHR", "LumSec",
                                                 "Fibroblasts", "Lymphatic", "Vascular",
                                                 "Pericytes", "Myeloid",
                                                 "T-cells",  "Mast", "Adipocytes"))
# horizontal mean heatmap for fig 1
df_fibro1 = run_mean_htmaps_fig1_horizontal(df = all_markers_sel_fig, seu = tum, cell = "all_cells", genes_df = tum_markers_top_rna, #sel_genes
                                            genes_df_cann = main_markers, genes_sel_df = genes_sel_df, 
                                            clus_levels = c("Basal", "LumHR", "LumSec",
                                                            "Fibroblasts", "Lymphatic", "Vascular",
                                                            "Pericytes", "Myeloid",
                                                            "T-cells",  "Mast", "Adipocytes"))
# supp mean heatmap for fig 1
df_fibro1 = run_mean_htmaps_fig1_supp(df = all_markers_sel_fig, seu = tum, cell = "all_cells", genes_df = tum_markers_top_rna, #sel_genes
                                            genes_df_cann = main_markers, genes_sel_df = genes_sel_df, 
                                            clus_levels = c("Basal", "LumHR", "LumSec",
                                                            "Fibroblasts", "Lymphatic", "Vascular",
                                                            "Pericytes", "Myeloid",
                                                            "T-cells",  "Mast", "Adipocytes"))
# supp horizontal mean heatmap for fig 1
df_fibro1 = run_mean_htmaps_fig1_horizontal_supp(df = all_markers_sel_fig, seu = tum, cell = "all_cells", genes_df = tum_markers_top_rna, #sel_genes
                                            genes_df_cann = main_markers, genes_sel_df = genes_sel_df, 
                                            clus_levels = c("Basal", "LumHR", "LumSec",
                                                            "Fibroblasts", "Lymphatic", "Vascular",
                                                            "Pericytes", "Myeloid",
                                                            "T-cells",  "Mast", "Adipocytes"))

# pathway ------
genes = read.csv("/volumes/lab/users/tapsi/projects/bcap/analysis/20200210_final1/figure_1_v2/allcell_markers.csv")
genes = genes[,-1]
genes$cluster = as.character(genes$cluster)

gsea_resultsc2hall = fun_gsea(genes = genes, set = "H", sub_cat = NULL, num_avg_logFC = 0.25)
gsea_resultsc2C5 = fun_gsea(genes = genes, set = "C5", sub_cat = "BP", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "fibro")
gsea_resultsc2C3 = fun_gsea(genes = genes, set = "C3", sub_cat = "TFT:GTRD", num_avg_logFC = 0.25, celltype = "lumsec", doplot = FALSE, ret_collap = FALSE)
gsea_resultsc2C8 = fun_gsea(genes = genes, set = "C8", sub_cat = NULL, num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "fibro")
gsea_resultsc2C2_reac = fun_gsea(genes = genes, set = "C2", sub_cat = "CP:REACTOME", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "fibro")
gsea_resultsc2C2_kegg = fun_gsea(genes = genes, set = "C2", sub_cat = "CP:KEGG", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "fibro")

clus = c("fibro-1", "fibro-2", "fibro-3")

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
colnames(p1.1) = c("pathway", "fibro-1", "fibro-2", "fibro-3")

# rearrange the datafram
test = p1

# for basal
test1 = test %>% dplyr::arrange(-`fibro-1`) %>% dplyr::filter(`fibro-1`>0) %>% pull(pathway)
test2 = test %>% dplyr::filter(!pathway %in% test1) %>% arrange(-`fibro-2`) %>% dplyr::filter(`fibro-1`<0, `fibro-2`>0, `fibro-3`<0) %>% dplyr::top_n(25, `fibro-2`) %>% pull(pathway)
test3 = test %>% dplyr::filter(!pathway %in% c(test1, test2)) %>% arrange(-`fibro-3`) %>% dplyr::filter(`fibro-1`<0, `fibro-2`<0, `fibro-3`>0) %>% dplyr::top_n(25, `fibro-3`) %>% pull(pathway)
seq_path = c(test1, test2, test3) #, test4, test5


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

png("/htmaps_pathways_kegg_reac_biocarta.png", res = 300, units = "in", width = 12, height = 12)
print(ha)
dev.off()

png("fibro_htmaps_pathways_GO_unique_ones_NES.png", res = 300, units = "in", width = 10, height = 24)
print(ht_hall)
dev.off()

png("htmaps_pathways_hallmark_NES.png", res = 300, units = "in", width = 6, height = 6)
print(ht_hall)
dev.off()

png("lum_hr_htmaps_pathways_kegg_NES.png", res = 300, units = "in", width = 6, height = 4)
print(ht_hall)
dev.off()

png("lum_hr_htmaps_pathways_reactome_NES_unique_top25.png", res = 300, units = "in", width = 10, height = 10)
print(ht_hall)
dev.off()

png("htmaps_pathways_immune_NES.png", res = 300, units = "in", width = 10, height = 20)
print(ht_hall)
dev.off()

png("lum_sec_htmaps_pathways_c3_NES.png", res = 300, units = "in", width = 10, height = 6)
print(ht_hall)
dev.off()

png("htmaps_pathways_atlas_NES.png", res = 300, units = "in", width = 10, height = 20)
print(ht_hall)
dev.off()

# Figure supp 7 violin plot --------

v1 = VlnPlot(tum, features = "LIF", assay = "RNA", group.by = "final_group1", pt.size = 0, cols = fig1_colors)
ggsave(plot = v1, filename = "periviolin_LIF_celltypes.pdf", device = "pdf", path = save_path, width = 5, height = 3, dpi = 300, units = "in")
DefaultAssay(tum) <- "RNA"
v1 = FeaturePlot(tum, features = "LIF", order = T)
ggsave(plot = v1, filename = "perifeature_LIF_celltypes.png", device = "png", path = save_path, width = 5, height = 3, dpi = 300, units = "in")
v1 = FeaturePlot(tum, features = "LIF")
ggsave(plot = v1, filename = "perifeature_LIF_celltypesNoorder.png", device = "png", path = save_path, width = 5, height = 3, dpi = 300, units = "in")


# Figure 1 stacked barplots --------

# cell
tum = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_1_v1/final_use.rds")
Idents(tum) = "final_group1"
df_meta = data.frame(tum@meta.data, 
                      celltype = tum$final_group1, 
                      stringsAsFactors = F) %>% 
  as_tibble() %>% clean_names()
write_rds(df_meta, "df_meta_cell.rds")
df_meta = readRDS("df_meta_cell.rds")

# nuc
nuc = read_rds("final_nuclei_use.rds")

df_meta_nuc = data.frame(nuc@meta.data, 
                     celltype = nuc$final_group1, 
                     stringsAsFactors = F) %>% 
  as_tibble() %>% clean_names()
write_rds(df_meta_nuc, "df_meta_nuc.rds")
df_meta_nuc = readRDS("df_meta_nuc.rds")

# merge 
df_meta = df_meta %>% dplyr::select(sheet_patient_id, celltype, final_sample_id) %>% mutate(type = "cell")
colnames(df_meta) = c("patient_id", "celltype", "sample_id", "type")
df_meta_nuc = df_meta_nuc %>% dplyr::select(patient_id, celltype, sample_id) %>% mutate(type = "nuclei")
df_meta = rbind(df_meta, df_meta_nuc)
df_meta  =  df_meta %>% mutate(major_group = factor(celltype, levels = c("Basal", "LumHR", "LumSec","Fibroblasts", "Lymphatic", "Vascular","Pericytes", "Myeloid", "T-cells", "B-cells", "Mast", "Adipocytes"), 
                                                    labels = c("Non-immune", "Non-immune", "Non-immune","Non-immune", "Non-immune", "Non-immune","Non-immune", "Immune", "Immune", "Immune", "Immune", "Non-immune")))
levels(df_meta$major_group) = c("Immune", "Non-immune")

# plots
func_bar(var_df = df_meta, var_x = "type", var_fill = "celltype", var_pos = "fill", var_cols = fig1_colors, celltype = "cells")
func_bar(var_df = df_meta, var_x = "type", var_fill = "major_group", var_pos = "fill", var_cols = c("green", "purple"), celltype = "cells")

# plot for reduction ipsi etc
df_meta = data.frame(tum@meta.data, 
                     celltype = tum$final_group1, 
                     stringsAsFactors = F) %>% 
  as_tibble() %>% clean_names()
df_meta = df_meta %>% dplyr::select(sheet_patient_id, celltype, final_sample_id, sheet_figure1_grp_updated) %>% mutate(type = "cell")
df_meta  =  df_meta %>% mutate(major_group = factor(celltype, levels = c("Basal", "LumHR", "LumSec","Fibroblasts", "Lymphatic", "Vascular","Pericytes", "Myeloid", "T-cells", "B-cells"), 
                                                    labels = c("Non-immune", "Non-immune", "Non-immune","Non-immune", "Non-immune", "Non-immune","Non-immune", "Immune", "Immune", "Immune")))
levels(df_meta$major_group) = c("Immune", "Non-immune")
func_bar(var_df = df_meta, var_x = "sheet_figure1_grp_updated", var_fill = "major_group", var_pos = "fill", var_cols = c("green", "purple"), celltype = "fig1grp")

