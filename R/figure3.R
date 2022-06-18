# geo ----------
setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/")
odir = "figure_3_v1";dir.create(odir)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")

p_load(future)
plan("multiprocess", workers = 40)
options(future.globals.maxSize = 90000 * 1024^2)

data_folder = "03_subset_lumhr_v2"
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/figure_3_v1/"
setwd(save_path)
col_m = c("#f9f871", "#00c9ae", "#ff8863", "#6e3452", "#00aafd", "#b178de", "#ff86da", "#7d9300","red2","#C70E7B", "#FC6882", "#A6E000", "#1BB6AF", "#6C6C9D", "#172869", "#D64358", "#EAFB88", "#3C8C4D", "#DFCEE0","orange", "red", "darkorchid1","darkorchid4", colors_dark)


# core --------
setwd("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/")
odir = "figure_3_v1";dir.create(odir)
setwd(paste0("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")

p_load(future)
plan("multiprocess", workers = 100)
options(future.globals.maxSize = Inf)

data_folder = "figure_3_v1"
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_3_v1/"
setwd(save_path)
col_m = c("#f9f871", "#00c9ae", "#ff8863", "#6e3452", "#00aafd", "#b178de", "#ff86da", "#7d9300","red2","#C70E7B", "#FC6882", "#A6E000", "#1BB6AF", "#6C6C9D", "#172869", "#D64358", "#EAFB88", "#3C8C4D", "#DFCEE0","orange", "red", "darkorchid1","darkorchid4", colors_dark)

# ************************* -----
# 1. MAIN figure -----
# ** epithelial cell figure -----

# **read data ----
lumhr = readRDS(file = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumhr_v2/clus_out.rds")
lumhr = readRDS(file = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumhr_v2/clus_out.rds")

lumsec = readRDS(file = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_v3/lumsec_clus_out.rds")
lumsec = readRDS(file = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_lumsec_v3/lumsec_clus_out.rds")

basal = readRDS(file = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_basal_v2/basal_clus_out.rds")
basal = readRDS(file = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_basal_v2/basal_clus_out.rds")
basal$sheet_patient_id[basal$sheet_ts_sampleid %in% c("hbca25_c")] = "pt37"

# **annotate clusters ----
lumhr$epi_state = lumhr$cellstate_tk2
lumsec$epi_state = lumsec$cellstate_tk1
basal$epi_state = "basal"

DimPlot(lumhr, group.by = "epi_state", cols = colors_bottle2)
DimPlot(lumsec, group.by = "epi_state", cols = col_m)
DimPlot(basal, group.by = "epi_state", cols = "#C00000FF")

# merge the objects
epi_merged = merge(lumhr, lumsec)
epi_merged = merge(epi_merged, basal)

epi_merged$celltype = "NULL"
epi_merged$celltype[epi_merged$epi_state %in% c("lhr_EGLN3", "lhr_FASN", "lhr_PIP")] = "lumhr"
epi_merged$celltype[epi_merged$epi_state %in% c('ls_cd74',"ls_HMOX1","ls_KIT","ls_KRT23", "ls_lactation","ls_MT1X","ls_prol")] = "lumsec"
epi_merged$celltype[epi_merged$epi_state %in% c("basal")] = "basal"


DefaultAssay(epi_merged) <- "RNA"
write_rds(epi_merged, "epi_merged.rds")
epi_merged1 = FindVariableFeatures(epi_merged)
epi_merged1 = ScaleData(epi_merged1)
epi_merged1 = RunPCA(epi_merged1, dims = 1:40)
epi_merged1 = RunUMAP(epi_merged1, dims = 1:20)
DimPlot(epi_merged1, group.by = "epi_state", cols = col_m)



# **re-integrate -------
seu_object = read_rds("epi_merged.rds")
dims = 20; k.param = 30; cols = colors_dark;var_split = "sheet_exp_proc"
fts1 = c("nFeature_RNA","nCount_RNA", "percent.mito","percent.rb", "percent.sg1", "percent.sg2", "percent.sg3")
fts2 = c("rna_EPCAM","rna_PTPRC","rna_PECAM1","rna_DCN","rna_KRT5",  "rna_LTF", "rna_KRT17", "rna_HPGD")
fts3 = c("rna_PDE3B","rna_ACACB","rna_ACSM1","rna_ACTA2","rna_CD4",  "rna_NKG7", "rna_CD19", "rna_MS4A1", "rna_RGS5")
fts4 = c("rna_VWF","rna_PRSS1","rna_CCL21","rna_FAP","COL1A1",  "rna_COL6A1", "rna_ESR2", "rna_KRT8", "rna_SERHL2")
clus_out = func_cluster_rpca(tum_int = seu_object, var_split = var_split, var_scale = "NULL", dims = dims, k.param = k.param, cols = colors_dark,  
                             markers1 = fts1, markers2 = fts2, markers3 = fts3, markers4 = fts4, save_path = save_path)


clus_out = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/figure_3_v1/dim_20k_30_NULL_tum.integrated.ns2_method_9rcp_0.1_1.75.rds")

DimPlot(clus_out, group.by = "epi_state", cols = colors_dark)

clus_out$final_group = "NULL"
clus_out$final_group[clus_out$epi_state %in% c("basal")] = "basal"
clus_out$final_group[clus_out$epi_state %in% c("lhr_FASN")] = "lumhr_fasn"
clus_out$final_group[clus_out$epi_state %in% c("lhr_PIP")] = "lumhr_pip"
clus_out$final_group[clus_out$epi_state %in% c("ls_cd74")] = "lumsec_hla"
clus_out$final_group[clus_out$epi_state %in% c("ls_HMOX1", "ls_KRT23", "ls_MT1X")] = "lumsec_krt23"
clus_out$final_group[clus_out$epi_state %in% c("ls_lactation")] = "lumsec_lactation"
clus_out$final_group[clus_out$epi_state %in% c("ls_prol")] = "lumsec_prol"
clus_out$final_group[clus_out$epi_state %in% c("ls_KIT")] = "lumsec_kit"

clus_out$final_group = factor(clus_out$final_group, levels = c("basal", "lumhr_egln3", "lumhr_fasn", "lumhr_pip", "lumsec_krt23", "lumsec_kit", "lumsec_lactation", 
                                                               "lumsec_hla", "lumsec_prol"))

clus_out$final_group1[clus_out$final_group %in% c("lumhr_egln3", "lumhr_fasn", "lumhr_pip")] = "LumHR"
clus_out$final_group1[clus_out$final_group %in% c("lumsec_krt23", "lumsec_kit", "lumsec_lactation", 
                                                  "lumsec_hla", "lumsec_prol")] = "LumSec"
clus_out$final_group1[clus_out$final_group %in% c("basal")] = "Basal"
clus_out$final_group1 = factor(clus_out$final_group1, levels = c("Basal", "LumHR", "LumSec"))

# ^^ Final object -------
write_rds(clus_out, "epi_final_use.rds")

#clus_out = read_rds("epi_final_use.rds")

fibro = clus_out

DimPlot(clus_out, group.by = "final_group", cols = c("#C00000FF", "#FF0000FF", "#FFC000FF"), label = T)
DimPlot(clus_out, group.by = "final_group1", cols = c("#C00000FF", "#FF0000FF", "#FFC000FF"), label = T)
Idents(b) <- "final_group"
mks = FindAllMarkers(b, assay = "RNA")


top_markers = mks %>% arrange(cluster) %>% 
  dplyr::group_by(cluster) %>%
  top_n(7, avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()

DoHeatmap(b, top_markers, group.colors = colors_dark)

# **Umap ------
tum = clus_out
umap_tx = tum@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(tum_meta = tum@meta.data)

umap_tx  =  umap_tx %>% mutate(major_group = factor(tum_meta.final_group, levels = c("basal", "lumhr_egln3", "lumhr_fasn", "lumhr_pip", "lumsec_krt23", "lumsec_kit","lumsec_lactation", 
                                                                                     "lumsec_hla", "lumsec_prol")))
umap_tx  =  umap_tx %>% mutate(major_group1 = factor(tum_meta.final_group1, levels = c("Basal", "LumHR", "LumSec")))

fill_col = c("#C00000FF", "#FF0000FF", "#FFC000FF")
fill_col = c("#C00000FF", "#FFBCD9", "#FF5159", "magenta1", "#FFD125", "#A26765",  "#FF8000", "#CA9300", "#B8BA00")
library(ggpubr);library(grid)
grob <- grobTree(textGrob("6977 cells", x=0.7,  y=0.95, hjust=0,
                          gp=gpar(col="black", fontsize=10, fontface="italic")))
umap_tx1 <- umap_tx[sample(x = 1:nrow(x = umap_tx)), ]
p1 = ggplot(umap_tx1, aes(x=UMAP_1, y=UMAP_2)) + #annotation_custom(grob) + 
  geom_jitter(shape = 21, size = 1, aes(fill = major_group),  color = "black", alpha = 0.9, stroke = 0.04) + 
  scale_fill_manual(values = fill_col) +
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

ggsave(plot = p_fixed, filename = "epi_major_umap.png", device = "png", path = save_path, width = 3, height = 3, dpi = 300, units = "in")
ggsave(plot = p_fixed, filename = "epi_major_umap.pdf", device = "pdf", path = save_path, width = 3, height = 3, dpi = 300, units = "in")

p03 = DimPlot(object = tum, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = fill_col, group.by = "final_group1", shuffle = T)+  NoAxes()
pdf(glue("{save_path}/epi_major_uamp_1.pdf"), width=6, height=5)
print(plot_grid(p03))
dev.off()


# **Umap basal ------

basal$final_group[basal$epi_state %in% c("basal")] = "basal"

umap_tx_basal = basal@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(tum_meta = basal@meta.data)

umap_tx_basal  =  umap_tx_basal %>% mutate(major_group = factor(tum_meta.final_group, levels = c("basal")))

fill_col = c("#C00000FF")
library(ggpubr);library(grid)
# grob <- grobTree(textGrob("6977 cells", x=0.7,  y=0.95, hjust=0,
#                           gp=gpar(col="black", fontsize=10, fontface="italic")))
umap_tx_basal1 <- umap_tx_basal[sample(x = 1:nrow(x = umap_tx_basal)), ]
p1 = ggplot(umap_tx_basal, aes(x=UMAP_1, y=UMAP_2)) + #annotation_custom(grob) + 
  geom_jitter(shape = 21, size = 1, aes(fill = major_group),  color = "black", alpha = 0.9, stroke = 0.04) + 
  scale_fill_manual(values = fill_col) +
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

ggsave(plot = p_fixed, filename = "basal_epi_major_umap.png", device = "png", path = save_path, width = 3, height = 3, dpi = 300, units = "in")
ggsave(plot = p_fixed, filename = "basal_epi_major_umap.pdf", device = "pdf", path = save_path, width = 3, height = 3, dpi = 300, units = "in")

p03 = DimPlot(object = basal, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = fill_col, group.by = "final_group", shuffle = T)+  NoAxes()
pdf(glue("{save_path}/basal_epi_major_uamp_1.pdf"), width=6, height=5)
print(plot_grid(p03))
dev.off()

write_rds(basal, "basal_final.rds")
basal = read_rds("basal_final.rds")

# **Umap lumHR ------

lumhr$final_group[lumhr$epi_state %in% c("lhr_EGLN3")] = "lumhr_egln3"
lumhr$final_group[lumhr$epi_state %in% c("lhr_FASN")] = "lumhr_fasn"
lumhr$final_group[lumhr$epi_state %in% c("lhr_PIP")] = "lumhr_pip"
lumhr$final_group = factor(lumhr$final_group, levels = c("lumhr_egln3", "lumhr_fasn", "lumhr_pip"))
Idents(lumhr) <- "final_group"

umap_tx_lumhr = lumhr@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(tum_meta = lumhr@meta.data)

umap_tx_lumhr  =  umap_tx_lumhr %>% mutate(major_group = factor(tum_meta.final_group, levels = c("lumhr_egln3", "lumhr_fasn", "lumhr_pip")))

fill_col = c("#FFBCD9", "#FF5159", "magenta1")
library(ggpubr);library(grid)
# grob <- grobTree(textGrob("6977 cells", x=0.7,  y=0.95, hjust=0,
#                           gp=gpar(col="black", fontsize=10, fontface="italic")))
#umap_tx_lumhr1 <- umap_tx_lumhr[sample(x = 1:nrow(x = umap_tx_lumhr)), ]
p1 = ggplot(umap_tx_lumhr, aes(x=UMAP_1, y=UMAP_2)) + #annotation_custom(grob) + 
  geom_jitter(shape = 21, size = 1, aes(fill = major_group),  color = "black", alpha = 0.9, stroke = 0.04) + 
  scale_fill_manual(values = fill_col) +
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

ggsave(plot = p_fixed, filename = "lumhr_epi_major_umap.png", device = "png", path = save_path, width = 3, height = 3, dpi = 300, units = "in")
ggsave(plot = p_fixed, filename = "lumhr_epi_major_umap.pdf", device = "pdf", path = save_path, width = 3, height = 3, dpi = 300, units = "in")

p03 = DimPlot(object = lumhr, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = fill_col, group.by = "final_group", shuffle = T)+  NoAxes()
pdf(glue("{save_path}/lumhr_epi_major_uamp_1.pdf"), width=6, height=5)
print(plot_grid(p03))
dev.off()

write_rds(lumhr, "lumhr_final.rds")
lumhr = read_rds("lumhr_final.rds")

all_lumhr = FindAllMarkers(object = lumhr, assay = "RNA", logfc.threshold = 0, only.pos = FALSE)
write.csv(all_lumhr, paste0(save_path, "lumhr_markers.csv"))


# **Umap LumSec ------
lumsec$final_group[lumsec$epi_state %in% c("ls_cd74")] = "lumsec_hla"
lumsec$final_group[lumsec$epi_state %in% c("ls_HMOX1", "ls_KRT23", "ls_MT1X")] = "lumsec_krt23"
lumsec$final_group[lumsec$epi_state %in% c("ls_lactation")] = "lumsec_lactation"
lumsec$final_group[lumsec$epi_state %in% c("ls_prol")] = "lumsec_prol"
lumsec$final_group[lumsec$epi_state %in% c("ls_KIT")] = "lumsec_kit"
lumsec$final_group  =  factor(lumsec$final_group, levels = c("lumsec_krt23", "lumsec_kit", "lumsec_lactation", 
                                                                                                   "lumsec_hla", "lumsec_prol"))

Idents(lumsec) <- "final_group"

umap_tx_lumsec= lumsec@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(tum_meta = lumsec@meta.data)

umap_tx_lumsec  =  umap_tx_lumsec %>% mutate(major_group = factor(tum_meta.final_group, levels = c("lumsec_KRT23", "lumsec_kit", "lumsec_lactation", 
                                                                                                   "lumsec_hla", "lumsec_prol")))
fill_col = c("#FFD125", "#A26765",  "#FF8000", "#CA9300", "#B8BA00")
library(ggpubr);library(grid)
grob <- grobTree(textGrob("6977 cells", x=0.7,  y=0.95, hjust=0,
                          gp=gpar(col="black", fontsize=10, fontface="italic")))
umap_tx_lumsec1 <- umap_tx_lumsec[sample(x = 1:nrow(x = umap_tx)), ]
p1 = ggplot(umap_tx_lumsec, aes(x=UMAP_1, y=UMAP_2)) + #annotation_custom(grob) + 
  geom_jitter(shape = 21, size = 1, aes(fill = major_group),  color = "black", alpha = 0.9, stroke = 0.04) + 
  scale_fill_manual(values = fill_col) +
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

ggsave(plot = p_fixed, filename = "lumsec_epi_major_umap.png", device = "png", path = save_path, width = 3, height = 3, dpi = 300, units = "in")
ggsave(plot = p_fixed, filename = "lumsec_epi_major_umap.pdf", device = "pdf", path = save_path, width = 3, height = 3, dpi = 300, units = "in")

p03 = DimPlot(object = lumsec, reduction = "umap", label = TRUE, pt.size = 0.000001, cols = fill_col, group.by = "final_group", shuffle = T)+  NoAxes()
pdf(glue("{save_path}/lumsec_epi_major_uamp_1.pdf"), width=6, height=5)
print(plot_grid(p03))
dev.off()

write_rds(lumsec, "lumsec_final.rds")
lumsec = read_rds("lumsec_final.rds")

all_lumsec = FindAllMarkers(object = lumsec, assay = "RNA", logfc.threshold = 0, only.pos = FALSE)
write.csv(all_lumsec, paste0(save_path, "lumsec_markers.csv"))


# **markers ----
# FINAL OBJECT 
# tum = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/figure_6_v1/mye_final_use.rds")
# tum = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_6_v1/mye_final_use.rds")

table(clus_out$final_group)

Idents(clus_out) = "final_group"

all_markers = FindAllMarkers(object = clus_out, assay = "RNA", logfc.threshold = 0.25, only.pos = TRUE)
write.csv(all_markers, paste0(save_path, "epi_fig_markers_final_group.csv"))
# all_markers = read.csv("./epi_fig_markers_final_group.csv")
# all_markers = all_markers[,-1]

allcell_markers = FindAllMarkers(object = clus_out, assay = "RNA", logfc.threshold = 0, only.pos = FALSE)
write.csv(allcell_markers, paste0(save_path, "epi_allcell_markers.csv"))
# allcell_markers = read.csv("epi_allcell_markers.csv")
# allcell_markers = allcell_markers[,-1]

# **Heatmap KRTs ------
epi_merged1 = clus_out
vec_krts = epi_merged1@assays$RNA@var.features[str_detect(epi_merged1@assays$RNA@var.features, '^KRT')]
vec_krts = vec_krts[!str_detect(vec_krts, '^KRTA|^KRTD')]
vec_krts1 = c("KRT34", "KRT86", "KRT75", "KRT17", "KRT5", "KRT33B", "KRT14", "KRT6A", "KRT1", "KRT6B", "KRT10", "KRT8", "KRT18", "KRT20", "KRT19",
             "KRT80")
vec_krts2 = vec_krts[!vec_krts %in% vec_krts1]
final_vec_krt = c(vec_krts1,vec_krts2) %>% unique()  

epi_merged1$celltype = "NULL"
epi_merged1$celltype[epi_merged1$epi_state %in% c("lhr_EGLN3", "lhr_FASN", "lhr_PIP")] = "lumhr"
epi_merged1$celltype[epi_merged1$epi_state %in% c('ls_cd74',"ls_HMOX1","ls_KIT","ls_KRT23", "ls_lactation","ls_MT1X","ls_prol")] = "lumsec"
epi_merged1$celltype[epi_merged1$epi_state %in% c("basal")] = "basal"
epi_merged1$celltype = factor(epi_merged1$celltype, levels = c("basal", "lumhr", "lumsec"))

Idents(epi_merged1) <- "celltype"
krt_obj = AverageExpression(epi_merged1, features = final_vec_krt, return.seurat = T, assays = "RNA")

test = krt_obj@assays$RNA@scale.data
test = cbind(rownames(test), test)
test1 = test[order(test[,2], decreasing = T),]
test2 = test[order(test[,3], decreasing = T),]


p1 = DoHeatmap(krt_obj, features = final_vec_krt, assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = c("#C00000FF", "#FF0000FF", "#FFC000FF"), raster = T)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
ggsave(plot = p1, filename = "heatmap_krts_celltype_v2.pdf", width = 4, height = 6)

# by cell state
Idents(epi_merged1) <- "final_group"
krt_obj = AverageExpression(epi_merged1, features = final_vec_krt, return.seurat = T, assays = "RNA")
# 
# test = krt_obj@assays$RNA@scale.data
# test = cbind(rownames(test), test)
# test1 = test[order(test[,2], decreasing = T),]
# test2 = test[order(test[,3], decreasing = T),]


p1 = DoHeatmap(krt_obj, features = final_vec_krt, assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = c("#C00000FF", "#FFBCD9", "#FF5159", "#A26765", "#CA9300", "#583000", "#FFD125", "#A58D65", "#B8BA00"), raster = T)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
ggsave(plot = p1, filename = "heatmap_krts_cellstate_v2.pdf", width = 4, height = 6)

# ** HLA's heatmap ----------
final_hlas = c("CD74", "HLA-DRA", "HLA-DPA1", "HLA-DRB1", "HLA-DPB1", "HLA-DQB1", 
"HLA-DQA1", "HLA-DMA", "HLA-DRB5", "HLA-DMB", "HLA-DQA2", "HLA-DOA", "HLA-A", "HLA-B", "HLA-C", 
"HLA-E", "HLA-F", "HLA-G")
hla_obj = AverageExpression(epi_merged1, features = c(final_hlas), return.seurat = T, assays = "RNA")
p1 = DoHeatmap(hla_obj, features = final_hlas, assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = c("#C00000FF", "#FF0000FF", "#FFC000FF"), raster = F)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
ggsave(plot = p1, filename = "heatmap_hlass_celltype_v2.pdf", width = 4, height = 4)

final_hlas = c("CD74", "HLA-DRA", "HLA-DPA1", "HLA-DRB1", "HLA-DPB1", "HLA-DQB1", 
               "HLA-DQA1", "HLA-DMA", "HLA-DRB5", "HLA-DMB", "HLA-DQA2", "HLA-DOA", "HLA-A", "HLA-B", "HLA-C", 
               "HLA-E", "HLA-F", "HLA-G")
Idents(epi_merged1) <- "final_group"
hla_obj2 = AverageExpression(epi_merged1, features = c(final_hlas), return.seurat = T, assays = "RNA")
p2 = DoHeatmap(hla_obj2, features = final_hlas, assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = c("#C00000FF", "#FFBCD9", "#FF5159", "magenta1", "#FFD125", "#A26765",  "#FF8000", "#CA9300", "#B8BA00")
, raster = F)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
ggsave(plot = p2, filename = "heatmap_hlass_celltype_v3.pdf", width = 4, height = 4)

p3 = VlnPlot(clus_out, features = final_hlas, pt.size = 0, cols = c("#C00000FF", "#FFBCD9", "#FF5159", "magenta1", "#FFD125", "#A26765",  "#FF8000", "#CA9300", "#B8BA00"))
ggsave(plot = p3, filename = "violin_hlass_cellstate_v3.pdf", width = 14, height = 14)

# ** Secretoglobins heatmap ----------
DefaultAssay(clus_out) <- "RNA"
epi_merged1 = ScaleData(clus_out)
vec_scgb = rownames(epi_merged1@assays$RNA@scale.data)[str_detect(rownames(epi_merged1@assays$RNA@scale.data), '^SCGB')]


hla_obj = AverageExpression(epi_merged1, features = c(vec_scgb), return.seurat = T, assays = "RNA")
p1 = DoHeatmap(hla_obj, features = vec_scgb, assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = c("#C00000FF", "#FFBCD9", "#FF5159", "magenta1", "#FFD125", "#A26765",  "#FF8000", "#CA9300", "#B8BA00"), raster = F)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
ggsave(plot = p1, filename = "heatmap_scgb_celltype_v2.pdf", width = 4, height = 4)

final_sec = c("SCGB1C1", "SCGB1C2","SCGB3A2" , "SCGB1A1","SCGB1B2P", "SCGB2B2", "SCGB1D2",  "SCGB2A2", "SCGB3A1","SCGB2A1" , "SCGB1D1",
              "SCGB1D4")
p2 = DoHeatmap(hla_obj, features = final_sec, assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = c("#C00000FF", "#FFBCD9", "#FF5159", "magenta1", "#FFD125", "#A26765",  "#FF8000", "#CA9300", "#B8BA00")
               , raster = F)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
ggsave(plot = p2, filename = "heatmap_scgb_celltype_v2.pdf", width = 4, height = 3)

p3 = VlnPlot(clus_out, features = final_sec, pt.size = 0, cols = c("#C00000FF", "#FFBCD9", "#FF5159", "magenta1", "#FFD125", "#A26765",  "#FF8000", "#CA9300", "#B8BA00"))
ggsave(plot = p3, filename = "violin_SCGBs_cellstate_v3.pdf", width = 8, height = 8)

# test avg ------
# testing avg calc
m1 = matrix(data = c(1,2,3,4,5,6,7,8,1,2,3,4), nrow = 3)
colnames(m1) = c("col1", "col1", "col2", "col3")
m1avg = vapply(unique(colnames(m1)), function(x) 
  rowMeans(m1[,colnames(m1)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(m1)) )

final_sec = c("SCGB1C1", "SCGB1C2","SCGB3A2" , "SCGB1A1","SCGB1B2P", "SCGB2B2", "SCGB1D2",  "SCGB2A2", "SCGB3A1","SCGB2A1" , "SCGB1D1",
              "SCGB1D4")
hla_obj = AverageExpression(epi_merged1, features = c(vec_scgb), return.seurat = T, assays = "RNA")

mat_exp_rna1 = epi_merged1@assays$RNA@scale.data[final_sec,] # avg the scaled data
colnames(mat_exp_rna1) = epi_merged1$final_group
hla_obj1 <- vapply(unique(colnames(mat_exp_rna1)), function(x) 
  rowMeans(mat_exp_rna1[,colnames(mat_exp_rna1)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(mat_exp_rna1)) )
hla_obj1_seu = hla_obj
hla_obj1_seu@assays$RNA@scale.data = hla_obj1

mat_exp_rna2 = epi_merged1@assays$RNA@counts[final_sec,] # avg the raw counts and then scaling
colnames(mat_exp_rna2) = epi_merged1$final_group
hla_obj2 <- vapply(unique(colnames(mat_exp_rna2)), function(x) 
  rowMeans(mat_exp_rna2[,colnames(mat_exp_rna2)== x,drop=FALSE], na.rm=TRUE),
  numeric(nrow(mat_exp_rna2)) )
hla_obj2_seu = hla_obj
hla_obj2_seu@assays$RNA@counts = hla_obj2
hla_obj2_seu = ScaleData(hla_obj2_seu, )
p1 = DoHeatmap(hla_obj, features = final_sec, assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = c("#C00000FF", "#FFBCD9", "#FF5159", "magenta1", "#FFD125", "#A26765",  "#FF8000", "#CA9300", "#B8BA00")
               , raster = F)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
p2 = DoHeatmap(hla_obj1_seu, features = final_sec, assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1, disp.max = 0.5, draw.lines = F, group.colors = c("#C00000FF", "#FFBCD9", "#FF5159", "magenta1", "#FFD125", "#A26765",  "#FF8000", "#CA9300", "#B8BA00")
               , raster = F)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
p1 + p2

# nuclei hormone signalling -----
load("epi_nucs_1.RObj")
Idents(epithelial.nucs) <- "epi_nuc_idents"
epithelial.nucs$epi_new = epithelial.nucs$epi_nuc_idents
epithelial.nucs$epi_new[epithelial.nucs$epi_nuc_idents == "Luminal_HR_Prolif"] = "Luminal_HR"
epithelial.nucs$epi_new[epithelial.nucs$epi_nuc_idents == "Luminal_Sec_Prolif"] = "Luminal_Sec"

epithelial.nucs$epi_new = factor(epithelial.nucs$epi_new, levels = c("Basal", "Luminal_HR", "Luminal_Sec"))
Idents(epithelial.nucs) <- "epi_new"
Idents(epithelial.nucs) <- "epi_nuc_idents"

hrs = c("ESR1", "AR", "PGR")
epithelial.nucs = ScaleData(epithelial.nucs)
HR_obj = AverageExpression(epithelial.nucs, features = hrs, return.seurat = T, assays = "RNA")
p1 = DoHeatmap(HR_obj, features = hrs, assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = c("#C00000FF", "#FF0000FF", "#FFC000FF"), raster = F)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
ggsave(plot = p1, filename = "heatmap_HR_nuctype_v2.pdf", width = 4, height = 4)

library(patchwork)
#fix.sc <- scale_color_gradientn(colours = c('beige', 'hotpink')) #,  limits = c(0, 4)
fix.sc <- scale_color_viridis_c(direction = -1)
p11 = FeaturePlot(epithelial.nucs, features = hrs, raster = F, combine = FALSE )
p11 + fix.sc
p2 <- lapply(p11, function(x) x + fix.sc)
p3 = DimPlot(epithelial.nucs, group.by = "epi_nuc_idents", cols = c("#C00000FF", "#FF0000FF", "purple", "#FFC000FF", "#B8BA00"))
p2c1 = p2[[1]] + p2[[2]] + p2[[3]]
ggsave(plot = p2c1, filename = "umap_HR_nuctype_v2.pdf", width = 8.5, height = 2.5)
ggsave(plot = , filename = "umap_HR_nuctype_v2.png", width = 8.5, height = 2.5)
ggsave(plot = p2[[1]], filename = "umap_ESR1_nuctype_v2.png", width = 3, height = 2.5)
ggsave(plot = p2[[2]], filename = "umap_AR_nuctype_v2.png", width = 3, height = 2.5)
ggsave(plot = p2[[3]], filename = "umap_PGR_nuctype_v2.png", width = 3, height = 2.5)
ggsave(plot = p3, filename = "umap_dimplot_nuctype_v2.png", width = 4.5, height = 2.5)

p12 = VlnPlot(epithelial.nucs, features = hrs, group.by = "epi_nuc_idents", cols = c("#C00000FF", "#FF0000FF", "purple", "#FFC000FF", "#B8BA00"), pt.size = 0)  +
ggsave(plot = p11, filename = "umap_HR_nuctype_v2.pdf", width = 4, height = 4)


# **Heatmap ------
all_markers_sel2 = all_markers
lvl = levels(umap_tx$major_group)
all_markers_sel2$cluster = factor(all_markers_sel2$cluster, levels = lvl)


top_markers = all_markers_sel2 %>% arrange(cluster) %>% 
  dplyr::group_by(cluster) %>%
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

main_markers = c("KRT5", "EGLN3", "FASN", "PIP", "KRT23", "KIT", "LALBA", "CD74", "TOP2A")

DefaultAssay(tum) <- "integrated"
Idents(tum) <- "final_group"

tum_sub1 = AverageExpression(tum, features = unique(c(main_markers, top_markers)), return.seurat = T, assays = "integrated")
levels(tum_sub1$orig.ident) <- levels(Idents(tum))

p1 = DoHeatmap(tum_sub1, features = c(main_markers, top_markers), assay = "integrated", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = fill_col, raster = F)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
ggsave(plot = p1, filename = "heatmap_cellstate_EPI_allcells_v1.pdf", width = 4, height = 10)


main_markers = c("KRT5", "EGLN3", "FASN", "PIP", "KRT23", "KIT", "LALBA", "CD74", "TOP2A")
Groups_set <- c("basal", "lumhr_egln3", "lumhr_fasn", "lumhr_pip", "lumsec_KRT23","lumsec_kit", "lumsec_lactation", 
                "lumsec_hla", "lumsec_prol")
col_grp <- fill_col[1:9]
names(col_grp) <- Groups_set

df_fibro1 = run_mean_htmaps_fig6_horizontal(df = all_markers_sel2, seu = tum, cell = "b_cells", genes_df = top_markers, #sel_genes
                                            genes_df_cann = main_markers, genes_sel_df = genes_sel_df, 
                                            clus_levels = c("b_naive", "bmem_switched", "bmem_unswitched", "plasma_IgA", "plasma_IgG"))

# **Heatmap indi------
Idents(lumhr) = "final_group"
Idents(lumsec) = "final_group"

lumhr_markers = read.csv("lumhr_markers.csv")
lumsec_markers1 = read.csv("lumsec_markers.csv")
lumsec_markers1$cluster = as.character(lumsec_markers1$cluster)
lumsec_markers1$cluster = factor(lumsec_markers1$cluster, levels = c("lumsec_krt23", "lumsec_kit", "lumsec_lactation",
                                "lumsec_hla", "lumsec_prol"))
basal_markers = all_markers %>% dplyr::filter(cluster == "basal")
  
all_markers_sel2 = all_markers
lvl = levels(clus_out$final_group)
all_markers_sel2$cluster = factor(all_markers_sel2$cluster, levels = lvl)


top_basal_markers = basal_markers %>% 
  filter(!str_detect(gene, "^RP|^HSP")) %>% 
  dplyr::slice_max(n = 7, order_by = avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()
top_lumhr_markers = lumhr_markers %>% arrange(cluster) %>% 
  dplyr::group_by(cluster) %>%
  filter(!str_detect(gene, "^RP|^HSP")) %>% 
  dplyr::slice_max(n = 7, order_by = avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()
top_lumsec_markers = lumsec_markers1 %>% arrange(cluster) %>% 
  dplyr::group_by(cluster) %>%
  filter(!str_detect(gene, "^RP|^HSP")) %>% 
  dplyr::slice_max(n = 7, order_by = avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()

main_markers = c("KRT5", "EGLN3", "FASN", "PIP", "KRT23", "KIT", "LALBA", "CD74", "TOP2A")

DefaultAssay(clus_out) <- "integrated"
Idents(clus_out) <- "final_group"

tum_sub1 = AverageExpression(clus_out, features = unique(c(main_markers, top_basal_markers, top_lumhr_markers, top_lumsec_markers)), return.seurat = T, assays = "RNA")
Idents(tum_sub1)
# tum_sub1$orig.ident = 
# levels(tum_sub1$orig.ident) <- levels(Idents(clus_out))

p1 = DoHeatmap(tum_sub1, features = c(main_markers, top_basal_markers, top_lumhr_markers, top_lumsec_markers), assay = "RNA", slot = "scale.data", label = TRUE,
               disp.min = -1.5, disp.max = 2, draw.lines = F, group.colors = fill_col, raster = F)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5), axis.text.x = element_text(size=6.5))
ggsave(plot = p1, filename = "heatmap_cellstate_EPI_allcells_v2.pdf", width = 4, height = 10)

# **Barplots ---------------------------------------------------------------
fill_col = c("#C00000FF", "#FFBCD9", "#FF5159", "magenta1", "#FFD125", "#A26765",  "#FF8000", "#CA9300", "#B8BA00")
tum = clus_out
df_meta = data.frame(tum@meta.data, 
                     celltype = tum$final_group, 
                     stringsAsFactors = F) %>% 
  as_tibble() %>% clean_names()

df_meta1 = df_meta
# create a summary table
func_bar(var_df = df_meta1, var_x = "sheet_patient_id", var_fill = "celltype", var_pos = "fill", var_cols = fill_col, celltype = "epi_cells")
func_bar(var_df = df_meta1, var_x = "sheet_figure1_grp_updated", var_fill = "celltype", var_pos = "fill", var_cols = fill_col, celltype = "epi_cells")
func_bar(var_df = df_meta1, var_x = "celltype", var_fill = "phase", var_pos = "fill", var_cols = colors_bottle2, celltype = "epi_cells")

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
ggsave(p1, file = paste0(save_path, "/figure6_b_tissue_location_patient.pdf"), width = 7, height = 2.5)

func_bar(var_df = df_meta3, var_x = "key", var_fill = "celltype", var_pos = "fill", var_cols = fill_col, celltype = "b_cells_leftright")
func_bar(var_df = df_meta3, var_x = "sheet_tissue_location", var_fill = "celltype", var_pos = "fill", var_cols = fill_col, celltype = "b_cells_leftright")

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
df_meta2_sel = df_meta1 %>% dplyr::select(final_group1, final_group, final_sample_id, ethnicity, age, sheet_patient_id, sheet_source, sheet_exp_proc, sheet_age_updated, sheet_brca_updated, sheet_parity, sheet_gravida, sheet_menopause, sheet_density, sheet_bmi_final, sheet_cancer_risk,sheet_figure1_grp_updated)
colnames(df_meta2_sel) = c("cellgroup", "celltype", "sample_id", "ethnicity", "age_orig", "patient_id",  "source", "exp_proc", "age", "brca", "parity", "gravida", "menopause", "density", "bmi", "cancer_risk", "surgery")
# df_mrg = df_meta2_sel %>% mutate(pt_x = paste0(patient_id, "_", type))
write_rds(df_meta2_sel, "epi_df_cell_meta2_sel.rds")
df_meta2_sel = read_rds("epi_df_cell_meta2_sel.rds")
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

df_mrg_nuc2$celltype = factor(df_mrg_nuc2$celltype, levels = c("basal", "lumhr_egln3", "lumhr_fasn", "lumhr_pip", "lumsec_krt23", "lumsec_kit","lumsec_lactation", 
                                                               "lumsec_hla", "lumsec_prol"))

write.csv(df_mrg_nuc1, "df_mrg_nuc1_samp_epi.csv")
write.csv(df_mrg_nuc2, "df_mrg_nuc2_pt_epi.csv")

save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/figure_3_v1/epi_all/"
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/figure_3_v1/epi_all/"
col_pal = fill_col

#df_mrg_nuc2$celltype = as.character(df_mrg_nuc2$celltype)
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

# cell cycle analysis Cell---------
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
DefaultAssay(clus_out) <- "RNA"
clus_out <- CellCycleScoring(clus_out, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
clus_out <- CellCycleScoring(lumsec, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
clus_out <- CellCycleScoring(clus_out, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

RidgePlot(clus_out, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2, group.by = "final_group", cols = fill_col)
DimPlot(clus_out, group.by = "Phase", cols = colors_bottle2)
DimPlot(clus_out, group.by = "final_group", cols = fill_col)
FeaturePlot(clus_out, features = c("S.Score", "G2M.Score"))

FeaturePlot(clus_out, features = c("rna_ERBB2", "rna_ERBB2"))

v1 = VlnPlot(clus_out, features = c("S.Score", "G2M.Score"), group.by = "final_group", cols = fill_col, pt.size = 0)
ggsave(v1, file =  "cellcycle_cellstate.pdf", width = 4, height = 3)

# cell cycle analysis Nuclei-------
load("epi_nucs_1.RObj")
epithelial.nucs$epi_nuc_idents = factor(epithelial.nucs$epi_nuc_idents, levels = c("Basal", "Luminal_HR", "Luminal_HR_Prolif",
                                                                                   "Luminal_Sec", "Luminal_Sec_Prolif"))
Idents(epithelial.nucs) <- "epi_nuc_idents"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
DefaultAssay(epithelial.nucs) <- "RNA"
epithelial.nucs <- CellCycleScoring(epithelial.nucs, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

RidgePlot(epithelial.nucs, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2, group.by = "epi_nuc_idents", cols = fill_col)
DimPlot(epithelial.nucs, group.by = "Phase", cols = colors_bottle2)
DimPlot(epithelial.nucs, group.by = "final_group", cols = fill_col)
FeaturePlot(epithelial.nucs, features = c("S.Score", "G2M.Score"))

FeaturePlot(epithelial.nucs, features = c("rna_ERBB2", "rna_ERBB2"))

v1 = VlnPlot(epithelial.nucs, features = c("S.Score", "G2M.Score"), group.by = "epi_nuc_idents", cols = fill_col, pt.size = 0)
ggsave(v1, file =  "cellcycle_cellstate_nuclei.pdf", width = 4, height = 3)
r1 = RidgePlot(epithelial.nucs, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2, group.by = "epi_nuc_idents", cols = fill_col)
ggsave(r1, file =  "genes_prolif_cellstate_nuclei.pdf", width = 9, height = 5)

# cycottrace -----------
install.packages("devtools")
devtools::install_local("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/CytoTRACE/") #stupid pacakge

CytoTRACE()

# plots for Kai --------
## feature/gene plot
f1 = FeatureScatter(clus_out, feature1 = "rna_PGR", 
                    feature2 = "rna_KRT8", group.by = "final_group", cols = c("#C00000FF", "#FFBCD9", "#FF5159", "magenta1", "#FFD125", "#A26765",  "#FF8000", "#CA9300", "#B8BA00"))
f2 = FeatureScatter(basal, feature1 = "rna_PGR", 
                    feature2 = "rna_KRT8", group.by = "final_group", cols = c("#C00000FF"))
f3 = FeatureScatter(lumhr, feature1 = "rna_PGR", 
                    feature2 = "rna_KRT8", group.by = "final_group", cols = c("#FFBCD9", "#FF5159", "magenta1"))
f4 = FeatureScatter(lumsec, feature1 = "rna_PGR", 
                    feature2 = "rna_KRT8", group.by = "final_group", cols = c("#FFD125", "#A26765",  "#FF8000", "#CA9300", "#B8BA00"))
f = cowplot::plot_grid(f1, f2, f3, f4, ncol = 2)
ggsave(f, file =  "krt8_pgr_featureplot.pdf", width = 12, height = 8)

v1 = VlnPlot(basal, features = "rna_KRT14", pt.size = 0, cols = "#C00000FF", group.by = "final_group")
ggsave(v1, file =  "krt14_basal_violinplot.pdf", width = 2.5, height = 3)


v1 = VlnPlot(clus_out, features = "rna_EREG", pt.size = 0, cols = c("#C00000FF", "#FFBCD9", "#FF5159", "magenta1", "#FFD125", "#A26765",  "#FF8000", "#CA9300", "#B8BA00"), group.by = "final_group")
ggsave(v1, file =  "Ereg_epi_violinplot.pdf", width = 5, height = 3)

# % of proliferative cells
df_lumsec = data.frame(lumsec@meta.data, 
                       celltype = lumsec$final_group, 
                       stringsAsFactors = F) %>% 
  as_tibble() %>% clean_names()


df_lumsec1 = df_lumsec %>% dplyr::group_by(sheet_patient_id) %>% 
  mutate(total_lumseccells = n()) %>% ungroup() %>% group_by(sheet_patient_id, final_group) %>% 
  mutate(n_lumsec = n(), prop_lumsec = n_lumsec/total_lumseccells) %>% 
  distinct(sheet_patient_id, final_group, prop_lumsec, n_lumsec, total_lumseccells)
df_lumsec2 = df_lumsec1 %>% dplyr::filter(final_group == "lumsec_prol")
p1 = ggplot(df_lumsec2, aes(x = sheet_patient_id, y = prop_lumsec)) +
  geom_bar(aes(fill = final_group), stat="identity") + 
  geom_text(aes(label = round(prop_lumsec, 2), angle = 90), nudge_y = 0.01,size = 3) + 
  xlab("") + ylab("") + scale_fill_manual(values = c("#B8BA00")) +
  theme(axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90), 
        axis.title.x = element_blank(), 
       axis.ticks = element_blank(), 
        legend.title = element_blank()) + 
  theme(strip.background = element_blank(), panel.spacing.x=unit(0.2, "lines"), strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(2,0,0,0, "pt")))
ggsave(p1, file = paste0(save_path, "/patient_lumsecprol.pdf"), width = 7, height = 2.5)

# % of proliferative cells
df_epi = data.frame(clus_out@meta.data, 
                       celltype = clus_out$final_group, 
                       stringsAsFactors = F) %>% 
  as_tibble() %>% clean_names()


df_epi1 = df_epi %>% dplyr::group_by(sheet_patient_id) %>% 
  mutate(total_epicells = n()) %>% ungroup() %>% group_by(sheet_patient_id, final_group) %>% 
  mutate(n_lumsec = n(), prop_lumsec = n_lumsec/total_lumseccells) %>% 
  distinct(sheet_patient_id, final_group, prop_lumsec, n_lumsec, total_lumseccells)
df_epi12 = df_epi1 %>% dplyr::filter(final_group == "lumsec_prol")
p1 = ggplot(df_epi12, aes(x = sheet_patient_id, y = prop_lumsec)) +
  geom_bar(aes(fill = final_group), stat="identity") + 
  geom_text(aes(label = round(prop_lumsec, 2), angle = 90), nudge_y = 0.002,size = 3) + 
  xlab("") + ylab("") + scale_fill_manual(values = c("#B8BA00")) +
  theme(axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90), 
        axis.title.x = element_blank(), 
        axis.ticks = element_blank(), 
        legend.title = element_blank()) + 
  theme(strip.background = element_blank(), panel.spacing.x=unit(0.2, "lines"), 
        strip.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(2,0,0,0, "pt")));p1
ggsave(p1, file = paste0(save_path, "/patient_lumsecprol_allPi.pdf"), width = 7, height = 2.5)

# HR signatures ---------
fill_col = c("#C00000FF", "#FFBCD9", "#FF5159", "magenta1", "#FFD125", "#A26765",  "#FF8000", "#CA9300", "#B8BA00")

hr_gene = c("ESR1", "PGR", "AR")
DefaultAssay(clus_out) <- "RNA"
clus_out = ScaleData(clus_out)
hr_gene_lst = list(hr_gene)

FeaturePlot(clus_out, features = c("ESR1", "PGR", "AR"))


clus_out <- AddModuleScore(clus_out, features = hr_gene_lst[[1]], assay = "RNA", name = "hormone_score", nbin = 20)
FeaturePlot(clus_out, features = c("hormone_score1"))

v1 = VlnPlot(clus_out, features = c("hormone_score1"), group.by = "final_group", cols = fill_col, pt.size = 0)
ggsave(v1, file =  "hormone_score.pdf", width = 4, height = 3)

# lactation signatures ---------
lac_genes3 = c("CSN2", "LALBA", "CSN3", "CSN1S1")

lac_genes3_lst = list(lac_genes3)
clus_out <- AddModuleScore(clus_out, features = lac_genes3_lst, assay = "RNA", name = "lactation_score_4", nbin = 20)
v3 = VlnPlot(clus_out, features = c("lactation_score_41"), group.by = "final_group", cols = fill_col, pt.size = 0)
ggsave(v3, file =  "lactation_score_4genes.pdf", width = 4, height = 3)


VlnPlot(clus_out, features = c("nCount_RNA"), group.by = "final_group", cols = fill_col, pt.size = 0)

VlnPlot(clus_out, "rna_KRT8", group.by = "final_group", pt.size = 0)

# lumSEC signatures ---------
lumsec3 = c("SLPI", "KRT23", "LTF", "KRT15", "CCL28","PIGR","ALDH1A3")
lumsec3 = c("SLPI", "KRT23", "LTF", "KRT15")
lumsec3 = read.csv("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/figure_1_v1/fig_markers_major_group1.csv")
lumsec3 = lumsec3 %>% dplyr::filter(cluster == "LumSec", avg_logFC>1, p_val_adj<0.05, pct.2<0.10) %>% pull(gene) %>% as.character()
lumsec3_lst = list(lumsec3)
clus_out <- AddModuleScore(clus_out, features = lumsec3_lst, assay = "RNA", name = "lumsec_score_1", nbin = 20)
v6 = VlnPlot(clus_out, features = c("lumsec_score_11"), group.by = "final_group", cols = fill_col, pt.size = 0) + ggtitle("logFC>1pct2<10-19genes")
ggsave(v3, file =  "lumsec_score_logFC>1pct2<30genes.pdf", width = 6, height = 3)
ggsave(v3, file =  "lumsec_score_4genes.pdf", width = 4, height = 3)
p = plot_grid(v3, v4, v5, v6)
ggsave(p, file =  "lumsec_score_mult genes.pdf", width = 12, height = 6)


# basal signatures ---------
basal = c("KRT5", "KRT14", "ACTG2", "KRT17", "LAMA3")
basal = read.csv("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/figure_1_v1/fig_markers_major_group1.csv")
basal = basal %>% dplyr::filter(cluster == "Basal", avg_logFC>1, p_val_adj<0.05, pct.2<0.20) %>% pull(gene) %>% as.character()
basal_lst = list(basal)

clus_out <- AddModuleScore(clus_out, features = basal_lst, assay = "RNA", name = "BASAL_score_1", nbin = 20)
v8 = VlnPlot(clus_out, features = c("BASAL_score_11"), group.by = "final_group", cols = fill_col, pt.size = 0) + ggtitle("logFC>1pct2<20-17genes")
ggsave(v4, file =  "basal_score_5genes.pdf", width = 4, height = 3)
p = plot_grid(v4, v5, v7, v8)
ggsave(p, file =  "basal_score_mult genes.pdf", width = 12, height = 8)

# LIM_MAMMARY_STEM_CELL_UP signature -----------
mammy_genes = read.csv("./kevin_epithelial/LIM_MAMMARY_STEM_CELL_UP.txt")
mammy_genes = mammy_genes$LIM_MAMMARY_STEM_CELL_UP %>% as.character()
mammy_genes = mammy_genes[-1]

mammy_genes3_lst = list(mammy_genes)
clus_out <- AddModuleScore(clus_out, features = mammy_genes3_lst, assay = "RNA", name = "LIM_MAMMARY_STEM_CELL_UP_score_4", nbin = 20)
v3 = VlnPlot(clus_out, features = c("LIM_MAMMARY_STEM_CELL_UP_score_41"), group.by = "final_group", cols = fill_col, pt.size = 0)
ggsave(v3, file =  "LIM_MAMMARY_STEM_CELL_UP.pdf", width = 4, height = 3)

VlnPlot(clus_out, features = c("HLA-C"), group.by = "final_group", cols = fill_col, pt.size = 0)
FeaturePlot(clus_out, "rna_HLA-C")

# PECE_MAMMARY_STEM_CELL_UP signature -----------
mammy_genes = read.csv("./kevin_epithelial/PECE_MAMMARY_STEM_CELL_UP.txt")
mammy_genes = mammy_genes$PECE_MAMMARY_STEM_CELL_UP %>% as.character()
mammy_genes = mammy_genes[-1]

mammy_genes3_lst = list(mammy_genes)
clus_out <- AddModuleScore(clus_out, features = mammy_genes3_lst, assay = "RNA", name = "PECE_MAMMARY_STEM_CELL_UP_score_4", nbin = 20)
v3 = VlnPlot(clus_out, features = c("PECE_MAMMARY_STEM_CELL_UP_score_41"), group.by = "final_group", cols = fill_col, pt.size = 0)
ggsave(v3, file =  "PECE_MAMMARY_STEM_CELL_UP_score_41.pdf", width = 4, height = 3)


# gobp_response_to_progesterone signature -----------
mammy_genes = read.csv("./kevin_epithelial/gobp_response_to_progesterone.txt")
mammy_genes = mammy_genes$GOBP_RESPONSE_TO_PROGESTERONE %>% as.character()
mammy_genes = mammy_genes[-1]

mammy_genes3_lst = list(mammy_genes)
clus_out <- AddModuleScore(clus_out, features = mammy_genes3_lst, assay = "RNA", name = "gobp_response_to_progesterone_score_4", nbin = 20)
v3 = VlnPlot(clus_out, features = c("gobp_response_to_progesterone_score_41"), group.by = "final_group", cols = fill_col, pt.size = 0)
ggsave(v3, file =  "gobp_response_to_progesterone_score_41.pdf", width = 4, height = 3)


# gobp_response_to_estrogen signature -----------
mammy_genes = read.csv("./kevin_epithelial/gobp_response_to_estrogen.txt")
mammy_genes = mammy_genes$GOBP_RESPONSE_TO_ESTROGEN %>% as.character()
mammy_genes = mammy_genes[-1]

mammy_genes3_lst = list(mammy_genes)
clus_out <- AddModuleScore(clus_out, features = mammy_genes3_lst, assay = "RNA", name = "gobp_response_to_estrogen_score_4", nbin = 20)
v3 = VlnPlot(clus_out, features = c("gobp_response_to_estrogen_score_41"), group.by = "final_group", cols = fill_col, pt.size = 0)
ggsave(v3, file =  "gobp_response_to_estrogen_score_41.pdf", width = 4, height = 3)


# REACTOME_PROLACTIN_RECEPTOR_SIGNALING signature -----------
mammy_genes = read.csv("./kevin_epithelial/REACTOME_PROLACTIN_RECEPTOR_SIGNALING.txt")
mammy_genes = mammy_genes$REACTOME_PROLACTIN_RECEPTOR_SIGNALING %>% as.character()
mammy_genes = mammy_genes[-1]

mammy_genes3_lst = list(mammy_genes)
clus_out <- AddModuleScore(clus_out, features = mammy_genes3_lst, assay = "RNA", name = "REACTOME_PROLACTIN_RECEPTOR_SIGNALING_score_4", nbin = 20)
v3 = VlnPlot(clus_out, features = c("REACTOME_PROLACTIN_RECEPTOR_SIGNALING_score_41"), group.by = "final_group", cols = fill_col, pt.size = 0)
ggsave(v3, file =  "REACTOME_PROLACTIN_RECEPTOR_SIGNALING_score_41.pdf", width = 4, height = 3)

# REACTOME_PROLACTIN_RECEPTOR_SIGNALING signature -----------
mammy_genes = read.csv("./kevin_epithelial/WP_MAMMARY_GLAND_DEVELOPMENT_PATHWAY_PREGNANCY_AND_LACTATION_STAGE_3_OF_4.txt")
mammy_genes = mammy_genes$WP_MAMMARY_GLAND_DEVELOPMENT_PATHWAY_PREGNANCY_AND_LACTATION_STAGE_3_OF_4 %>% as.character()
mammy_genes = mammy_genes[-1]

mammy_genes3_lst = list(mammy_genes)
clus_out <- AddModuleScore(clus_out, features = mammy_genes3_lst, assay = "RNA", name = "WP_MAMMARY_GLAND_DEVELOPMENT_PATHWAY_PREGNANCY_AND_LACTATION_STAGE_3_OF_4_score_4", nbin = 20)
v3 = VlnPlot(clus_out, features = c("WP_MAMMARY_GLAND_DEVELOPMENT_PATHWAY_PREGNANCY_AND_LACTATION_STAGE_3_OF_4_score_41"), group.by = "final_group", cols = fill_col, pt.size = 0)
ggsave(v3, file =  "WP_MAMMARY_GLAND_DEVELOPMENT_PATHWAY_PREGNANCY_AND_LACTATION_STAGE_3_OF_4_score_41.pdf", width = 4, height = 3)


# myoep signature -----------
mammy_genes = read.csv("./kevin_epithelial/myoepithelial_markers.xlsx")
mammy_genes = mammy_genes$WP_MAMMARY_GLAND_DEVELOPMENT_PATHWAY_PREGNANCY_AND_LACTATION_STAGE_3_OF_4 %>% as.character()
mammy_genes = mammy_genes[-1]

mammy_genes3_lst = list(mammy_genes)
clus_out <- AddModuleScore(clus_out, features = mammy_genes3_lst, assay = "RNA", name = "WP_MAMMARY_GLAND_DEVELOPMENT_PATHWAY_PREGNANCY_AND_LACTATION_STAGE_3_OF_4_score_4", nbin = 20)
v3 = VlnPlot(clus_out, features = c("WP_MAMMARY_GLAND_DEVELOPMENT_PATHWAY_PREGNANCY_AND_LACTATION_STAGE_3_OF_4_score_41"), group.by = "final_group", cols = fill_col, pt.size = 0)
ggsave(v3, file =  "WP_MAMMARY_GLAND_DEVELOPMENT_PATHWAY_PREGNANCY_AND_LACTATION_STAGE_3_OF_4_score_41.pdf", width = 4, height = 3)


# pathway ------
allcell_markers = read.csv("epi_allcell_markers.csv")
allcell_markers = allcell_markers[,-1]
genes = allcell_markers
genes$cluster = as.character(genes$cluster)

gsea_resultsc2hall = fun_gsea(genes = genes, set = "H", sub_cat = NULL, num_avg_logFC = 0.25)
gsea_resultsc2C5 = fun_gsea(genes = genes, set = "C5", sub_cat = "BP", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "lumsec")
gsea_resultsc2C3 = fun_gsea(genes = genes, set = "C3", sub_cat = "TFT:GTRD", num_avg_logFC = 0.25, celltype = "lumsec", doplot = FALSE, ret_collap = FALSE)
gsea_resultsc2C8 = fun_gsea(genes = genes, set = "C8", sub_cat = NULL, num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "lumsec")
gsea_resultsc2C2_reac = fun_gsea(genes = genes, set = "C2", sub_cat = "CP:REACTOME", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "lumsec")
gsea_resultsc2C2_kegg = fun_gsea(genes = genes, set = "C2", sub_cat = "CP:KEGG", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "lumsec")

clus = c("basal","lumhr_egln3","lumhr_fasn","lumhr_pip","lumsec_hla","lumsec_kit","lumsec_KRT23","lumsec_lactation","lumsec_prol")

gs = gsea_resultsc2hall
gs = gsea_resultsc2C2;gs_kegg = gsea_resultsc2C2_kegg;gs_reac = gsea_resultsc2C2_reac
gs = gs_kegg
gs = gs_reac
gs = gsea_resultsc2C5
gs = gsea_resultsc2C3
gs = gsea_resultsc2C7
gs = gsea_resultsc2C8

gs = rbind(gsea_resultsc2hall, gs_kegg)
gs = rbind(gs, gs_reac)
gs = rbind(gs, gsea_resultsc2C5)
gs_copy = gs
gs = gs_copy
gs = gs %>% dplyr::filter(cluster %in% c("lumhr_egln3","lumhr_fasn","lumhr_pip"))

gs1 = gs %>% dplyr::select(-leadingEdge) %>% dplyr::group_by(cluster) %>% 
  dplyr::filter(padj <= 0.05) %>%
  top_n(10, wt = NES) %>% 
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
colnames(p1.1) = c("pathway", "basal","lumhr_egln3","lumhr_fasn","lumhr_pip","lumsec_hla","lumsec_kit","lumsec_KRT23","lumsec_lactation","lumsec_prol")

# rearrange the datafram
test = p1

# for basal
test1 = test %>% dplyr::arrange(-`basal`) %>% dplyr::filter(basal>0, lumhr_egln3<0, lumhr_fasn<0, lumhr_pip<0, lumsec_hla<0, lumsec_kit<0, lumsec_KRT23<0, lumsec_lactation<0, lumsec_prol<0) %>% pull(pathway)
test2 = test %>% dplyr::filter(!pathway %in% test1) %>% arrange(-lumhr_egln3) %>% dplyr::filter(basal<0, lumhr_egln3>0, lumhr_fasn<0, lumhr_pip<0, lumsec_hla<0, lumsec_kit<0, lumsec_KRT23<0, lumsec_lactation<0, lumsec_prol<0) %>% dplyr::top_n(25, lumhr_egln3) %>% pull(pathway)
test3 = test %>% dplyr::filter(!pathway %in% c(test1, test2)) %>% arrange(-lumhr_fasn) %>% dplyr::filter(basal<0, lumhr_egln3<0, lumhr_fasn>0, lumhr_pip<0, lumsec_hla<0, lumsec_kit<0, lumsec_KRT23<0, lumsec_lactation<0, lumsec_prol<0) %>% dplyr::top_n(25, lumhr_fasn) %>% pull(pathway)
test4 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3)) %>% arrange(-lumhr_pip) %>% dplyr::filter(basal<0, lumhr_egln3<0, lumhr_fasn<0, lumhr_pip>0, lumsec_hla<0, lumsec_kit<0, lumsec_KRT23<0, lumsec_lactation<0, lumsec_prol<0) %>% dplyr::top_n(25, lumhr_pip) %>% pull(pathway)
test5 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3, test4)) %>% arrange(-lumsec_hla) %>% dplyr::filter(basal<0, lumhr_egln3<0, lumhr_fasn<0, lumhr_pip<0, lumsec_hla>0, lumsec_kit<0, lumsec_KRT23<0, lumsec_lactation<0, lumsec_prol<0) %>% dplyr::top_n(25, lumsec_hla) %>% pull(pathway)
test6 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3, test4, test5)) %>% arrange(-lumsec_kit) %>% dplyr::filter(basal<0, lumhr_egln3<0, lumhr_fasn<0, lumhr_pip<0, lumsec_hla<0, lumsec_kit>0, lumsec_KRT23<0, lumsec_lactation<0, lumsec_prol<0) %>% dplyr::top_n(25, lumsec_kit) %>% pull(pathway)
test7 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3, test4, test5, test6)) %>% arrange(-lumsec_KRT23) %>% dplyr::filter(basal<0, lumhr_egln3<0, lumhr_fasn<0, lumhr_pip<0, lumsec_hla<0, lumsec_kit<0, lumsec_KRT23>0, lumsec_lactation<0, lumsec_prol<0) %>% dplyr::top_n(25, lumsec_KRT23) %>% pull(pathway)
test8 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3, test4, test5, test6, test7)) %>% arrange(-lumsec_lactation) %>% dplyr::filter(basal<0, lumhr_egln3<0, lumhr_fasn<0, lumhr_pip<0, lumsec_hla<0, lumsec_kit<0, lumsec_KRT23<0, lumsec_lactation>0, lumsec_prol<0) %>% dplyr::top_n(25, lumsec_lactation) %>% pull(pathway)
test9 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3, test4, test5, test6, test7, test8)) %>% arrange(-lumsec_prol) %>% dplyr::filter(basal<0, lumhr_egln3<0, lumhr_fasn<0, lumhr_pip<0, lumsec_hla<0, lumsec_kit<0, lumsec_KRT23<0, lumsec_lactation<0, lumsec_prol>0) %>% dplyr::top_n(25, lumsec_prol) %>% pull(pathway)
seq_path = c(test1, test2, test3, test4, test5, test6, test7, test8, test9) #

seq_path = test$pathway
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
                  cluster_columns = T,
                  #column_dend_side = "top",
                  #show_column_dend = T,
                  row_dend_side = "right",
                  row_names_side = "left",
                  width = unit(5, "cm"), 
                  row_names_max_width = unit(20, "cm"),
                  column_names_gp = gpar(cex = 0.8));ht_hall

png("epi_htmaps_pathways_unique_onespval0.01.png", res = 300, units = "in", width = 10, height = 38)
print(ht_hall)
dev.off()

png("epi_htmaps_pathways_unique_v1.png", res = 300, units = "in", width = 10, height = 14)
print(ht_hall)
dev.off()

png("epi_htmaps_pathways_top30_padj0.2_v1.png", res = 300, units = "in", width = 10, height = 18)
print(ht_hall)
dev.off()


# pathway individual------
lumsec = read.csv("lumsec_markers.csv")
lumsec = lumsec[,-1]
genes = lumsec
genes$cluster = as.character(genes$cluster)

gsea_resultsc2hall = fun_gsea(genes = genes, set = "H", sub_cat = NULL, num_avg_logFC = 0.25)
gsea_resultsc2C5 = fun_gsea(genes = genes, set = "C5", sub_cat = "BP", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "lumsec")
gsea_resultsc2C3 = fun_gsea(genes = genes, set = "C3", sub_cat = "TFT:GTRD", num_avg_logFC = 0.25, celltype = "lumsec", doplot = FALSE, ret_collap = FALSE)
gsea_resultsc2C8 = fun_gsea(genes = genes, set = "C8", sub_cat = NULL, num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "lumsec")
gsea_resultsc2C2_reac = fun_gsea(genes = genes, set = "C2", sub_cat = "CP:REACTOME", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "lumsec")
gsea_resultsc2C2_kegg = fun_gsea(genes = genes, set = "C2", sub_cat = "CP:KEGG", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "lumsec")

clus = c("lumsec_hla","lumsec_kit","lumsec_krt23","lumsec_lactation","lumsec_prol") #"basal","lumhr_egln3","lumhr_fasn","lumhr_pip",

gs = gsea_resultsc2hall
gs = gsea_resultsc2C2;gs_kegg = gsea_resultsc2C2_kegg;gs_reac = gsea_resultsc2C2_reac
gs = gs_kegg
gs = gs_reac
gs = gsea_resultsc2C5
gs = gsea_resultsc2C3
gs = gsea_resultsc2C7
gs = gsea_resultsc2C8

gs = rbind(gsea_resultsc2hall, gs_kegg)
gs = rbind(gs, gs_reac)
gs = rbind(gs, gsea_resultsc2C5)
gs_copy = gs
gs = gs_copy
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
colnames(p1.1) = c("pathway", "basal","lumhr_egln3","lumhr_fasn","lumhr_pip","lumsec_hla","lumsec_kit","lumsec_KRT23","lumsec_lactation","lumsec_prol")

# rearrange the datafram
test = p1

# for basal
test1 = test %>% dplyr::arrange(-`basal`) %>% dplyr::filter(basal>0, lumhr_egln3<0, lumhr_fasn<0, lumhr_pip<0, lumsec_hla<0, lumsec_kit<0, lumsec_KRT23<0, lumsec_lactation<0, lumsec_prol<0) %>% pull(pathway)
test2 = test %>% dplyr::filter(!pathway %in% test1) %>% arrange(-lumhr_egln3) %>% dplyr::filter(basal<0, lumhr_egln3>0, lumhr_fasn<0, lumhr_pip<0, lumsec_hla<0, lumsec_kit<0, lumsec_KRT23<0, lumsec_lactation<0, lumsec_prol<0) %>% dplyr::top_n(25, lumhr_egln3) %>% pull(pathway)
test3 = test %>% dplyr::filter(!pathway %in% c(test1, test2)) %>% arrange(-lumhr_fasn) %>% dplyr::filter(basal<0, lumhr_egln3<0, lumhr_fasn>0, lumhr_pip<0, lumsec_hla<0, lumsec_kit<0, lumsec_KRT23<0, lumsec_lactation<0, lumsec_prol<0) %>% dplyr::top_n(25, lumhr_fasn) %>% pull(pathway)
test4 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3)) %>% arrange(-lumhr_pip) %>% dplyr::filter(basal<0, lumhr_egln3<0, lumhr_fasn<0, lumhr_pip>0, lumsec_hla<0, lumsec_kit<0, lumsec_KRT23<0, lumsec_lactation<0, lumsec_prol<0) %>% dplyr::top_n(25, lumhr_pip) %>% pull(pathway)
test5 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3, test4)) %>% arrange(-lumsec_hla) %>% dplyr::filter(basal<0, lumhr_egln3<0, lumhr_fasn<0, lumhr_pip<0, lumsec_hla>0, lumsec_kit<0, lumsec_KRT23<0, lumsec_lactation<0, lumsec_prol<0) %>% dplyr::top_n(25, lumsec_hla) %>% pull(pathway)
test6 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3, test4, test5)) %>% arrange(-lumsec_kit) %>% dplyr::filter(basal<0, lumhr_egln3<0, lumhr_fasn<0, lumhr_pip<0, lumsec_hla<0, lumsec_kit>0, lumsec_KRT23<0, lumsec_lactation<0, lumsec_prol<0) %>% dplyr::top_n(25, lumsec_kit) %>% pull(pathway)
test7 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3, test4, test5, test6)) %>% arrange(-lumsec_KRT23) %>% dplyr::filter(basal<0, lumhr_egln3<0, lumhr_fasn<0, lumhr_pip<0, lumsec_hla<0, lumsec_kit<0, lumsec_KRT23>0, lumsec_lactation<0, lumsec_prol<0) %>% dplyr::top_n(25, lumsec_KRT23) %>% pull(pathway)
test8 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3, test4, test5, test6, test7)) %>% arrange(-lumsec_lactation) %>% dplyr::filter(basal<0, lumhr_egln3<0, lumhr_fasn<0, lumhr_pip<0, lumsec_hla<0, lumsec_kit<0, lumsec_KRT23<0, lumsec_lactation>0, lumsec_prol<0) %>% dplyr::top_n(25, lumsec_lactation) %>% pull(pathway)
test9 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3, test4, test5, test6, test7, test8)) %>% arrange(-lumsec_prol) %>% dplyr::filter(basal<0, lumhr_egln3<0, lumhr_fasn<0, lumhr_pip<0, lumsec_hla<0, lumsec_kit<0, lumsec_KRT23<0, lumsec_lactation<0, lumsec_prol>0) %>% dplyr::top_n(25, lumsec_prol) %>% pull(pathway)
seq_path = c(test1, test2, test3, test4, test5, test6, test7, test8, test9) #

seq_path = test$pathway
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
                  cluster_columns = T,
                  #column_dend_side = "top",
                  #show_column_dend = T,
                  row_dend_side = "right",
                  row_names_side = "left",
                  width = unit(5, "cm"), 
                  row_names_max_width = unit(20, "cm"),
                  column_names_gp = gpar(cex = 0.8));ht_hall


png("lumsecepi_htmaps_pathways_top20padj0.05_v1.png", res = 300, units = "in", width = 10, height = 14)
print(ht_hall)
dev.off()

# pathway lumhr------
lumhr = read.csv("lumhr_markers.csv")
lumhr = lumhr[,-1]
genes = lumhr
genes$cluster = as.character(genes$cluster)

gsea_resultsc2hall = fun_gsea(genes = genes, set = "H", sub_cat = NULL, num_avg_logFC = 0.25)
gsea_resultsc2C5 = fun_gsea(genes = genes, set = "C5", sub_cat = "BP", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "lumsec")
gsea_resultsc2C3 = fun_gsea(genes = genes, set = "C3", sub_cat = "TFT:GTRD", num_avg_logFC = 0.25, celltype = "lumsec", doplot = FALSE, ret_collap = FALSE)
gsea_resultsc2C8 = fun_gsea(genes = genes, set = "C8", sub_cat = NULL, num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "lumsec")
gsea_resultsc2C2_reac = fun_gsea(genes = genes, set = "C2", sub_cat = "CP:REACTOME", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "lumsec")
gsea_resultsc2C2_kegg = fun_gsea(genes = genes, set = "C2", sub_cat = "CP:KEGG", num_avg_logFC = 0.25, doplot = FALSE, ret_collap = FALSE, celltype = "lumsec")

clus = c("lumhr_egln3","lumhr_fasn","lumhr_pip") #"basal",,

# gs = gsea_resultsc2hall
gs_kegg = gsea_resultsc2C2_kegg;gs_reac = gsea_resultsc2C2_reac
# gs = gs_kegg
# gs = gs_reac
# gs = gsea_resultsc2C5
# gs = gsea_resultsc2C3
# gs = gsea_resultsc2C7
# gs = gsea_resultsc2C8

gs = rbind(gsea_resultsc2hall, gs_kegg)
gs = rbind(gs, gs_reac)
gs = rbind(gs, gsea_resultsc2C5)
# gs_copy = gs
# gs = gs_copy
#gs = gs %>% dplyr::filter(cluster %in% c("lumhr_egln3","lumhr_fasn","lumhr_pip"))
write_rds(gs, "lumhr_gs.rds")
gs = read_rds("lumhr_gs.rds")

gs1 = gs %>% dplyr::select(-leadingEdge) %>% dplyr::group_by(cluster) %>% 
  dplyr::filter(padj <= 0.05) %>%
  top_n(10, wt = NES) %>% 
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
colnames(p1.1) = c("pathway", "basal","lumhr_egln3","lumhr_fasn","lumhr_pip","lumsec_hla","lumsec_kit","lumsec_KRT23","lumsec_lactation","lumsec_prol")

# rearrange the datafram
test = p1

# # for basal
# test1 = test %>% dplyr::arrange(-`basal`) %>% dplyr::filter(basal>0, lumhr_egln3<0, lumhr_fasn<0, lumhr_pip<0, lumsec_hla<0, lumsec_kit<0, lumsec_KRT23<0, lumsec_lactation<0, lumsec_prol<0) %>% pull(pathway)
# test2 = test %>% dplyr::filter(!pathway %in% test1) %>% arrange(-lumhr_egln3) %>% dplyr::filter(basal<0, lumhr_egln3>0, lumhr_fasn<0, lumhr_pip<0, lumsec_hla<0, lumsec_kit<0, lumsec_KRT23<0, lumsec_lactation<0, lumsec_prol<0) %>% dplyr::top_n(25, lumhr_egln3) %>% pull(pathway)
# test3 = test %>% dplyr::filter(!pathway %in% c(test1, test2)) %>% arrange(-lumhr_fasn) %>% dplyr::filter(basal<0, lumhr_egln3<0, lumhr_fasn>0, lumhr_pip<0, lumsec_hla<0, lumsec_kit<0, lumsec_KRT23<0, lumsec_lactation<0, lumsec_prol<0) %>% dplyr::top_n(25, lumhr_fasn) %>% pull(pathway)
# test4 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3)) %>% arrange(-lumhr_pip) %>% dplyr::filter(basal<0, lumhr_egln3<0, lumhr_fasn<0, lumhr_pip>0, lumsec_hla<0, lumsec_kit<0, lumsec_KRT23<0, lumsec_lactation<0, lumsec_prol<0) %>% dplyr::top_n(25, lumhr_pip) %>% pull(pathway)
# test5 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3, test4)) %>% arrange(-lumsec_hla) %>% dplyr::filter(basal<0, lumhr_egln3<0, lumhr_fasn<0, lumhr_pip<0, lumsec_hla>0, lumsec_kit<0, lumsec_KRT23<0, lumsec_lactation<0, lumsec_prol<0) %>% dplyr::top_n(25, lumsec_hla) %>% pull(pathway)
# test6 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3, test4, test5)) %>% arrange(-lumsec_kit) %>% dplyr::filter(basal<0, lumhr_egln3<0, lumhr_fasn<0, lumhr_pip<0, lumsec_hla<0, lumsec_kit>0, lumsec_KRT23<0, lumsec_lactation<0, lumsec_prol<0) %>% dplyr::top_n(25, lumsec_kit) %>% pull(pathway)
# test7 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3, test4, test5, test6)) %>% arrange(-lumsec_KRT23) %>% dplyr::filter(basal<0, lumhr_egln3<0, lumhr_fasn<0, lumhr_pip<0, lumsec_hla<0, lumsec_kit<0, lumsec_KRT23>0, lumsec_lactation<0, lumsec_prol<0) %>% dplyr::top_n(25, lumsec_KRT23) %>% pull(pathway)
# test8 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3, test4, test5, test6, test7)) %>% arrange(-lumsec_lactation) %>% dplyr::filter(basal<0, lumhr_egln3<0, lumhr_fasn<0, lumhr_pip<0, lumsec_hla<0, lumsec_kit<0, lumsec_KRT23<0, lumsec_lactation>0, lumsec_prol<0) %>% dplyr::top_n(25, lumsec_lactation) %>% pull(pathway)
# test9 = test %>% dplyr::filter(!pathway %in% c(test1, test2, test3, test4, test5, test6, test7, test8)) %>% arrange(-lumsec_prol) %>% dplyr::filter(basal<0, lumhr_egln3<0, lumhr_fasn<0, lumhr_pip<0, lumsec_hla<0, lumsec_kit<0, lumsec_KRT23<0, lumsec_lactation<0, lumsec_prol>0) %>% dplyr::top_n(25, lumsec_prol) %>% pull(pathway)
# seq_path = c(test1, test2, test3, test4, test5, test6, test7, test8, test9) #

seq_path = test$pathway
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
                  cluster_columns = T,
                  #column_dend_side = "top",
                  #show_column_dend = T,
                  row_dend_side = "right",
                  row_names_side = "left",
                  width = unit(5, "cm"), 
                  row_names_max_width = unit(20, "cm"),
                  column_names_gp = gpar(cex = 0.8));ht_hall


png("lumHRepi_htmaps_pathways_top20padj0.05_v1.png", res = 300, units = "in", width = 10, height = 14)
print(ht_hall)
dev.off()

# ***************************** -------------------
# 2. SUPP figures ------------

# ** Full heatmap -----------
# supp mean heatmap for fig 1
DefaultAssay(clus_out) <- "RNA"
clus_out = ScaleData(clus_out)
all_markers = read.csv("./epi_fig_markers_final_group.csv")
all_markers = all_markers[,-1]
lumhr_markers = read.csv("lumhr_markers.csv")
lumsec_markers1 = read.csv("lumsec_markers.csv")
lumsec_markers1$cluster = as.character(lumsec_markers1$cluster)
lumsec_markers1$cluster = factor(lumsec_markers1$cluster, levels = c("lumsec_krt23", "lumsec_kit", "lumsec_lactation",
                                                                     "lumsec_hla", "lumsec_prol"))
basal_markers = all_markers %>% dplyr::filter(cluster == "basal")

all_markers_sel2 = all_markers
lvl = levels(clus_out$final_group)
all_markers_sel2$cluster = factor(all_markers_sel2$cluster, levels = lvl)


top_basal_markers = basal_markers %>% 
  filter(!str_detect(gene, "^RP|^HSP")) %>% 
  dplyr::slice_max(n = 7, order_by = avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()
top_lumhr_markers = lumhr_markers %>% arrange(cluster) %>% 
  dplyr::group_by(cluster) %>%
  filter(!str_detect(gene, "^RP|^HSP")) %>% 
  dplyr::slice_max(n = 7, order_by = avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()
top_lumsec_markers = lumsec_markers1 %>% arrange(cluster) %>% 
  dplyr::group_by(cluster) %>%
  filter(!str_detect(gene, "^RP|^HSP")) %>% 
  dplyr::slice_max(n = 7, order_by = avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()

main_markers = c("KRT5", "EGLN3", "FASN", "PIP", "KRT23", "KIT", "LALBA", "CD74", "TOP2A")

all_markers_sel_fig = data.frame(gene = c(top_basal_markers, top_lumhr_markers, top_lumsec_markers), cluster1 = c(rep("basal", 7), rep("lumhr_egln3", 7),
                                                                                                                   rep("lumhr_fasn", 7), rep("lumhr_pip", 7),
                                                                                                                   rep("lumsec_krt23", 7), rep("lumsec_kit", 7),
                                                                                                                   rep("lumsec_lactation", 7), rep("lumsec_hla", 7),
                                                                                                                   rep("lumsec_prol", 7)))
DefaultAssay(clus_out) <- "integrated"
Idents(clus_out) <- "final_group"
all_markers_sel_fig_genes = unique(c(main_markers, top_basal_markers, top_lumhr_markers, top_lumsec_markers))
col_grp = c("basal" = "#C00000FF", "#FFBCD9", "#FF5159", "magenta1", "#FFD125", "#A26765",  "#FF8000", "#CA9300", "#B8BA00")
names(col_grp) = lvl
df_fibro1 = run_mean_htmaps_fig1_supp(df = all_markers_sel_fig, seu = clus_out, cell = "epi_cells", genes_df = all_markers_sel_fig_genes, #sel_genes
                                      genes_df_cann = main_markers, genes_sel_df = genes_sel_df, 
                                      clus_levels = lvl)
 
# *** violin plots -----------------------
v1 = VlnPlot(clus_out, features = c("LGR5", "PROCR", "ALDH1A3", "ELF5", "KIT"), pt.size = 0, assay = "RNA", group.by = "final_group", cols = fill_col)
ggsave(v1, file =  "vlnplot_progenitor_genes.pdf", width = 8, height = 7)
DefaultAssay(clus_out) <- "RNA"
v2 = FeaturePlot(clus_out, features = c("LGR5", "PROCR", "ALDH1A3"))
ggsave(v2, file =  "feature_progenitor_genes.pdf", width = 8, height = 7)

v3 = VlnPlot(clus_out, features = c("EGLN3", "FASN", "PIP"), pt.size = 0, assay = "RNA", group.by = "final_group", cols = fill_col)
ggsave(v3, file =  "vlnplot_lumhr_genes.pdf", width = 8, height = 7)

# *** feature plots -----------------------
# The genes to plot.
to_plots <- c("EPCAM", "KRT8", "KRT19", "TP63")
to_plots <- c("PROCR", "LGR5", "ALDH1A3")

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
ggsave("feature_canonical_genes_canonical.pdf", plot = plt, width = 2, height = 6, limitsize = FALSE)  
ggsave("feature_progenitor.pdf", plot = plt, width = 2, height = 6, limitsize = FALSE)  

# ** > lumhr ----
to_plots <- c("EGLN3", "PIP", "FASN")
DefaultAssay(lumhr) <- "RNA"
plots <- lapply(to_plots, function(to_plot) {               
  plt <- FeaturePlot(lumhr, to_plot, order = TRUE) + 
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
ggsave("feature_lumhr_canonical.pdf", plot = plt, width = 2, height = 6, limitsize = FALSE)  

# ** > lumsec ----
to_plots <- c("KRT23", "KIT", "LALBA", "CCL20", "PCNA")
DefaultAssay(lumsec) <- "RNA"
plots <- lapply(to_plots, function(to_plot) {               
  plt <- FeaturePlot(lumsec, to_plot, order = TRUE) + 
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
ggsave("feature_lumsec_canonical.pdf", plot = plt, width = 2, height = 8, limitsize = FALSE)  

# ** > basal ----
to_plots <- c("KRT5", "KRT14", "ACTA2", "TP63")
DefaultAssay(basal) <- "RNA"
plots <- lapply(to_plots, function(to_plot) {               
  plt <- FeaturePlot(basal, to_plot, order = TRUE) + 
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
ggsave("feature_basal_canonical.pdf", plot = plt, width = 2, height = 8, limitsize = FALSE)  


# *** barplots for prol cells -------
# cell
tum = clus_out
df_meta = data.frame(tum@meta.data, 
                     celltype = tum$epi_state, 
                     stringsAsFactors = F) %>% 
  as_tibble() %>% clean_names()

df_meta1 = df_meta
df_meta1$final_group1 = as.character(df_meta1$final_group1)
df_meta1 = df_meta1 %>% dplyr::mutate(new_state = if_else(epi_state == "ls_prol", "prol", final_group1))
df_meta1$new_state = factor(df_meta1$new_state, levels = c( "Basal", "LumHR", "LumSec", "prol"))
table(df_meta1$new_state)
table(df_meta1$final_group1)

# create a summary table
fill_col = c("#C00000FF", "#FF0000FF", "#FFC000FF","black")
func_bar(var_df = df_meta1, var_x = "final_group1", var_fill = "new_state", var_pos = "fill", var_cols = fill_col, celltype = "epi_cells")

p1 = ggplot(df_meta1) +
  geom_bar(aes(x = final_group1, fill = new_state), position = "fill") + 
  xlab("") + ylab("") + scale_fill_manual(values = fill_col) + 
  scale_y_continuous(name = "Y-axis",labels = scales::percent) +
  # set axis limits in coord_cartesian
  coord_cartesian(ylim = c(0, 0.10))+ #0.90, 1
  theme(axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
        panel.background = element_blank(),axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(2,0,0,0, "pt")), axis.ticks = element_blank(), 
        legend.title = element_blank())
p1

pdf("p1_inverted.pdf", width = 2.8, height = 3);p1;dev.off()

# nuclei
load("epi_nucs_1.RObj")
Idents(epithelial.nucs) <- "epi_nuc_idents"
df_meta11 = data.frame(epithelial.nucs@meta.data, 
                     celltype = epithelial.nucs$epi_nuc_idents, 
                     stringsAsFactors = F) %>% 
  as_tibble() %>% clean_names()
df_meta11 = df_meta11

# df_meta1$final_group1 = as.character(df_meta1$final_group1)
#df_meta11$epi_nuc_idents = factor(df_meta11$epi_nuc_idents, levels = c("Basal", "Luminal_HR_Prolif", "Luminal_HR", "Luminal_Sec_Prolif", "Luminal_Sec"))
df_meta11$epi_nuc_idents = factor(df_meta11$epi_nuc_idents, levels = c("Basal", "Luminal_HR","Luminal_HR_Prolif",   "Luminal_Sec","Luminal_Sec_Prolif"))
table(df_meta11$epi_nuc_idents)
#fill_col = c("#C00000FF","black", "#FF0000FF", "black", "#FFC000FF")
fill_col = c("#C00000FF", "#FF0000FF","black",  "#FFC000FF","black")
func_bar(var_df = df_meta11, var_x = "epi_new", var_fill = "epi_nuc_idents", var_pos = "fill", var_cols = fill_col, celltype = "epi_nucs")

p11 = ggplot(df_meta11) +
  geom_bar(aes(x = epi_new, fill = epi_nuc_idents), position = "fill") + 
  xlab("") + ylab("") + scale_fill_manual(values = fill_col) + 
  theme(axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
        panel.background = element_blank(),axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(2,0,0,0, "pt")), axis.ticks = element_blank(), 
        legend.title = element_blank()) +
  coord_cartesian(ylim = c(0.0, 0.1))
# ylim(0.5, 1.2)
p11
pdf("p2_inverted.pdf", width = 2.8, height = 3);p11;dev.off()
png("p2_inverted.png", width = 2.8, height = 3, units = "in", res = 300);p11;dev.off()
