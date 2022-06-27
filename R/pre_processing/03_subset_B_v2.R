# geo ----
setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_B_v2";dir.create(odir)
options(future.globals.maxSize = 90000 * 1024^2)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "03_subset_B_v1"
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_B_v2/"
# stress_sig = read.delim("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_subset_fibroblasts_final1/stromal_stress_signature_032720.txt")
# stress_sig = as.character(stress_sig$stromal.stress.sig)

# core ----
setwd("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_subset_B_v2"; dir.create(odir)
library(future)
plan("multiprocess", workers = 50)
options(future.globals.maxSize = Inf)
setwd(paste0("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
data_folder = "03_subset_B_v1"
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_B_v2/"
# stress_sig = read.delim("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis//03_subset_fibroblasts_final1/stromal_stress_signature_032720.txt")
# stress_sig = as.character(stress_sig$stromal.stress.sig)


# read original object --------------
tum = readRDS(file = paste0(wd, data_folder,"/b_cell_final.rds"))

# subset ---------
Idents(tum) <- "cellstate"
tum_sub = subset(tum, idents = c("b_plasma_hsp"), invert = TRUE) #endo_vas
tum_sub1 = tum_sub
tum_sub@meta.data$cellstate = tum_sub$cellstate

write_rds(tum_sub, path = glue("{save_path}/b_subset_v2.rds"))

tum_sub = RunUMAP(tum_sub, dims = 1:12)

# method 9 ----------
# Log T integration with no scaling
library(future)
plan("multiprocess", workers = 30)
#plan("sequential")

# ** Clus function --------
# run clustering function
seu_object = read_rds(glue("{save_path}/b_subset_v2.rds"))
DefaultAssay(seu_object) <- "RNA"
dims = 12; k.param = 20; cols = colors_dark; var_split = "sheet_exp_proc"; var_scale = "nCount_RNA"
fts1 = c("nFeature_RNA","nCount_RNA", "percent.mito","percent.rb", "percent.stress")
fts2 = c("rna_EPCAM","rna_PTPRC","rna_PECAM1","rna_DCN","rna_KRT5",  "rna_CD86", "rna_COL1A1", "rna_CLU")
fts3 = c("rna_PTPRC","rna_CD79A","rna_CD19","rna_IGKC","rna_IGLC3","rna_IGHG3", "rna_CD79B", "rna_IGHM")
fts4 = c("rna_IRF8","rna_MS4A1","rna_IGHD","rna_CD27","rna_MS4A1", "rna_IGHG1","rna_JCHAIN", "rna_IGHA1","rna_IGHM", "rna_IGLC2", "rna_IGLC3")

#ref_set = c("pt01", "pt06","pt09", "pt20", "pt21","pt25", "pt27", "pt28", "pt31")
ref_set = NULL

clus_out = func_cluster_vst(tum_int = seu_object, var_scale = var_scale, var_split = var_split, dims = dims, k.param = k.param, 
                            markers1 = fts1, markers2 = fts2, markers3 = fts3, markers4 = fts4, save_path = save_path,
                            cols = colors_dark, ref_set = ref_set)

char_ident = "integrated_snn_res.0.4"


# ** Markers function (9)--------
dims = 12; k.param = 20;
char_ident = "integrated_snn_res.0.3"  #integrated_snn_res.0.1
ht_fts1 = c("MS4A1", "JCHAIN", "IGHA1")
#clus_out = read_rds(path = glue("{save_path}/dim_{dims}k_{k.param}_tum.integrated.ns_SCT_v4_intgeratedassay_{char_ident}_vst_method9.rds")) 
mark_out = func_markers_vst(clus_out = clus_out, char_ident = char_ident, fts1 = ht_fts1, dims = dims, k.param = k.param, cols = colors_dark, 
                            save_path = save_path)

clus_out = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_B_v2/dim_12k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")
Idents(clus_out) <- "cellstate"
clus_out$cellstate = factor(clus_out$cellstate, levels = c("b_plasma_igl", "b_plasma_igh_igk", 
                                                 "b_naive_activated_mhc2", "b_naive_cd27", "b_naive_tcl1a", "b_naive_myc"))
markers = FindAllMarkers(clus_out, assay = "RNA", logfc.threshold = 0, only.pos = FALSE)
write.csv(markers, "markers.csv")
#markers = read.csv("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_B_v2/markers.csv")
tum = clus_out

markers_v2 = markers %>% dplyr::filter(gene %in% markers$gene)

markers_v2 %>% 
  mutate(cluster1 = factor(cluster, levels = as.character(levels(tum$cellstate)))) %>%
  arrange(cluster1) %>% 
  filter(!str_detect(gene, "^RP")) %>% 
  filter(!str_detect(gene, "^MT-")) %>% 
  filter(!str_detect(gene, "^HBB")) %>% 
  group_by(cluster1) %>% 
  #filter(pct.2 < 0.50) %>% 
  top_n(10,avg_logFC) %>% 
  pull(gene) %>% 
  unique() %>% as.character() ->
  top_markers

#Idents(tum.integrated.ns)
#tum.integrated.ns$clusters1 <- factor(tum.integrated.ns$clusters1, levels = levels(Idents(tum.integrated.ns)))
message("Making Heatmaps")
Idents(tum) <- "cellstate"
DefaultAssay(tum) <- "integrated"
p1 = DoHeatmap(tum, features = c(top_markers), size = 2,raster = FALSE,
               group.by = "cellstate", slot = "scale.data", disp.min = -1.5, disp.max = 1.5, group.colors = colors_dark) + 
  scale_color_manual(values = colors_dark) +
  scale_fill_gradient2(low="steelblue1", high = "tomato3", mid = "white", midpoint = 0) + 
  theme(axis.text.y.left = element_text(size=6))

#tum$cellstate = factor(tum$cellstate, levels = levels(markers_v2$cluster))
Idents(tum) <- "cellstate"
tavg = AverageExpression(tum, features = unique(top_markers), return.seurat = T)
a1 = DoHeatmap(tavg, features = c(top_markers), assay = "integrated", label = TRUE,
               slot = "scale.data", disp.min = -2, disp.max = 2, draw.lines = F, group.colors = colors_dark, raster = F)  +
  scale_fill_gradient2(low = "#29A7B8", mid = "white", high = "#E85041", na.value = "00000000") + 
  theme(axis.text.y = element_text(size=7.5)) -> tcell_heatmap

png("heatmap_all_labelled_integrated_snn_res.0.6_Tcells_v1.png", width=26, height=26, units = "in", res = 100)
p1
dev.off()
ggsave(plot = a1, filename = "figure_heatmap_cell_states_avg_integratedtop10.png", device = "png", width = 10, height = 14)

u1 = DimPlot(tum, group.by = "cellstate", cols = col_t, label = T)
ggsave(plot = u1, filename = "figure_uamp_cell_states.png", device = "png", width = 8, height = 6)

b = c("rna_MS4A1", "rna_BANK1", "rna_PAX5", "rna_CD24", "rna_CD19", "rna_CR2", "rna_TNFRSF13C")
plasma = c("rna_IGLC2", "rna_IGHA1", "rna_IGHG3", "rna_MZB1", "rna_SDC1", "rna_TNFRSF17", "rna_LEU1", "rna_PRDM1", "rna_CD1D")
prol = c("rna_TOP2A", "rna_MKI67")
naive = c("rna_IGHD", "rna_MS4A1", "rna_NT5E", "rna_RGS13", "rna_CD24", "rna_NEIL1") #NT5E is cd73
memory = c("rna_CD27","rna_IGHD")

igs = c("rna_IGHD","rna_IGHG1","rna_IGHG4","rna_IGHG3","rna_IGHG2","rna_IGHGP","rna_IGHD","rna_IGHA2","rna_IGHA1","rna_IGHD","rna_IGHM", 
        "rna_IGLL5", "rna_IGLC2")

p1=FeaturePlot(tum, reduction = reduction, features = c(b, plasma, prol, naive, memory), order = TRUE, pt.size = 0.001)
p2=FeaturePlot(tum, reduction = reduction, features = c("rna_ACTA2", "rna_PDGFRB", "rna_FAP", "rna_TGFB1", "rna_PDGFA","rna_PDGFB", "rna_FGF2", 
                                                        "rna_S100A8", "rna_VIM", "rna_DDR2", "rna_HGF", "rna_TIMP1", "rna_ADAM10"))
p3=FeaturePlot(tum, reduction = reduction, features = c("rna_ACTA2", "rna_MYH11", "rna_TAGLN", "rna_HHIP", "rna_ASPN",  "rna_JUNB", "rna_SCX"))
p4=FeaturePlot(tum, reduction = reduction, features = c("rna_COL1A1", "rna_CA4", "rna_CD36", "rna_ACKR1", "rna_CREM"))

# add module score
p5=VlnPlot(tum, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), group.by = "cellstate", pt.size = 0)
p6=VlnPlot(tum, features = igs, group.by = "cellstate", pt.size = 0, cols = colors_dark)

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

png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_igs_5violin_9a.png"), width=12, height=15, units = "in", res = 300)
print(p6)
dev.off()


write_rds(tum, "FINAL_B_OBJECT.rds")

# Rename clusters (old)-------
p6 = VlnPlot(tum, features = c("rna_IGLL5", "rna_IGLC2", 
                               "rna_HLA-DRA", "rna_HLA-DRB1",
                               "rna_CD27",
                               "rna_IGHG4", "rna_IGHG2",
                               "rna_TCL1A", "rna_SELL",
                               "rna_MYC", "rna_BCL2A1", "rna_EGR3"), group.by = "cellstate",
             pt.size = 0,cols = colors_dark) #c("#5a189a", "#f4a261","#2a9d8f", "#f15bb5","#fee440", "#00bbf9","#00f5d4")

png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"vln_all.png"), width=9, height=14, units = "in", res = 300)
print(p6)
dev.off()

p7 = DimPlot(tum, group.by = "cellstate", label = T, cols = c("#5a189a", "#f4a261","#2a9d8f", "#f15bb5","#fee440", "#00bbf9","#00f5d4"))
png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_final_9a.png"), width=6, height=3.5, units = "in", res = 300)
print(p7)
dev.off()

png("umap_B_v2.png", width=9, height=7, units = "in", res = 300) #
DimPlot(tum_sub, group.by = "cellstate", label = F, cols = c("#5a189a", "#f4a261","#2a9d8f", "#f15bb5","#fee440", "#00bbf9","#00f5d4"))
dev.off()
DimPlot(tum_sub, group.by = "integrated_snn_res.0.4")

# compare with kerrigan --------
tum = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_subset_B_v2/dim_12k_20_nCount_RNA_sheet_exp_proc_tum.integrated.ns2_method_9.rds")
meta_df = read.csv("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/HCA_immune_metadata_cell_labels.txt",sep = "\t")
meta_df$cellname = rownames(meta_df)
meta_df_b = meta_df %>% dplyr::filter(final_group == "b")

tum1 = tum
meta_df_orig = tum1@meta.data
head(meta_df_orig)

meta_df_mrg = left_join(meta_df_orig, meta_df)
tum1@meta.data = meta_df_mrg
rownames(tum1@meta.data) = tum1@meta.data$cellname
tum1@meta.data$second_pass_cluster_names = as.character(tum1@meta.data$second_pass_cluster_names)
DimPlot(tum1, group.by = "second_pass_cluster_names", cols = colors_dark, label = T)
DimPlot(tum1, group.by = "cellstate", cols = colors_dark, label = T)
DimPlot(tum1, group.by = "cellstate", cols = colors_dark, label = T, split.by = "sheet_exp_proc")
DimPlot(tum1, group.by = "integrated_snn_res.0.3", cols = colors_dark, label = T)
FeaturePlot(tum1, features = c("rna_MS4A1", "rna_JCHAIN", "rna_IGHA1", "rna_IGLL5", "rna_MYC", "rna_CD27", "rna_MKI67"))
table(tum1@meta.data$second_pass_cluster_names)

tum2 = tum
meta_df_orig = tum2@meta.data
head(meta_df_orig)

meta_df_mrg = left_join(meta_df_orig, meta_df_b)
tum2@meta.data = meta_df_mrg
rownames(tum2@meta.data) = tum2@meta.data$cellname
tum2@meta.data$second_pass_cluster_names = as.character(tum2@meta.data$second_pass_cluster_names)
DimPlot(tum2, group.by = "second_pass_cluster_names", cols = colors_dark, label = T)
DimPlot(tum2, group.by = "cellstate", cols = colors_dark, label = T)
DimPlot(tum2, group.by = "integrated_snn_res.0.3", cols = colors_dark, label = T)


# create a confusion matrix
cM <- confusionMatrix(data = tum2@meta.data, paste0(tum2$cellstate), paste0(tum2$second_pass_cluster_names))
print(cor(tum2$cellstate, tum2$second_pass_cluster_names, method = "kendall")) 

cc = table(tum2$cellstate, tum2$second_pass_cluster_names)
g1 = ggplot(as.data.frame(cc), aes(Var1,Var2, fill= Freq)) +
  geom_tile() + geom_text(aes(label=Freq)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  #scale_fill_viridis_c()+
  scale_fill_gradient2(low="white",high="#009194") +
  labs(x = "Tapsi's labels",y = "Kerri's labels")
png(paste0(save_path,"/tapsi_kerri_confusionmat.png"), width=8, height=6, units = "in", res = 300)
print(g1)
dev.off()
