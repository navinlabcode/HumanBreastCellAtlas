setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/")
odir = "03_remove_doublets_recluster"; dir.create(odir)
library(future)
plan("multiprocess", workers = 50)
options(future.globals.maxSize = 50000 * 1024^2)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_remove_doublets_recluster/"

# core ---
setwd("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/")
odir = "02_markers"; dir.create(odir)
#plan("multiprocess", workers = 100)
options(future.globals.maxSize = 50000 * 1024^2)
library(future)
plan("multiprocess", workers = 10)
plan(strategy = "multicore", workers = 20)
wd = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20200210_final1/00_load_packages.R")
source("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")
save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20201212_final2/03_remove_doublets_recluster/"


# set variables ----------------------------------------------------------------------------------------------------------------
dims = 20; k.param = 30; cols = colors_dark; reduction = "umap"
fts1 = c("EPCAM", "PTPRC", "PECAM1")

# read object ----------------------------------------------------------------------------------------------------------------
tum.integrated.ns = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/02_cluster_integration_core/dim_20k_30_NULL_tum.integrated.ns2_method_9rcp.rds")

p0 = DimPlot(tum.integrated.ns, group.by = "integrated_snn_res.0.2", cols = colors_dark, label = TRUE)
png(paste0(save_path,"/umap_integrated_snn_res.0.2_method_9.png"), width=6, height=5, units = "in", res = 300)
print(p0)
dev.off()

p00 = DimPlot(tum.integrated.ns, group.by = "patient_id", cols = colors_dark, label = F, shuffle = T)
png(paste0(save_path,"/umap_integrated_snn_res.0.2_pateint_method_9.png"), width=6, height=5, units = "in", res = 300)
print(p00)
dev.off()

# qc -------------------------------------------------------------------------------------------------------------------
p1 = VlnPlot(object = tum.integrated.ns, pt.size = 0, cols = colors_dark, group.by = "integrated_snn_res.0.2", features = c("nCount_RNA", "nFeature_RNA"))
message("aggregating data ...2/6")
tum_count_agg = aggregate(nCount_RNA~integrated_snn_res.0.2, tum.integrated.ns@meta.data, summary)
tum_RNA_agg = aggregate(nFeature_RNA~integrated_snn_res.0.2, tum.integrated.ns@meta.data, summary)
tum_percent.rb_agg = aggregate(percent.rb~integrated_snn_res.0.2, tum.integrated.ns@meta.data, summary)
tum_percent.mito_agg = aggregate(percent.mito~integrated_snn_res.0.2, tum.integrated.ns@meta.data, summary)
tum_percent.sg1_agg = aggregate(percent.sg1~integrated_snn_res.0.2, tum.integrated.ns@meta.data, summary)
tum_percent.sg2_agg = aggregate(percent.sg2~integrated_snn_res.0.2, tum.integrated.ns@meta.data, summary)
tum_percent.sg3_agg = aggregate(percent.sg3~integrated_snn_res.0.2, tum.integrated.ns@meta.data, summary)
tum_summary = cbind(tum_count_agg, tum_RNA_agg, tum_percent.rb_agg, tum_percent.mito_agg, tum_percent.sg1_agg, tum_percent.sg2_agg, tum_percent.sg3_agg)

# save the aggreate summary of each filter
message("saving csv ...3/6")
write.csv(tum_summary, glue("{save_path}/tum_summary_integrated_snn_res.0.2.csv"))


png(paste0(save_path,"/VLNplot_integrated_snn_res.0.2_method_9.png"), width=13, height=6, units = "in", res = 300)
print(p1)
dev.off()

VlnPlot(tum.integrated.ns, features = c("rna_JCHAIN", "rna_CD4"), pt.size = 0,group.by = "integrated_snn_res.0.2",cols = colors_dark)
VlnPlot(tum.integrated.ns, features = c("rna_COL1A1", "rna_EPCAM"), pt.size = 0,group.by = "integrated_snn_res.0.2",cols = colors_dark)
VlnPlot(tum.integrated.ns, features = c("rna_COL1A1", "rna_PTPRC"), pt.size = 0,group.by = "integrated_snn_res.0.2",cols = colors_dark)
VlnPlot(tum.integrated.ns, features = c("rna_KRT5", "rna_PTPRC"), pt.size = 0,group.by = "integrated_snn_res.0.2",cols = colors_dark)
VlnPlot(tum.integrated.ns, features = c("rna_EPCAM", "rna_PTPRC"), pt.size = 0,group.by = "integrated_snn_res.0.2",cols = colors_dark)

tum.integrated.ns$clusters1 = as.character(tum.integrated.ns$integrated_snn_res.0.2)
tum.integrated.ns$clusters1[tum.integrated.ns$integrated_snn_res.0.2 %in% c(0)] = "fibro"
tum.integrated.ns$clusters1[tum.integrated.ns$integrated_snn_res.0.2 %in% c(1)] = "t"
tum.integrated.ns$clusters1[tum.integrated.ns$integrated_snn_res.0.2 %in% c(2)] = "endo_vas"
tum.integrated.ns$clusters1[tum.integrated.ns$integrated_snn_res.0.2 %in% c(3)] = "epi_basal"
tum.integrated.ns$clusters1[tum.integrated.ns$integrated_snn_res.0.2 %in% c(4)] = "epi_lumsec"
tum.integrated.ns$clusters1[tum.integrated.ns$integrated_snn_res.0.2 %in% c(5)] = "epi_lumhr"
tum.integrated.ns$clusters1[tum.integrated.ns$integrated_snn_res.0.2 %in% c(6)] = "pericytes"
tum.integrated.ns$clusters1[tum.integrated.ns$integrated_snn_res.0.2 %in% c(7)] = "fibro"
tum.integrated.ns$clusters1[tum.integrated.ns$integrated_snn_res.0.2 %in% c(8)] = "b"
tum.integrated.ns$clusters1[tum.integrated.ns$integrated_snn_res.0.2 %in% c(9)] = "macro"
tum.integrated.ns$clusters1[tum.integrated.ns$integrated_snn_res.0.2 %in% c(10)] = "fibro_low"
tum.integrated.ns$clusters1[tum.integrated.ns$integrated_snn_res.0.2 %in% c(11)] = "macro_low"
tum.integrated.ns$clusters1[tum.integrated.ns$integrated_snn_res.0.2 %in% c(12)] = "macro"
tum.integrated.ns$clusters1[tum.integrated.ns$integrated_snn_res.0.2 %in% c(13)] = "epi_lumhr_low"
tum.integrated.ns$clusters1[tum.integrated.ns$integrated_snn_res.0.2 %in% c(14)] = "endo_lymph"
tum.integrated.ns$clusters1[tum.integrated.ns$integrated_snn_res.0.2 %in% c(14)] = "endo_lymph"
write_rds(tum.integrated.ns, glue("{save_path}/tum.integrated.ns.rds"))

# try removing bad clusters -------

#tum.integrated.ns_1 = DietSeurat(tum.integrated.ns, counts = TRUE, scale.data = FALSE, dimreducs = TRUE, assays = TRUE)
Idents(tum.integrated.ns) = "integrated_snn_res.0.2"
tum.integrated.ns_1 = subset(tum.integrated.ns, idents = c("15", "16", "17", "18",  "19", "20", "21"), invert = TRUE)
write_rds(tum.integrated.ns_1, glue("{save_path}/tum.integrated.ns_subsetted_1.rds"))

dims = 20; k.param = 31; cols = colors_dark; reduction = "umap"
fts1 = c("nFeature_RNA","nCount_RNA", "percent.mito","percent.rb", "percent.sg1", "percent.sg2", "percent.sg3")
fts2 = c("rna_EPCAM","rna_PTPRC","rna_PECAM1","rna_DCN","rna_KRT5",  "rna_LTF", "rna_KRT17", "rna_HPGD", "rna_TPM2")
fts3 = c("rna_CD69","rna_CD3D","rna_CD68","rna_CD8A","rna_CD4",  "rna_NKG7", "rna_CD19", "rna_MS4A1", "rna_RGS5")
fts4 = c("rna_VWF","rna_ACTA2","rna_CCL21","rna_FAP","rna_COL1A1",  "rna_COL6A1", "rna_ESR1", "rna_KRT8", "rna_SERHL2")

# recluster
plan("multiprocess", workers = 5)
plan(strategy = "multicore", workers = 5)
tum.integrated.ns_1 = FindNeighbors(object = tum.integrated.ns_1, dims = 1:dims, reduction = "pca", k.param = k.param)
tum.integrated.ns_1 = FindClusters(object = tum.integrated.ns_1, resolution = c(1, 0.8, 0.6,0.4,0.3,0.2, 0.1, 0.75, 0.05, 0.02), verbose = TRUE)
write_rds(tum.integrated.ns_1, glue("{save_path}/tum.integrated.ns_1.rds"))
tum.integrated.ns_1 = RunUMAP(object = tum.integrated.ns_1, reduction = "pca", dims = 1:20) #, umap.method = 'umap-learn', metric = 'correlation'
write_rds(tum.integrated.ns_1, glue("{save_path}/dim_{dims}_k_{k.param}_tum.integrated.ns_1.rds"))
seu_object = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_remove_doublets_recluster/tum.integrated.ns_subsetted_1.rds")

DimPlot(tum.integrated.ns_1, cols = colors_dark)

# plots ----
# make plots
p1 = DimPlot(object = tum.integrated.ns_1, reduction = "umap", label = TRUE, pt.size = 0.001, cols = colors_dark, group.by = "clusters1")
#save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20200210_final1//"

png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_all_labelled_method_9.png"), width=6, height=5, units = "in", res = 300)
print(p1)
dev.off()

p1 = DimPlot(object = tum.integrated.ns_1, reduction = reduction, label = TRUE, pt.size = 0.0001, cols = cols, group.by = "clusters1")+  NoAxes()
#p2 = DimPlot(object = tum.integrated.ns_1, reduction = reduction, label = TRUE, pt.size = 0.001, cols = cols, group.by = "clusters2") +  NoAxes()

# system("pdftoppm input.pdf outputname -png")
# save plots
png(paste0(save_path,"/", reduction,dims,"k_", k.param, "umap_all_labelled_method_9.png"), width=6, height=5, units = "in", res = 300)
print(p1)
dev.off()


# find markers ----------------------------------------------------------------------------------------------------------------
# ** Markers function --------
fts1 = c("EPCAM", "PTPRC", "PECAM1")
save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20200210_final1/03_remove_doublets_recluster/"
#save_path = "/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20200210_final1/02_markers/"
char_ident = "clusters1"  #char_ident = "clusters1"
clus_out = read_rds("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_remove_doublets_recluster/tum.integrated.ns_subsetted_1.rds")
mark_out = func_markers_vst(clus_out = clus_out, char_ident = char_ident, fts1 = fts1, dims = dims, k.param = k.param, cols = colors_dark, 
                            save_path = save_path)

markers = read.csv("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/03_remove_doublets_recluster/dims_20k_31_tum.integrated.ns_RNA_markers_intgeratedassay_clusters1_vst_method_9.csv")
clus_out = ScaleData(clus_out, assay = "RNA", vars.to.regress = NULL)
top_markers = markers %>%
  mutate(cluster_1 = factor(cluster, levels = c("b", "endo_lymph", "endo_vas", "epi_basal",
                                                "epi_lumhr", "epi_lumhr_low", "epi_lumsec",
                                                "fibro", "fibro_low", "macro","macro_low",
                                                "pericytes", "t"))) %>%
  arrange(cluster_1) %>% 
  dplyr::group_by(cluster_1) %>%
  top_n(15, avg_logFC) %>%
  dplyr::filter(avg_logFC > 0) %>%
  pull(gene) %>% as.character()

data.frame(new_cluster_factor = clus_out$clusters1, 
           cell = names(clus_out$clusters1)) %>%
  #filter(!str_detect(cell, "untreated")) %>% 
  group_by(new_cluster_factor) %>% 
  sample_n(200) %>% 
  mutate(cell = as.character(cell)) %>% 
  pull(cell) ->
  tum_cells

DefaultAssay(clus_out) <- "RNA"
png(paste0(save_path,"/heatmap_dim_",dims,"k_",k.param,"_tum.integrated.ns_assaymarkers", set_assay, "_integratedassay_", char_ident,"vst_method_9.png"),units = "in", res = 200, width=10,height=10)
png(paste0(save_path,"htmap1.png"),units = "in", res = 200, width=10,height=20)
print(DoHeatmap(clus_out, features = top_markers, size = 2,cells = tum_cells, 
                group.by = "clusters1", slot = "scale.data", disp.min = -1.5, disp.max = 1.5, group.colors = colors_dark) + 
        scale_color_manual(values = colors_dark) +
        scale_fill_gradient2(low="steelblue1", high = "tomato3", mid = "white", midpoint = 0) + 
        theme(axis.text.y.left = element_text(size=8)))
dev.off()

# find de between macro clusters
macro_markers = FindMarkers(clus_out, ident.1 = "macro_low", ident.2 = "macro",assay = "RNA",group.by = "clusters1", logfc.threshold = 0)
write.csv(macro_markers, "macro_markers.csv")
fibro_markers = FindMarkers(clus_out, ident.1 = "fibro_low", ident.2 = "fibro",assay = "RNA",group.by = "clusters1", logfc.threshold = 0)
write.csv(macro_markers, "fibro_markers.csv")
epilumHR_markers = FindMarkers(clus_out, ident.1 = "epi_lumhr_low", ident.2 = "epi_lumhr",assay = "RNA",group.by = "clusters1", logfc.threshold = 0)
write.csv(epilumHR_markers, "epi_lumhr_markers.csv")

P1 = VlnPlot(clus_out, features = c("rna_APOD", "rna_DCN","rna_IGKC"), group.by = "clusters1", pt.size = 0, cols = colors_dark)
P2 = VlnPlot(clus_out, features = c("rna_JCHAIN", "rna_KRT14","rna_IGKC"), group.by = "clusters1", pt.size = 0, cols = colors_dark)
P3 = VlnPlot(clus_out, features = c("rna_JCHAIN", "rna_KRT19","rna_KRT7"), group.by = "clusters1", pt.size = 0, cols = colors_dark)

P = cowplot::plot_grid(P1, P2, P3)
png(paste0(save_path,"vlnplots_low_qc_p1.png"), width=11, height=4, units = "in", res = 300)
print(P1)
dev.off()
png(paste0(save_path,"vlnplots_low_qc_p2.png"), width=11, height=4, units = "in", res = 300)
print(P2)
dev.off()
png(paste0(save_path,"vlnplots_low_qc_p3.png"), width=11, height=4, units = "in", res = 300)
print(P3)
dev.off()

# **Barplots ---------------------------------------------------------------
clus_out = tum.integrated.ns_1
df_meta = data.frame(clus_out@meta.data, 
                     celltype = clus_out$clusters1, 
                     stringsAsFactors = F) %>% 
  tbl_df() %>% clean_names()

#df_meta = df_meta %>% mutate(orig_pt_ident = factor(orig_pt_ident, levels = c("HR+ Luminal Cells", "Secretory Luminal Cells", "Endothelial",  "Fibroblasts",  "Pericytes", "T cells", "Macrophages",
#                                                              "TumA", "TumB", "TumC", "TumD")))
#df_meta1 = df_meta %>% mutate(patient_id =  factor(patient_id, levels = unique(patient_id[order(source)]), ordered = TRUE))
# fn = factor(f, levels=unique(f[order(a,b,f)]), ordered=TRUE)
df_meta1 = df_meta

# create a summary table
func_bar(var_df = df_meta1, var_x = "patient_id", var_fill = "sampletype", var_pos = "fill", var_cols = colors_dark, celltype = "sampletype")
func_bar(var_df = df_meta1, var_x = "patient_id", var_fill = "exp_proc", var_pos = "fill", var_cols = colors_dark[10:17], celltype = "exp_proc")
func_bar(var_df = df_meta1, var_x = "patient_id", var_fill = "procedure", var_pos = "fill", var_cols = colors_dark[4:7], celltype = "procedure")

func_bar(var_df = df_meta1, var_x = "sample_id", var_fill = "sampletype", var_pos = "fill", var_cols = colors_dark, celltype = "source")
func_bar(var_df = df_meta1, var_x = "sample_id", var_fill = "exp_proc", var_pos = "fill", var_cols = colors_dark[10:17], celltype = "exp_proc")
func_bar(var_df = df_meta1, var_x = "sample_id", var_fill = "procedure", var_pos = "fill", var_cols = colors_dark[4:7], celltype = "procedure")

func_bar(var_df = df_meta1, var_x = "patient_id", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters1")
func_bar(var_df = df_meta1, var_x = "patient_id", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters1")
func_bar(var_df = df_meta1, var_x = "sample_id", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters1")
func_bar(var_df = df_meta1, var_x = "sample_id", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters1")
func_bar(var_df = df_meta1, var_x = "procedure", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters1")
func_bar(var_df = df_meta1, var_x = "procedure", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters1")
func_bar(var_df = df_meta1, var_x = "exp_proc", var_fill = "celltype", var_pos = "stack", var_cols = colors_dark, celltype = "clusters1")
func_bar(var_df = df_meta1, var_x = "exp_proc", var_fill = "celltype", var_pos = "fill", var_cols = colors_dark, celltype = "clusters1")

 # END -------
#
# read original object --------------
tum = readRDS(file = paste0(wd, odir,"/dim_30k_20_tum.integrated.ns_RNA_v4_intgeratedassay_clusters1_vst.rds"))
# add meta data ---

table(tum$clusters1)
DimPlot(tum, group.by = "clusters1", label = TRUE, reduction = "umap",cols = colors_dark)

meta_tracker = read_sheet("/volumes/lab/users/tapsi/projects/bcap/HBCA_sample_sheet_for_metadata.xlsx")
meta_tracker = meta_tracker %>% dplyr::select(samp_id, patient_id, ts_sampleid, ts_ptid, tissue_location,  brca_updated, cancer_history, age_updated, gravida.1, parity.1, menopause)

tum@meta.data = left_join(tum@meta.data, meta_tracker, by = c("sample_id" = "samp_id"))

table(is.na(tum@meta.data))

# FALSE 
# 3275259 

# clinial analysis -----
df_meta = data.frame(tum@meta.data,
                     stringsAsFactors = F) %>%
  tbl_df() %>% clean_names()


# # select samples with atleast 30 cells
# sel_samples = df_meta %>% 
#   group_by(sample_id, integrated_snn_res_0_1) %>%
#   summarise(ncells = n()) %>% 
#   dplyr::filter(ncells > 30) %>% pull(sample_id)

# randomly sample 100 cells
df_meta1 = df_meta %>% 
  group_by(sample_id) %>% sample_n(size = 100, replace = TRUE) %>% ungroup()

df_meta_boot = lapply(1:1000, function(i){
  df_meta1 = df_meta %>% 
    group_by(sample_id) %>% 
    sample_n(size = 200, replace = TRUE) %>% 
    ungroup() %>% 
    mutate(iter = i)
}) %>% bind_rows()

df_meta_boot2 = df_meta_boot %>% 
  #filter(new.clonotype_id %in% bottom80_clonos) %>%
  # grouo by sample id and cell state, bootstrap iteration
  group_by(sample_id, integrated_snn_res_0_1, iter) %>%
  summarise(ncells = n()) %>% ungroup() %>% 
  group_by(sample_id, iter) %>%
  mutate(freq = ncells/sum(ncells)) %>% ungroup() %>% 
  left_join(meta_tracker, by = c("sample_id" = "samp_id"))

# View(df_meta_boot)

# barplots ----
df_meta_boot2 %>% 
  dplyr::filter(integrated_snn_res_0_1 == "1") %>% 
  mutate(sample_id2 = fct_reorder(sample_id, as.numeric(factor(tissue_location)))) %>% 
  ggplot(aes(x=sample_id2, y = freq, fill= tissue_location)) +  
  geom_violin() +
  scale_fill_manual(values = colors_rcart) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(axis.ticks.x = element_blank(), 
        #axis.title.y = element_blank(), 
        panel.background = element_blank(), axis.title.x = element_blank(), 
        #axis.text.y = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(1.75,0,0,0, "pt")), axis.ticks = element_blank(), 
        legend.title = element_blank())

df_meta_boot2 %>% 
  # dplyr::filter(integrated_snn_res_0_1 == "1") %>% 
  # mutate(sample_id2 = fct_reorder(sample_id, as.numeric(factor(tissue_location)))) %>% 
  ggplot(aes(x=tissue_location, y = freq, fill= integrated_snn_res_0_1)) +  
  geom_boxplot() +
  scale_fill_manual(values = colors_rcart) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(axis.ticks.x = element_blank(), 
        #axis.title.y = element_blank(), 
        panel.background = element_blank(), axis.title.x = element_blank(), 
        #axis.text.y = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(1.75,0,0,0, "pt")), axis.ticks = element_blank(), 
        legend.title = element_blank())

df_meta_boot2 %>% 
  # dplyr::filter(integrated_snn_res_0_1 == "1") %>% 
  # mutate(sample_id2 = fct_reorder(sample_id, as.numeric(factor(tissue_location)))) %>% 
  ggplot(aes(x=tissue_location, y = freq, fill= integrated_snn_res_0_1)) +  
  geom_violin() +
  scale_fill_manual(values = colors_rcart) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(axis.ticks.x = element_blank(), 
        #axis.title.y = element_blank(), 
        panel.background = element_blank(), axis.title.x = element_blank(), 
        #axis.text.y = element_blank(),
        #axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(1.75,0,0,0, "pt")), axis.ticks = element_blank(), 
        legend.title = element_blank())

func_bar(var_df = df_meta, var_x = "menopause", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells")
func_bar(var_df = df_meta, var_x = "menopause", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells")

func_bar(var_df = df_meta, var_x = "gravida_1", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells")
func_bar(var_df = df_meta, var_x = "gravida_1", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells")

func_bar(var_df = df_meta, var_x = "parity_1", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells")
func_bar(var_df = df_meta, var_x = "parity_1", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells")

func_bar(var_df = df_meta, var_x = "cancer_history", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells")
func_bar(var_df = df_meta, var_x = "cancer_history", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells")

func_bar(var_df = df_meta, var_x = "tissue_location", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells")
func_bar(var_df = df_meta, var_x = "tissue_location", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells")

func_bar(var_df = df_meta, var_x = "brca_updated", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells")
func_bar(var_df = df_meta, var_x = "brca_updated", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells")

func_bar(var_df = df_meta, var_x = "age_updated", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells")
func_bar(var_df = df_meta, var_x = "age_updated", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells")

func_bar(var_df = df_meta, var_x = "exp_proc", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells")
func_bar(var_df = df_meta, var_x = "exp_proc", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells")

func_bar(var_df = df_meta, var_x = "sample_id", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells")
func_bar(var_df = df_meta, var_x = "sample_id", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells")

func_bar(var_df = df_meta, var_x = "sample_id", var_fill = "tissue_location", var_pos = "fill", var_cols = colors_rcart, celltype = "all_cells")
func_bar(var_df = df_meta, var_x = "sample_id", var_fill = "exp_proc", var_pos = "fill", var_cols = colors_set1[c(1,5,3)], celltype = "all_cells")
func_bar(var_df = df_meta, var_x = "sample_id", var_fill = "brca_updated", var_pos = "fill", var_cols = colors_bottle2, celltype = "all_cells")
func_bar(var_df = df_meta, var_x = "sample_id", var_fill = "cancer_history", var_pos = "fill", var_cols = colors_darj, celltype = "all_cells")
func_bar(var_df = df_meta, var_x = "sample_id", var_fill = "age_updated", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells")
func_bar(var_df = df_meta, var_x = "sample_id", var_fill = "gravida_1", var_pos = "fill", var_cols = colors_neon, celltype = "all_cells")
func_bar(var_df = df_meta, var_x = "sample_id", var_fill = "parity_1", var_pos = "fill", var_cols = colors_set1[c(1,6,9)], celltype = "all_cells")
func_bar(var_df = df_meta, var_x = "sample_id", var_fill = "menopause", var_pos = "fill", var_cols = colors_bottle2, celltype = "all_cells")

# clinial analysis digestion type -----
df_meta_short = df_meta %>% dplyr::filter(exp_proc == "short_digestion")
df_meta_long = df_meta %>% dplyr::filter(exp_proc == "overnight_digestion")
df_meta_contra = df_meta %>% dplyr::filter(procedure == "mastectomy")

# short
func_bar(var_df = df_meta_short, var_x = "menopause", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_short")
func_bar(var_df = df_meta_short, var_x = "menopause", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_short")

func_bar(var_df = df_meta_short, var_x = "gravida_1", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_short")
func_bar(var_df = df_meta_short, var_x = "gravida_1", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_short")

func_bar(var_df = df_meta_short, var_x = "parity_1", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_short")
func_bar(var_df = df_meta_short, var_x = "parity_1", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_short")

func_bar(var_df = df_meta_short, var_x = "cancer_history", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_short")
func_bar(var_df = df_meta_short, var_x = "cancer_history", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_short")

func_bar(var_df = df_meta_short, var_x = "tissue_location", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_short")
func_bar(var_df = df_meta_short, var_x = "tissue_location", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_short")

func_bar(var_df = df_meta_short, var_x = "brca_updated", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_short")
func_bar(var_df = df_meta_short, var_x = "brca_updated", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_short")

func_bar(var_df = df_meta_short, var_x = "age_updated", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_short")
func_bar(var_df = df_meta_short, var_x = "age_updated", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_short")

func_bar(var_df = df_meta_short, var_x = "exp_proc", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_short")
func_bar(var_df = df_meta_short, var_x = "exp_proc", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_short")

func_bar(var_df = df_meta_short, var_x = "sample_id", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_short")
func_bar(var_df = df_meta_short, var_x = "sample_id", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_short")

func_bar(var_df = df_meta_short, var_x = "sample_id", var_fill = "tissue_location", var_pos = "fill", var_cols = colors_rcart, celltype = "all_cells_short")
func_bar(var_df = df_meta_short, var_x = "sample_id", var_fill = "exp_proc", var_pos = "fill", var_cols = colors_set1[c(1,5,3)], celltype = "all_cells_short")
func_bar(var_df = df_meta_short, var_x = "sample_id", var_fill = "brca_updated", var_pos = "fill", var_cols = colors_bottle2, celltype = "all_cells_short")
func_bar(var_df = df_meta_short, var_x = "sample_id", var_fill = "cancer_history", var_pos = "fill", var_cols = colors_darj, celltype = "all_cells_short")
func_bar(var_df = df_meta_short, var_x = "sample_id", var_fill = "age_updated", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_short")
func_bar(var_df = df_meta_short, var_x = "sample_id", var_fill = "gravida_1", var_pos = "fill", var_cols = colors_neon, celltype = "all_cells_short")
func_bar(var_df = df_meta_short, var_x = "sample_id", var_fill = "parity_1", var_pos = "fill", var_cols = colors_set1[c(1,6,9)], celltype = "all_cells_short")
func_bar(var_df = df_meta_short, var_x = "sample_id", var_fill = "menopause", var_pos = "fill", var_cols = colors_bottle2, celltype = "all_cells_short")

# contra
func_bar(var_df = df_meta_contra, var_x = "menopause", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_tissue")
func_bar(var_df = df_meta_contra, var_x = "menopause", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_tissue")

func_bar(var_df = df_meta_contra, var_x = "gravida_1", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_tissue")
func_bar(var_df = df_meta_contra, var_x = "gravida_1", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_tissue")

func_bar(var_df = df_meta_contra, var_x = "parity_1", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_tissue")
func_bar(var_df = df_meta_contra, var_x = "parity_1", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_tissue")

func_bar(var_df = df_meta_contra, var_x = "cancer_history", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_tissue")
func_bar(var_df = df_meta_contra, var_x = "cancer_history", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_tissue")

func_bar(var_df = df_meta_contra, var_x = "tissue_location", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_tissue")
func_bar(var_df = df_meta_contra, var_x = "tissue_location", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_tissue")

func_bar(var_df = df_meta_contra, var_x = "brca_updated", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_tissue")
func_bar(var_df = df_meta_contra, var_x = "brca_updated", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_tissue")

func_bar(var_df = df_meta_contra, var_x = "exp_proc", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_tissue")
func_bar(var_df = df_meta_contra, var_x = "exp_proc", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_tissue")

func_bar(var_df = df_meta_contra, var_x = "sample_id", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_tissue")
func_bar(var_df = df_meta_contra, var_x = "sample_id", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_tissue")

func_bar(var_df = df_meta_contra, var_x = "sample_id", var_fill = "tissue_location", var_pos = "fill", var_cols = colors_rcart, celltype = "all_cells_tissue")
func_bar(var_df = df_meta_contra, var_x = "sample_id", var_fill = "exp_proc", var_pos = "fill", var_cols = colors_set1[c(1,5,3)], celltype = "all_cells_tissue")
func_bar(var_df = df_meta_contra, var_x = "sample_id", var_fill = "brca_updated", var_pos = "fill", var_cols = colors_bottle2, celltype = "all_cells_tissue")
func_bar(var_df = df_meta_contra, var_x = "sample_id", var_fill = "cancer_history", var_pos = "fill", var_cols = colors_darj, celltype = "all_cells_tissue")
func_bar(var_df = df_meta_contra, var_x = "sample_id", var_fill = "age_updated", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_tissue")
func_bar(var_df = df_meta_contra, var_x = "sample_id", var_fill = "gravida_1", var_pos = "fill", var_cols = colors_neon, celltype = "all_cells_tissue")
func_bar(var_df = df_meta_contra, var_x = "sample_id", var_fill = "parity_1", var_pos = "fill", var_cols = colors_set1[c(1,6,9)], celltype = "all_cells_tissue")
func_bar(var_df = df_meta_contra, var_x = "sample_id", var_fill = "menopause", var_pos = "fill", var_cols = colors_bottle2, celltype = "all_cells_tissue")

# long
func_bar(var_df = df_meta_long, var_x = "menopause", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_long")
func_bar(var_df = df_meta_long, var_x = "menopause", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_long")

func_bar(var_df = df_meta_long, var_x = "gravida_1", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_long")
func_bar(var_df = df_meta_long, var_x = "gravida_1", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_long")

func_bar(var_df = df_meta_long, var_x = "parity_1", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_long")
func_bar(var_df = df_meta_long, var_x = "parity_1", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_long")

func_bar(var_df = df_meta_long, var_x = "cancer_history", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_long")
func_bar(var_df = df_meta_long, var_x = "cancer_history", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_long")

func_bar(var_df = df_meta_long, var_x = "tissue_location", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_long")
func_bar(var_df = df_meta_long, var_x = "tissue_location", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_long")

func_bar(var_df = df_meta_long, var_x = "brca_updated", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_long")
func_bar(var_df = df_meta_long, var_x = "brca_updated", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_long")

func_bar(var_df = df_meta_long, var_x = "exp_proc", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_long")
func_bar(var_df = df_meta_long, var_x = "exp_proc", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_long")

func_bar(var_df = df_meta_long, var_x = "age_updated", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_long")
func_bar(var_df = df_meta_long, var_x = "age_updated", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_long")

func_bar(var_df = df_meta_long, var_x = "sample_id", var_fill = "clusters1", var_pos = "stack", var_cols = colors_dark, celltype = "all_cells_long")
func_bar(var_df = df_meta_long, var_x = "sample_id", var_fill = "clusters1", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_long")

func_bar(var_df = df_meta_long, var_x = "sample_id", var_fill = "tissue_location", var_pos = "fill", var_cols = colors_rcart, celltype = "all_cells_long")
func_bar(var_df = df_meta_long, var_x = "sample_id", var_fill = "exp_proc", var_pos = "fill", var_cols = colors_set1[c(1,5,3)], celltype = "all_cells_long")
func_bar(var_df = df_meta_long, var_x = "sample_id", var_fill = "brca_updated", var_pos = "fill", var_cols = colors_bottle2, celltype = "all_cells_long")
func_bar(var_df = df_meta_long, var_x = "sample_id", var_fill = "cancer_history", var_pos = "fill", var_cols = colors_darj, celltype = "all_cells_long")
func_bar(var_df = df_meta_long, var_x = "sample_id", var_fill = "age_updated", var_pos = "fill", var_cols = colors_dark, celltype = "all_cells_long")
func_bar(var_df = df_meta_long, var_x = "sample_id", var_fill = "gravida_1", var_pos = "fill", var_cols = colors_neon, celltype = "all_cells_long")
func_bar(var_df = df_meta_long, var_x = "sample_id", var_fill = "parity_1", var_pos = "fill", var_cols = colors_set1[c(1,6,9)], celltype = "all_cells_long")
func_bar(var_df = df_meta_long, var_x = "sample_id", var_fill = "menopause", var_pos = "fill", var_cols = colors_bottle2, celltype = "all_cells_long")
  
  
  