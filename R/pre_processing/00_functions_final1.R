#options(future.globals.maxSize = 50000 * 1024^2)
options(scipen = 999)
library(egg)
library(cowplot)

# ***Filtering -------

#' Function to apply filters to merged seurat object
#'
#' \code{func_filters} returns a filtered seurat object
#'
#' This function applied to merged seurat object and takes in cutoffs for
#' QC metrics. It saves plots pre and post filtering.
#'
#' @param seu_object: merged seurat object
#' @param setident: ident from seurat object
#' @param save_path: data path where files will be saved
#' @param rna_min: lower cutoff for nFeature_RNA, default is 100
#' @param rna_max: upper cutoff for nFeature_RNA, default is 8000
#' @param umi_min: lower cutoff for nCount_RNA, default is 200
#' @param umi_max: upper cutoff for nCount_RNA, default is 40000
#' @param mito_max: upper cutoff for percent of ribo genes, default is 0.50
#' @param rb_max: upper cutoff for percent of mito genes, default is 0.15
#' @param ug_min: upper cutoff for umi to gene ratio, default is 0.25
#' @param ug_max: upper cutoff for umi to gene ratio, default is 1.5
#' @return filtered seurat object
#' @examples
#' tum = func_filters(seu_object = tum, setident = "sample_id", save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20200210_final1/01_read_files/",
#'                   rna_min = 100, rna_max = 8000, umi_min = 200, umi_max = 40000, mito_max = 0.15, rb_max = 0.50, ug_min = 0.25, ug_max = 1.5)

func_filters = function(seu_object = tum, setident = "sample_id", save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis",
                        gene_min = 100, gene_max = 8000, umi_min = 200, umi_max = 40000, mito_max = 0.15, rb_max = 0.50, ug_min = 0.25, ug_max = 1.5){
  
  tum = SetIdent(seu_object, value = setident) # length(tum$sample_id) # 282122
  
  # save plot
  message("saving pre plots ...1/6")
  png(file = glue("{save_path}/pre_filter1.png"), width=20,height=14, units = "in", res = 200)
  print(VlnPlot(object = tum, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "ug.prop"), ncol = 2, pt.size = 0))
  dev.off() 
  
  png(file = glue("{save_path}/pre_filter2.png"), width=20,height=14, units = "in", res = 200)
  print(VlnPlot(object = tum, features = c("percent.rb", "percent.sg1", "percent.sg2", "percent.sg3"), ncol = 2, pt.size = 0))
  dev.off() 
  
  # aggreate the summary of each filter
  message("aggregating data ...2/6")
  tum_count_agg = aggregate(nCount_RNA~sample_id, tum@meta.data, summary)
  tum_RNA_agg = aggregate(nFeature_RNA~sample_id, tum@meta.data, summary)
  tum_percent.rb_agg = aggregate(percent.rb~sample_id, tum@meta.data, summary)
  tum_percent.mito_agg = aggregate(percent.mito~sample_id, tum@meta.data, summary)
  tum_percent.sg1_agg = aggregate(percent.sg1~sample_id, tum@meta.data, summary)
  tum_percent.sg2_agg = aggregate(percent.sg2~sample_id, tum@meta.data, summary)
  tum_percent.sg3_agg = aggregate(percent.sg3~sample_id, tum@meta.data, summary)
  tum_summary = cbind(tum_count_agg, tum_RNA_agg, tum_percent.rb_agg, tum_percent.mito_agg, tum_percent.sg1_agg, tum_percent.sg2_agg, tum_percent.sg3_agg)
  
  # save the aggreate summary of each filter
  message("saving csv ...3/6")
  write.csv(tum_summary, glue("{save_path}/tum_summary.csv"))
  
  # apply filters
  message("subsetting ...4/6")
  
  tum_filt <- subset(x = tum, subset = nCount_RNA > umi_min & nCount_RNA < umi_max) #length(tum_filt$sample_id) #342948
  tum_filt <- subset(x = tum_filt, subset = nFeature_RNA > gene_min & nFeature_RNA < gene_max & percent.mito < mito_max & ug.prop > ug_min & ug.prop < ug_max & percent.rb < rb_max) #262266
  
  #tum_filt <- subset(x = tum, subset = nCount_RNA > 500)
  tum_filt = SetIdent(tum_filt, value = setident) # length(tum_filt$sample_id) # 275621
  
  # save plot
  message("saving plots ...5/6")
  png(file = glue("{save_path}/post_filter1.png"), width=20,height=16, units = "in", res = 200)
  print(VlnPlot(object = tum_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "ug.prop"), ncol = 2, pt.size = 0))
  dev.off() 
  
  png(file = glue("{save_path}/post_filter2.png"), width=20,height=16, units = "in", res = 200)
  print(VlnPlot(object = tum_filt, features = c("percent.rb", "percent.sg1", "percent.sg2", "percent.sg3"), ncol = 2, pt.size = 0))
  dev.off() 
  
  # save file
  message("saving files ...6/6")
  write_rds(tum_filt, glue("{save_path}/tum_filt_v1.rds"))
  return(tum_filt)
}

# example usage
# tum = func_filters(seu_object = tum, setident = "sample_id", save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20200210_final1/01_read_files/",
#                    rna_min = 100, rna_max = 8000, umi_min = 200, umi_max = 40000, mito_max = 0.15, rb_max = 0.50, ug_min = 0.25, ug_max = 1.5)



## ***Clustering < 9 > -----------------
# LogT, intechgration

# initiliaze variables
func_cluster_vst = function(tum_int = seu_object, var_split = var_split, dims = dims, k.param = k.param, markers1 = fts1, markers2 = fts2, 
                            markers3 = fts3, markers4 = fts4,
                            var_scale = var_scale,
                            cols = colors_canv,
                            ref_set = NULL,
                            save_path = save_path){
  
  DefaultAssay(tum_int) <- "RNA"
  message("splitting object....1/8")
  tum.sep.list   <- SplitObject(object = tum_int, split.by = var_split)
  reference.list <- tum.sep.list[names(tum.sep.list)]
  
  #var_regex = '^HLA-|^IG[HJKL]|^RNA|^MT|^RP'
  var_regex = "TR[ABDG][VDJ]|^IG[HJKL]|PLCG2"
  stress_sig = list(c("JUNB", "JUN", "FOS","FOSB", "DNAJB1", "HSPA1A", "HSPA1B", "CXCL8","CXCL1","HMOX1","SAT1","SOD2","CXCL3","AKR1C1","GAPDH","TXN","CXCL2","EDNRB","MEDAG","CD59",    
                  "EIF5A","CLEC2B","TXNRD1","TPI1","SERPINE1", "SLC43A3","PRDX1","PGAM1","WTAP","ADRM1","EIF3I","PRELID1","RAD23A",  
                  "TUBB","PSMB2","ATP1A1","MPZL1","NEAT1","PIM3","AKR1C2","HSPA5"))
  
  message("normalizing object....2/8")
  for (i in 1:length(tum.sep.list)) {
    message(i)
    reference.list[[i]] <- NormalizeData(reference.list[[i]])
    reference.list[[i]] <- AddModuleScore(reference.list[[i]], features = stress_sig, assay = "RNA", name = "stress_mod_rna")
    reference.list[[i]] <- FindVariableFeatures(reference.list[[i]], selection.method = "vst", nfeatures = 5000)
    reference.list[[i]]@assays$RNA@var.features = reference.list[[i]]@assays$RNA@var.features[!str_detect(reference.list[[i]]@assays$RNA@var.features, var_regex)]
    reference.list[[i]]@assays$RNA@var.features = reference.list[[i]]@assays$RNA@var.features[!reference.list[[i]]@assays$RNA@var.features %in% stress_sig]
  }
  

  message("writing object....3/8")
  write_rds(tum.sep.list, glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_tum.sep.list_method_9.rds"))
  write_rds(reference.list, glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_reference.list_method_9.rds")) 
  
  
  #Next, identify anchors and integrate the datasets.
  # i replaced tum_int_50 with tum.sep.list so as to help with low numebr of cells issue
  message("finding anchors....4/8")
  #c("pt01", "pt06","pt09", "pt20", "pt24")) #"pt12", "pt17",
  #reference_dataset <- which(names(reference.list) %in% ref_set)
  #ref_set <- which(names(reference.list) %in%  c("pt41", "pt46","pt49", "pt50", "pt51","pt55", "pt57",  "pt61")) #"pt12", "pt17",
  #ref_set <- c("pt41", "pt46","pt49", "pt50", "pt51","pt55", "pt57",  "pt61") #"pt12", "pt17",
  
  tum.anchors       <- FindIntegrationAnchors(object.list = reference.list, dims = 1:dims, k.filter = k.param, anchor.features = 3000,
                                              k.score = 30, reference = ref_set) #reference = ref_set, 
  write_rds(tum.anchors, glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_tum.anchors_method_9.rds")) 
  to_integrate <- Reduce(intersect, lapply(tum.anchors@object.list, rownames))
  
  message("integrating data....5/8")
  tum.integrated.ns <- IntegrateData(anchorset = tum.anchors, dims = 1:dims, features.to.integrate = NULL) #, features.to.integrate = to_integrate
  write_rds(tum.integrated.ns, glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_{var_split}_tum.integrated.ns1_method_9.rds")) 
  
  # ** SCALE 2
  message("scaling data....6/8")
  DefaultAssay(tum.integrated.ns) <- "integrated"
  #tum.integrated.ns = Seurat::AddModuleScore(tum.integrated.ns, features = stress_sig, assay = "integrated", name = "stress_mod_int")
  #var_scale = "stress_mod_int1"
  #var_scale = "stress_mod_rna1"
  tum.integrated.ns = ScaleData(tum.integrated.ns, assay = "integrated", vars.to.regress = var_scale, features = rownames(tum.integrated.ns@assays$integrated@data)) # no scaling # no scaling
  tum.integrated.ns = RunPCA(tum.integrated.ns, npcs = 100)
  write_rds(tum.integrated.ns, glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_{var_split}_tum.integrated.ns11_method_9.rds")) #tum_int = read_rds(glue("{odir}/dim_{dims}k_{k.param}_tum.integrated.ns11_method_9.rds"))
  e1 = ElbowPlot(tum.integrated.ns, ndims = 30)
  png(paste0(save_path,"elbow", dims,"k_", k.param, "elbow100_method_9.png"), width=8, height=5, units = "in", res = 300)
  print(e1)
  dev.off()
  message("finding neighbors....7/8")
  #tum.integrated.ns = FindNeighbors(object = tum.integrated.ns, dims = 1:dims, reduction = "pca", k.param = k.param)
  #tum.integrated.ns = FindClusters(object = tum.integrated.ns, resolution = c(1, 0.8, 0.6,0.4,0.3,0.2, 0.1, 0.75, 0.05, 0.023), verbose = TRUE)
  tum.integrated.ns = RunUMAP(object = tum.integrated.ns, reduction = "pca", dims = 1:dims) #, umap.method = 'umap-learn', metric = 'manhattan'
  tum.integrated.ns = RunTSNE(object = tum.integrated.ns, reduction = "pca", dims = 1:dims)
  tum.integrated.ns = FindNeighbors(object = tum.integrated.ns, dims = 1:dims, reduction = "pca", k.param = k.param)
  tum.integrated.ns = FindClusters(object = tum.integrated.ns, resolution = c(1, 0.8, 0.6,0.4,0.3,0.2, 0.15,0.17, 0.1, 0.085, 0.090, 0.075, 0.05, 0.023), verbose = TRUE)
  write_rds(tum.integrated.ns, glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_{var_split}_tum.integrated.ns2_method_9.rds")) #tum_int = read_rds(glue("{odir}/dim_{dims}k_{k.param}_tum.integrated.ns_method_9.rds"))
  #tum.integrated.ns = read_rds(glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_tum.integrated.ns2_method_9.rds"))
  
  # ** PLOTS 
  tum.integrated.ns = SetIdent(tum.integrated.ns, value = "integrated_snn_res.0.1")
  message("saving plots....8/8")
  reduction = c("umap", "tsne")
  for (i in 1:2) {
    reduction = c("umap", "tsne")
    DefaultAssay(tum.integrated.ns) <- "RNA"
    reduction = reduction[i]
    message(reduction)
    p1 = DimPlot(object = tum.integrated.ns, reduction = reduction, label = TRUE, pt.size = 0.001, cols = cols, group.by = "integrated_snn_res.0.1")+  NoAxes()
    p2 = DimPlot(object = tum.integrated.ns, reduction = reduction, label = FALSE, pt.size = 0.001, cols = cols, group.by = "integrated_snn_res.0.2") +  NoAxes()
    p3 = DimPlot(object = tum.integrated.ns, reduction = reduction, group.by = "sheet_procedure",  pt.size = 0.00001, cols = c("tan1","steelblue1")) + NoAxes()
    p4 = DimPlot(object = tum.integrated.ns, reduction = reduction, group.by = "batch",  pt.size = 0.0001, cols = c("darkorange","darkorchid")) +  NoAxes()
    p5 = DimPlot(object = tum.integrated.ns, reduction = reduction, group.by = "sheet_patient_id", pt.size = 0.00013, cols = colors_dark) +  NoAxes()
    p6 = DimPlot(object = tum.integrated.ns, reduction = reduction, group.by = "sheet_exp_proc", pt.size = 0.00003, cols = colors_dark) +  NoAxes()
    
    # system("pdftoppm input.pdf outputname -png")
    # save plots
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale, var_split, "umap_all_labelled_method_9.png"), width=6, height=5, units = "in", res = 300)
    print(p1)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale, var_split, "umap_all_method_9.png"), width=6, height=5, units = "in", res = 300)
    print(p2)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,var_split, "umap_procedure_method_9.png"), width=6, height=5, units = "in", res = 300)
    print(p3)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,var_split, "umap_source_method_9.png"), width=6, height=5, units = "in", res = 300)
    print(p4)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,var_split, "umap_patient_method_9.png"), width=6, height=5, units = "in", res = 300)
    print(p5)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,var_split, "umap_exp_proc_method_9.png"), width=6, height=5, units = "in", res = 300)
    print(p6)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,var_split, "umap_markers1_method_9.png"), width=16, height=15, units = "in", res = 300)
    print(FeaturePlot(tum.integrated.ns, features = markers1, reduction = reduction, pt.size = 0.01) + NoAxes())
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,var_split, "umap_markers2_method_9.png"), width=16, height=15, units = "in", res = 300)
    print(FeaturePlot(tum.integrated.ns, features = markers2, reduction = reduction, pt.size = 0.01) + NoAxes())
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,var_split, "umap_markers3_method_9.png"), width=16, height=15, units = "in", res = 300)
    print(FeaturePlot(tum.integrated.ns, features = markers3, reduction = reduction, pt.size = 0.01) + NoAxes())
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,var_split, "umap_markers4_method_9.png"), width=16, height=15, units = "in", res = 300)
    print(FeaturePlot(tum.integrated.ns, features = markers4, reduction = reduction, pt.size = 0.01) + NoAxes())
    dev.off()
    
    # png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_markers4_method_9ki67.png"), width=6, height=5, units = "in", res = 300)
    # print(FeaturePlot(tum.integrated.ns, features = "rna_MKI67", reduction = reduction, pt.size = 0.01) + NoAxes())
    # dev.off()
    
  }
  return(tum.integrated.ns)
}


# initiliaze variables
func_cluster_intsct = function(tum_int = seu_object, var_split = "source", dims = dims, k.param = k.param, markers1 = fts1, markers2 = fts2, 
                            markers3 = fts3, markers4 = fts4,
                            var_scale = var_scale,
                            cols = colors_canv,
                            ref_set = NULL,
                            save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20190504_cellrangerv3_cell/04_subset_Tcells_v2/all"){
  
  DefaultAssay(seu_object) <- "RNA"
  message("splitting object....1/8")
  tum.sep.list   <- SplitObject(object = tum_int, split.by = var_split)
  reference.list <- tum.sep.list[names(tum.sep.list)]
  
  var_regex = '^HLA-|^IG[HJKL]|^RNA|^MT|^RP'
  stress_sig = list(c("CXCL8","CXCL1","HMOX1","SAT1","SOD2","CXCL3","AKR1C1","GAPDH","TXN","CXCL2","EDNRB","MEDAG","CD59",    
                      "EIF5A","CLEC2B","TXNRD1","TPI1","SERPINE1", "SLC43A3","PRDX1","PGAM1","WTAP","ADRM1","EIF3I","PRELID1","RAD23A",  
                      "TUBB","PSMB2","ATP1A1","MPZL1","NEAT1","PIM3","AKR1C2","HSPA5"))
  
  message("normalizing object....2/8")
  for (i in 1:length(tum.sep.list)) {
    message(i)
    reference.list[[i]] <- NormalizeData(reference.list[[i]])
    reference.list[[i]] <- AddModuleScore(reference.list[[i]], features = stress_sig, assay = "RNA", name = "stress_mod_rna")
    reference.list[[i]] <- FindVariableFeatures(reference.list[[i]], selection.method = "vst", nfeatures = 5000)
    reference.list[[i]]@assays$RNA@var.features = reference.list[[i]]@assays$RNA@var.features[!str_detect(reference.list[[i]]@assays$RNA@var.features, var_regex)]
    reference.list[[i]]@assays$RNA@var.features = reference.list[[i]]@assays$RNA@var.features[!reference.list[[i]]@assays$RNA@var.features %in% stress_sig]
  }
  
  
  message("writing object....3/8")
  write_rds(tum.sep.list, glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_tum.sep.list_method_9.rds"))
  write_rds(reference.list, glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_reference.list_method_9.rds")) 
  
  
  #Next, identify anchors and integrate the datasets.
  # i replaced tum_int_50 with tum.sep.list so as to help with low numebr of cells issue
  message("finding anchors....4/8")
  #c("pt01", "pt06","pt09", "pt20", "pt24")) #"pt12", "pt17",
  #reference_dataset <- which(names(reference.list) %in% ref_set)
  #ref_set <- which(names(reference.list) %in%  c("pt41", "pt46","pt49", "pt50", "pt51","pt55", "pt57",  "pt61")) #"pt12", "pt17",
  #ref_set <- c("pt41", "pt46","pt49", "pt50", "pt51","pt55", "pt57",  "pt61") #"pt12", "pt17",
  
  tum.anchors       <- FindIntegrationAnchors(object.list = reference.list, dims = 1:dims, k.filter = k.param, anchor.features = 2500,
                                              k.score = 5) #reference = ref_set, 
  write_rds(tum.anchors, glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_tum.anchors_method_9.rds")) 
  to_integrate <- Reduce(intersect, lapply(tum.anchors@object.list, rownames))
  message("integrating data....5/8")
  tum.integrated.ns <- IntegrateData(anchorset = tum.anchors, dims = 1:dims, features.to.integrate = to_integrate)
  write_rds(tum.integrated.ns, glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_tum.integrated.ns1_method_9.rds")) 
  
  # ** SCALE 2
  message("scaling data....6/8")
  DefaultAssay(tum.integrated.ns) <- "integrated"
  #tum.integrated.ns = Seurat::AddModuleScore(tum.integrated.ns, features = stress_sig, assay = "integrated", name = "stress_mod_int")
  #var_scale = "stress_mod_int1"
  #var_scale = "stress_mod_rna1"
  tum.integrated.ns = ScaleData(tum.integrated.ns, assay = "integrated", vars.to.regress = NULL, features = rownames(tum.integrated.ns@assays$integrated@data)) # no scaling # no scaling
  tum.integrated.ns = RunPCA(tum.integrated.ns, npcs = 100)
  write_rds(tum.integrated.ns, glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_tum.integrated.ns11_method_9.rds")) #tum_int = read_rds(glue("{odir}/dim_{dims}k_{k.param}_tum.integrated.ns11_method_9.rds"))
  e1 = ElbowPlot(tum.integrated.ns, ndims = 30)
  png(paste0(save_path,"elbow", dims,"k_", k.param, "elbow100_method_9.png"), width=8, height=5, units = "in", res = 300)
  print(e1)
  dev.off()
  message("finding neighbors....7/8")
  #tum.integrated.ns = FindNeighbors(object = tum.integrated.ns, dims = 1:dims, reduction = "pca", k.param = k.param)
  #tum.integrated.ns = FindClusters(object = tum.integrated.ns, resolution = c(1, 0.8, 0.6,0.4,0.3,0.2, 0.1, 0.75, 0.05, 0.023), verbose = TRUE)
  tum.integrated.ns = RunUMAP(object = tum.integrated.ns, reduction = "pca", dims = 1:dims) #, umap.method = 'umap-learn', metric = 'manhattan'
  tum.integrated.ns = RunTSNE(object = tum.integrated.ns, reduction = "pca", dims = 1:dims)
  tum.integrated.ns = FindNeighbors(object = tum.integrated.ns, dims = 1:dims, reduction = "pca", k.param = k.param)
  tum.integrated.ns = FindClusters(object = tum.integrated.ns, resolution = c(1, 0.8, 0.6,0.4,0.3,0.2, 0.15,0.17, 0.1, 0.085, 0.090, 0.075, 0.05, 0.023), verbose = TRUE)
  write_rds(tum.integrated.ns, glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_tum.integrated.ns2_method_9.rds")) #tum_int = read_rds(glue("{odir}/dim_{dims}k_{k.param}_tum.integrated.ns_method_9.rds"))
  #tum.integrated.ns = read_rds(glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_tum.integrated.ns2_method_9.rds"))
  
  # ** PLOTS 
  tum.integrated.ns = SetIdent(tum.integrated.ns, value = "integrated_snn_res.0.1")
  message("saving plots....8/8")
  reduction = c("umap", "tsne")
  for (i in 1:2) {
    reduction = c("umap", "tsne")
    DefaultAssay(tum.integrated.ns) <- "RNA"
    reduction = reduction[i]
    message(reduction)
    p1 = DimPlot(object = tum.integrated.ns, reduction = reduction, label = TRUE, pt.size = 0.001, cols = cols, group.by = "integrated_snn_res.0.1")+  NoAxes()
    p2 = DimPlot(object = tum.integrated.ns, reduction = reduction, label = FALSE, pt.size = 0.001, cols = cols, group.by = "integrated_snn_res.0.2") +  NoAxes()
    p3 = DimPlot(object = tum.integrated.ns, reduction = reduction, group.by = "sheet_procedure",  pt.size = 0.00001, cols = c("tan1","steelblue1")) + NoAxes()
    p4 = DimPlot(object = tum.integrated.ns, reduction = reduction, group.by = "batch",  pt.size = 0.0001, cols = c("darkorange","darkorchid")) +  NoAxes()
    p5 = DimPlot(object = tum.integrated.ns, reduction = reduction, group.by = "sheet_patient_id", pt.size = 0.00013, cols = colors_dark) +  NoAxes()
    p6 = DimPlot(object = tum.integrated.ns, reduction = reduction, group.by = "sheet_exp_proc", pt.size = 0.00003, cols = colors_dark) +  NoAxes()
    
    # system("pdftoppm input.pdf outputname -png")
    # save plots
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale, "umap_all_labelled_method_9.png"), width=6, height=5, units = "in", res = 300)
    print(p1)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale, "umap_all_method_9.png"), width=6, height=5, units = "in", res = 300)
    print(p2)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_procedure_method_9.png"), width=6, height=5, units = "in", res = 300)
    print(p3)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_source_method_9.png"), width=6, height=5, units = "in", res = 300)
    print(p4)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_patient_method_9.png"), width=6, height=5, units = "in", res = 300)
    print(p5)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_exp_proc_method_9.png"), width=6, height=5, units = "in", res = 300)
    print(p6)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_markers1_method_9.png"), width=16, height=15, units = "in", res = 300)
    print(FeaturePlot(tum.integrated.ns, features = markers1, reduction = reduction, pt.size = 0.01) + NoAxes())
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_markers2_method_9.png"), width=16, height=15, units = "in", res = 300)
    print(FeaturePlot(tum.integrated.ns, features = markers2, reduction = reduction, pt.size = 0.01) + NoAxes())
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_markers3_method_9.png"), width=16, height=15, units = "in", res = 300)
    print(FeaturePlot(tum.integrated.ns, features = markers3, reduction = reduction, pt.size = 0.01) + NoAxes())
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_markers4_method_9.png"), width=16, height=15, units = "in", res = 300)
    print(FeaturePlot(tum.integrated.ns, features = markers4, reduction = reduction, pt.size = 0.01) + NoAxes())
    dev.off()
    
    # png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_markers4_method_9ki67.png"), width=6, height=5, units = "in", res = 300)
    # print(FeaturePlot(tum.integrated.ns, features = "rna_MKI67", reduction = reduction, pt.size = 0.01) + NoAxes())
    # dev.off()
    
  }
  return(tum.integrated.ns)
}

func_cluster_rpca = function(tum_int = seu_object, var_split = "source", dims = dims, k.param = k.param, markers1 = fts1, markers2 = fts2, 
                            markers3 = fts3, markers4 = fts4,
                            var_scale = var_scale,
                            cols = colors_canv,
                            ref_set = NULL,
                            save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20190504_cellrangerv3_cell/04_subset_Tcells_v2/all"){
  
  message("splitting object....1/8")
  tum.sep.list   <- SplitObject(object = tum_int, split.by = var_split)
  reference.list <- tum.sep.list[names(tum.sep.list)]
  
  var_regex = '^HLA-|^IG[HJKL]|^RNA|^MT|^RP'
  stress_sig = list(c("JUNB", "JUN", "FOS","FOSB", "DNAJB1", "HSPA1A", "HSPA1B", "CXCL8","CXCL1","HMOX1","SAT1","SOD2","CXCL3","AKR1C1","GAPDH","TXN","CXCL2","EDNRB","MEDAG","CD59",    
                      "EIF5A","CLEC2B","TXNRD1","TPI1","SERPINE1", "SLC43A3","PRDX1","PGAM1","WTAP","ADRM1","EIF3I","PRELID1","RAD23A",  
                      "TUBB","PSMB2","ATP1A1","MPZL1","NEAT1","PIM3","AKR1C2","HSPA5"))
  # 
  message("normalizing object....2/8")
  for (i in 1:length(tum.sep.list)) {
    message(i)
    reference.list[[i]] <- NormalizeData(reference.list[[i]])
    reference.list[[i]] <- AddModuleScore(reference.list[[i]], features = stress_sig, assay = "RNA", name = "stress_mod_rna")
    reference.list[[i]] <- FindVariableFeatures(reference.list[[i]], selection.method = "vst", nfeatures = 5000)
    reference.list[[i]]@assays$RNA@var.features = reference.list[[i]]@assays$RNA@var.features[!str_detect(reference.list[[i]]@assays$RNA@var.features, var_regex)]
    reference.list[[i]]@assays$RNA@var.features = reference.list[[i]]@assays$RNA@var.features[!reference.list[[i]]@assays$RNA@var.features %in% stress_sig]
  }
  
  
  message("writing object....3/8")
  write_rds(tum.sep.list, glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_tum.sep.list_method_9rcp.rds"))
  write_rds(reference.list, glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_reference.list_method_9rcp.rds")) 
  
  
  #Next, identify anchors and integrate the datasets.
  # i replaced tum_int_50 with tum.sep.list so as to help with low numebr of cells issue
  message("finding anchors....4/8")
  #c("pt01", "pt06","pt09", "pt20", "pt24")) #"pt12", "pt17",
  #reference_dataset <- which(names(reference.list) %in% ref_set)
  #reference_dataset <- which(names(reference.list) %in%  c("pt01", "pt06","pt09", "pt20", "pt21","pt25", "pt27", "pt28", "pt31")) #"pt12", "pt17",
  
  # for rpca integration
  features <- SelectIntegrationFeatures(object.list = reference.list)
  reference.list <- lapply(X = reference.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = TRUE)
    x <- RunPCA(x, features = features, verbose = TRUE)
  })
  tum.anchors <- FindIntegrationAnchors(object.list = reference.list, reduction = "cca", 
                                        anchor.features = 1500, dims = 1:dims, k.filter = k.param, k.score = 20, reference = ref_set) #,reference = ref_set
  # for rpca integration end
  
  #tum.anchors       <- FindIntegrationAnchors(object.list = reference.list, dims = 1:dims, k.filter = k.param, anchor.features = 2500,
  #                                            k.score = 25) # 
  write_rds(tum.anchors, glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_tum.anchors_method_9rcp.rds")) 
  to_integrate <- Reduce(intersect, lapply(tum.anchors@object.list, rownames))
  message("integrating data....5/8")
  tum.integrated.ns <- IntegrateData(anchorset = tum.anchors, dims = 1:dims, eps = 0.2, features.to.integrate = to_integrate) #, features.to.integrate = to_integrate
  write_rds(tum.integrated.ns, glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_tum.integrated.ns1_method_9rcp.rds")) 
  
  # ** SCALE 2
  message("scaling data....6/8")
  DefaultAssay(tum.integrated.ns) <- "integrated"
  #tum.integrated.ns = Seurat::AddModuleScore(tum.integrated.ns, features = stress_sig, assay = "integrated", name = "stress_mod_int")
  #var_scale = "stress_mod_int1"
  #var_scale = "stress_mod_rna1"
  tum.integrated.ns = ScaleData(tum.integrated.ns, assay = "integrated", vars.to.regress = NULL, features = rownames(tum.integrated.ns@assays$integrated@data)) # no scaling # no scaling
  tum.integrated.ns = RunPCA(tum.integrated.ns, npcs = 100)
  write_rds(tum.integrated.ns, glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_tum.integrated.ns11_method_9rcp.rds")) #tum_int = read_rds(glue("{odir}/dim_{dims}k_{k.param}_tum.integrated.ns11_method_9.rds"))
  e1 = ElbowPlot(tum.integrated.ns, ndims = 30)
  png(paste0(save_path,"elbow", dims,"k_", k.param, "elbow100_method_9rcp.png"), width=8, height=5, units = "in", res = 300)
  print(e1)
  dev.off()
  message("finding neighbors....7/8")
  #tum.integrated.ns = FindNeighbors(object = tum.integrated.ns, dims = 1:dims, reduction = "pca", k.param = k.param)
  #tum.integrated.ns = FindClusters(object = tum.integrated.ns, resolution = c(1, 0.8, 0.6,0.4,0.3,0.2, 0.1, 0.75, 0.05, 0.023), verbose = TRUE)
  tum.integrated.ns = RunUMAP(object = tum.integrated.ns, reduction = "pca", dims = 1:dims, min.dist = 0.1, spread = 1.5) #, umap.method = 'umap-learn', metric = 'manhattan'
  tum.integrated.ns = RunTSNE(object = tum.integrated.ns, reduction = "pca", dims = 1:dims)
  tum.integrated.ns = FindNeighbors(object = tum.integrated.ns, dims = 1:dims, reduction = "pca", k.param = k.param)
  tum.integrated.ns = FindClusters(object = tum.integrated.ns, resolution = c(1, 0.8, 0.6,0.4,0.3,0.2, 0.15,0.17, 0.1, 0.085, 0.090, 0.075, 0.05, 0.023), verbose = TRUE)
  #tum.integrated.ns = FindClusters(object = tum.integrated.ns, resolution = c(1, 1.8, 1.6,1.4,1.3,1.2, 2.15,2.17, 2, 3, 2.5), verbose = TRUE)
  #tum.integrated.ns = FindClusters(object = tum.integrated.ns, resolution = c(3, 3.8, 3.6,3.4,4,5, 2.15,2.17, 2, 3, 2.5, 1.5,1), verbose = TRUE)
  write_rds(tum.integrated.ns, glue("{save_path}/dim_{dims}k_{k.param}_{var_scale}_tum.integrated.ns2_method_9rcp_0.1_1.75.rds")) #tum_int = read_rds(glue("{odir}/dim_{dims}k_{k.param}_tum.integrated.ns_method_9.rds"))
  
  # ** PLOTS 
  tum.integrated.ns = SetIdent(tum.integrated.ns, value = "integrated_snn_res.0.1")
  message("saving plots....8/8")
  reduction = c("umap", "tsne")
  for (i in 1:2) {
    reduction = c("umap", "tsne")
    DefaultAssay(tum.integrated.ns) <- "RNA"
    reduction = reduction[i]
    message(reduction)
    p1 = DimPlot(object = tum.integrated.ns, reduction = reduction, label = TRUE, pt.size = 0.001, cols = cols, group.by = "integrated_snn_res.0.1")+  NoAxes()
    p2 = DimPlot(object = tum.integrated.ns, reduction = reduction, label = FALSE, pt.size = 0.001, cols = cols, group.by = "integrated_snn_res.0.2") +  NoAxes()
    p3 = DimPlot(object = tum.integrated.ns, reduction = reduction, group.by = "procedure",  pt.size = 0.00001, cols = c("tan1","steelblue1")) + NoAxes()
    p4 = DimPlot(object = tum.integrated.ns, reduction = reduction, group.by = "batch",  pt.size = 0.0001, cols = c("darkorange","darkorchid")) +  NoAxes()
    p5 = DimPlot(object = tum.integrated.ns, reduction = reduction, group.by = "patient_id", pt.size = 0.00013, cols = colors_dark) +  NoAxes()
    p6 = DimPlot(object = tum.integrated.ns, reduction = reduction, group.by = "exp_proc", pt.size = 0.00003, cols = colors_dark) +  NoAxes()
    
    # system("pdftoppm input.pdf outputname -png")
    # save plots
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale, "umap_all_labelled_method_9rcp.png"), width=6, height=5, units = "in", res = 300)
    print(p1)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale, "umap_all_method_9rcp.png"), width=6, height=5, units = "in", res = 300)
    print(p2)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_procedure_method_9rcp.png"), width=6, height=5, units = "in", res = 300)
    print(p3)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_source_method_9rcp.png"), width=6, height=5, units = "in", res = 300)
    print(p4)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_patient_method_9rcp.png"), width=10, height=5, units = "in", res = 300)
    print(p5)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_exp_proc_method_9rcp.png"), width=8, height=5, units = "in", res = 300)
    print(p6)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_markers1_method_9rcp.png"), width=18, height=15, units = "in", res = 300)
    print(FeaturePlot(tum.integrated.ns, features = markers1, reduction = reduction, pt.size = 0.01) + NoAxes())
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_markers2_method_9rcp.png"), width=16, height=15, units = "in", res = 300)
    print(FeaturePlot(tum.integrated.ns, features = markers2, reduction = reduction, pt.size = 0.01) + NoAxes())
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_markers3_method_9rcp.png"), width=16, height=15, units = "in", res = 300)
    print(FeaturePlot(tum.integrated.ns, features = markers3, reduction = reduction, pt.size = 0.01) + NoAxes())
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_markers4_method_9rcp.png"), width=16, height=15, units = "in", res = 300)
    print(FeaturePlot(tum.integrated.ns, features = markers4, reduction = reduction, pt.size = 0.01) + NoAxes())
    dev.off()
    
    # png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_markers4_method_9ki67.png"), width=6, height=5, units = "in", res = 300)
    # print(FeaturePlot(tum.integrated.ns, features = "rna_MKI67", reduction = reduction, pt.size = 0.01) + NoAxes())
    # dev.off()
    
  }
  return(tum.integrated.ns)
}


# example run
# dims = 30; k.param = 30
# clus_out = read_rds({save_path}/dim_{dims}k_{k.param}_tum.integrated.ns_SEL_v4.rds"))
# mark_out = func_markers_vst(clus_out = clus_out, char_ident = "integrated_snn_res.0.1", dims = dims, k.param = k.param, save_path = "/volumes/lab/users/tapsi/projects/pdac/analysis/20190410_seu3_all26samp_v1/06_2_macro/all")

# ******markers ----
func_markers_vst = function(clus_out = clus_out, char_ident = char_ident, dims = dims, fts1 = fts1, 
                            k.param = k.param, 
                            cols = colors_canv, 
                            save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20190504_cellrangerv3_cell/04_subset_Tcells_v2/all"){
  
  tum.integrated.ns = SetIdent(clus_out, value = char_ident)
  tum.integrated.ns = ScaleData(tum.integrated.ns, assay = "RNA", vars.to.regress = NULL) # no scaling
  
  char_assays = c("RNA", "integrated")
  for (i in 1:length(char_assays)) {
    set_assay = char_assays[i]
    DefaultAssay(tum.integrated.ns) <- set_assay
    markers = FindAllMarkers(tum.integrated.ns, logfc.threshold = 0)
    write.csv(markers, glue("{save_path}/dims_{dims}k_{k.param}_tum.integrated.ns_{set_assay}_markers_intgeratedassay_{char_ident}_vst_method_9.csv"))
    write_rds(tum.integrated.ns, glue("{save_path}/dim_{dims}k_{k.param}_tum.integrated.ns_{set_assay}_v4_intgeratedassay_{char_ident}_vst_method_9.rds"))
    
   # tum.integrated.ns = read_rds("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20200210_final1/03_subset_epi_v1/dim_25k_31_tum.integrated.ns_RNA_v4_intgeratedassay_integrated_snn_res.0.8_vst_method_9.rds")
  # markers = read.csv("/volumes/lab/users/tapsi/projects/bcap/analysis/20200210_final1/03_remove_doublets_recluster/dims_30k_20_tum.integrated.ns_RNA_markers_intgeratedassay_clusters1_vst.csv")
    #markers = read.csv("/volumes/core/users/tapsi/geo_tapsi/projects/bcap/analysis/20200210_final1/03_subset_epi_v1/dims_25k_31_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.8_vst_method_9.csv")
    
    # heatmap
    # FindAllMarkers usually uses data slot in the RNA assay to find differential genes.
    # For a heatmap or dotplot of markers, the scale.data in the RNA assay should be used.
    # markers_v1 = markers$gene[!str_detect(markers$gene, "^MT-")]
    # markers_v1 = markers_v1[!str_detect(markers_v1, "^(RP)|(MRP)")]
    # markers_v1 = markers_v1[!str_detect(markers_v1, "^HBB")]
    markers_v2 = markers %>% dplyr::filter(gene %in% markers$gene)
    
    #tum.integrated.ns = NormalizeData(tum.integrated.ns, assay = "RNA")
    #levels(markers_v2$cluster) = c("M1", "M2", "M3")
    #markers_v2 = markers_v2 %>% dplyr::mutate(cluster_1 = factor(x = cluster, levels = c("M1", "M2", "M3"), ordered = TRUE))
    #markers_v2 = markers_v2 %>% dplyr::mutate(cluster = factor(x = cluster, levels = unique(tum.integrated.ns$clusters1), ordered = TRUE))
    top_markers = markers_v2 %>%
      #mutate(cluster_1 = factor(cluster, levels = c("M1", "M2", "M3"))) %>%
      dplyr::group_by(cluster) %>%
      #arrange(avg_logFC) %>% 
      top_n(7, avg_logFC) %>%
      dplyr::filter(avg_logFC > 0) %>%
      pull(gene) %>% as.character()
    
    #Idents(tum.integrated.ns)
    #tum.integrated.ns$clusters1 <- factor(tum.integrated.ns$clusters1, levels = levels(Idents(tum.integrated.ns)))
    message("Making Heatmaps")
    DefaultAssay(tum.integrated.ns) <- "RNA"
    png(paste0(save_path,"/heatmap_dim_",dims,"k_",k.param,"_tum.integrated.ns_assaymarkers", set_assay, "_RNAassay_", char_ident,"vst_method_9.png"),  width=20, height=22, units = "in", res = 200) #units = "in", res = 200,
    print(DoHeatmap(tum.integrated.ns, features = c(fts1, top_markers), size = 2,
                    group.by = char_ident, slot = "scale.data", disp.min = -1.5, disp.max = 1.5, group.colors = colors_dark) + 
            scale_color_manual(values = colors_dark) +
            scale_fill_gradient2(low="steelblue1", high = "tomato3", mid = "white", midpoint = 0) + 
            theme(axis.text.y.left = element_text(size=6)))
    dev.off()
    
    DefaultAssay(tum.integrated.ns) <- "integrated"
    png(paste0(save_path,"/heatmap_dim_",dims,"k_",k.param,"_tum.integrated.ns_assaymarkers", set_assay, "_integratedassay_", char_ident,"vst_method_9.png"),units = "in", res = 200, width=20, height=22)
    print(DoHeatmap(tum.integrated.ns, features = c(fts1, top_markers), size = 2,
                    group.by = char_ident, slot = "scale.data", disp.min = -1.5, disp.max = 1.5, group.colors = colors_dark) + 
            scale_color_manual(values = colors_dark) +
            scale_fill_gradient2(low="steelblue1", high = "tomato3", mid = "white", midpoint = 0) + 
            theme(axis.text.y.left = element_text(size=6)))
    dev.off()
    
    message("Making featureplots")
    reduction = c("umap", "tsne")
    for (i in 1:2) {
      reduction = c("umap", "tsne")
      DefaultAssay(tum.integrated.ns) <- set_assay
      reduction = reduction[i]
      message(reduction)
      png(paste0(save_path,"/", reduction, dims,"k_", k.param, "_tum.integrated.ns", set_assay, "_integratedassay_", char_ident,"vst_method_9.png"), width=6,height=5,units = "in", res = 300)
      print(DimPlot(object = tum.integrated.ns, reduction = reduction, label = TRUE, pt.size = 0.002, cols = cols, group.by = "integrated_snn_res.0.2", shuffle = T) +  NoAxes())
      print(DimPlot(object = tum.integrated.ns, reduction = reduction, label = FALSE, pt.size = 0.001, cols = cols, group.by = char_ident, shuffle = T) +  NoAxes())
      dev.off()
    }
    
  }
}

# example run
# dims = 9; k.param = 20
# clus_out = read_rds({save_path}/dim_{dims}k_{k.param}_tum.integrated.ns_SEL_v4.rds"))
# mark_out = func_markers_vst(clus_out = clus_out, char_ident = "integrated_snn_res.0.1", dims = dims, k.param = k.param, save_path = "/volumes/lab/users/tapsi/projects/pdac/analysis/20190410_seu3_all26samp_v1/06_2_macro/all")

## ***Clustering < 13 >-----------------
# SCTransform, scale by a variable, integration

# initiliaze variables
func_cluster_scale = function(tum_int = seu_object, var_split = "source", dims = dims, k.param = k.param, markers1 = fts1, markers2 = fts2, 
                              markers3 = fts3, markers4 = fts4,cols = colors_canv,
                              var_scale = "nFeature_RNA",
                              save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20190504_cellrangerv3_cell/04_subset_Tcells_v2/all")
{
  
  message("splitting object....1/8")
  tum.sep.list   <- SplitObject(object = seu_object, split.by = var_split)
  reference.list <- tum.sep.list[names(tum.sep.list)]
  
  var_regex = '^HLA-|^IG[HJKL]|^RNA|^MT|^RP'
  stress_sig = list(c("CXCL8","CXCL1","HMOX1","SAT1","SOD2","CXCL3","AKR1C1","GAPDH","TXN","CXCL2","EDNRB","MEDAG","CD59",    
                 "EIF5A","CLEC2B","TXNRD1","TPI1","SERPINE1", "SLC43A3","PRDX1","PGAM1","WTAP","ADRM1","EIF3I","PRELID1","RAD23A",  
                 "TUBB","PSMB2","ATP1A1","MPZL1","NEAT1","PIM3","AKR1C2","HSPA5"))
  #var_scale = "stress_mod_rna1"
  for (i in 1:length(tum.sep.list)) {
    message("SCT object....2/8")
    tum.sep.list[[i]] <- AddModuleScore(tum.sep.list[[i]], features = stress_sig, assay = "RNA", name = "stress_mod_rna", ctrl = 5)
    tum.sep.list[[i]] <- SCTransform(tum.sep.list[[i]], verbose = TRUE, vars.to.regress = var_scale)
    tum.sep.list[[i]]@assays$SCT@var.features = tum.sep.list[[i]]@assays$SCT@var.features[!str_detect(tum.sep.list[[i]]@assays$SCT@var.features, var_regex)] # remove vtcr variable genes
  }
  
  message("writing object....3/8")
  write_rds(reference.list, glue("{save_path}/dim_{dims}k_{k.param}_split_by_{var_split}_scaled_by_{var_scale}_reference.list_method_13.rds"))
  write_rds(tum.sep.list, glue("{save_path}/dim_{dims}k_{k.param}_split_by_{var_split}_scaled_by_{var_scale}_tum.sep.list_method_13.rds"))
  
  #Next, select features for downstream integration, and run PrepSCTIntegration, which ensures that all necessary Pearson residuals have been calculated.
  pancreas.features <- SelectIntegrationFeatures(object.list = tum.sep.list, nfeatures = 5000)
  pancreas.features = pancreas.features[!str_detect(pancreas.features, var_regex)]
  #pancreas.features = pancreas.features[!(pancreas.features %in% stress_sig)]
  tum.sep.list      <- PrepSCTIntegration(object.list = tum.sep.list, anchor.features = pancreas.features, 
                                          verbose = TRUE)
  write_rds(tum.sep.list, glue("{save_path}/dim_{dims}k_{k.param}_split_by_{var_split}_scaled_by_{var_scale}_tum.sep.list1_method_13.rds"))
  # sel_cells = names(table(tum_int$patient_id)[table(tum_int$patient_id) > 5])
  # tum_int_50 = tum.sep.list[names(tum.sep.list) %in% sel_cells]
  # 
  #Next, identify anchors and integrate the datasets.
  # i replaced tum_int_50 with tum.sep.list so as to help with low numebr of cells issue
  message("finding anchors....4/8")
  #reference_dataset <- which(names(reference.list) %in%  c("pt01", "pt06","pt09", "pt12", "pt17", "pt20", "pt24"))
  tum.anchors       <- FindIntegrationAnchors(object.list = tum.sep.list, dims = 1:dims, k.filter = k.param, normalization.method = "SCT", anchor.features = pancreas.features, reference = NULL)
  write_rds(tum.anchors, glue("{save_path}/dim_{dims}k_{k.param}_split_by_{var_split}_scaled_by_{var_scale}_tum.anchors_method_13.rds"))
  to_integrate <- Reduce(intersect, lapply(tum.anchors@object.list, rownames))
  
  message("integrating data....5/8")
  tum.integrated.ns <- IntegrateData(anchorset = tum.anchors, normalization.method = "SCT", dims = 1:dims, features.to.integrate = to_integrate)
  
  write_rds(tum.integrated.ns, glue("{save_path}/dim_{dims}k_{k.param}_split_by_{var_split}_scaled_by_{var_scale}_tum.integrated.ns1_method_13.rds"))
  
  # ** SCALE 2 
  message("scaling data....6/8")
  DefaultAssay(tum.integrated.ns) <- "integrated"
  tum.integrated.ns = RunPCA(tum.integrated.ns, npcs = 100)
  png(glue("{save_path}/dim_{dims}k_{k.param}_split_by_{var_split}_scaled_by_{var_scale}_elbow100_method_13.rds"), width=8, height=5, units = "in", res = 300)
  ep = ElbowPlot(tum.integrated.ns, ndims = 75)
  dev.off()
  message("finding neighbors....7/8")
  tum.integrated.ns = RunUMAP(object = tum.integrated.ns, reduction = "pca", dims = 1:dims)
  tum.integrated.ns = RunTSNE(object = tum.integrated.ns, reduction = "pca", nthreads=100, dims = 1:dims)
  tum.integrated.ns = FindNeighbors(object = tum.integrated.ns, dims = 1:dims, reduction = "pca", k.param = k.param)
  tum.integrated.ns = FindClusters(object = tum.integrated.ns, resolution = c(1, 0.8, 0.6,0.4,0.3,0.25, 0.2, 0.15, 0.1, 0.75, 0.05, 0.02), verbose = TRUE)
  write_rds(tum.integrated.ns, glue("{save_path}/dim_{dims}k_{k.param}_split_by_{var_split}_scaled_by_{var_scale}_tum.integrated.ns2_method_13.rds")) #tum_int = read_rds(glue("{save_path}/dim_{dims}k_{k.param}_tum.integrated.ns2_method_13.rds"))
  
  # ** PLOTS 
  tum.integrated.ns = SetIdent(tum.integrated.ns, value = "integrated_snn_res.1")
  message("saving plots....8/8")
  reduction = c("umap", "tsne")
  for (i in 1:2) {
    reduction = c("umap", "tsne")
    DefaultAssay(tum.integrated.ns) <- "RNA"
    reduction = reduction[i]
    message(reduction)
    p1 = DimPlot(object = tum.integrated.ns, reduction = reduction, label = TRUE, pt.size = 0.001, cols = cols, group.by = "integrated_snn_res.0.4")+  NoAxes()
    p2 = DimPlot(object = tum.integrated.ns, reduction = reduction, label = FALSE, pt.size = 0.001, cols = cols, group.by = "integrated_snn_res.0.6") +  NoAxes()
    p3 = DimPlot(object = tum.integrated.ns, reduction = reduction, group.by = "sheet_procedure",  pt.size = 0.00001, cols = c("tan1","steelblue1", "pink", "seagreen")) + NoAxes()
    p4 = DimPlot(object = tum.integrated.ns, reduction = reduction, group.by = "batch",  pt.size = 0.0001, cols = c("darkorange","darkorchid")) +  NoAxes()
    p5 = DimPlot(object = tum.integrated.ns, reduction = reduction, group.by = "sheet_patient_id", pt.size = 0.00013, cols = colors_dark) +  NoAxes()
    p6 = DimPlot(object = tum.integrated.ns, reduction = reduction, group.by = "sheet_exp_proc", pt.size = 0.00003, cols = colors_dark) +  NoAxes()
    p7 = DimPlot(object = tum.integrated.ns, reduction = reduction, group.by = "integrated_snn_res.0.4", split.by = "sheet_exp_proc", pt.size = 0.00003, cols = colors_dark) +  NoAxes()
    
    # system("pdftoppm input.pdf outputname -png")
    # save plots
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, "split_by_",var_split, "_scaled_by_", var_scale,"umap_all_labelled_method_13.png"), width=6, height=5, units = "in", res = 300)
    print(p1)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, "split_by_",var_split, "_scaled_by_", var_scale,"umap_all_method_13.png"), width=6, height=5, units = "in", res = 300)
    print(p2)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, "split_by_",var_split, "_scaled_by_", var_scale, "umap_procedure_method_13.png"), width=6, height=5, units = "in", res = 300)
    print(p3)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, "split_by_",var_split, "_scaled_by_", var_scale,"umap_batch_method_13.png"), width=6, height=5, units = "in", res = 300)
    print(p4)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, "split_by_",var_split, "_scaled_by_", var_scale, "umap_patient_method_13.png"), width=6, height=5, units = "in", res = 300)
    print(p5)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, "split_by_",var_split, "_scaled_by_", var_scale, "umap_exp_proc_method_13.png"), width=6, height=5, units = "in", res = 300)
    print(p6)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, "split_by_",var_split, "_scaled_by_", var_scale, "umap_exp_proc_SPLIT_method_13.png"), width=20, height=5, units = "in", res = 300)
    print(p7)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, "split_by_",var_split, "_scaled_by_", var_scale, "umap_markers1_method_13.png"), width=16, height=15, units = "in", res = 300)
    print(FeaturePlot(tum.integrated.ns, features = markers1, reduction = reduction, pt.size = 0.01) + NoAxes())
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, "split_by_",var_split, "_scaled_by_", var_scale,"umap_markers2_method_13.png"), width=16, height=15, units = "in", res = 300)
    print(FeaturePlot(tum.integrated.ns, features = markers2, reduction = reduction, pt.size = 0.01) + NoAxes())
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, "split_by_",var_split, "_scaled_by_", var_scale,"umap_markers3_method_13.png"), width=16, height=15, units = "in", res = 300)
    print(FeaturePlot(tum.integrated.ns, features = markers3, reduction = reduction, pt.size = 0.01) + NoAxes())
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, "split_by_",var_split, "_scaled_by_", var_scale, "umap_markers4_method_13.png"), width=16, height=15, units = "in", res = 300)
    print(FeaturePlot(tum.integrated.ns, features = markers4, reduction = reduction, pt.size = 0.01) + NoAxes())
    dev.off()
  }
  return(tum.integrated.ns)
}

# ***markers -------
func_markers_scale = function(clus_out = clus_out, var_split = var_split, var_scale = NULL, char_ident = char_ident, dims = dims, 
                              fts1 = fts1, k.param = k.param, cols = colors_canv, 
                              save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20190504_cellrangerv3_cell/04_subset_Tcells_v2/all"){
  
  tum.integrated.ns = SetIdent(clus_out, value = char_ident)
  tum.integrated.ns = ScaleData(tum.integrated.ns, assay = "RNA", vars.to.regress = NULL)# no scaling
  
  tum.integrated.ns = NormalizeData(tum.integrated.ns, assay = "RNA")
  var_scale = "NULL"
  
  char_assays = c("RNA", "SCT", "integrated")
  for (i in 1:length(char_assays)) {
    set_assay = char_assays[i]
    DefaultAssay(tum.integrated.ns) <- set_assay
    markers = FindAllMarkers(tum.integrated.ns, logfc.threshold = 0)
    write.csv(markers, glue("{save_path}/dims_{dims}k_{k.param}_tum.integrated.ns_{set_assay}_markers_intgeratedassay_{char_ident}_split_by_{var_split}_scaled_by_{var_scale}_method_13.csv"))
    write_rds(tum.integrated.ns, glue("{save_path}/dim_{dims}k_{k.param}_tum.integrated.ns_{set_assay}_v4_intgeratedassay_{char_ident}_split_by_{var_split}_scaled_by_{var_scale}_method_13.rds"))
    #markers = read.csv(glue("{save_path}/dims_17k_25_tum.integrated.ns_RNA_markers_intgeratedassay_integrated_snn_res.0.3_split_by_sheet_exp_proc_scaled_by_NULL_method_13.csv"))

    # heatmap
    # FindAllMarkers usually uses data slot in the RNA assay to find differential genes.
    # For a heatmap or dotplot of markers, the scale.data in the RNA assay should be used.
    #markers_v1 = markers$gene[!str_detect(markers$gene, "^MT-")]
    markers_v2 = markers %>% dplyr::filter(gene %in% markers$gene)
    
    #tum.integrated.ns = NormalizeData(tum.integrated.ns, assay = "RNA")
    #levels(markers_v2$cluster) = c("M1", "M2", "M3")
    #markers_v2 = markers_v2 %>% dplyr::mutate(cluster_1 = factor(x = cluster, levels = c("M1", "M2", "M3"), ordered = TRUE))
    top_markers = markers_v2 %>%
      #mutate(cluster_1 = factor(cluster, levels = c("M1", "M2", "M3"))) %>%
      dplyr::group_by(cluster) %>%
      top_n(20, avg_logFC) %>%
      dplyr::filter(avg_logFC > 0) %>%
      pull(gene) %>% as.character()
    
    message("Making Heatmaps")
    DefaultAssay(tum.integrated.ns) <- "RNA"
    png(paste0(save_path,"/heatmap_dim_",dims,"k_",k.param,"_tum.integrated.ns_assaymarkers", set_assay, "_RNAassay_", char_ident, "_split_by_", var_split, "_scaled_by_",var_scale,"_method_13.png"),  width=10,height=12, units = "in", res = 200)
    print(DoHeatmap(tum.integrated.ns, features = c(fts1, top_markers), size = 4,
                    group.by = char_ident, slot = "scale.data", disp.min = -1.5, disp.max = 1.5, group.colors = colors_dark) +
            scale_color_manual(values = colors_dark) +
            scale_fill_gradient2(low="steelblue1", high = "tomato3", mid = "white", midpoint = 0) +
            theme(axis.text.y.left = element_text(size=6)))
    dev.off()
    
    DefaultAssay(tum.integrated.ns) <- "SCT"
    png(paste0(save_path,"/heatmap_dim_",dims,"k_",k.param,"_tum.integrated.ns_assaymarkers", set_assay, "_SCTassay_", char_ident, "_split_by_", var_split, "_scaled_by_",var_scale,"_method_13.png"),  width=10,height=12, units = "in", res = 200)
    print(DoHeatmap(tum.integrated.ns, features = c(fts1, top_markers), size = 4,
                    group.by = char_ident, slot = "scale.data",disp.min = -1.5, disp.max = 1.5, group.colors = colors_dark) +
            scale_color_manual(values = colors_dark) +
            scale_fill_gradient2(low="steelblue1", high = "tomato3", mid = "white", midpoint = 0) +
            theme(axis.text.y.left = element_text(size=6)))

    dev.off()
    
    DefaultAssay(tum.integrated.ns) <- "integrated"
    png(paste0(save_path,"/heatmap_dim_",dims,"k_",k.param,"_tum.integrated.ns_assaymarkers", set_assay, "_intassay_", char_ident, "_split_by_", var_split, "_scaled_by_",var_scale,"_method_13.png"),  width=10,height=12, units = "in", res = 200)
    print(DoHeatmap(tum.integrated.ns, features = c(fts1, top_markers), size = 4,
                    group.by = char_ident, slot = "scale.data",disp.min = -1.5, disp.max = 1.5, group.colors = colors_dark) +
            scale_color_manual(values = colors_dark) +
            scale_fill_gradient2(low="steelblue1", high = "tomato3", mid = "white", midpoint = 0) +
            theme(axis.text.y.left = element_text(size=6)))
    
    dev.off()
    
    # cols.use <- list(sheet_patient_id = colors_dark[1:60])
    # names(cols.use[['sheet_patient_id']]) <- unique(tum@meta.data$sheet_patient_id)
    # DoMultiBarHeatmap(tum, features=top_markers, group.by=char_ident, additional.group.by = c('sheet_patient_id'), assay = "RNA", cols.use = cols.use)
    # DoMultiBarHeatmap(tum, features=top_markers, group.by='sheet_patient_id', additional.group.by = c(char_ident), assay = "RNA", cols.use = cols.use, draw.lines = F)
    # 
    # 
    # Idents(tum) <- "integrated_snn_res.0.3"
    # Idents(tum) <- "integrated_snn_res.0.3"
    # tum1 = subset(tum, idents = "0")
    # VlnPlot(tum1, "SOD2")
    # tum2<-WhichCells(object=tum1, expression = SOD2 > 1)
    # Idents(tum1) <- "sheet_patient_id"
    # tum2<-WhichCells(object=tum1,  idents = c("pt12", "pt20", "pt25", "pt26", "pt29", "pt30", "pt31", "pt34"))
    # 
    # tum$cellToDisplay <- row.names(tum@meta.data) %in% tum2
    # FeaturePlot(tum, features = c("cellToDisplay", "rna_SOD2"))
    
                                                           
    message("Making featureplots")
    reduction = c("umap", "tsne")
    for (i in 1:2) {
      reduction = c("umap", "tsne")
      DefaultAssay(tum.integrated.ns) <- set_assay
      reduction = reduction[i]
      message(reduction)
      png(paste0(save_path,"/", reduction, dims,"k_", k.param, "_tum.integrated.ns", set_assay, "_integratedassay_", char_ident,"vst_method_13.png"), width=6,height=5,units = "in", res = 300)
      print(DimPlot(object = tum.integrated.ns, reduction = reduction, label = TRUE, pt.size = 0.002, cols = cols, group.by = "integrated_snn_res.0.3", shuffle = T) +  NoAxes())
      print(DimPlot(object = tum.integrated.ns, reduction = reduction, label = FALSE, pt.size = 0.001, cols = cols, group.by = char_ident, shuffle = T) +  NoAxes())
      dev.off()
    }
    
  }
}


## ***Clustering <10> -----------------
# SCTransform, no integration

# initialize variables
func_cluster_sct = function(tum_int = seu_object, var_scale = "source", dims = dims, k.param = k.param, 
                            markers1 = fts1, markers2 = fts2, 
                            markers3 = fts3, markers4 = fts4,cols = colors_canv,
                            save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20190504_cellrangerv3_cell/04_subset_Tcells_v2/all"){
  
  seu_object <- SCTransform(object = seu_object, vars.to.regress = var_scale, verbose = TRUE)
  
  # ** SCALE 2
  tum.integrated.ns = seu_object
  #tum.integrated.ns = ScaleData(tum.integrated.ns, assay = "integrated") # no scaling
  tum.integrated.ns = RunPCA(tum.integrated.ns, npcs = 100)
  #tum.integrated.ns = RunUMAP(object = tum.integrated.ns, reduction = "pca", nthreads=100, dims = 1:dims,min.dist = 0.003, spread = 0.5)
  
  tum.integrated.ns = FindNeighbors(object = tum.integrated.ns, dims = 1:dims, reduction = "pca", k.param = k.param)
  tum.integrated.ns = FindClusters(object = tum.integrated.ns, resolution = c(1, 0.8, 0.6,0.4,0.3,0.2, 0.1, 0.05, 0.02), verbose = TRUE)
  tum.integrated.ns = RunUMAP(object = tum.integrated.ns, reduction = "pca", nthreads=100, dims = 1:dims)
  tum.integrated.ns = RunTSNE(object = tum.integrated.ns, reduction = "pca", nthreads=100, dims = 1:dims)
  write_rds(tum.integrated.ns, glue("{save_path}/dim_{dims}k_{k.param}_split_by_{var_scale}_tum.sct_SEL_v4_method_10.rds")) #tum_int = read_rds(glue("{save_path}/dim_{dims}k_{k.param}_tum.sct_SEL_v4.rds"))
  tum.integrated.ns = read_rds(glue("{save_path}/dim_{dims}k_{k.param}_split_by_{var_scale}_tum.sct_SEL_v4_method_10.rds"))
  
  # ** PLOTS 
  tum.integrated.ns = SetIdent(tum.integrated.ns, value = "SCT_snn_res.0.1")
  
  reduction = c("umap", "tsne")
  for (i in 1:2) {
    reduction = c("umap", "tsne")
    #DefaultAssay(tum.integrated.ns) <- "RNA"
    reduction = reduction[i]
    message(reduction)
    p1 = DimPlot(object = tum.integrated.ns, reduction = reduction, label = TRUE, pt.size = 0.25, cols = cols, group.by = "SCT_snn_res.0.1")+  NoAxes()
    p2 = DimPlot(object = tum.integrated.ns, reduction = reduction, label = FALSE, pt.size = 0.25, cols = cols, group.by = "SCT_snn_res.0.2")+  NoAxes()
    p3 = DimPlot(object = tum.integrated.ns, reduction = reduction, group.by = "procedure",  pt.size = 0.1, cols = c("cornsilk3","lightblue"))+  NoAxes()
    p4 = DimPlot(object = tum.integrated.ns, reduction = reduction, group.by = "source",  pt.size = 0.1, cols = c("pink","palegreen"))+ NoAxes()
    p5 = DimPlot(object = tum.integrated.ns, reduction = reduction, group.by = "patient_id", pt.size = 0.03,  cols = colors_dark)+  NoAxes()
    p6 = DimPlot(object = tum.integrated.ns, reduction = reduction, group.by = "exp_proc", pt.size = 0.03, cols = colors_dark)+  NoAxes()
    
    # system("pdftoppm input.pdf outputname -png")
    # save plots
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale, "umap_all_labelled_method_10.png"), width=6, height=5, units = "in", res = 300)
    print(p1)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale, "umap_all_method_10.png"), width=6, height=5, units = "in", res = 300)
    print(p2)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_procedure_method_10.png"), width=6, height=5, units = "in", res = 300)
    print(p3)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_source_method_10.png"), width=6, height=5, units = "in", res = 300)
    print(p4)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_patient_method_10.png"), width=6, height=5, units = "in", res = 300)
    print(p5)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_exp_proc_method_10.png"), width=6, height=5, units = "in", res = 300)
    print(p6)
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_markers1_method_10.png"), width=16, height=15, units = "in", res = 300)
    print(FeaturePlot(tum.integrated.ns, features = markers1, reduction = reduction, pt.size = 0.01) + NoAxes())
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_markers2_method_10.png"), width=16, height=15, units = "in", res = 300)
    print(FeaturePlot(tum.integrated.ns, features = markers2, reduction = reduction, pt.size = 0.01) + NoAxes())
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_markers3_method_10.png"), width=16, height=15, units = "in", res = 300)
    print(FeaturePlot(tum.integrated.ns, features = markers3, reduction = reduction, pt.size = 0.01) + NoAxes())
    dev.off()
    
    png(paste0(save_path,"/", reduction,dims,"k_", k.param, var_scale,"umap_markers4_method_10.png"), width=16, height=15, units = "in", res = 300)
    print(FeaturePlot(tum.integrated.ns, features = markers4, reduction = reduction, pt.size = 0.01) + NoAxes())
    dev.off()
    
  }
  return(tum.integrated.ns)
}

# ***markers -------
func_markers_sct = function(clus_out = clus_out, char_ident = char_ident, var_scale = var_scale, dims = dims, fts1 = fts1, k.param = k.param, cols = colors_canv, save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20190504_cellrangerv3_cell/04_subset_Tcells_v2/all"){
  
  tum.integrated.ns = SetIdent(clus_out, value = char_ident)
  #tum.integrated.ns = ScaleData(tum.integrated.ns, assay = "SCT")
  tum.integrated.ns = ScaleData(tum.integrated.ns, assay = "RNA")# no scaling
  
  tum.integrated.ns = NormalizeData(tum.integrated.ns, assay = "RNA")
  
  char_assays = c("RNA", "SCT")
  for (i in 1:length(char_assays)) {
    set_assay = char_assays[i]
    DefaultAssay(tum.integrated.ns) <- set_assay
    markers = FindAllMarkers(tum.integrated.ns)
    write.csv(markers, glue("{save_path}/dims_{dims}k_{k.param}_tum.sct_{set_assay}_splitby_{var_scale}_markers_intgeratedassay_{char_ident}_method_10.csv"))
    write_rds(tum.integrated.ns, glue("{save_path}/dim_{dims}k_{k.param}_tum.sct_{set_assay}_v4_splitby_{var_scale}_intgeratedassay_{char_ident}_method_10.rds"))
    
    # heatmap
    # FindAllMarkers usually uses data slot in the RNA assay to find differential genes.
    # For a heatmap or dotplot of markers, the scale.data in the RNA assay should be used.
    markers_v1 = markers$gene[!str_detect(markers$gene, "^MT-")]
    markers_v1 = markers_v1[!str_detect(markers_v1, "^RP")]
    markers_v2 = markers %>% dplyr::filter(gene %in% markers_v1)
    
    #tum.integrated.ns = NormalizeData(tum.integrated.ns, assay = "RNA")
    #levels(markers_v2$cluster) = c("M1", "M2", "M3")
    #markers_v2 = markers_v2 %>% dplyr::mutate(cluster_1 = factor(x = cluster, levels = c("M1", "M2", "M3"), ordered = TRUE))
    top_markers = markers_v2 %>%
      #mutate(cluster_1 = factor(cluster, levels = c("M1", "M2", "M3"))) %>%
      dplyr::group_by(cluster) %>%
      top_n(20, avg_logFC) %>%
      dplyr::filter(avg_logFC > 0) %>%
      pull(gene) %>% as.character()
    
    message("Making Heatmaps")
    DefaultAssay(tum.integrated.ns) <- "RNA"
    pdf(paste0(save_path,"/heatmap_dim_",dims,"k_",k.param,"_tum.sct_assaymarkers_split_by",var_scale,"_", set_assay, "_RNAassay_", char_ident,"_method_10.pdf"),width=10,height=8)
    print(DoHeatmap(tum.integrated.ns, features = c(fts1, top_markers), size = 4,
                    group.by = char_ident, disp.min = -1.5, disp.max = 1.5) +
            scale_fill_gradient2(low="steelblue1",high = "tomato3", mid = "white", midpoint = 0) + theme(axis.text.y.left = element_text(size=8)))
    dev.off()
    
    DefaultAssay(tum.integrated.ns) <- "SCT"
    pdf(paste0(save_path,"/heatmap_dim_",dims,"k_",k.param,"_tum.sct_assaymarkers_split_by",var_scale,"_", set_assay, "_SCTassay_", char_ident,"_method_10.pdf"),width=10,height=8)
    print(DoHeatmap(tum.integrated.ns, features = c(fts1, top_markers), size = 4,
                    group.by = char_ident, slot = "scale.data",disp.min = -1.5, disp.max = 1.5) +
            scale_fill_gradient2(low="steelblue1",high = "tomato3", mid = "white", midpoint = 0) + theme(axis.text.y.left = element_text(size=8)))
    dev.off()
    
    
    message("Making featureplots")
    reduction = c("umap", "tsne")
    for (i in 1:2) {
      reduction = c("umap", "tsne")
      DefaultAssay(tum.integrated.ns) <- set_assay
      reduction = reduction[i]
      message(reduction)
      pdf(paste0(save_path,"/", reduction, dims,"k_", k.param, "_tum.sct_split_by",var_scale,"_", set_assay, "_sct_", char_ident,"_method_10.pdf"), width=6,height=5)
      print(DimPlot(object = tum.integrated.ns, reduction = reduction, label = TRUE, pt.size = 0.25, cols = cols,  group.by = char_ident) +  NoAxes())
      print(DimPlot(object = tum.integrated.ns, reduction = reduction, label = FALSE, pt.size = 0.25, cols = cols, group.by = char_ident) +  NoAxes())
      dev.off()
    }
    
  }
}
## ***Barplots -----------------
# A function to create bar plots
func_bar = function(var_df, var_x, var_fill, var_pos, var_cols, celltype) {
  # create barplot
  p1 = ggplot(var_df) +
    geom_bar(aes_string(x = var_x, fill = var_fill), position = var_pos) + 
    xlab("") + ylab("") + scale_fill_manual(values = var_cols) +
    theme(axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
          panel.background = element_blank(),axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(2,0,0,0, "pt")), axis.ticks = element_blank(), 
          legend.title = element_blank())
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
  #         legend.text = element_text(),
  #         panel.background = element_blank())
  
  # save the plot
  ggsave(p1, file = paste0(save_path, "/p_bar", celltype,"_", var_pos, "_", var_fill, "_", var_x, ".pdf"), width = 4, height = 3) #12, 6.5
  ggsave(p1, file = paste0(save_path, "/p_bar", celltype,"_", var_pos, "_", var_fill, "_", var_x, ".png"), width = 4, height = 3) #12, 6.5
  return(p1)
}

# example run
#func_bar(var_df = df_meta, var_x = "trt_type", var_fill = "t_cells_v2", var_pos = "stack", var_cols = colors_darj, celltype = "T_cells")
#func_bar(var_df = df_meta, var_x = "trt_type", var_fill = "t_cells_v2", var_pos = "fill", var_cols = colors_darj, celltype = "T_cells")
# df_meta %>% 
#   #filter(new.clonotype_id %in% bottom80_clonos) %>% 
#   group_by(sample_id, integrated_snn_res_0_2) %>%
#   dplyr::summarise(ncells = n()) %>% 
#   mutate(freq = ncells/sum(ncells)) %>% 
#   ggplot() + aes(x=sample_id, y=freq, fill=integrated_snn_res_0_2) + geom_bar(stat="identity") +
#   scale_fill_manual(values = colors_rcart) +
#   scale_y_continuous(expand = c(0,0)) +
#   theme(axis.ticks.x = element_blank(), axis.title.y = element_blank(), 
#         panel.background = element_blank(),axis.title.x = element_blank(), axis.text.y = element_blank(),
#         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(2,0,0,0, "pt")), axis.ticks = element_blank(), 
#         legend.title = element_blank())


## ***Feature plots ------------
# functions
get_features = function(seu_obj, genes) {
  mat = seu_obj@assays$RNA@data
  mat = as.matrix(mat)
  mat[genes, ]
}

# get umap object
get_umap.seu <- function(x, genes = NULL){
  df_umap = x@reductions$umap@cell.embeddings
  df_umap = cbind(df_umap, x@meta.data[rownames(df_umap), ])
  
  # if genes supplied, extract them
  if(!is.null(genes)){
    mat_feat = get_features(x, genes)
    mat_feat = as.matrix(mat_feat)
    #mat_feat = t(mat_feat)
    cellnms = rownames(df_umap)
    df_umap_ann = cbind(df_umap[cellnms, ], mat_feat[cellnms,])
    return(df_umap_ann)
  }
  
  df_umap
}

# get tsne object
get_tsne.seu <- function(x, genes = NULL){
  df_tsne = x@reductions$tsne@cell.embeddings
  df_tsne = cbind(df_tsne, x@meta.data[rownames(df_tsne), ])
  
  # if genes supplied, extract them
  if(!is.null(genes)){
    mat_feat = get_features(x, genes)
    mat_feat = as.matrix(mat_feat)
    #mat_feat = t(mat_feat)
    cellnms = rownames(df_tsne)
    df_tsne_ann = cbind(df_tsne[cellnms, ], mat_feat[cellnms,])
    return(df_tsne_ann)
  }
  
  df_tsne
} 

# ** ex plot 
gg_feature_plot <- function(x, genes, reduction = "umap", save_folder, celltype = celltype){
  
  if(reduction %in% c("umap")){
    dim1 = "UMAP_1"
    dim2 = "UMAP_2"
    df_dim_ann = get_umap.seu(x, genes)
  } else{
    dim1 = "tSNE_1"
    dim2 = "tSNE_2"
    df_dim_ann = get_tsne.seu(x, genes)
  }
  
  
  
  # fix range of genes
  df_dim_ann[, genes] %>% head(100) %>% range()
  
  
  # plot every gene
  # gene = "ESR1"
  tmp = lapply(genes, function(gene){
    message(gene)
    # sort df, by gene value
    df_dim_ann = dplyr::arrange_at(df_dim_ann, .vars = gene)
    vals = df_dim_ann[, gene]
    df_dim_ann[, gene] = ifelse(vals == 0, NA, vals)
    
    tail(df_dim_ann, 30)
    p = ggplot(df_dim_ann, aes_string(dim1, dim2, color = gene)) +
      geom_point(alpha = 0.9, size = 0.01) +
      #scale_fill_material("red")
      # scale_colour_gradient(low = "grey", high = "darkblue")
      #scale_fill_gradientn(colours = (10))
      scale_color_gradient2(low = "white", mid = "bisque2", high = "red", na.value = "grey90") +
      theme_cowplot()
    p
    #save_plot(filename = paste0(wd, save_folder, "/gene_plots_clus/", celltype, "_", gene, "_", reduction, ".pdf"), plot = p, base_height = 5, base_width = 6)
    #save_plot(filename = paste0(wd, save_folder, "/all/gene_plots/", celltype, "_", gene, "_", reduction, ".png"), plot = p, base_height = 5, base_width = 6)
    
    #save_plot(filename = paste0(wd, save_folder, "/gene_plots_clus_all/", gene, "_", reduction, ".pdf"), plot = p, base_height = 5, base_width = 6)
    save_plot(filename = paste0(wd, save_folder, "/gene_plots_all/", gene, "_", reduction, ".png"), plot = p, base_height = 5, base_width = 6)
  })
}
#df_umap = get_tsne.seu(tum)

# example run
#gg_feature_plot(tum, genes = genes, reduction = "umap", save_folder = "02_gene_markers")

# ***DoHeatmap -------------

run_htmaps_state <- function(df, seu, cell, genes_df) {
  # de scores ------
  #mat_de = df %>% dplyr::filter(abs(avg_logFC)>=1, p_val_adj <= 0.05)
  #genes_df = genes_df %>% dplyr::mutate(gene = sapply(strsplit(gene, "[.]"), "[", 1))
  #df = df %>% dplyr::mutate(gene = sapply(strsplit(gene, "[.]"), "[", 1))
  mat_de = df %>% dplyr::filter(gene %in% genes_df)
  #mat_de = mat_de %>% tbl_df() %>% dplyr::mutate(gene = sapply(strsplit(gene, "[.]"), "[", 1))
  mat_de = mat_de %>% dplyr::select(gene, avg_logFC, cluster)
  mat_de_wd = reshape2::dcast(mat_de, gene ~ cluster, value.var = "avg_logFC", fill = 0)
  mat_de_fil = data.matrix(mat_de_wd[,-1])
  
  rownames(mat_de_fil) = mat_de_wd$gene
  states = as.character(unique(seu$minor_groups))
  mat_de_fil = mat_de_fil[genes_df,]
  mat_de_fil = mat_de_fil[,states]
  message("de scores done")
  
  # heatmap of de
  ht_de = Heatmap(data.matrix(mat_de_fil),
                  column_title = "DE scores", column_title_gp = gpar(fontsize = 8),
                  col = colorRamp2(c(-1.5,0, 1.5), c("purple","black", "yellow")),
                  name = "Avg Log FC",
                  #top_annotation = ha_top,
                  show_column_names = T, show_row_dend = F, cluster_rows = F,
                  use_raster = T, cluster_columns = FALSE,show_row_names = F, 
                  row_names_gp = gpar(fontsize = 8),  width = unit(2, "cm"))
  message("htmap of de done")
  #ht_de
  
  # exp scores -------
  # get the normalized data matrix and extract subset of de genes
  # normalized and log-transformed single cell expression
  # Setting center to TRUE will center the expression for each gene by subtracting the 
  # average expression for that gene. Setting scale to TRUE will scale the expression level 
  # for each gene by dividing the centered gene expression levels by their standard deviations 
  # if center is TRUE and by their root mean square otherwise. 
  mat_exp = seu@assays$integrated@scale.data
  #mat_exp = seu@assays$SCT_vasc@scale.data
  mat_exp = data.matrix(mat_exp)
  #dim(mat_exp)
  message("exp scores done")
  
  # filter for genes that are de from wilcox test, abs(avg FC) >=1, adj p val <=0.05
  d1 = mat_exp %>% tbl_df() %>% dplyr::mutate(genes = rownames(mat_exp))
  #d1 = mat_exp %>% tbl_df() %>% dplyr::mutate(genes = sapply(strsplit(rownames(mat_exp), "[.]"), "[", 1))
  d2 = d1 %>% dplyr::filter(genes %in% mat_de_wd$gene)
  mat_exp_fil = data.matrix(d2)
  mat_exp_fil = suppressWarnings(mat_exp_fil[,-ncol(mat_exp_fil)])
  rownames(mat_exp_fil) = d2$genes
  mat_exp_fil = mat_exp_fil[rownames(mat_de_fil),]
  table(rownames(mat_de_fil) == rownames(mat_exp_fil)) # check the order of genes matches to that in de list from wilcox
  message("filter done")
  
  # Annotations -------
  # get df for annotation
  coldata = data.frame(Cluster = seu@meta.data$minor_groups, 
                       Procedure = seu@meta.data$procedure,
                       Source = seu@meta.data$source,
                       Protocol = seu@meta.data$exp_proc,
                       Patient = seu@meta.data$patient_id,
                       #meno = seu$,
                       #brca = seu$brca_status,
                       nFeature_RNA = seu$nFeature_RNA,
                       nCount_RNA = seu$nCount_RNA,
                       sample_id = colnames(mat_exp_fil),
                       stringsAsFactors = F) %>% dplyr::arrange(Cluster)
  
  # top annotation
  ha_top = dplyr::select(coldata, Cluster, Procedure, Source, Protocol) %>% 
    HeatmapAnnotation(df = ., show_annotation_name = TRUE, 
                      col = list(Cluster = c("L-1" = endo_colors[1], "L-2" = endo_colors[2], "V-1" = endo_colors[3],"V-2" = endo_colors[4],"V-3" = endo_colors[5]),
                                 Procedure = c("reduction" = "cyan", "mastectomy" = "lightpink"),
                                 Source = c("uci" = "red", "mda" = "blue"),
                                 Protocol = c("short_digestion" = colors_darjee[2], 
                                              "6hr_digestion" = colors_darjee[3], "overnight_digestion" = colors_darjee[4])))
  
  # get a sorted list of samples
  samps_sorted = coldata$sample_id
  
  # rearrnge the matrix
  mat = mat_exp_fil[, samps_sorted]
  
  # heatmap of exp
  ht_exp = Heatmap(data.matrix(mat),
                   column_title = "Scaled Expression",
                   column_title_gp = gpar(fontsize = 6),
                   name = "Scaled Expression",
                   top_annotation = ha_top,
                   show_row_names = T,
                   show_column_names = F, cluster_rows = F,
                   use_raster = T, cluster_columns = FALSE,
                   col = colorRamp2(c(-2, 0, 2), c("black","black", "cornsilk")),
                   row_names_gp = gpar(fontsize = 8),  width = unit(4, "cm"))
  ht_exp_notrow = Heatmap(data.matrix(mat),
                          column_title = "Scaled Expression",
                          column_title_gp = gpar(fontsize = 6),
                          name = "Scaled Expression",
                          top_annotation = ha_top,
                          col = colorRamp2(c(-2, 0, 2), c("steelblue2","white", "deeppink4")),
                          show_column_names = F, cluster_rows = F,
                          use_raster = T, cluster_columns = FALSE,
                          row_names_gp = gpar(fontsize = 8),  width = unit(12, "cm"))
  #ht_exp
  message("ht exp done")
  
  # Mean ht ------
  # get the mean of exp 
  scores = cbind(coldata, t(mat))#mat_exp_fil
  scores = scores[,colnames(unique(as.matrix(scores),MARGIN=2))] 
  scores_m = suppressWarnings(scores %>% group_by(Cluster) %>% unique() %>% summarise_all(funs(mean)))
  mat_exp_mean = data.matrix(t(scores_m))
  mat_exp_mean = mat_exp_mean[rownames(mat_de_fil),]
  #mat_exp_mean = mat_exp_mean[-c(1,2), ]
  colnames(mat_exp_mean) = levels(scores$Cluster)
  mat_exp_mean = apply(mat_exp_mean, 2,as.numeric)
  rownames(mat_exp_mean) = rownames(mat_de_fil)
  
  # heatmap of exp
  ha_top_mn =  data.frame(Cluster = levels(scores$Cluster)) %>% HeatmapAnnotation(df = ., show_annotation_name = TRUE,
                                                                                  col = list(Cluster = c("Secretory" = colors_darj[1], "Capillary" = colors_darj[2], "Venous" = colors_darj[3])))
  #%>% HeatmapAnnotation(col = list(Condition = c("0" = "orange", "1" = "cyan", "2" = "grey")))
  ht_mean = Heatmap(data.matrix(mat_exp_mean),
                    column_title = "Avg Exp",column_title_gp = gpar(fontsize = 6),
                    name = "Avg Expression",
                    top_annotation = ha_top_mn,
                    show_row_names = T, row_names_side = c("left"), 
                    show_column_names = F, show_row_dend = F, cluster_rows = F,
                    use_raster = T, cluster_columns = FALSE,
                    row_names_gp = gpar(fontsize = 8),  width = unit(1.5, "cm"))
  #ht_mean
  message("ht means done")
  
  pdf(glue("{odir}/htmaps_de_{cell}_mean_exp_descores.pdf"), width = 12, height = 10)
  print(ht_de + ht_exp_notrow)
  dev.off()
  
  pdf(glue("{odir}/htmaps_de_{cell}_mean_exp_descores_rowclust.pdf"), width = 16, height = 10)
  print(ht_mean + ht_exp + ht_de)
  dev.off()
  
  
}

# ***fig_heatmap -------------

run_htmaps_fig1 <- function(df, seu, cell, genes_df) {
  # de scores ------
  #mat_de = df %>% dplyr::filter(abs(avg_logFC)>=1, p_val_adj <= 0.05)
  #genes_df = genes_df %>% dplyr::mutate(gene = sapply(strsplit(gene, "[.]"), "[", 1))
  #df = df %>% dplyr::mutate(gene = sapply(strsplit(gene, "[.]"), "[", 1))
  mat_de = df %>% dplyr::filter(gene %in% genes_df)
  #mat_de = mat_de %>% tbl_df() %>% dplyr::mutate(gene = sapply(strsplit(gene, "[.]"), "[", 1))
  mat_de = mat_de %>% dplyr::select(gene, avg_logFC, cluster)
  mat_de_wd = reshape2::dcast(mat_de, gene ~ cluster, value.var = "avg_logFC", fill = 0)
  mat_de_fil = data.matrix(mat_de_wd[,-1])
  
  rownames(mat_de_fil) = mat_de_wd$gene
  states = as.character(unique(seu$minor_groups))
  mat_de_fil = mat_de_fil[genes_df,]
  mat_de_fil = mat_de_fil[,states]
  message("de scores done")
  
  # heatmap of de
  ht_de = Heatmap(data.matrix(mat_de_fil),
                  column_title = "DE scores", column_title_gp = gpar(fontsize = 8),
                  col = colorRamp2(c(-1.5,0, 1.5), c("purple","black", "yellow")),
                  name = "Avg Log FC",
                  #top_annotation = ha_top,
                  show_column_names = T, show_row_dend = F, cluster_rows = F,
                  use_raster = T, cluster_columns = FALSE,show_row_names = F, 
                  row_names_gp = gpar(fontsize = 8),  width = unit(2, "cm"))
  message("htmap of de done")
  #ht_de
  
  # exp scores -------
  # get the normalized data matrix and extract subset of de genes
  # normalized and log-transformed single cell expression
  # Setting center to TRUE will center the expression for each gene by subtracting the 
  # average expression for that gene. Setting scale to TRUE will scale the expression level 
  # for each gene by dividing the centered gene expression levels by their standard deviations 
  # if center is TRUE and by their root mean square otherwise. 
  mat_exp = seu@assays$RNA@scale.data
  #mat_exp = seu@assays$SCT_vasc@scale.data
  mat_exp = data.matrix(mat_exp)
  #dim(mat_exp)
  message("exp scores done")
  
  # filter for genes that are de from wilcox test, abs(avg FC) >=1, adj p val <=0.05
  d1 = mat_exp %>% tbl_df() %>% dplyr::mutate(genes = rownames(mat_exp))
  #d1 = mat_exp %>% tbl_df() %>% dplyr::mutate(genes = sapply(strsplit(rownames(mat_exp), "[.]"), "[", 1))
  d2 = d1 %>% dplyr::filter(genes %in% mat_de_wd$gene)
  mat_exp_fil = data.matrix(d2)
  mat_exp_fil = suppressWarnings(mat_exp_fil[,-ncol(mat_exp_fil)])
  rownames(mat_exp_fil) = d2$genes
  mat_exp_fil = mat_exp_fil[rownames(mat_de_fil),]
  table(rownames(mat_de_fil) == rownames(mat_exp_fil)) # check the order of genes matches to that in de list from wilcox
  message("filter done")
  
  # Annotations -------
  # get df for annotation
  coldata = data.frame(Cluster = seu@meta.data$minor_groups, 
                       Procedure = seu@meta.data$procedure,
                       Source = seu@meta.data$source,
                       Protocol = seu@meta.data$exp_proc,
                       Patient = seu@meta.data$patient_id,
                       #meno = seu$,
                       #brca = seu$brca_status,
                       nFeature_RNA = seu$nFeature_RNA,
                       nCount_RNA = seu$nCount_RNA,
                       sample_id = colnames(mat_exp_fil),
                       stringsAsFactors = F) %>% dplyr::arrange(Cluster)
  
  # top annotation
  ha_top = dplyr::select(coldata, Cluster, Procedure, Source, Protocol) %>% 
    HeatmapAnnotation(df = ., show_annotation_name = TRUE, 
                      col = list(Cluster = c("L-1" = endo_colors[1], "L-2" = endo_colors[2], "V-1" = endo_colors[3],"V-2" = endo_colors[4],"V-3" = endo_colors[5]),
                                 Procedure = c("reduction" = "cyan", "mastectomy" = "lightpink"),
                                 Source = c("uci" = "red", "mda" = "blue"),
                                 Protocol = c("short_digestion" = colors_darjee[2], 
                                              "6hr_digestion" = colors_darjee[3], "overnight_digestion" = colors_darjee[4])))
  
  # get a sorted list of samples
  samps_sorted = coldata$sample_id
  
  # rearrnge the matrix
  mat = mat_exp_fil[, samps_sorted]
  
  # heatmap of exp
  ht_exp = Heatmap(data.matrix(mat),
                   column_title = "Scaled Expression",
                   column_title_gp = gpar(fontsize = 6),
                   name = "Scaled Expression",
                   top_annotation = ha_top,
                   show_row_names = T,
                   show_column_names = F, cluster_rows = F,
                   use_raster = T, cluster_columns = FALSE,
                   col = colorRamp2(c(-2, 0, 2), c("black","black", "cornsilk")),
                   row_names_gp = gpar(fontsize = 8),  width = unit(4, "cm"))
  ht_exp_notrow = Heatmap(data.matrix(mat),
                          column_title = "Scaled Expression",
                          column_title_gp = gpar(fontsize = 6),
                          name = "Scaled Expression",
                          top_annotation = ha_top,
                          col = colorRamp2(c(-2, 0, 2), c("steelblue2","white", "deeppink4")),
                          show_column_names = F, cluster_rows = F,
                          use_raster = T, cluster_columns = FALSE,
                          row_names_gp = gpar(fontsize = 8),  width = unit(12, "cm"))
  #ht_exp
  message("ht exp done")
  
  # Mean ht ------
  # get the mean of exp 
  scores = cbind(coldata, t(mat))#mat_exp_fil
  scores = scores[,colnames(unique(as.matrix(scores),MARGIN=2))] 
  scores_m = suppressWarnings(scores %>% group_by(Cluster) %>% unique() %>% summarise_all(funs(mean)))
  mat_exp_mean = data.matrix(t(scores_m))
  mat_exp_mean = mat_exp_mean[rownames(mat_de_fil),]
  #mat_exp_mean = mat_exp_mean[-c(1,2), ]
  colnames(mat_exp_mean) = levels(scores$Cluster)
  mat_exp_mean = apply(mat_exp_mean, 2,as.numeric)
  rownames(mat_exp_mean) = rownames(mat_de_fil)
  
  # heatmap of exp
  ha_top_mn =  data.frame(Cluster = levels(scores$Cluster)) %>% HeatmapAnnotation(df = ., show_annotation_name = TRUE,
                                                                                  col = list(Cluster = c("Secretory" = colors_darj[1], "Capillary" = colors_darj[2], "Venous" = colors_darj[3])))
  #%>% HeatmapAnnotation(col = list(Condition = c("0" = "orange", "1" = "cyan", "2" = "grey")))
  ht_mean = Heatmap(data.matrix(mat_exp_mean),
                    column_title = "Avg Exp",column_title_gp = gpar(fontsize = 6),
                    name = "Avg Expression",
                    top_annotation = ha_top_mn,
                    show_row_names = T, row_names_side = c("left"), 
                    show_column_names = F, show_row_dend = F, cluster_rows = F,
                    use_raster = T, cluster_columns = FALSE,
                    row_names_gp = gpar(fontsize = 8),  width = unit(1.5, "cm"))
  #ht_mean
  message("ht means done")
  
  pdf(glue("{odir}/htmaps_de_{cell}_mean_exp_descores.pdf"), width = 12, height = 10)
  print(ht_de + ht_exp_notrow)
  dev.off()
  
  pdf(glue("{odir}/htmaps_de_{cell}_mean_exp_descores_rowclust.pdf"), width = 16, height = 10)
  print(ht_mean + ht_exp + ht_de)
  dev.off()
  
  
}

# heatmaps cell type ---------------
run_htmaps <- function(df, seu, cell, genes_df, clus_levels) {
  # de scores ------
  mat_de = df %>% dplyr::filter(gene %in% genes_df)
  mat_de = mat_de %>% dplyr::select(gene, avg_logFC, cluster)
  mat_de_wd = reshape2::dcast(mat_de, gene ~ cluster, value.var = "avg_logFC" )
  mat_de_fil = data.matrix(mat_de_wd[,-1])
  
  rownames(mat_de_fil) = mat_de_wd$gene
  mat_de_fil = mat_de_fil[genes_df,]
  mat_de_fil = mat_de_fil[,clus_levels]
  message("de scores done")
# 
#   
#   ha_top_mn =  data.frame(Cluster = Groups_set) %>% HeatmapAnnotation(df = ., show_annotation_name = FALSE,
#                                                                               col = list(Cluster = col))
#   
#   # heatmap of de
#   ht_de = Heatmap(data.matrix(mat_de_fil),
#                   column_title = "DE scores", column_title_gp = gpar(fontsize = 8),
#                   col = colorRamp2(c(-1.5,0, 1.5), c("steelblue","black", "tomato3")),
#                   name = "Avg Log FC",
#                   #top_annotation = ha_top_mn,
#                   show_column_names = F, show_row_dend = F, cluster_rows = F,
#                   use_raster = T, cluster_columns = FALSE,show_row_names = F, 
#                   row_names_gp = gpar(fontsize = 8),  width = unit(1, "cm"))
  #message("htmap of de done")
  #ht_de
  
  # exp scores -------
  mat_exp = seu@assays$RNA@scale.data
  mat_exp = data.matrix(mat_exp)
  #dim(mat_exp)
  message("exp scores done")
  
  # filter for genes that are de from wilcox test, abs(avg FC) >=1, adj p val <=0.05
  d1 = mat_exp %>% tbl_df() %>% dplyr::mutate(genes = rownames(mat_exp))
  #d1 = mat_exp %>% tbl_df() %>% dplyr::mutate(genes = sapply(strsplit(rownames(mat_exp), "[.]"), "[", 1))
  d2 = d1 %>% dplyr::filter(genes %in% mat_de_wd$gene)
  mat_exp_fil = data.matrix(d2)
  mat_exp_fil = suppressWarnings(mat_exp_fil[,-ncol(mat_exp_fil)])
  rownames(mat_exp_fil) = d2$genes
  mat_exp_fil = mat_exp_fil[rownames(mat_de_fil),]
  table(rownames(mat_de_fil) == rownames(mat_exp_fil)) # check the order of genes matches to that in de list from wilcox
  message("filter done")
  
  # Annotations -------
  # get df for annotation
  coldata = data.frame(Cluster = seu$major_group_factor, 
                       #Source = seu@meta.data$source,
                       Protocol = seu@meta.data$exp_proc,
                       #Patient = seu@meta.data$patient_id.x,
                       Menopause = seu@meta.data$menopause.y,
                       Paridy = seu@meta.data$parity.1.y,
                       Procedure = seu@meta.data$figure1_grp_updated,
                       #cancer_risk = seu@meta.data$cancer_risk,
                       #nFeature_RNA = seu$nFeature_RNA,
                       #nCount_RNA = seu$nCount_RNA,
                       
                       sample_id = colnames(mat_exp_fil),
                       stringsAsFactors = F) %>% dplyr::arrange(Cluster)
  
  # top annotation
  ha_top = dplyr::select(coldata, Cluster, Protocol, Menopause, Paridy, Procedure) %>% 
    HeatmapAnnotation(df = ., show_annotation_name = TRUE, 
                      col = list(Cluster = col,
                                 Procedure = col_procedure,
                                 Menopause = col_meno,
                                 Paridy = col_paridy,
                                 Protocol = col_Protocol))
  
  # get a sorted list of samples
  samps_sorted = coldata$sample_id
  
  # rearrnge the matrix
  mat = mat_exp_fil[, samps_sorted]
  
  # heatmap of exp
  ht_exp = Heatmap(data.matrix(mat),
                   column_title = "Scaled Expression",
                   column_title_gp = gpar(fontsize = 6),
                   name = "Scaled Expression",
                   top_annotation = ha_top,show_row_names = T,
                   show_column_names = F, cluster_rows = F,
                   use_raster = F, cluster_columns = FALSE,
                   col = colorRamp2(c(-1.5, 0, 1.5), c("steelblue","white", "tomato3")),
                   row_names_gp = gpar(fontsize = 8),  width = unit(8, "cm"))
 
  png(glue("{save_path}/htmaps_de_{cell}_mean_exp_descores_rowclust.png"), res = 300, units = "in", width = 13, height = 11)
  print(ht_exp)
  dev.off()
  
}

# example run
## Read Files
# df_fibro = read.csv(file = glue("/volumes/lab/users/tapsi/projects/bcap/analysis/20190504_cellrangerv3_cell/03_subset_fibroblasts/all/dims_17k_20_tum.integrated.ns_markers_RNA_filtered.csv"), stringsAsFactors = FALSE)
# df_fibro <- df_fibro[order(df_fibro$avg_logFC, decreasing=T),]
# ## get genes from all
# genes_fibro <- df_fibro %>% dplyr::group_by(cluster) %>% top_n(15, avg_logFC) %>% dplyr::select(gene)

# heatmaps mean ---------------
run_mean_htmaps <- function(df, seu, cell, genes_df, clus_levels) {
  # # de scores ------
  # mat_de = df %>% dplyr::filter(gene %in% genes_df)
  # mat_de = mat_de %>% dplyr::select(gene, avg_logFC, cluster)
  # mat_de_wd = reshape2::dcast(mat_de, gene ~ cluster, value.var = "avg_logFC" )
  # mat_de_fil = data.matrix(mat_de_wd[,-1])
  # 
  # rownames(mat_de_fil) = mat_de_wd$gene
  # mat_de_fil = mat_de_fil[genes_df,]
  # mat_de_fil = mat_de_fil[,clus_levels]
  # message("de scores done")
  # 
  
  # exp scores -------
  mat_exp = seu@assays$integrated@scale.data
  mat_exp = data.matrix(mat_exp)
  message("exp scores done")
  
  # filter for genes that are de from wilcox test, abs(avg FC) >=1, adj p val <=0.05
  d1 = mat_exp %>% tbl_df() %>% dplyr::mutate(genes = rownames(mat_exp))
  #d1 = mat_exp %>% tbl_df() %>% dplyr::mutate(genes = sapply(strsplit(rownames(mat_exp), "[.]"), "[", 1))
  d2 = d1 %>% dplyr::filter(genes %in% genes_df)
  mat_exp_fil = data.matrix(d2)
  mat_exp_fil = suppressWarnings(mat_exp_fil[,-ncol(mat_exp_fil)])
  rownames(mat_exp_fil) = d2$genes
  mat_exp_fil = mat_exp_fil[genes_df,] #match the order of orginal genes
  #table(rownames(mat_de_fil) == rownames(mat_exp_fil)) # check the order of genes matches to that in de list from wilcox
  message("filter done")
  
  # Annotations -------
  # get df for annotation
  coldata = data.frame(Cluster = seu$minor_group_factor, 
                        Protocol = seu@meta.data$exp_proc,
                        Menopause = seu@meta.data$menopause.y,
                        Paridy = seu@meta.data$parity.1.y,
                        Procedure = seu@meta.data$figure1_grp_updated,
                        sample_id = colnames(mat_exp_fil),
                        stringsAsFactors = F) %>% dplyr::arrange(Cluster)
   
  # top annotation
  ha_top = dplyr::select(coldata, Cluster, Protocol, Menopause, Paridy, Procedure) %>% 
    HeatmapAnnotation(df = ., show_annotation_name = TRUE, 
                      col = list(Cluster = col_grp,
                                 Procedure = col_procedure,
                                 Menopause = col_meno,
                                 Paridy = col_paridy,
                                 Protocol = col_Protocol))
  
  # get a sorted list of samples
  samps_sorted = coldata$sample_id
  
  # rearrnge the matrix
  mat = mat_exp_fil[, samps_sorted]
  
  # # heatmap of exp
  # ht_exp = Heatmap(data.matrix(mat),
  #             column_title = "Scaled Expression",
  #             column_title_gp = gpar(fontsize = 6),
  #             name = "Scaled Expression",
  #             top_annotation = ha_top,show_row_names = T,
  #             show_column_names = F, cluster_rows = F,
  #             use_raster = F, cluster_columns = FALSE,
  #             col = colorRamp2(c(-1.5, 0, 1.5), c("steelblue","white", "tomato3")),
  #             row_names_gp = gpar(fontsize = 8),  width = unit(8, "cm"))
              
  # Mean ht ------
  # get the mean of exp 
  scores = cbind(coldata, t(mat))#ma
  scores = scores[,colnames(unique(as.matrix(scores),MARGIN=2))] 
  scores_m = suppressWarnings(scores %>% group_by(Cluster) %>% unique() %>% summarise_all(funs(mean)))
  mat_exp_mean = data.matrix(t(scores_m))
  mat_exp_mean = mat_exp_mean[genes_df,]
  #mat_exp_mean = mat_exp_mean[-c(1,2), ]
  colnames(mat_exp_mean) = levels(scores$Cluster)
  mat_exp_mean = apply(mat_exp_mean, 2,as.numeric)
  rownames(mat_exp_mean) = genes_df
  
  # heatmap of exp
  ha = HeatmapAnnotation(Protocol = anno_barplot(t(df_meta_dig2), border = FALSE, gp = gpar(fill = col_Protocol, border = NA, lty = "blank"),axis = FALSE,bar_width = 0.95),
                         Procedure = anno_barplot(t(df_meta_grp2), border = FALSE, gp = gpar(fill = col_procedure, border = NA, lty = "blank"),axis = FALSE,bar_width = 0.95), 
                         Menopause = anno_barplot(t(df_meta_meno2), border = FALSE, gp = gpar(fill = col_meno, border = NA, lty = "blank"), axis = FALSE, bar_width = 0.95),#size = unit(0.85, "mm"), 
                         Paridy = anno_barplot(t(df_meta_paridy2), border = FALSE, gp = gpar(fill = col_paridy, border = NA, lty = "blank"), axis = FALSE,bar_width = 0.95),
                         Cluster = levels(scores$Cluster), col = list(Cluster = col_grp), show_annotation_name = TRUE)
  ht_mean = Heatmap(data.matrix(mat_exp_mean),
                    
                    #column_title = "Avg Exp",column_title_gp = gpar(fontsize = 6),
                    name = "Avg Expression",
                    top_annotation = ha,
                    show_row_names = T, row_names_side = c("left"), 
                    show_column_names = T,column_names_side = c("bottom"), 
                    show_row_dend = F, cluster_rows = F,
                    use_raster = T, cluster_columns = FALSE,
                    col = colorRamp2(c(-1, 0, 2), c("steelblue","white", "tomato3")),
                    #col = colorRamp2(c(-1, 0, 2), c("cadetblue4","white", "coral3")),
                    #col = colorRamp2(c(-1.5, 0, 1.5), c("black","black", "grey99")),
                    row_names_gp = gpar(fontsize = 7),  width = unit(4.5, "cm"))
  ht1 = draw(ht_mean, show_heatmap_legend = FALSE, show_annotation_legend = FALSE)
  ht2 = draw(ht_mean, show_heatmap_legend = TRUE, show_annotation_legend = TRUE)
  message("ht means done")
  
  png("htmaps_mean_big_bars_usual_nolegend.png", res = 300, units = "in", width = 8, height = 7)
  print(ht1)
  dev.off()
  
  png("htmaps_mean_big_bars_usual.png", res = 300, units = "in", width = 8, height = 7)
  print(ht2)
  dev.off()
  
  pdf("htmaps_mean_big_bars.pdf", width = 8, height = 7)
  print(ht1)
  dev.off()
  
}

# heatmaps mean main figure---------------
library(ComplexHeatmap)
run_mean_htmaps_fig1 <- function(df, seu, cell, genes_df, genes_df_cann, clus_levels, genes_sel_df) {
  
  # exp scores -------
  mat_exp = seu@assays$integrated@scale.data
  mat_exp = data.matrix(mat_exp)
  message("exp scores done")
  
  # exp scores -------
  mat_exp_rna = seu@assays$RNA@scale.data
  mat_exp_rna = data.matrix(mat_exp_rna)
  message("exp scores done")
  
  # cannonical genes
  d1_cann = mat_exp_rna %>% as_tibble() %>% dplyr::mutate(genes = rownames(mat_exp_rna))
  d2_cann = d1_cann %>% dplyr::filter(genes %in% genes_df_cann)
  mat_exp_fil_cann = data.matrix(d2_cann)
  mat_exp_fil_cann = suppressWarnings(mat_exp_fil_cann[,-ncol(mat_exp_fil_cann)])
  rownames(mat_exp_fil_cann) = d2_cann$genes
  mat_exp_fil_cann = mat_exp_fil_cann[genes_df_cann,] #match the order of orginal genes
  message("cannonical done")
  
  # main genes body ht
  d1 = mat_exp_rna %>% as_tibble() %>% dplyr::mutate(genes = rownames(mat_exp_rna))
  d2 = d1 %>% dplyr::filter(genes %in% genes_df)
  mat_exp_fil = data.matrix(d2)
  mat_exp_fil = suppressWarnings(mat_exp_fil[,-ncol(mat_exp_fil)])
  rownames(mat_exp_fil) = d2$genes
  mat_exp_fil = mat_exp_fil[genes_df,] #match the order of orginal genes
  message("main done")
  
  # Annotations -------
  # get df for annotation
  coldata = data.frame(Cluster = seu$final_group1, #minor_group_factor
                       # Protocol = seu@meta.data$exp_proc,
                       # Menopause = seu@meta.data$menopause.y,
                       # Paridy = seu@meta.data$parity.1.y,
                       # Procedure = seu@meta.data$figure1_grp_updated,
                       sample_id = colnames(mat_exp_fil),
                       stringsAsFactors = F) %>% dplyr::arrange(Cluster)
  
  # top annotation
  ha_top = dplyr::select(coldata, Cluster) %>% 
    HeatmapAnnotation(df = ., show_annotation_name = FALSE, 
                      col = list(Cluster = col_grp))
                                 # Procedure = col_procedure,
                                 # Menopause = col_meno,
                                 # Paridy = col_paridy,
                                 # Protocol = col_Protocol))
  
  # get a sorted list of samples
  samps_sorted = coldata$sample_id
  
  # rearrnge the matrix
  mat = mat_exp_fil[, samps_sorted]
  mat0 = mat_exp_fil_cann[, samps_sorted]
  
  # # heatmap of exp
  # ht_exp = Heatmap(data.matrix(mat),
  #                  column_title = "Scaled Expression",
  #                  column_title_gp = gpar(fontsize = 6),
  #                  name = "Scaled Expression",
  #                  top_annotation = ha_top,show_row_names = T,
  #                  show_column_names = F, cluster_rows = F,
  #                  use_raster = F, cluster_columns = FALSE,
  #                  col = colorRamp2(c(-1.5, 0, 1.5), c("steelblue","white", "tomato3")),
  #                  row_names_gp = gpar(fontsize = 8),  width = unit(8, "cm"))
  
  # Mean cann ht ------
  # get the mean of exp 
  scores_cann = cbind(coldata, t(mat0))#ma
  scores_cann = scores_cann[,colnames(unique(as.matrix(scores_cann),MARGIN=2))] 
  scores_m_cann = suppressWarnings(scores_cann %>% group_by(Cluster) %>% unique() %>% summarise_all(funs(mean)))
  mat_exp_mean_cann = data.matrix(t(scores_m_cann))
  mat_exp_mean_cann = mat_exp_mean_cann[genes_df_cann,]
  #mat_exp_mean = mat_exp_mean[-c(1,2), ]
  colnames(mat_exp_mean_cann) = levels(scores_cann$Cluster)
  mat_exp_mean_cann = apply(mat_exp_mean_cann, 2,as.numeric)
  rownames(mat_exp_mean_cann) = genes_df_cann
  
  # Mean ht ------
  # get the mean of exp 
  scores = cbind(coldata, t(mat))#ma
  #scores = t(mat)#ma
  scores = scores[,colnames(unique(as.matrix(scores),MARGIN=2))] 
  scores_m = suppressWarnings(scores %>% group_by(Cluster) %>% unique() %>% summarise_all(funs(mean)))
  mat_exp_mean = data.matrix(t(scores_m))
  mat_exp_mean = mat_exp_mean[genes_df,]
  #mat_exp_mean = mat_exp_mean[-c(1,2), ]
  colnames(mat_exp_mean) = levels(scores$Cluster)
  mat_exp_mean = apply(mat_exp_mean, 2,as.numeric)
  rownames(mat_exp_mean) = genes_df
  
  # mean ht row annotation ------
  df_genes_sel = df %>% dplyr::filter(gene %in% rownames(mat_exp_mean)) %>% dplyr::select(gene, cluster1) %>% mutate(id = paste0(gene, "-", cluster1))
  
  # for figure1
  #tff3 tcim stc2 saa1 pdk4 mt1a mcam myve1 krt23 igfbp5 hladqa cldn5 cd93 cd69 cavin1adamst4d
  # df_genes_sel1 = df_genes_sel %>% dplyr::filter(id %in% genes_sel_df$id)
  # seq_clusters = dplyr::select(coldata, Cluster) %>% pull(Cluster) %>% unique() %>% as.character()
  # df_genes_sel2 = df_genes_sel1 %>% mutate(cluster1 = forcats::fct_relevel(cluster, levels = seq_clusters)) %>% arrange(cluster1)
  # df_genes_sel2 = df_genes_sel2[-c(8,9,16,17,24, 33,34,35,49, 56, 63,64,65,78, 79),]
  # df_genes_sel2 = df_genes_sel2[-c(50,43),]
  
  
  ha_row = df_genes_sel %>% dplyr::select(cluster1) %>% 
       rowAnnotation(df = ., show_annotation_name = FALSE, width =  unit(0.1, "cm"),
                         col = list(cluster1 = col_grp))
 
  # for endothelial
  # df_genes_sel1 = df_genes_sel[-6,] 
  # ha_row = df_genes_sel1 %>% dplyr::select(cluster) %>% 
  #   rowAnnotation(df = ., show_annotation_name = FALSE, width =  unit(0.1, "cm"),
  #                     col = list(cluster = col_grp))
   
  # heatmap of exp
   ha = HeatmapAnnotation(#Menopause = anno_barplot(t(df_meta_meno2), border = FALSE, gp = gpar(fill = col_meno, border = NA, lty = "blank"), axis = FALSE, size = unit(0.85, "mm"), bar_width = 0.75),
  #                        Paridy = anno_barplot(t(df_meta_paridy2), border = FALSE, gp = gpar(fill = col_paridy, border = NA, lty = "blank"), axis = FALSE,bar_width = 0.85),
  #                        Procedure = anno_barplot(t(df_meta_grp2), border = FALSE, gp = gpar(fill = col_procedure, border = NA, lty = "blank"),axis = FALSE,bar_width = 0.85), 
  #                        Protocol = anno_barplot(t(df_meta_dig2), border = FALSE, gp = gpar(fill = col_Protocol, border = NA, lty = "blank"),axis = FALSE,bar_width = 0.85),
                          #pt = anno_points(t(mat_exp_mean_cann)[,1], pch = 21, col = c("red", "blue")),
                          ggplot = anno_empty(height = unit(1, "cm")),show_annotation_name = FALSE)
                          #Cluster = levels(scores$Cluster), col = list(Cluster = col_grp), )
  #ha_row = HeatmapAnnotation(Cluster = levels(scores$Cluster), col = list(Cluster = col), show_annotation_name = TRUE)
   
   ha1 = HeatmapAnnotation(Cluster = levels(scores$Cluster), col = list(Cluster = col_grp), show_annotation_name = FALSE)
   
   # df_exp_mean_cann = mat_exp_mean_cann %>% as.data.frame()
   # df_exp_mean_cann$gene = rownames(df_exp_mean_cann)
   # df_exp_mean_cann1= reshape::melt(df_exp_mean_cann, id.vars=c("gene"))
   # df_exp_mean_cann1$cluster = df_exp_mean_cann1$variable
   # g<-ggplot(df_exp_mean_cann1, aes(x=cluster, y=gene, color=value, size=value)) + scale_radius(range=c(0, 3)) +
   #   geom_point() + theme_void() + theme(legend.position="none") + scale_color_gradient2(low = "steelblue", mid = "white", high = "tomato3")  
   #   
   # 
   
  ht_mean_cann = Heatmap(data.matrix(mat_exp_mean_cann),
                    #column_title = "Avg Exp",column_title_gp = gpar(fontsize = 6),
                    name = "Avg Expression",
                    top_annotation = ha1,
                    #left_annotation = ha_row1,
                    show_row_names = T, row_names_side = c("right"),
                    show_column_names = F, show_row_dend = F, cluster_rows = F,
                    use_raster = T, cluster_columns = FALSE,
                    #col = colorRamp2(c(0, 0.5, 1), c("steelblue","white", "tomato3")),
                    col = colorRamp2(c(-1.5, 0, 1.5), c("steelblue","white", "tomato3")),
                    #col = colorRamp2(c(-1.5, 0, 1.5), c("black","black", "grey99")),
                    row_names_gp = gpar(fontsize = 4),  width = unit(4.5, "cm"), height = unit(1, "cm"))
  ht_mean_cann
  
  ht_mean = Heatmap(data.matrix(mat_exp_mean),
                    name = "Avg Expression",
                    top_annotation = ha1, #ha for endothelial
                    left_annotation = ha_row,
                    show_row_names = T, row_names_side = c("right"), 
                    show_column_names = F, show_row_dend = F, cluster_rows = F,
                    use_raster = T, cluster_columns = FALSE,
                    # heatmap_legend_param = list(legend_direction = "horizontal",
                    #                             legend_width = unit(3, "cm"), title_position = "lefttop"),
                    col = colorRamp2(c(-1.5, 0, 2), c("steelblue","white", "tomato3")),
                    #col = colorRamp2(c(-1, 0, 2), c("cadetblue4","white", "coral3")),
                    #col = colorRamp2(c(-1.5, 0, 1.5), c("black","black", "grey99")),
                    row_names_gp = gpar(fontsize = 4),  width = unit(5, "cm"), height = unit(11, "cm"))
  ht_mean
  message("ht means done")
  
  #ht = draw(ht_mean_cann %v% ht_mean, heatmap_legend_side = "bottom")
  ht = draw(ht_mean, heatmap_legend_side = "bottom")
  png("htmaps_nucmean_v1_40_40_30_pct225.png", res = 300, units = "in", width = 8, height = 8.5)
  #ht = draw(ha1 %v% ht_mean, heatmap_legend_side = "bottom")
  print(ht)
  dev.off()
  
  pdf("htmaps_nuc_mean_v1_40_40_30_pct225_with10more.pdf", width = 8, height = 8.5)
  #ht = draw(ha1 %v% ht_mean, heatmap_legend_side = "bottom")
  print(ht)
  dev.off()
  
  
  # pdf("htmaps_mean_small_bars.pdf", width = 8, height = 7)
  # print(ht)
  # dev.off()
  
}

run_mean_htmaps_fig1_supp <- function(df, seu, cell, genes_df, genes_df_cann, clus_levels, genes_sel_df) {
  
  
  # exp scores -------
  mat_exp_rna = seu@assays$integrated@scale.data
  mat_exp_rna = data.matrix(mat_exp_rna)
  message("exp scores done")
  
  # cannonical genes
  d1_cann = mat_exp_rna %>% as_tibble() %>% dplyr::mutate(genes = rownames(mat_exp_rna))
  d2_cann = d1_cann %>% dplyr::filter(genes %in% genes_df_cann)
  mat_exp_fil_cann = data.matrix(d2_cann)
  mat_exp_fil_cann = suppressWarnings(mat_exp_fil_cann[,-ncol(mat_exp_fil_cann)])
  rownames(mat_exp_fil_cann) = d2_cann$genes
  mat_exp_fil_cann = mat_exp_fil_cann[genes_df_cann,] #match the order of orginal genes
  message("cannonical done")
  
  # main genes body ht
  d1 = mat_exp_rna %>% as_tibble() %>% dplyr::mutate(genes = rownames(mat_exp_rna))
  d2 = d1 %>% dplyr::filter(genes %in% genes_df)
  mat_exp_fil = data.matrix(d2)
  mat_exp_fil = suppressWarnings(mat_exp_fil[,-ncol(mat_exp_fil)])
  rownames(mat_exp_fil) = d2$genes
  mat_exp_fil = mat_exp_fil[genes_df,] #match the order of orginal genes
  message("main done")
  
  # Annotations -------
  # get df for annotation
  coldata = data.frame(Cluster = seu$final_group, #minor_group_factor
                       # Protocol = seu@meta.data$exp_proc,
                       # Menopause = seu@meta.data$menopause.y,
                       # Paridy = seu@meta.data$parity.1.y,
                       # Procedure = seu@meta.data$figure1_grp_updated,
                       sample_id = colnames(mat_exp_fil),
                       stringsAsFactors = F) %>% dplyr::arrange(Cluster)
  
  # top annotation
  ha_top = dplyr::select(coldata, Cluster) %>% 
    HeatmapAnnotation(df = ., show_annotation_name = FALSE, 
                      col = list(Cluster = col_grp))
  
  # get a sorted list of samples
  samps_sorted = coldata$sample_id
  
  # rearrnge the matrix
  mat = mat_exp_fil[, samps_sorted]
  mat0 = mat_exp_fil_cann[, samps_sorted]
  
  # Mean cann ht ------
  # get the mean of exp 
  scores_cann = cbind(coldata, t(mat0))#ma
  scores_cann = scores_cann[,colnames(unique(as.matrix(scores_cann),MARGIN=2))] 
  scores_m_cann = suppressWarnings(scores_cann %>% group_by(Cluster) %>% unique() %>% summarise_all(funs(mean)))
  mat_exp_mean_cann = data.matrix(t(scores_m_cann))
  mat_exp_mean_cann = mat_exp_mean_cann[genes_df_cann,]
  #mat_exp_mean = mat_exp_mean[-c(1,2), ]
  colnames(mat_exp_mean_cann) = levels(scores_cann$Cluster)
  mat_exp_mean_cann = apply(mat_exp_mean_cann, 2,as.numeric)
  rownames(mat_exp_mean_cann) = genes_df_cann
  
  # Mean ht ------
  # only get 100 cells
  seu$cellname = rownames(seu@meta.data)
  data.frame(new_cluster_factor = seu@meta.data$final_group, cell = seu$cellname) %>%
    #filter(!str_detect(cell, "untreated")) %>% 
    group_by(new_cluster_factor) %>% 
    sample_n(500) %>% 
    mutate(cell = as.character(cell)) %>% 
    pull(cell) ->
    tum_cells
  scores = cbind(coldata, t(mat))#ma
  mat_exp_mean = data.matrix(t(scores))
  mat_exp_mean = mat_exp_mean[genes_df,]
  mat_exp_mean = mat_exp_mean[,tum_cells]
  mat_exp_mean = apply(mat_exp_mean, 2,as.numeric)
  rownames(mat_exp_mean) = genes_df
  
  
  
  # mean ht row annotation ------
  df_genes_sel = df %>% dplyr::filter(gene %in% rownames(mat_exp_mean)) %>% dplyr::select(gene, cluster1) %>% mutate(id = paste0(gene, "-", cluster1))
  
  ha_row = df_genes_sel %>% dplyr::select(cluster1) %>% 
    rowAnnotation(df = ., show_annotation_name = FALSE, width =  unit(0.1, "cm"),
                  col = list(cluster1 = col_grp))

  coldata1 = coldata %>% dplyr::filter(sample_id %in% colnames(mat_exp_mean))
  ha1 = HeatmapAnnotation(Cluster = coldata1$Cluster, col = list(Cluster = col_grp), show_annotation_name = FALSE)
  
  # ht_mean_cann = Heatmap(data.matrix(mat_exp_mean_cann),
  #                        #column_title = "Avg Exp",column_title_gp = gpar(fontsize = 6),
  #                        name = "Avg Expression",
  #                        top_annotation = ha1,
  #                        #left_annotation = ha_row1,
  #                        show_row_names = T, row_names_side = c("right"),
  #                        show_column_names = F, show_row_dend = F, cluster_rows = F,
  #                        use_raster = T, cluster_columns = FALSE,
  #                        #col = colorRamp2(c(0, 0.5, 1), c("steelblue","white", "tomato3")),
  #                        col = colorRamp2(c(-1.5, 0, 1.5), c("steelblue","white", "tomato3")),
  #                        #col = colorRamp2(c(-1.5, 0, 1.5), c("black","black", "grey99")),
  #                        row_names_gp = gpar(fontsize = 4),  width = unit(4.5, "cm"), height = unit(1, "cm"))
  # ht_mean_cann
  
  ht_mean = Heatmap(data.matrix(mat_exp_mean),
                    name = "Avg Expression",
                    top_annotation = ha1, #ha for endothelial
                    #left_annotation = ha_row,
                    show_row_names = T, row_names_side = c("right"), 
                    show_column_names = F, show_row_dend = F, cluster_rows = F,
                    use_raster = F, cluster_columns = FALSE,
                    # heatmap_legend_param = list(legend_direction = "horizontal",
                    #                             legend_width = unit(3, "cm"), title_position = "lefttop"),
                    col = colorRamp2(c(-1.5, 0, 1), c("steelblue4","white", "tomato3")),
                    #col = colorRamp2(c(-1, 0, 2), c("cadetblue4","white", "coral3")),
                    #col = colorRamp2(c(-1.5, 0, 1.5), c("black","black", "grey99")),
                    row_names_gp = gpar(fontsize = 5),  width = unit(6, "cm"), height = unit(16, "cm"))
  ht_mean
  message("ht means done")
  
  #ht = draw(ht_mean_cann %v% ht_mean, heatmap_legend_side = "bottom")
  ht = draw(ht_mean, heatmap_legend_side = "bottom")
  png("supp_htmaps_nucmean_v1_40_40_30_pct225.png", res = 300, units = "in", width = 8, height = 8.5)
  #ht = draw(ha1 %v% ht_mean, heatmap_legend_side = "bottom")
  print(ht)
  dev.off()
  
  pdf("supp_htmaps_cell_mean_v1_40_40_30_pct225_with10more.pdf", width = 8, height = 12)
  #ht = draw(ha1 %v% ht_mean, heatmap_legend_side = "bottom")
  print(ht)
  dev.off()
  
  
  # pdf("htmaps_mean_small_bars.pdf", width = 8, height = 7)
  # print(ht)
  # dev.off()
  
}


run_mean_htmaps_fig1_horizontal <- function(df, seu, cell, genes_df, genes_df_cann, clus_levels, genes_sel_df) {
  
  # exp scores -------
  mat_exp_rna = seu@assays$RNA@scale.data
  mat_exp_rna = data.matrix(mat_exp_rna)
  message("exp scores done")
  
  # main genes body ht
  d1 = mat_exp_rna %>% as_tibble() %>% dplyr::mutate(genes = rownames(mat_exp_rna))
  d2 = d1 %>% dplyr::filter(genes %in% genes_df)
  mat_exp_fil = data.matrix(d2)
  mat_exp_fil = suppressWarnings(mat_exp_fil[,-ncol(mat_exp_fil)])
  rownames(mat_exp_fil) = d2$genes
  mat_exp_fil = mat_exp_fil[genes_df,] #match the order of orginal genes
  message("main done")
  
  # Annotations -------
  # get df for annotation
  coldata = data.frame(Cluster = seu$final_group1, #minor_group_factor
                       # Protocol = seu@meta.data$exp_proc,
                       # Menopause = seu@meta.data$menopause.y,
                       # Paridy = seu@meta.data$parity.1.y,
                       # Procedure = seu@meta.data$figure1_grp_updated,
                       sample_id = colnames(mat_exp_fil),
                       stringsAsFactors = F) %>% dplyr::arrange(Cluster)
  
  # top annotation
  ha_top = dplyr::select(coldata, Cluster) %>% 
    HeatmapAnnotation(df = ., show_annotation_name = FALSE, 
                      col = list(Cluster = col_grp))

  
  # get a sorted list of samples
  samps_sorted = coldata$sample_id
  
  # rearrnge the matrix
  mat = mat_exp_fil[, samps_sorted]
  
  # Mean ht ------
  # get the mean of exp 
  scores = cbind(coldata, t(mat))#ma
  #scores = t(mat)#ma
  scores = scores[,colnames(unique(as.matrix(scores),MARGIN=2))] 
  scores_m = suppressWarnings(scores %>% group_by(Cluster) %>% unique() %>% summarise_all(funs(mean)))
  mat_exp_mean = data.matrix(t(scores_m))
  mat_exp_mean = mat_exp_mean[genes_df,]
  #mat_exp_mean = mat_exp_mean[-c(1,2), ]
  colnames(mat_exp_mean) = levels(scores$Cluster)
  mat_exp_mean = apply(mat_exp_mean, 2,as.numeric)
  rownames(mat_exp_mean) = genes_df
  
  # mean ht row annotation ------
  df_genes_sel = df %>% dplyr::filter(gene %in% rownames(mat_exp_mean)) %>% dplyr::select(gene, cluster1) %>% mutate(id = paste0(gene, "-", cluster1))
  
  ha_row = df_genes_sel %>% dplyr::select(cluster1) %>% 
    HeatmapAnnotation(df = ., show_annotation_name = FALSE, 
                  col = list(cluster1 = col_grp))
  
  ha1 = rowAnnotation(Cluster = levels(scores$Cluster), col = list(Cluster = col_grp), show_annotation_name = FALSE)
  
  #dim(t(data.matrix(mat_exp_mean)))
  ht_mean = Heatmap(t(data.matrix(mat_exp_mean)),
                    name = "Avg Expression",
                    # row params
                    #left_annotation = ha1, 
                    show_row_dend = F, cluster_rows = F,
                    show_row_names = F, #row_names_side = c("right"), 
                    # row params
                    top_annotation = ha_row,
                    cluster_columns = FALSE, 
                    column_names_gp = gpar(fontsize = 5),
                    column_names_rot = 65,
                    show_column_names = T, 
                    
                    use_raster = T, 
                    col = colorRamp2(c(-1.5, 0, 2), c("steelblue","white", "tomato3")),
                    width = unit(15, "cm"), height = unit(4, "cm"))
  ht_mean
  message("ht means done")
  
  ht = draw(ht_mean, heatmap_legend_side = "bottom")
  png("htmaps_nucmean_v1_40_40_30_pct225_horizontal.png", res = 300, units = "in", width = 12.5, height = 8.5)
  print(ht)
  dev.off()
  
  pdf("htmaps_nuc_mean_v1_40_40_30_pct225_with10more_horizontal_65.pdf", width = 12.5, height = 8)
  print(ht)
  dev.off()

  
}

run_mean_htmaps_fig1_horizontal_supp <- function(df, seu, cell, genes_df, genes_df_cann, clus_levels, genes_sel_df) {
  
  # exp scores -------
  mat_exp_rna = seu@assays$RNA@scale.data
  mat_exp_rna = data.matrix(mat_exp_rna)
  message("exp scores done")
  
  # main genes body ht
  d1 = mat_exp_rna %>% as_tibble() %>% dplyr::mutate(genes = rownames(mat_exp_rna))
  d2 = d1 %>% dplyr::filter(genes %in% genes_df)
  mat_exp_fil = data.matrix(d2)
  mat_exp_fil = suppressWarnings(mat_exp_fil[,-ncol(mat_exp_fil)])
  rownames(mat_exp_fil) = d2$genes
  mat_exp_fil = mat_exp_fil[genes_df,] #match the order of orginal genes
  message("main done")
  
  # Annotations -------
  # get df for annotation
  coldata = data.frame(Cluster = seu$final_group1, #minor_group_factor
                       # Protocol = seu@meta.data$exp_proc,
                       # Menopause = seu@meta.data$menopause.y,
                       # Paridy = seu@meta.data$parity.1.y,
                       # Procedure = seu@meta.data$figure1_grp_updated,
                       sample_id = colnames(mat_exp_fil),
                       stringsAsFactors = F) %>% dplyr::arrange(Cluster)
  
  # top annotation
  ha_top = dplyr::select(coldata, Cluster) %>% 
    HeatmapAnnotation(df = ., show_annotation_name = FALSE, 
                      col = list(Cluster = col_grp))
  
  
  # get a sorted list of samples
  samps_sorted = coldata$sample_id
  
  # rearrnge the matrix
  mat = mat_exp_fil[, samps_sorted]
  
  # Mean ht ------
  # only get 100 cells
  seu$cellname = rownames(seu@meta.data)
  data.frame(new_cluster_factor = seu@meta.data$final_group1, cell = seu$cellname) %>%
    #filter(!str_detect(cell, "untreated")) %>% 
    group_by(new_cluster_factor) %>% 
    sample_n(500) %>% 
    mutate(cell = as.character(cell)) %>% 
    pull(cell) ->
    tum_cells
  scores = cbind(coldata, t(mat))#ma
  mat_exp_mean = data.matrix(t(scores))
  mat_exp_mean = mat_exp_mean[genes_df,]
  mat_exp_mean = mat_exp_mean[,tum_cells]
  mat_exp_mean = apply(mat_exp_mean, 2,as.numeric)
  rownames(mat_exp_mean) = genes_df
  
  # mean ht row annotation ------
  df_genes_sel = df %>% dplyr::filter(gene %in% rownames(mat_exp_mean)) %>% dplyr::select(gene, cluster1) %>% mutate(id = paste0(gene, "-", cluster1))
  
  ha_row = df_genes_sel %>% dplyr::select(cluster1) %>% 
    HeatmapAnnotation(df = ., show_annotation_name = FALSE, 
                      col = list(cluster1 = col_grp))
  
  ha1 = rowAnnotation(Cluster = levels(scores$Cluster), col = list(Cluster = col_grp), show_annotation_name = FALSE)
  
  #dim(t(data.matrix(mat_exp_mean)))
  ht_mean = Heatmap(t(data.matrix(mat_exp_mean)),
                    name = "Avg Expression",
                    # row params
                    #left_annotation = ha1, 
                    show_row_dend = F, cluster_rows = F,
                    show_row_names = F, #row_names_side = c("right"), 
                    # row params
                    top_annotation = ha_row,
                    cluster_columns = FALSE, 
                    column_names_gp = gpar(fontsize = 5),
                    column_names_rot = 65,
                    show_column_names = T, 
                    
                    use_raster = T, 
                    col = colorRamp2(c(-1.5, 0, 2), c("steelblue","white", "tomato3")),
                    width = unit(18, "cm"), height = unit(5, "cm"))
  ht_mean
  message("ht means done")
  
  ht = draw(ht_mean, heatmap_legend_side = "bottom")
  png("supp_htmaps_cellmean_v1_40_40_30_pct225_horizontal.png", res = 300, units = "in", width = 12.5, height = 8.5)
  print(ht)
  dev.off()
  
  pdf("supp_htmaps_cell_mean_v1_40_40_30_pct225_with10more_horizontal_65.pdf", width = 12.5, height = 8)
  print(ht)
  dev.off()
  
  
}

run_mean_htmaps_fig_endo <- function(df, seu, cell, genes_df, genes_df_cann, clus_levels, genes_sel_df) {
  
  # exp scores -------
  mat_exp = seu@assays$integrated@scale.data
  mat_exp = data.matrix(mat_exp)
  message("exp scores done")
  
  # exp scores -------
  mat_exp_rna = seu@assays$RNA@scale.data
  mat_exp_rna = data.matrix(mat_exp_rna)
  message("exp scores done")
  
  # cannonical genes
  d1_cann = mat_exp_rna %>% as_tibble() %>% dplyr::mutate(genes = rownames(mat_exp_rna))
  d2_cann = d1_cann %>% dplyr::filter(genes %in% genes_df_cann)
  mat_exp_fil_cann = data.matrix(d2_cann)
  mat_exp_fil_cann = suppressWarnings(mat_exp_fil_cann[,-ncol(mat_exp_fil_cann)])
  rownames(mat_exp_fil_cann) = d2_cann$genes
  mat_exp_fil_cann = mat_exp_fil_cann[genes_df_cann,] #match the order of orginal genes
  message("cannonical done")
  
  # main genes body ht
  d1 = mat_exp_rna %>% as_tibble() %>% dplyr::mutate(genes = rownames(mat_exp_rna))
  d2 = d1 %>% dplyr::filter(genes %in% genes_df)
  mat_exp_fil = data.matrix(d2)
  mat_exp_fil = suppressWarnings(mat_exp_fil[,-ncol(mat_exp_fil)])
  rownames(mat_exp_fil) = d2$genes
  mat_exp_fil = mat_exp_fil[genes_df,] #match the order of orginal genes
  message("main done")
  
  # Annotations -------
  # get df for annotation
  coldata = data.frame(Cluster = seu$minor_group_factor, #minor_group_factor
                       # Protocol = seu@meta.data$exp_proc,
                       # Menopause = seu@meta.data$menopause.y,
                       # Paridy = seu@meta.data$parity.1.y,
                       # Procedure = seu@meta.data$figure1_grp_updated,
                       sample_id = colnames(mat_exp_fil),
                       stringsAsFactors = F) %>% dplyr::arrange(Cluster)
  
  # top annotation
  ha_top = dplyr::select(coldata, Cluster) %>% 
    HeatmapAnnotation(df = ., show_annotation_name = FALSE, 
                      col = list(Cluster = col_grp))
  # Procedure = col_procedure,
  # Menopause = col_meno,
  # Paridy = col_paridy,
  # Protocol = col_Protocol))
  
  # get a sorted list of samples
  samps_sorted = coldata$sample_id
  
  # rearrnge the matrix
  mat = mat_exp_fil[, samps_sorted]
  mat0 = mat_exp_fil_cann[, samps_sorted]
  
  # # heatmap of exp
  # ht_exp = Heatmap(data.matrix(mat),
  #                  column_title = "Scaled Expression",
  #                  column_title_gp = gpar(fontsize = 6),
  #                  name = "Scaled Expression",
  #                  top_annotation = ha_top,show_row_names = T,
  #                  show_column_names = F, cluster_rows = F,
  #                  use_raster = F, cluster_columns = FALSE,
  #                  col = colorRamp2(c(-1.5, 0, 1.5), c("steelblue","white", "tomato3")),
  #                  row_names_gp = gpar(fontsize = 8),  width = unit(8, "cm"))
  
  # Mean cann ht ------
  # get the mean of exp 
  scores_cann = cbind(coldata, t(mat0))#ma
  scores_cann = scores_cann[,colnames(unique(as.matrix(scores_cann),MARGIN=2))] 
  scores_m_cann = suppressWarnings(scores_cann %>% group_by(Cluster) %>% unique() %>% summarise_all(funs(mean)))
  mat_exp_mean_cann = data.matrix(t(scores_m_cann))
  mat_exp_mean_cann = mat_exp_mean_cann[genes_df_cann,]
  #mat_exp_mean = mat_exp_mean[-c(1,2), ]
  colnames(mat_exp_mean_cann) = levels(scores_cann$Cluster)
  mat_exp_mean_cann = apply(mat_exp_mean_cann, 2,as.numeric)
  rownames(mat_exp_mean_cann) = genes_df_cann
  
  # Mean ht ------
  # get the mean of exp 
  scores = cbind(coldata, t(mat))#ma
  #scores = t(mat)#ma
  scores = scores[,colnames(unique(as.matrix(scores),MARGIN=2))] 
  scores_m = suppressWarnings(scores %>% group_by(Cluster) %>% unique() %>% summarise_all(funs(mean)))
  mat_exp_mean = data.matrix(t(scores_m))
  mat_exp_mean = mat_exp_mean[genes_df,]
  #mat_exp_mean = mat_exp_mean[-c(1,2), ]
  colnames(mat_exp_mean) = levels(scores$Cluster)
  mat_exp_mean = apply(mat_exp_mean, 2,as.numeric)
  rownames(mat_exp_mean) = genes_df
  
  # mean ht row annotation ------
  df_genes_sel = df %>% dplyr::filter(gene %in% rownames(mat_exp_mean)) %>% dplyr::select(gene, cluster) %>% mutate(id = paste0(gene, "-", cluster))
  
  # for figure1
  df_genes_sel1 = df_genes_sel %>% dplyr::filter(id %in% genes_sel_df$id)
  ha_row = df_genes_sel1 %>% dplyr::select(cluster) %>% 
    rowAnnotation(df = ., show_annotation_name = FALSE, width =  unit(0.1, "cm"),
                  col = list(cluster = col_grp))
  
  # for endothelial
  # df_genes_sel1 = df_genes_sel[-6,] 
   ha_row = df_genes_sel %>% dplyr::select(cluster) %>% 
     rowAnnotation(df = ., show_annotation_name = FALSE, width =  unit(0.1, "cm"),
                       col = list(cluster = col_grp))
  
  # heatmap of exp
  ha = HeatmapAnnotation(#Menopause = anno_barplot(t(df_meta_meno2), border = FALSE, gp = gpar(fill = col_meno, border = NA, lty = "blank"), axis = FALSE, size = unit(0.85, "mm"), bar_width = 0.75),
    #                        Paridy = anno_barplot(t(df_meta_paridy2), border = FALSE, gp = gpar(fill = col_paridy, border = NA, lty = "blank"), axis = FALSE,bar_width = 0.85),
    #                        Procedure = anno_barplot(t(df_meta_grp2), border = FALSE, gp = gpar(fill = col_procedure, border = NA, lty = "blank"),axis = FALSE,bar_width = 0.85), 
    #                        Protocol = anno_barplot(t(df_meta_dig2), border = FALSE, gp = gpar(fill = col_Protocol, border = NA, lty = "blank"),axis = FALSE,bar_width = 0.85),
    #pt = anno_points(t(mat_exp_mean_cann)[,1], pch = 21, col = c("red", "blue")),
    ggplot = anno_empty(height = unit(1, "cm")),show_annotation_name = FALSE)
  #Cluster = levels(scores$Cluster), col = list(Cluster = col_grp), )
  #ha_row = HeatmapAnnotation(Cluster = levels(scores$Cluster), col = list(Cluster = col), show_annotation_name = TRUE)
  
  ha1 = HeatmapAnnotation(Cluster = levels(scores$Cluster), col = list(Cluster = col_grp), show_annotation_name = FALSE)
  
  df_exp_mean_cann = mat_exp_mean_cann %>% as.data.frame()
  df_exp_mean_cann$gene = rownames(df_exp_mean_cann)
  df_exp_mean_cann1= reshape::melt(df_exp_mean_cann, id.vars=c("gene"))
  df_exp_mean_cann1$cluster = df_exp_mean_cann1$variable
  g<-ggplot(df_exp_mean_cann1, aes(x=cluster, y=gene, color=value, size=value)) + scale_radius(range=c(0, 3)) +
    geom_point() + theme_void() + theme(legend.position="none") + scale_color_gradient2(low = "steelblue", mid = "white", high = "tomato3")  
  
  
  
  ht_mean_cann = Heatmap(data.matrix(mat_exp_mean_cann),
                    #column_title = "Avg Exp",column_title_gp = gpar(fontsize = 6),
                    name = "Avg Expression",
                    top_annotation = ha1,

                    show_row_names = T, row_names_side = c("right"),
                    show_column_names = F, show_row_dend = F, cluster_rows = F,
                    use_raster = T, cluster_columns = FALSE,
                    col = colorRamp2(c(-1, 0, 1), c("steelblue","white", "tomato3")),
                    #col = colorRamp2(c(-1, 0, 2), c("cadetblue4","white", "coral3")),
                    #col = colorRamp2(c(-1.5, 0, 1.5), c("black","black", "grey99")),
                    row_names_gp = gpar(fontsize = 5),  width = unit(2.5, "cm"), height = unit(0.75, "cm"))
  ht_mean_cann
  
  ht_mean = Heatmap(data.matrix(mat_exp_mean),
                    #column_title = "Avg Exp",column_title_gp = gpar(fontsize = 6),
                    name = "Avg Expression",
                    top_annotation = ha1, #ha for endothelial
                    left_annotation = ha_row,
                    show_row_names = T, row_names_side = c("right"), 
                    show_column_names = F, show_row_dend = F, cluster_rows = F,
                    use_raster = T, cluster_columns = FALSE,
                    # heatmap_legend_param = list(legend_direction = "horizontal",
                    #                             legend_width = unit(3, "cm"), title_position = "lefttop"),
                    col = colorRamp2(c(-1, 0, 1), c("steelblue","white", "tomato3")),
                    #col = colorRamp2(c(-1, 0, 2), c("cadetblue4","white", "coral3")),
                    #col = colorRamp2(c(-1.5, 0, 1.5), c("black","black", "grey99")),
                    row_names_gp = gpar(fontsize = 5),  width = unit(2.5, "cm"), height = unit(4, "cm"))
  ht_mean
  message("ht means done")
  
  # png("htmaps_mean_no_bars_canno_integrated_dot_leftbar.png", res = 300, units = "in", width = 5, height = 5.5)
  # ht = draw(ha1 %v% ht_mean, heatmap_legend_side = "bottom")
  # decorate_annotation("ggplot", {
  #   vp = current.viewport()$name
  #   print(g, vp = vp)
  #   
  # })
  # dev.off()
  # 
  #ht = draw(ht_mean_cann %v% ht_mean, heatmap_legend_side = "bottom")
  ht = draw(ht_mean, heatmap_legend_side = "bottom")
  png("htmaps_mean_v1_updated.png", res = 300, units = "in", width = 4, height = 4)
  #ht = draw(ha1 %v% ht_mean, heatmap_legend_side = "bottom")
  print(ht)
  dev.off()
  
  
  
  # pdf("htmaps_mean_small_bars.pdf", width = 8, height = 7)
  # print(ht)
  # dev.off()
  
}

run_mean_htmaps_fig_perio <- function(df, seu, cell, genes_df, genes_df_cann, clus_levels, genes_sel_df) {
  
  # exp scores -------
  mat_exp = seu@assays$integrated@scale.data
  mat_exp = data.matrix(mat_exp)
  message("exp scores done")
  
  # exp scores -------
  mat_exp_rna = seu@assays$RNA@scale.data
  mat_exp_rna = data.matrix(mat_exp_rna)
  message("exp scores done")
  
  # cannonical genes
  d1_cann = mat_exp_rna %>% as_tibble() %>% dplyr::mutate(genes = rownames(mat_exp_rna))
  d2_cann = d1_cann %>% dplyr::filter(genes %in% genes_df_cann)
  mat_exp_fil_cann = data.matrix(d2_cann)
  mat_exp_fil_cann = suppressWarnings(mat_exp_fil_cann[,-ncol(mat_exp_fil_cann)])
  rownames(mat_exp_fil_cann) = d2_cann$genes
  mat_exp_fil_cann = mat_exp_fil_cann[genes_df_cann,] #match the order of orginal genes
  message("cannonical done")
  
  # main genes body ht
  d1 = mat_exp_rna %>% as_tibble() %>% dplyr::mutate(genes = rownames(mat_exp_rna))
  d2 = d1 %>% dplyr::filter(genes %in% genes_df)
  mat_exp_fil = data.matrix(d2)
  mat_exp_fil = suppressWarnings(mat_exp_fil[,-ncol(mat_exp_fil)])
  rownames(mat_exp_fil) = d2$genes
  mat_exp_fil = mat_exp_fil[genes_df,] #match the order of orginal genes
  message("main done")
  
  # Annotations -------
  # get df for annotation
  coldata = data.frame(Cluster = seu$minor_group_factor, #minor_group_factor
                       # Protocol = seu@meta.data$exp_proc,
                       # Menopause = seu@meta.data$menopause.y,
                       # Paridy = seu@meta.data$parity.1.y,
                       # Procedure = seu@meta.data$figure1_grp_updated,
                       sample_id = colnames(mat_exp_fil),
                       stringsAsFactors = F) %>% dplyr::arrange(Cluster)
  
  # top annotation
  ha_top = dplyr::select(coldata, Cluster) %>% 
    HeatmapAnnotation(df = ., show_annotation_name = FALSE, 
                      col = list(Cluster = col_grp))
  # Procedure = col_procedure,
  # Menopause = col_meno,
  # Paridy = col_paridy,
  # Protocol = col_Protocol))
  
  # get a sorted list of samples
  samps_sorted = coldata$sample_id
  
  # rearrnge the matrix
  mat = mat_exp_fil[, samps_sorted]
  mat0 = mat_exp_fil_cann[, samps_sorted]
  
  # # heatmap of exp
  # ht_exp = Heatmap(data.matrix(mat),
  #                  column_title = "Scaled Expression",
  #                  column_title_gp = gpar(fontsize = 6),
  #                  name = "Scaled Expression",
  #                  top_annotation = ha_top,show_row_names = T,
  #                  show_column_names = F, cluster_rows = F,
  #                  use_raster = F, cluster_columns = FALSE,
  #                  col = colorRamp2(c(-1.5, 0, 1.5), c("steelblue","white", "tomato3")),
  #                  row_names_gp = gpar(fontsize = 8),  width = unit(8, "cm"))
  
  # Mean cann ht ------
  # get the mean of exp 
  scores_cann = cbind(coldata, t(mat0))#ma
  scores_cann = scores_cann[,colnames(unique(as.matrix(scores_cann),MARGIN=2))] 
  scores_m_cann = suppressWarnings(scores_cann %>% group_by(Cluster) %>% unique() %>% summarise_all(funs(mean)))
  mat_exp_mean_cann = data.matrix(t(scores_m_cann))
  mat_exp_mean_cann = mat_exp_mean_cann[genes_df_cann,]
  #mat_exp_mean = mat_exp_mean[-c(1,2), ]
  colnames(mat_exp_mean_cann) = levels(scores_cann$Cluster)
  mat_exp_mean_cann = apply(mat_exp_mean_cann, 2,as.numeric)
  rownames(mat_exp_mean_cann) = genes_df_cann
  
  # Mean ht ------
  # get the mean of exp 
  scores = cbind(coldata, t(mat))#ma
  #scores = t(mat)#ma
  scores = scores[,colnames(unique(as.matrix(scores),MARGIN=2))] 
  scores_m = suppressWarnings(scores %>% group_by(Cluster) %>% unique() %>% summarise_all(funs(mean)))
  mat_exp_mean = data.matrix(t(scores_m))
  mat_exp_mean = mat_exp_mean[genes_df,]
  #mat_exp_mean = mat_exp_mean[-c(1,2), ]
  colnames(mat_exp_mean) = levels(scores$Cluster)
  mat_exp_mean = apply(mat_exp_mean, 2,as.numeric)
  rownames(mat_exp_mean) = genes_df
  
  # mean ht row annotation ------
  df_genes_sel = df %>% dplyr::filter(gene %in% rownames(mat_exp_mean)) %>% dplyr::select(gene, cluster) %>% mutate(id = paste0(gene, "-", cluster))
  
  ha_row = df_genes_sel %>% dplyr::select(cluster) %>% 
    rowAnnotation(df = ., show_annotation_name = FALSE, width =  unit(0.05, "mm"),
                  col = list(cluster = col_grp))
  
  # heatmap of exp
  ha = HeatmapAnnotation(#Menopause = anno_barplot(t(df_meta_meno2), border = FALSE, gp = gpar(fill = col_meno, border = NA, lty = "blank"), axis = FALSE, size = unit(0.85, "mm"), bar_width = 0.75),
    #                        Paridy = anno_barplot(t(df_meta_paridy2), border = FALSE, gp = gpar(fill = col_paridy, border = NA, lty = "blank"), axis = FALSE,bar_width = 0.85),
    #                        Procedure = anno_barplot(t(df_meta_grp2), border = FALSE, gp = gpar(fill = col_procedure, border = NA, lty = "blank"),axis = FALSE,bar_width = 0.85), 
    #                        Protocol = anno_barplot(t(df_meta_dig2), border = FALSE, gp = gpar(fill = col_Protocol, border = NA, lty = "blank"),axis = FALSE,bar_width = 0.85),
    #pt = anno_points(t(mat_exp_mean_cann)[,1], pch = 21, col = c("red", "blue")),
    ggplot = anno_empty(height = unit(1, "cm")),show_annotation_name = FALSE)
  #Cluster = levels(scores$Cluster), col = list(Cluster = col_grp), )
  #ha_row = HeatmapAnnotation(Cluster = levels(scores$Cluster), col = list(Cluster = col), show_annotation_name = TRUE)
  
  #ha1 = HeatmapAnnotation(Cluster = levels(scores$Cluster), col = list(Cluster = col_grp), show_annotation_name = FALSE)
  
  df_exp_mean_cann = mat_exp_mean_cann %>% as.data.frame()
  df_exp_mean_cann$gene = rownames(df_exp_mean_cann)
  df_exp_mean_cann1= reshape::melt(df_exp_mean_cann, id.vars=c("gene"))
  df_exp_mean_cann1$cluster = df_exp_mean_cann1$variable
  g<-ggplot(df_exp_mean_cann1, aes(x=cluster, y=gene, color=value, size=value)) + scale_radius(range=c(0, 3)) +
    geom_point() + theme_void() + theme(legend.position="none") + scale_color_gradient2(low = "steelblue", mid = "white", high = "tomato3")  
  
  
  
  ht_mean_cann = Heatmap(data.matrix(mat_exp_mean_cann),
                         #column_title = "Avg Exp",column_title_gp = gpar(fontsize = 6),
                         name = "Avg Expression",
                         top_annotation = ha1,
                         
                         show_row_names = T, row_names_side = c("right"),
                         show_column_names = F, show_row_dend = F, cluster_rows = F,
                         use_raster = T, cluster_columns = FALSE,
                         col = colorRamp2(c(-1, 0, 1), c("steelblue","white", "tomato3")),
                         #col = colorRamp2(c(-1, 0, 2), c("cadetblue4","white", "coral3")),
                         #col = colorRamp2(c(-1.5, 0, 1.5), c("black","black", "grey99")),
                         row_names_gp = gpar(fontsize = 5),  width = unit(2.5, "cm"), height = unit(0.75, "cm"))
  ht_mean_cann
  
  ht_mean = Heatmap(data.matrix(mat_exp_mean),
                    #column_title = "Avg Exp",column_title_gp = gpar(fontsize = 6),
                    name = "Avg Expression",
                    top_annotation = ha1, #ha for endothelial
                    left_annotation = ha_row,
                    show_row_names = T, row_names_side = c("right"), 
                    show_column_names = F, show_row_dend = F, cluster_rows = F,
                    use_raster = T, cluster_columns = FALSE,
                    #heatmap_legend_param = list(legend_direction = "horizontal",
                      #                          legend_width = unit(3, "cm"), 
                      #                          title_position = "lefttop"),
                    col = colorRamp2(c(-1, 0, 1), c("steelblue","white", "tomato3")),
                    #col = colorRamp2(c(-1, 0, 2), c("cadetblue4","white", "coral3")),
                    #col = colorRamp2(c(-1.5, 0, 1.5), c("black","black", "grey99")),
                    row_names_gp = gpar(fontsize = 5),  width = unit(2.5, "cm"), height = unit(4, "cm"))
  ht_mean
  message("ht means done")
  
  # png("htmaps_mean_no_bars_canno_integrated_dot_leftbar.png", res = 300, units = "in", width = 5, height = 5.5)
  # ht = draw(ha1 %v% ht_mean, heatmap_legend_side = "bottom")
  # decorate_annotation("ggplot", {
  #   vp = current.viewport()$name
  #   print(g, vp = vp)
  #   
  # })
  # dev.off()
  # 
  #ht = draw(ht_mean_cann %v% ht_mean, heatmap_legend_side = "bottom")
  ht = draw(ht_mean, heatmap_legend_side = "bottom")
  png("htmaps_mean_v1.png", res = 300, units = "in", width = 4, height = 4)
  #ht = draw(ha1 %v% ht_mean, heatmap_legend_side = "bottom")
  print(ht)
  dev.off()
  
  
  
  # pdf("htmaps_mean_small_bars.pdf", width = 8, height = 7)
  # print(ht)
  # dev.off()
  
}

run_mean_htmaps_fig_fibro <- function(df, seu, cell, genes_df, genes_df_cann, clus_levels, genes_sel_df) {
  
  # exp scores -------
  mat_exp = seu@assays$integrated@scale.data
  mat_exp = data.matrix(mat_exp)
  message("exp scores done")

  # exp scores -------
  mat_exp_rna = seu@assays$RNA@scale.data
  mat_exp_rna = data.matrix(mat_exp_rna)
  message("exp scores done")

  # cannonical genes
  # d1_cann = mat_exp_rna %>% as_tibble() %>% dplyr::mutate(genes = rownames(mat_exp_rna))
  # d2_cann = d1_cann %>% dplyr::filter(genes %in% genes_df_cann)
  # mat_exp_fil_cann = data.matrix(d2_cann)
  # mat_exp_fil_cann = suppressWarnings(mat_exp_fil_cann[,-ncol(mat_exp_fil_cann)])
  # rownames(mat_exp_fil_cann) = d2_cann$genes
  # mat_exp_fil_cann = mat_exp_fil_cann[genes_df_cann,] #match the order of orginal genes
  # message("cannonical done")
  # 
  
  # main genes body ht
  d1 = mat_exp_rna %>% as_tibble() %>% dplyr::mutate(genes = rownames(mat_exp_rna))#_rna
  d2 = d1 %>% dplyr::filter(genes %in% genes_df)
  mat_exp_fil = data.matrix(d2)
  mat_exp_fil = suppressWarnings(mat_exp_fil[,-ncol(mat_exp_fil)])
  rownames(mat_exp_fil) = d2$genes
  mat_exp_fil = mat_exp_fil[genes_df,] #match the order of orginal genes
  message("main done")
  
  # Annotations -------
  # get df for annotation
  coldata = data.frame(Cluster = seu$minor_group_factor, #minor_group_factor
                       # Protocol = seu@meta.data$exp_proc,
                       # Menopause = seu@meta.data$menopause.y,
                       # Paridy = seu@meta.data$parity.1.y,
                       # Procedure = seu@meta.data$figure1_grp_updated,
                       sample_id = colnames(mat_exp_fil),
                       stringsAsFactors = F) %>% dplyr::arrange(Cluster)
  
  # top annotation
  ha_top = dplyr::select(coldata, Cluster) %>% 
    HeatmapAnnotation(df = ., show_annotation_name = FALSE, 
                      col = list(Cluster = col_grp))
  # Procedure = col_procedure,
  # Menopause = col_meno,
  # Paridy = col_paridy,
  # Protocol = col_Protocol))
  
  # get a sorted list of samples
  samps_sorted = coldata$sample_id
  
  # rearrnge the matrix
  mat = mat_exp_fil[, samps_sorted]
  #mat0 = mat_exp_fil_cann[, samps_sorted]
  
  # # heatmap of exp
  # ht_exp = Heatmap(data.matrix(mat),
  #                  column_title = "Scaled Expression",
  #                  column_title_gp = gpar(fontsize = 6),
  #                  name = "Scaled Expression",
  #                  top_annotation = ha_top,show_row_names = T,
  #                  show_column_names = F, cluster_rows = F,
  #                  use_raster = F, cluster_columns = FALSE,
  #                  col = colorRamp2(c(-1.5, 0, 1.5), c("steelblue","white", "tomato3")),
  #                  row_names_gp = gpar(fontsize = 8),  width = unit(8, "cm"))
  
  # Mean cann ht ------
  # get the mean of exp 
  # scores_cann = cbind(coldata, t(mat0))#ma
  # scores_cann = scores_cann[,colnames(unique(as.matrix(scores_cann),MARGIN=2))] 
  # scores_m_cann = suppressWarnings(scores_cann %>% group_by(Cluster) %>% unique() %>% summarise_all(funs(mean)))
  # mat_exp_mean_cann = data.matrix(t(scores_m_cann))
  # mat_exp_mean_cann = mat_exp_mean_cann[genes_df_cann,]
  # colnames(mat_exp_mean_cann) = levels(scores_cann$Cluster)
  # mat_exp_mean_cann = apply(mat_exp_mean_cann, 2,as.numeric)
  # rownames(mat_exp_mean_cann) = genes_df_cann
  # 
  # Mean ht ------
  # get the mean of exp 
  scores = cbind(coldata, t(mat))#ma
  #scores = t(mat)#ma
  scores = scores[,colnames(unique(as.matrix(scores),MARGIN=2))] 
  scores_m = suppressWarnings(scores %>% group_by(Cluster) %>% unique() %>% summarise_all(funs(mean)))
  mat_exp_mean = data.matrix(t(scores_m))
  mat_exp_mean = mat_exp_mean[genes_df,]
  colnames(mat_exp_mean) = levels(scores$Cluster)
  mat_exp_mean = apply(mat_exp_mean, 2,as.numeric)
  rownames(mat_exp_mean) = genes_df
  
  # mean ht row annotation ------
  df_genes_sel = df %>% dplyr::filter(gene %in% rownames(mat_exp_mean)) %>% dplyr::select(gene, cluster) %>% mutate(id = paste0(gene, "-", cluster))
  #df_genes_sel = df_genes_sel[-c(29, 33, 54,32,61,36,38,35,34,37,59,57,58,56,60),] 
  df_genes_sel = df_genes_sel[-c(11, 15, 18, 17, 26),] 
  
  ha_row = df_genes_sel %>% dplyr::select(cluster) %>% 
    rowAnnotation(df = ., show_annotation_name = FALSE, width =  unit(0.05, "mm"),
                  col = list(cluster = col_grp))
  
  # heatmap of exp
  #ha = HeatmapAnnotation(#Menopause = anno_barplot(t(df_meta_meno2), border = FALSE, gp = gpar(fill = col_meno, border = NA, lty = "blank"), axis = FALSE, size = unit(0.85, "mm"), bar_width = 0.75),
    #                        Paridy = anno_barplot(t(df_meta_paridy2), border = FALSE, gp = gpar(fill = col_paridy, border = NA, lty = "blank"), axis = FALSE,bar_width = 0.85),
    #                        Procedure = anno_barplot(t(df_meta_grp2), border = FALSE, gp = gpar(fill = col_procedure, border = NA, lty = "blank"),axis = FALSE,bar_width = 0.85), 
    #                        Protocol = anno_barplot(t(df_meta_dig2), border = FALSE, gp = gpar(fill = col_Protocol, border = NA, lty = "blank"),axis = FALSE,bar_width = 0.85),
    #pt = anno_points(t(mat_exp_mean_cann)[,1], pch = 21, col = c("red", "blue")),
   # ggplot = anno_empty(height = unit(1, "cm")),show_annotation_name = FALSE)
  #Cluster = levels(scores$Cluster), col = list(Cluster = col_grp), )
  #ha_row = HeatmapAnnotation(Cluster = levels(scores$Cluster), col = list(Cluster = col), show_annotation_name = TRUE)
  
  ha1 = HeatmapAnnotation(Cluster = levels(scores$Cluster), col = list(Cluster = col_grp), show_annotation_name = FALSE)
  
  # df_exp_mean_cann = mat_exp_mean_cann %>% as.data.frame()
  # df_exp_mean_cann$gene = rownames(df_exp_mean_cann)
  # df_exp_mean_cann1= reshape::melt(df_exp_mean_cann, id.vars=c("gene"))
  # df_exp_mean_cann1$cluster = df_exp_mean_cann1$variable
  # g<-ggplot(df_exp_mean_cann1, aes(x=cluster, y=gene, color=value, size=value)) + scale_radius(range=c(0, 3)) +
  #   geom_point() + theme_void() + theme(legend.position="none") + scale_color_gradient2(low = "steelblue", mid = "white", high = "tomato3")  
  # 
  # 
  # 
  # ht_mean_cann = Heatmap(data.matrix(mat_exp_mean_cann),
  #                        #column_title = "Avg Exp",column_title_gp = gpar(fontsize = 6),
  #                        name = "Avg Expression",
  #                        top_annotation = ha1,
  #                        
  #                        show_row_names = T, row_names_side = c("right"),
  #                        show_column_names = F, show_row_dend = F, cluster_rows = F,
  #                        use_raster = T, cluster_columns = FALSE,
  #                        col = colorRamp2(c(-1, 0, 1), c("steelblue","white", "tomato3")),
  #                        #col = colorRamp2(c(-1, 0, 2), c("cadetblue4","white", "coral3")),
  #                        #col = colorRamp2(c(-1.5, 0, 1.5), c("black","black", "grey99")),
  #                        row_names_gp = gpar(fontsize = 5),  width = unit(2.5, "cm"), height = unit(0.75, "cm"))
  # ht_mean_cann
  # 
  ht_mean = Heatmap(data.matrix(mat_exp_mean),
                    #column_title = "Avg Exp",column_title_gp = gpar(fontsize = 6),
                    name = "Avg Expression",
                    top_annotation = ha1, #ha for endothelial
                    left_annotation = ha_row,
                    show_row_names = T, row_names_side = c("right"), 
                    show_column_names = F, show_row_dend = F, cluster_rows = F,
                    use_raster = T, cluster_columns = FALSE,
                    #heatmap_legend_param = list(legend_direction = "horizontal",
                    #                          legend_width = unit(3, "cm"), 
                    #                          title_position = "lefttop"),
                    col = colorRamp2(c(-0.5, 0, 0.5), c("steelblue","white", "tomato3")),
                    #col = colorRamp2(c(-1, 0, 2), c("cadetblue4","white", "coral3")),
                    #col = colorRamp2(c(-1.5, 0, 1.5), c("black","black", "grey99")),
                    row_names_gp = gpar(fontsize = 5),  width = unit(2.5, "cm"), height = unit(4, "cm"))
  ht_mean
  message("ht means done")
  
  # png("htmaps_mean_no_bars_canno_integrated_dot_leftbar.png", res = 300, units = "in", width = 5, height = 5.5)
  # ht = draw(ha1 %v% ht_mean, heatmap_legend_side = "bottom")
  # decorate_annotation("ggplot", {
  #   vp = current.viewport()$name
  #   print(g, vp = vp)
  #   
  # })
  # dev.off()
  # 
  #ht = draw(ht_mean_cann %v% ht_mean, heatmap_legend_side = "bottom")
  ht = draw(ht_mean, heatmap_legend_side = "bottom")
  png("htmaps_mean_v1_RNa_kevinshortlist.png", res = 300, units = "in", width = 4, height = 4)
  #ht = draw(ha1 %v% ht_mean, heatmap_legend_side = "bottom")
  print(ht)
  dev.off()
  
  
  
  # pdf("htmaps_mean_small_bars.pdf", width = 8, height = 7)
  # print(ht)
  # dev.off()
  
}

run_mean_htmaps_fig_epi <- function(df, seu, cell, genes_df, genes_df_cann, clus_levels, genes_sel_df) {
  
  # exp scores -------
  mat_exp = seu@assays$integrated@scale.data
  mat_exp = data.matrix(mat_exp)
  message("exp scores done")
  
  # exp scores -------
  mat_exp_rna = seu@assays$RNA@scale.data
  mat_exp_rna = data.matrix(mat_exp_rna)
  message("exp scores done")
  
  # cannonical genes
  # d1_cann = mat_exp_rna %>% as_tibble() %>% dplyr::mutate(genes = rownames(mat_exp_rna))
  # d2_cann = d1_cann %>% dplyr::filter(genes %in% genes_df_cann)
  # mat_exp_fil_cann = data.matrix(d2_cann)
  # mat_exp_fil_cann = suppressWarnings(mat_exp_fil_cann[,-ncol(mat_exp_fil_cann)])
  # rownames(mat_exp_fil_cann) = d2_cann$genes
  # mat_exp_fil_cann = mat_exp_fil_cann[genes_df_cann,] #match the order of orginal genes
  # message("cannonical done")
  # 
  
  # main genes body ht
  d1 = mat_exp_rna %>% as_tibble() %>% dplyr::mutate(genes = rownames(mat_exp_rna))#_rna
  #d1 = mat_exp %>% as_tibble() %>% dplyr::mutate(genes = rownames(mat_exp))#_rna
  d2 = d1 %>% dplyr::filter(genes %in% genes_df)
  mat_exp_fil = data.matrix(d2)
  mat_exp_fil = suppressWarnings(mat_exp_fil[,-ncol(mat_exp_fil)])
  rownames(mat_exp_fil) = d2$genes
  mat_exp_fil = mat_exp_fil[genes_df,] #match the order of orginal genes
  message("main done")
  
  # Annotations -------
  # get df for annotation
  coldata = data.frame(Cluster = seu$minor_group_factor, #minor_group_factor
                       # Protocol = seu@meta.data$exp_proc,
                       # Menopause = seu@meta.data$menopause.y,
                       # Paridy = seu@meta.data$parity.1.y,
                       # Procedure = seu@meta.data$figure1_grp_updated,
                       sample_id = colnames(mat_exp_fil),
                       stringsAsFactors = F) %>% dplyr::arrange(Cluster)
  
  # top annotation
  ha_top = dplyr::select(coldata, Cluster) %>% 
    HeatmapAnnotation(df = ., show_annotation_name = FALSE, 
                      col = list(Cluster = col_grp))
  # Procedure = col_procedure,
  # Menopause = col_meno,
  # Paridy = col_paridy,
  # Protocol = col_Protocol))
  
  # get a sorted list of samples
  samps_sorted = coldata$sample_id
  
  # rearrnge the matrix
  mat = mat_exp_fil[, samps_sorted]
  #mat0 = mat_exp_fil_cann[, samps_sorted]
  
  # # heatmap of exp
  # ht_exp = Heatmap(data.matrix(mat),
  #                  column_title = "Scaled Expression",
  #                  column_title_gp = gpar(fontsize = 6),
  #                  name = "Scaled Expression",
  #                  top_annotation = ha_top,show_row_names = T,
  #                  show_column_names = F, cluster_rows = F,
  #                  use_raster = F, cluster_columns = FALSE,
  #                  col = colorRamp2(c(-1.5, 0, 1.5), c("steelblue","white", "tomato3")),
  #                  row_names_gp = gpar(fontsize = 8),  width = unit(8, "cm"))
  
  # Mean cann ht ------
  # get the mean of exp 
  # scores_cann = cbind(coldata, t(mat0))#ma
  # scores_cann = scores_cann[,colnames(unique(as.matrix(scores_cann),MARGIN=2))] 
  # scores_m_cann = suppressWarnings(scores_cann %>% group_by(Cluster) %>% unique() %>% summarise_all(funs(mean)))
  # mat_exp_mean_cann = data.matrix(t(scores_m_cann))
  # mat_exp_mean_cann = mat_exp_mean_cann[genes_df_cann,]
  # colnames(mat_exp_mean_cann) = levels(scores_cann$Cluster)
  # mat_exp_mean_cann = apply(mat_exp_mean_cann, 2,as.numeric)
  # rownames(mat_exp_mean_cann) = genes_df_cann
  # 
  # Mean ht ------
  # get the mean of exp 
  scores = cbind(coldata, t(mat))#ma
  #scores = t(mat)#ma
  scores = scores[,colnames(unique(as.matrix(scores),MARGIN=2))] 
  scores_m = suppressWarnings(scores %>% group_by(Cluster) %>% unique() %>% summarise_all(funs(mean)))
  mat_exp_mean = data.matrix(t(scores_m))
  mat_exp_mean = mat_exp_mean[genes_df,]
  colnames(mat_exp_mean) = levels(scores$Cluster)
  mat_exp_mean = apply(mat_exp_mean, 2,as.numeric)
  rownames(mat_exp_mean) = genes_df
  
  # mean ht row annotation ------
  df_genes_sel = df %>% dplyr::filter(gene %in% rownames(mat_exp_mean)) %>% dplyr::select(gene, cluster) %>% mutate(id = paste0(gene, "-", cluster))
  #df_genes_sel = df_genes_sel[-c(8, 15, 24),] 
  
  ha_row = df_genes_sel %>% dplyr::select(cluster) %>% 
    rowAnnotation(df = ., show_annotation_name = FALSE, width =  unit(0.05, "mm"),
                  col = list(cluster = col_grp))
  
  # heatmap of exp
  #ha = HeatmapAnnotation(#Menopause = anno_barplot(t(df_meta_meno2), border = FALSE, gp = gpar(fill = col_meno, border = NA, lty = "blank"), axis = FALSE, size = unit(0.85, "mm"), bar_width = 0.75),
  #                        Paridy = anno_barplot(t(df_meta_paridy2), border = FALSE, gp = gpar(fill = col_paridy, border = NA, lty = "blank"), axis = FALSE,bar_width = 0.85),
  #                        Procedure = anno_barplot(t(df_meta_grp2), border = FALSE, gp = gpar(fill = col_procedure, border = NA, lty = "blank"),axis = FALSE,bar_width = 0.85), 
  #                        Protocol = anno_barplot(t(df_meta_dig2), border = FALSE, gp = gpar(fill = col_Protocol, border = NA, lty = "blank"),axis = FALSE,bar_width = 0.85),
  #pt = anno_points(t(mat_exp_mean_cann)[,1], pch = 21, col = c("red", "blue")),
  # ggplot = anno_empty(height = unit(1, "cm")),show_annotation_name = FALSE)
  #Cluster = levels(scores$Cluster), col = list(Cluster = col_grp), )
  #ha_row = HeatmapAnnotation(Cluster = levels(scores$Cluster), col = list(Cluster = col), show_annotation_name = TRUE)
  
  ha1 = HeatmapAnnotation(Cluster = levels(scores$Cluster), col = list(Cluster = col_grp), show_annotation_name = FALSE)
  
  # df_exp_mean_cann = mat_exp_mean_cann %>% as.data.frame()
  # df_exp_mean_cann$gene = rownames(df_exp_mean_cann)
  # df_exp_mean_cann1= reshape::melt(df_exp_mean_cann, id.vars=c("gene"))
  # df_exp_mean_cann1$cluster = df_exp_mean_cann1$variable
  # g<-ggplot(df_exp_mean_cann1, aes(x=cluster, y=gene, color=value, size=value)) + scale_radius(range=c(0, 3)) +
  #   geom_point() + theme_void() + theme(legend.position="none") + scale_color_gradient2(low = "steelblue", mid = "white", high = "tomato3")  
  # 
  # 
  # 
  # ht_mean_cann = Heatmap(data.matrix(mat_exp_mean_cann),
  #                        #column_title = "Avg Exp",column_title_gp = gpar(fontsize = 6),
  #                        name = "Avg Expression",
  #                        top_annotation = ha1,
  #                        
  #                        show_row_names = T, row_names_side = c("right"),
  #                        show_column_names = F, show_row_dend = F, cluster_rows = F,
  #                        use_raster = T, cluster_columns = FALSE,
  #                        col = colorRamp2(c(-1, 0, 1), c("steelblue","white", "tomato3")),
  #                        #col = colorRamp2(c(-1, 0, 2), c("cadetblue4","white", "coral3")),
  #                        #col = colorRamp2(c(-1.5, 0, 1.5), c("black","black", "grey99")),
  #                        row_names_gp = gpar(fontsize = 5),  width = unit(2.5, "cm"), height = unit(0.75, "cm"))
  # ht_mean_cann
  # 
  ht_mean = Heatmap(data.matrix(mat_exp_mean),
                    #column_title = "Avg Exp",column_title_gp = gpar(fontsize = 6),
                    name = "Avg Expression",
                    top_annotation = ha1, #ha for endothelial
                    #left_annotation = ha_row,
                    show_row_names = T, row_names_side = c("right"), 
                    show_column_names = F, show_row_dend = F, cluster_rows = T,
                    use_raster = T, cluster_columns = FALSE,
                    #heatmap_legend_param = list(legend_direction = "horizontal",
                    #                          legend_width = unit(3, "cm"), 
                    #                          title_position = "lefttop"),
                    col = colorRamp2(c(-1, 0, 1), c("steelblue","white", "tomato3")),
                    #col = colorRamp2(c(-1, 0, 2), c("cadetblue4","white", "coral3")),
                    #col = colorRamp2(c(-1.5, 0, 1.5), c("black","black", "grey99")),
                    row_names_gp = gpar(fontsize = 5),  width = unit(2.5, "cm"), height = unit(4, "cm"))
  ht_mean
  message("ht means done")
  
  # png("htmaps_mean_no_bars_canno_integrated_dot_leftbar.png", res = 300, units = "in", width = 5, height = 5.5)
  # ht = draw(ha1 %v% ht_mean, heatmap_legend_side = "bottom")
  # decorate_annotation("ggplot", {
  #   vp = current.viewport()$name
  #   print(g, vp = vp)
  #   
  # })
  # dev.off()
  # 
  #ht = draw(ht_mean_cann %v% ht_mean, heatmap_legend_side = "bottom")
  ht = draw(ht_mean, heatmap_legend_side = "bottom")
  #png("htmaps_mean_v1_basal.png", res = 300, units = "in", width = 4, height = 4)
  png("htmaps_mean_v1_lumhr_keratins.png", res = 300, units = "in", width = 4, height = 4)
  #ht = draw(ha1 %v% ht_mean, heatmap_legend_side = "bottom")
  print(ht)
  dev.off()
  
  
  
  # pdf("htmaps_mean_small_bars.pdf", width = 8, height = 7)
  # print(ht)
  # dev.off()
  
}


# ***Pathway FGSEA ----------------------
#BiocManager::install("GSVAdata")
#BiocManager::install("fgsea")
#BiocManager::install("GSVA")
#BiocManager::install("GO.db")
#install.packages("msigdbr")
#install.packages("virdis") #BiocManager::install("goseq")
library(GSVA);library(GSVAdata);library(GO.db); library(msigdbr)
library(fgsea)
#library(fgsea,lib.loc = "/volumes/core/users/jerome/R/x86_64-redhat-linux-gnu-library/3.6/")
library(viridis)
library(GSEABase);library(org.Hs.eg.db);library(goseq)
data(c2BroadSets)
data(leukemia)
library(stringr)
fun_gsea = function(genes, fc = 0.10, set = "C2", min_cat=5, sub_cat = NULL, num_avg_logFC = 0.25, doplot = TRUE, ret_collap = TRUE, celltype = "lumsec"){
  
  out = list()
  out_collapsed = list()
  
  for(i in 1:length(unique(genes$cluster))){
    
    #message("running on cluster - ", i)
    j = unique(genes$cluster)[i]
    message("running on cluster - ", j)
    
    final_genes = genes %>% dplyr::filter(cluster == j) %>% 
      #dplyr::filter(as.numeric(p_val_adj) <= 0.05, avg_logFC > num_avg_logFC) %>% 
      dplyr::arrange(-avg_logFC) %>% dplyr::select(gene, avg_logFC)
    
    gene_vec = dplyr::pull(final_genes, var = avg_logFC)
    names(gene_vec) = final_genes$gene
    
    if(set %in% c("H","C6","C7")){
      msig_df = msigdbr(species = "Homo sapiens", category = set, subcategory = NULL)
    }else{
      msig_df = msigdbr(species = "Homo sapiens", category = set, subcategory = sub_cat)
    }
    msig_list = msig_df %>% split(x = .$gene_symbol, f = .$gs_name)
    fgsea_run = fgsea(pathways = msig_list, stats = gene_vec, nperm = 1000)
    
    if (doplot == TRUE) {
      
      for (k in 1:length(msig_list)) {
        gset_name = names(msig_list)[k]
        pp = plotEnrichment(msig_list[[gset_name]], + gene_vec) + labs(title=gset_name)
        
        #message("working on - ", k)
        #png(paste0(celltype, "_", gset_name, "_", j, ".png"), res = 300, units = "in", width = 5, height = 3.5)
        #print(pp)
        #dev.off()
      }  
    }
    
    
    fgsea_out = fgsea_run %>% dplyr::filter(size > min_cat) %>% dplyr::arrange(pval)
    fgsea_out_collapsed = fgsea::collapsePathways(fgsea_run, pathways = msig_list, 
                                                  stats = gene_vec,  nperm = 1000, pval.threshold = 0.05)
    
    # select the main pathways ---
    fgsea_out_collapsed =  fgsea_out[pathway %in% fgsea_out_collapsed$mainPathways,]
    
    out[[length(out)+1]] = dplyr::tibble(fgsea_out)
    names(out)[length(out)] = j
    
    out_collapsed[[length(out_collapsed)+1]] = dplyr::tibble(fgsea_out_collapsed)
    names(out_collapsed)[length(out_collapsed)] = j
  }
  
  res_out = do.call(rbind, out)
  res_out$cluster = str_split(row.names(res_out), "\\.", simplify = TRUE)[,1]
  
  res_out_collapsed = do.call(rbind, out_collapsed)
  res_out_collapsed$cluster = str_split(row.names(res_out_collapsed), "\\.", simplify = TRUE)[,1]
  
  
  write_rds(res_out, paste0(celltype, "_", set, "_", sub_cat, "_min_size_",min_cat,".rds"))
  write_rds(res_out_collapsed, paste0(celltype, "_", set, "_", sub_cat, "_min_size_",min_cat,"collapsed.rds"))
  
  if(ret_collap == TRUE){
    
    return(res_out_collapsed)
  }
  
  return(res_out) 
  
}
# example run
# gsea_resultsc2reac = fun_gsea(genes = markers, set = "H", sub_cat = "NULL", num_avg_logFC = 0.25)

fun_gsea_results = function(df_gsea = gsea_resultsc2, num_pval = 0.05, num_top = 20){
  
  char_genesets = df_gsea %>%    
    #mutate(treatment = str_sub(cluster,1,1)) %>% 
    arrange(cluster) %>% 
    group_by(cluster) %>% 
    filter(pval < num_pval) %>% 
    #top_n(num_top, -pval) %>% 
    pull(pathway) %>% 
    unique()
  return(char_genesets)
}

# example run
# char_genesets = fun_gsea_results(df_gsea = gsea_resultsc2, num_pval = 0.05, num_top = 10)

fun_gsea_df = function(df_gsea_out = gsea_resultsc2reac, char_genesets = char_genesets){
  
  df = df_gsea_out %>% 
    dplyr::filter(pathway %in% char_genesets) %>%
    #mutate(group = ("DCIS")) %>% 
    dplyr::select(pathway, NES, cluster) %>% 
    spread(key = cluster, value = NES, fill=0) %>% 
    gather(key = "cluster", value="NES", -pathway) %>% 
    mutate(pathway = factor(pathway, levels = char_genesets))
}

# example run
# df_htmap_input = fun_gsea_df(df_gsea_out = gsea_resultsc2, char_genesets = char_genesets)

# *** Get top genes for heatmap ---------



func_select_markers_thresh = function(all_markers_sel, pct_var = 0.25){
  
  # df with all values, filter based on pct 1,2 for all cluster
  #cls = "Fibroblasts"
  check = lapply(df_thresh$cluster, function(cls){
    
    all_markers_sel %>% dplyr::filter(cluster == cls) %>% 
      dplyr::arrange(pct.1) %>% 
      tidylog::filter(pct.1 > df_thresh[cls, "thresh_pct1"]) %>% 
      tidylog::filter(pct.2 < pct_var)
    
  }) %>% bind_rows()
  
  # filter(check, cluster1 == "Fibroblasts")
  
  
  df_top_genes_cluster = check %>% as_tibble() %>% 
    dplyr::select(gene, cluster_top = cluster) %>% 
    mutate(in_top = 1) %>% 
    left_join(all_markers_sel, by = "gene") %>% 
    filter(avg_logFC > 0 ) %>% 
    group_by(gene, cluster_top) %>% tally(name = "n_top_cluster") %>% 
    ungroup() %>% filter(n_top_cluster == 1)
  
  df_top_genes_cluster_key = df_top_genes_cluster %>% mutate(key = paste0(gene, "_", cluster_top))
  df_top_genes_cluster_sel = all_markers_sel %>% mutate(key = paste0(gene, "_", cluster)) %>% 
    dplyr::filter(key %in% df_top_genes_cluster_key$key, avg_logFC > 0)
  return(df_top_genes_cluster_sel)
  
} 



func_select_markers = function(all_markers_sel, pct_var = 0.25){
  
  all_markers_sel %>% 
    filter(pct.2 < 0.25) %>% 
    filter(pct.1 > pct_var) -> check # %>% filter(avg_logFC > 1.2) 
  
  df_top_genes_cluster = check %>% as_tibble() %>% 
    dplyr::select(gene, cluster_top = cluster1) %>% 
    mutate(in_top = 1) %>% 
    left_join(all_markers_sel, by = "gene") %>% 
    filter(avg_logFC > 0 ) %>% 
    group_by(gene, cluster_top) %>% tally(name = "n_top_cluster") %>% 
    ungroup() %>% filter(n_top_cluster == 1)
  
  df_top_genes_cluster_key = df_top_genes_cluster %>% mutate(key = paste0(gene, "_", cluster_top))
  df_top_genes_cluster_sel = all_markers_sel %>% mutate(key = paste0(gene, "_", cluster1)) %>% 
    dplyr::filter(key %in% df_top_genes_cluster_key$key, avg_logFC > 0)
  return(df_top_genes_cluster_sel)
  
} 

run_mean_htmaps_fig6_horizontal <- function(df, seu, cell, genes_df, genes_df_cann, clus_levels, genes_sel_df) {
  
  # exp scores -------
  mat_exp_rna = seu@assays$integrated@data
  mat_exp_rna = data.matrix(mat_exp_rna)
  mat_exp_rna =  expm1(mat_exp_rna)
  message("exp scores done")
  
  # main genes body ht
  d1 = mat_exp_rna %>% as_tibble() %>% dplyr::mutate(genes = rownames(mat_exp_rna))
  d2 = d1 %>% dplyr::filter(genes %in% genes_df)
  mat_exp_fil = data.matrix(d2)
  mat_exp_fil = suppressWarnings(mat_exp_fil[,-ncol(mat_exp_fil)])
  rownames(mat_exp_fil) = d2$genes
  mat_exp_fil = mat_exp_fil[genes_df,] #match the order of orginal genes
  message("main done")
  
  # Annotations -------
  # get df for annotation
  coldata = data.frame(Cluster = seu$final_group1, #minor_group_factor
                       # Protocol = seu@meta.data$exp_proc,
                       # Menopause = seu@meta.data$menopause.y,
                       # Paridy = seu@meta.data$parity.1.y,
                       # Procedure = seu@meta.data$figure1_grp_updated,
                       sample_id = colnames(mat_exp_fil),
                       stringsAsFactors = F) %>% dplyr::arrange(Cluster)
  
  # top annotation
  ha_top = dplyr::select(coldata, Cluster) %>% 
    HeatmapAnnotation(df = ., show_annotation_name = FALSE, 
                      col = list(Cluster = col_grp))
  
  
  # get a sorted list of samples
  samps_sorted = coldata$sample_id
  
  # rearrnge the matrix
  mat = mat_exp_fil[, samps_sorted]
  
  # Mean ht ------
  # get the mean of exp 
  scores = cbind(coldata, t(mat))#ma
  #scores = t(mat)#ma
  scores = scores[,colnames(unique(as.matrix(scores),MARGIN=2))] 
  scores_m = suppressWarnings(scores %>% group_by(Cluster) %>% unique() %>% summarise_all(funs(mean)))
  mat_exp_mean = data.matrix(t(scores_m))
  mat_exp_mean = mat_exp_mean[genes_df,]
  #mat_exp_mean = mat_exp_mean[-c(1,2), ]
  colnames(mat_exp_mean) = levels(scores$Cluster)
  mat_exp_mean = apply(mat_exp_mean, 2,as.numeric)
  rownames(mat_exp_mean) = genes_df
  
  # scale and center
  mat_exp_mean = ScaleData(mat_exp_mean)
  
  # mean ht row annotation ------
  df_genes_sel = df %>% dplyr::filter(gene %in% rownames(mat_exp_mean)) %>% dplyr::select(gene, cluster) %>% mutate(id = paste0(gene, "-", cluster))
  df_genes_sel = df_genes_sel %>% dplyr::filter(id %in% genes_sel_df$id)
  
  ha_row = df_genes_sel %>% dplyr::select(cluster) %>% 
    HeatmapAnnotation(df = ., show_annotation_name = FALSE, 
                      col = list(cluster = col_grp))
  
  ha1 = rowAnnotation(Cluster = levels(scores$Cluster), col = list(Cluster = col_grp), show_annotation_name = FALSE)
  
  #dim(t(data.matrix(mat_exp_mean)))
  ht_mean = Heatmap(t(data.matrix(mat_exp_mean)),
                    name = "Avg Expression",
                    # row params
                    #left_annotation = ha1, 
                    show_row_dend = F, cluster_rows = F,
                    show_row_names = F, #row_names_side = c("right"), 
                    # row params
                    top_annotation = ha_row,
                    cluster_columns = FALSE, 
                    column_names_gp = gpar(fontsize = 5),
                    column_names_rot = 65,
                    show_column_names = T, 
                    
                    use_raster = T, 
                    col = colorRamp2(c(-1, 0, 2), c("steelblue","white", "tomato3")),
                    width = unit(15, "cm"), height = unit(4, "cm"))
  ht_mean
  message("ht means done")
  
  ht = draw(ht_mean, heatmap_legend_side = "bottom")
  png("htmaps_mean_v1_horizontal.png", res = 300, units = "in", width = 12.5, height = 8.5)
  print(ht)
  dev.off()
  
  pdf("htmaps_mean_v1_horizontal.pdf", width = 12.5, height = 8)
  print(ht)
  dev.off()
  
  
}

library(ggpubr)
do_clinical_comp = function(df_mrg_nuc2 = df_mrg_nuc2){
  
  # ....... // patient BMI ----
  df_mrg_nuc_filt = df_mrg_nuc2 %>% distinct(patient_id, bmi2, exp_proc, exp_proc_updated,celltype, cellgroup, prop_cell_per_pt, prop_cell_per_grp_pt, cell_per_grp_pt, cell_per_pt, total_cells_per_pt)
  df_mrg_nuc_filt = df_mrg_nuc_filt %>% dplyr::filter(!bmi2 %in% c("unknown"))
  df_mrg_nuc_filt$bmi2 = factor(df_mrg_nuc_filt$bmi2, levels = c("normal", "overweight"), ordered = T) 
  
  my_comparisons <- list(c("normal", "overweight"))
  theme_sel = theme(legend.position = "none",
                    strip.background = element_blank(),
                    panel.spacing.x = unit(0,"line"),
                    legend.background = element_blank(),
                    legend.title = element_text(size=6),
                    legend.text = element_text(size=6),
                    legend.key = element_blank(),
                    panel.grid = element_blank(),
                    #panel.border = element_blank(),
                    panel.border = element_rect(colour = "grey", fill=NA, size=0.3),
                    axis.line = element_blank(),
                    panel.background = element_blank(),
                    axis.text = element_text(size = 10),
                    axis.text.x = element_text(angle = 90),
                    axis.ticks = element_blank(),
                    axis.title=element_text(size=10))
  
  stat.test = compare_means(formula = prop_cell_per_pt ~ bmi2, group.by = "celltype", data = df_mrg_nuc_filt %>% dplyr::filter(!bmi2 == 0), 
                            p.adjust.method = "fdr") %>% mutate(y_pos = max(df_mrg_nuc_filt$prop_cell_per_pt) + 0.05, p.adj = format.pval(p.adj, digits = 2))
  
  p0 = ggboxplot(df_mrg_nuc_filt, x = "bmi2", y = "prop_cell_per_pt",title = "prop_by_all_cell",
                 fill = "celltype", palette = col_pal,
                 add = "jitter", repel = T, outlier.shape = NA, size = 0.3, alpha = 0.9, add.params = list(size  = 0.5, alpha = 0.7)) +  theme_sel +
    facet_wrap(~celltype, ncol = 14) + geom_signif(
      data=stat.test, 
      aes(xmin = group1, xmax = group2, y_position = y_pos, annotations = p.signif), # 
      manual= TRUE);p0 #+ 
  
  png(paste0(save_path,"cell_proportions_bmi_paper_patient_v2.png"), width=8, height=5, units = "in", res = 300)
  print(p0)
  dev.off()
  
  pdf(paste0(save_path,"cell_proportions_bmi_paper_patient_v2.pdf"), width=6, height=4)
  print(p0)
  dev.off()
  
  # ....... // patient meno ----
  df_mrg_nuc_filt = df_mrg_nuc2 %>% distinct(patient_id, menopause, exp_proc, exp_proc_updated,celltype, cellgroup, prop_cell_per_pt, prop_cell_per_grp_pt, cell_per_grp_pt, cell_per_pt, total_cells_per_pt)
  df_mrg_nuc_filt = df_mrg_nuc_filt %>% dplyr::filter(!menopause %in% c("unknown"))
  df_mrg_nuc_filt$menopause = factor(df_mrg_nuc_filt$menopause, levels = c("pre", "post"), ordered = T) 
  
  my_comparisons <- list( c("pre", "post"))
  stat.test = compare_means(formula = prop_cell_per_pt ~ menopause, group.by = "celltype", data = df_mrg_nuc_filt %>% dplyr::filter(!menopause == 0), 
                            p.adjust.method = "fdr") %>% mutate(y_pos = max(df_mrg_nuc_filt$prop_cell_per_pt) + 0.05, p.adj = format.pval(p.adj, digits = 2))
  
  p0 = ggboxplot(df_mrg_nuc_filt, x = "menopause", y = "prop_cell_per_pt",title = "prop_by_all_cell",
                 fill = "celltype", palette = col_pal,
                 add = "jitter", outlier.shape = NA, size = 0.3, alpha = 0.9, add.params = list(size  = 0.5, alpha = 0.7)) +  theme_sel +
    facet_wrap(~celltype, ncol = 14) + geom_signif(
      data=stat.test, 
      aes(xmin = group1, xmax = group2, y_position = y_pos, annotations = p.signif), # 
      manual= TRUE);p0
  png(paste0(save_path,"cell_proportions_menopause_paper_patient.png"), width=8, height=5, units = "in", res = 300)
  print(p0)
  dev.off()
  
  pdf(paste0(save_path,"cell_proportions_menopause_paper_patient.pdf"), width=6, height=4)
  print(p0)
  dev.off()
  
  # ....... // patient age ----
  df_mrg_nuc_filt = df_mrg_nuc2 %>% distinct(patient_id, age2, exp_proc, exp_proc_updated,celltype, cellgroup, prop_cell_per_pt, prop_cell_per_grp_pt, cell_per_grp_pt, cell_per_pt, total_cells_per_pt)
  df_mrg_nuc_filt = df_mrg_nuc_filt %>% dplyr::filter(!age2 %in% c("unknown"), !is.na(age2))
  df_mrg_nuc_filt$age2 = factor(df_mrg_nuc_filt$age2, levels = c("young", "old"), ordered = T) 
  
  my_comparisons <- list(c("young", "old"))
  stat.test = compare_means(formula = prop_cell_per_pt ~ age2, group.by = "celltype", data = df_mrg_nuc_filt %>% dplyr::filter(!age2 == 0), 
                            p.adjust.method = "fdr") %>% mutate(y_pos = max(df_mrg_nuc_filt$prop_cell_per_pt) + 0.05, p.adj = format.pval(p.adj, digits = 2))
  
  p0 = ggboxplot(df_mrg_nuc_filt, x = "age2", y = "prop_cell_per_pt",title = "prop_by_all_cell",
                 fill = "celltype", palette = col_pal,
                 add = "jitter", outlier.shape = NA, size = 0.3, alpha = 0.9, add.params = list(size  = 0.5, alpha = 0.7)) +  theme_sel +
    facet_wrap(~celltype, ncol = 14) + geom_signif(
      data=stat.test, 
      aes(xmin = group1, xmax = group2, y_position = y_pos, annotations = p.signif), # 
      manual= TRUE);p0
  
  png(paste0(save_path,"cell_proportions_age2_paper_patient.png"), width=8, height=5, units = "in", res = 300)
  print(p0)
  dev.off()
  
  pdf(paste0(save_path,"cell_proportions_age2_paper_patient.pdf"), width=6, height=4)
  print(p0)
  dev.off() 
  
  # ....... // paridy ----
  df_mrg_nuc_filt = df_mrg_nuc2 %>% distinct(patient_id, parity, exp_proc, exp_proc_updated,celltype, cellgroup, prop_cell_per_pt, prop_cell_per_grp_pt, cell_per_grp_pt, cell_per_pt, total_cells_per_pt)
  df_mrg_nuc_filt = df_mrg_nuc_filt %>% dplyr::filter(!parity %in% c("unknown"), !is.na(parity))
  
  my_comparisons <- list( c("0", "1"))
  stat.test = compare_means(formula = prop_cell_per_pt ~ parity, group.by = "celltype", data = df_mrg_nuc_filt, 
                            p.adjust.method = "fdr") %>% mutate(y_pos = max(df_mrg_nuc_filt$prop_cell_per_pt) + 0.05, p.adj = format.pval(p.adj, digits = 2))
  
  p0 = ggboxplot(df_mrg_nuc_filt, x = "parity", y = "prop_cell_per_pt",title = "prop_by_all_cell",
                 fill = "celltype", palette = col_pal,
                 add = "jitter", outlier.shape = NA, size = 0.3, alpha = 0.9, add.params = list(size  = 0.5, alpha = 0.7)) +  theme_sel +
    facet_wrap(~celltype, ncol = 14) + geom_signif(
      data=stat.test, 
      aes(xmin = group1, xmax = group2, y_position = y_pos, annotations = p.signif), # 
      manual= TRUE);p0
  
  
  png(paste0(save_path,"cell_proportions_parity_paper_patient.png"), width=8, height=5, units = "in", res = 300)
  print(p0)
  dev.off()
  
  pdf(paste0(save_path,"cell_proportions_parity_paper_patient.pdf"), width=6, height=4)
  print(p0)
  dev.off() 
  
  # ....... // gravida ----
  df_mrg_nuc_filt = df_mrg_nuc2 %>% distinct(patient_id, gravida, exp_proc, exp_proc_updated,celltype, cellgroup, prop_cell_per_pt, prop_cell_per_grp_pt, cell_per_grp_pt, cell_per_pt, total_cells_per_pt)
  df_mrg_nuc_filt = df_mrg_nuc_filt %>% dplyr::filter(!gravida %in% c("unknown"), !is.na(gravida))
  
  my_comparisons <- list( c("0", "1"))
  stat.test = compare_means(formula = prop_cell_per_pt ~ gravida, group.by = "celltype", data = df_mrg_nuc_filt, 
                            p.adjust.method = "fdr") %>% mutate(y_pos = max(df_mrg_nuc_filt$prop_cell_per_pt) + 0.05, p.adj = format.pval(p.adj, digits = 2))
  
  p0 = ggboxplot(df_mrg_nuc_filt, x = "gravida", y = "prop_cell_per_pt",title = "prop_by_all_cell",
                 fill = "celltype", palette = col_pal,
                 add = "jitter", outlier.shape = NA, size = 0.3, alpha = 0.9, add.params = list(size  = 0.5, alpha = 0.7)) +  theme_sel +
    facet_wrap(~celltype, ncol = 14) + geom_signif(
      data=stat.test, 
      aes(xmin = group1, xmax = group2, y_position = y_pos, annotations = p.signif), # 
      manual= TRUE);p0
  
  png(paste0(save_path,"cell_proportions_gravida_paper_patient.png"), width=8, height=5, units = "in", res = 300)
  print(p0)
  dev.off()
  
  pdf(paste0(save_path,"cell_proportions_gravida_paper_patient.pdf"), width=6, height=4)
  print(p0)
  dev.off()
  
  
  # ....... // density ----
  df_mrg_nuc_filt = df_mrg_nuc2 %>% distinct(patient_id, density2, exp_proc, exp_proc_updated,celltype, cellgroup, prop_cell_per_pt, prop_cell_per_grp_pt, cell_per_grp_pt, cell_per_pt, total_cells_per_pt)
  df_mrg_nuc_filt = df_mrg_nuc_filt %>% dplyr::filter(!density2 %in% c("unknown"), !is.na(density2))
  df_mrg_nuc_filt$density2 = factor(df_mrg_nuc_filt$density2, levels = c("low", "high"), ordered = T)
  
  my_comparisons <- list( c("low", "high"))
  stat.test = compare_means(formula = prop_cell_per_pt ~ density2, group.by = "celltype", data = df_mrg_nuc_filt, 
                            p.adjust.method = "fdr") %>% mutate(y_pos = max(df_mrg_nuc_filt$prop_cell_per_pt) + 0.05, p.adj = format.pval(p.adj, digits = 2))
  
  p0 = ggboxplot(df_mrg_nuc_filt, x = "density2", y = "prop_cell_per_pt",title = "prop_by_all_cell",
                 fill = "celltype", palette = col_pal,
                 add = "jitter", outlier.shape = NA, size = 0.3, alpha = 0.9, add.params = list(size  = 0.5, alpha = 0.7)) +  theme_sel +
    facet_wrap(~celltype, ncol = 14) + geom_signif(
      data=stat.test, 
      aes(xmin = group1, xmax = group2, y_position = y_pos, annotations = p.signif), # 
      manual= TRUE);p0
  
  png(paste0(save_path,"cell_proportions_density_paper_patient.png"), width=8, height=5, units = "in", res = 300)
  print(p0)
  dev.off()
  
  pdf(paste0(save_path,"cell_proportions_density_paper_patient.pdf"), width=6, height=4)
  print(p0)
  dev.off()
  
  # ....... // germline_cancer_risk ----
  df_mrg_nuc_filt = df_mrg_nuc2 %>% distinct(patient_id, germline_cancer_risk, exp_proc, exp_proc_updated,celltype, cellgroup, prop_cell_per_pt, prop_cell_per_grp_pt, cell_per_grp_pt, cell_per_pt, total_cells_per_pt)
  df_mrg_nuc_filt = df_mrg_nuc_filt %>% dplyr::filter(!germline_cancer_risk %in% c("unknown"), !is.na(germline_cancer_risk))
  df_mrg_nuc_filt$germline_cancer_risk = factor(df_mrg_nuc_filt$germline_cancer_risk, levels = c("normal", "high"), ordered = T)
  
  my_comparisons <- list(c("normal", "high"))
  stat.test = compare_means(formula = prop_cell_per_pt ~ germline_cancer_risk, group.by = "celltype", data = df_mrg_nuc_filt, 
                            p.adjust.method = "fdr") %>% mutate(y_pos = max(df_mrg_nuc_filt$prop_cell_per_pt) + 0.05, p.adj = format.pval(p.adj, digits = 2))
  
  p0 = ggboxplot(df_mrg_nuc_filt, x = "germline_cancer_risk", y = "prop_cell_per_pt",title = "prop_by_all_cell",
                 fill = "celltype", palette = col_pal,
                 add = "jitter", outlier.shape = NA, size = 0.3, alpha = 0.9, add.params = list(size  = 0.5, alpha = 0.7)) +  theme_sel +
    facet_wrap(~celltype, ncol = 14) + geom_signif(
      data=stat.test, 
      aes(xmin = group1, xmax = group2, y_position = y_pos, annotations = p.signif), # 
      manual= TRUE);p0
  
  png(paste0(save_path,"cell_proportions_germline_cancer_risk_paper_patient.png"), width=8, height=5, units = "in", res = 300)
  print(p0)
  dev.off()
  
  pdf(paste0(save_path,"cell_proportions_germline_cancer_risk_paper_patient.pdf"), width=6, height=4)
  print(p0)
  dev.off()
  
  # ....... // surgery ----
  df_mrg_nuc_filt = df_mrg_nuc2 %>% distinct(patient_id, surgery, exp_proc, exp_proc_updated,celltype, cellgroup, prop_cell_per_pt, prop_cell_per_grp_pt, cell_per_grp_pt, cell_per_pt, total_cells_per_pt)
  df_mrg_nuc_filt = df_mrg_nuc_filt %>% dplyr::filter(!surgery %in% c("unknown"), !is.na(surgery))
  df_mrg_nuc_filt$surgery = factor(df_mrg_nuc_filt$surgery, levels = c("Reduction Mammoplasty", "Prophylatic Mastectomy", "Cancer Mastectomy"), ordered = T)
  
  my_comparisons <- list(c("Reduction Mammoplasty", "Prophylatic Mastectomy"), c("Prophylatic Mastectomy", "Cancer Mastectomy"), c("Reduction Mammoplasty", "Cancer Mastectomy"))
  stat.test = compare_means(formula = prop_cell_per_pt ~ surgery, group.by = "celltype", data = df_mrg_nuc_filt, 
                            p.adjust.method = "fdr") %>% mutate(y_pos = max(df_mrg_nuc_filt$prop_cell_per_pt) + 0.05, p.adj = format.pval(p.adj, digits = 2))
  
  p0 = ggboxplot(df_mrg_nuc_filt, x = "surgery", y = "prop_cell_per_pt",title = "prop_by_all_cell",
                 fill = "celltype", palette = col_pal,
                 add = "jitter", outlier.shape = NA, size = 0.3, alpha = 0.9, add.params = list(size  = 0.5, alpha = 0.7)) +  theme_sel +
    facet_wrap(~celltype, ncol = 14) + geom_signif(
      data=stat.test, 
      aes(xmin = group1, xmax = group2, y_position = y_pos, annotations = p.signif), # 
      manual= TRUE);p0
  
  png(paste0(save_path,"cell_proportions_surgery_paper_patient.png"), width=8, height=5, units = "in", res = 300)
  print(p0)
  dev.off()
  
  pdf(paste0(save_path,"cell_proportions_surgery_risk_paper_patient.pdf"), width=9, height=4)
  print(p0)
  dev.off()
  
  # ....... // exp proc ----
  df_mrg_nuc_filt = df_mrg_nuc2 %>% distinct(patient_id, surgery, exp_proc, exp_proc_updated,celltype, cellgroup, prop_cell_per_pt, prop_cell_per_grp_pt, cell_per_grp_pt, cell_per_pt, total_cells_per_pt)
  df_mrg_nuc_filt = df_mrg_nuc_filt %>% dplyr::filter(!exp_proc_updated %in% c("unknown"), !is.na(exp_proc_updated))
  df_mrg_nuc_filt$exp_proc_updated = factor(df_mrg_nuc_filt$exp_proc_updated, levels = c("short", "overnight"), ordered = T)
  
  my_comparisons <- list(c("short", "overnight"))
  stat.test = compare_means(formula = prop_cell_per_pt ~ exp_proc_updated, group.by = "celltype", data = df_mrg_nuc_filt, 
                            p.adjust.method = "fdr") %>% mutate(y_pos = max(df_mrg_nuc_filt$prop_cell_per_pt) + 0.05, p.adj = format.pval(p.adj, digits = 2))
  
  p0 = ggboxplot(df_mrg_nuc_filt, x = "exp_proc_updated", y = "prop_cell_per_pt",title = "prop_by_all_cell",
                 fill = "celltype", palette = col_pal,
                 add = "jitter", outlier.shape = NA, size = 0.3, alpha = 0.9, add.params = list(size  = 0.5, alpha = 0.7)) +  theme_sel +
    facet_wrap(~celltype, ncol = 14) + geom_signif(
      data=stat.test, 
      aes(xmin = group1, xmax = group2, y_position = y_pos, annotations = p.signif), # 
      manual= TRUE);p0
  
  png(paste0(save_path,"cell_proportions_exp_proc_updated_paper_patient.png"), width=8, height=5, units = "in", res = 300)
  print(p0)
  dev.off()
  
  pdf(paste0(save_path,"cell_proportions_exp_proc_updated_risk_paper_patient.pdf"), width=6, height=4)
  print(p0)
  dev.off()
  
  
  # ....... // ethnicity2 ----
  df_mrg_nuc_filt = df_mrg_nuc2 %>% distinct(patient_id, surgery, exp_proc, ethnicity2,celltype, cellgroup, prop_cell_per_pt, prop_cell_per_grp_pt, cell_per_grp_pt, cell_per_pt, total_cells_per_pt)
  df_mrg_nuc_filt = df_mrg_nuc_filt %>% dplyr::filter(!ethnicity2 %in% c("unknown"), !is.na(ethnicity2))
  df_mrg_nuc_filt$ethnicity2 = factor(df_mrg_nuc_filt$ethnicity2, levels = c("african-american","caucasian","latino"), ordered = T)
  
  my_comparisons <- list(c("african-american","caucasian"), c("caucasian","latino"), c("african-american", "latino"))
  stat.test = compare_means(formula = prop_cell_per_pt ~ ethnicity2, group.by = "celltype", data = df_mrg_nuc_filt, 
                            p.adjust.method = "fdr") %>% mutate(y_pos = max(df_mrg_nuc_filt$prop_cell_per_pt) + 0.05, p.adj = format.pval(p.adj, digits = 2))
  
  p0 = ggboxplot(df_mrg_nuc_filt, x = "ethnicity2", y = "prop_cell_per_pt",title = "prop_by_all_cell",
                 fill = "celltype", palette = col_pal,
                 add = "jitter", outlier.shape = NA, size = 0.3, alpha = 0.9, add.params = list(size  = 0.5, alpha = 0.7)) +  theme_sel +
    facet_wrap(~celltype, ncol = 14) + geom_signif(
      data=stat.test, 
      aes(xmin = group1, xmax = group2, y_position = y_pos, annotations = p.signif), # 
      manual= TRUE);p0
  
  png(paste0(save_path,"cell_proportions_ethnicity2_paper_patient.png"), width=8, height=5, units = "in", res = 300)
  print(p0)
  dev.off()
  
  pdf(paste0(save_path,"cell_proportions_ethnicity2_paper_paper_patient.pdf"), width=6, height=4)
  print(p0)
  dev.off()
}


library(cowplot)
theme_sel = theme(legend.position = "none",
                  strip.background = element_blank(),
                  panel.spacing.x = unit(0,"line"),
                  legend.background = element_blank(),
                  legend.title = element_text(size=8),
                  legend.text = element_text(size=8),
                  legend.key = element_blank(),
                  panel.grid = element_blank(),
                  #panel.border = element_blank(),
                  panel.border = element_rect(colour = "grey", fill=NA, size=0.3),
                  axis.line = element_blank(),
                  panel.background = element_blank(),
                  axis.text = element_text(size = 10),
                  axis.text.x = element_text(angle = 90),
                  axis.ticks = element_blank(),
                  axis.title=element_text(size=10),
                  title = element_text(size=8))

func_metadata_new = function(df_input = df_mrg_nuc3,
                         save_path = save_path,
                         pal_new = c("orange3", "hotpink3"),
                         celltype = "epithelial",
                         var = "menopause",
                         var_lvls = c("pre", "post"),
                         my_comparisons = list(c("pre", "post"))){
  
  # new test -----
  df_mrg_nuc_filt = df_input %>% 
    dplyr::select(patient_id, !!var, exp_proc_updated, surgery,  
                  celltype,  
                  # pt level summary
                  prop_celltype, n_celltype, n_tot)
  df_mrg_nuc_filt$var = df_mrg_nuc_filt[,var] %>% pull()
  
  # rm unknown menopuase
  df_mrg_nuc_filt = df_mrg_nuc_filt %>% dplyr::filter(!var %in% c("unknown"))
  #df_mrg_nuc_filt$var
  # conv. menopause to 2 levels
  df_mrg_nuc_filt$var2 = factor(df_mrg_nuc_filt$var, levels = var_lvls, ordered = T) 
  
  # wilcox test
  stat.test <- df_mrg_nuc_filt %>%
    dplyr::group_by(celltype) %>%
    wilcox_test(prop_celltype ~ var2) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj")
  #stat.test
  
  stat.test <- stat.test %>% add_xy_position(x = "celltype") %>% 
    mutate(p.adj = format.pval(p.adj, digits = 2))
  stat.test$custom.label <- ifelse(stat.test$p.adj <= 0.05, stat.test$p.adj, "ns")
  df_stat.test = data.frame(stat.test)
  df_stat.test <- apply(df_stat.test,2,as.character)
  
  write_rds(stat.test, paste0(save_path, "stattest_", celltype, "_", var, ".rds"))
  write.csv(df_stat.test, paste0(save_path, "stattest_", celltype, "_", var, ".csv"))
  
  # side by side boxplot
  g1 = ggboxplot(df_mrg_nuc_filt, x = "celltype", y = "prop_celltype", color = "var2", add = "jitter",  ylab = paste0("Prop of ", celltype," states"),
                 palette = pal_new)+
    stat_pvalue_manual(stat.test, label = "custom.label") + 
    rotate_x_text(angle = 45, hjust = 1, vjust = 1)
  ggsave(plot = g1,filename =  paste0(save_path, "boxplot_", celltype, "_", var, ".pdf"), width = 8.5, height = 5)
  
  
  # do_clinical_comp = function(df_mrg_nuc2 = df_mrg_nuc2){
  # stat.test0 = compare_means(formula = prop_celltype ~ var2, group.by = "celltype", 
  #                           data = df_mrg_nuc_filt, #%>% dplyr::filter(!var2 == 0), 
  #                           p.adjust.method = "fdr") %>% 
  #   mutate(y_pos = max(df_mrg_nuc_filt$prop_celltype) + 0.05, 
  #          p.adj = format.pval(p.adj, digits = 2))
  
  p0 = ggboxplot(df_mrg_nuc_filt, x = "var2", y = "prop_celltype",
                 title = "prop_by_all_cell",
                 #add = "jitter", 
                 fill = "celltype", palette = col_pal,
                 repel = T, outlier.shape = NA, size = 0.3, alpha = 0.9, add.params = list(size  = 0.5, alpha = 0.7)) +  
    theme_sel +
    facet_wrap(~celltype, ncol = 9) + 
    geom_signif(
      data=stat.test, 
      aes(xmin = group1, xmax = group2, y_position = y.position, annotations = custom.label), # 
      manual= TRUE);p0#+
  p1 = p0 + geom_jitter(color = "black", size = 0.3, alpha = 0.9, width = 0.1)
  save_plot(paste0(save_path,"p_boxplot1_", var, "-", celltype, ".pdf"), p0)
  save_plot(paste0(save_path,"p_boxplot2_", var, "-", celltype, ".pdf"), p1)
  
  ## check distributions of proportions
  p_load(ggridges, cowplot)
  
  # p1 = df_mrg_nuc_filt %>% 
  #   ggplot(aes(prop_celltype, cellgroup)) +
  #   geom_density_ridges() +
  #   theme_cowplot()
  
  p2 = df_mrg_nuc_filt %>% 
    ggplot(aes(prop_celltype, celltype)) +
    geom_density_ridges() +
    theme_cowplot()
  
  p = plot_grid(p2)
  save_plot(paste0(save_path,"p_ridge_", var, "-",celltype, ".pdf"), p)
  p_load(broom)
  
  # Y is binomial
  # lst_fit = group_by(df_mrg_nuc_filt, cellgroup) %>% 
  #   do(tidy(glm(bmi3 ~ prop_cell_per_pt + exp_proc_updated, family=binomial, data=.)))
  
  lst_fit = group_by(df_mrg_nuc_filt, celltype) %>% 
    do(generics::tidy(glm(var2 ~ prop_celltype + exp_proc_updated, family=binomial, data=.)))
  
  ### add CI
  
  lst_fit %<>% 
    mutate(ci_low = estimate - 1.96*std.error,
           ci_high = estimate + 1.96*std.error)
  
  write_tsv(lst_fit, glue({save_path},"lst_fit_", {var}, "_", {celltype}, ".csv"))
  
  lst_fit1 = lst_fit %>% dplyr::filter(term == "prop_celltype") %>% mutate(clinical_var = var)
  write_tsv(lst_fit1, glue({save_path},"lst_fit_", {var}, "_", {celltype}, "_supptable.csv"))
  
  ## create forest plot
  # devtools::install_github("rdboyes/forester")
  # p_load(forester)
  p_forest = lst_fit %>% 
    filter(term == "prop_celltype") %>% 
    ggplot(aes(estimate, celltype)) +
    geom_segment(aes(x = ci_low, xend = ci_high, y = celltype, yend = celltype)) +
    geom_vline(xintercept = 0, linetype = 2, alpha = 0.3) +
    geom_point(color = "red") +
    theme_cowplot() +
    theme_sel +
    ggtitle("Point estimate of cell type proportions on menopause", "Adjusted for protocol")
  p_forest
  
  save_plot(paste0(save_path, "p_forest_", var, "-", celltype, ".pdf"), p_forest)
}

func_metadata = function(df_input = df_mrg_nuc3,
                         save_path = save_path,
                         celltype = "epithelial",
                         var = "menopause",
                         var_lvls = c("pre", "post"),
                         my_comparisons = list(c("pre", "post"))){
  
  # wilcox test
  stat.test <- df_input %>%
    dplyr::group_by(celltype) %>%
    wilcox_test(prop_celltype ~ !!var) %>%
    adjust_pvalue() %>%
    add_significance("p.adj")
  #stat.test
  
  stat.test <- stat.test %>% add_xy_position(x = "celltype")
  stat.test$custom.label <- ifelse(stat.test$p.adj <= 0.05, stat.test$p.adj, "ns")
  
  # side by side boxplot
  g1 = ggboxplot(df_meta1, x = "celltype", y = "prop_celltype", color = var, add = "jitter",  ylab = paste0("Prop of ", celltype," states"),
                 palette = c("orange3", "hotpink3"))+
    stat_pvalue_manual(stat.test, label = "p.adj") + 
    rotate_x_text(angle = 45, hjust = 1, vjust = 1)
  ggsave(plot = g1,filename =  paste0("boxplot_", celltype, "_", var, ".pdf"), width = 8.5, height = 5)
  
  # new test -----
  df_mrg_nuc_filt = df_input %>% 
    dplyr::select(patient_id, !!var, exp_proc_updated,surgery,  
                  celltype, cellgroup, prop_cell_per_pt, 
                  surgery,
                  # pt level summary
                  prop_cell_per_grp_pt, cell_per_grp_pt, cell_per_pt, total_cells_per_pt) %>% 
    distinct()
  df_mrg_nuc_filt$var = df_mrg_nuc_filt[,var] %>% pull()
  
  # rm unknown menopuase
  df_mrg_nuc_filt = df_mrg_nuc_filt %>% dplyr::filter(!var %in% c("unknown"))
  df_mrg_nuc_filt$var
  # conv. menopause to 2 levels
  df_mrg_nuc_filt$var2 = factor(df_mrg_nuc_filt$var, levels = var_lvls, ordered = T) 
  
  # do_clinical_comp = function(df_mrg_nuc2 = df_mrg_nuc2){
  stat.test = compare_means(formula = prop_cell_per_pt ~ var2, group.by = "celltype", 
                            data = df_mrg_nuc_filt, #%>% dplyr::filter(!var2 == 0), 
                            p.adjust.method = "fdr") %>% 
    mutate(y_pos = max(df_mrg_nuc_filt$prop_cell_per_pt) + 0.05, 
           p.adj = format.pval(p.adj, digits = 2))
  
  p0 = ggboxplot(df_mrg_nuc_filt, x = "var2", y = "prop_cell_per_pt",
                 title = "prop_by_all_cell",
                 #add = "jitter", 
                 fill = "celltype", palette = col_pal,
                 repel = T, outlier.shape = NA, size = 0.3, alpha = 0.9, add.params = list(size  = 0.5, alpha = 0.7)) +  
    theme_sel +
    facet_wrap(~celltype, ncol = 9) + 
    geom_signif(
      data=stat.test, 
      aes(xmin = group1, xmax = group2, y_position = y_pos, annotations = p.signif), # 
      manual= TRUE);p0#+
  p1 = p0 + geom_jitter(color = "black", size = 0.3, alpha = 0.9, width = 0.1)
  save_plot(paste0(save_path,"p_boxplot_", var, "-", celltype, ".pdf"), p0)
  save_plot(paste0(save_path,"p_boxplot_", var, "-", celltype, ".pdf"), p1)
  
  ## check distributions of proportions
  p_load(ggridges, cowplot)
  
  p1 = df_mrg_nuc_filt %>% 
    ggplot(aes(prop_cell_per_pt, cellgroup)) +
    geom_density_ridges() +
    theme_cowplot()
  
  p2 = df_mrg_nuc_filt %>% 
    ggplot(aes(prop_cell_per_pt, celltype)) +
    geom_density_ridges() +
    theme_cowplot()
  
  p = plot_grid(p1, p2)
  save_plot(paste0(save_path,"p_ridge_", var, "-",celltype, ".pdf"), p)
  p_load(broom)
  
  # Y is binomial
  # lst_fit = group_by(df_mrg_nuc_filt, cellgroup) %>% 
  #   do(tidy(glm(bmi3 ~ prop_cell_per_pt + exp_proc_updated, family=binomial, data=.)))
  
  lst_fit = group_by(df_mrg_nuc_filt, celltype) %>% 
    do(generics::tidy(glm(var2 ~ prop_cell_per_pt + exp_proc_updated, family=binomial, data=.)))
  
  ### add CI
  
  lst_fit %<>% 
    mutate(ci_low = estimate - 1.96*std.error,
           ci_high = estimate + 1.96*std.error)
  
  write_tsv(lst_fit, glue({save_path},"lst_fit_", {var}, "_", {celltype}, ".csv"))
  
  lst_fit1 = lst_fit %>% dplyr::filter(term == "prop_cell_per_pt") %>% mutate(clinical_var = var)
  write_tsv(lst_fit1, glue({save_path},"lst_fit_", {var}, "_", {celltype}, "_supptable.csv"))
  
  ## create forest plot
  # devtools::install_github("rdboyes/forester")
  # p_load(forester)
  p_forest = lst_fit %>% 
    filter(term == "prop_cell_per_pt") %>% 
    ggplot(aes(estimate, celltype)) +
    geom_segment(aes(x = ci_low, xend = ci_high, y = celltype, yend = celltype)) +
    geom_vline(xintercept = 0, linetype = 2, alpha = 0.3) +
    geom_point(color = "red") +
    theme_cowplot() +
    theme_sel +
    ggtitle("Point estimate of cell type proportions on menopause", "Adjusted for protocol, surgery")
  p_forest
  
  save_plot(paste0(save_path, "p_forest_", var, "-", celltype, ".pdf"), p_forest)
}







