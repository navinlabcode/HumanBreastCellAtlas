## the scripts is for 10X visium ST analysis

HBCA40_ST <- Load10X_Spatial('./dir/outs/')
HBCA12_ST <- Load10X_Spatial('./dir/outs/')
BCMHBCA03L_ST <- Load10X_Spatial('./dir/outs/')
BCMHBCA04L_ST <- Load10X_Spatial('./dir/outs/')

HBCA_ST_int_list <- list(HBCA40_ST, HBCA12_ST, BCMHBCA03L_ST, BCMHBCA04L_ST)

for (i in 1:length(HBCA_ST_int_list)) {
  HBCA_ST_int_list[[i]] <- SCTransform(HBCA_ST_int_list[[i]], assay = 'Spatial')
}

HBCA_ST_int_features_all <- Reduce(intersect, lapply(HBCA_ST_int_list, function(x) rownames(x@assays$SCT@data)))
HBCA_ST_int_features <- SelectIntegrationFeatures(object.list = HBCA_ST_int_list, nfeatures = 3000)
HBCA_ST_int_list <- PrepSCTIntegration(object.list = HBCA_ST_int_list, anchor.features = HBCA_ST_int_features)

HBCA_ST_int_anchors <- FindIntegrationAnchors(object.list = HBCA_ST_int_list, normalization.method = "SCT", anchor.features = HBCA_ST_int_features, verbose=T)
HBCA_ST_int <- IntegrateData(anchorset = HBCA_ST_int_anchors, normalization.method = "SCT", features.to.integrate = HBCA_ST_int_features_all)

# after integration

HBCA_ST_int <- RunPCA(HBCA_ST_int)
HBCA_ST_int <- RunUMAP(HBCA_ST_int, dims = 1:20)

# Filter-1 ----------------------------------------------------------------
HBCA_ST_int_sub <- subset(HBCA_ST_int, nFeature_Spatial>180)
HBCA_ST_int_sub <- ScaleData(HBCA_ST_int_sub)
HBCA_ST_int_sub <- RunPCA(HBCA_ST_int_sub)
HBCA_ST_int_sub <- RunUMAP(HBCA_ST_int_sub, dims = 1:20)
HBCA_ST_int_sub <- FindNeighbors(HBCA_ST_int_sub, reduction='pca', dims=1:20)
HBCA_ST_int_sub <- FindClusters(HBCA_ST_int_sub, resolution = 2)

HBCA_ST_sub_markers_0 <- FindAllMarkers(HBCA_ST_int_sub, assay = 'Spatial')
HBCA_ST_sub_top_0 <- HBCA_ST_sub_markers_0 %>% group_by(cluster) %>%
  dplyr::filter(avg_logFC>0 & p_val_adj<0.1 & gene %!in% grep('^RPL|^RPS|^MT-', HBCA_ST_markers$gene, value=T)) %>%
  top_n(n=7, wt=avg_logFC)

Idents(HBCA_ST_int_sub) <- HBCA_ST_int_sub$seurat_clusters

# Filter-2 ----------------------------------------------------------------
idx_rm <- c('1', '2', '11', '13', '20')

####
HBCA_ST_int_sub <- subset(HBCA_ST_int_sub, seurat_clusters %!in% idx_rm)
HBCA_ST_int_sub$merge_clusters <- as.character(HBCA_ST_int_sub$seurat_clusters)
HBCA_ST_int_sub$merge_clusters[HBCA_ST_int_sub$merge_clusters %in% c('0', '5', '9')] <- 'ST05_fib'
HBCA_ST_int_sub$merge_clusters[HBCA_ST_int_sub$merge_clusters %in% c('3', '6', '18', '23', '24')] <- 'ST10_b'
HBCA_ST_int_sub$merge_clusters[HBCA_ST_int_sub$merge_clusters %in% c('4', '19')] <- 'ST03_lob'
HBCA_ST_int_sub$merge_clusters[HBCA_ST_int_sub$merge_clusters %in% c('7')] <- 'ST04_duct'
HBCA_ST_int_sub$merge_clusters[HBCA_ST_int_sub$merge_clusters %in% c('8', '17', '21')] <- 'ST02_lob'
HBCA_ST_int_sub$merge_clusters[HBCA_ST_int_sub$merge_clusters %in% c('10', '25')] <- 'ST08_vas'
HBCA_ST_int_sub$merge_clusters[HBCA_ST_int_sub$merge_clusters %in% c('12')] <- 'ST07_lym'
HBCA_ST_int_sub$merge_clusters[HBCA_ST_int_sub$merge_clusters %in% c('16')] <- 'ST09_peri'
HBCA_ST_int_sub$merge_clusters[HBCA_ST_int_sub$merge_clusters %in% c('14', '22')] <- 'ST06_adp'
HBCA_ST_int_sub$merge_clusters[HBCA_ST_int_sub$merge_clusters %in% c('15')] <- 'ST01_myo'
# HBCA_ST_int_sub$merge_clusters %<>% 
#   factor(., levels=c('ST10_b', 'ST09_peri', 'ST08_vas', 'ST07_lym', 'ST06_adp', 'ST05_fib', 'ST04_duct', 'ST03_lob', 'ST02_lob', 'ST01_myo'))
####

DimPlot(HBCA_ST_int_sub, group.by='merge_clusters', label=T) + NoLegend()
DimPlot(HBCA_ST_int_sub, group.by='sample') + scale_color_brewer(palette = 'Set1') + theme_void()
SpatialDimPlot(HBCA_ST_int_sub, group.by='merge_clusters')

Idents(HBCA_ST_int_sub) <- factor(HBCA_ST_int_sub$merge_clusters, sort(unique(HBCA_ST_int_sub$merge_clusters)))
HBCA_ST_sub_markers_1 <- FindAllMarkers(HBCA_ST_int_sub, assay = 'Spatial', logfc.threshold=0, return.thresh=0.05)
HBCA_ST_top_7 <- HBCA_ST_sub_markers_1 %>% group_by(cluster) %>%
  dplyr::filter(avg_logFC>0 & p_val_adj<0.1 & gene %!in% grep('^RPL|^RPS|^MT-', HBCA_ST_sub_markers_1$gene, value=T)) %>%
  top_n(n=7, wt=avg_logFC)
HBCA_ST_top_3 <- HBCA_ST_sub_markers_1 %>% group_by(cluster) %>%
  dplyr::filter(avg_logFC>0 & p_val_adj<0.1 & gene %!in% grep('^RPL|^RPS|^MT-', HBCA_ST_sub_markers_1$gene, value=T)) %>%
  top_n(n=3, wt=avg_logFC)

get_palette('Paired', 15)
plot(c(1:15), col=get_palette('Paired', 15), pch=16, cex=10)
text(x=c(1:15), y=c(1:15), labels = c(1:15))
## Epi_myo:15; Lob_lumhr:10; Lob_lummix:7; Duct_lumsec:12; ConTis_fib:1; Adp_adp:14; Lym_endlym:3; Vas_endvas:6; Vas_peri:11; ConTis_b ##
col_order <- c(15, 10, 7, 12, 1, 14, 3, 6, 11, 13)

DimPlot(HBCA_ST_int_sub, cols=get_palette('Paired', 15)[col_order], label=T) + theme_void() + NoLegend()
# SpatialDimPlot(HBCA_ST_int_sub, images = 'slice1', pt.size.factor = 0) + NoLegend()
SpatialDimPlot(HBCA_ST_int_sub, images = 'slice1', pt.size.factor=1.2, stroke=.1) + 
  scale_fill_manual(values = get_palette('Paired', 15)[col_order]) + rremove('legend.title') + 
  theme(legend.key.size = unit(1, 'mm'), legend.position = c(.2, .1)) + 
  font('legend.text', size=8) + guides(fill=guide_legend(ncol=2))

HBCA_ST_int_avg <- AverageExpression(HBCA_ST_int_sub, assays='Spatial', return.seurat=T)
DoHeatmap(HBCA_ST_int_avg, features=HBCA_ST_top_7$gene, size=3, group.bar.height=.01, 
          group.colors=get_palette('Paired', 15)[col_order], disp.min = -2, disp.max = 2, draw.lines=F) + 
  scale_fill_gradientn(colours = rev(brewer.pal(10, 'RdBu'))) + NoLegend()

DoHeatmap(HBCA_ST_int_avg, features=HBCA_ST_top_3$gene, size=3, group.bar.height=.01, 
          group.colors=get_palette('Paired', 15)[col_order], disp.min = -2, disp.max = 2, draw.lines=F) + 
  scale_fill_gradientn(colours = rev(brewer.pal(10, 'RdBu'))) + NoLegend()

# SaveData ----------------------------------------------------------------
saveRDS(HBCA_ST_int_sub, 'HBCA_ST_int_sub.rds')