# Loading packages ---------------------------------------------------------
setwd("your/working/director")
library(Seurat)
library(SeuratWrappers)
require(Hmisc)
require(zoo)
require(RColorBrewer)
require(plotly)
library(future)
library(GeneOverlap)

# Preparation -------------------------------------------------------------
color_temp <- c('#B15928', '#FB820F', '#E93E3F', '#8C66AF', '#A6CEE3', '#EEDB80', '#72B29C', '#EC9A91', '#D1AAB7', '#A99099')

'%!in%' <- function(x,y) {!('%in%'(x,y))}

clust_gene_ovlp <- function(top_de, marker_list, bk_gene) {
  top_de$cluster %<>% as.character()
  clust_query <- unique(top_de$cluster)
  clust_pval_df <- matrix(NA, length(marker_list), length(clust_query)) %>%
    set_rownames(names(marker_list)) %>% set_colnames(clust_query) %>% data.frame
  clust_or_df <- matrix(NA, length(marker_list), length(clust_query)) %>%
    set_rownames(names(marker_list)) %>% set_colnames(clust_query) %>% data.frame
  for (i in 1:length(clust_query)) {
    gene_i <- top_de$gene[top_de$cluster==clust_query[i]]
    ovlp_res <- sapply(marker_list, function(marker_x) {
      ovlp_test <- GeneOverlap::newGeneOverlap(listA = gene_i, listB = marker_x, genome.size = length(unique(bk_gene))) %>% 
        GeneOverlap::testGeneOverlap()
      # if (is.infinite(fshexact_test$estimate)) fshexact_test <- fisher.test(matrix(c(x11, x21, x12, x22)+1, nrow = 2, ncol = 2))
      c(pval=ovlp_test@pval, or=ovlp_test@odds.ratio)
    })
    clust_pval_df[, i] <- ovlp_res[1, rownames(clust_pval_df)]
    clust_or_df[, i] <- ovlp_res[2, rownames(clust_pval_df)]
  }
  return(list(pval=clust_pval_df, or=clust_or_df))
}

# 1. ST Analysis -------------------------------------------------------------
# 1.1 Load data ---------------------------------------------------------------
HBCA_ST_int_sub <- readRDS('HBCA_ST_int_sub_06212022.rds')
HBCA_ST_markers <- readRDS('HBCA_ST_markers_06212022.rds')
HBCA_ST_sub_markers_1 <- readRDS('HBCA_ST_sub_markers_1_06212022.rds')
HBCA_marker_cmbn_lis <- readRDS('HBCA_marker_cmbn_lis.rds')
HBCA_marker_df_filter <- read.csv('../HBCA_cells_filtered_markers_tps_022321.csv', sep=',', header=T, row.names = 1)
HBCA_sn_marker_df_adpmast <- readRDS('HBCA_sn_marker_df_adpmast_06220222.rds')
HBCA_SC_int_sub <- readRDS('~/Analysis/HBCA/HBCA_new/HBCA_SC_int_sub_50k_013121.rds')
Idents(HBCA_SC_int_sub) <- HBCA_SC_int_sub$final_group1
HBCA_SN_int_sub <- readRDS('~/Analysis/HBCA/HBCA_new/HBCA_SN_int_sub_5k_022321.rds')
HBCA_ST_int_epi_sub <- readRDS('HBCA_ST_int_epi_sub_06222022.rds')
HBCA_sn_marker_df_all <- read.csv('../HBCA_allnuc_markers_022621.csv', sep=',', header=T, row.names = 1)

# 1.2 Figure2 ---------------------------------------------------------------
Figure2_a <- SpatialDimPlot(HBCA_ST_int_sub, group.by='merge_clusters', images='slice1') + scale_fill_manual(values = color_temp)

## DE Analysis ##
HBCA_ST_sub_markers_1 <- FindAllMarkers(HBCA_ST_int_sub, assay = 'Spatial', logfc.threshold=0, return.thresh=0.05)

HBCA_ST_top_7 <- HBCA_ST_sub_markers_1 %>% group_by(cluster) %>%
  dplyr::filter(avg_logFC>0 & p_val_adj<0.1 & gene %!in% grep('^RPL|^RPS|^MT-', HBCA_ST_sub_markers_1$gene, value=T)) %>%
  top_n(n=7, wt=avg_logFC)
HBCA_ST_top_3 <- HBCA_ST_sub_markers_1 %>% group_by(cluster) %>%
  dplyr::filter(avg_logFC>0 & p_val_adj<0.1 & gene %!in% grep('^RPL|^RPS|^MT-', HBCA_ST_sub_markers_1$gene, value=T)) %>%
  top_n(n=3, wt=avg_logFC)

HBCA_ST_int_avg <- AverageExpression(HBCA_ST_int_sub, assays='Spatial', return.seurat=T)

DoHeatmap(HBCA_ST_int_avg, features=HBCA_ST_top_7$gene, size=3, group.bar.height=.01, 
          group.colors=color_temp, disp.min = -2, disp.max = 2, draw.lines=F) + 
  scale_fill_gradientn(colours = rev(brewer.pal(10, 'RdBu'))) + NoLegend()

Figure2_b <- DoHeatmap(HBCA_ST_int_avg, features=HBCA_ST_top_3$gene, size=3, group.bar.height=.01, group.colors=color_temp, disp.min = -2, disp.max = 2, draw.lines=F) + 
  scale_fill_gradientn(colours = rev(brewer.pal(10, 'RdBu'))) + NoLegend()


# 1.3 SI Figure 2 -------------------------------------------------------------
SF2_a <- DimPlot(HBCA_ST_int_sub, cols=color_temp, label=T) + theme_void() + NoLegend()


## Gene Overlap Test ##
HBCA_ST_top20 <- HBCA_ST_sub_markers_1 %>% group_by(cluster) %>%
  dplyr::filter(avg_logFC>0, p_val_adj<0.1, 
                # pct.2<.5, pct.1>0.3, pct.1>pct.2,
                gene %!in% grep('^RPL|^RPS|^MT-', HBCA_ST_sub_markers_1$gene, value=T)) %>%
  top_n(n=20, wt=avg_logFC)

HBCA_ST_marker_ovlp_res <- clust_gene_ovlp(data.frame(HBCA_ST_top20), HBCA_marker_cmbn_lis, bk_gene=unique(HBCA_ST_markers$gene))

HBCA_ST_marker_ovlp_pval <- -log(HBCA_ST_marker_ovlp_res$pval)
HBCA_ST_marker_ovlp_pval <- HBCA_ST_marker_ovlp_pval[apply(HBCA_ST_marker_ovlp_pval, 1, sum)!=0, ]
HBCA_ST_marker_ovlp_pval$Cell <- rownames(HBCA_ST_marker_ovlp_pval)
HBCA_ST_marker_ovlp_pval_melt <- tidyr::pivot_longer(HBCA_ST_marker_ovlp_pval, cols=-Cell, names_to='Group', values_to='neglogP')
HBCA_ST_marker_ovlp_pval_melt$Group %<>% factor(., levels=c('ST10_b', 'ST09_peri', 'ST08_vas', 'ST07_lym', 'ST06_adp', 'ST05_fib', 'ST04_duct', 'ST03_lob', 'ST02_lob', 'ST01_myo'))
HBCA_ST_marker_ovlp_pval_melt$Cell %<>% factor(., levels=c('Basal', 'LumHR', 'LumSec', 'Fibroblasts', 'Adipocytes', 'Lymphatic','Vascular', 'Pericytes',
                                                           'Myeloid', 'B.cells', 'Mast'))

HBCA_ST_marker_ovlp_or <- log1p(HBCA_ST_marker_ovlp_res$or)
HBCA_ST_marker_ovlp_or <- HBCA_ST_marker_ovlp_or[apply(HBCA_ST_marker_ovlp_or, 1, sum)!=0, ]
HBCA_ST_marker_ovlp_or$Cell <- rownames(HBCA_ST_marker_ovlp_or)
HBCA_ST_marker_ovlp_or_melt <- tidyr::pivot_longer(HBCA_ST_marker_ovlp_or, cols=-Cell, names_to='Group', values_to='logOR')
HBCA_ST_marker_ovlp_or_melt$Group %<>% factor(., levels=c('ST10_b', 'ST09_peri', 'ST08_vas', 'ST07_lym', 'ST06_adp', 'ST05_fib', 'ST04_duct', 'ST03_lob', 'ST02_lob', 'ST01_myo'))
HBCA_ST_marker_ovlp_or_melt$Cell %<>% factor(., levels=c('Basal', 'LumHR', 'LumSec', 'Fibroblasts', 'Adipocytes', 'Lymphatic','Vascular', 'Pericytes',
                                                         'Myeloid', 'B.cells', 'Mast'))
HBCA_ST_marker_ovlp <- dplyr::full_join(HBCA_ST_marker_ovlp_or_melt, HBCA_ST_marker_ovlp_pval_melt)

SF2_b <- ggplot(HBCA_ST_marker_ovlp, aes(Cell, Group)) + geom_point(aes(color=logOR, size=neglogP)) + 
  theme_bw() + theme(text = element_text(size=10)) + 
  scale_color_gradient2(low = 'white',  high = 'red') + 
  labs(color="logOR", size="-logP", y='ST Group') + rotate_x_text(angle = 45)

SF2_c <- SpatialDimPlot(HBCA_ST_int_sub, group.by='merge_clusters', images='slice1.1') + scale_fill_manual(values = color_temp)
SF2_d <- SpatialDimPlot(HBCA_ST_int_sub, group.by='merge_clusters', images='slice1.3') + scale_fill_manual(values = color_temp)
SF2_e <- SpatialDimPlot(HBCA_ST_int_sub, group.by='merge_clusters', images='slice1.2') + scale_fill_manual(values = color_temp)

## ST-SC marker ##
HBCA_SCST_marker_top7 <- HBCA_marker_df_filter %>% group_by(cluster) %>% top_n(7, avg_logFC) %>% bind_rows(., HBCA_ST_top_7) %>% bind_rows(., HBCA_sn_marker_df_adpmast)
HBCA_SCST_marker_top7 <- HBCA_SCST_marker_top7[HBCA_SCST_marker_top7$gene %in% intersect(rownames(HBCA_ST_int_sub@assays$Spatial@data), rownames(HBCA_SC_int_sub@assays$RNA@data)), ]

STSC_group <- list(c('Basal', 'ST01_myo'), c('LumHR', 'ST02_lob'), c('LumSec', 'ST03_lob'), c('LumSec', 'ST04_duct'), 
                   c('Fibroblasts', 'ST05_fib'), c('Adipocytes', 'ST06_adp'), 
                   c('Lymphatic', 'ST07_lym'), c('Vascular', 'ST08_vas'), c('Pericytes', 'ST09_peri'), c('B-cells', 'ST10_b'))

ST_SC_data_temp <- bind_rows(
  data.frame(HBCA_ST_int_sub@assays$Spatial@data[unique(HBCA_SCST_marker_top7$gene), ] %>% t, type='ST', group=HBCA_ST_int_sub$merge_clusters), 
  data.frame(HBCA_SC_int_sub@assays$RNA@data[unique(HBCA_SCST_marker_top7$gene), ] %>% t, type='SC', group=HBCA_SC_int_sub$final_group1),
  data.frame(HBCA_SN_int_sub@assays$RNA@data[unique(HBCA_SCST_marker_top7$gene), ] %>% t, type='SN', group=HBCA_SN_int_sub$final_group1)
)

## DotPlot ##
for (i in 1:length(STSC_group)) {
  print(i)
  STSC_i <- STSC_group[[i]]
  feat_temp <- HBCA_SCST_marker_top7$gene[HBCA_SCST_marker_top7$cluster %in% c(STSC_i)]
  data_temp <- ST_SC_data_temp %>% dplyr::select(contains(feat_temp), type, group) %>% 
    dplyr::filter(group %in% STSC_i) %>% 
    tidyr::pivot_longer(cols=contains(feat_temp), names_to='Gene', values_to='Expression')
  data_temp$group %<>% factor(., levels = STSC_i)
  feat_temp_sub <- data_temp %>% group_by(type, Gene) %>% summarise(AvgExp = mean(Expression)) %>% ungroup() %>% 
    group_by(Gene) %>% summarise(AvgExp = mean(AvgExp)) %>% ungroup %>% top_n(7, AvgExp) %>% extract2('Gene')
  data_temp <- data_temp %>% dplyr::filter(Gene %in% feat_temp_sub)
  data_temp <- data_temp %>% group_by(group, Gene) %>% summarise(mean_exp=mean(Expression), prop_exp=mean(Expression>0))
  data_temp$mean_exp <- oob_squish(data_temp$mean_exp, range=c(0, 3))
  data_temp$prop_exp <- oob_squish(data_temp$prop_exp, range=c(.1, .8))
  p_i <- ggplot(data_temp, aes(group, Gene)) + 
    geom_point(aes(color=mean_exp, size=prop_exp)) + 
    scale_color_viridis(limits = c(0, 3)) + 
    scale_size_continuous(limits = c(.1, .8)) + 
    theme_bw() + NoLegend() + 
    rremove('grid') + rremove('xlab') + rremove('ylab') + rotate_x_text(45) + 
    theme(strip.text=element_text(size=6), 
          panel.border = element_rect(fill=NA, size=0.5, linetype="solid"),
          text=element_text(size=6), axis.text = element_text(size = 6),
          plot.margin = unit(c(0,0,0,0), "mm"))
  plot_str <- paste0('p_', i, ' <- p_i')
  eval(parse(text=plot_str))
}

SF2_e <- p_1 + p_2 + p_3 + p_4 + p_5 + p_6 + p_7 + p_8 + p_9 + p_10 + plot_layout(ncol=5)


# 1.4 Figure 3 ---------------------------------------------------------
SpatialDimPlot(HBCA_ST_int_epi_sub, images = 'slice1', pt.size.factor=1.5) + 
  scale_fill_brewer(palette = 'Set1') + theme(legend.position = c(.1, .1))
SpatialDimPlot(HBCA_ST_int_epi_sub, images = 'slice1.1', pt.size.factor=1.5) + 
  scale_fill_brewer(palette = 'Set1') + theme(legend.position = c(.1, .1))
SpatialDimPlot(HBCA_ST_int_epi_sub, images = 'slice1.2', pt.size.factor=1.5) + 
  scale_fill_brewer(palette = 'Set1')
SpatialDimPlot(HBCA_ST_int_epi_sub, images = 'slice1.3', pt.size.factor=1) + 
  scale_fill_brewer(palette = 'Set1') + theme(legend.position = c(.1, .8))

## DE between Lobule vs Duct ##
HBCA_ST_int_lobduct_marker <- FindAllMarkers(HBCA_ST_int_epi_sub, assay='Spatial')
HBCA_ST_int_lobduct_marker_top <- HBCA_ST_int_lobduct_marker %>% 
  dplyr::filter(avg_logFC>0, pct.1>.3, p_val_adj<.1, !grepl('^RPL|^RPS|^MT-', gene)) %>% 
  group_by(cluster) %>% top_n(15, avg_logFC) %>% group_by(gene) %>% top_n(1, avg_logFC)

HBCA_ST_int_epi_sub_avg_exp <- AverageExpression(HBCA_ST_int_epi_sub, features = HBCA_ST_int_lobduct_marker_top$gene, return.seurat=TRUE, assays='Spatial')

## Recalculate Average Expression after Scale ##
HBCA_ST_int_epi_sub_avg_exp_duct <- apply(HBCA_ST_int_epi_sub@assays$Spatial@scale.data[, Cells(HBCA_ST_int_epi_sub)[HBCA_ST_int_epi_sub$LobDuct=='Duct']], 1, mean)
HBCA_ST_int_epi_sub_avg_exp_lob <- apply(HBCA_ST_int_epi_sub@assays$Spatial@scale.data[, Cells(HBCA_ST_int_epi_sub)[HBCA_ST_int_epi_sub$LobDuct=='Lobule']], 1, mean)
HBCA_ST_int_epi_sub_avg_exp@assays$Spatial@scale.data[, 1] <- HBCA_ST_int_epi_sub_avg_exp_duct[HBCA_ST_int_lobduct_marker_top$gene]
HBCA_ST_int_epi_sub_avg_exp@assays$Spatial@scale.data[, 2] <- HBCA_ST_int_epi_sub_avg_exp_lob[HBCA_ST_int_lobduct_marker_top$gene]

Figure3_q <- DoHeatmap(HBCA_ST_int_epi_sub_avg_exp, features = HBCA_ST_int_lobduct_marker_top$gene, assay="Spatial", slot="scale.data", 
          label = TRUE, size=4, draw.lines = F, group.colors=get_palette('Set1', 3)[c(1, 2)], raster = F)  +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato3", na.value = "00000000", midpoint = 0)

# 1.5 Figure 5 ------------------------------------------------------------
## Adipocytes ##
HBCA_ST_top_adp <- HBCA_ST_sub_markers_1 %>% dplyr::filter(cluster=='ST06_adp') %>% 
  dplyr::filter(avg_logFC>0, p_val_adj<0.1, pct.1>.25, pct.2<.4, gene %!in% grep('^RPL|^RPS|^MT-', HBCA_ST_sub_markers_1$gene, value=T)) %>% 
  top_n(n=15, wt=avg_logFC)

HBCA_ST_int_adp <- subset(HBCA_ST_int_sub, subset=merge_clusters=='ST06_adp')
DefaultAssay(HBCA_ST_int_adp) <- 'integrated'

SpatialDimPlot(HBCA_ST_int_adp, images = 'slice1.1', group.by='merge_clusters', pt.size.factor = 1) + 
  scale_fill_manual(values=get_palette('Set3', 5)[2]) + theme(legend.position = c(.2, .2))

SpatialDimPlot(HBCA_ST_int_adp, images = 'slice1.2', group.by='merge_clusters', pt.size.factor = 1) + 
  scale_fill_manual(values=get_palette('Set3', 5)[2]) + theme(legend.position = c(.2, .2))


## Adipocytes ##
HBCA_SN_adp_marker_top15 <- HBCA_sn_marker_df_all %>% dplyr::filter(cluster=='Adipocytes') %>% top_n(15, avg_logFC)
HBCA_ST_adp_marker_top15 <- HBCA_ST_sub_markers_1 %>% dplyr::filter(cluster=='ST06_adp') %>% top_n(15, avg_logFC)
HBCA_SNST_adp_marker_top15 <- unique(c(HBCA_SN_adp_marker_top15$gene, HBCA_ST_adp_marker_top15$gene))

adp_marker <- c('ADH1B', 'CD36', 'PLIN1', 'PLIN4', 
                'ADIPOQ', 'FABP4', 'LEP', 'LPL', 
                'GK', 'PDE4A', 'PRDM16', 'UCP1')

ST_SC_data_temp <- bind_rows(
  data.frame(HBCA_ST_int_sub@assays$Spatial@data[adp_marker, ] %>% t, type='ST', group=HBCA_ST_int_sub$merge_clusters, check.names = FALSE), 
  data.frame(HBCA_SN_int_sub@assays$RNA@data[adp_marker, ] %>% t, type='SN', group=HBCA_SN_int_sub$final_group1, check.names = FALSE)
) %>% dplyr::filter(group %in% c('Adipocytes', 'ST06_adp'))

data_temp <- ST_SC_data_temp %>% tidyr::pivot_longer(cols=-c(type, group), names_to='Gene', values_to='Expression')
data_temp$Gene %<>% factor(., levels=adp_marker)

Figure5_l <- ggplot(data_temp, aes(group, Expression)) + geom_violin(aes(fill=group)) + 
  facet_wrap(.~Gene, scales='free', ncol=3, dir='v') + 
  theme_bw() + 
  rremove('grid') + rremove('xlab') + rremove('x.text') + rremove('x.ticks') + NoLegend() + 
  scale_fill_brewer(palette = 'Set1') +
  theme(strip.text=element_text(size=6), 
        panel.border = element_rect(fill=NA, size=0.5, linetype="solid"),
        text=element_text(size=6), plot.margin = unit(c(0,0,0,0), "mm"), 
        legend.spacing = unit(.1,"cm"), legend.spacing.x = unit(.1,"cm"), legend.spacing.y = unit(.1,"cm"))

# 2. Resolve Analysis -------------------------------------------------------------
# 2.1 Load data -----------------------------------------------------------
resolve_hbca_A2_batch1_srt <- readRDS('resolve_hbca_A2_batch1_srt_06222022.rds')
resolve_hbca_B1_batch1_srt <- readRDS('resolve_hbca_B1_batch1_srt_06222022.rds')
resolve_hbca_C1_batch1_srt <- readRDS('resolve_hbca_C1_batch1_srt_06222022.rds')
hbca_resolve_A1_1_batch2 <- readRDS('hbca_resolve_A1_1_batch2_06222022.rds')
hbca_resolve_A2_1_batch2 <- readRDS('hbca_resolve_A2_1_batch2_06222022.rds')
hbca_resolve_B1_1_batch2 <- readRDS('hbca_resolve_B1_1_batch2_06222022.rds')
hbca_resolve_B2_1_batch2 <- readRDS('hbca_resolve_B2_1_batch2_06222022.rds')
hbca_resolve_D2_1_batch2 <- readRDS('hbca_resolve_D2_1_batch2_06222022.rds')
hbca_resolve_A1_batch3 <- readRDS('hbca_resolve_A1_batch3_06222022.rds')
hbca_resolve_B1_batch3 <- readRDS('hbca_resolve_B1_batch3_06222022.rds')
hbca_resolve_B2_batch3 <- readRDS('hbca_resolve_B2_batch3_06222022.rds')
hbca_resolve_D1_batch3 <- readRDS('hbca_resolve_D1_batch3_06222022.rds')
hbca_resolve_D2_batch3 <- readRDS('hbca_resolve_D2_batch3_06222022.rds')
resolve_hbca_merge <- readRDS('resolve_hbca_merge_06222022.rds')
resolve_hbca_merge_sp_df <- readRDS('resolve_hbca_merge_sp_df_06222022.rds')
resolve_hbca_merge_mannot_roi_sum <- readRDS('resolve_hbca_merge_mannot_roi_sum_06222022.rds')
resolve_hbca_srt_merge_sp_df <- readRDS('resolve_hbca_srt_merge_sp_df_06222022.rds')

DimPlot(subset(hbca_resolve_D2_batch3, subset=celltype!='lowconf'), group.by = 'celltype', reduction='SP') + 
  scale_color_manual(values = color_temp) + theme_bw() + rremove('xylab') + rremove('xy.text') + rremove('ticks')

# 2.2 Figure2 -----------------------------------------------------------------
Figure2_f <- DimPlot(subset(resolve_hbca_A2_batch1_srt, subset=celltype!='lowconf'), group.by = 'celltype', reduction='SP') + 
  scale_color_manual(values = color_temp) + theme_bw() + rremove('xylab') + rremove('xy.text') + rremove('ticks')

## DE Analysis ##
resolve_hbca_merge_deg <- FindAllMarkers(resolve_hbca_merge)
resolve_hbca_merge_avg <- AverageExpression(resolve_hbca_merge, return.seurat = TRUE)
resolve_hbca_merge_deg_top <- resolve_hbca_merge_deg %>% dplyr::filter(avg_log2FC>0, pct.2<.5) %>% group_by(cluster) %>% top_n(5, avg_log2FC)

Figure2_e <- DoHeatmap(resolve_hbca_merge_avg, rev(resolve_hbca_merge_deg_top$gene), disp.min=-2, disp.max = 2, draw.lines=F, group.colors=color_temp, angle=90) + 
  scale_fill_gradientn(colours=viridis(20))

## Spatial cell proximity ##
resolve_hbca_srt_merge_sp_df_dt <- CellTrek:::DT_boot_mst(resolve_hbca_srt_merge_sp_df[, c(3, 4)], 
                                                          coord_df=resolve_hbca_srt_merge_sp_df[, c(1, 2)], col_cell='cell_type', boot_n=20, dist_cutoff=1000)

resolve_hbca_srt_merge_sp_df_meta <- data.frame(id=rownames(resolve_hbca_srt_merge_sp_df_dt$mst_cons), 
                                                type=c('epi', 'strm', 'epi', 'epi', 'endo', 'imm', 'strm', 'imm', 'endo', 'imm'), 
                                                freq=c(table(resolve_hbca_merge$celltype)[rownames(resolve_hbca_srt_merge_sp_df_dt$mst_cons)]))

resolve_hbca_srt_merge_sp_df_dt_boot <- apply(resolve_hbca_srt_merge_sp_df_dt$boot_array, c(1, 2), mean)
diag(resolve_hbca_srt_merge_sp_df_dt_boot) <- NA
resolve_hbca_srt_merge_sp_df_dt_boot <- max(resolve_hbca_srt_merge_sp_df_dt_boot, na.rm = T) - resolve_hbca_srt_merge_sp_df_dt_boot
CellTrek::scoloc_vis(resolve_hbca_srt_merge_sp_df_dt_boot, meta_data = resolve_hbca_srt_merge_sp_df_meta)

## Area calculation ##
p0 <- ggplot(data=resolve_hbca_merge_mannot_roi_sum %>% na.omit %>% dplyr::filter(celltype!='lowconf' & roi_area_type!='Lobule_abnm'), aes(x=roi_type_abbr, y=cell_n, fill=celltype)) + 
  geom_bar(color='black', stat='identity', position='fill') + scale_fill_manual(values = color_temp) + 
  theme_bw(base_size = 12) + ylab('Proportion') + NoLegend() + coord_flip() + rremove('grid')
p1 <- ggplot(data=resolve_hbca_merge_mannot_roi_sum %>% na.omit %>% dplyr::filter(celltype!='lowconf' & roi_area_type!='Lobule_abnm'), aes(x=roi_type_abbr, y=cell_dens, fill=celltype)) + 
  geom_bar(color='black', stat='identity') + scale_fill_manual(values = color_temp) + 
  theme_bw() + ylab('Density') + NoLegend() + coord_flip() + rremove('grid')
Figure2_h <- plot_grid(p0, p1, nrow=2)

# 2.3 SI Figure 3 ---------------------------------------------------------
SF3_a <- ggplot(resolve_hbca_merge_sp_df %>% na.omit, aes(SP_1, SP_2)) + geom_point(aes(color=celltype), size=.01) + 
  facet_wrap(.~sample, scales = 'free', nrow=3) + scale_color_manual(values = color_temp) + 
  theme_bw() + rremove('xylab') + rremove('xy.text') + rremove('ticks') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(plot.title = element_text(hjust = 0.5)) + guides(colour = guide_legend(override.aes = list(size=2)))

SF3_b <- ggplot(resolve_hbca_merge_sp_df %>% na.omit, aes(SP_1, SP_2)) + geom_point(aes(color=roi_area_type), size=.01) + 
  facet_wrap(.~sample, scales = 'free', nrow=3) + scale_color_brewer(palette = 'Set1') + 
  theme_bw() + rremove('xylab') + rremove('xy.text') + rremove('ticks') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(plot.title = element_text(hjust = 0.5)) + guides(colour = guide_legend(override.aes = list(size=2)))

# 2.4 Figure3 -------------------------------------------------------------
## Lob vs Duct ##
resolve_hbca_merge_epi <- subset(resolve_hbca_merge, roi_area_type %in% c('Duct', 'Lobule'))
Idents(resolve_hbca_merge_epi) <- factor(resolve_hbca_merge_epi$roi_area_type, levels=c('Duct', 'Lobule'))
resolve_hbca_merge_area_deg <- FindAllMarkers(resolve_hbca_merge_epi)
resolve_hbca_merge_area_deg_top <- resolve_hbca_merge_area_deg %>% dplyr::filter(pct.1>.3) %>% group_by(cluster) %>% top_n(5, avg_log2FC)

p0 <- StackedVlnPlot(resolve_hbca_merge_epi, resolve_hbca_merge_area_deg_top$gene[1:5], pt.size = 0, 
                     plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"))
p1 <- StackedVlnPlot(resolve_hbca_merge_epi, resolve_hbca_merge_area_deg_top$gene[6:10], pt.size = 0, 
                     plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"))

Figure3_s <- p0 | p1


# 2.5 SI Figure5 ----------------------------------------------------------
## KI67 ##
resolve_hbca_merge_sp_df_sub <- resolve_hbca_merge_sp_df %>% na.omit %>% 
  dplyr::filter(celltype %in% c('basal', 'lumhr', 'lumsec') & sample %in% c('B1_batch1', 'A1_1_batch2', 'A2_1_batch2', 'A1_batch3'))

SF5_d <- ggplot(data=resolve_hbca_merge_sp_df_sub, aes(SP_1, SP_2)) + 
  geom_point(data=resolve_hbca_merge_sp_df_sub %>% dplyr::filter(MKI67_bin=='neg'), size=.1, color='grey') + 
  geom_point(data=resolve_hbca_merge_sp_df_sub %>% dplyr::filter(MKI67_bin=='pos'), size=.5, color='red') + 
  facet_wrap(.~sample, scales = 'free', nrow=3) + 
  theme_bw() + rremove('xylab') + rremove('xy.text') + rremove('ticks') + 
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(plot.title = element_text(hjust = 0.5))

# 2.6 Figure 4 ------------------------------------------------------------
## Immune ##
resolve_hbca_merge_sp_imm_df <- resolve_hbca_merge_sp_df[, 1:5]
resolve_hbca_merge_sp_imm_df$imm <- 'other'
resolve_hbca_merge_sp_imm_df$imm[resolve_hbca_merge_sp_imm_df$celltype %in% c('myeloid', 'tcell', 'bcell')] <- 'immune'
resolve_hbca_merge_imm_summary <- table(resolve_hbca_merge_sp_imm_df$imm) %>% data.frame(data='Resolve', .) %>% set_colnames(c('data', 'imm', 'cell_n'))
resolve_hbca_merge_imm_summary$imm %<>% factor(., levels=c('other', 'immune'))

ggplot(data=resolve_hbca_merge_imm_summary, aes(x=data, y=cell_n, fill=imm)) + 
  geom_bar(color='black', stat='identity', position='fill') + 
  theme_bw(base_size = 12) + ylab('Proportion') + rremove('grid')

# 2.7 SI Figure 7 ---------------------------------------------------------
p0 <- ggplot(data=resolve_hbca_merge_mannot_roi_sum %>% na.omit %>% dplyr::filter(celltype %in% c('myeloid', 'tcell', 'bcell') & roi_area_type!='Lobule_abnm'), 
             aes(x=roi_type_abbr, y=cell_n, fill=celltype)) + 
  geom_bar(color='black', stat='identity', position='fill') + scale_fill_manual(values = immune_color_temp) + 
  theme_bw(base_size = 12) + ylab('Proportion') + rremove('grid')
p1 <- ggplot(data=resolve_hbca_merge_mannot_roi_sum %>% na.omit %>% dplyr::filter(celltype %in% c('myeloid', 'tcell', 'bcell') & roi_area_type!='Lobule_abnm'), 
             aes(x=roi_type_abbr, y=cell_dens, fill=celltype)) + 
  geom_bar(color='black', stat='identity') + scale_fill_manual(values = immune_color_temp) + 
  theme_bw(base_size = 12) + ylab('Density') + rremove('grid')
plot_grid(p1, p0, ncol=2)

