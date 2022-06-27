
# setwd -------

setwd("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/")
odir = "01_read_files"; dir.create(odir)
#plan("multiprocess", workers = 100)
options(future.globals.maxSize = Inf)
setwd(paste0("/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/", odir))
wd = "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/"
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_load_packages.R")
source("/volumes/lab/users/tapsi/projects/bcap/.wd/20201212_final2/00_functions_final1.R")


# Read TRACKER ---------------------------------------------------------------
snames = c(paste0("hbca_c",seq(46,100)))
tracker = read_sheet(x = "/volumes/lab/users/tapsi/projects/bcap/hbca_sample_tracker.xlsx", sheet = 2) %>% as_tibble() %>% clean_names()
length(snames) #55


# Read 10X ---------------------------------------------------------------
# read 10x data
data_path = tracker %>% dplyr::filter(samp_id %in% snames) %>% pull(path_v31)
tum_lst = lapply(data_path, function(x){return(Read10X(data.dir = x))})
names(tum_lst) = snames
write_rds(tum_lst, "tum_lst.rds")


# Create seurat and add meta ---------------------------------------------
patient_id = tracker %>% dplyr::filter(samp_id %in% snames) %>% pull(patient_id)
source = tracker %>% dplyr::filter(samp_id %in% snames) %>% pull(source)
exp_proc = tracker %>% dplyr::filter(samp_id %in% snames) %>% pull(exp_proc)
procedure = tracker %>% dplyr::filter(samp_id %in% snames) %>% pull(procedure)
sampletype = tracker %>% dplyr::filter(samp_id %in% snames) %>% pull(sampletype)


tum = lapply(tum_lst, function(x){return(CreateSeuratObject(counts = x, project = "HBCA"))})
tum = mapply(function(x,y) {x@meta.data$sample_id = y;return(x)}, tum, snames)
tum = mapply(function(x,y) {x@meta.data$patient_id = y;return(x)}, tum, patient_id)
tum = mapply(function(x,y) {x@meta.data$source = y;return(x)}, tum, source)
tum = mapply(function(x,y) {x@meta.data$exp_proc = y;return(x)}, tum, exp_proc)
tum = mapply(function(x,y) {x@meta.data$procedure = y;return(x)}, tum, procedure)
tum = mapply(function(x,y) {x@meta.data$sampletype = y;return(x)}, tum, sampletype)

table(is.na(tum[[1]]@meta.data))
write_rds(tum, path = "tum.rds")

# MERGE Seurat ---------------------------------------------------------------
# merger sample 1 and 2
for (i in names(tum)) {
  tum[[i]] <- RenameCells(tum[[i]], add.cell.id = i)
}

# merge all the objects in the list
tum_merged <- purrr::reduce(tum, merge)

# write file
write_rds(tum_merged, path = "tum_merged_v1.rds")
tum_merged = read_rds(path = "/volumes/lab/users/tapsi/projects/bcap/analysis//01_read_files/tum_merged_v1.rds")

# FILTER ----------------

tum = tum_merged

# mito genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = GetAssayData(object = tum)), value = TRUE)
percent.mito <- Matrix::colSums(GetAssayData(object = tum)[mito.genes, ])/Matrix::colSums(GetAssayData(object = tum))
# rb genes
rb.genes <- grep(pattern = "^(RP)|(MRP)", x = rownames(x = GetAssayData(object = tum)), value = TRUE)
percent.rb <- Matrix::colSums(GetAssayData(object = tum)[rb.genes, ])/Matrix::colSums(GetAssayData(object = tum))
# ug genes
ug.prop = log10(tum@meta.data[,"nCount_RNA"])/log10(tum@meta.data[,"nFeature_RNA"])
# stress genes
stressgenes1 =  c("FOS", "CXCL12", "ZFP36", "FOSB", "DUSP1", "ATF3", "CXCL8", "CXCL3", "NR4A1", "JUNB", "HSP1A", "HSP1B") # from core list
stressgenes2 = c("ICAM1", "IRF", "ATF3", "HSPA6", "NFKBIA", "TNFAIP3", "DLX2", "NFATC2", "ZCTH12A") # https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1830-0/figures/4 b
stressgenes3 = c("DUSP", "FOS", "FOSB", "JUN", "JUNB", "IER2", "SOCS3", "ZFP36", "IGFBP7") # https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1830-0/figures/4 c
stressgenes1 = stressgenes1[stressgenes1 %in% rownames(x = GetAssayData(object = tum))]
stressgenes2 = stressgenes2[stressgenes2 %in% rownames(x = GetAssayData(object = tum))]
stressgenes3 = stressgenes3[stressgenes3 %in% rownames(x = GetAssayData(object = tum))]
percent.sg1 <- Matrix::colSums(GetAssayData(object = tum)[stressgenes1, ])/Matrix::colSums(GetAssayData(object = tum))
percent.sg2 <- Matrix::colSums(GetAssayData(object = tum)[stressgenes2, ])/Matrix::colSums(GetAssayData(object = tum))
percent.sg3 <- Matrix::colSums(GetAssayData(object = tum)[stressgenes3, ])/Matrix::colSums(GetAssayData(object = tum))


# add meta
tum <- AddMetaData(tum, metadata = percent.mito, col.name = "percent.mito")
tum <- AddMetaData(tum, metadata = ug.prop, col.name = "ug.prop")
tum <- AddMetaData(tum, metadata = percent.rb, col.name = "percent.rb")
tum <- AddMetaData(tum, metadata = percent.sg1, col.name = "percent.sg1")
tum <- AddMetaData(tum, metadata = percent.sg2, col.name = "percent.sg2")
tum <- AddMetaData(tum, metadata = percent.sg3, col.name = "percent.sg3")

tum$orig.ident = tum$sample_id
Idents(tum) <- "sample_id"

write_rds(tum, "tum_v0.rds")


message("saving plots ...4/6")
save_path =  "/volumes/lab/users/tapsi/projects/bcap/analysis/20201212_final2/01_read_files/"
png(file = glue("{save_path}/pre_filter1.png"), width=24,height=16, units = "in", res = 200)
print(VlnPlot(object = tum, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "ug.prop"), ncol = 2, pt.size = 0))
dev.off() 

png(file = glue("{save_path}/pre_filter2.png"), width=24,height=16, units = "in", res = 200)
print(VlnPlot(object = tum, features = c("percent.rb", "percent.sg1", "percent.sg2", "percent.sg3"), ncol = 2, pt.size = 0))
dev.off() 

# does not run in functions ------
length(tum$orig.ident) # 581872
# tum2 = func_filters(seu_object = tum, setident = "sample_id", save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20200210_final1/01_read_filesv2",
#                         gene_min = 150, gene_max = 7000, umi_min = 200, umi_max = 40000, mito_max = 0.15, rb_max = 0.50, ug_min = 0.25, ug_max = 1.5) # useful
# 
# tum0 = func_filters(seu_object = tum, setident = "sample_id", save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20200210_final1/01_read_filesv0",
#                     gene_min = 150, gene_max = 100000, umi_min = 200, umi_max = 10000000, mito_max = 0.15, rb_max = 0.50, ug_min = 0.25, ug_max = 1.5) # least stringent

tum3 = func_filters(seu_object = tum, setident = "sample_id", save_path = "/volumes/lab/users/tapsi/projects/bcap/analysis/20200210_final1/01_read_files/",
                    gene_min = 200, gene_max = 5000, umi_min = 500, umi_max = 20000, mito_max = 0.10, rb_max = 0.50, ug_min = 0.50, ug_max = 1.5) # most stringent

tum3 = subset(x = tum, subset = nCount_RNA > 500)
tum3 = subset(x = tum3, subset = nCount_RNA < 20000)
tum3 = subset(x = tum3, subset = nFeature_RNA > 200)
tum3 = subset(x = tum3, subset = nFeature_RNA < 5000)
tum3 = subset(x = tum3, subset = percent.mito < 0.10)
tum3 = subset(x = tum3, subset = percent.rb < 0.50)
tum3 = subset(x = tum3, subset = ug.prop < 1.5)
tum3 = subset(x = tum3, subset = ug.prop > 0.50)
write_rds(tum3, "tum_v1.rds")


message("saving plots ...5/6")
png(file = glue("{save_path}/post_filter1.png"), width=24,height=16, units = "in", res = 200)
print(VlnPlot(object = tum3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "ug.prop"), ncol = 2, pt.size = 0))
dev.off() 

png(file = glue("{save_path}/post_filter2.png"), width=24,height=16, units = "in", res = 200)
print(VlnPlot(object = tum3, features = c("percent.rb", "percent.sg1", "percent.sg2", "percent.sg3"), ncol = 2, pt.size = 0))
dev.off()


# qc -------------------------------------------------------------------------------------------------------------------
message("aggregating data ...2/6")
#tum3 = tum
tum_count_agg = aggregate(nCount_RNA~sample_id, tum3@meta.data, summary)
tum_RNA_agg = aggregate(nFeature_RNA~sample_id, tum3@meta.data, summary)
tum_percent.rb_agg = aggregate(percent.rb~sample_id, tum3@meta.data, summary)
tum_percent.mito_agg = aggregate(percent.mito~sample_id, tum3@meta.data, summary)
tum_percent.sg1_agg = aggregate(percent.sg1~sample_id, tum3@meta.data, summary)
tum_percent.sg2_agg = aggregate(percent.sg2~sample_id, tum3@meta.data, summary)
tum_percent.sg3_agg = aggregate(percent.sg3~sample_id, tum3@meta.data, summary)
tum_summary = cbind(tum_count_agg, tum_RNA_agg, tum_percent.rb_agg, tum_percent.mito_agg, tum_percent.sg1_agg, tum_percent.sg2_agg, tum_percent.sg3_agg)

# save the aggreate summary of each filter
message("saving csv ...3/6")
write.csv(tum_summary, glue("{save_path}/tum3_summary_sample_id.csv"))


# end -----------------



