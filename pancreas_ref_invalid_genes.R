

# Load required libraries 
required_packages <- c( 'Matrix', 'dplyr','rhdf5', 'Seurat', 'patchwork', 'SeuratData', 'ggplot2', 'SeuratDisk')

invisible(lapply(required_packages, library, character.only = T))

# obtain directories
wp <- "/data/cephfs-2/unmirrored/groups/krieger/users/spiese_c/visium/"
# wp <- "C:/Users/Eli/Desktop/Patho/"

ref_dir <- paste0(wp, "spot_deconvolution/celltype_deconvolution/References/")
peng_nr <- "PRJCA001063_besca"
peng_dir <- paste0(ref_dir,peng_nr, "/")
ig_dir <- paste0(ref_dir, "invalid_genes/")
gse_dir <- paste0(ig_dir, "gse111672/")

peng_dir <- "folder where peng seurat file  (StdWf1_PRJCA001063_CRC_besca2.annotated.h5ad) can be loaded"
ig_dir <- "folder where pk_all.rds can be loaded"
gse_dir <- paste0(ig_dir, "gse111672/")


########################
### pk_all reference ###
########################

# load data for all cells
# wget "https://zenodo.org/records/6024273/files/pk_all.rds"
pk_ref <- readRDS(paste0(ig_dir, "pk_all.rds"))
DefaultAssay(pk_ref)
# [1] "RNA"

dim(pk_ref)
# [1]  34746 136163

# subset to counts
pk_counts <- pk_ref@assays$RNA@counts
str(pk_counts)
# Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
# ..@ i       : int [1:305589707] 3 6 7 14 15 16 26 27 28 29 ...
# ..@ p       : int [1:136164] 0 3342 4316 5806 6907 10096 11536 12799 13971 15290 ...
# ..@ Dim     : int [1:2] 34746 136163
# ..@ Dimnames:List of 2
# .. ..$ : chr [1:34746] "LINC00115" "FAM41C" "SAMD11" "NOC2L" ...
# .. ..$ : chr [1:136163] "T1_AAACCTGAGATGTCGG" "T1_AAACGGGGTCATGCAT" "T1_AAAGATGCATGTTGAC" "T1_AAAGATGGTCGAGTTT" ...
# ..@ x       : num [1:305589707] 0.613 1.264 1.264 0.613 0.991 ...
# ..@ factors : list()

# check invalid features
invalid_features <- pk_ref[c(grep("-Mar$", rownames(pk_ref)), grep("-Sep$", rownames(pk_ref)), grep("-Dec$", rownames(pk_ref)), grep("<class 'str'>", rownames(pk_ref))),]
dim(invalid_features)
# [1]     29 136163

# exclude cells with no counts
invalid_features <- subset(invalid_features, subset=nFeature_RNA >0 & nCount_RNA >0)
dim(invalid_features)
# [1]    29 10674

# add info about cells with invalid counts to pk_ref
pk_ref[["cells_with_invalid_counts"]] <- colnames(pk_ref) %in% colnames(invalid_features)

table(pk_ref@meta.data$Project, pk_ref@meta.data$cells_with_invalid_counts)
#            FALSE  TRUE
# CA001063   50130  7293
# GSE111672    244  3381
# GSE154778   8000     0
# GSE155698  50943     0
# GSM4293555  4874     0
# OUGS       11298     0

library(RColorBrewer)
# pal_set <- colorRampPalette(brewer.pal(8,"Paired"))
png(paste0(ig_dir, "cells_with_invalid_counts.png"), 
    width = 35, height = 20, units = "cm", res = 600)
par(mfrow=c(2,1), mar=c(8,4,3,4))
barplot(table(pk_ref@meta.data$cells_with_invalid_counts, pk_ref@meta.data$Project), col= c("orange", "blue"),legend.text=T, horiz=F, args.legend=list(x=7.3, y=55000, bty = "n"), 
        las=1, main="Cells with invalid Gene Counts per Project")
barplot(table(pk_ref@meta.data$cells_with_invalid_counts[pk_ref@meta.data$Project %in% c("GSE111672", "CA001063")], pk_ref@meta.data$Patient2[pk_ref@meta.data$Project %in% c("GSE111672", "CA001063")]), col= c("orange", "blue"),legend.text=T, horiz=F, args.legend=list(x=45, y=4000, bty = "n"), 
        las=2, main="Cells in Studies with invalid Gene Counts per Sample")
dev.off()


DefaultAssay(pk_ref)
# # [1] "RNA"

# create dim plots
library(patchwork)
p1 <- DimPlot(pk_ref, group.by = "Cell_type", label = TRUE, label.size=3, repel = T, pt.size = 0.2) + ggtitle("Cell Type")
p2 <-  DimPlot(pk_ref, group.by = "Type",label = TRUE, label.size=3, repel = T, pt.size = 0.2) + ggtitle("Disease State")
p3 <-  DimPlot(pk_ref, group.by = "Project",label = TRUE, label.size=3, repel = T, pt.size = 0.2) + ggtitle("Project ID")
p4 <-  DimPlot(pk_ref, group.by = "cells_with_invalid_counts",label = TRUE, label.size=3, repel = T, pt.size = 0.2) + ggtitle("Cells with invalid gene counts")
plot_patched <- wrap_plots(p1,p2,p3,p4, ncol=2) + plot_annotation("UMAP for pk_all.rds") + theme(plot.title = element_text(size = 18))
ggsave(paste0(ig_dir, "dim_plot_umap_integrated.png"), plot = plot_patched, width = 35, height =25, units = "cm")

table(pk_ref@meta.data$predicted.id)
# Acinar cell     B cell    Ductal cell type 1    Ductal cell type 2  Endocrine cell   Endothelial cell    Fibroblast cell    Macrophage cell  Stellate cell    T cell
# 3857            5838      11605                 33607               1298              10059              11003              23703            7624             27569

table(pk_ref@meta.data$Project)
# CA001063  GSE111672  GSE154778  GSE155698 GSM4293555   OUGS
# 57423       3625       8000      50943       4874      11298
# hg19        hg38       hg19      hg19        hg19      hg19          


# replace class genes from peng
# using raw data would be advantageous to avoid subsetting & normalizing twice (but class genes are also present in raw data)
class_gene_dict <- read.csv(paste0(ig_dir, "class_gene_dict.csv"))
class_gene_dict
#           ENSEMBL       SYMBOL      class_gene
# 1 ENSG00000249201       TERLR1   <class 'str'>
# 2 ENSG00000233937 CTC-338M12.4 <class 'str'>-1
# 3 ENSG00000237310    LINC03011 <class 'str'>-2
# 4 ENSG00000229180 GS1-124K5.11 <class 'str'>-3
# 5 ENSG00000249825    THBS4-AS1 <class 'str'>-4
# 6 ENSG00000245857   GS1-24F4.2 <class 'str'>-5
# 7 ENSG00000253978    TENM2-AS1 <class 'str'>-6
# 8 ENSG00000261729 GS1-204I12.4 <class 'str'>-7
# 9 ENSG00000225833      ACE2-DT <class 'str'>-8

# replace month genes from gse111672
# using raw data would be advantageous to avoid invalid genes as well as subsetting & normalizing twice
month_gene_dict <- read.csv(paste0(ig_dir, "month_gene_dict.csv"))
month_gene_dict
#     Month    Gene
# 1   1-Dec    DEC1
# 2   1-Mar   MARC1
# 3   2-Mar   MARC2
# 4   1-Mar  MARCH1
# 5  10-Mar MARCH10
# 6  11-Mar MARCH11
# 7   2-Mar  MARCH2
# 8   3-Mar  MARCH3
# 9   4-Mar  MARCH4
# 10  5-Mar  MARCH5
# 11  6-Mar  MARCH6
# 12  7-Mar  MARCH7
# 13  8-Mar  MARCH8
# 14  9-Mar  MARCH9
# 15  1-Sep   SEPT1
# 16 10-Sep  SEPT10
# 17 11-Sep  SEPT11
# 18 12-Sep  SEPT12
# 19 14-Sep  SEPT14
# 20  2-Sep   SEPT2
# 21  3-Sep   SEPT3
# 22  4-Sep   SEPT4
# 23  5-Sep   SEPT5
# 24  6-Sep   SEPT6
# 25  7-Sep   SEPT7
# 26  8-Sep   SEPT8
# 27  9-Sep   SEPT9

# to replace class genes in pk_ref data
class_gene_dict$SYMBOL[class_gene_dict$SYMBOL %in% rownames(pk_ref)]
# [1] "CTC-338M12.4" "GS1-124K5.11" "GS1-24F4.2"   "GS1-204I12.4"

setdiff(class_gene_dict$class_gene, rownames(pk_ref))
# character(0)
setdiff(class_gene_dict$class_gene, rownames(invalid_features))
# character(0)

# replace month genes in pk_ref data
month_gene_dict$Gene[month_gene_dict$Gene %in% rownames(pk_ref)]
# [1] "DEC1"    "MARC1"   "MARC2"   "MARCH1"  "MARCH10" "MARCH11" "MARCH2"
# [8] "MARCH3"  "MARCH4"  "MARCH5"  "MARCH6"  "MARCH7"  "MARCH8"  "MARCH9"
# [15] "SEPT1"   "SEPT10"  "SEPT11"  "SEPT12"  "SEPT14"  "SEPT2"   "SEPT3"
# [22] "SEPT4"   "SEPT5"   "SEPT6"   "SEPT7"   "SEPT8"   "SEPT9"

all(month_gene_dict$Gene %in% rownames(pk_ref))
# [1] TRUE

setdiff(month_gene_dict$Month, rownames(pk_ref))
# [1] "1-Dec"  "10-Mar" "11-Mar" "12-Sep" "5-Sep"
setdiff(month_gene_dict$Month, rownames(invalid_features))
# [1] "1-Dec"  "10-Mar" "11-Mar" "12-Sep" "5-Sep"

any(duplicated(rownames(pk_ref))) || any(duplicated(rownames(invalid_features)))
# [1] FALSE


# gene_map: named vector where names are invalid features, values are valid gene symbols
gene_map <- c(class_gene_dict$SYMBOL, month_gene_dict$Gene)
names(gene_map) <-  c(class_gene_dict$class_gene, month_gene_dict$Month)
gene_map
# <class 'str'> <class 'str'>-1 <class 'str'>-2 <class 'str'>-3 <class 'str'>-4
# "TERLR1"  "CTC-338M12.4"     "LINC03011"  "GS1-124K5.11"     "THBS4-AS1"
# <class 'str'>-5 <class 'str'>-6 <class 'str'>-7 <class 'str'>-8           1-Dec
# "GS1-24F4.2"     "TENM2-AS1"  "GS1-204I12.4"       "ACE2-DT"          "DEC1"
# 1-Mar           2-Mar           1-Mar          10-Mar          11-Mar
# "MARC1"         "MARC2"        "MARCH1"       "MARCH10"       "MARCH11"
# 2-Mar           3-Mar           4-Mar           5-Mar           6-Mar
# "MARCH2"        "MARCH3"        "MARCH4"        "MARCH5"        "MARCH6"
# 7-Mar           8-Mar           9-Mar           1-Sep          10-Sep
# "MARCH7"        "MARCH8"        "MARCH9"         "SEPT1"        "SEPT10"
# 11-Sep          12-Sep          14-Sep           2-Sep           3-Sep
# "SEPT11"        "SEPT12"        "SEPT14"         "SEPT2"         "SEPT3"
# 4-Sep           5-Sep           6-Sep           7-Sep           8-Sep
# "SEPT4"         "SEPT5"         "SEPT6"         "SEPT7"         "SEPT8"
# 9-Sep
# "SEPT9"


# replace invalid genes in invalid features
rownames(pk_ref)[c(grep("-Mar", rownames(pk_ref)), grep("-Sep", rownames(pk_ref)), grep("-Dec", rownames(pk_ref)), grep("<class 'str'>", rownames(pk_ref)))]
# [1] "1-Mar"           "2-Mar"           "3-Mar"           "4-Mar"
# [5] "5-Mar"           "6-Mar"           "7-Mar"           "8-Mar"
# [9] "9-Mar"           "1-Sep"           "10-Sep"          "11-Sep"
# [13] "14-Sep"          "2-Sep"           "3-Sep"           "4-Sep"
# [17] "6-Sep"           "7-Sep"           "8-Sep"           "9-Sep"
# [21] "<class 'str'>"   "<class 'str'>-1" "<class 'str'>-2" "<class 'str'>-3"
# [25] "<class 'str'>-4" "<class 'str'>-5" "<class 'str'>-6" "<class 'str'>-7"
# [29] "<class 'str'>-8"

invalid_features <- pk_ref[c(grep("-Mar$", rownames(pk_ref)), grep("-Sep$", rownames(pk_ref)), grep("-Dec$", rownames(pk_ref)), grep("<class 'str'>", rownames(pk_ref))),]
dim(invalid_features)
# [1]     29 136163


##########################################################################################################
rownames(GSE111672_A) <- make.names(GSE111672_A$Genes)
# function (names, unique = FALSE, allow_ = TRUE) 
# {
#   names <- as.character(names)
#   names2 <- .Internal(make.names(names, allow_))
#   if (unique) {
#     o <- order(names != names2)
#     names2[o] <- make.unique(names2[o])
#   }
#   names2
# }
# by default make.names(unique=FALSE)
GSE111672_A <- distinct(GSE111672_A,Genes ,.keep_all = TRUE)
# distinct() keeps the first occurrence of each duplicate by default.
# That means:
# When "1-Mar" was encountered, distinct() kept the first raw gene that mapped to "1-Mar" = MARC1, and dropped MARCH1.
# Same for "2-Mar" = MARC2, dropping MARCH2.
# genes_filtered[duplicated(genes_filtered)] gave "1-Mar", "2-Mar", showing both MARC and MARCH collapsed.
# only_raw shows MARC1 and MARC2 occur before MARCH1 and MARCH2 in the original order.
# Therefore, "1-Mar" and "2-Mar" must correspond to MARC1 and MARC2, since distinct() kept the first.
##########################################################################################################

# subset gene map
# exclude MARCH1 & MARCH2 (duplicated 1-Mar & 2-Mar)
gene_map <- gene_map[!duplicated(names(gene_map))]
# only keep entries which can be found in pk_ref
gene_map <- gene_map[names(gene_map) %in% rownames(pk_ref)]
gene_map
# <class 'str'> <class 'str'>-1 <class 'str'>-2 <class 'str'>-3 <class 'str'>-4
# "TERLR1"  "CTC-338M12.4"     "LINC03011"  "GS1-124K5.11"     "THBS4-AS1"
# <class 'str'>-5 <class 'str'>-6 <class 'str'>-7 <class 'str'>-8           1-Mar
# "GS1-24F4.2"     "TENM2-AS1"  "GS1-204I12.4"       "ACE2-DT"         "MARC1"
# 2-Mar           3-Mar           4-Mar           5-Mar           6-Mar
# "MARC2"        "MARCH3"        "MARCH4"        "MARCH5"        "MARCH6"
# 7-Mar           8-Mar           9-Mar           1-Sep          10-Sep
# "MARCH7"        "MARCH8"        "MARCH9"         "SEPT1"        "SEPT10"
# 11-Sep          14-Sep           2-Sep           3-Sep           4-Sep
# "SEPT11"        "SEPT14"         "SEPT2"         "SEPT3"         "SEPT4"
# 6-Sep           7-Sep           8-Sep           9-Sep
# "SEPT6"         "SEPT7"         "SEPT8"         "SEPT9"

# to replace invalid genes in invalid_features
invalid_names <- rownames(invalid_features)
rownames(invalid_features@assays$RNA@counts) <- gene_map[match(invalid_names, names(gene_map))]
rownames(invalid_features@assays$RNA@data) <- gene_map[match(invalid_names, names(gene_map))]
rownames(invalid_features@assays$RNA@meta.features) <- gene_map[match(invalid_names, names(gene_map))]
dim(invalid_features)
# [1]     29 136163

# exclude invalid features in pk_ref
valid_features <- pk_ref[-c(grep("-Mar$", rownames(pk_ref)), grep("-Sep$", rownames(pk_ref)), grep("-Dec$", rownames(pk_ref)), grep("<class 'str'>", rownames(pk_ref))),]
dim(valid_features)
# [1]  34717 136163

# 24 genes are already present in valid features -> counts (embedding) would change
intersect(rownames(invalid_features), rownames(valid_features))
# [1] "MARC1"        "MARC2"        "MARCH3"       "MARCH4"       "MARCH5"
# [6] "MARCH6"       "MARCH7"       "MARCH8"       "MARCH9"       "SEPT1"
# [11] "SEPT10"       "SEPT11"       "SEPT14"       "SEPT2"        "SEPT3"
# [16] "SEPT4"        "SEPT6"        "SEPT7"        "SEPT8"        "SEPT9"
# [21] "CTC-338M12.4" "GS1-124K5.11" "GS1-24F4.2"   "GS1-204I12.4"
length(intersect(rownames(invalid_features), rownames(valid_features)))
# [1] 24
length(intersect(rownames(invalid_features), rownames(valid_features)))/length(rownames(invalid_features))
# [1] 0.8275862

# only 5 new genes in invalid features
setdiff(rownames(invalid_features), rownames(valid_features))
# [1] "TERLR1"    "LINC03011" "THBS4-AS1" "TENM2-AS1" "ACE2-DT"

# -> re-run using raw matrices???

################################################################################################################################################################################

###########################
### PRJCA001063 | BESCA ###
###########################

# replace class genes from peng besca
# using raw data would be advantageous to avoid subsetting & normalizing twice (but class genes are also present in raw data)
# peng raw
# wget "https://zenodo.org/records/3969339/files/StdWf1_PRJCA001063_CRC_besca2.raw.h5ad"
# Convert(paste0(peng_dir, "StdWf1_PRJCA001063_CRC_besca2.raw.h5ad"), dest="h5seurat", overwrite = TRUE)
# peng_raw <- LoadH5Seurat(paste0(peng_dir, "StdWf1_PRJCA001063_CRC_besca2.raw.h5seurat"))

# https://www.nature.com/articles/s41422-019-0195-y
# https://pmc.ncbi.nlm.nih.gov/articles/PMC6796938/
# https://zenodo.org/record/3969339#.YgORVfhUtaY
# https://www.nature.com/articles/s41422-019-0195-y#Sec11

# get, read & adjust data
library(rhdf5)
# wget "https://zenodo.org/records/3969339/files/StdWf1_PRJCA001063_CRC_besca2.annotated.h5ad"

list.files(paste0(peng_dir, ""))

#Convert(paste0(peng_dir, "StdWf1_PRJCA001063_CRC_besca2.annotated.h5ad"), dest="h5seurat", overwrite = TRUE)
peng_seurat <- LoadH5Seurat(paste0(peng_dir, "StdWf1_PRJCA001063_CRC_besca2.annotated.h5seurat"))
h5ls(paste0(peng_dir, "StdWf1_PRJCA001063_CRC_besca2.annotated.h5ad"))
peng_X <- h5read(paste0(peng_dir, "StdWf1_PRJCA001063_CRC_besca2.annotated.h5ad"), name = "raw")
peng_var <- peng_X$var

class_genes <- peng_var$ENSEMBL[grep("class", peng_var$index)]
# replace class genes by proper gene symbols in peng seurat objects

# all(rownames(peng_seurat) %in% rownames(peng_raw))
# [1] TRUE

# no proper names available in mart37
# (request on 08/30/2025)
library(biomaRt)
# options(biomaRt.cache = "/data/cephfs-1/scratch/groups/krieger/users/spiese_c/.cache/biomaRt/")

# load hg19 / grch37 mart data
# file.remove(list.files("/data/cephfs-1/scratch/groups/krieger/users/spiese_c/.cache/biomaRt/", full.names = TRUE))
ensembl_hg19 <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://grch37.ensembl.org")
class_gene_names <- biomaRt::getBM(attributes = c("chromosome_name", "ensembl_gene_id", "hgnc_symbol", "band"), filters = "ensembl_gene_id", values=class_genes, mart = ensembl_hg19)
class_gene_names
# chromosome_name ensembl_gene_id hgnc_symbol   band
# 1               X ENSG00000225833          NA  p22.2
# 2               7 ENSG00000229180          NA q11.21
# 3               5 ENSG00000233937          NA  q35.3
# 4               7 ENSG00000237310          NA q11.21
# 5               8 ENSG00000245857          NA  p23.1
# 6               5 ENSG00000249201          NA p15.33
# 7               5 ENSG00000249825          NA  q14.1
# 8               5 ENSG00000253978          NA    q34
# 9               1 ENSG00000261729          NA  q25.3


library(org.Hs.eg.db)  # use org.Hs.eg.db since mart37 provided not any symbol
class_gene_names <- AnnotationDbi::select(org.Hs.eg.db, keys = class_genes, column = c("ENTREZID", "MAP", "SYMBOL"), keytype = c("ENSEMBL"))
# 'select()' returned 1:1 mapping between keys and columns
# (request on 08/30/2025)
# replace missing symbol manually
class_gene_names
# 'select()' returned 1:1 mapping between keys and columns
# ENSEMBL  ENTREZID     MAP       SYMBOL
# 1 ENSG00000249201 101928857 5p15.33       TERLR1
# 2 ENSG00000233937      <NA>    <NA>         <NA> 
# (request on 08/30/2025)
# CTC-338M12.4 (https://www.bgee.org/gene/ENSG00000233937?pageNumber=2)
# CTC-338M12.4 (https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000233937;r=5:180673523-180699168)
# 3 ENSG00000237310 100289098 7q11.21    LINC03011
# 4 ENSG00000229180      <NA>    <NA>         <NA> 
# (request on 08/30/2025)
# GS1-124K5.11 (https://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000229180;r=7:65991075-66007536;t=ENST00000449307)
# GS1-124K5.11 (http://angiogenes.uni-frankfurt.de/gene/ENSG00000229180)
# 5 ENSG00000249825 101929215  5q14.1    THBS4-AS1
# 6 ENSG00000245857 100652791  8p23.1   GS1-24F4.2
# 7 ENSG00000253978 101927862    5q34    TENM2-AS1
# 8 ENSG00000261729 105371653  1q25.3 GS1-204I12.4
# 9 ENSG00000225833 104798195  Xp22.2      ACE2-DT

class_gene_names$SYMBOL[class_gene_names$ENSEMBL == "ENSG00000233937"] <- "CTC-338M12.4"
class_gene_names$SYMBOL[class_gene_names$ENSEMBL == "ENSG00000229180"] <- "GS1-124K5.11"
class_gene_dict <- class_gene_names[,c("ENSEMBL", "SYMBOL")]
class_gene_dict$class_gene <- peng_var$index[grepl("class", peng_var$index)][match(class_gene_dict$ENSEMBL, peng_var$ENSEMBL[grepl("class", peng_var$index)])]
print(class_gene_dict)
# ENSEMBL       SYMBOL      class_gene
# 1 ENSG00000249201       TERLR1   <class 'str'>
# 2 ENSG00000233937 CTC-338M12.4 <class 'str'>-1
# 3 ENSG00000237310    LINC03011 <class 'str'>-2
# 4 ENSG00000229180 GS1-124K5.11 <class 'str'>-3
# 5 ENSG00000249825    THBS4-AS1 <class 'str'>-4
# 6 ENSG00000245857   GS1-24F4.2 <class 'str'>-5
# 7 ENSG00000253978    TENM2-AS1 <class 'str'>-6
# 8 ENSG00000261729 GS1-204I12.4 <class 'str'>-7
# 9 ENSG00000225833      ACE2-DT <class 'str'>-8

# any(class_gene_dict$SYMBOL %in% rownames(peng_raw))
# [1] FALSE

write.csv(class_gene_dict, paste0(ig_dir, "class_gene_dict.csv"), row.names=FALSE)

################################################################################################################################################################################

is_integer_0L <- function(x){
  all(is.integer(x) & length(x) == 0L)}

# to replace class genes in peng_seurat data
if (all(class_gene_dict$ENSEMBL == peng_var$ENSEMBL[grep("class", peng_var$index)])){
  print("Replacing class genes by proper gene symbols")
  # add to gene names
  peng_var$index[grepl("class", peng_var$index)] <- class_gene_names$SYMBOL[match(peng_var$ENSEMBL[grepl("class", peng_var$index)], class_gene_names$ENSEMBL)]
  peng_var$index <- make.unique(peng_var$index)
  rownames(peng_seurat@assays$RNA@counts) <- peng_var$index
  rownames(peng_seurat@assays$RNA@data) <- peng_var$index
  if (all(peng_seurat@assays$RNA@meta.features$ENSEMBL == peng_var$ENSEMBL)){
    peng_seurat@assays$RNA@meta.features$Gene_Name <- peng_var$index
    rownames(peng_seurat@assays$RNA@meta.features) <- peng_seurat@assays$RNA@meta.features$Gene_Name}
  if (is_integer_0L(grep("class", rownames(peng_seurat))) & all(rownames(peng_seurat) == peng_var$index)==TRUE){
    print("Replaced all class gene names.")}}

dim(peng_seurat)
# [1] 17004 57423

any(grepl("class", rownames(peng_seurat)))
# [1] FALSE

################################################################################################################################################################################

#################
### GSE111672 ###
#################

# replace month genes from gse111672
# using raw matrix data would be advantageous to avoid subsetting & normalizing twice, as well as no month genes are in raw matrix files

# get, read & adjust data
# use wget or download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111672
# wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE111672&format=file&file=GSE111672%5FPDAC%2DA%2Dindrop%2Dfiltered%2DexpMat%2Etxt%2Egz"
# wget "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE111672&format=file&file=GSE111672%5FPDAC%2DB%2Dindrop%2Dfiltered%2DexpMat%2Etxt%2Egz"

GSE111672_A <- read.delim(paste0(gse_dir, "GSE111672_PDAC-A-indrop-filtered-expMat.txt.gz"), sep = "\t")
dim(GSE111672_A)
# [1] 19738  1927
GSE111672_B <- read.delim(paste0(gse_dir,"GSE111672_PDAC-B-indrop-filtered-expMat.txt.gz"), sep = "\t")
dim(GSE111672_B)
# [1] 19738  1734

GSE111672_A$Genes[c(grep("-Dec$", GSE111672_A$Genes), grep("-Mar$", GSE111672_A$Genes), grep("-Sep$", GSE111672_A$Genes), grep("<class 'str'>", GSE111672_A$Genes))]
# [1] "1-Dec"  "1-Mar"  "2-Mar"  "1-Mar"  "10-Mar" "11-Mar" "2-Mar"  "3-Mar"  "4-Mar"  "5-Mar"  "6-Mar"  "7-Mar"  "8-Mar"  "9-Mar"  "1-Sep"  "10-Sep" "11-Sep" "12-Sep" "14-Sep" "2-Sep"  "3-Sep"  "4-Sep"  "5-Sep"  "6-Sep"  "7-Sep"  "8-Sep"  "9-Sep" 

c(grep("-Dec$", GSE111672_A$Genes), grep("-Mar$", GSE111672_A$Genes), grep("-Sep$", GSE111672_A$Genes), grep("<class 'str'>", GSE111672_A$Genes))
# [1]  4772  9945  9946  9947  9948  9949  9950  9951  9952  9953  9954  9955  9956  9957 15064 15065 15066 15067 15068 15069 15070 15071 15072 15073 15074 15075 15076

all(GSE111672_A$Genes == GSE111672_B$Genes)
# [1] TRUE

length(unique(GSE111672_A$Genes))
# [1] 19736

GSE111672_A$Genes[duplicated(GSE111672_A$Genes)]
# [1] "1-Mar" "2-Mar"

genes_filtered <- GSE111672_A$Genes
str(genes_filtered)
# chr [1:19738] "A1BG" "A1CF" "A2M" "A2ML1" "A3GALT2" "A4GALT" "A4GNT" ...

################################################################################################################################################################################

# load raw matrix files
library(GEOquery)
# load gse info
geo_nr <- "GSE111672"
# gse <- getGEO(GEO = geo_nr, GSEMatrix = FALSE)
# alternative on the cluster, if ssl refuses to connect
# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111672/soft/GSE111672_family.soft.gz
gse <- getGEO(filename = paste0(gse_dir, "GSE111672_family.soft.gz"))
  
names(gse@gsms)
# [1] "GSM3036909" "GSM3036910" "GSM3036911" "GSM3405527" "GSM3405528" "GSM3405529" "GSM3405530" "GSM3405531" "GSM3405532" "GSM3405533" "GSM3405534" "GSM4100717" "GSM4100718" "GSM4100719" "GSM4100720"
#[16] "GSM4100721" "GSM4100722" "GSM4100723" "GSM4100724" "GSM4100725" "GSM4100726" "GSM4100727" "GSM4100728"

# 23 data sets
# 13 scRNA data sets
# 10 st data sets
# Obtain GEO-Sample Characteristics
sample_characteristics <- list()
for (x in names(gse@gsms)){
  index_genome <- grep("Genome_build", gse@gsms[[x]]@header$data_processing)
  sample_characteristics[[x]] <-  c(gse@gsms[[x]]@header$title, gse@gsms[[x]]@header$characteristics_ch1, gse@gsms[[x]]@header$data_processing[index_genome])}

df <- data.frame(matrix(nrow = length(sample_characteristics), ncol=4))
rownames(df) <- unlist(names(sample_characteristics))
colnames(df) <- c( "sample", "tissue", "tissue preparation", "Genome_build")
for (x in 1:length(sample_characteristics)){
  for (y in 1:ncol(df)){
    df[x,y] <- gsub(paste0(colnames(df)[y], ": "), "", sample_characteristics[[x]][y])}}

df$tissue <- "PDAC"
sample_info <- gsub("PDAC-", "", df$sample)
tech_info <- gsub("^[A-G]\ ", "", sample_info )
sample_info <- gsub("\ .*$", "", sample_info )
df$sample <- tech_info
colnames(df)[1] <- "technology"
df <- cbind("sample"=sample_info, df)

df
#            sample technology tissue               tissue preparation Genome_build
# GSM3036909      A    inDrop1   PDAC single-cell suspension of tissue         hg38
# GSM3036910      A    inDrop2   PDAC single-cell suspension of tissue         hg38
# GSM3036911      A        ST1   PDAC              tissue cryosections         hg38
# GSM3405527      A    inDrop3   PDAC single-cell suspension of tissue         hg38
# GSM3405528      A    inDrop4   PDAC single-cell suspension of tissue         hg38
# GSM3405529      A    inDrop5   PDAC single-cell suspension of tissue         hg38
# GSM3405530      A    inDrop6   PDAC single-cell suspension of tissue         hg38
# GSM3405531      B    inDrop1   PDAC single-cell suspension of tissue         hg38
# GSM3405532      B    inDrop2   PDAC single-cell suspension of tissue         hg38
# GSM3405533      B    inDrop3   PDAC single-cell suspension of tissue         hg38
# GSM3405534      B        ST1   PDAC              tissue cryosections         hg38
# GSM4100717      C    inDrop1   PDAC single-cell suspension of tissue         hg38
# GSM4100718      C    inDrop2   PDAC single-cell suspension of tissue         hg38
# GSM4100719      C    inDrop3   PDAC single-cell suspension of tissue         hg38
# GSM4100720      C    inDrop4   PDAC single-cell suspension of tissue         hg38
# GSM4100721      A        ST2   PDAC              tissue cryosections         hg38
# GSM4100722      A        ST3   PDAC              tissue cryosections         hg38
# GSM4100723      B        ST2   PDAC              tissue cryosections         hg38
# GSM4100724      B        ST3   PDAC              tissue cryosections         hg38
# GSM4100725      D        ST1   PDAC              tissue cryosections         hg38
# GSM4100726      E        ST1   PDAC              tissue cryosections         hg38
# GSM4100727      F        ST1   PDAC              tissue cryosections         hg38
# GSM4100728      G        ST1   PDAC              tissue cryosections         hg38


# subset to scrna data
in_drop <- list()
for (x in 1:length(names(gse@gsms))){
  in_drop[[names(gse@gsms)[x]]] <-  gse@gsms[[x]]@header$title}

# exclude ST / include inDrop only
in_drop <- in_drop[grep("inDrop", in_drop)]

# Obtain GEO-Sample-URLs
geo_url <- list()
for (x in names(in_drop)){
  geo_url[[x]] <- gse@gsms[[x]]@header$supplementary_file_1}

ifelse(!dir.exists(paste0(gse_dir, geo_nr, "_data/")),dir.create(paste0(gse_dir, geo_nr, "_data/")), "Directory already exists.")

# download data
for (x in 1:length(geo_url)){
  destfile <- paste0(gse_dir, geo_nr, "_data/", names(geo_url)[x], ".tsv.gz")
  download.file(geo_url[[x]], method="wget", destfile= destfile, mode = "wb")}

# read files
files_to_read <- list.files(path = paste0(gse_dir, geo_nr, "_data/"),pattern = "\\.tsv.gz$",full.names = T)
fm <- lapply(files_to_read,function(x) {
  read.delim(file = x, sep = '\t', header = TRUE,stringsAsFactors = FALSE)})
unique(lapply(fm, dim))
# [[1]]
# [1] 19831  2001
# 
# [[2]]
# [1] 19831  2501
# 
# [[3]]
# [1] 19936  2001

names(fm) <- names(in_drop)

length(fm)
# [1] 13

# no invalid gene names
for (i in seq_along(fm)){
 print(fm[[i]]$Genes[c(grep("-Mar$", fm[[i]]$Genes), grep("-Sep$", fm[[i]]$Genes), grep("<class 'str'>", fm[[i]]$Genes))])}
# character(0)
# character(0)
# character(0)
# character(0)
# character(0)
# character(0)
# character(0)
# character(0)
# character(0)
# character(0)
# character(0)
# character(0)
# character(0)

# all gene names are unique
for (i in seq_along(fm)){
  print(length(fm[[i]]$Genes) == length(unique(fm[[i]]$Genes)))}
# [1] TRUE
# [1] TRUE
# [1] TRUE
# [1] TRUE
# [1] TRUE
# [1] TRUE
# [1] TRUE
# [1] TRUE
# [1] TRUE
# [1] TRUE
# [1] TRUE
# [1] TRUE
# [1] TRUE

# get all gene names across all samples
all_genes <- c()
for (i in seq_along(fm)){
  all_genes <- unique(c(all_genes, fm[[i]]$Genes))}
str(all_genes)
# chr [1:20330] "A1BG" "A1CF" "A2M" "A2ML1" "A3GALT2" "A4GALT" "A4GNT" "AAAS" "AACS" "AADAC" "AADACL2" "AADACL3" "AADACL4" "AADAT" "AAED1" "AAGAB" "AAK1" "AAMDC" ...

all_genes[c(grep("-Mar$", all_genes), grep("-Sep$", all_genes), grep("<class 'str'>", all_genes))]
# character(0)

for (i in seq_along(fm)){
  print(all(fm[[i]]$Genes %in% all_genes))}
# [1] TRUE
# [1] TRUE
# [1] TRUE
# [1] TRUE
# [1] TRUE
# [1] TRUE
# [1] TRUE
# [1] TRUE
# [1] TRUE
# [1] TRUE
# [1] TRUE
# [1] TRUE
# [1] TRUE

############################################################################################################################################################

# check genes between filtered matrices (GSE111672_A / GSE111672_B -> genes_filtered) and raw matrices (fm -> all_genes)
only_filtered <- setdiff(genes_filtered, all_genes)
only_filtered
# [1] "1-Dec"  "1-Mar"  "2-Mar"  "10-Mar" "11-Mar" "3-Mar"  "4-Mar"  "5-Mar"  "6-Mar"  "7-Mar"  "8-Mar"  "9-Mar"  "1-Sep"  "10-Sep" "11-Sep" "12-Sep" "14-Sep" "2-Sep"  "3-Sep"  "4-Sep"  "5-Sep"  "6-Sep"  "7-Sep"  "8-Sep"  "9-Sep"

only_raw <- setdiff(all_genes, genes_filtered)
str(only_raw)
# chr [1:619] "DEC1" "FAU" "MARC1" "MARC2" "MARCH1" "MARCH10" "MARCH11" "MARCH2" "MARCH3" "MARCH4" "MARCH5" "MARCH6" "MARCH7" "MARCH8" "MARCH9" "MT-ATP6" "MT-ATP8" "MT-CO1" "MT-CO2" "MT-CO3" "MT-CYB" "MT-ND1" "MT-ND2" "MT-ND3" "MT-ND4" "MT-ND4L" ...
only_raw[grepl("DEC|MARC[[:digit:]]|MARCH|SEP", only_raw)]
# [1] "DEC1"    "MARC1"   "MARC2"   "MARCH1"  "MARCH10" "MARCH11" "MARCH2"  "MARCH3"  "MARCH4"  "MARCH5"  "MARCH6"  "MARCH7"  "MARCH8"  "MARCH9"  "SEPT1"   "SEPT10"  "SEPT11"  "SEPT12"  "SEPT14"  "SEPT2"   "SEPT3"   "SEPT4"   "SEPT5"   "SEPT6"  
# [25] "SEPT7"   "SEPT8"   "SEPT9"

# with MARC1 & MARCH1 -> 1-Mar
# and MARC2 & MARCH2 -> 2-Mar
genes_filtered[duplicated(genes_filtered)]
# [1] "1-Mar" "2-Mar"

genes_filtered[genes_filtered %in% only_filtered]
# [1] "1-Dec"  "1-Mar"  "2-Mar"  "1-Mar"  "10-Mar" "11-Mar" "2-Mar"  "3-Mar"  "4-Mar"  "5-Mar"  "6-Mar"  "7-Mar"  "8-Mar"  "9-Mar"  "1-Sep"  "10-Sep" "11-Sep" "12-Sep" "14-Sep" "2-Sep"  "3-Sep"  "4-Sep"  "5-Sep"  "6-Sep"  "7-Sep"  "8-Sep"  "9-Sep" 

month_gene_dict <- data.frame("Month" = genes_filtered[genes_filtered %in% only_filtered], "Gene" = only_raw[grepl("DEC|MARC[[:digit:]]|MARCH|SEP", only_raw)])
month_gene_dict
#     Month    Gene
# 1   1-Dec    DEC1
# 2   1-Mar   MARC1
# 3   2-Mar   MARC2
# 4   1-Mar  MARCH1
# 5  10-Mar MARCH10
# 6  11-Mar MARCH11
# 7   2-Mar  MARCH2
# 8   3-Mar  MARCH3
# 9   4-Mar  MARCH4
# 10  5-Mar  MARCH5
# 11  6-Mar  MARCH6
# 12  7-Mar  MARCH7
# 13  8-Mar  MARCH8
# 14  9-Mar  MARCH9
# 15  1-Sep   SEPT1
# 16 10-Sep  SEPT10
# 17 11-Sep  SEPT11
# 18 12-Sep  SEPT12
# 19 14-Sep  SEPT14
# 20  2-Sep   SEPT2
# 21  3-Sep   SEPT3
# 22  4-Sep   SEPT4
# 23  5-Sep   SEPT5
# 24  6-Sep   SEPT6
# 25  7-Sep   SEPT7
# 26  8-Sep   SEPT8
# 27  9-Sep   SEPT9

write.csv(month_gene_dict, paste0(ig_dir, "month_gene_dict.csv"), row.names = FALSE)

############################################################################################################################################################

# to replace month genes in GSE filtered matrix data
GSE111672_A$Genes[GSE111672_A$Genes %in% only_filtered] <- only_raw[grepl("DEC|MARC[[:digit:]]|MARCH|SEP", only_raw)]
GSE111672_B$Genes[GSE111672_B$Genes %in% only_filtered] <- only_raw[grepl("DEC|MARC[[:digit:]]|MARCH|SEP", only_raw)]



