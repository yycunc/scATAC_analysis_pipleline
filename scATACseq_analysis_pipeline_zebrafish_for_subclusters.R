library(SnapATAC)
snap.files = c("day0_scATAC/day0_scATAC/day0_with_subtype.snap", "day2_scATAC/day2_scATAC/day2_with_subtype.snap",
               "day7_scATAC/day7_scATAC/day7_with_subtype.snap", "day14_scATAC/day14_scATAC/day14_with_subtype.snap")
sample.names = c("day0", "day2", "day7", "day14")
barcode.files = c("day0_scATAC/day0_scATAC/singlecell.csv", "day2_scATAC/day2_scATAC/singlecell.csv",
                  "day7_scATAC/day7_scATAC/singlecell.csv", "day14_scATAC/day14_scATAC/singlecell.csv")
ZF_with_subtype.sp.ls = lapply(seq(snap.files), function(i){
  createSnap(file=snap.files[i], sample=sample.names[i])
})
names(ZF_with_subtype.sp.ls) = sample.names

barcode.ls = lapply(seq(snap.files), function(i){
  barcodes = read.csv(barcode.files[i], head=TRUE)
  # remove NO BAROCDE line
  barcodes = barcodes[2:nrow(barcodes),];
  barcodes$logUMI = log10(barcodes$passed_filters + 1);
  barcodes$promoter_ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1);
  barcodes
})

# Barcodes filtering folowing the same thresholds
cutoff.logUMI.low = c(3.5, 3.5, 3.5, 3.5)
cutoff.logUMI.high = c(5, 5, 5, 5)
barcode.ls = lapply(seq(snap.files), function(i){
  barcodes = barcode.ls[[i]]
  idx = which(
    barcodes$logUMI >= cutoff.logUMI.low[i] & 
    barcodes$logUMI <= cutoff.logUMI.high[i]
  )
  barcodes[idx,]
})
ZF_with_subtype.sp.ls = lapply(seq(snap.files), function(i){
  barcodes = barcode.ls[[i]]
  ZF_with_subtype.sp = ZF_with_subtype.sp.ls[[i]]
  barcode.shared = intersect(ZF_with_subtype.sp@barcode, barcodes$barcode)
  ZF_with_subtype.sp = ZF_with_subtype.sp[match(barcode.shared, ZF_with_subtype.sp@barcode),]
  barcodes = barcodes[match(barcode.shared, barcodes$barcode),]
  ZF_with_subtype.sp@metaData = barcodes
  ZF_with_subtype.sp
})
names(ZF_with_subtype.sp.ls) = sample.names
ZF_with_subtype.sp.ls

# combine two snap object
ZF_with_subtype.sp_raw = Reduce(snapRbind, ZF_with_subtype.sp.ls)
ZF_with_subtype.sp_raw@metaData["sample"] = ZF_with_subtype.sp_raw@sample

load("ZF.sp_overall_tsne.RData")
load("ZF.sp_with_subtype_raw.RData")
ZF_with_subtype.sp = ZF.sp
ZF_with_subtype.sp@file = ZF_with_subtype.sp_raw@file

# Only focusing on fibroblasts, endothelial cells and macrophages
ZF_with_subtype.sp = ZF_with_subtype.sp[ZF_with_subtype.sp@cluster %in% c(3, 8, 12),]

# loading scRNA-seq data for annotation of sub-clusters in scATAC-seq data
library(Seurat)
wt_all = readRDS("wt_nonCM_liger_as_seurat_38factors.RDS")
MC <- readRDS("MC_HM_updated.RDS")
wt_all <- SetIdent(object = wt_all, cells = rownames(MC@meta.data)[which(MC@meta.data$combined_cluster_HM_updated == "c1_init" & MC@meta.data$genotype == "wt")], value = 'MC1')
wt_all <- SetIdent(object = wt_all, cells = rownames(MC@meta.data)[which(MC@meta.data$combined_cluster_HM_updated == "c5_ap" & MC@meta.data$genotype == "wt")], value = 'MC2')
wt_all <- SetIdent(object = wt_all, cells = rownames(MC@meta.data)[which(MC@meta.data$combined_cluster_HM_updated == "c2_phag" & MC@meta.data$genotype == "wt")], value = 'MC3')
wt_all <- SetIdent(object = wt_all, cells = rownames(MC@meta.data)[which(MC@meta.data$combined_cluster_HM_updated == "c3_prol" & MC@meta.data$genotype == "wt")], value = 'MC4')
wt_all <- SetIdent(object = wt_all, cells = rownames(MC@meta.data)[which(MC@meta.data$combined_cluster_HM_updated == "c4_stressed" & MC@meta.data$genotype == "wt")], value = 'MC5')
FB <- readRDS("FB_liger_as_seurat.RDS")
wt_all <- SetIdent(object = wt_all, cells = rownames(FB@meta.data)[which(FB@meta.data$cell_subtype == "FB1" & FB@meta.data$genotype == "wt")], value = 'FB1')
wt_all <- SetIdent(object = wt_all, cells = rownames(FB@meta.data)[which(FB@meta.data$cell_subtype == "FB2" & FB@meta.data$genotype == "wt")], value = 'FB2')
wt_all <- SetIdent(object = wt_all, cells = rownames(FB@meta.data)[which(FB@meta.data$cell_subtype == "FB3" & FB@meta.data$genotype == "wt")], value = 'FB3')
wt_all <- SetIdent(object = wt_all, cells = rownames(FB@meta.data)[which(FB@meta.data$cell_subtype == "FB4" & FB@meta.data$genotype == "wt")], value = 'FB4')
EC <- readRDS("EC_liger_as_seurat.RDS")
wt_all <- SetIdent(object = wt_all, cells = rownames(EC@meta.data)[which(EC@meta.data$genotype_subtype == "wt_EC1" & EC@meta.data$genotype == "wt")], value = 'EC1')
wt_all <- SetIdent(object = wt_all, cells = rownames(EC@meta.data)[which(EC@meta.data$genotype_subtype == "wt_EC2" & EC@meta.data$genotype == "wt")], value = 'EC2')
wt_all <- SetIdent(object = wt_all, cells = rownames(EC@meta.data)[which(EC@meta.data$genotype_subtype == "wt_EC3" & EC@meta.data$genotype == "wt")], value = 'EC3')
wt_all <- SetIdent(object = wt_all, cells = rownames(EC@meta.data)[which(EC@meta.data$genotype_subtype == "wt_EC4" & EC@meta.data$genotype == "wt")], value = 'EC4')

wt_all = SubsetData(wt_all, cells = names(wt_all@active.ident[which(wt_all@active.ident %in% c("MC1", "MC2", "MC3", "MC4", "MC5", "FB1", "FB2", "FB3", "FB4", "EC1", "EC2", "EC3", "EC4"))]))
wt_all$tech = "rna"
variable.genes = rownames(wt_all@assays$RNA@scale.data)

# Subcluster annotation based on scRNA-seq
library(GenomicRanges)
genes.df = read.table("Danio_rerio.GRCz11.98.filtered.bed")
genes.gr = GRanges(genes.df[,1], IRanges(genes.df[,2], genes.df[,3]), name=genes.df[,4])
genes.sel.gr = genes.gr[which(genes.gr$name %in% variable.genes)]
ZF_with_subtype.sp = addBmatToSnap(ZF_with_subtype.sp)
ZF_with_subtype.sp = createGmatFromMat(obj=ZF_with_subtype.sp, input.mat="bmat", genes=genes.sel.gr, do.par=TRUE, num.cores=10)
ZF_with_subtype.atac <- snapToSeurat(obj= ZF_with_subtype.sp, eigs.dims=1:20, norm=TRUE, scale=TRUE)
transfer.anchors <- FindTransferAnchors(reference = wt_all, query = ZF_with_subtype.atac, features = variable.genes, 
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = wt_all@active.ident,
                                     weight.reduction = ZF_with_subtype.atac[["SnapATAC"]], dims = 1:15)
ZF_with_subtype.sp@metaData$predicted.id = celltype.predictions$predicted.id
ZF_with_subtype.sp@metaData$predict.max.score = apply(celltype.predictions[,-1], 1, max)
ZF_with_subtype.sp@cluster = as.factor(ZF_with_subtype.sp@metaData$predicted.id)
refdata <- GetAssayData(object = wt_all, assay = "RNA", slot = "data")
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = ZF_with_subtype.atac[["SnapATAC"]], dims = 1:20)
ZF_with_subtype.sp@gmat = t(imputation@data)

# Only retaining the cells correctly predicted for downstream analysis
id.correctly_predicted = which((ZF_with_subtype.sp@metaData$cluster == 8 & ZF_with_subtype.sp@metaData$predicted.id %in% c("EC1", "EC2", "EC3", "EC4")) | 
                               (ZF_with_subtype.sp@metaData$cluster == 3 & ZF_with_subtype.sp@metaData$predicted.id %in% c("FB1", "FB2", "FB3", ",FB4")) | 
                               (ZF_with_subtype.sp@metaData$cluster == 12 & ZF_with_subtype.sp@metaData$predicted.id %in% c("MC1", "MC2", "MC3", "MC4")))
ZF_with_subtype.sp = ZF_with_subtype.sp[id.correctly_predicted,]
id.correctly_predicted_eEC = which((ZF_with_subtype.sp@metaData$cluster == 8 & ZF_with_subtype.sp@metaData$predicted.id %in% c("EC1", "EC2", "EC3", "EC4")))
ZF_with_subtype.sp_eEC = ZF_with_subtype.sp[id.correctly_predicted_eEC,]
id.correctly_predicted_FB = which((ZF_with_subtype.sp@metaData$cluster == 3 & ZF_with_subtype.sp@metaData$predicted.id %in% c("FB1", "FB2", "FB3", ",FB4")))
ZF_with_subtype.sp_FB = ZF_with_subtype.sp[id.correctly_predicted_FB,]

# Only the prediction with high score were retained for downstream analysis
ZF_with_subtype.sp = ZF_with_subtype.sp[ZF_with_subtype.sp@metaData$predict.max.score > 0.3,]


# Running MACS for peak calling in each subcluster of fibroblasts, endothelial cells and macrophages
library(parallel)
source("runMACS_modified.R")
clusters.sel = names(table(ZF_with_subtype.sp@cluster))[which(table(ZF_with_subtype.sp@cluster) > 100)]
peaks.ls = mclapply(seq(clusters.sel), function(i){
    print(clusters.sel[i])
    runMACS.modified(
        obj=ZF_with_subtype.sp[which(ZF_with_subtype.sp@cluster == clusters.sel[i]),],
        output.prefix=paste0("ZF.", clusters.sel[i]),
        path.to.snaptools="/nas/longleaf/home/yyuchen/.pyenv/shims/snaptools",
        path.to.macs.shell="~/scATAC",
        gsize=1.5e9,
        buffer.size=500,
        num.cores=10,
        macs.options="--nomodel --shift 100 --ext 200 --qval 1e-2 -B --SPMR --call-summits",
        tmp.folder=tempdir()
    )
}, mc.cores = 2)

# assuming all .narrowPeak files in the current folder are generated from the clusters
library(GenomicRanges)
peaks.names = system("ls /proj/liulab/users/yanhan/scATAC_yanhan/results/narrowPeak_with_subtype/ | grep narrowPeak", intern=TRUE)
peak.gr.ls = lapply(peaks.names, function(x){
  peak.df = read.table(paste0("results/narrowPeak_with_subtype/", x))
  GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
})
peak.gr = reduce(Reduce(c, peak.gr.ls))


# Create a cell-by-peak matrix
peaks.df = as.data.frame(peak.gr)[,1:3];
write.table(peaks.df,file = "ZF.peaks.combined_with_subtype.bed",append=FALSE,
            quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".",
            row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

# Create cell-by-peak matrix and add to the snap file in linux
$ module add bedtools
$ snaptools snap-add-pmat --snap-file day0_scATAC/day0_with_subtype.snap --peak-file ZF.peaks.combined_with_subtype.bed
$ snaptools snap-add-pmat --snap-file day2_scATAC/day2_with_subtype.snap --peak-file ZF.peaks.combined_with_subtype.bed
$ snaptools snap-add-pmat --snap-file day7_scATAC/day7_with_subtype.snap --peak-file ZF.peaks.combined_with_subtype.bed
$ snaptools snap-add-pmat --snap-file day14_scATAC/day14_with_subtype.snap --peak-file ZF.peaks.combined_with_subtype.bed

# Add the cell-by-peak matrix to the existing snap object in R
ZF_with_subtype.sp = addPmatToSnap(ZF_with_subtype.sp)
ZF_with_subtype.sp = makeBinary(ZF_with_subtype.sp, "pmat")

library(umap)
ZF_with_subtype.sp = runViz(obj=ZF_with_subtype.sp, tmp.folder=tempdir(), dims=2, eigs.dims=1:15, method="umap", seed.use=10)

# Motif variability analysis using chromVAR 
library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(BSgenome.Drerio.UCSC.danRer11)
library(GenomicRanges)
library(TFBSTools)
library(ggplot2)
source("runChromVAR_modified.R")
ZF_with_subtype.sp@mmat = runChromVAR_modified(
  obj=ZF_with_subtype.sp,
  input.mat="pmat",
  genome=BSgenome.Drerio.UCSC.danRer11,
  min.count=10
)

#  Identify differentially accessible peaks in, for example, FB3
source("/proj/liulab/users/ycyang/scATAC/mouse_lacZ_MGT/lacZ_MGT_aggr/findDAR_new.R")
DARs = findDAR_new(obj=ZF_with_subtype.sp, input.mat="pmat", cluster.pos="FB3", cluster.neg=c("FB1", "FB2"), 
  cluster.neg.method="knn", test.method="exactTest", bcv=0.1, seed.use=10)
DARs$FDR = p.adjust(DARs$PValue, method="BH")
idy = which(DARs$FDR < 5e-2 & DARs$logFC > 0)
if ((x=length(idy)) < 2000L){
   PValues = DARs$PValue;
   PValues[DARs$logFC < 0] = 1;
   idy = order(PValues, decreasing=FALSE)[1:2000];
   rm(PValues); # free memory
}

# Motif analysis identifies master regulators using HOMER
motifs = runHomer(ZF_with_subtype.sp[,idy,"pmat"], mat = "pmat", path.to.homer = "~/bin/findMotifsGenome.pl", result.dir = "homer_with_subtype/FB/FB3", num.cores=5,
  genome = 'danRer11', motif.length = 10, scan.size = 300, optimize.count = 2, background = 'automatic', local.background = FALSE, only.known = TRUE,
  only.denovo = FALSE, fdr.num = 5, cache = 100, overwrite = TRUE, keep.minimal = FALSE)



# UMAP plot for feature genes
library(GenomicRanges)
genes = read.table("/proj/liulab/users/yfxie/Yuchen/scATAC/scATAC_zebrafish/Danio_rerio.GRCz11.98.filtered.bed")
genes.gr = GRanges(genes[,1], IRanges(genes[,2], genes[,3]), name=genes[,4])
marker.genes = c("fosl2", "bach2a", "bach2b", "batf3", "fn1a", "fn1b", "nppc", "tnfa", "nr4a1")
genes.sel.gr <- genes.gr[which(genes.gr$name %in% marker.genes)]
ZF_with_subtype.sp = createGmatFromMat(obj=ZF_with_subtype.sp, input.mat="bmat", genes=genes.sel.gr, do.par=TRUE, num.cores=10)
ZF_with_subtype.sp = scaleCountMatrix(obj=ZF_with_subtype.sp, cov=ZF_with_subtype.sp@metaData$passed_filters + 1, mat="gmat", method = "RPM")
ZF_with_subtype.sp = runMagic(obj=ZF_with_subtype.sp, input.mat="gmat", step.size=3)

library(ggplot2)
cols = c("#90BFF9", "#cfe3fc", "#e8f0fa", "white", "orange", "red", "firebrick")
df <- data.frame(xval=c(1:100), yval=c(1:100))
test.plot = ggplot(df, aes(x=xval, y=yval, colour=yval)) + scale_color_gradientn(colors = cols, guide = "colorbar")
col.used = ggplot_build(test.plot)$data[[1]][,1]
source("plotFeatureSingle_modified.R")
par(mfrow = c(3, 3));
for(i in 1:9){
  plotFeatureSingle_modified(obj=ZF_with_subtype.sp, feature.value=ZF_with_subtype.sp@gmat[,i], method="umap",
    main=colnames(ZF_with_subtype.sp@gmat)[i], point.size=0.1, point.shape=19, down.sample=10000, quantiles=c(0.01, 0.99))
}


# Predict gene-enhancer pairs
library(GenomicRanges)
### fosl2
TSS.loci = GRanges("17", IRanges(41302605, 41302606))
pairs = predictGenePeakPair(ZF_with_subtype.sp, input.mat = "pmat", gene.name = "fosl2", gene.loci=resize(TSS.loci, width = 500000, fix="center"), do.par = F)
pairs.df = as.data.frame(pairs)
pairs.df_all = data.frame(chr1 = paste0("chr", pairs.df[,"seqnames"]), start1=pairs.df[,"start"], end1=pairs.df[,"end"], chr = paste0("chr17:41302605-41302606,", pairs.df[,"logPval"]), number = rownames(pairs.df), plus = "-")
write.table(pairs.df_all, file = "fosl2_all_interaction.txt", sep = "\t", quote = F, row.names = F, col.names = F)
pairs.df_sig = pairs.df[pairs.df$Pval < 0.05,]
pairs.df_sig = data.frame(chr1 = paste0("chr", pairs.df_sig[,"seqnames"]), start1=pairs.df_sig[,"start"], end1=pairs.df_sig[,"end"], chr = paste0("chr17:41302605-41302606,", pairs.df_sig[,"logPval"]), number = rownames(pairs.df_sig), plus = "-")
write.table(pairs.df_sig, file = "fosl2_significant_interaction.txt", sep = "\t", quote = F, row.names = F, col.names = F)






