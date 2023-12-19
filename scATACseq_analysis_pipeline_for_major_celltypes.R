library(SnapATAC)
snap.files = c("day0_scATAC/day0.snap", "day2_scATAC/day2.snap",
               "day7_scATAC/day7.snap", "day14_scATAC/day14.snap")
sample.names = c("day0", "day2", "day7", "day14")
barcode.files = c("day0_scATAC/singlecell.csv", "day2_scATAC/singlecell.csv",
                  "day7_scATAC/singlecell.csv", "day14_scATAC/singlecell.csv")
ZF.sp.ls = lapply(seq(snap.files), function(i){
  createSnap(file=snap.files[i], sample=sample.names[i])
})
names(ZF.sp.ls) = sample.names

barcode.ls = lapply(seq(snap.files), function(i){
  barcodes = read.csv(barcode.files[i], head=TRUE)
  # remove NO BAROCDE line
  barcodes = barcodes[2:nrow(barcodes),];
  barcodes$logUMI = log10(barcodes$passed_filters + 1);
  barcodes$promoter_ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1);
  barcodes
})

library(ggplot2)
plot_list = list()
for (i in seq(snap.files)) {
  p = ggplot(barcode.ls[[i]], aes(x=logUMI, y=promoter_ratio)) + 
      geom_point(size=0.3, col="grey") + theme_classic() +
      ggtitle(sample.names[[i]]) + ylim(0, 1) + xlim(0, 6) + 
      labs(x = "log10(UMI)", y="promoter ratio")
  plot_list[[i]] = p
}

# We identify usable barcodes using [3.5-5] for log10(UMI) as cutoff.
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
ZF.sp.ls = lapply(seq(snap.files), function(i){
  barcodes = barcode.ls[[i]]
  ZF.sp = ZF.sp.ls[[i]]
  barcode.shared = intersect(ZF.sp@barcode, barcodes$barcode)
  ZF.sp = ZF.sp[match(barcode.shared, ZF.sp@barcode),]
  barcodes = barcodes[match(barcode.shared, barcodes$barcode),]
  ZF.sp@metaData = barcodes
  ZF.sp
})
names(ZF.sp.ls) = sample.names

# combine two snap object
ZF.sp = Reduce(snapRbind, ZF.sp.ls)
ZF.sp@metaData["sample"] = ZF.sp@sample

# Add cell-by-bin matrix
ZF.sp = addBmatToSnap(ZF.sp, bin.size=5000, num.cores=1)

# Matrix binarization. we next remove 0.1% items of the highest coverage in the count matrix and then convert the remaining non-zero items to 1.
ZF.sp = makeBinary(ZF.sp, mat="bmat")

# The top 5% bins that overlap with invariant features such as promoters of the house keeping genes were removed.
bin.cov = log10(Matrix::colSums(ZF.sp@bmat)+1)
hist(bin.cov[bin.cov > 0],  xlab="log10(bin cov)",  main="log10(Bin Cov)",  col="lightblue",  xlim=c(0, 5))
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95)
idy = which(bin.cov <= bin.cutoff & bin.cov > 0)
#idy = which(bin.cov > 0)
ZF.sp = ZF.sp[, idy, mat="bmat"]

# Dimentionality reduction.
ZF.sp = runDiffusionMaps(obj=ZF.sp, input.mat="bmat", num.eigs=50)

# PCA plot
plotDimReductPW(obj=ZF.sp, eigs.dims=1:50, point.size=0.3, point.color="grey", point.shape=19,
    point.alpha=0.6, down.sample=5000, pdf.file.name=NULL, pdf.height=7, pdf.width=7)

# Graph-based clustering
ZF.sp = runKNN(obj=ZF.sp, eigs.dims=1:15, k=15)
ZF.sp=runCluster(obj=ZF.sp, tmp.folder=tempdir(), louvain.lib="R-igraph", seed.use=10)
ZF.sp@metaData$cluster = ZF.sp@cluster

# Visualization
library(umap)
ZF.sp = runViz(obj=ZF.sp, tmp.folder=tempdir(), dims=2, eigs.dims=1:15, method="umap", seed.use=10)
plotViz(obj=ZF.sp, method="umap", main="ZF clustering", point.color=ZF.sp@cluster,
  point.size=0.2, point.shape=19, point.alpha=0.8,
  text.add=F, text.size=1, text.color="black", down.sample=10000, legend.add=TRUE
)

# Expression of marker genes
library(GenomicRanges)
genes = read.table("Danio_rerio.GRCz11.98.filtered.bed")
genes.gr = GRanges(genes[,1], IRanges(genes[,2], genes[,3]), name=genes[,4])
marker.genes = c("nppa", "tnnt2a", "itga2b", "slc4a1a", "cahz", "tcf21", "fn1b", "kdrl", "fli1a", "rca2.2", "cldn5b", "spock3", "f8", "angptl7", "aif1l", "tnfaip6", "repo1", "mgp", "mfap4", "cd74a", "lyz", "mpx", "lect2l", "npsn", "pnp5a", "cxcr4a", "irf4b", "rgs5b", "ccl25b", "tagln")
genes.sel.gr <- genes.gr[which(genes.gr$name %in% marker.genes)]
ZF.sp = createGmatFromMat(obj=ZF.sp, input.mat="bmat", genes=genes.sel.gr, do.par=TRUE, num.cores=10)
ZF.sp = scaleCountMatrix(obj=ZF.sp, cov=ZF.sp@metaData$passed_filters + 1, mat="gmat", method = "RPM")
ZF.sp = runMagic(obj=ZF.sp, input.mat="gmat", step.size=3)

for(i in 1:29){
  plotFeatureSingle(obj=ZF.sp, feature.value=ZF.sp@gmat[,i], method="umap",
      main=colnames(ZF.sp@gmat)[i], point.size=0.1, point.shape=19, down.sample=10000, quantiles=c(0.01, 0.99))
}


### cluster 11 is cardiomyocyte and cluster 1, 3, 4, 5, 6, 12, 14, 15 and 16 are the clusters we analysis in the previous study.
ZF.sp = ZF.sp[ZF.sp@cluster %in% c(1, 3:6, 12, 14, 15, 16),]
celltype.label = c()
for(i in ZF.sp@cluster){
   if (i == 1){
      celltype.label = c(celltype.label, "eEC")
   } else if (i == 3){
      celltype.label = c(celltype.label, "Neutro")
   } else if (i == 4){
      celltype.label = c(celltype.label, "MC")
   } else if (i == 5){
      celltype.label = c(celltype.label, "FB")
   } else if (i == 6){
      celltype.label = c(celltype.label, "Eryth")
   } else if (i == 12){
      celltype.label = c(celltype.label, "Mes")
   } else if (i == 14){
      celltype.label = c(celltype.label, "cEC")
   } else if (i == 15){
      celltype.label = c(celltype.label, "Throm")
   } else if (i == 16){
      celltype.label = c(celltype.label, "T_NK_B")
   }
}
ZF.sp@cluster = factor(celltype.label)


# Call peaks in each cell type
library(parallel)
source("runMACS_modified.R")
clusters.sel = names(table(ZF.sp@cluster))[which(table(ZF.sp@cluster) > 100)]
peaks.ls = mclapply(seq(clusters.sel), function(i){
    print(clusters.sel[i])
    runMACS.modified(
        obj=ZF.sp[which(ZF.sp@cluster == clusters.sel[i]),],
        output.prefix=paste0("ZF.", clusters.sel[i]),
        path.to.snaptools="/nas/longleaf/home/yyuchen/.pyenv/shims/snaptools",
        path.to.macs.shell="/proj/liulab/users/ycyang/scATAC",
        gsize=1.5e9,
        buffer.size=500,
        num.cores=10,
        macs.options="--nomodel --shift 100 --ext 200 --qval 1e-2 -B --SPMR --call-summits",
        tmp.folder=tempdir()
    )
}, mc.cores = 2)

library(GenomicRanges)
peaks.names = system("ls results_overall/narrowPeak_overall/ | grep narrowPeak", intern=TRUE)
peak.gr.ls = lapply(peaks.names, function(x){
  peak.df = read.table(paste0("results_overall/narrowPeak_overall/", x))
  GRanges(peak.df[,1], IRanges(peak.df[,2], peak.df[,3]))
})
peak.gr = reduce(Reduce(c, peak.gr.ls))
peak.gr

# Create a cell-by-peak matrix
peaks.df = as.data.frame(peak.gr)[,1:3];
write.table(peaks.df,file = "ZF.peaks.combined_overall.bed",append=FALSE,
              quote= FALSE,sep="\t", eol = "\n", na = "NA", dec = ".",
              row.names = FALSE, col.names = FALSE, qmethod = c("escape", "double"),
              fileEncoding = "")

# Create cell-by-peak matrix and add to the snap file in linux
$ module add bedtools
$ snaptools snap-add-pmat --snap-file day0_scATAC/day0.snap --peak-file ZF.peaks.combined_overall.bed
$ snaptools snap-add-pmat --snap-file day2_scATAC/day2.snap --peak-file ZF.peaks.combined_overall.bed
$ snaptools snap-add-pmat --snap-file day7_scATAC/day7.snap --peak-file ZF.peaks.combined_overall.bed
$ snaptools snap-add-pmat --snap-file day14_scATAC/day14.snap --peak-file ZF.peaks.combined_overall.bed

# Add the cell-by-peak matrix to the existing snap object in R
ZF.sp = addPmatToSnap(ZF.sp)
ZF.sp = makeBinary(ZF.sp, "pmat")

# Motif enrichment analysis using chromVAR
library(chromVAR)
library(motifmatchr)
library(SummarizedExperiment)
library(BSgenome.Drerio.UCSC.danRer11)
library(GenomicRanges)
library(TFBSTools)
library(ggplot2)
source("runChromVAR_modified.R")
ZF.sp@mmat = runChromVAR_modified(
  obj=ZF.sp,
  input.mat="pmat",
  genome=BSgenome.Drerio.UCSC.danRer11,
  min.count=10
)


#  Identify differentially accessible peaks
source("findDAR_new.R")
idy.ls = lapply(levels(ZF.sp@cluster), function(cluster_i){
  DARs = findDAR_new(
    obj=ZF.sp,
    input.mat="pmat",
    cluster.pos=cluster_i,
    cluster.neg=levels(ZF.sp@cluster)[which(levels(ZF.sp@cluster) != cluster_i)],
    cluster.neg.method="knn",
    test.method="exactTest",
    bcv=0.1,
    seed.use=10
  )
  DARs$FDR = p.adjust(DARs$PValue, method="BH")
  idy = which(DARs$FDR < 5e-2 & DARs$logFC > 0)
  if ((x=length(idy)) < 2000L){
     PValues = DARs$PValue;
     PValues[DARs$logFC < 0] = 1;
     idy = order(PValues, decreasing=FALSE)[1:2000];
     rm(PValues); # free memory
  }
  idy
})
names(idy.ls) = levels(ZF.sp@cluster)
covs = Matrix::rowSums(ZF.sp@pmat)
library(rlist)
summary(idy.ls)


# Motif analysis identifies master regulators
motif.ls = lapply(levels(ZF.sp@cluster), function(cluster_i){
  motifs = runHomer(
    ZF.sp[,idy.ls[[cluster_i]],"pmat"],
    mat = "pmat",
    path.to.homer = "~/bin/findMotifsGenome.pl",
    result.dir = paste0("homer_overall/", cluster_i),
    num.cores=5,
    genome = 'danRer11',
    motif.length = 10,
    scan.size = 300,
    optimize.count = 2,
    background = 'automatic',
    local.background = FALSE,
    only.known = TRUE,
    only.denovo = FALSE,
    fdr.num = 5,
    cache = 100,
    overwrite = TRUE,
    keep.minimal = FALSE
  )
  motifs
})

