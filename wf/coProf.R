library("ape")
library("ArchR")
library("biovizBase")
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome.Mmusculus.UCSC.mm10")
library("cowplot")
library("dplyr")
library("EnsDb.Hsapiens.v86")
library("EnsDb.Mmusculus.v79")
library("GenomicRanges")
library("getopt")
library("ggplot2")
library("patchwork")
library("pdftools")
library("plyr")
library("readr")
library("Seurat")
library("Signac")
library("qdap")

source("/root/addgeneintegrationmatrix.R")


# Compute before filter stats
qc_stats_df <- function(seurat_object, row.name) {
  df <- data.frame(
    Num_Tixels = round(ncol(seurat_object)),
    UMI_Avg = round(mean(seurat_object@meta.data$nCount_Spatial)),
    Genes_Avg = round(mean(seurat_object@meta.data$nFeature_Spatial)),
    percent_mito_Avg = round(mean(seurat_object@meta.data$percent_mito)),
    percent_hb_Avg = round(min(seurat_object@meta.data$percent_hb)),
    row.names = row.name
  ) %>% t()
  return(df)
}

dir.create(file.path("reports"), recursive = TRUE)
setwd("/root/reports")

args <- commandArgs(trailingOnly = TRUE)

project_name <- args[1]
species <- args[2]
point_size <- as.integer(args[3])
resolution <- as.numeric(args[4])

tile_size <- as.integer(args[5])
min_tss <- as.numeric(args[6])
min_frags <- as.integer(args[7])
lsi_iterations <- as.integer(args[8])
lsi_resolution <- as.numeric(args[9])
lsi_varfeatures <- as.integer(args[10])

runs <- strsplit(args[11:length(args)], ",")
runs

obj.lst_off <- list()
for (run in runs) {
  data <- Read10X(data.dir = run[2])
  colnames(data) <- paste0(
    run[1], "#", colnames(Read10X(data.dir = run[2])), "-1"
  )
  obj.lst_off[run[1]] <- CreateSeuratObject(counts = data, assay = "Spatial")
  obj.lst_off[run[1]][[1]]$Sample <- rep(
    run[1], ncol(obj.lst_off[run[1]][[1]])
  )
  image <- Read10X_Image(
    image.dir = run[4],
    filter.matrix = FALSE
  )
  rownames(image@coordinates) <- paste0(
    run[1], "#", rownames(image@coordinates), "-1"
  )
  image <- image[Cells(x = obj.lst_off[run[1]][[1]])]
  DefaultAssay(object = image) <- "Spatial"

  obj.lst_off[run[1]][[1]][[run[1]]] <- image
  on_tissue_cells <- rownames(obj.lst_off[run[1]][[1]]@images[[1]]@coordinates)
  obj.lst_off[run[1]] <- obj.lst_off[run[1]][[1]][, on_tissue_cells]
  obj.lst_off[run[1]][[1]]@meta.data$row <- obj.lst_off[run[1]][[1]]@images[[1]]@coordinates$row
  obj.lst_off[run[1]][[1]]@meta.data$col <- obj.lst_off[run[1]][[1]]@images[[1]]@coordinates$col
}

obj.lst_off
bar_plot_r <- list()

for (i in seq_along(obj.lst_off)) {

  nbc <- sqrt(dim(obj.lst_off[[i]])[2])
  if (nbc == 96) {
    ncol <- 1
  } else {
    ncol <- 2
  }

  # Compute umi by row/col before filtering
  meta.data_plots <- obj.lst_off[[i]]@meta.data

  # since number of row/col start with 0 we change it to start with 1
  meta.data_plots[, c("col", "row")] <- meta.data_plots[, c("col", "row")] + 1
  meta.data_plots$col <- as.character(meta.data_plots$col)
  meta.data_plots$row <- as.character(meta.data_plots$row)
  meta.data_plots$col <- factor(meta.data_plots$col, levels = 0:nbc)
  meta.data_plots$row <- factor(meta.data_plots$row, levels = 0:nbc)

  umi_row <- ggplot(
    data = meta.data_plots, aes(x = col, y = nCount_Spatial)
  ) +
    geom_bar(stat = "identity", width = 0.75) +
    ylab("#UMIs / row") +
    xlab(paste("Row (", project_name, ")", sep = "")) +
    theme(text = element_text(size = 10))

  # umi_row
  umi_column <- ggplot(
    data = meta.data_plots,
    aes(x = row, y = nCount_Spatial)
  ) +
    geom_bar(stat = "identity", width = 0.75) +
    ylab("#UMIs / column") +
    xlab(paste("Column (", project_name, ")", sep = "")) +
    theme(text = element_text(size = 10))

  options(repr.plot.width = 10, repr.plot.height = 5)
  bar_plot_r[[i]] <- wrap_plots(umi_row, umi_column, ncol = ncol)
}

names(bar_plot_r) <- names(obj.lst_off)

obj.lst_on <- list()
for (run in runs) {
  data <- Read10X(data.dir = run[2])
  colnames(data) <- paste0(
    run[1], "#", colnames(Read10X(data.dir = run[2])), "-1"
  )
  obj.lst_on[run[1]] <- CreateSeuratObject(counts = data, assay = "Spatial")
  obj.lst_on[run[1]][[1]]$Sample <- rep(run[1], ncol(obj.lst_on[run[1]][[1]]))
  image <- Read10X_Image(
    image.dir = run[4],
    filter.matrix = TRUE
  )
  rownames(image@coordinates) <- paste0(
    run[1], "#", rownames(image@coordinates), "-1"
  )
  image <- image[Cells(x = obj.lst_on[run[1]][[1]])]
  DefaultAssay(object = image) <- "Spatial"

  obj.lst_on[run[1]][[1]][[run[1]]] <- image


  on_tissue_cells <- rownames(obj.lst_on[run[1]][[1]]@images[[1]]@coordinates)
  obj.lst_on[run[1]] <- obj.lst_on[run[1]][[1]][, on_tissue_cells]
  obj.lst_on[run[1]][[1]]@meta.data$row <- obj.lst_on[run[1]][[1]]@images[[1]]@coordinates$row
  obj.lst_on[run[1]][[1]]@meta.data$col <- obj.lst_on[run[1]][[1]]@images[[1]]@coordinates$col
}

saveRDS(obj.lst_on, "RNA_lst.rds")
cat(paste0("RDS saved in ", getwd(), "/", "RNA_lst.rds"))

ctm <- c()
meta.data <- c()
for (i in names(obj.lst_on)){
  tmp <- obj.lst_on[[i]]@assays$Spatial@layers$counts
  colnames(tmp) <- colnames(obj.lst_on[[i]])
  rownames(tmp) <- rownames(obj.lst_on[[i]])
  tmp_met <- obj.lst_on[[i]]@meta.data
  ctm <- cbind(ctm, tmp)
  meta.data <- rbind(meta.data, tmp_met)
}

# Use Reduce to merge the data frames sequentially
brain <- CreateSeuratObject(
  counts = ctm, meta.data = meta.data, assay = "Spatial"
)
for (i in names(obj.lst_on)){
  brain@images[[i]] <- obj.lst_on[[i]]@images[[1]]
}
brain <- PercentageFeatureSet(brain, "^mt-", col.name = "percent_mito")
brain <- PercentageFeatureSet(brain, "^Hb.*-", col.name = "percent_hb")
brain

on_tiss_before_filter <- qc_stats_df(brain, row.name = "Before Filtering")

# filtering
min_nFeature_non_zero <- min(
  brain@meta.data$nFeature_Spatial[brain@meta.data$nFeature_Spatial > 0]
)
max_percent_mito <- max(brain@meta.data$percent_mito)
max_percent_hb <- max(brain@meta.data$percent_hb)

res <- try(
  brain_filt <- brain[
    ,
    brain$nFeature_Spatial >= min_nFeature_non_zero &
      brain$percent_mito <= max_percent_mito &
      brain$percent_hb <= max_percent_hb
  ]
)

if (inherits(res, "try-error")) {
  #error handling code, maybe just skip this iteration using
  brain_filt <- brain
}

# Filter Bl1
brain_filt <- brain_filt[!grepl("Bc1", rownames(brain_filt)), ]

# Filter Mitocondrial
brain_filt <- brain_filt[!grepl("^mt-", rownames(brain_filt)), ]

# # Filter Hemoglobin gene (optional if that is a problem on your data)
brain_filt <- brain_filt[!grepl("^Hb.*-", rownames(brain_filt)), ]

brain_filt

on_tiss_after_filter <- qc_stats_df(brain_filt, row.name = "After Filtering")

on_tiss_stats_r <- cbind(on_tiss_before_filter, on_tiss_after_filter)

qc_rna <- VlnPlot(
  brain_filt,
  features = c(
    "nCount_Spatial", "nFeature_Spatial", "percent_mito", "percent_hb"
  ),
  pt.size = 0.1,
  ncol = 4
) +
  NoLegend()

cat("RNA-seq data analysis starting ...")

rna_processing <- function(resolution) {

  res <- try( # without this filtering we can't use SCTransform
    RNA <- SCTransform(
      subset(brain_filt, subset = nCount_Spatial > 0),
      assay = "Spatial",
      verbose = TRUE,
      method = "poisson"
    )
  )

  if (inherits(res, "try-error")) {
    #error handling code, maybe just skip this iteration using
    print(
      "SCT can not be used because there is not enough cells, so we use the \
      object with out nCount_Spatial > 0"
    )

    RNA <- brain_filt

    print("Running Clustering and UMAP Perform standard analysis not SCT! ")

    npc <- min(nrow(RNA), ncol(RNA))

    if (npc > 50) {
      npcs <- 50
      n.neighbors <- 30
    } else {
      npcs <- npc - 1
      n.neighbors <- npcs
    }

    print("number of pca is")
    print(npcs)
    print("number of neighbors is")
    print(n.neighbors)

    RNA <- NormalizeData(
      RNA, normalization.method = "LogNormalize", scale.factor = 10000
    )
    RNA <- FindVariableFeatures(RNA, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(RNA)
    RNA <- ScaleData(RNA, features = all.genes)

    RNA <- RunPCA(RNA, npcs = npcs, verbose = TRUE)

    RNA <- FindNeighbors(RNA, reduction = "pca", dims = 1:5)
    RNA <- FindClusters(
      RNA, resolution = resolution, verbose = TRUE, random.seed = 101
    )
    RNA <- RunUMAP(
      RNA,
      reduction = "pca",
      n.neighbors = n.neighbors,
      dims = 1:5,
      seed.use = 101
    )

  } else {

    print("Running Clustering and UMAP Perform SCT! ")

    RNA <- RunPCA(RNA, assay = "SCT", verbose = FALSE)
    RNA <- FindNeighbors(RNA, reduction = "pca", dims = 1:30)
    RNA <- FindClusters(RNA, verbose = FALSE, resolution = resolution)
    RNA <- RunUMAP(RNA, reduction = "pca", dims = 1:30)
  }
}

RNA <- rna_processing(resolution)

print("check number of clusters")

if (length(unique(RNA$seurat_clusters)) < 2) {
  print(
    "since number of clusters is 1, i will change the resolution to 1 to get \
    more clusters"
  )
  seuratObj <- rna_processing(1)
} else {
  print("will not change the resolution")
  seuratObj <- RNA
}

seuratObj$seurat_clusters <- paste0("R", seuratObj$seurat_clusters)
Idents(seuratObj) <- seuratObj$seurat_clusters

saveRDS(seuratObj, "RNA_all.rds")
cat(paste0("RDS saved in ", getwd(), "/", "RNA_all.rds"))

req_DF <- as.data.frame(seuratObj@meta.data)

req_table <- melt(table(req_DF$seurat_clusters, req_DF$Sample))

colnames(req_table) <- c("Cluster", "Sample", "%cells in knn clusters")

req_table$Cluster <- factor(
  req_table$Cluster,
  levels = unique(
    req_table[order(as.numeric(gsub("C", "", req_table$Cluster))), ]$Cluster
  )
)

pal <- paletteDiscrete(values = seuratObj$Sample)

sample_bar_plots_r <- ggplot(
  req_table,
  aes(fill = Sample, y = "%cells in knn clusters", x = Cluster)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_classic() +
    scale_fill_manual(values = pal) +
    theme(text = element_text(size = 30)) +
    theme(axis.title.x = element_blank())

n_clusters <- length(unique(seuratObj$seurat_clusters))
if (n_clusters >= 5) {
  ncol <- 5
} else {
  ncol <- n_clusters
}

color <- colorRampPalette(
  c(
    "blue",
    "red",
    "dark green",
    "orange",
    "purple",
    "brown",
    "cyan",
    "yellow"
  )
)(n_clusters)
cols1 <- ArchRPalettes$stallion2[seq_len(n_clusters)]
names(cols1) <- sort(unique(seuratObj$seurat_clusters))
n_samples <- length(unique(seuratObj$Sample))

cols2 <- pal
names(cols2) <- unique(seuratObj$Sample)

cols <- c(cols1, cols2)

dimrd3_r <- SpatialDimPlot(
  seuratObj,
  group.by = "seurat_clusters",
  pt.size.factor = point_size,
  cols = cols
)

dimrd1_r <- DimPlot(
  seuratObj,
  reduction = "umap",
  group.by = c("seurat_clusters", "Sample"),
  cols = cols
)

dimrd4_r <- list()
Idents(seuratObj) <- seuratObj$seurat_clusters
cells.highlight <- CellsByIdentities(
  object = seuratObj
)[order(names(CellsByIdentities(object = seuratObj)))]

for (i in names(seuratObj@images)){

  dimrd4_r[[i]] <- SpatialDimPlot(
    seuratObj,
    cells.highlight = cells.highlight,
    facet.highlight = TRUE,
    ncol = ncol,
    images = i
  )
}

cat("FindAllMarkers starting ...")

n_clusters <- length(unique(seuratObj$seurat_clusters))
n_samples <- length(unique(seuratObj$Sample))
samples <- unique(seuratObj$Sample)

# Identification of Spatially Variable Features
de_markers <- try(
  FindAllMarkers(
    seuratObj,
    only.pos = TRUE,
    logfc.threshold = 0.01,
    min.pct = 0.01
  )
)

if (nrow(de_markers) < n_clusters) {
  print(
    "Not Enough features pass min.pct=0.01 and threshold=0.1; returning empty \
    data.frame"
  )
  top_spatial_plots_r <- list()
} else {
  print("there we gooooooo!")
  top1 <- de_markers %>% group_by(cluster) %>% top_n(n = 1)
  brain_ss <- list()
  for (i in samples){

    brain_ss[[i]] <- seuratObj[, which(seuratObj$Sample == i)]
    brain_ss[[i]]@images <- seuratObj@images[i]
  }

  topGenes <- top1$gene
  print(topGenes)
  if (length(topGenes) >= 5) {
    ncol <- 5
  } else {
    ncol <- length(topGenes)
  }

  top_spatial_plots_r <- list()
  top_list <- list()

  for (i in seq_along(samples)) {
    for (j in topGenes){
      top_list[[j]] <- SpatialFeaturePlot(
        brain_ss[[i]],
        features = j,
        ncol = 1,
        alpha = c(.1, 1),
        pt.size.factor = point_size) +
          theme(
            legend.position = "top",
            legend.text = element_text(size = 5),
            legend.title = element_text(size = 20)
          )

      top_list[[j]]$layers[[1]]$aes_params <- c(
        top_list[[j]]$layers[[1]]$aes_params, shape = 22
      )

      top_spatial_plots_r[[i]] <- wrap_plots(top_list, ncol = ncol)
    }
  }
}

names(top_spatial_plots_r) <- samples

if (rlang::is_empty(top_spatial_plots_r) == TRUE) {
  cat("Not Enough features pass min.pct=0.01 and threshold=0.1")
} else {
  top_spatial_plots_r
}

# Identification of Spatially Variable Features
seuratObj <- FindSpatiallyVariableFeatures(
  seuratObj,
  assay = "SCT",
  features = VariableFeatures(seuratObj)[1:1000],
  selection.method = "moransi"
)

SpatiallyVariableFeatures_workaround <- function(
  object, assay = "SCT", selection.method = "moransi"
) {
  #' This is work around function to replace
  #' SeuratObject::SpatiallyVariableFeatures function. Returns ranked list of
  #' Spatially Variable Features

  # Check if object is a Seurat object
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }

  # Check if assay is a valid assay
  if (!assay %in% names(object@assays)) {
    stop("assay must be a valid assay")
  }

  # Extract meta.features from the specified object and assay
  data <- object@assays[[assay]]@meta.features

  # Select columns starting with the provided col_prefix
  moransi_cols <- grep(
    paste0("^", selection.method), colnames(data), value = TRUE
  )

  # Filter rows where "moransi.spatially.variable" is TRUE
  filtered_data <- data[
    data[[paste0(selection.method, ".spatially.variable")]],
    moransi_cols
  ]

  # Sort filtered data by "moransi.spatially.variable.rank" column in
  # ascending order
  sorted_data <- filtered_data[order(
    filtered_data[[paste0(selection.method, ".spatially.variable.rank")]]
  ), ]

  # Return row names of the sorted data frame
  rownames(sorted_data)
}

top.features <- head(
  SpatiallyVariableFeatures_workaround(
    seuratObj,
    selection.method = "moransi"
  ),
  12
)

top.features

if (length(top.features) >= 5) {
  ncol <- 5
} else {
  ncol <- length(top.features)
}

top_spatial_plots_r2 <- list()
top_list <- list()

for (i in seq_along(samples)) {
  for (j in top.features) {
    top_list[[j]] <- SpatialFeaturePlot(
      brain_ss[[i]],
      features = j,
      ncol = 1,
      alpha = c(.1, 1),
      pt.size.factor = point_size
    ) +
      theme(
        legend.position = "top",
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 20)
      )

    top_list[[j]]$layers[[1]]$aes_params <- c(
      top_list[[j]]$layers[[1]]$aes_params, shape = 22
    )

    top_spatial_plots_r2[[i]] <- wrap_plots(top_list, ncol = ncol)
  }
}

names(top_spatial_plots_r2) <- samples
if (rlang::is_empty(top_spatial_plots_r2) == TRUE) {
  cat("Not Enough features by moransi method")
} else {
  top_spatial_plots_r2
}

##############################################################################
cat("Creation and Processing of ATAC-seq")

genome <- args[2]

if (genome == "hg38") {
  addArchRGenome(genome)
  geneAnnotation <- getGeneAnnotation()
  genomeAnnotation <- getGenomeAnnotation()
  genomeSize <- 3.3e+09
} else if (genome == "mm10") {
  addArchRGenome(genome)
  geneAnnotation <- getGeneAnnotation()
  genomeAnnotation <- getGenomeAnnotation()
  genomeSize <- 3.0e+09
} else if (genome == "rnor6") {
  geneAnnotation <- readRDS(file.path(ArchRref, genome, "geneAnnotation.rds"))
  genomeAnnotation <- readRDS(
    file.path(ArchRref, genome, "genomeAnnotation.rds")
  )
  genomeSize <- 3.15e+09
} else {
  stop("Error : Organisms other than Mouse/Rat/Human not supported currently")
}

inputFiles <- list()

for (run in runs) {
  inputFiles[run[1]] <- run[3]
  names(inputFiles[run[1]]) <- run[1]
}
inputFiles <- unlist(inputFiles)
inputFiles

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  minTSS = min_tss,
  minFrags = min_frags,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  maxFrags = 1e+07,
  offsetPlus = 0,
  offsetMinus = 0,
  TileMatParams = list(tileSize = tile_size),
  force = TRUE
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "./Results",
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  copyArrows = FALSE
)
proj

getAvailableMatrices(proj)

samples <- sort(unique(proj$Sample))

bar_plot_a <- list()
nbc <- sqrt(dim(obj.lst_off[[1]])[2])

if (nbc == 96) {
  ncol <- 1
} else {
  ncol <- 2
}

for (i in seq_along(samples)){

  # Compute umi by row/col before filtering
  req_coord <- obj.lst_off[[i]]@meta.data
  meta.data_plots <- as.data.frame(
    getCellColData(proj)
  )[rownames(req_coord), ]
  meta.data_plots$col <- req_coord$col
  meta.data_plots$row <- req_coord$row

  # since number of row/col start with 0 we change it to start with 1
  meta.data_plots[, c("col", "row")] <- meta.data_plots[, c("col", "row")] + 1
  meta.data_plots$col <- as.character(meta.data_plots$col)
  meta.data_plots$row <- as.character(meta.data_plots$row)
  meta.data_plots$col <- factor(meta.data_plots$col, levels = 0:nbc)
  meta.data_plots$row <- factor(meta.data_plots$row, levels = 0:nbc)

  umi_row <- ggplot(
    data = meta.data_plots,
    aes(x = col, y = nFrags)) +
      geom_bar(stat = "identity", width = 0.75) +
      ylab("#nFrags / row") +
      xlab(paste("Row (", project_name, ")", sep = "")) +
      theme(text = element_text(size = 10))

  # umi_row
  umi_column <- ggplot(data = meta.data_plots, aes(x = row, y = nFrags)) +
    geom_bar(stat = "identity", width = 0.75) +
    ylab("#nFrags / column") +
    xlab(paste("Column (", project_name, ")", sep = "")) +
    theme(text = element_text(size = 10))

  options(repr.plot.width = 10, repr.plot.height = 5)
  bar_plot_a[[i]] <- wrap_plots(umi_row, umi_column, ncol = ncol)
}

names(bar_plot_a) <- names(obj.lst_off)

samples <- sort(unique(proj$Sample))

data <- list()

for (i in seq_along(samples)) {
  data[[i]] <- read.csv(
    paste0(run[4], "/tissue_positions_list.csv"),
    header = FALSE
  )
  colnames(data[[i]]) <- c(
    "barcode", "tissue", "row", "col", "imagerow", "imagecol"
  )
  data[[i]]$barcode <- paste0(samples[i], "#", data[[i]]$barcode, "-1")
}

all_tissue_positions <- do.call(rbind, data)

off_tissues_to_remove <- all_tissue_positions[
  which(all_tissue_positions$tissue == 0),
]$barcode

proj$barcodes1 <- rownames(getCellColData(proj))
proj$barcodes2 <- unname(unlist(genXtract(proj$barcodes1, "#", "-")))

proj2 <- proj[which(proj$barcodes1 %ni% off_tissues_to_remove)]

# remove barcodes with TSS < 1 and nFrags < 1000
proj2 <- proj2[which(proj2$TSSEnrichment > 1)]
proj2 <- proj2[which(proj2$nFrags > 1000)]

proj2

############### TSS and fragment distributions Plots

plot_frag_size <- plotFragmentSizes(ArchRProj = proj2) +
  labs(x = "Fragment Size (bp)") +
  theme(text = element_text(size = 21))

plot_tss_enrich <- plotTSSEnrichment(ArchRProj = proj2) +
  theme(text = element_text(size = 21))

df <- getCellColData(proj2, select = c("log10(nFrags)", "TSSEnrichment"))
plot_uniq_frag <- ggPoint(
  x = df[, 1],
  y = df[, 2],
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  baseSize = 12
)
layout <- "ABC"
frag_tss_Plot <- wrap_plots(
  A = plot_frag_size, B = plot_tss_enrich, C = plot_uniq_frag, design = layout
)

# dimension reduction and clustering
proj2 <- addIterativeLSI(
  ArchRProj = proj2,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = lsi_iterations,
  clusterParams = list(
    resolution = lsi_resolution,
    sampleCells = 20000,
    n.start = 10
  ),
  varFeatures = lsi_varfeatures,
  verbose = FALSE,
  dimsToUse = 1:30,
  force = TRUE
)

if (length(runs) > 1) {
  proj2 <- addHarmony(
    ArchRProj = proj2,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = TRUE
  )

  name <- "Harmony"
} else {
  name <- "IterativeLSI"
}

proj2 <- addClusters(
  input = proj2,
  reducedDims = name,
  method = "Seurat",
  name = "Clusters",
  resolution = resolution,
  force = TRUE
)

proj2 <- addUMAP(
  ArchRProj = proj2,
  reducedDims = name,
  name = "UMAP",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  force = TRUE
)

proj2 <- addImputeWeights(proj2)

proj2$Clusters <- gsub("C", "A", proj2$Clusters)

proj2 <- addGroupCoverages(
  ArchRProj = proj2,
  groupBy = "Clusters",
  maxCells = 1500,
  force = TRUE
)
pathToMacs2 <- findMacs2()

proj2 <- addReproduciblePeakSet(
  ArchRProj = proj2,
  groupBy = "Clusters",
  pathToMacs2 = pathToMacs2,
  genomeSize = genomeSize,
  maxPeaks = 300000,
  force = TRUE
)

getAvailableMatrices(proj2)
proj2 <- addPeakMatrix(proj2, force = TRUE)

# Save Archr project
saveArchRProject(
  ArchRProj = proj2, outputDirectory = "ArchrProject", load = FALSE
)

req_DF <- as.data.frame(getCellColData(proj2))
req_table <- melt(table(req_DF$Clusters, req_DF$Sample))
colnames(req_table) <- c("Cluster", "Sample", "%cells in knn clusters")
req_table$Cluster <- factor(
  req_table$Cluster,
  levels = (
    unique(
      req_table[order(as.numeric(gsub("C", "", req_table$Cluster))), ]$Cluster
    )
  )
)

pal <- paletteDiscrete(values = seuratObj$Sample)

sample_bar_plots_a <- ggplot(
  req_table, aes(fill = Sample, y = `%cells in knn clusters`, x = Cluster)
) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic() +
  scale_fill_manual(values = pal) +
  theme(text = element_text(size = 30)) +
  theme(axis.title.x = element_blank())


# get gene score matrix
seGS <- getMatrixFromProject(proj2)
rownames(seGS) <- rowData(seGS)$name

proj2 <- addImputeWeights(proj2)

# impute required gene score matrix
matGS <- imputeMatrix(assay(seGS), getImputeWeights(proj2))

meta.data <- as.data.frame(getCellColData(proj2))

samples <- sort(unique(proj2$Sample))

GS_mat <- list()
metaData <- list()
for (i in seq_along(samples)) {
  GS_mat[[i]] <- matGS[, which(
    colnames(matGS) %in% colnames(matGS)[grep(samples[i], colnames(matGS))]
  )]
  metaData[[i]] <- as.data.frame(
    getCellColData(proj2[which(proj2$Sample == samples[i])])
  )
}
names(GS_mat) <- samples
names(metaData) <- samples

spatial_in_tissue.obj <- list()
for (run in runs) {

  spatial_in_tissue.obj[[run[1]]] <- CreateSeuratObject(
    counts = as.data.frame(GS_mat[[run[1]]]),
    assay = "Spatial",
    meta.data = metaData[[run[1]]]
  )
  Idents(spatial_in_tissue.obj[[run[1]]]) <- factor(
    spatial_in_tissue.obj[[run[1]]]@meta.data$Clusters,
    levels = c(
      sort(unique(spatial_in_tissue.obj[[run[1]]]@meta.data$Clusters))
    )
  )

  image <- Read10X_Image(image.dir = run[4], filter.matrix = TRUE)
  rownames(image@coordinates) <- paste0(
    run[1], "#", rownames(image@coordinates), "-1"
  )
  sequenced_tixels <- row.names(meta.data)
  image <- image[sequenced_tixels, ]
  image <- image[Cells(x = spatial_in_tissue.obj[[run[1]]])]
  DefaultAssay(spatial_in_tissue.obj = image) <- "Spetial"
  spatial_in_tissue.obj[[run[1]]][[paste0(run[1], "_spatial")]] <- image
}

saveRDS(spatial_in_tissue.obj, "ATAC_lst.rds")
cat(paste0("RDS saved in ", getwd(), "/", "ATAC_lst.rds"))

# merge them
ctm <- c()
meta.data <- c()
for (i in names(spatial_in_tissue.obj)){
  tmp <- spatial_in_tissue.obj[[i]]@assays$Spatial@layers$counts
  colnames(tmp) <- colnames(spatial_in_tissue.obj[[i]])
  rownames(tmp) <- rownames(spatial_in_tissue.obj[[i]])
  tmp_met <- spatial_in_tissue.obj[[i]]@meta.data
  ctm <- cbind(ctm, tmp)
  meta.data <- rbind(meta.data, tmp_met)
}

spatial.obj <- CreateSeuratObject(
  counts = ctm, meta.data = meta.data, assay = "Spatial"
)
for (i in names(spatial_in_tissue.obj)){
  spatial.obj@images[[i]] <- spatial_in_tissue.obj[[i]]@images[[1]]
}
#### we run this but finally use the umap in ArchR project
spatial.obj <- NormalizeData(
  spatial.obj, normalization.method = "LogNormalize", verbose = FALSE
)
spatial.obj <- FindVariableFeatures(
  spatial.obj, selection.method = "vst", verbose = FALSE
)
all.genes <- rownames(spatial.obj)
spatial.obj <- ScaleData(spatial.obj, features = all.genes, verbose = FALSE)
spatial.obj <- RunPCA(spatial.obj, npcs = 5, verbose = FALSE)
spatial.obj <- FindNeighbors(spatial.obj, dims = 1:5, verbose = FALSE)
spatial.obj <- FindClusters(spatial.obj, resolution = 0.5, verbose = FALSE)
spatial.obj <- RunUMAP(
  spatial.obj, reduction = "pca", dims = 1:5, verbose = FALSE
)

##### now add umap cooridinations of ArchRproject to seurat object
spatial.obj@reductions$umap@cell.embeddings <- as.matrix(getEmbedding(proj2))
colnames(spatial.obj@reductions$umap@cell.embeddings) <- c("umap_1", "umap_2")
Idents(spatial.obj) <- spatial.obj$Clusters
spatial.obj

spatial.obj@meta.data$Clusters <- factor(
  spatial.obj@meta.data$Clusters,
  levels = sort(unique(spatial.obj@meta.data$Clusters))
)

Idents(spatial.obj) <- spatial.obj$Clusters

saveRDS(spatial.obj, "ATAC_all.rds")
cat(paste0("RDS saved in ", getwd(), "/", "ATAC_all.rds"))

n_clusters <- length(unique(spatial.obj$Clusters))
cols1 <- ArchRPalettes$stallion2[seq_len(n_clusters)]
names(cols1) <- unique(spatial.obj$Clusters)
n_samples <- length(unique(spatial.obj$Sample))

cols2 <- pal
names(cols2) <- unique(spatial.obj$Sample)

cols <- c(cols1, cols2)

# cols
dimrd3_a <- SpatialDimPlot(
  spatial.obj,
  group.by = "Clusters",
  pt.size.factor = point_size,
  cols = cols
)

dimrd1_a <- DimPlot(
  spatial.obj,
  reduction = "umap",
  group.by = c("Clusters", "Sample"),
  cols = cols
)

dimrd4_a <- list()
if (n_clusters >= 5) {
  ncol <- 5
} else {
  ncol <- n_clusters
}

for (i in names(spatial.obj@images)) {

  dimrd4_a[[i]] <- SpatialDimPlot(
    spatial.obj,
    cells.highlight = CellsByIdentities(object = spatial.obj),
    facet.highlight = TRUE,
    images = i,
    ncol = ncol
  )
}

# per cluster
markersGS <- getMarkerFeatures(
  ArchRProj = proj2,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "ttest"
)

threshold <- 0.1

if (threshold > 0) {
  markerList <- getMarkers(
    markersGS, cutOff = "Pval <= 0.05 & Log2FC >= threshold"
  )
} else {
  markerList <- getMarkers(
    markersGS, cutOff = "Pval <= 0.05 & Log2FC <= threshold"
  )
}

markerGenes <- list()
for (i in seq_len(length(markerList))) {
  markerGenes <- c(markerGenes, markerList[[i]]$name)
}

markerGenes <- unlist(markerGenes)

top_spatial_plots_a <- list()
for (i in samples){

  brain_ss[[i]] <- spatial.obj[, which(spatial.obj$Sample == i)]
  brain_ss[[i]]@images <- spatial.obj@images[i]
}

topGenes <- markerGenes[1:5]
topGenes <- topGenes[!is.na(topGenes)]

if (length(topGenes) >= 5) {
  ncol <- 5
} else {
  ncol <- length(topGenes)
}

top_list <- list()
for (i in seq_along(samples)) {
  for (j in topGenes) {
    top_list[[j]] <- SpatialFeaturePlot(
      brain_ss[[i]],
      features = j,
      ncol = 1,
      alpha = c(.8, 2),
      pt.size.factor = point_size) +
        theme(
          legend.position = "top",
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 20)
        )
    top_list[[j]]$layers[[1]]$aes_params <- c(
      top_list[[j]]$layers[[1]]$aes_params, shape = 22
    )

    top_spatial_plots_a[[i]] <- wrap_plots(top_list, ncol = ncol)
  }
}

names(top_spatial_plots_a) <- samples

######### linking peaks to genes
# for having peak to gene links we have to make a new atac object
# like what signac does

RNA <- seuratObj

for (run in runs) {
  # first add run id to fragments barcodes

  rawdata <- run[3]
  dataPath <- "/root"
  file.copy(rawdata, dataPath, overwrite = TRUE)

  frag <- fread("/root/fragments.tsv.gz")
  frag$V4 <- paste0(run[1], "#", frag$V4)

  write.table(
    frag,
    paste0(run[1], "_fragments.tsv"),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )

  system("rm /root/*fragments.tsv.gz")

  system(
    paste(
      "sort -k 1,1V -k2,2n",
      paste0(run[1], "_fragments.tsv"),
      ">",
      paste0(run[1], "_sorted_fragments.tsv")
    ),
    wait = TRUE
  )

  system(paste("rm", paste0(run[1], "_fragments.tsv")), wait = TRUE)
}

# merge and sort all fragments together
system("sort -m -k 1,1V -k2,2n *_sorted_fragments.tsv > merged_fragments.tsv")

system("bgzip -@ 4 merged_fragments.tsv")     # -@ 4 uses 4 threads

system("tabix -p bed merged_fragments.tsv.gz")

# now create signac object
fragpath <- "/root/reports/merged_fragments.tsv.gz"
counts <- CountFragments(fragments = fragpath, verbose = FALSE)
frags <- CreateFragmentObject(path = fragpath, cells = counts$CB)

species <- args[2]

if (species == "mm10") {
  genome1 <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)
} else if (species == "hg38") {
  genome1 <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
} else {
  genome1 <- seqlengths(BSgenome.Rnorvegicus.UCSC.rn6)
}

bin_matrix <- GenomeBinMatrix(
  fragments = frags,
  genome = genome1,
  binsize = 5000,
  verbose = TRUE
)

if (species == "mm10") {
  genome2 <- "mm10"
}else if (species == "hg38") {
  genome2 <= "hg38"
} else {
  genome2 <- "rn6"
}

chrom <- CreateChromatinAssay(
  counts = bin_matrix,
  genome = genome2,
  fragments = fragpath
)

signac_obj <- CreateSeuratObject(
  counts = chrom,
  assay = "ATAC"
)

# call peaks using MACS2
peaks <- CallPeaks(signac_obj)

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(signac_obj),
  features = peaks,
  cells = colnames(signac_obj)
)

if (species == "mm10") {
  ensdb <- EnsDb.Mmusculus.v79
} else {
  ensdb <- EnsDb.Hsapiens.v86
}

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = ensdb)
seqlevels(annotation) <- paste0("chr", seqlevels(annotation))

# create a new assay using the MACS2 peak set and add it to the Seurat object
signac_obj[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)
signac_obj
# keep only samples that include in RNA seq
signac_obj <- signac_obj[, colnames(signac_obj) %in% colnames(RNA)]

# no QC required
# Normalization and linear dimensional reduction
DefaultAssay(signac_obj) <- "peaks"

signac_obj <- RunTFIDF(signac_obj)
signac_obj <- FindTopFeatures(signac_obj, min.cutoff = "q0")
signac_obj <- RunSVD(object = signac_obj)

# we don't remove any lsi because there might be some information

signac_obj <- RunUMAP(
  object = signac_obj,
  reduction = "lsi",
  dims = 1:30
)
signac_obj <- FindNeighbors(
  object = signac_obj,
  reduction = "lsi",
  dims = 1:30
)

signac_obj <- FindClusters(
  object = signac_obj,
  algorithm = 1,
  resolution = 0.5,
  verbose = FALSE
)

# we don't add gene activity score but use our rna data
counts <- RNA@assays$SCT@counts

# create a new assay to store ADT information
rna_assay <- CreateAssay5Object(counts = counts)

# add this assay to the previously created Seurat object
signac_obj[["RNA"]] <- rna_assay

DefaultAssay(signac_obj) <- "RNA"
signac_obj <- SCTransform(signac_obj)
signac_obj <- RunPCA(signac_obj)

# since rna clustering usually is better than ATAC one, add rna clusters to
# comapping signac object
df <- signac_obj@meta.data

df <- df[match(rownames(df), colnames(RNA)), ]

df$RNA_clusters <- RNA@meta.data$seurat_clusters

df <- df[match(colnames(signac_obj), rownames(df)), ]

signac_obj@meta.data <- df

# linking peaks to genes
genome <- args[2]

if (genome == "hg38") {
  bsgenome <- BSgenome.Hsapiens.UCSC.hg38
} else {
  bsgenome <- BSgenome.Mmusculus.UCSC.mm10
}
DefaultAssay(signac_obj) <- "peaks"

# first compute the GC content for each peak
signac_obj <- RegionStats(signac_obj, genome = bsgenome)
genes.use <- top.features

# link peaks to genes
signac_obj <- LinkPeaks(
  object = signac_obj,
  peak.assay = "peaks",
  expression.assay = "SCT",
  genes.use = genes.use
)

saveRDS(signac_obj, "signac_obj.rds")

Idents(signac_obj) <- "RNA_clusters"

clusters2 <- sort(unique(signac_obj$RNA_clusters))
colors2 <- ArchRPalettes$stallion2[seq_len(length(clusters2))]
names(colors2) <- clusters2

idents.plot <- clusters2
top_coverage_plots_r <- list()
topGenes <- unique(signac_obj@assays[["peaks"]]@links$gene)

for (j in topGenes){
  top_coverage_plots_r[[j]] <- CoveragePlot(
    object = signac_obj,
    region = j,
    features = j,
    expression.assay = "SCT",
    idents = idents.plot,
    extend.upstream = 10000,
    extend.downstream = 10000) & scale_fill_manual(values = colors2)
}

pdf("linkpeakstogenes.pdf", width = 15, height = 9)
print(top_coverage_plots_r)
dev.off()


# add heatmap of integration
RNA <- readRDS("/root/reports/RNA_all.rds")
cells <- intersect(rownames(getCellColData(proj2)), colnames(RNA))

proj3 <- subsetArchRProject(
  ArchRProj = proj2,
  cells = cells,
  outputDirectory = "ArchRSubset",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)

RNA2 <- RNA[, colnames(RNA) %in% rownames(getCellColData(proj3))]
ct <- as.matrix(RNA2@assays$SCT@counts)
colnames(ct) <- colnames(RNA2)
rownames(ct) <- rownames(RNA2)

ct <- ct[, match(rownames(getCellColData(proj3)), colnames(ct))]
meta.data <- RNA2@meta.data[match(
  rownames(getCellColData(proj3)), rownames(RNA2@meta.data)
), ]

RNA_reorder <- CreateSeuratObject(
  counts = as.data.frame(ct),
  assay = "RNA",
  meta.data = meta.data
)
RNA_reorder

# testIntegration() from /root/addgeneintegrationmatrix.R
proj3 <- testIntegration(
  ArchRProj = proj3,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = RNA_reorder,
  addToArrow = TRUE,
  force = TRUE,
  groupRNA = "seurat_clusters",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)

getAvailableMatrices(proj3)

proj3 <- addPeak2GeneLinks(
  ArchRProj = proj3,
  useMatrix = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI"
)

heatmap <- plotPeak2GeneHeatmap(
  ArchRProj = proj3,
  groupBy = "predictedGroup",
  k = 10
)

##########################################
report_name <- paste0(project_name, "_report.Rmd")

file.copy("/root/HTML_report_knit.Rmd", report_name, overwrite = TRUE)

knitr::opts_chunk$set(echo = FALSE, warning = FALSE)
rmarkdown::render(report_name)
