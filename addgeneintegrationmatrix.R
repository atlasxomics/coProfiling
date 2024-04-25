library("ArchR")
library("dplyr")
library("ggplot2")
library("Matrix")
library("matrixStats")
library("parallel")
library("pdftools")
library("rhdf5")
library("Seurat")
library("stringr")
library("SummarizedExperiment")
library("S4Vectors")
library("utils")
library("uwot")

testIntegration <- function(
  ArchRProj = NULL,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = NULL,
  groupATAC = NULL,
  groupRNA = NULL,
  groupList = NULL,
  sampleCellsATAC = 10000,
  sampleCellsRNA = 10000,
  embeddingATAC = NULL,
  embeddingRNA = NULL,
  dimsToUse = 1:30,
  scaleDims = NULL,
  corCutOff = 0.75,
  plotUMAP = TRUE,
  UMAPParams = list(
    n_neighbors = 40,
    min_dist = 0.4,
    metric = "cosine",
    verbose = FALSE
  ),
  nGenes = 2000,
  useImputation = TRUE,
  reduction = "cca",
  addToArrow = TRUE,
  scaleTo = 10000,
  genesUse = NULL,
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  transferParams = list(),
  threads = getArchRThreads(),
  verbose = TRUE,
  force = FALSE,
  logFile = createLogFile("addGeneIntegrationMatrix"),
  ...
) {
  if (is.null(groupList)) {
    groupList <- SimpleList()
    groupList[[1]] <- SimpleList(
      ATAC = ArchRProj$cellNames,
      RNA = colnames(seRNA)
    )
  }

  if (useMatrix %ni% getAvailableMatrices(ArchRProj)) {
    stop("Matrix name provided to useMatrix does not exist in ArchRProject!")
  }

  if (!is.null(groupATAC)) {
    dfATAC <- getCellColData(
      ArchRProj = ArchRProj, select = groupATAC, drop = FALSE
    )
  }

  nCell <- rep(0, length(ArchRProj$cellNames))
  names(nCell) <- ArchRProj$cellNames
  groupList <- lapply(seq_along(groupList), function(x) {
    ATAC <- groupList[[x]]$ATAC
    if (!is.null(groupATAC)) {
      if (any(ATAC %in% dfATAC[, 1])) {
        idx <- which(ATAC %in% dfATAC[, 1])
        ATAC2 <- rownames(dfATAC)[which(dfATAC[, 1] %in% ATAC[idx])]
        if (length(idx) == length(ATAC)) {
          ATAC <- ATAC2
        }
        else {
          ATAC <- c(ATAC[-idx], ATAC2)
        }
      }
    }
    SimpleList(ATAC = ATAC, RNA = groupList[[x]]$RNA)
  }) %>% SimpleList

  for (i in seq_along(groupList)) {
    nCell[groupList[[i]]$ATAC] <- nCell[groupList[[i]]$ATAC] + 1
  }

  if (!all(nCell == 1)) {
    stop(
      "Missing ",
      length(which(nCell == 0)),
      " cells. Found ",
      length(which(nCell > 1)),
      " overlapping cells from ArchRProj in groupList! Cannot have \
      overlapping/missing cells in ATAC input, check 'groupList' argument!"
    )
  }

  if (inherits(seRNA, "SummarizedExperiment")) {
    seuratRNA <- CreateSeuratObject(counts = assay(seRNA))
    if (groupRNA %ni% colnames(colData(seRNA))) {
      stop("groupRNA not in colData of seRNA")
    }
    seuratRNA$Group <- paste0(colData(seRNA)[, groupRNA, drop = TRUE])
    rm(seRNA)
  } else {
    if (groupRNA %ni% colnames(seRNA@meta.data)) {
      stop("groupRNA not in meta.data of Seurat Object")
    }
    seuratRNA <- seRNA
    seuratRNA$Group <- paste0(seRNA@meta.data[, groupRNA])
    rm(seRNA)
  }

  if ("RNA" %in% names(seuratRNA@assays)) {
    DefaultAssay(seuratRNA) <- "RNA"
  } else {
    stop(
      "'RNA' is not present in Seurat Object's Assays! Please make sure \
      that this assay is present!"
    )
  }

  gc()
  if (!is.null(groupRNA)) {
    dfRNA <- DataFrame(row.names = colnames(seuratRNA), Group = seuratRNA$Group)
  }

  groupList <- lapply(seq_along(groupList), function(x) {
    RNA <- groupList[[x]]$RNA
    if (!is.null(groupRNA)) {
      if (any(RNA %in% dfRNA[, 1])) {
        idx <- which(RNA %in% dfRNA[, 1])
        RNA2 <- rownames(dfRNA)[which(dfRNA[, 1] %in% RNA[idx])]
        if (length(idx) == length(RNA)) {
          RNA <- RNA2
        }
        else {
          RNA <- c(RNA[-idx], RNA2)
        }
      }
    }
    SimpleList(ATAC = groupList[[x]]$ATAC, RNA = RNA)
  }) %>% SimpleList

  cellRNA <- unlist(lapply(groupList, function(x) x$RNA))
  if (!all(cellRNA %in% colnames(seuratRNA))) {
    stop("Found cells for RNA not in colnames(seRNA)! Please retry your input!")
  }

  seuratRNA <- seuratRNA[, unique(cellRNA)]
  seuratRNA <- NormalizeData(object = seuratRNA, verbose = FALSE)
  geneDF <- ArchR:::.getFeatureDF(getArrowFiles(ArchRProj), useMatrix)
  sumOverlap <- sum(unique(geneDF$name) %in% unique(rownames(seuratRNA)))
  if (sumOverlap < 5) {
    stop("Error not enough overlaps (",
    sumOverlap,
    ") between gene names from gene scores (ArchR) and rna matrix (seRNA)!")
  }

  blockList <- SimpleList()
  for (i in seq_along(groupList)) {
    gLi <- groupList[[i]]
    if (length(gLi$ATAC) > sampleCellsATAC) {
      if (!is.null(embeddingATAC)) {
        probATAC <- ArchR:::.getDensity(
          embeddingATAC[gLi$ATAC, 1],
          embeddingATAC[gLi$ATAC, 2]
        )$density
        probATAC <- probATAC / max(probATAC)
        cellsATAC <- gLi$ATAC[order(probATAC, decreasing = TRUE)]
      } else {
        cellsATAC <- sample(gLi$ATAC, length(gLi$ATAC))
      }
      cutoffs <- lapply(seq_len(1000), function(x) length(gLi$ATAC) / x) %>%
        unlist
      blockSize <- ceiling(
        min(cutoffs[order(abs(cutoffs - sampleCellsATAC))[1]] + 1,
        length(gLi$ATAC))
      )
      nBlocks <- ceiling(length(gLi$ATAC) / blockSize)
      blocks <- lapply(seq_len(nBlocks), function(x) {
        cellsATAC[seq(x, length(cellsATAC), nBlocks)]
      }) %>% SimpleList
    } else {
      blocks <- list(gLi$ATAC)
    }
    if (!is.null(embeddingRNA)) {
      probRNA <- ArchR:::.getDensity(
        embeddingRNA[gLi$RNA, 1],
        embeddingRNA[gLi$RNA, 2]
      )$density
      probRNA <- probRNA / max(probRNA)
    } else {
      probRNA <- rep(1, length(gLi$RNA))
    }
    blockListi <- lapply(
      seq_along(blocks), function(x) {
        SimpleList(
          ATAC = blocks[[x]],
          RNA = sample(x = gLi$RNA, size = min(sampleCellsRNA, length(gLi$RNA)),
          prob = probRNA)
        )
      }
    ) %>% SimpleList
    blockList <- c(blockList, blockListi)
  }

  rm(groupList)
  subProj <- ArchRProj
  subProj@imputeWeights <- SimpleList()
  geneDF <- ArchR:::.getFeatureDF(getArrowFiles(subProj), useMatrix)
  geneDF <- geneDF[geneDF$name %in% rownames(seuratRNA), , drop = FALSE]
  splitGeneDF <- S4Vectors::split(geneDF, geneDF$seqnames)
  featureDF <- lapply(splitGeneDF, function(x) {
    x$idx <- seq_len(nrow(x))
    return(x)
  }) %>% Reduce("rbind", .)

  dfParams <- data.frame(reduction = reduction)
  allChr <- unique(featureDF$seqnames)
  tmpFile <- ArchR:::.tempfile()
  o <- suppressWarnings(
    file.remove(paste0(
      tmpFile,
      "-IntegrationBlock-",
      seq_along(blockList),
      ".h5")
    )
  )

  h5lock <- addArchRLogging(useLogs = TRUE)
  if (h5lock) {
    if (threads > 1) {
      message(
        "subThreading Disabled since ArchRLocking is TRUE see `addArchRLocking`"
      )
      threads <- 1
    }
  } else {
    if (threads > 1) {
      message(
        "subThreading Enabled since ArchRLocking is FALSE see `addArchRLocking`"
      )
    }
  }

  rD <- getReducedDims(ArchRProj = ArchRProj, reducedDims = reducedDims,
                       corCutOff = corCutOff, dimsToUse = dimsToUse)
  outDir1 <- getOutputDirectory(ArchRProj)
  outDir2 <- file.path(outDir1, "RNAIntegration")
  outDir3 <- file.path(outDir2, matrixName)
  dir.create(outDir1, showWarnings = FALSE)
  dir.create(outDir2, showWarnings = FALSE)
  dir.create(outDir3, showWarnings = FALSE)
  prevFiles <- list.files(outDir3, full.names = TRUE)
  prevFiles <- ArchR:::.suppressAll(file.remove(prevFiles))

  threads2 <- max(ceiling(threads * 0.75), 1)
  dfAll <- ArchR:::.safelapply(seq_along(blockList), function(i) {
    prefix <- sprintf("Block (%s of %s) :", i, length(blockList))
    blocki <- blockList[[i]]
    subProj@cellColData <- subProj@cellColData[blocki$ATAC, ]
    subProj@sampleColData <- subProj@sampleColData[
      unique(subProj$Sample), , drop = FALSE
    ]
    subRNA <- seuratRNA[, blocki$RNA]
    subRNA <- subRNA[rownames(subRNA) %in% geneDF$name, ]
    subRNA <- FindVariableFeatures(
      object = subRNA, nfeatures = nGenes, verbose = FALSE
    )
    subRNA <- ScaleData(object = subRNA, verbose = FALSE)
    if (is.null(genesUse)) {
      genesUse <- VariableFeatures(object = subRNA)
    }

    mat <- ArchR:::.getPartialMatrix(
      getArrowFiles(subProj),
      featureDF = geneDF[geneDF$name %in% genesUse, ],
      threads = 2,
      cellNames = subProj$cellNames,
      useMatrix = useMatrix,
      verbose = FALSE
    )
    rownames(mat) <- geneDF[geneDF$name %in% genesUse, "name"]

    if (useImputation) {
      imputeParams <- list()
      imputeParams$ArchRProj <- subProj
      imputeParams$randomSuffix <- TRUE
      imputeParams$reducedDims <- reducedDims
      imputeParams$dimsToUse <- dimsToUse
      imputeParams$scaleDims <- scaleDims
      imputeParams$corCutOff <- corCutOff
      imputeParams$threads <- 1
      subProj <- suppressMessages(do.call(addImputeWeights, imputeParams))
      mat <- suppressMessages(
        imputeMatrix(mat = mat,
        imputeWeights = getImputeWeights(subProj),
        verbose = FALSE)
      )
      o <- suppressWarnings(file.remove(unlist(getImputeWeights(subProj)[[1]])))
    }

    mat <- log(mat + 1)
    colnames(mat) <- as.character(colnames(mat))
    rownames(mat) <- as.character(rownames(mat))

    seuratATAC <- Seurat::CreateSeuratObject(
      counts = mat[head(seq_len(nrow(mat)), 5), , drop = FALSE]
    )
    seuratATAC[["GeneScore"]] <- Seurat::CreateAssayObject(counts = mat)
    rm(mat)
    DefaultAssay(seuratATAC) <- "GeneScore"
    seuratATAC <- Seurat::ScaleData(seuratATAC, verbose = FALSE)
    transferAnchors <- ArchR:::.retryCatch(
      {
        gc()
        Seurat::FindTransferAnchors(
          reference = subRNA,
          query = seuratATAC,
          reduction = reduction,
          features = genesUse,
          verbose = FALSE,
          ...
        )
      },
      maxAttempts = 2,
      logFile = logFile
    )

    rDSub <- rD[colnames(seuratATAC), , drop = FALSE]

    transferParams$anchorset <- transferAnchors
    transferParams$weight.reduction <- CreateDimReducObject(
      embeddings = rDSub, key = "LSI_", assay = DefaultAssay(seuratATAC)
    )
    transferParams$verbose <- FALSE
    transferParams$dims <- seq_len(ncol(rDSub))

    transferParams$refdata <- subRNA$Group
    rnaLabels <- do.call(Seurat::TransferData, transferParams)

    transferParams$refdata <- colnames(subRNA)
    rnaLabels2 <- do.call(Seurat::TransferData, transferParams)[, 1]
    if (addToArrow) {
      transferParams$refdata <- GetAssayData(
        subRNA, assay = "RNA", slot = "data"
      )
      gc()
      matchedRNA <- do.call(Seurat::TransferData, transferParams)
      matchedRNA <- matchedRNA@data
    }
    matchDF <- DataFrame(
      cellNames = colnames(seuratATAC),
      predictionScore = rnaLabels$prediction.score.max,
      predictedGroup = rnaLabels$predicted.id,
      predictedCell = rnaLabels2
    )
    rownames(matchDF) <- matchDF$cellNames

    jointCCA <- DataFrame(
      transferAnchors@object.list[[1]]@reductions$cca@cell.embeddings
    )
    jointCCA$Assay <- ifelse(
      endsWith(rownames(jointCCA), "_reference"), "RNA", "ATAC"
    )
    jointCCA$Group <- NA
    jointCCA$Score <- NA
    jointCCA[paste0(colnames(subRNA), "_reference"), "Group"] <- subRNA$Group
    jointCCA[paste0(matchDF$cellNames, "_query"), "Group"] <- matchDF$predictedGroup
    jointCCA[paste0(matchDF$cellNames, "_query"), "Score"] <- matchDF$predictionScore
    ArchR:::.safeSaveRDS(
      object = jointCCA,
      file = file.path(
        outDir3,
        paste0("Save-Block", i, "-JointCCA.rds")
      )
    )
    rm(transferParams, transferAnchors)
    gc()
    if (addToArrow) {

      tmpFilei <- paste0(tmpFile, "-IntegrationBlock-", i, ".h5")
      o <- h5createFile(tmpFilei)
      sampleNames <- getCellColData(subProj, "Sample")[matchDF$cellNames, ]
      uniqueSamples <- unique(sampleNames)
      matchedRNA <- ArchR:::.safeSubset(
        mat = matchedRNA,
        subsetRows = paste0(featureDF$name),
        subsetCols = matchDF$cellNames
      )
      for (z in seq_along(uniqueSamples)) {
        mat <- matchedRNA[
          ,
          which(sampleNames == uniqueSamples[z]),
          drop = FALSE
        ]
        Group <- uniqueSamples[z]
        o <- tryCatch(
          {
            h5delete(tmpFilei, paste0(Group))
          },
          error = function(x) {
          }
        )

        o <- h5createGroup(tmpFilei, paste0(Group))
        j <- Rle(findInterval(seq(mat@x) - 1, mat@p[-1]) + 1)
        lengthRle <- length(j@lengths)
        lengthI <- length(mat@i)
        o <- ArchR:::.suppressAll(h5createDataset(
          tmpFilei,
          paste0(Group, "/i"),
          storage.mode = "integer",
          dims = c(lengthI, 1), level = 0)
        )
        o <- ArchR:::.suppressAll(
          h5createDataset(
            tmpFilei,
            paste0(Group, "/jLengths"),
            storage.mode = "integer",
            dims = c(lengthRle, 1),
            level = 0
          )
        )
        o <- ArchR:::.suppressAll(
          h5createDataset(
            tmpFilei,
            paste0(Group, "/jValues"),
            storage.mode = "integer",
            dims = c(lengthRle, 1),
            level = 0
          )
        )
        o <- ArchR:::.suppressAll(
          h5createDataset(
            tmpFilei,
            paste0(Group, "/x"),
            storage.mode = "double",
            dims = c(lengthI, 1),
            level = 0
          )
        )
        o <- ArchR:::.suppressAll(
          h5write(obj = mat@i + 1, file = tmpFilei, name = paste0(Group, "/i"))
        )
        o <- ArchR:::.suppressAll(
          h5write(
            obj = j@lengths, file = tmpFilei, name = paste0(Group, "/jLengths")
          )
        )
        o <- ArchR:::.suppressAll(
          h5write(
            obj = j@values, file = tmpFilei, name = paste0(Group, "/jValues")
          )
        )
        o <- ArchR:::.suppressAll(
          h5write(obj = mat@x, file = tmpFilei, name = paste0(Group, "/x"))
        )
        o <- ArchR:::.suppressAll(
          h5write(
            obj = colnames(mat),
            file = tmpFilei,
            name = paste0(Group, "/cellNames")
          )
        )
      }
      rm(matchedRNA, mat, j)
    }

    gc()
    matchDF$Block <- Rle(i)
    matchDF
  }, threads = threads2) %>% Reduce("rbind", .)

  if (plotUMAP) {
    for (i in seq_along(blockList)) {
      o <- tryCatch({
        prefix <- sprintf("Block (%s of %s) :", i, length(blockList))

        jointCCA <- readRDS(
          file.path(outDir3, paste0("Save-Block", i, "-JointCCA.rds"))
        )
        set.seed(1)
        UMAPParams <- .mergeParams(
          UMAPParams,
          list(
            n_neighbors = 40,
            min_dist = 0.4,
            metric = "cosine",
            verbose = FALSE
          )
        )
        UMAPParams$X <- as.data.frame(
          jointCCA[, grep("CC_", colnames(jointCCA))]
        )
        UMAPParams$ret_nn <- FALSE
        UMAPParams$ret_model <- FALSE
        UMAPParams$n_threads <- 1
        uwotUmap <- tryCatch({
          do.call(uwot::umap, UMAPParams)
        }, error = function(e) {
          errorList <- UMAPParams
        })

        jointCCA$UMAP1 <- uwotUmap[, 1]
        jointCCA$UMAP2 <- uwotUmap[, 2]
        .safeSaveRDS(
          object = jointCCA,
          file = file.path(outDir3, paste0("Save-Block", i, "-JointCCA.rds"))
        )
        p1 <- ggPoint(
          x = uwotUmap[, 1],
          y = uwotUmap[, 2],
          color = jointCCA$Assay,
          randomize = TRUE,
          size = 0.2,
          title = paste0(prefix, " colored by Assay"),
          xlabel = "UMAP Dimension 1",
          ylabel = "UMAP Dimension 2",
          rastr = TRUE
        ) +
          theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
          )
        p2 <- ggPoint(
          x = uwotUmap[, 1],
          y = uwotUmap[, 2],
          color = jointCCA$Group,
          randomize = TRUE,
          size = 0.2,
          title = paste0(prefix, " colored by scRNA Group"),
          xlabel = "UMAP Dimension 1",
          ylabel = "UMAP Dimension 2",
          rastr = TRUE
        ) +
          theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
          )
        pdf(
          file.path(outDir3, paste0("Save-Block", i, "-JointCCA-UMAP.pdf")),
          width = 12,
          height = 6,
          useDingbats = FALSE
        )
        ggAlignPlots(p1, p2, type = "h")
        dev.off()
      },
      error = function(e) {
      }
      )
    }
  }
  if (addToArrow) {
    matrixName <- ArchR:::.isProtectedArray(matrixName)
    integrationFiles <- paste0(
      tmpFile,
      "-IntegrationBlock-",
      seq_along(blockList), ".h5"
    )
    if (!all(file.exists(integrationFiles))) {
      stop(
        "Something went wrong with integration as not all temporary files \
        containing integrated RNA exist!"
      )
    }

    h5list <- ArchR:::.safelapply(
      seq_along(integrationFiles), function(x) {
        h5ls(integrationFiles[x])
      },
      threads = threads
    )

    ArrowFiles <- getArrowFiles(ArchRProj)
    allSamples <- names(ArrowFiles)
    o <- .safelapply(
      seq_along(allSamples), function(y) {
        sample <- allSamples[y]
        prefix <- sprintf("%s (%s of %s)", sample, y, length(ArrowFiles))

        sampleIF <- lapply(seq_along(h5list), function(x) {
          if (any(h5list[[x]]$group == paste0("/", sample))) {
            integrationFiles[x]
          } else {
            NULL
          }
        }) %>% unlist
        sampleMat <- lapply(
          seq_along(sampleIF), function(x) {
            cellNames <- ArchR:::.h5read(
              sampleIF[x], paste0(sample, "/cellNames")
            )
            mat <- sparseMatrix(
              i = ArchR:::.h5read(sampleIF[x], paste0(sample, "/i"))[, 1],
              j = as.vector(
                Rle(
                  ArchR:::.h5read(sampleIF[x], paste0(sample, "/jValues"))[, 1],
                  ArchR:::.h5read(sampleIF[x], paste0(sample, "/jLengths"))[, 1]
                )
              ),
              x = ArchR:::.h5read(
                sampleIF[x],
                paste0(sample, "/x")
              )[, 1],
              dims = c(nrow(featureDF), length(cellNames))
            )
            colnames(mat) <- cellNames
            mat
          }
        ) %>% Reduce("cbind", .)

        sampleMat@x <- exp(sampleMat@x) - 1
        sampleMat <- ArchR:::.normalizeCols(sampleMat, scaleTo = scaleTo)
        sampleMat <- drop0(sampleMat)
        rownames(sampleMat) <- paste0(featureDF$name)
        sampleMat <- sampleMat[
          ,
          ArchRProj$cellNames[BiocGenerics::which(ArchRProj$Sample == sample)],
          drop = FALSE
        ]
        o <- ArchR:::.createArrowGroup(
          ArrowFile = ArrowFiles[sample], group = matrixName, force = force
        )
        o <- ArchR:::.initializeMat(
          ArrowFile = ArrowFiles[sample],
          Group = matrixName,
          Class = "double",
          Units = "NormCounts",
          cellNames = colnames(sampleMat),
          params = dfParams,
          featureDF = featureDF,
          force = force
        )
        o <- h5write(
          obj = dfAll[colnames(sampleMat), "predictionScore"],
          file = ArrowFiles[sample],
          name = paste0(
            matrixName,
            "/Info/predictionScore"
          )
        )
        o <- h5write(
          obj = dfAll[colnames(sampleMat), "predictedGroup"],
          file = ArrowFiles[sample],
          name = paste0(
            matrixName,
            "/Info/predictedGroup"
          )
        )
        o <- h5write(
          obj = dfAll[colnames(sampleMat), "predictedCell"],
          file = ArrowFiles[sample],
          name = paste0(
            matrixName,
            "/Info/predictedCell"
          )
        )

        for (z in seq_along(allChr)) {
          chrz <- allChr[z]

          idz <- BiocGenerics::which(featureDF$seqnames %bcin% chrz)
          matz <- sampleMat[idz, , drop = FALSE]
          stopifnot(identical(paste0(featureDF$name[idz]),
                              paste0(rownames(matz))))
          o <- ArchR:::.addMatToArrow(
            mat = matz,
            ArrowFile = ArrowFiles[sample],
            Group = paste0(matrixName, "/", chrz),
            binarize = FALSE,
            addColSums = TRUE,
            addRowSums = TRUE,
            addRowVarsLog2 = TRUE,
            logFile = logFile
          )
          rm(matz)
          if (z %% 3 == 0 | z == length(allChr)) {
            gc()
          }
        }
        0
      },
      threads = threads
    )
    o <- suppressWarnings(file.remove(integrationFiles))
  }
  ArchRProj <- addCellColData(
    ArchRProj = ArchRProj,
    cells = dfAll$cellNames,
    data = dfAll$predictedCell,
    name = nameCell,
    force = TRUE
  )
  ArchRProj <- addCellColData(
    ArchRProj = ArchRProj,
    cells = dfAll$cellNames,
    data = dfAll$predictedGroup,
    name = nameGroup,
    force = TRUE
  )
  ArchRProj <- addCellColData(
    ArchRProj = ArchRProj,
    cells = dfAll$cellNames,
    data = dfAll$predictionScore,
    name = nameScore,
    force = TRUE
  )

  return(ArchRProj)
}

###############################################################################
# Helper Intermediate Methods
###############################################################################

.mergeParams <- function(paramInput = NULL, paramDefault = NULL) {
  for (i in seq_along(paramDefault)){
    if (!(names(paramDefault)[i] %in% names(paramInput))) {
      paramInput[[names(paramDefault)[i]]] <- paramDefault[[i]]
    }
  }
  return(paramInput)
}

.requirePackage <- function(
  x = NULL, load = TRUE, installInfo = NULL, source = NULL
) {
  if (x %in% rownames(installed.packages())) {
    if (load) {
      suppressPackageStartupMessages(require(x, character.only = TRUE))
    } else {
      return(0)
    }
  } else {
    if (!is.null(source) & is.null(installInfo)) {
      if (tolower(source) == "cran") {
        installInfo <- paste0('install.packages("', x, '")')
      }else if (tolower(source) == "bioc") {
        installInfo <- paste0('BiocManager::install("', x, '")')
      } else {
        stop("Unrecognized package source, available are cran/bioc!")
      }
    }
    if (!is.null(installInfo)) {
      stop(
        paste0(
          "Required package : ",
          x,
          " is not installed/found!\n  Package Can Be Installed : ",
          installInfo
        )
      )
    } else {
      stop(paste0("Required package : ", x, " is not installed/found!"))
    }
  }
}

###############################################################################
# Safe saveRDS check
###############################################################################

.safeSaveRDS <- function(
  object = NULL,
  file = "",
  ascii = FALSE,
  version = NULL,
  compress = TRUE,
  refhook = NULL
) {
  #Try to save a test data.frame in location
  testDF <- data.frame(a = 1, b = 2)
  canSave <- suppressWarnings(
    tryCatch(
      {
        saveRDS(
          object = testDF,
          file = file,
          ascii = ascii,
          version = version,
          compress = compress,
          refhook = refhook
        )
        TRUE
      },
      error = function(x) {
      FALSE
      }
    )
  )
  if (!canSave) {
    dirExists <- dir.exists(dirname(file))
    if (dirExists) {
      stop("Cannot saveRDS. File Path : ", file)
    } else {
      stop(
        "Cannot saveRDS because directory does not exist (",
        dirname(file),
        "). File Path : ",
        file
      )
    }
  } else {
    saveRDS(
      object = object,
      file = file,
      ascii = ascii,
      version = version,
      compress = compress,
      refhook = refhook
    )
  }
}

###############################################################################
# Stat/Summary Methods
###############################################################################

.computeKNN <- function(
  data = NULL,
  query = NULL,
  k = 50,
  includeSelf = FALSE,
  ...
) {
  .validInput(input = data, name = "data", valid = c("dataframe", "matrix"))
  .validInput(input = query, name = "query", valid = c("dataframe", "matrix"))
  .validInput(input = k, name = "k", valid = c("integer"))
  .validInput(input = includeSelf, name = "includeSelf", valid = c("boolean"))
  if (is.null(query)) {
    query <- data
    searchSelf <- TRUE
  } else {
    searchSelf <- FALSE
  }
  .requirePackage("nabor", source = "cran")
  if (searchSelf & !includeSelf) {
    knnIdx <- nabor::knn(data = data, query = query, k = k + 1, ...)$nn.idx
    knnIdx <- knnIdx[, -1, drop = FALSE]
  } else {
    knnIdx <- nabor::knn(data = data, query = query, k = k, ...)$nn.idx
  }
  knnIdx
}

.rowZscores <- function(m = NULL, min = -2, max = 2, limit = FALSE) {
  z <- sweep(m - rowMeans(m), 1, matrixStats::rowSds(m), `/`)
  if (limit) {
    z[z > max] <- max
    z[z < min] <- min
  }
  return(z)
}

.computeROC <- function(labels = NULL, scores = NULL, name = "ROC") {
  .calcAUC <- function(TPR = NULL, FPR = NULL) {
    # http://blog.revolutionanalytics.com/2016/11/calculating-auc.html
    dFPR <- c(diff(FPR), 0)
    dTPR <- c(diff(TPR), 0)
    out <- sum(TPR * dFPR) + sum(dTPR * dFPR) / 2
    return(out)
  }
  labels <- labels[order(scores, decreasing = TRUE)]
  df <- data.frame(
    False_Positive_Rate = cumsum(!labels) / sum(!labels),
    True_Positive_Rate =  cumsum(labels) / sum(labels)
  )
  df$AUC <- round(.calcAUC(df$True_Positive_Rate, df$False_Positive_Rate), 3)
  df$name <- name
  return(df)
}

.getQuantiles <- function(v = NULL, len = length(v)) {
  if (length(v) < len) {
    v2 <- rep(0, len)
    v2[seq_along(v)] <- v
  } else {
    v2 <- v
  }
  p <- trunc(rank(v2)) / length(v2)
  if (length(v) < len) {
    p <- p[seq_along(v)]
  }
  return(p)
}

.rowScale <- function(mat = NULL, min = NULL, max = NULL) {
  if (!is.null(min)) {
    rMin <- min
  } else {
    rMin <- matrixStats::rowMins(mat)
  }
  if (!is.null(max)) {
    rMax <- max
  } else {
    rMax <- matrixStats::rowMaxs(mat)
  }
  rScale <- rMax - rMin
  matDiff <- mat - rMin
  matScale <- matDiff / rScale
  out <- list(mat = matScale, min = rMin, max = rMax)
  return(out)
}

.quantileCut <- function(x = NULL, lo = 0.025, hi = 0.975, maxIf0 = TRUE) {
  q <- quantile(x, probs = c(lo, hi))
  if (q[2] == 0) {
    if (maxIf0) {
      q[2] <- max(x)
    }
  }
  x[x < q[1]] <- q[1]
  x[x > q[2]] <- q[2]
  return(x)
}

.normalizeCols <- function(mat = NULL, colSm = NULL, scaleTo = NULL) {

  if (is.null(colSm)) {
      colSm <- Matrix::colSums(mat)
  }
  if (!is.null(scaleTo)) {
    mat@x <- scaleTo * mat@x / rep.int(colSm, Matrix::diff(mat@p))
  } else {
    mat@x <- mat@x / rep.int(colSm, Matrix::diff(mat@p))
  }
  return(mat)
}

.safeSubset <- function(mat = NULL, subsetRows = NULL, subsetCols = NULL) {

  if (!is.null(subsetRows)) {
    idxNotIn <- which(subsetRows %ni% rownames(mat))
    if (length(idxNotIn) > 0) {
      subsetNamesNotIn <- subsetRows[idxNotIn]
      matNotIn <- Matrix::sparseMatrix(
        i = 1, j = 1, x = 0, dims = c(length(idxNotIn), ncol = ncol(mat))
      )
      rownames(matNotIn) <- subsetNamesNotIn
      mat <- rbind(mat, matNotIn)
    }
    mat <- mat[subsetRows, ]
  }

  if (!is.null(subsetCols)) {

    idxNotIn <- which(subsetCols %ni% colnames(mat))
    if (length(idxNotIn) > 0) {
      subsetNamesNotIn <- subsetCols[idxNotIn]
      matNotIn <- Matrix::sparseMatrix(
        i = 1, j = 1, x = 0, dims = c(nrow(mat), ncol = length(idxNotIn))
      )
      colnames(matNotIn) <- subsetNamesNotIn
      mat <- cbind(mat, matNotIn)
    }
    mat <- mat[, subsetCols]
  }
  mat
}

.groupMeans <- function(
  mat = NULL, groups = NULL, na.rm = TRUE, sparse = FALSE
) {
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == ncol(mat))
  gm <- lapply(unique(groups), function(x) {
    if (sparse) {
      Matrix::rowMeans(mat[, which(groups == x), drop = FALSE], na.rm = na.rm)
    } else {
      rowMeans(mat[, which(groups == x), drop = FALSE], na.rm = na.rm)
    }
  }) %>% Reduce("cbind", .)
  colnames(gm) <- unique(groups)
  return(gm)
}

.groupSums <- function(
  mat = NULL, groups = NULL, na.rm = TRUE, sparse = FALSE
) {
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == ncol(mat))
  gm <- lapply(unique(groups), function(x) {
    if (sparse) {
      Matrix::rowSums(mat[, which(groups == x), drop = FALSE], na.rm = na.rm)
    } else {
      rowSums(mat[, which(groups == x),drop=FALSE], na.rm = na.rm)
    }
  }) %>% Reduce("cbind", .)
  colnames(gm) <- unique(groups)
  return(gm)
}

.groupSds <- function(mat = NULL, groups = NULL, na.rm = TRUE, sparse = FALSE) {
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == ncol(mat))
  gs <- lapply(unique(groups), function(x) {
    if (sparse) {
      matrixStats::rowSds(
        as.matrix(mat[, which(groups == x), drop = FALSE]), na.rm = na.rm
      )
    } else {
      matrixStats::rowSds(
        mat[, which(groups == x), drop = FALSE], na.rm = na.rm
      )
    }
  }) %>% Reduce("cbind", .)
  colnames(gs) <- unique(groups)
  return(gs)
}

.centerRollMean <- function(v = NULL, k = NULL) {
  o1 <- data.table::frollmean(v, k, align = "right", na.rm = FALSE)
  if (k %% 2 == 0) {
    o2 <- c(
      rep(o1[k], floor(k / 2) - 1),
      o1[-seq_len(k - 1)],
      rep(o1[length(o1)], floor(k / 2))
    )
  }else if (k %% 2 == 1) {
    o2 <- c(
      rep(o1[k], floor(k / 2)),
      o1[-seq_len(k - 1)],
      rep(o1[length(o1)], floor(k / 2))
    )
  } else {
    stop("Error!")
  }
  o2
}

##############################################################################
# Miscellaneous Methods
###############################################################################

.splitEvery <- function(x = NULL, n = NULL) {
  #https://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
  if (is.atomic(x)) {
    split(x, ceiling(seq_along(x) / n))
  } else {
    split(x, ceiling(seq_len(nrow(x)) / n))
  }
}

.suppressAll <- function(expr = NULL) {
  suppressPackageStartupMessages(suppressMessages(suppressWarnings(expr)))
}

.getAssay <- function(se = NULL, assayName = NULL) {
  .assayNames <- function(se) {
    names(SummarizedExperiment::assays(se))
  }
  if (is.null(assayName)) {
    o <- SummarizedExperiment::assay(se)
  }else if (assayName %in% .assayNames(se)) {
    o <- SummarizedExperiment::assays(se)[[assayName]]
  } else {
    stop(
      sprintf(
        "assayName '%s' is not in assayNames of se : %s",
        assayName,
        paste(.assayNames(se), collapse = ", ")
      )
    )
  }
  return(o)
}

.fileExtension <- function(x = NULL) {
  pos <- regexpr("\\.([[:alnum:]]+)$", x)
  ifelse(pos > -1L, substring(x, pos + 1L), "")
}

.checkPath <- function(u = NULL, path = NULL, throwError = TRUE) {
  if (is.null(u)) {
    out <- TRUE
  }
  out <- lapply(u, function(x, error = TRUE) {
    if (Sys.which(x) == "") {
      if (!is.null(path) && file.exists(file.path(path, x))) {
        o <- TRUE
      } else {
        if (throwError) {
          stop(x, " not found in path, please add ", x, " to path!")
        } else {
          o <- FALSE
        }
      }
    } else {
      o <- TRUE
    }
    return(o)
  }) %>% unlist %>% all
  return(out)
}

.tempfile <- function(
  pattern = "tmp", tmpdir = "tmp", fileext = "", addDOC = TRUE
) {
  #if the directory doesnt already exist and file.exists evaluates to true,
  # then a file exists with that name

  if (!dir.exists(tmpdir)) {
    if (file.exists(tmpdir)) {
      stop(
        paste0(
          "Attempted to create temporary directory ",
          tmpdir,
          " but a file already exists with this name. Please remove this file \
          and try again!"
        )
      )
    }
  }

  dir.create(tmpdir, showWarnings = FALSE)

  if (!dir.exists(tmpdir)) {
    stop(
      paste0(
        "Unable to create temporary directory ",
        tmpdir,
        ". Check file permissions!"
      )
    )
  }

  if (addDOC) {
    doc <- paste0(
      "-Date-",
      Sys.Date(),
      "_Time-",
      gsub(
        ":",
        "-",
        stringr::str_split(Sys.time(),
        pattern = " ",
        simplify = TRUE)[1, 2]
      )
    )
  } else {
    doc <- ""
  }
  tempfile(
    pattern = paste0(pattern, "-"),
    tmpdir = tmpdir,
    fileext = paste0(doc, fileext)
  )
}

.ArchRLogo <- function(ascii = "Logo", messageLogo = TRUE) {
  Ascii <- list(
    Package = c("
           ___      .______        ______  __    __  .______      
          /   \\\     |   _  \\\      /      ||  |  |  | |   _  \\\     
         /  ^  \\\    |  |_)  |    |  ,----'|  |__|  | |  |_)  |    
        /  /_\\\  \\\   |      /     |  |     |   __   | |      /     
       /  _____  \\\  |  |\\\  \\\\___ |  `----.|  |  |  | |  |\\\  \\\\___.
      /__/     \\__\\ | _| `._____| \\______||__|  |__| | _| `._____|
    "),

    #modified from cyu@athena.mit.edu
    Logo = c("
                                                   / |
                                                 /    \\\
            .                                  /      |.
            \\\\\\                              /        |.
              \\\\\\                          /           `|.
                \\\\\\                      /              |.
                  \\\                    /                |\\\
                  \\\\#####\\\           /                  ||
                ==###########>      /                   ||
                 \\\\##==......\\\    /                     ||
            ______ =       =|__ /__                     ||      \\\\\\\
        ,--' ,----`-,__ ___/'  --,-`-===================##========>
       \\\               '        ##_______ _____ ,--,__,=##,__   ///
        ,    __==    ___,-,__,--'#'  ==='      `-'    | ##,-/
        -,____,---'       \\\\####\\\\________________,--\\\\_##,/
           ___      .______        ______  __    __  .______      
          /   \\\     |   _  \\\      /      ||  |  |  | |   _  \\\     
         /  ^  \\\    |  |_)  |    |  ,----'|  |__|  | |  |_)  |    
        /  /_\\\  \\\   |      /     |  |     |   __   | |      /     
       /  _____  \\\  |  |\\\  \\\\___ |  `----.|  |  |  | |  |\\\  \\\\___.
      /__/     \\__\\ | _| `._____| \\______||__|  |__| | _| `._____|
    ")
  )

  if(messageLogo) {
    message(Ascii[[ascii]])
  } else {
    Ascii[[ascii]]
  }

}

###############################################################################
# Batch Methods
###############################################################################

.safelapply <- function(..., threads = 1, preschedule = FALSE) {

  if (tolower(.Platform$OS.type) == "windows") {
    threads <- 1
  }

  if (threads > 1) {
    .requirePackage("parallel", source = "cran")
    o <- mclapply(..., mc.cores = threads, mc.preschedule = preschedule)

    errorMsg <- list()

    for (i in seq_along(o)) { #Make Sure this doesnt explode!
      if (inherits(o[[i]], "try-error")) {
        capOut <- utils::capture.output(o[[i]])
        capOut <- capOut[!grepl("attr\\(\\,|try-error", capOut)]
        capOut <- head(capOut, 10)
        capOut <- unlist(lapply(capOut, function(x) substr(x, 1, 250)))
        capOut <- paste0("\t", capOut)
        errorMsg[[length(errorMsg) + 1]] <- paste0(
          c(paste0("Error Found Iteration ", i, " : "), capOut), "\n"
        )
      }
    }
    if (length(errorMsg) != 0) {

      errorMsg <- unlist(errorMsg)
      errorMsg <- head(errorMsg, 50)
      errorMsg[1] <- paste0("\n", errorMsg[1])
      stop(errorMsg)

    }
  } else {
    o <- lapply(...)
  }
  o
}

.batchlapply <- function(args = NULL, sequential = FALSE) {

  if (is.null(args$tstart)) {
    args$tstart <- Sys.time()
  }

  #Determine Parallel Backend
  if (inherits(args$parallelParam, "BatchtoolsParam")) {

    .logStop(
      "Batchtools not yet fully supported please use local parallel threading!",
      logFile = args$logFile
    )

    .logDiffTime(
      "Batch Execution w/ BatchTools through BiocParallel!",
      t1 = args$tstart,
      verbose = TRUE,
      logFile = args$logFile
    )

    require(BiocParallel)

    args$parallelParam <- btParam
    #Unlink registry Directory
    if (dir.exists(args$registryDir)) {
      #Clean Up Registry
      unlink(args$registryDir, recursive = TRUE)# Delete registry directory
    }

    #Set Up Registry For Runnning
    args$parallelParam$registryargs <- batchtoolsRegistryargs(
      file.dir = args$registryDir,
      work.dir = getwd(),
      packages = character(0L),
      namespaces = character(0L),
      source = character(0L),
      load = character(0L)
    )

    #Register
    BPPARAM <- args$parallelParam
    register(BPPARAM)

    #Add To Args
    args$BPPARAM <- BPPARAM

    if ("..." %in% names(args)) {
      args["..."] <- NULL
    }

    #Run
    args <- args[names(args) %ni% c("threads", "parallelParam", "subThreading")]
    outlist <- do.call(bplapply, args)

  } else {

    .logDiffTime(
      "Batch Execution w/ safelapply!",
      t1 = args$tstart,
      verbose = TRUE,
      logFile = args$logFile
    )
    if (sequential) {
      args$subThreads <- args$threads
      args$threads <- 1
    } else {
      if (args$threads > length(args$X)) {
        args$subThreads <- floor(args$threads / length(args$X))
        args$threads <- length(args$X)
      } else {
        args$subThreads <- 1
      }
    }

    args <- args[
      names(args) %ni% c("registryDir", "parallelParam", "subThreading")
    ]
    outlist <- do.call(.safelapply, args)

  }

  return(outlist)

}

.retryCatch <- function(
  expr,
  ...,
  maxAttempts = 3,
  warnAttempts = FALSE,
  nameFN = "FN",
  printInfo = NULL,
  logFile = NULL
) {
  currentAttempt <- 0
  completed <- FALSE
  while (!completed & currentAttempt <= maxAttempts) {
    currentAttempt <- currentAttempt + 1
    if (currentAttempt > 1) {
      .logMessage(
        nameFN,
        " : Error occured, attempting again (",
        currentAttempt - 1,
        " of ",
        maxAttempts,
        ")",
        logFile = logFile
      )
    }
    ###########################################################
    tryResult <- tryCatch({
      #########################################################
      #Try Catch Statement Here
      if (warnAttempts) {
        out <- return(expr)
      } else {
        out <- suppressWarnings(return(expr))
      }
      #########################################################
      list(out = out, completed = TRUE)
    }, error = function(e) {
      list(out = e, completed = FALSE)
    }, ...)
    ###########################################################
    completed <- tryResult$completed
  }
  if (!completed) {
    .logMessage(
      nameFN,
      " : Error occured and could not be resolved after ",
      maxAttempts,
      " additional attempts!",
      logFile = logFile
    )
    if (!is.null(printInfo)) {
      .logMessage("Error occured at ", printInfo, logFile = logFile)
    }
    print(tryResult[[1]])
    stop()
  }

  tryResult[[1]]

}

###############################################################################
# Developer Utils
###############################################################################

.devMode <- function(package = "ArchR") {
  fn <- unclass(lsf.str(envir = asNamespace("ArchR"), all = TRUE))
  for (i in seq_along(fn)) {
    tryCatch({
      eval(parse(text = paste0(fn[i], "<-ArchR:::", fn[i])))
    }, error = function(x) {
    })
  }
}

.convertToPNG <- function(
  ArchRProj = NULL,
  paths = c("QualityControl"),
  recursive = TRUE,
  outDir = "Figures",
  command = "mv"
) {

  .requirePackage("pdftools", source = "cran")

  if (!is.null(ArchRProj)) {
    paths <- c(paths, file.path(getOutputDirectory(ArchRProj), "Plots"))
  }

  pdfFiles <- lapply(seq_along(paths), function(i) {
    if (recursive) {
      dirs <- list.dirs(paths[i], recursive = FALSE, full.names = FALSE)
      if (length(dirs) > 0) {
        pdfs <- lapply(seq_along(dirs), function(j) {
          list.files(
            file.path(paths[i], dirs[j]), full.names = TRUE, pattern = "\\.pdf"
          )
        }) %>% unlist
      } else {
        pdfs <- c()
      }
      pdfs <- c(
        list.files(paths[i], full.names = TRUE, pattern = "\\.pdf"), pdfs
      )
    } else {
      pdfs <- list.files(paths[i], full.names = TRUE, pattern = "\\.pdf")
    }
    pdfs
  }) %>% unlist

  dir.create(outDir, showWarnings = FALSE)

  for (i in seq_along(pdfFiles)) {
    print(i)
    tryCatch({
      pdf_convert(
        pdfFiles[i],
        format = "png",
        pages = NULL,
        filenames = file.path(
          outDir, gsub("\\.pdf", "_%d.png", basename(pdfFiles[i]))
        ),
        dpi = 300,
        opw = "",
        upw = "",
        verbose = TRUE
      )
      system(
        paste0(
          command,
          " ",
          pdfFiles[i],
          " ",
          file.path(outDir, basename(pdfFiles[i]))
        )
      )
    },
    error = function(x) {
      0
    })
  }
}
