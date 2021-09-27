estimateCellCounts2Mod <- function (rgSet,
                                 compositeCellType = "Blood",
                                 processMethod = "preprocessNoob",
                                 probeSelect = c("auto", "any", "IDOL"),
                                 cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"),
                                 referencePlatform = c("IlluminaHumanMethylation450k",
                                                       "IlluminaHumanMethylationEPIC",
                                                       "IlluminaHumanMethylation27k"),
                                 referenceset = NULL,
                                 IDOLOptimizedCpGs = NULL,
                                 returnAll = FALSE,
                                 meanPlot = FALSE,
                                 verbose = TRUE,
                                 keepProbes = NULL,
                                 ...)
{
  if ((!is(rgSet, "RGChannelSet")) && (!is(rgSet, "MethylSet")))
    stop(strwrap(sprintf("object is of class '%s', but needs to be of \n                                class 'RGChannelSet' 'RGChannelSetExtended' or \n                                'MethylSet' to use this function",
                         class(rgSet)), width = 80, prefix = " ", initial = ""))
  if (!is(rgSet, "RGChannelSet") && (processMethod[1] != "preprocessQuantile"))
    stop(strwrap(sprintf("object is of class '%s', but needs to be of \n                                class 'RGChannelSet' or 'RGChannelSetExtended' \n                                to use other methods different to \n                                'preprocessQuantile'",
                         class(rgSet)), width = 80, prefix = " ", initial = ""))
  if (is(rgSet, "MethylSet") && (processMethod[1] == "preprocessQuantile"))
    message(strwrap("[estimateCellCounts2] The function will assume that\n                            no preprocessing has been performed. Using \n                            'preprocessQuantile' in prenormalized data is \n                            experimental and it should only be run under the \n                            user responsibility",
                    width = 80, prefix = " ", initial = ""))
  if (is(rgSet, "RGChannelSetExtended"))
    rgSet <- as(rgSet, "RGChannelSet")
  referencePlatform <- match.arg(referencePlatform)
  rgPlatform <- sub("IlluminaHumanMethylation", "", annotation(rgSet)[which(names(annotation(rgSet)) ==
                                                                              "array")])
  platform <- sub("IlluminaHumanMethylation", "", referencePlatform)
  if ((compositeCellType == "CordBlood") && (!"nRBC" %in% cellTypes))
    message(strwrap("[estimateCellCounts2] Consider including 'nRBC' in \n                        argument 'cellTypes' for cord blood estimation.\n",
                    width = 80, prefix = " ", initial = ""))
  if ((compositeCellType == "Blood") && (referencePlatform ==
                                         "IlluminaHumanMethylationEPIC") && ("Gran" %in% cellTypes))
    message(strwrap("[estimateCellCounts2] Replace 'Gran' for 'Neu' in \n                        argument 'cellTypes' for EPIC blood estimation.\n",
                    width = 80, prefix = " ", initial = ""))
  referencePkg <- sprintf("FlowSorted.%s.%s", compositeCellType,
                          platform)
  subverbose <- max(as.integer(verbose) - 1L, 0L)
  if (!is.null(referenceset)) {
    referenceRGset <- get(referenceset)
    if (!is(rgSet, "RGChannelSet"))
      referenceRGset <- preprocessRaw(referenceRGset)
  }
  else {
    if (!require(referencePkg, character.only = TRUE))
      stop(strwrap(sprintf("Could not find reference data package for \n                                compositeCellType '%s' and referencePlatform \n                                '%s' (inferred package name is '%s')",
                           compositeCellType, platform, referencePkg), width = 80,
                   prefix = " ", initial = ""))
    if (referencePkg != "FlowSorted.Blood.EPIC") {
      referenceRGset <- get(referencePkg)
    }
    else {
      hub <- ExperimentHub()
      referenceRGset <- hub[["EH1136"]]
    }
    if (!is(rgSet, "RGChannelSet"))
      referenceRGset <- preprocessRaw(referenceRGset)
  }
  if (rgPlatform != platform) {
    rgSet <- convertArray(rgSet, outType = referencePlatform,
                          verbose = TRUE)
  }
  if (!"CellType" %in% names(colData(referenceRGset)))
    stop(strwrap(sprintf("the reference sorted dataset (in this case '%s') \n                            needs to have a phenoData column called \n                            'CellType'"),
                 names(referencePkg), width = 80, prefix = " ", initial = ""))
  if (sum(colnames(rgSet) %in% colnames(referenceRGset)) >
      0)
    stop(strwrap("the sample/column names in the user set must not be in \n                    the reference data ",
                 width = 80, prefix = " ", initial = ""))
  if (!all(cellTypes %in% referenceRGset$CellType))
    stop(strwrap(sprintf("all elements of argument 'cellTypes' needs to be \n                            part of the reference phenoData columns 'CellType' \n                            (containg the following elements: '%s')",
                         paste(unique(referenceRGset$cellType), collapse = "', '")),
                 width = 80, prefix = " ", initial = ""))
  if (length(unique(cellTypes)) < 2)
    stop("At least 2 cell types must be provided.")
  if ((processMethod == "auto") && (compositeCellType %in%
                                    c("Blood", "DLPFC")))
    processMethod <- "preprocessQuantile"
  if ((processMethod == "auto") && (!compositeCellType %in%
                                    c("Blood", "DLPFC")) && (is(rgSet, "RGChannelSet")))
    processMethod <- "preprocessNoob"
  processMethod <- get(processMethod)
  if ((probeSelect == "auto") && (compositeCellType == "CordBlood")) {
    probeSelect <- "any"
  }
  if ((probeSelect == "auto") && (compositeCellType != "CordBlood")) {
    probeSelect <- "both"
  }
  if (verbose)
    message(strwrap("[estimateCellCounts2] Combining user data with \n                        reference (flow sorted) data.\n",
                    width = 80, prefix = " ", initial = ""))
  newpd <- DataFrame(sampleNames = c(colnames(rgSet), colnames(referenceRGset)),
                     studyIndex = rep(c("user", "reference"), times = c(ncol(rgSet),
                                                                        ncol(referenceRGset))), stringsAsFactors = FALSE)
  if (is.null(rgSet$CellType))
    rgSet$CellType <- rep("NA", dim(rgSet)[2])
  commoncolumn <- intersect(names(colData(rgSet)), names(colData(referenceRGset)))
  colData(referenceRGset)[commoncolumn] <- mapply(FUN = as,
                                                  colData(referenceRGset)[commoncolumn], vapply(colData(rgSet)[commoncolumn],
                                                                                                class, FUN.VALUE = character(1)), SIMPLIFY = FALSE)
  colData(referenceRGset) <- colData(referenceRGset)[commoncolumn]
  colData(rgSet) <- colData(rgSet)[commoncolumn]
  referencePd <- colData(referenceRGset)
  combinedRGset <- combineArrays(rgSet, referenceRGset, outType = referencePlatform)
  colData(combinedRGset) <- newpd
  colnames(combinedRGset) <- newpd$sampleNames
  rm(referenceRGset)
  if (verbose)
    message(strwrap("[estimateCellCounts2] Processing user and reference \n                        data together.\n",
                    width = 80, prefix = " ", initial = ""))
  if (compositeCellType == "CordBlood") {
    if (!is(combinedRGset, "RGChannelSet"))
      combinedRGset@preprocessMethod["rg.norm"] <- "Raw (no normalization or bg correction)"
    combinedMset <- processMethod(combinedRGset, verbose = subverbose)
    rm(combinedRGset)
    gc()
    compTable <- get(paste0(referencePkg, ".compTable"))
    combinedMset <- combinedMset[which(rownames(combinedMset) %in%
                                         rownames(compTable)), ]
  }
  else {
    if (!is(combinedRGset, "RGChannelSet"))
      combinedRGset@preprocessMethod["rg.norm"] <- "Raw (no normalization or bg correction)"
    combinedMset <- processMethod(combinedRGset)
    rm(combinedRGset)
    gc()
  }
  referenceMset <- combinedMset[, combinedMset$studyIndex ==
                                  "reference"]
  if(!is.null(keepProbes)){
    referenceMset <- referenceMset[featureNames(referenceMset) %in% keepProbes,]
  }
  colData(referenceMset) <- as(referencePd, "DataFrame")
  mSet <- combinedMset[, combinedMset$studyIndex == "user"]
  colData(mSet) <- as(colData(rgSet), "DataFrame")
  rm(combinedMset)
  if (probeSelect != "IDOL") {
    if (verbose)
      message(strwrap("[estimateCellCounts2] Picking probes for \n                            composition estimation.\n",
                      width = 80, prefix = " ", initial = ""))
    compData <- minfi:::pickCompProbes(referenceMset, cellTypes = cellTypes,
                                       compositeCellType = compositeCellType, probeSelect = probeSelect)
    coefs <- compData$coefEsts
    if (verbose)
      message("[estimateCellCounts2] Estimating composition.\n")
    counts <- minfi:::projectCellType(getBeta(mSet)[rownames(coefs),
                                                    ], coefs)
    rownames(counts) <- colnames(rgSet)
    if (meanPlot) {
      smeans <- compData$sampleMeans
      smeans <- smeans[order(names(smeans))]
      sampleMeans <- c(colMeans(minfi::getBeta(mSet)[rownames(coefs),
                                                     ]), smeans)
      sampleColors <- c(rep(1, ncol(mSet)), 1 + as.numeric(factor(names(smeans))))
      plot(sampleMeans, pch = 21, bg = sampleColors)
      legend("bottomleft", c("blood", levels(factor(names(smeans)))),
             col = seq_len(7), pch = 15)
    }
    if (returnAll) {
      list(counts = counts, compTable = compData$compTable,
           normalizedData = mSet)
    }
    else {
      list(counts = counts)
    }
  }
  else {
    if (verbose)
      message(strwrap("[estimateCellCounts2] Using IDOL L-DMR probes for \n                            composition estimation.\n",
                      width = 80, prefix = " ", initial = ""))
    p <- getBeta(referenceMset)
    pd <- as.data.frame(colData(referenceMset))
    rm(referenceMset)
    if (!is.null(cellTypes)) {
      if (!all(cellTypes %in% pd$CellType))
        stop(strwrap("elements of argument 'cellTypes' is not part of \n                            'referenceMset$CellType'",
                     width = 80, prefix = " ", initial = ""))
      keep <- which(pd$CellType %in% cellTypes)
      pd <- pd[keep, ]
      p <- p[, keep]
    }
    pd$CellType <- factor(pd$CellType, levels = cellTypes)
    ffComp <- rowFtests(p, pd$CellType)
    tIndexes <- split(seq(along = pd$CellType), pd$CellType)

    prof <- vapply(tIndexes, function(i) rowMeans(p[, i]),
                   FUN.VALUE = numeric(dim(p)[1]))
    r <- rowRanges(p)
    compTable <- cbind(ffComp, prof, r, abs(r[, 1] - r[,
                                                       2]))
    names(compTable)[1] <- "Fstat"
    names(compTable)[c(-2, -1, 0) + ncol(compTable)] <- c("low",
                                                          "high", "range")
    tstatList <- lapply(tIndexes, function(i) {
      x <- rep(0, ncol(p))
      x[i] <- 1
      return(rowttests(p, factor(x)))
    })
    trainingProbes <- IDOLOptimizedCpGs
    trainingProbes <- trainingProbes[trainingProbes %in%
                                       rownames(p)]
    p <- p[trainingProbes, ]
    pMeans <- colMeans(p)
    names(pMeans) <- pd$CellType
    form <- as.formula(sprintf("y ~ %s - 1", paste(levels(pd$CellType),
                                                   collapse = "+")))
    phenoDF <- as.data.frame(model.matrix(~pd$CellType -
                                            1))
    colnames(phenoDF) <- sub("^pd\\$CellType", "", colnames(phenoDF))
    if (ncol(phenoDF) == 2) {
      X <- as.matrix(phenoDF)
      coefEsts <- t(solve(t(X) %*% X) %*% t(X) %*% t(p))
      coefs <- coefEsts
    }
    else {
      tmp <- validationCellType(Y = p, pheno = phenoDF,
                                modelFix = form)
      coefEsts <- tmp$coefEsts
      coefs <- coefEsts
    }
    compData <- list(coefEsts = coefEsts, compTable = compTable,
                     sampleMeans = pMeans)
    if (verbose)
      message("[estimateCellCounts2] Estimating composition.\n")
    counts <- projectCellType(getBeta(mSet)[rownames(coefs),
                                            ], coefs)
    rownames(counts) <- colnames(rgSet)
    if (meanPlot) {
      smeans <- compData$sampleMeans
      smeans <- smeans[order(names(smeans))]
      sampleMeans <- c(colMeans(getBeta(mSet)[rownames(coefs),
                                              ]), smeans)
      sampleColors <- c(rep(1, ncol(mSet)), 1 + as.numeric(factor(names(smeans))))
      plot(sampleMeans, pch = 21, bg = sampleColors)
      legend("bottomleft", c("blood", levels(factor(names(smeans)))),
             col = seq_len(7), pch = 15)
    }
    if (returnAll) {
      list(counts = counts, compTable = compTable, normalizedData = mSet)
    }
    else {
      list(counts = counts)
    }
  }
}


