args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
    library(mutscan)
    library(SummarizedExperiment)
})

print(dataset)
print(inrds)
print(outrds)

mutscan_config <- yaml::read_yaml(file.path(dataset, "mutscan_config.yml"))
(replicateCol <- mutscan_config$replicateCol)
(conditionCol <- mutscan_config$conditionCol)
(numerator <- mutscan_config$numerator)
(denominator <- mutscan_config$denominator)
(odCol <- mutscan_config$odCol)
wtrows <- mutscan_config$wtrows
(normMethod <- mutscan_config$normMethod)

sec <- readRDS(inrds)

(wtrows <- intersect(wtrows, rownames(sec)))

sec$replicate <- sapply(strsplit(sec[[replicateCol]], ","), .subset, 1)
sec$condition <- sapply(strsplit(sec[[conditionCol]], ","), .subset, 1)
if (!(odCol %in% colnames(colData(sec)))) {
    sec[[odCol]] <- 1
}
if (is.character(sec[[odCol]])) {
    sec[[odCol]] <- sapply(strsplit(sec[[odCol]], ","), .subset, 1)
    sec[[odCol]] <- as.numeric(sec[[odCol]])
}

## Keep only features where all input samples have at least 50 reads
sec <- sec[rowSums(assay(sec[, sec$condition == denominator], "counts") > 50) == 
               length(which(sec$condition == denominator))]
dim(sec)

## Generate all subsets of ceiling(N/2) replicates. In each step, 
## use one such subset as the artificial "group1" and all other replicates 
## as the artificial "group2". Fit a model to test whether there are 
## differences between the log fold changes in the two groups. 
(subsets <- combn(unique(sec$replicate), 
                  ceiling(length(unique(sec$replicate))/2), 
                  simplify = FALSE))

pvals <- do.call(dplyr::bind_rows, lapply(subsets, function(s) {
    sec$group <- ifelse(sec$replicate %in% s, "group1", "group2")
    design <- model.matrix(~ replicate, data = colData(sec))
    design <- cbind(design, 
                    sel1 = as.numeric(sec$group == "group1" & 
                                          sec$condition == numerator),
                    sel2 = as.numeric(sec$group == "group2" & 
                                          sec$condition == numerator))
    contrast <- (colnames(design) == "sel2") - (colnames(design) == "sel1")
    resfc_edgeR <- calculateRelativeFC(
        sec, 
        design = design,
        contrast = contrast,
        WTrows = wtrows, normMethod = normMethod,
        method = "edgeR")
    resfc_limma <- calculateRelativeFC(
        sec, 
        design = design,
        contrast = contrast,
        WTrows = wtrows, normMethod = normMethod,
        method = "limma")
    
    dplyr::bind_rows(
        data.frame(g1 = paste(s, collapse = "_"),
                   method = "edgeR",
                   pvalue = resfc_edgeR$PValue),
        data.frame(g1 = paste(s, collapse = "_"),
                   method = "limma",
                   pvalue = resfc_limma$P.Value)
    )
}))

saveRDS(pvals, file = outrds)

date()
sessionInfo()
