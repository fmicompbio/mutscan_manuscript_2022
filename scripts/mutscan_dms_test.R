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
resppi <- calculateFitnessScore(sec, pairingCol = "replicate",
                                ODCols = odCol, WTrows = wtrows,
                                comparison = c("condition", numerator, denominator))
design <- model.matrix(~ replicate + condition, data = colData(sec))
resfc_edgeR <- calculateRelativeFC(
    sec, 
    design = design,
    contrast = (colnames(design) == paste0("condition", numerator)) - 
        (colnames(design) == paste0("condition", denominator)),
    WTrows = wtrows, normMethod = normMethod,
    method = "edgeR")
resfc_limma <- calculateRelativeFC(
    sec, 
    design = design,
    contrast = (colnames(design) == paste0("condition", numerator)) - 
        (colnames(design) == paste0("condition", denominator)),
    WTrows = wtrows, normMethod = normMethod,
    method = "limma")

saveRDS(list(ppi = resppi, fc_edgeR = resfc_edgeR, fc_limma = resfc_limma),
        file = outrds)

date()
sessionInfo()
