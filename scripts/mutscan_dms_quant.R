args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
    library(mutscan)
    library(yaml)
})

print(dataset)
print(sample)
print(nthreads)
print(nreads)
print(outrds)

samplefile <- read.delim(file.path(dataset, "samples.txt"), header = TRUE)
idx <- which(samplefile$mutscan_name == sample)
fastqforward <- file.path(dataset, "FASTQ", samplefile$pair1[idx])
fastqreverse <- file.path(dataset, "FASTQ", samplefile$pair2[idx])
if (length(fastqreverse) == 0) {
    fastqreverse <- NULL
}
dataset_config <- yaml::read_yaml(file.path(dataset, "dataset_config.yml"))
mutscan_config <- yaml::read_yaml(file.path(dataset, "mutscan_config.yml"))

print(fastqforward)
print(fastqreverse)
print(dataset_config)
print(mutscan_config)

print(unlist(dataset_config$wildTypeForward))
print(unlist(dataset_config$wildTypeReverse))
print(any(grepl(mutscan_config$mutNameDelimiter, c(names(unlist(dataset_config$wildTypeForward)), 
                              names(unlist(dataset_config$wildTypeReverse))), fixed = TRUE)))
print(names(unlist(dataset_config$wildTypeForward)))
print(names(unlist(dataset_config$wildTypeReverse)))

res <- mutscan::digestFastqs(
    fastqForward = fastqforward,
    fastqReverse = fastqreverse,
    mergeForwardReverse = mutscan_config$mergeForwardReverse,
    minOverlap = mutscan_config$minOverlap,
    maxOverlap = mutscan_config$maxOverlap,
    minMergedLength = mutscan_config$minMergedLength,
    maxMergedLength = mutscan_config$maxMergedLength,
    maxFracMismatchOverlap = mutscan_config$maxFracMismatchOverlap,
    greedyOverlap = mutscan_config$greedyOverlap,
    revComplForward = mutscan_config$revComplForward,
    revComplReverse = mutscan_config$revComplReverse,
    adapterForward = dataset_config$adapterForward,
    adapterReverse = dataset_config$adapterReverse,
    elementsForward = dataset_config$elementsForward,
    elementLengthsForward = dataset_config$elementLengthsForward,
    elementsReverse = dataset_config$elementsReverse,
    elementLengthsReverse = dataset_config$elementLengthsReverse,
    primerForward = dataset_config$primerForward,
    primerReverse = dataset_config$primerReverse,
    wildTypeForward = unlist(dataset_config$wildTypeForward),
    wildTypeReverse = unlist(dataset_config$wildTypeReverse),
    constantForward = dataset_config$constantForward,
    constantReverse = dataset_config$constantReverse,
    avePhredMinForward = mutscan_config$avePhredMinForward,
    avePhredMinReverse = mutscan_config$avePhredMinReverse,
    variableNMaxForward = mutscan_config$variableNMaxForward,
    variableNMaxReverse = mutscan_config$variableNMaxReverse,
    umiNMax = mutscan_config$umiNMax,
    nbrMutatedCodonsMaxForward = mutscan_config$nbrMutatedCodonsMaxForward,
    nbrMutatedCodonsMaxReverse = mutscan_config$nbrMutatedCodonsMaxReverse,
    nbrMutatedBasesMaxForward = mutscan_config$nbrMutatedBasesMaxForward,
    nbrMutatedBasesMaxReverse = mutscan_config$nbrMutatedBasesMaxReverse,
    forbiddenMutatedCodonsForward = mutscan_config$forbiddenMutatedCodonsForward,
    forbiddenMutatedCodonsReverse = mutscan_config$forbiddenMutatedCodonsReverse,
    useTreeWTmatch = mutscan_config$useTreeWTmatch,
    mutatedPhredMinForward = mutscan_config$mutatedPhredMinForward,
    mutatedPhredMinReverse = mutscan_config$mutatedPhredMinReverse,
    mutNameDelimiter = mutscan_config$mutNameDelimiter,
    constantMaxDistForward = mutscan_config$constantMaxDistForward,
    constantMaxDistReverse = mutscan_config$constantMaxDistReverse,
    variableCollapseMaxDist = mutscan_config$variableCollapseMaxDist,
    variableCollapseMinReads = mutscan_config$variableCollapseMinReads,
    variableCollapseMinRatio = mutscan_config$variableCollapseMinRatio,
    umiCollapseMaxDist = mutscan_config$umiCollapseMaxDist,
    filteredReadsFastqForward = mutscan_config$filteredReadsFastqForward,
    filteredReadsFastqReverse = mutscan_config$filteredReadsFastqReverse,
    maxNReads = nreads,
    verbose = mutscan_config$verbose,
    nThreads = nthreads,
    chunkSize = mutscan_config$chunkSize
)

saveRDS(res, file = outrds)

date()
sessionInfo()
