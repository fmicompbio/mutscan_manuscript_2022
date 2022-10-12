args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
    library(mutscan)
    library(yaml)
    library(dplyr)
    library(Matrix)
    library(sparseMatrixStats)
    library(SummarizedExperiment)
})

print(dataset)
print(outrds)

samplefile <- read.delim(file.path(dataset, "samples.txt"), header = TRUE)
digests <- file.path(dataset, "output", "mutscan", "digestFastqs", 
                     paste0(unique(samplefile$mutscan_name), ".rds"))
names(digests) <- sub("\\.rds$", "", basename(digests))
print(digests)
stopifnot(all(file.exists(digests)))
digests <- lapply(digests, readRDS)

mutscan_config <- yaml::read_yaml(file.path(dataset, "mutscan_config.yml"))
print(mutscan_config)

samplefile <- samplefile %>%
    dplyr::rename(Name = mutscan_name) %>%
    dplyr::select(-any_of(c("pair1", "pair2"))) %>%
    dplyr::group_by(Name) %>%
    dplyr::summarize(across(everything(), paste, collapse = ",")) %>%
    dplyr::distinct()
samplefile <- samplefile[match(names(digests), samplefile$Name), ]

res <- mutscan::summarizeExperiment(digests, coldata = samplefile, 
                                    countType = mutscan_config$countType)

if (mutscan_config$collapseToAA) {
    ## For consistency with DiMSum (for comparison purposes), don't allow 
    ## sequences with both synonymous and non-synonymous mutations
    res <- res[!(grepl(",", rowData(res)$mutationTypes) & 
                     grepl("silent", rowData(res)$mutationTypes)), ]
    res <- collapseMutantsByAA(res)
}

dim(res)

## Filter
colData(res)$selsub <- sapply(strsplit(res$selection_id, ","), .subset, 1)
countInput <- as(assay(res, "counts")[, res$selsub == "0"], "dgCMatrix")
countOutput <- as(assay(res, "counts")[, res$selsub == "1"], "dgCMatrix")
resFilt <- res[rowMins(countInput) >= mutscan_config$minInputCountAll & 
                   rowMaxs(countInput) >= mutscan_config$minInputCountAny & 
                   rowMins(countOutput) >= mutscan_config$minOutputCountAll & 
                   rowMaxs(countOutput) >= mutscan_config$minOutputCountAny, ]
dim(resFilt)

mutscan::generateQCReport(resFilt, outFile = sub("_se.rds", "_qc.html", outrds),
                          reportTitle = dataset, forceOverwrite = TRUE)

saveRDS(res, file = outrds)
saveRDS(resFilt, file = sub("_se.rds", "_se_filt.rds", outrds))

date()
sessionInfo()
