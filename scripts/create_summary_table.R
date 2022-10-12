args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(ComplexUpset)
    library(tibble)
    library(mutscan)
    library(SummarizedExperiment)
})
Rcpp::sourceCpp("scripts/seqToName.cpp")
source("scripts/enrich2_to_mutscan_ids.R")

print(dataset)
print(outrds)

## -----------------------------------------------------------------------------
## Define input files
## -----------------------------------------------------------------------------
## mutscan
mutscanCountFile <- file.path(dataset, "output", "mutscan", 
                              paste0(dataset, "_mutscan_se.rds"))
mutscanResFile <- file.path(dataset, "output", "mutscan", 
                            paste0(dataset, "_mutscan_testresults.rds"))

## DiMSum
dimsumCountFile <- file.path(dataset, "output", "dimsum", 
                             "dimsum_variant_data_merge.tsv")
dimsumResFile <- file.path(dataset, "output", "dimsum", 
                           "dimsum_fitness_replicates.RData")

## Enrich2
enrich2CountFile <- file.path(dataset, "output", "enrich2", "tsv", 
                              paste0(dataset, "_exp"), "main_variants_counts.tsv")
enrich2ResFile <- file.path(dataset, "output", "enrich2", "tsv", 
                            paste0(dataset, "_exp"), "main_variants_scores_shared_full.tsv")

if (!file.exists(enrich2CountFile) | !file.exists(enrich2ResFile)) {
    include_enrich2 <- FALSE
} else {
    include_enrich2 <- TRUE
}

## -----------------------------------------------------------------------------
## Set some parameters
## -----------------------------------------------------------------------------
mutscan_config <- yaml::read_yaml(file.path(dataset, "mutscan_config.yml"))
if (mutscan_config$collapseToAA) {
    nmutColDimsum <- "Nham_aa"
    nmutColMutscan <- "minNbrMutAAs"
    mutType <- "amino acid"
    seqColDimsum <- "aa_seq"
} else {
    nmutColDimsum <- "Nham_nt"
    nmutColMutscan <- "minNbrMutBases"
    mutType <- "base"
    seqColDimsum <- "nt_seq"
}

## -----------------------------------------------------------------------------
## Read mutscan data
## -----------------------------------------------------------------------------
## Read mutscan counts (SE object)
mutscan_se <- readRDS(mutscanCountFile)
mutscan_se

## Read mutscan fitness values
mutscan_fitness <- readRDS(mutscanResFile)
head(mutscan_fitness$fc_edgeR, 4)

mutscan <- cbind(variant = rownames(mutscan_se), 
                 as.data.frame(as.matrix(assay(mutscan_se, "counts"))) %>%
                     setNames(paste0("mutscan_", colnames(mutscan_se), "_counts"))) %>%
    dplyr::full_join(data.frame(variant = rownames(mutscan_se), 
                                mutscan_nbr_mutations = rowData(mutscan_se)[[nmutColMutscan]], 
                                mutscan_detected = TRUE), 
                     by = "variant") %>%
    dplyr::full_join(data.frame(variant = rownames(mutscan_fitness$fc_edgeR), 
                                mutscan_edgeR_fitness = mutscan_fitness$fc_edgeR$logFC_shrunk), 
                     by = "variant") %>%
    dplyr::full_join(data.frame(variant = rownames(mutscan_fitness$fc_limma), 
                                mutscan_limma_fitness = mutscan_fitness$fc_limma$logFC), 
                     by = "variant")
head(mutscan, 4)

## -----------------------------------------------------------------------------
## Read DiMSum data
## -----------------------------------------------------------------------------
## Read DiMSum counts and fix column names
dimsum_counts <- read.delim(dimsumCountFile)
(idx <- which(grepl("_e.*_s.*_b.*_count", colnames(dimsum_counts))))
colnames(dimsum_counts)[idx] <- paste0(vapply(strsplit(colnames(dimsum_counts)[idx], "_"), 
                                              .subset, 1, FUN.VALUE = ""), "_counts")
if (dataset == "Bolognesi_TDP43_290_331") {
    for (cn in colnames(mutscan_se)) {
        g <- grep(cn, colnames(dimsum_counts), value = TRUE)
        if (length(g) > 1) {
            dimsum_counts[[paste0(cn, "_counts")]] <- rowSums(dimsum_counts[, g, drop = FALSE])
            dimsum_counts <- dimsum_counts[, !(colnames(dimsum_counts) %in% g)]
        }
    }
}
head(dimsum_counts, 4)
dim(dimsum_counts)
length(unique(dimsum_counts[[seqColDimsum]]))

## If summarizing on the aa level, aggregate the counts
if (mutscan_config$collapseToAA) {
    dimsum_counts <- dimsum_counts %>%
        dplyr::select(aa_seq, Nham_aa, paste0(colnames(mutscan_se), "_counts")) %>%
        dplyr::group_by(aa_seq, Nham_aa) %>%
        dplyr::summarize(across(everything(), sum), .groups = "drop")
}

## Read DiMSum fitness values and fix column names
load(dimsumResFile)
dimsum_fitness <- all_variants
samples <- read.delim(file.path(dataset, "samples.txt"), 
                      header = TRUE) %>%
    dplyr::mutate(matchName = paste0("count_e", experiment, "_s", selection_id)) %>%
    dplyr::select(matchName, sample_name)
samples
(idx <- match(samples$matchName, colnames(dimsum_fitness)))
colnames(dimsum_fitness)[idx] <- paste0(samples$sample_name, "_counts")
head(dimsum_fitness, 4)
dim(dimsum_fitness)
length(unique(dimsum_fitness[[seqColDimsum]]))

## Add variant name column
if (mutscan_config$collapseToAA) {
    wtseq <- rowData(mutscan_se)$sequenceAA[rowData(mutscan_se)$nbrMutAAs == "0"]
} else {
    wtseq <- rowData(mutscan_se)$sequence[rowData(mutscan_se)$nbrMutBases == "0"]
}

if (grepl("_", wtseq)) {
    wtseqf <- sapply(strsplit(wtseq, "_"), .subset, 1)
    wtseqr <- sapply(strsplit(wtseq, "_"), .subset, 2)
    ## counts
    dimsum_counts$mutname <- paste0(
        vapply(substr(toupper(dimsum_counts[[seqColDimsum]]), 1, 
                      nchar(wtseqf)), function(w) {
            seqToName(w, wtseqf, "base", "f")
        }, ""), 
        "_", 
        vapply(substr(toupper(dimsum_counts[[seqColDimsum]]), nchar(wtseqf) + 1, 
                      nchar(wtseqf) + nchar(wtseqr)), function(w) {
            seqToName(w, wtseqr, "base", "r")
        }, "")
    )
    ## fitness
    dimsum_fitness$mutname <- paste0(
        vapply(substr(toupper(dimsum_fitness[[seqColDimsum]]), 1, 
                      nchar(wtseqf)), function(w) {
                          seqToName(w, wtseqf, "base", "f")
                      }, ""), 
        "_", 
        vapply(substr(toupper(dimsum_fitness[[seqColDimsum]]), nchar(wtseqf) + 1, 
                      nchar(wtseqf) + nchar(wtseqr)), function(w) {
                          seqToName(w, wtseqr, "base", "r")
                      }, "")
    )
} else {
    ## counts
    dimsum_counts$mutname <- vapply(toupper(dimsum_counts[[seqColDimsum]]), function(w) {
        seqToName(w, wtseq, "base", "f")
    }, "")
    ## fitness
    dimsum_fitness$mutname <- vapply(toupper(dimsum_fitness[[seqColDimsum]]), function(w) {
        seqToName(w, wtseq, "base", "f")
    }, "")
}
    
## For Diss_FOS_JUN, only one mutated amino acid per protein was allowed with 
## mutscan - filter out variants with more than one mutated amino acid per 
## protein in the DiMSum results
if (dataset == "Diss_FOS_JUN") {
    dimsum_counts <- dimsum_counts %>% 
        dplyr::filter(
            stringr::str_count(mutname, "f\\.") == 1 & 
                stringr::str_count(mutname, "r\\.") == 1
        )
    dimsum_fitness <- dimsum_fitness %>% 
        dplyr::filter(
            stringr::str_count(mutname, "f\\.") == 1 & 
                stringr::str_count(mutname, "r\\.") == 1
        )
}


dimsum_counts <- as.data.frame(dimsum_counts)
dimsum_fitness <- as.data.frame(dimsum_fitness)

dimsum <- dimsum_counts %>%
    dplyr::select(mutname, .data[[nmutColDimsum]], contains("_counts")) %>%
    dplyr::full_join(dimsum_fitness %>% 
                         dplyr::select(mutname, contains("fitness")), 
                     by = "mutname") %>%
    dplyr::rename(nbr_mutations = .data[[nmutColDimsum]]) %>%
    dplyr::rename_with(.fn = function(x) paste0("DiMSum_", x)) %>%
    dplyr::rename(variant = DiMSum_mutname) %>%
    dplyr::mutate(DiMSum_detected = TRUE)

head(dimsum, 4)

## -----------------------------------------------------------------------------
## Read Enrich2 data
## -----------------------------------------------------------------------------
if (include_enrich2) {
    ## Read Enrich2 counts (if applicable)
    ## This will only be used for Diss_FOS - the code below works for that 
    ## dataset, but may not be generally applicable for other datasets
    enrich2_counts <- read.delim(enrich2CountFile, skip = 3, header = FALSE)
    enrich2header <- read.delim(enrich2CountFile, nrows = 3, header = FALSE)
    colnames(enrich2_counts) <- c("variant", apply(enrich2header[-1, -1], 2, 
                                                   paste, collapse = "_"))
    enrich2_counts <- enrich2_to_mutscan(df = enrich2_counts, 
                                         variantCol = "variant", 
                                         codonPrefix = "f")
    head(enrich2_counts, 4)
    
    ## Read Enrich2 fitness values
    enrich2_fitness <- read.delim(enrich2ResFile, skip = 3, header = FALSE)
    enrich2header <- read.delim(enrich2ResFile, nrows = 3, header = FALSE)
    colnames(enrich2_fitness) <- c("variant", apply(enrich2header[-1, -1], 2, 
                                                    paste, collapse = "_"))
    enrich2_fitness <- enrich2_to_mutscan(df = enrich2_fitness, 
                                          variantCol = "variant", 
                                          codonPrefix = "f")
    head(enrich2_fitness, 4)
    
    colnames(enrich2_counts) <- sub("([0-9])_c_1$", paste0("OU", "\\1"),
                                    sub("([0-9])_c_0$", paste0("IN", "\\1"), 
                                        colnames(enrich2_counts)))
    stopifnot(all(colnames(mutscan_se) %in% colnames(enrich2_counts)))
    for (cn in colnames(mutscan_se)) {
        enrich2_counts[[cn]][is.na(enrich2_counts[[cn]])] <- 0
    }
    enrich2_counts$Enrich2_nbr_mutations = vapply(enrich2_counts$mutname, function(w) {
        if (grepl("WT", w)) {
            0
        } else {
            length(strsplit(w, "_")[[1]])
        }
    }, NA_real_)
    enrich2_counts$Enrich2_detected <- TRUE
    countCols <- which(colnames(enrich2_counts) %in% colnames(mutscan_se))
    colnames(enrich2_counts)[countCols] <- paste0("Enrich2_", colnames(enrich2_counts)[countCols], "_counts")
    
    colnames(enrich2_fitness) <- gsub("(.+)_score", paste0("Enrich2_fitness_\\1"), colnames(enrich2_fitness))
    colnames(enrich2_fitness) <- gsub("(.+)_SE", paste0("Enrich2_se_\\1"), colnames(enrich2_fitness))
    (fitnessCols <- grep("_fitness_", colnames(enrich2_fitness), value = TRUE))
    enrich2_fitness$Enrich2_fitness <- rowMeans(enrich2_fitness[, fitnessCols, drop = FALSE], na.rm = TRUE)
    
    enrich2 <- dplyr::full_join(enrich2_counts %>% dplyr::select(-variant), 
                                enrich2_fitness %>% dplyr::select(-variant), 
                                by = "mutname") %>%
        dplyr::rename(variant = mutname)
} else {
    enrich2 <- NULL
}

head(enrich2, 4)

## -------------------------------------------------------------------------- ##
## Merge
## -------------------------------------------------------------------------- ##
summary_table <- dplyr::full_join(mutscan, dimsum, by = "variant")
if (include_enrich2) {
    summary_table <- summary_table %>%
        dplyr::full_join(enrich2, by = "variant")
}

## Number of mutations
(mutCols <- grep("_nbr_mutations", colnames(summary_table), value = TRUE))
summary_table$nbr_mutations <- rowMeans(
    summary_table[, mutCols, drop = FALSE], 
    na.rm = TRUE)

## Detected in all
(detCols <- grep("_detected", colnames(summary_table), value = TRUE))
summary_table <- summary_table %>% 
    dplyr::mutate(across(all_of(detCols), .fns = function(x) replace(x, is.na(x), FALSE)))
summary_table$detected_in_all <- rowSums(summary_table[, detCols, drop = FALSE]) == length(detCols)

## Average count where detected
(countCols <- grep("_counts", colnames(summary_table), value = TRUE))
summary_table$ave_count <- rowMeans(summary_table[, countCols, drop = FALSE], 
                                    na.rm = TRUE)
summary_table$ave_log10count <- rowMeans(log10(summary_table[, countCols, drop = FALSE] + 1), 
                                         na.rm = TRUE)

## Replace NA counts with 0
summary_table <- summary_table %>%
    dplyr::mutate(across(all_of(countCols), .fns = function(x) replace(x, is.na(x), 0)))

## -------------------------------------------------------------------------- ##
## Save
## -------------------------------------------------------------------------- ##
saveRDS(list(summary_table = summary_table, mutation_type = mutType, 
             dataset = dataset), 
        file = outrds)

date()
sessionInfo()
