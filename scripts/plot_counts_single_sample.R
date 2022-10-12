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

print(sumtbl)
print(sname)
print(outbase)

## -------------------------------------------------------------------------- ##
## Read data
## -------------------------------------------------------------------------- ##
dat <- readRDS(sumtbl)

## -------------------------------------------------------------------------- ##
## Pairs plot for each sample - shared variants
## -------------------------------------------------------------------------- ##
df <- dat$summary_table %>% 
    dplyr::filter(detected_in_all) %>%
    dplyr::select(matches(paste0("_", sname, "_counts")), "nbr_mutations")
colnames(df) <- gsub(paste0("_", sname, "_counts"), "", colnames(df))

png(paste0(outbase, "_", sname, "_pairs_sharedvariants.png"), width = 8, height = 8, units = "in", res = 200)
a <- mutscan::plotPairs(
    se = SummarizedExperiment(assays = list(counts = df[, colnames(df) != "nbr_mutations"])), 
    selAssay = "counts", corrColorRange = c(0.9, 1),
    addIdentityLine = TRUE) + 
    theme(strip.text = element_text(size = 14))
a$title <- paste0(sname, " (", dat$dataset, "), ", a$title, "\nShared variants")
print(a)
dev.off()

## -------------------------------------------------------------------------- ##
## Pairs plot for each sample - all variants
## -------------------------------------------------------------------------- ##
df <- dat$summary_table %>% 
    dplyr::select(matches(paste0("_", sname, "_counts")), "nbr_mutations")
colnames(df) <- gsub(paste0("_", sname, "_counts"), "", colnames(df))
df[is.na(df)] <- 0

png(paste0(outbase, "_", sname, "_pairs_allvariants.png"), width = 8, height = 8, units = "in", res = 200)
a <- mutscan::plotPairs(
    se = SummarizedExperiment(assays = list(counts = df[, colnames(df) != "nbr_mutations"])), 
    selAssay = "counts", corrColorRange = c(0.9, 1),
    addIdentityLine = TRUE) + 
    theme(strip.text = element_text(size = 14))
a$title <- paste0(sname, " (", dat$dataset, "), ", a$title, "\nAll variants")
print(a)
dev.off()

for (nmb in setdiff(unique(df$nbr_mutations[df$nbr_mutations <= 6]), 0)) {
    png(paste0(outbase, "_", sname, "_pairs_allvariants_", nmb, "mutations.png"),
        width = 8, height = 8, units = "in", res = 200)
    a <- mutscan::plotPairs(
        se = SummarizedExperiment(assays = list(counts = df[df$nbr_mutations == nmb, 
                                                            colnames(df) != "nbr_mutations"])), 
        selAssay = "counts", corrColorRange = c(0.9, 1),
        addIdentityLine = TRUE) + 
        theme(strip.text = element_text(size = 14))
    a$title <- paste0(sname, " (", dat$dataset, "), ", a$title, " \n", nmb, 
                      " mutated ", dat$mutation_type, ifelse(nmb == 1, "", "s"), " (all variants)")
    print(a)
    dev.off()
}


date()
sessionInfo()
