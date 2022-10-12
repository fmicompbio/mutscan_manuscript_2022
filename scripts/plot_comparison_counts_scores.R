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
print(outbase)

dat <- readRDS(sumtbl)

## -----------------------------------------------------------------------------
## Summary of number of variants
## -----------------------------------------------------------------------------
dim(dat$summary_table)
summary(dat$summary_table %>% dplyr::select(contains("detected")))

## -----------------------------------------------------------------------------
## UpSet plots
## -----------------------------------------------------------------------------
df <- dat$summary_table %>%
    dplyr::select(contains("_detected"), "ave_log10count", "nbr_mutations")
colnames(df) <- gsub("_detected", "", colnames(df))
png(paste0(outbase, "_upset.png"), width = 8, height = 7, units = "in", res = 200)
print(ComplexUpset::upset(df, 
                          intersect = setdiff(colnames(df), c("ave_log10count", 
                                                              "nbr_mutations")), 
                          annotations = list("Abundance" = (
                              ggplot(mapping = aes(y = ave_log10count)) + 
                                  geom_boxplot() + scale_y_sqrt()
                          ))))
dev.off()

for (nmb in setdiff(unique(df$nbr_mutations), 0)) {
    png(paste0(outbase, "_upset_", nmb, "_mutations.png"), width = 8,
        height = 7, units = "in", res = 200)
    print(ComplexUpset::upset(df[df$nbr_mutations == nmb, ], 
                              intersect = setdiff(colnames(df), c("ave_log10count", 
                                                                  "nbr_mutations")), 
                              annotations = list("Abundance" = (
                                  ggplot(mapping = aes(y = ave_log10count)) + 
                                      geom_boxplot() + scale_y_sqrt()
                              )),
                              wrap = TRUE) + 
              ggtitle(paste0(nmb, " mutated ", dat$mutation_type, 
                             ifelse(nmb == 1, "", "s"))))
    dev.off()
}

## -----------------------------------------------------------------------------
## Pairs plots fitness scores
## -----------------------------------------------------------------------------
df <- dat$summary_table %>% 
    dplyr::filter(detected_in_all) %>%
    dplyr::select(matches("_fitness$"), "nbr_mutations")
df <- df[rowSums(is.na(df[, grep("_fitness$", colnames(df))])) == 0, ]
colnames(df) <- gsub("_fitness", "", colnames(df))
dim(df)
head(df)

png(paste0(outbase, "_pairs.png"), width = 8, height = 8, units = "in", res = 200)
a <- mutscan::plotPairs(se = SummarizedExperiment(
    assays = list(logFC_fitness = df[, colnames(df) != "nbr_mutations"])), 
    selAssay = "logFC_fitness", doLog = FALSE, 
    corrColorRange = c(0.9, 1)) + 
    theme(strip.text = element_text(size = 14))
a$title <- paste0(dat$dataset, ", ", sub("FC_fit", "FC/fit", a$title))
print(a)
dev.off()

for (nmb in setdiff(unique(df$nbr_mutations), 0)) {
    png(paste0(outbase, paste0("_pairs_", nmb, "mutations.png")), width = 8, 
        height = 8, units = "in", res = 200)
    a <- mutscan::plotPairs(se = SummarizedExperiment(
        assays = list(logFC_fitness = df[df$nbr_mutations == nmb, 
                                         colnames(df) != "nbr_mutations"])), 
        selAssay = "logFC_fitness", doLog = FALSE, 
        corrColorRange = c(0.9, 1)) + 
        theme(strip.text = element_text(size = 14))
    a$title <- paste0(dat$dataset, ", ", sub("FC_fit", "FC/fit", a$title),
                      " \n", nmb, " mutated ", dat$mutation_type, 
                      ifelse(nmb == 1, "", "s"))
    print(a)
    dev.off()
}

## -----------------------------------------------------------------------------
## Pairs plots fitness scores, between replicates
## -----------------------------------------------------------------------------
df <- dat$summary_table %>% 
    dplyr::filter(detected_in_all)
for (m in c("DiMSum", "Enrich2")) {
    cols <- grep(paste0(m, "_fitness[0-9]*_"), colnames(df), 
                 value = TRUE)
    print(cols)
    if (length(cols) > 1) {
        png(paste0(outbase, "_", m, "_pairs_replicates.png"), width = 8, 
            height = 8, units = "in", res = 200)
        dfsub <- df[, cols, drop = FALSE]
        dfsub <- dfsub[rowSums(is.na(dfsub)) == 0, ]
        a <- mutscan::plotPairs(se = SummarizedExperiment(
            assays = list(logFC_fitness = dfsub)), 
            selAssay = "logFC_fitness", doLog = FALSE, 
            corrColorRange = c(0.7, 1)) + 
            theme(strip.text = element_text(size = 10))
        a$title <- paste0(dat$dataset, ", ", sub("FC_fit", "FC/fit", a$title), 
                          "\nAcross replicates, ", m)
        print(a)
        dev.off()
    }
}

date()
sessionInfo()
