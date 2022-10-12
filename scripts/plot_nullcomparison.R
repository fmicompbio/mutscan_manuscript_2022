args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2, lib.loc = "rlib_devel")
})

print(fileglob)
print(outpdf)

colvec <- c(edgeR = "#3c8f15", limma = "#643a9e")

infiles <- Sys.glob(fileglob)

dataset <- sapply(strsplit(infiles, "/"), .subset, 1)
names(infiles) <- dataset
df <- do.call(rbind, lapply(dataset, function(ds) {
    data.frame(dataset = ds, readRDS(infiles[ds]))
})) %>%
    dplyr::mutate(dataset = factor(dataset, levels = c("Diss_FOS", "Diss_FOS_JUN",
                                                       "Bolognesi_TDP43_290_331", 
                                                       "Li_tRNA_sel30")))
head(df)

## Number of variants per method/dataset/subset
nvars <- df %>%
    dplyr::group_by(dataset, method, g1) %>%
    dplyr::summarize(n = length(pvalue), .groups = "drop") %>%
    dplyr::select(dataset, n) %>%
    dplyr::distinct() %>%
    dplyr::mutate(n = paste0("p = ", n))

pdf(outpdf, width = 10, height = 6)
ggplot(df) + 
    geom_density(aes(x = pvalue, fill = method, color = method, group = g1),
                 alpha = 0.05, bounds = c(0, 1)) + 
    geom_density(aes(x = pvalue), bounds = c(0, 1), color = "grey20",
                 size = 1.25) + 
    facet_grid(method ~ dataset) + 
    scale_fill_manual(values = colvec) +
    scale_color_manual(values = colvec) + 
    coord_cartesian(xlim = c(0, 1)) + 
    theme_bw() + theme(legend.position = "none") + 
    labs(x = "Nominal p-value", y = "Density") + 
    geom_text(data = nvars, aes(x = 1, y = 0, label = n), 
              hjust = "right", vjust = "bottom")
dev.off()

date()
sessionInfo()
