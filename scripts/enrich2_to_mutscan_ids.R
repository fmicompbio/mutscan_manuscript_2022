suppressPackageStartupMessages({
    library(rlang)
    library(dplyr)
    library(tibble)
})

## Convert Enrich2 variant IDs to mutscan-style IDs
enrich2_to_mutscan <- function(df, variantCol = "variant", codonPrefix = "f") {
    ## Remove undetermined variants (containing >X) from Enrich2 results and rename variants
    df <- df %>%
        dplyr::filter(!grepl(">X", .data[[variantCol]], fixed = TRUE)) %>%
        dplyr::mutate(mutname = vapply(.data[[variantCol]], function(v) {
            s <- unique(vapply(strsplit(v, ", ")[[1]], 
                               function(w) strsplit(w, " ")[[1]][1], ""))
            s <- s[s != ""]
            paste(vapply(s, function(w) {
                paste0(codonPrefix, ".", 
                       paste(strsplit(strsplit(w, "\\.")[[1]][2], 
                                      "[A,C,G,T]>")[[1]], collapse = "."))
            }, ""), collapse = "_")
        }, "")) %>%
        dplyr::mutate(mutname = replace(mutname, mutname == paste0(codonPrefix, ".NA"), 
                                        paste0(codonPrefix, ".0.WT"))) %>%
        as.data.frame()

    df
}
