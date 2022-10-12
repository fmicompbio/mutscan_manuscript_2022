args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
    library(dplyr)
})

print(fileglob)
print(outrds)

infiles <- Sys.glob(fileglob)

temp <- sub(".txt$", "", basename(infiles))
dataset <- sapply(strsplit(infiles, "/"), .subset, 1)
df1 <- strcapture(pattern = "^([^_]+)_?([^_]*)_?([^_]*)_?([^_]*)_?([^_]*)$",
                  x = temp,
                  proto = data.frame(tool = "", scaling = "", step = "", 
                                     sample = "", reads = ""))
df1$dataset <- dataset
df1$file = temp
df2 <- do.call(dplyr::bind_rows, lapply(infiles, function(f) {
    read.delim(f) %>%
        dplyr::mutate(file = sub(".txt$", "", basename(f)))
}))
df <- dplyr::full_join(df1, df2, by = "file")

## Some clean-up/clarifications
df <- df %>%
    dplyr::mutate(reads = as.numeric(sub("reads", "", reads)))

df

saveRDS(df, file = outrds)

date()
sessionInfo()
