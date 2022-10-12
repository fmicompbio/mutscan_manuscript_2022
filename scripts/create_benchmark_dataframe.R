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

## Don't consider the 'scaling' ones here
infiles <- grep("_scaling_", infiles, invert = TRUE, value = TRUE)

temp <- sub(".txt$", "", basename(infiles))
dataset <- sapply(strsplit(infiles, "/"), .subset, 1)
df1 <- strcapture(pattern = "^([^_]+)_?([^_]*)_?([^_]*)$",
                  x = temp,
                  proto = data.frame(tool = "", step = "", sample = ""))
df1$dataset <- dataset
df1$file = temp
df2 <- do.call(dplyr::bind_rows, lapply(infiles, function(f) {
    read.delim(f) %>%
        dplyr::mutate(file = sub(".txt$", "", basename(f)),
                      dataset = sapply(strsplit(f, "/"), .subset, 1))
}))
df <- dplyr::full_join(df1, df2, by = c("dataset", "file"))

## Some clean-up/clarifications
df <- df %>%
    dplyr::mutate(step = replace(step, tool == "dimsum" & step == "wrap", "01-wrap"),
                  step = replace(step, tool == "dimsum" & step == "steam", "02-steam"),
                  step = replace(step, tool == "mutscan" & step == "digestFastqs", "01-digestFastqs"),
                  step = replace(step, tool == "mutscan" & step == "summarize", "02-summarize"),
                  step = replace(step, tool == "mutscan" & step == "test", "03-test"),
                  step = replace(step, tool == "enrich2" & step == "", "Enrich2")) %>%
    dplyr::mutate(tool = replace(tool, tool == "dimsum", "DiMSum"),
                  tool = replace(tool, tool == "enrich2", "Enrich2")) %>%
    dplyr::filter(step != "nullcomparison")

df
saveRDS(df, file = outrds)

date()
sessionInfo()
