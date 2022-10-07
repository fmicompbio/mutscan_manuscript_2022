args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

## List the R version and all packages used in the analyses together with the 
## version, by parsing the files in the "Routdir" and "Reportdir" directories.
## The results are written to the "outtxt" text file.

Routdirs <- strsplit(Routdirs, ",")[[1]]
Reportdirs <- strsplit(Reportdirs, ",")[[1]]
print(Routdirs)
print(Reportdirs)
print(outtxt)
print(outrds)

all_packages <- list()

## process Rout
for (Routdir in Routdirs) {
    if (file.exists(Routdir)) {
        lf <- list.files(Routdir, full.names = TRUE)
        for (f in lf) {
            x <- readLines(f)
            idx1 <- which(x == "> sessionInfo()")
            idx2 <- which(x == "other attached packages:")
            idx3 <- which(x == "loaded via a namespace (and not attached):")
            if (length(idx1) != 0 & length(idx2) != 0 & length(idx3) != 0) {
                thisR <- x[idx1 + 1]
                if (!(thisR %in% names(all_packages))) {
                    all_packages[[thisR]] <- list(files = character(0), pkgs = character(0))
                }
                all_packages[[thisR]][["files"]] <- c(all_packages[[thisR]]$files, f)
                all_packages[[thisR]][["pkgs"]] <- 
                    unique(c(all_packages[[thisR]][["pkgs"]],
                             do.call(c, lapply((idx2 + 1):(idx3 - 2), function(i) {
                                 grep("\\[", setdiff(setdiff(strsplit(x[i], " ")[[1]], " "), ""), 
                                      value = TRUE, invert = TRUE)
                             }))))
            }
        }
    }
}

## process Reportdir
for (Reportdir in Reportdirs) {
    if (file.exists(Reportdir)) {
        lf <- list.files(Reportdir, pattern = "*.knit.md$", full.names = TRUE, recursive = TRUE)
        for (f in lf) {
            x <- readLines(f)
            idx1 <- which(x == "sessionInfo()")
            idx2 <- which(x == "## other attached packages:")
            idx3 <- which(x == "## loaded via a namespace (and not attached):")
            if (length(idx1) != 0 & length(idx2) != 0 & length(idx3) != 0) {
                thisR <- sub("^## ","",x[idx1 + 4])
                if (!(thisR %in% names(all_packages))) {
                    all_packages[[thisR]] <- list(files = character(0), pkgs = character(0))
                }
                all_packages[[thisR]][["files"]] <- c(all_packages[[thisR]]$files, f)
                all_packages[[thisR]][["pkgs"]] <- 
                    unique(c(all_packages[[thisR]][["pkgs"]],
                             do.call(c, lapply((idx2 + 1):(idx3 - 2), function(i) {
                                 grep("^\\[[0-9]+\\]$", setdiff(strsplit(x[i], " ")[[1]], c("", " ", "##")), 
                                      value = TRUE, invert = TRUE)
                             }))))
            }
        }
    }
}

## sort
all_packages <- lapply(all_packages, function(x) {
    x$files <- sort(x$files)
    x$pkgs <- sort(x$pkgs)
    x
})

## output
## ... text
if (exists("outtxt")) {
    if (file.exists(outtxt)) {
        unlink(outtxt)
    }
    fh <- file(outtxt, open = "wt")
    for (thisR in names(all_packages)) {
        cat(thisR, "\n## used in files:\n", file = fh, append = TRUE, sep = "")
        cat(paste(all_packages[[thisR]]$files, collapse = "\n"), "\n", 
            file = fh, append = TRUE, sep = "")
        cat("## with packages:\n", paste(all_packages[[thisR]]$pkgs, collapse = "\n"), 
            "\n\n", file = fh, append = TRUE, sep = "")
    }
    close(fh)
}
## ... rds
if (exists("outrds")) {
    saveRDS(all_packages, file = outrds)
}

date()
sessionInfo()
