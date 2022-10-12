args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

print(inrds)
print(outpdf)

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(cowplot)
    library(ggforce)
    library(RColorBrewer)
})

df <- readRDS(inrds)
df

sharedPlotArgs <- list(
    geom_line(size = 1, linetype = "dashed"),
    geom_point(size = 4, alpha = 1),
    theme_bw(),
    theme(panel.grid.major.x = element_blank()),
    expand_limits(y = 0)
)

pdf(outpdf, width = 10, height = 3)

g1 <- ggplot(df %>% dplyr::group_by(cores, step, dataset) %>%
                 dplyr::summarize(s = mean(s)), aes(x = cores, y = s)) + 
    sharedPlotArgs + 
    geom_jitter(data = df, width = 0.1, height = 0, size = 1.5, alpha = 0.5, color = "red") + 
    labs(y = "Total execution time (s)", x = "Cores") 

g2 <- ggplot(df %>% dplyr::group_by(cores, step, dataset) %>%
                 dplyr::summarize(max_rss = mean(max_rss)), 
             aes(x = cores, y = max_rss/1024)) + 
    sharedPlotArgs + 
    geom_jitter(data = df, width = 0.1, height = 0, size = 1.5, alpha = 0.5, color = "red") + 
    labs(y = "Max RSS (GB)", x = "Cores")

g3a <- ggplot(df, aes(x = cores, y = io_char_in/1024)) + 
    sharedPlotArgs + 
    labs(y = "Total I/O in (GB)", x = "Cores")

g3b <- ggplot(df, aes(x = cores, y = io_char_out/1024)) + 
    sharedPlotArgs + 
    labs(y = "Total I/O out (GB)", x = "Cores")

g4 <- ggplot(df %>% dplyr::group_by(cores, step, dataset) %>%
                 dplyr::summarize(mean_load = mean(mean_load)), 
             aes(x = cores, y = mean_load)) + 
    sharedPlotArgs + 
    geom_jitter(data = df, width = 0.1, height = 0, size = 1.5, alpha = 0.5, color = "red") + 
    labs(y = "Mean load", x = "Cores")

g5 <- ggplot(df %>% dplyr::group_by(cores, step, dataset) %>%
                 dplyr::summarize(cpu_time = mean(cpu_time)), 
             aes(x = cores, y = cpu_time)) + 
    sharedPlotArgs + 
    geom_jitter(data = df, width = 0.1, height = 0, size = 1.5, alpha = 0.5, color = "red") + 
    labs(y = "Total CPU time (s)", x = "Cores")

cowplot::plot_grid(
    g1 + theme(legend.position = "none"), 
    g2 + theme(legend.position = "none"), 
    g4 + theme(legend.position = "none"),
    nrow = 1, align = "hv", labels = c("A", "B", "C"))

dev.off()

date()
sessionInfo()
