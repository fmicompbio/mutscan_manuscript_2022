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
    geom_point(size = 4, alpha = 1),
    theme_bw(),
    theme(panel.grid.major.x = element_blank()),
    expand_limits(y = 0)
)

pdf(outpdf, width = 9, height = 4)

g1 <- ggplot(df %>% dplyr::group_by(reads, step, dataset) %>%
                 dplyr::summarize(s = mean(s)), aes(x = reads, y = s)) + 
    sharedPlotArgs + 
    geom_smooth(method = "lm", linetype = "dashed", color = "black", se = FALSE, size = 1) + 
    geom_jitter(data = df, width = 0.1, height = 0, size = 1.5, alpha = 0.5, color = "red") + 
    labs(y = "Total execution time (s)", x = "Number of reads") 

g2 <- ggplot(df %>% dplyr::group_by(reads, step, dataset) %>%
                 dplyr::summarize(max_rss = mean(max_rss)), 
             aes(x = reads, y = max_rss/1024)) + 
    sharedPlotArgs + 
    geom_line(linetype = "dashed", size = 1) + 
    geom_jitter(data = df, width = 0.1, height = 0, size = 1.5, alpha = 0.5, color = "red") + 
    labs(y = "Max RSS (GB)", x = "Number of reads")

g3a <- ggplot(df, aes(x = reads, y = io_char_in/1024)) +
    sharedPlotArgs + 
    geom_smooth(method = "lm") + 
    labs(y = "Total I/O in (GB)", x = "Number of reads")

g3b <- ggplot(df, aes(x = reads, y = io_char_out/1024)) +
    sharedPlotArgs + 
    geom_smooth(method = "lm") + 
    labs(y = "Total I/O out (GB)", x = "Number of reads")

g4 <- ggplot(df %>% dplyr::group_by(reads, step, dataset) %>%
                 dplyr::summarize(mean_load = mean(mean_load)), 
             aes(x = reads, y = mean_load)) +
    sharedPlotArgs + 
    geom_smooth(method = "lm") + 
    geom_jitter(data = df, width = 0.1, height = 0, size = 1.5, alpha = 0.5, color = "red") + 
    labs(y = "Mean load", x = "Number of reads")

g5 <- ggplot(df %>% dplyr::group_by(reads, step, dataset) %>%
                 dplyr::summarize(cpu_time = mean(cpu_time)), 
             aes(x = reads, y = cpu_time)) + 
    sharedPlotArgs + 
    geom_smooth(method = "lm") + 
    geom_jitter(data = df, width = 0.1, height = 0, size = 1.5, alpha = 0.5, color = "red") + 
    labs(y = "Total CPU time (s)", x = "Number of reads")

cowplot::plot_grid(
    g1 + theme(legend.position = "none"), 
    g2 + theme(legend.position = "none"), 
    ncol = 2, align = "hv", labels = c("A", "B"))

dev.off()

date()
sessionInfo()
