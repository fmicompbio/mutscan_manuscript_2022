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
})

df <- readRDS(inrds)

## Change level order, define colors
df <- df %>% 
    dplyr::mutate(step = factor(step, levels = c("01-wrap", "02-steam",
                                                 "Enrich2", "01-digestFastqs",
                                                 "02-summarize", "03-test"))) %>%
    dplyr::select(-h.m.s, -file) %>%
    dplyr::mutate(dataset = factor(dataset, levels = c("Diss_FOS", "Diss_FOS_JUN",
                                                       "Bolognesi_TDP43_290_331",
                                                       "Li_tRNA_sel30")))

stepCols <- c("01-wrap" = "#0b03fc",
              "02-steam" = "#b1aef2",
              "Enrich2" = "#d1a104",
              "01-digestFastqs" = "#f7202f",
              "02-summarize" = "#f5828a",
              "03-test" = "#fad4d7")

## Summarize digestFastqs performance by averaging across the five input samples
## in the Li data set (only need to run this once with mutscan)
df0 <- df %>%
    dplyr::filter(tool == "mutscan" & 
                      step == "01-digestFastqs" &
                      dataset == "Li_tRNA_sel30" & 
                      grepl("input", sample)) %>%
    dplyr::mutate(sample = "input") %>%
    dplyr::group_by(tool, step, sample, dataset) %>%
    dplyr::summarize(across(everything(), mean), 
                     .groups = "drop")
df <- df %>%
    dplyr::filter(!(tool == "mutscan" & 
                        step == "01-digestFastqs" &
                        dataset == "Li_tRNA_sel30" & 
                        grepl("input", sample))) %>%
    dplyr::bind_rows(df0)

pdf(outpdf, width = 10, height = 12)

g1 <- ggplot(df, aes(x = tool, y = s, fill = step)) + 
    geom_bar(stat = "identity") + 
    theme_bw() + theme(panel.grid.major.x = element_blank()) + 
    facet_row(~ dataset, scales = "free", space = "free") + 
    labs(y = "Total execution time (s)", x = "") + 
    scale_fill_manual(values = stepCols)
    
g2 <- ggplot(df %>% dplyr::group_by(tool, step, dataset) %>%
                 dplyr::summarize(max_rss = max(max_rss)), 
             aes(x = tool, y = max_rss/1024, fill = step)) + 
    geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) + 
    theme_bw() + theme(panel.grid.major.x = element_blank()) + 
    facet_row(~ dataset, scales = "free", space = "free") + 
    labs(y = "Max RSS (GB)", x = "") + 
    scale_fill_manual(values = stepCols)

g3 <- ggplot() + 
    geom_bar(data = df %>% dplyr::select(tool, step, sample, dataset, io_char_in, io_char_out) %>%
                 tidyr::gather(key = "io", value = "value", io_char_in, io_char_out) %>%
                 dplyr::filter(io == "io_char_in"), 
             aes(x = tool, y = value/1024, fill = step), stat = "identity",
             position = position_nudge(x = -0.21), width = 0.4) +
    geom_text(data = df %>% dplyr::select(tool, step, sample, dataset, io_char_in, io_char_out) %>%
                  tidyr::gather(key = "io", value = "value", io_char_in, io_char_out) %>%
                  dplyr::filter(io == "io_char_in"), label = "I",
              aes(x = tool, y = 0), position = position_nudge(x = -0.21), vjust = -2) + 
    geom_bar(data = df %>% dplyr::select(tool, step, sample, dataset, io_char_in, io_char_out) %>%
                 tidyr::gather(key = "io", value = "value", io_char_in, io_char_out) %>%
                 dplyr::filter(io == "io_char_out"), 
             aes(x = tool, y = value/1024, fill = step), stat = "identity",
             position = position_nudge(x = 0.21), width = 0.4) + 
    geom_text(data = df %>% dplyr::select(tool, step, sample, dataset, io_char_in, io_char_out) %>%
                  tidyr::gather(key = "io", value = "value", io_char_in, io_char_out) %>%
                  dplyr::filter(io == "io_char_out"), label = "O",
              aes(x = tool, y = 0), position = position_nudge(x = 0.21), vjust = -2) + 
    theme_bw() + theme(panel.grid.major.x = element_blank()) + 
    facet_row(~ dataset, scales = "free", space = "free") + 
    labs(y = "Total I/O (GB)", x = "") + 
    scale_fill_manual(values = stepCols)

g3a <- ggplot(df, aes(x = tool, y = io_char_in/1024, fill = step)) + 
    geom_bar(stat = "identity") + 
    theme_bw() + theme(panel.grid.major.x = element_blank()) + 
    facet_row(~ dataset, scales = "free", space = "free") + 
    labs(y = "Total I/O in (GB)", x = "") + 
    scale_fill_manual(values = stepCols)

g3b <- ggplot(df, aes(x = tool, y = io_char_out/1024, fill = step)) + 
    geom_bar(stat = "identity") + 
    theme_bw() + theme(panel.grid.major.x = element_blank()) + 
    facet_row(~ dataset, scales = "free", space = "free") + 
    labs(y = "Total I/O out (GB)", x = "") + 
    scale_fill_manual(values = stepCols)

g4 <- ggplot(df %>% dplyr::group_by(tool, step, dataset) %>%
                 dplyr::summarize(mean_load = mean(mean_load)), 
             aes(x = tool, y = mean_load, fill = step)) + 
    geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) + 
    theme_bw() + theme(panel.grid.major.x = element_blank()) + 
    facet_row(~ dataset, scales = "free", space = "free") + 
    labs(y = "Mean load (with 10 cores)", x = "") + 
    scale_fill_manual(values = stepCols)

g5 <- ggplot(df, aes(x = tool, y = cpu_time, fill = step)) + 
    geom_bar(stat = "identity") + 
    theme_bw() + theme(panel.grid.major.x = element_blank()) + 
    facet_row(~ dataset, scales = "free", space = "free") + 
    labs(y = "Total CPU time (s)", x = "") + 
    scale_fill_manual(values = stepCols)

cowplot::plot_grid(
    g1 + theme(legend.position = "none"), 
    g2 + theme(legend.position = "none"), 
    g3 + theme(legend.position = "none"),
    g4 + theme(legend.position = "none"),
    get_legend(g1 + theme(legend.position = "bottom")),
    ncol = 1, align = "hl",
    rel_heights = c(1, 1, 1, 1, 0.3)
)

dev.off()

date()
sessionInfo()
