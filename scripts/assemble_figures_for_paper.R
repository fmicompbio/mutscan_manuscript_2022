args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
    eval(parse(text = args[[i]]))
}

suppressPackageStartupMessages({
    library(magick)
    library(cowplot)
})

print(outdir)

## -----------------------------------------------------------------------------
## Fig 1 - overview
## -----------------------------------------------------------------------------
fig1 <- ggdraw() + draw_image("paper/figures/mutscan-overview/mutscan-overview.003.png")
pdf(file.path(outdir, "Fig1.pdf"), width = 7, height = 12)
print(fig1)
dev.off()

png(file.path(outdir, "Fig1.png"), width = 7, height = 12, units = "in", res = 200)
print(fig1)
dev.off()

## -----------------------------------------------------------------------------
## Fig 2 - case study
## -----------------------------------------------------------------------------
fig2 <- ggdraw() + draw_image("Diss_FOS_JUN/CaseStudy/Diss_FOS_JUN_CaseStudy_files/figure-html/summary-figure-1.png")
pdf(file.path(outdir, "Fig2.pdf"), width = 17, height = 17)
print(fig2)
dev.off()

png(file.path(outdir, "Fig2.png"), width = 17, height = 17, units = "in", res = 200)
print(fig2)
dev.off()

## -----------------------------------------------------------------------------
## Fig 3 - p-value distributions
## -----------------------------------------------------------------------------
fig3 <- ggdraw() + draw_image(magick::image_read_pdf("benchmark/pdf/nullcomparison.pdf"))
pdf(file.path(outdir, "Fig3.pdf"), width = 10, height = 6)
print(fig3)
dev.off()

png(file.path(outdir, "Fig3.png"), width = 10, height = 6, units = "in", res = 200)
print(fig3)
dev.off()

## -----------------------------------------------------------------------------
## Fig 4 - computational performance
## -----------------------------------------------------------------------------
fig4 <- ggdraw() + draw_image(magick::image_read_pdf("benchmark/pdf/benchmark_results.pdf"))
pdf(file.path(outdir, "Fig4.pdf"), width = 10, height = 12)
print(fig4)
dev.off()

png(file.path(outdir, "Fig4.png"), width = 10, height = 12, units = "in", res = 200)
print(fig4)
dev.off()

## -----------------------------------------------------------------------------
## Fig 5 - comparison of detected variants
## -----------------------------------------------------------------------------
fig5 <- cowplot::plot_grid(
    ggdraw() + draw_image("Diss_FOS/plots/compare_counts_CISOU1_pairs_allvariants.png"),
    ggdraw() + draw_image("Diss_FOS_JUN/plots/compare_counts_TRANSOU1_pairs_allvariants.png"),
    ggdraw() + draw_image("Bolognesi_TDP43_290_331/plots/compare_counts_bTDP1OU7_pairs_allvariants.png"),
    ggdraw() + draw_image("Li_tRNA_sel30/plots/compare_counts_sel30A_pairs_allvariants.png"),
    nrow = 2
)
pdf(file.path(outdir, "Fig5.pdf"), width = 10, height = 10)
print(fig5)
dev.off()

png(file.path(outdir, "Fig5.png"), width = 10, height = 10, units = "in", res = 200)
print(fig5)
dev.off()

## -----------------------------------------------------------------------------
## Fig 6 - comparison of scores
## -----------------------------------------------------------------------------
fig6 <- cowplot::plot_grid(
    ggdraw() + draw_image("Diss_FOS/plots/compare_output_pairs.png"),
    ggdraw() + draw_image("Diss_FOS_JUN/plots/compare_output_pairs.png"),
    ggdraw() + draw_image("Bolognesi_TDP43_290_331/plots/compare_output_pairs.png"),
    ggdraw() + draw_image("Li_tRNA_sel30/plots/compare_output_pairs.png"),
    nrow = 2
)
pdf(file.path(outdir, "Fig6.pdf"), width = 10, height = 10)
print(fig6)
dev.off()

png(file.path(outdir, "Fig6.png"), width = 10, height = 10, units = "in", res = 200)
print(fig6)
dev.off()
