#################################################################################################################
#################################################################################################################
### G.Discovery_DIA_FP_MaxLFQ/R6.DIA_overall_view.R
### Project: jcalvet_202110_marato
#################################################################################################################
#################################################################################################################
# Libraries, paths and options

rm(list = ls())
.res <- paste0(.res0, "G.Discovery_DIA_FP_MaxLFQ/R6.DIA_overall_view/")
dir.create(.res, F, T)

#################################################################################################################
# Merge differential expression results for PRM

rm(list = ls())
nd <- readRDS(paste0(.dat, "1_DIA/nd.fp.rds"))
# source("code/functions/assoc.pc.R")
source("code/functions/make_transparent.R")

pc <- prcomp(t(exprs(nd)))
pcvar <- pc$sdev^2 / sum(pc$sdev^2) * 100


nv <- c("1")
nvr <- c("quant.batch")
form <- paste("~ ", paste(c(nv, nvr), collapse = " + "))
design <- model.matrix(as.formula(form), nd)
formr <- paste("~ ", paste(nvr, collapse = " + "))
designr <- model.matrix(as.formula(formr), nd)[, -1, drop = F]
lmf <- lmFit(exprs(nd)[, rownames(design), drop = F], design)
betas <- coef(lmf)
exprs(nd) <- exprs(nd) - betas[, colnames(designr), drop = F] %*% t(designr)


### PC1 vs PC2 ORIGINAL COMPARISON

cols <- c("deepskyblue3", "orange3", "blueviolet", "red3")
pdf(paste0(.res, "pc_1_2_protein_corr.pdf"), width = 8.5, height = 8)
par(mar = c(5, 6, 5, 5))
plot(pc$x[, 1], pc$x[, 2],
  col = makeTransparent(cols[as.numeric(nd$group)], 0.7), pch = 16, cex = 2.25,
  xlim = c(-24, 24), ylim = c(-24, 24),
  main = "Principal components 1 and 2\nProtein level - corrected",
  xlab = paste0("PC1 (", formatC(pcvar[1], format = "f", dig = 1), "%)"),
  ylab = paste0("PC2 (", formatC(pcvar[2], format = "f", dig = 1), "%)"),
  axes = F, cex.lab = 1.75
)
axis(1, lwd = 2.5, cex.axis = 1.25)
axis(2, lwd = 2.5, las = 2, cex.axis = 1.25)
abline(v = 0, lty = 1, lwd = 4, col = makeTransparent("gray80", 0.1))
abline(h = 0, lty = 1, lwd = 4, col = makeTransparent("gray80", 0.1))
box(lwd = 2.5)
legend(
  x = "topleft", inset = 0.01, pch = 16,
  col = rev(cols), pt.cex = 2,
  legend = rev(levels(nd$group)), cex = 1.25
)
dev.off()

# PC1 vs PC2 with ellipses for hospitalized and non-hospitalized clusters
# Added for the draft, to illustrate the separation of hospitalized vs non-hospitalized patients in the PCA space
# parametric, rotated, stretched ellipse helper
ellipse_custom <- function(cx, cy,
                           a, b,
                           angle = 0, # radians
                           npoints = 200) {
  t <- seq(0, 2 * pi, length.out = npoints)
  # basic axis-aligned ellipse
  x0 <- a * cos(t)
  y0 <- b * sin(t)
  # rotate by 'angle'
  xr <- x0 * cos(angle) - y0 * sin(angle)
  yr <- x0 * sin(angle) + y0 * cos(angle)
  # translate to centre
  list(x = xr + cx, y = yr + cy)
}

cols <- c("deepskyblue3", "orange3", "blueviolet", "red3")

# indices for hospitalized and non-hospitalized clinical groups
hosp_idx <- nd$group %in% c("Severe", "Critical")
nonhosp_idx <- nd$group %in% c("Mild", "Asymptomatic")

# centres from the data
cx_hosp <- mean(pc$x[hosp_idx, 1])
cy_hosp <- mean(pc$x[hosp_idx, 2])
cx_nonhosp <- mean(pc$x[nonhosp_idx, 1])
cy_nonhosp <- mean(pc$x[nonhosp_idx, 2])

# Horizontal stretch + left/right shift
# Adjustment of the ellipse parameters.
stretch_x <- 12 # semi-major axis (horizontal)
stretch_y <- 18 # semi-minor axis (vertical)
angle_hosp <- 2.5 # rotation in radians
angle_non <- -2.5

# explicit left/right separation
shift_x <- 6 # how far to push centres apart

ell_hosp <- ellipse_custom(cx_hosp + shift_x, cy_hosp,
  a = stretch_x, b = stretch_y,
  angle = angle_hosp
)

ell_nonhosp <- ellipse_custom(cx_nonhosp - shift_x, cy_nonhosp,
  a = stretch_x, b = stretch_y,
  angle = angle_non
)

pdf(paste0(.res, "pc_1_2_protein_corr_ellipsis.pdf"), width = 8.5, height = 8)
par(mar = c(5, 6, 5, 5))
plot(pc$x[, 1], pc$x[, 2],
  col = makeTransparent(cols[as.numeric(nd$group)], 0.7),
  pch = 16, cex = 2.25,
  xlim = c(-24, 24), ylim = c(-24, 24),
  main = "Principal components 1 and 2\nProtein level - corrected",
  xlab = paste0("PC1 (", formatC(pcvar[1], format = "f", dig = 1), "%)"),
  ylab = paste0("PC2 (", formatC(pcvar[2], format = "f", dig = 1), "%)"),
  axes = F, cex.lab = 1.75
)

## dashed Venn-style ellipses
lines(ell_hosp$x, ell_hosp$y, col = "darkred", lty = 2, lwd = 2)
lines(ell_nonhosp$x, ell_nonhosp$y, col = "darkgreen", lty = 2, lwd = 2)

axis(1, lwd = 2.5, cex.axis = 1.25)
axis(2, lwd = 2.5, las = 2, cex.axis = 1.25)
abline(v = 0, lty = 1, lwd = 4, col = makeTransparent("gray80", 0.1))
abline(h = 0, lty = 1, lwd = 4, col = makeTransparent("gray80", 0.1))
box(lwd = 2.5)

legend(
  x = "topleft", inset = 0.01,
  pch = c(rep(16, length(cols)), NA, NA),
  col = c(rev(cols), "darkred", "darkgreen"),
  pt.cex = c(rep(2, length(cols)), NA, NA),
  lty = c(rep(NA, length(cols)), 2, 2),
  lwd = c(rep(NA, length(cols)), 2, 2),
  legend = c(rev(levels(nd$group)), "Hospitalized", "Non-hospitalized"),
  cex = 0.75 # reduced the size of the legend
)

dev.off()
################################################################################################################
################################################################################################################
################################################################################################################
