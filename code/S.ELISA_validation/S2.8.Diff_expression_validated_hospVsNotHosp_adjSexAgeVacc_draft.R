#################################################################################################################
#################################################################################################################
### S2.8.Diff_expression_validated_hospVsNotHosp_adjSexAgeVacc_draft.R
### Project: jcalvet_202110_marato
#################################################################################################################
#################################################################################################################
  # Libraries, paths and options

rm(list=ls())
.res <- paste0(.res0, "S.ELISA_validation/S2.8.Diff_expression_validated_hospVsNotHosp_adjSexAgeVacc_draft.R/")
dir.create(.res, F, T)


# Select and prepare data for analysis 2-group design
rm(list=ls())
da <- readRDS(paste0(.dat,"3_ELISA/ELISA.res.long.average.withClinicalData.rds"))

prots <- c("LBP","FGL-1","SAA-2","APOC-1","HPRG","PGLYRP2","SAA-1",
           "SERPINA-3","ORM-1")

d <- da[da$protein %in% prots, ]
d$group <- as.character(d$group)

# Collapse Severe + Critical into Hospitalized; Mild into Not.hospitalized
d$group[d$group %in% c("Severe","Critical")] <- "Hospitalized"
d$group[d$group %in% c("Mild")]             <- "Not.hospitalized"

d$group <- factor(d$group, levels = c("Not.hospitalized","Hospitalized"))
d.list <- split(d, d$protein)
d.list <- d.list[prots]

conts <- c("Hosp-Not.hosp" = "groupHospitalized - groupNot.hospitalized = 0")

# Remove PL5 from SERPINA-3 (old assay, conducted in a first batch of failed experiments)
s <- d.list$`SERPINA-3`$plate == "PL5"
d.list$`SERPINA-3` <- d.list$`SERPINA-3`[!s, ]

# Remove PL1 and PL3 from ORM-1 (failed assay)
s <- d.list$`ORM-1`$plate %in% c("PL1","PL3")
d.list$`ORM-1` <- d.list$`ORM-1`[!s, ]

saveRDS(d.list, paste0(.res,"dataForAnalysis.rds"))
saveRDS(conts,  paste0(.res,"conts.rds"))

#################################################################################################################
# Comparisons
rm(list=ls())

conts  <- readRDS(paste0(.res, "conts.rds"))
d.list <- readRDS(paste0(.res, "dataForAnalysis.rds"))
source(paste0(.wd,"code/functions/make_contrast_list.R"))

levg <- levels(d.list[[1]]$group)
cont_list <- make_contrast_list(conts, levg)

lr <- mclapply(seq_along(d.list), function(i) {
  
  prot <- names(d.list)[i]
  d <- d.list[[i]]
  if (nrow(d) == 0) return(NA)
  
  d$sample.id <- as.character(d$sample.id)
  rownames(d) <- d$sample.id
  d$y <- log2(d$conc.dilCor)
  
  d <- d[complete.cases(d[, c("y","group","center","plate","sex","age","vaccinated")]), ]
  d <- droplevels(d)
  if (nrow(d) == 0 || nlevels(d$group) < 2) return(NA)
  
  # n per group
  n_by_group <- table(d$group)
  n_df <- data.frame(group = names(n_by_group),
                     n     = as.integer(n_by_group),
                     stringsAsFactors = FALSE)
  
  # model
  mm <- try(lm(y ~ group + center + plate + sex + age + vaccinated, data = d), silent = TRUE)
  if (inherits(mm, "try-error")) return(NA)
  
  # emmeans group means (proportional)
  emm_g <- try(emmeans(mm, ~ group, weights = "proportional"), silent = TRUE)
  if (inherits(emm_g, "try-error")) return(NA)
  
  sm_g <- as.data.frame(summary(emm_g))
  dd1 <- data.frame(
    group    = as.character(sm_g$group),
    log2mean = sm_g$emmean,
    log2se   = sm_g$SE,
    stringsAsFactors = FALSE
  )
  
  # global mean (equal-weight average of group means)
  ng <- length(levels(d$group))
  w  <- rep(1/ng, ng)
  gm <- contrast(emm_g, list("global.mean" = w), adjust = "none")
  sm_gm <- as.data.frame(summary(gm))
  dd0 <- c("global.mean" = sm_gm$estimate, "global.se" = sm_gm$SE)
  
  # contrasts (only those applicable; here always Hosp-Not.hosp)
  lv_present <- levels(emm_g@grid$group)
  cont_list_2 <- lapply(cont_list, function(wv) {
    w_named <- setNames(as.numeric(wv), levg)
    w_named[lv_present]
  })
  
  ctr <- try(contrast(emm_g, method = cont_list_2, adjust = "none"), silent = TRUE)
  if (inherits(ctr, "try-error")) return(NA)
  
  sm_c <- as.data.frame(summary(ctr, infer = c(TRUE, TRUE)))
  dd2 <- data.frame(
    log2fc  = sm_c$estimate,
    log2se  = sm_c$SE,
    log2_lo = sm_c$lower.CL,
    log2_hi = sm_c$upper.CL,
    t       = sm_c$t.ratio,
    pv      = sm_c$p.value,
    stringsAsFactors = FALSE,
    row.names = sm_c$contrast
  )
  
  # ycorr for visualization: subtract center+plate term contributions
  tm <- predict(mm, type = "terms")
  batch_terms <- intersect(colnames(tm), c("center", "plate"))
  batch_hat <- if (length(batch_terms)) rowSums(tm[, batch_terms, drop = FALSE]) else 0
  
  ycorr <- d$y - batch_hat
  names(ycorr) <- d$sample.id
  
  list(
    ycorr = ycorr,
    sm1   = dd1,
    sm2   = dd2,
    sm0   = dd0,
    n     = n_df,
    form  = "y ~ group + center + plate + sex + age + vaccinated",
    m     = mm
  )
  
}, mc.cores = 24)

names(lr) <- names(d.list)
lr <- lr[!is.na(lr)]
saveRDS(lr, paste0(.res, "lmod.group_emm.rds"))
saveRDS(cont_list, paste0(.res, "cont_list.rds"))

################################################################################
## Extract n per group
rm(list=ls())
li <- readRDS(paste0(.res, "lmod.group_emm.rds"))

gr_levels <- c("Not.hospitalized","Hospitalized")

dn <- t(sapply(li, function(o) {
  nn <- o$n
  out <- unlist(lapply(gr_levels, function(lv) {
    val <- nn$n[nn$group == lv]
    if (length(val) == 0) val <- NA
    setNames(as.numeric(val), paste0(lv, ".n"))
  }), use.names = TRUE)
  out
}))

saveRDS(dn, paste0(.res, "dcounts_emm.rds"))

################################################################################
## Extract group means
rm(list=ls())
li <- readRDS(paste0(.res, "lmod.group_emm.rds"))

gr_levels <- c("Not.hospitalized","Hospitalized")

dm <- t(sapply(li, function(o) {
  o <- as.data.frame(o$sm1)
  out <- unlist(lapply(gr_levels, function(lv) {
    row <- o[o$group == lv, , drop = FALSE]
    if (nrow(row) == 0) {
      c(setNames(NA_real_, paste0(lv, ".log2mean")),
        setNames(NA_real_, paste0(lv, ".log2se")))
    } else {
      c(setNames(as.numeric(row$log2mean[1]), paste0(lv, ".log2mean")),
        setNames(as.numeric(row$log2se[1]),   paste0(lv, ".log2se")))
    }
  }), use.names = TRUE)
  out
}))

saveRDS(dm, paste0(.res, "dmeans_emm.rds"))

################################################################################
## Extract contrasts
rm(list=ls())
li <- readRDS(paste0(.res, "lmod.group_emm.rds"))
cont_list <- readRDS(paste0(.res, "cont_list.rds"))

contrast_names <- names(cont_list)  # "Hosp-Not.hosp"

dt <- t(sapply(li, function(o) {
  o <- as.data.frame(o$sm2)
  
  # FC + CI
  o$fc    <- 2^(o$log2fc)
  o$fc_lo <- 2^(o$log2_lo)
  o$fc_hi <- 2^(o$log2_hi)
  
  out <- unlist(lapply(contrast_names, function(nm) {
    rr <- o[nm, c("log2fc","log2se","log2_lo","log2_hi","t","pv","fc","fc_lo","fc_hi")]
    setNames(as.numeric(rr), paste0(nm, ".", names(rr)))
  }), use.names = TRUE)
  
  out
}))

dt <- as.data.frame(dt, check.names = FALSE)

# BH across proteins for this contrast
pv_col <- "Hosp-Not.hosp.pv"
dt[["Hosp-Not.hosp.adj.pv"]] <- p.adjust(dt[[pv_col]], method = "BH")

# reorder
ord <- c("Hosp-Not.hosp.log2fc","Hosp-Not.hosp.log2se","Hosp-Not.hosp.log2_lo","Hosp-Not.hosp.log2_hi",
         "Hosp-Not.hosp.t","Hosp-Not.hosp.fc","Hosp-Not.hosp.fc_lo","Hosp-Not.hosp.fc_hi",
         "Hosp-Not.hosp.pv","Hosp-Not.hosp.adj.pv")
dt <- dt[, ord, drop = FALSE]

saveRDS(dt, paste0(.res, "dtests_emm.rds"))


################################################################################
## Extract global means
rm(list=ls())
li <- readRDS(paste0(.res, "lmod.group_emm.rds"))
dmg <- t(sapply(li, function(o) o$sm0))
saveRDS(dmg, paste0(.res, "dmeans.glob_emm.rds"))

################################################################################
## Build dataForAnalysis with ycorr merged
rm(list=ls())
d.list <- readRDS(paste0(.res, "dataForAnalysis.rds"))
li <- readRDS(paste0(.res, "lmod.group_emm.rds"))

lr.corr <- lapply(names(li), function(prot) {
  conc <- data.frame(sample.id = names(li[[prot]]$ycorr),
                     ycorr = as.numeric(li[[prot]]$ycorr),
                     stringsAsFactors = FALSE)
  
  d <- d.list[[prot]]
  d$sample.id <- as.character(d$sample.id)
  d$y <- log2(d$conc.dilCor)
  d$ycorr <- conc$ycorr[match(d$sample.id, conc$sample.id)]
  d
})
names(lr.corr) <- names(li)

saveRDS(lr.corr, paste0(.res, "dataForAnalysis.ycorr_emm.rds"))


################################################################################
# Merge results
rm(list=ls())
dm  <- readRDS(paste0(.res, "dmeans_emm.rds"))
dmg <- readRDS(paste0(.res, "dmeans.glob_emm.rds"))
dn  <- readRDS(paste0(.res, "dcounts_emm.rds"))
dt  <- readRDS(paste0(.res, "dtests_emm.rds"))

nd <- readRDS(paste0(.dat, "1_DIA/nd.fp.rds"))
df <- fData(nd)

dr <- data.frame(gene.symbol=c("LBP","FGL1","SAA2","APOC1","HRG","PGLYRP2","SAA1","SERPINA3","ORM1"),
                 stringsAsFactors=FALSE)
rownames(dr) <- dr$gene.symbol
rownames(dr)[rownames(dr)=="HRG"] <- "HPRG"

dr$protein.id <- df[match(dr$gene.symbol, df$gene.symbol), "protein.group"]
dr$gene.name  <- df[match(dr$gene.symbol, df$gene.symbol), "gene.name"]
dr <- dr[, c("protein.id","gene.symbol","gene.name")]

# rownames to protein labels
rownames(dr)[rownames(dr) %in% c("FGL1","SAA2","APOC1","SAA1","SERPINA3","ORM1")] <-
  c("FGL-1","SAA-2","APOC-1","SAA-1","SERPINA-3","ORM-1")

dr <- cbind(dr, dmg[match(rownames(dr), rownames(dmg)), ])
dr <- cbind(dr, dm [match(rownames(dr), rownames(dm)),  ])
dr <- cbind(dr, dn [match(rownames(dr), rownames(dn)),  ])
dr <- cbind(dr, dt [match(rownames(dr), rownames(dt)),  ])

saveRDS(dr, paste0(.res, "dcomp.rds"))


##################################################################################
# Reporting table
rm(list=ls())

dr <- readRDS(paste0(.res, "dcomp.rds"))


tab <- data.frame(
  Protein      = rownames(dr),
  `UniProt ID` = dr$protein.id,
  `Gene symbol`= dr$gene.symbol,
  `Gene name`  = dr$gene.name,
  `n (Not.hosp / Hosp)` =
    paste0(dr$`Not.hospitalized.n`, " / ", dr$Hospitalized.n),
  
  `Not.hosp mean (SE)` =
    sprintf("%.2f (%.2f)", dr$`Not.hospitalized.log2mean`, dr$`Not.hospitalized.log2se`),
  
  `Hosp mean (SE)` =
    sprintf("%.2f (%.2f)", dr$Hospitalized.log2mean, dr$Hospitalized.log2se),
  
  `Δlog2 (95% CI)` =
    sprintf("%+.2f [%.2f, %.2f]",
            dr$`Hosp-Not.hosp.log2fc`,
            dr$`Hosp-Not.hosp.log2_lo`,
            dr$`Hosp-Not.hosp.log2_hi`),
  
  `FC (95% CI)` =
    sprintf("%.2f [%.2f, %.2f]",
            dr$`Hosp-Not.hosp.fc`,
            dr$`Hosp-Not.hosp.fc_lo`,
            dr$`Hosp-Not.hosp.fc_hi`),
  
  `p` =
    formatC(dr$`Hosp-Not.hosp.pv`, format = "e", digits = 2),
  
  `p (BH)` =
    formatC(dr$`Hosp-Not.hosp.adj.pv`, format = "e", digits = 2),
  
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# numeric vectors used only for bolding
p_nom_num <- dr$`Hosp-Not.hosp.pv`
p_bh_num  <- dr$`Hosp-Not.hosp.adj.pv`

#  GT formatting
gt_tab <- gt(tab) |>
  # tab_header(
  #   title = "ELISA validation by hospitalization status (adjusted for plate and center)",
  #   subtitle = paste0(
  #     "Adjusted means on log2 scale; contrast Hospitalized − Not.hospitalized; ",
  #     "FC = 2^Δlog2; BH-FDR across all proteins in this table"
  #   )
  # ) |>
  cols_align(
    align = "left",
    columns = c(Protein, `UniProt ID`, `Gene symbol`, `Gene name`)
  ) |>
  cols_align(
    align = "right",
    columns = c(`n (Not.hosp / Hosp)`, `Not.hosp mean (SE)`, `Hosp mean (SE)`,
                `Δlog2 (95% CI)`, `FC (95% CI)`, `p`, `p (BH)`)
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      rows = which(!is.na(p_nom_num) & p_nom_num < 0.05),
      columns = `p`
    )
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      rows = which(!is.na(p_bh_num) & p_bh_num < 0.05),
      columns = `p (BH)`
    )
  ) |>
  # tab_footnote(
  #   footnote = "Group means are adjusted marginal means from emmeans (weights='proportional').",
  #   locations = cells_column_labels(columns = c(`Not.hosp mean (SE)`, `Hosp mean (SE)`))
  # ) |>
  # tab_footnote(
  #   footnote = "Δlog2 is the adjusted difference (Hospitalized − Not.hospitalized). FC = 2^Δlog2. 95% CIs are on the log2 scale and transformed to FC.",
  #   locations = cells_column_labels(columns = c(`Δlog2 (95% CI)`, `FC (95% CI)`))
  # ) |>
  # tab_footnote(
  #   footnote = "Per protein: linear model log2(conc.dilCor) ~ group + center + plate + sex + age + vaccinated. BH-FDR is applied across proteins for this contrast.",
  #   locations = cells_column_labels(columns = c(`p (BH)`))
  # ) |>
  tab_options(
    table.font.size = px(12),
    heading.title.font.size = px(14),
    heading.subtitle.font.size = px(12),
    data_row.padding = px(4)
  )

gt_tab
gtsave(gt_tab, filename = paste0(.res, "ELISA_table_hospComparison_adjSexAgeVacc.html"))
gtsave(gt_tab, filename = paste0(.res, "ELISA_table_hospComparison_adjSexAgeVacc.png"),
       vwidth = 3000, vheight = 1600, zoom = 2
)





################################################################################################################
# Plots for proteins
################################################################################################################

rm(list=ls())

dr    <- readRDS(paste0(.res, "dcomp.rds"))
ycorr <- readRDS(paste0(.res, "dataForAnalysis.ycorr_emm.rds"))

grp_levels <- c("Not.hospitalized","Hospitalized")
cols <- c("Not.hospitalized" = "orange3", "Hospitalized" = "red4")

unit_map <- c(
  "SERPINA-3" = "ng/ml",
  "ORM-1"     = "ng/ml",
  "SAA-1"     = "ng/ml",
  "LBP"       = "ug/ml",
  "FGL-1"     = "ng/ml",
  "SAA-2"     = "ng/ml",
  "PGLYRP2"   = "ng/ml",
  "APOC-1"    = "ug/ml",
  "HPRG"      = "ng/ml",
  "MSN"       = "pg/ml"
)
facet_lab <- function(p) {
  u <- unit_map[p]
  if (is.null(u) || !nzchar(u)) p else paste0(p, "\n(", u, ")")
}
stars <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01,  "**",
                ifelse(p < 0.05,  "*", "")))
}

#  long plot data 
df_plot <- do.call(rbind, lapply(names(ycorr), function(p) {
  d <- ycorr[[p]]
  d <- d[!is.na(d$ycorr) & !is.na(d$group), , drop=FALSE]
  d <- d[d$group %in% grp_levels, , drop=FALSE]
  if (nrow(d) == 0) return(NULL)
  d$group <- factor(as.character(d$group), levels = grp_levels)
  d$protein <- p
  d$facet <- facet_lab(p)
  d$y_plot <- 2^(d$ycorr)
  d[, c("protein","facet","group","y_plot","sample.id")]
}))
df_plot <- as.data.frame(df_plot, stringsAsFactors = FALSE)
df_plot$group <- factor(df_plot$group, levels = grp_levels)

facet_order <- vapply(unique(df_plot$protein), facet_lab, character(1))
df_plot$facet <- factor(df_plot$facet, levels = facet_order)

#  mean overlay 
mean_df <- do.call(rbind, lapply(rownames(dr), function(p) {
  if (!p %in% unique(df_plot$protein)) return(NULL)
  r <- dr[p, , drop=FALSE]
  mean_log2 <- c(r$`Not.hospitalized.log2mean`, r$Hospitalized.log2mean)
  se_log2   <- c(r$`Not.hospitalized.log2se`,   r$Hospitalized.log2se)
  
  lo <- mean_log2 - 1.96*se_log2
  hi <- mean_log2 + 1.96*se_log2
  
  data.frame(
    facet = facet_lab(p),
    group = factor(grp_levels, levels = grp_levels),
    mean  = 2^(mean_log2),
    lo    = 2^(lo),
    hi    = 2^(hi),
    stringsAsFactors = FALSE
  )
}))
mean_df$facet <- factor(mean_df$facet, levels = facet_order)

#  brackets rule: BH if sig else raw only if nominal sig 
p_nom <- dr$`Hosp-Not.hosp.pv`;      names(p_nom) <- rownames(dr)
p_bh  <- dr$`Hosp-Not.hosp.adj.pv`;  names(p_bh)  <- rownames(dr)

ymax <- tapply(df_plot$y_plot, df_plot$facet, max, na.rm=TRUE)

bh_prots  <- names(p_bh)[!is.na(p_bh) & p_bh < 0.05]
nom_prots <- setdiff(names(p_nom)[!is.na(p_nom) & p_nom < 0.05], bh_prots)

make_br_df <- function(prots, mult0, mult1, mult_lab, label_fun) {
  if (!length(prots)) {
    return(data.frame(facet=factor(character(0), levels=facet_order),
                      x1=numeric(0), x2=numeric(0), y0=numeric(0), y=numeric(0),
                      ylab=numeric(0), lab=character(0), stringsAsFactors=FALSE))
  }
  facs <- vapply(prots, facet_lab, character(1))
  ytop <- as.numeric(ymax[facs])
  data.frame(
    facet = factor(facs, levels=facet_order),
    x1=1, x2=2,
    y0=ytop*mult0,
    y =ytop*mult1,
    ylab=ytop*mult_lab,
    lab=label_fun(prots),
    stringsAsFactors=FALSE
  )
}

bh_df <- make_br_df(bh_prots, 1.08, 1.12, 1.15, function(prots) {
  pv <- as.numeric(p_bh[prots])
  paste0(stars(pv), "  BH p=", formatC(pv, format = "e", digits = 2))
})

nom_df <- make_br_df(nom_prots, 1.02, 1.06, 1.09, function(prots) {
  pv <- as.numeric(p_nom[prots])
  paste0(stars(pv), "  p=", formatC(pv, format = "e", digits = 2))
})

cap_raw <- paste0(
  "Each panel shows batch-corrected ELISA concentrations (points) by hospitalization status. ",
  "Batch correction was performed on log2 scale by subtracting fitted center and plate contributions; values are displayed as 2^ycorr. ",
  "Adjusted group means (points) and 95% CIs (bars) are from y ~ group + center + plate + sex + age + vaccinated using emmeans (weights='proportional'). ",
  "Brackets show BH significance when BH p<0.05; otherwise raw p is shown only if nominally significant."
)
cap <- paste(strwrap(cap_raw, width=150), collapse="\n")

p <- ggplot(df_plot, aes(x=group, y=y_plot)) +
  geom_violin(aes(fill=group), alpha=0.06, color="grey35",
              width=0.72, trim=FALSE, linewidth=0.25) +
  geom_jitter(aes(color=group), width=0.10, height=0,
              alpha=0.75, size=1.15) +
  geom_errorbar(data=mean_df, aes(x=group, y=mean, ymin=lo, ymax=hi),
                inherit.aes=FALSE, position=position_nudge(x=0.25),
                width=0.05, linewidth=0.65, color="black") +
  geom_point(data=mean_df, aes(x=group, y=mean, fill=group),
             inherit.aes=FALSE, position=position_nudge(x=0.25),
             shape=21, size=2.8, stroke=0.45, color="black") +
  # nominal bracket
  geom_segment(data=nom_df, aes(x=x1, xend=x2, y=y0, yend=y0),
               inherit.aes=FALSE, linewidth=0.7, color="black") +
  geom_segment(data=nom_df, aes(x=x1, xend=x1, y=y0, yend=y),
               inherit.aes=FALSE, linewidth=0.7, color="black") +
  geom_segment(data=nom_df, aes(x=x2, xend=x2, y=y0, yend=y),
               inherit.aes=FALSE, linewidth=0.7, color="black") +
  geom_text(data=nom_df, aes(x=1.5, y=ylab, label=lab),
            inherit.aes=FALSE, size=3.1, color="black") +
  # BH bracket
  geom_segment(data=bh_df, aes(x=x1, xend=x2, y=y0, yend=y0),
               inherit.aes=FALSE, linewidth=0.7, color="black") +
  geom_segment(data=bh_df, aes(x=x1, xend=x1, y=y0, yend=y),
               inherit.aes=FALSE, linewidth=0.7, color="black") +
  geom_segment(data=bh_df, aes(x=x2, xend=x2, y=y0, yend=y),
               inherit.aes=FALSE, linewidth=0.7, color="black") +
  geom_text(data=bh_df, aes(x=1.5, y=ylab, label=lab),
            inherit.aes=FALSE, size=3.1, color="black") +
  facet_wrap(~facet, scales="free_y") +
  expand_limits(y=0) +
  coord_cartesian(ylim=c(0, NA)) +
  scale_color_manual(values=cols) +
  scale_fill_manual(values=cols) +
  labs(x=NULL, y="Batch-corrected concentration (original scale)", caption=cap) +
  theme_classic(base_size=12) +
  theme(
    legend.position="none",
    strip.background=element_rect(fill="grey92", color="black", linewidth=0.4),
    strip.text=element_text(face="bold", size=10),
    panel.border=element_rect(fill=NA, color="black", linewidth=0.5),
    axis.line=element_line(color="black", linewidth=0.4),
    axis.ticks=element_line(color="black", linewidth=0.4),
    axis.text.x=element_text(size=9),
    plot.margin=ggplot2::margin(10,10,10,10),
    plot.caption=element_text(hjust=0, size=9, margin=ggplot2::margin(t=8))
  )

p

ggsave(
  filename = paste0(.res, "ELISA_Hosp_vs_NotHosp_ycorr_originalScale_facets_adjSexAgeVacc.png"),
  plot = p,
  width = 1000, height = 600, units = "px", dpi = 96
)








################################################################################################################
################################################################################################################
################################################################################################################
