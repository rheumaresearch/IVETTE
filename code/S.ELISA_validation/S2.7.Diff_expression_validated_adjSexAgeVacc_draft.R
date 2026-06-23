#################################################################################################################
#################################################################################################################
### S2.7.Diff_expression_validated_adjSexAgeVacc_draft.R
### Project: jcalvet_202110_marato
#################################################################################################################
#################################################################################################################
  # Libraries, paths and options

rm(list=ls())
.res <- paste0(.res0, "S.ELISA_validation/S2.7.Diff_expression_validated_adjSexAgeVacc_draft/")
dir.create(.res, F, T)

#################################################################################################################
# Select and prepare data for analysis 2-group design
rm(list=ls())
da <- readRDS(paste0(.dat,"3_ELISA/ELISA.res.long.average.withClinicalData.rds"))

prots <- c("LBP","FGL-1","SAA-2","APOC-1","HPRG","PGLYRP2","SAA-1",
           "SERPINA-3","ORM-1","MSN")

d <- da[da$protein %in% prots,]
d$group <- droplevels(d$group)
d.list <- split(d,d$protein)
d.list <- d.list[prots]


### Contrasts

conts <- c("Severe-Mild"="groupSevere - groupMild = 0",
           "Critical-Mild"="groupCritical - groupMild = 0",
           "Critical-Severe"="groupCritical - groupSevere = 0")


# Remove PL5 from SERPINA-3
s <- d.list$`SERPINA-3`$plate=="PL5"
d.list$`SERPINA-3` <- d.list$`SERPINA-3`[!s,]

# Remove PL1 and PL3 from ORM-1 (very different standard curves - failed assay)
s <- d.list$`ORM-1`$plate%in%c("PL1","PL3")
d.list$`ORM-1` <- d.list$`ORM-1`[!s,]

# Sensitivity - check whether results change substantialy when removing outliers
# despite being ok in the standard curve
# SAA2 - Remove outliers
# s <- d.list$`SAA-2`$conc.dilCor<8912
# d.list$`SAA-2` <- d.list$`SAA-2`[!s,]

# HPRG - Remove outliers
# s <- d.list$HPRG$conc.dilCor<4096
# d.list$HPRG <- d.list$HPRG[!s,]


saveRDS(d.list,paste0(.res,"dataForAnalysis.rds"))
saveRDS(conts, paste0(.res, "conts.rds"))


#################################################################################################################
# Comparisons
rm(list=ls())
conts <- readRDS( paste0(.res, "conts.rds"))
d.list <- readRDS(paste0(.res, "dataForAnalysis.rds"))
source(paste0(.wd,"code/functions/make_contrast_list.R"))
# levg <- levels(d.list$`SAA-2`$group)
levg <- levels(d.list$`SERPINA-3`$group)

cont_list <- make_contrast_list(conts, levg)


lr <- mclapply(seq_along(d.list), function(i) {
  
  prot <- names(d.list)[i]
  d <- d.list[[i]]
  d <- d[d$protein == prot, , drop = FALSE]
  if (nrow(d) == 0) return(NA)
  
  d$sample.id <- as.character(d$sample.id)
  rownames(d) <- d$sample.id
  d$y <- log2(d$conc.dilCor)
  
  d <- d[complete.cases(d[, c("y","group","center","plate","sex","age","vaccinated")]), ]
  d <- droplevels(d)
  if (nrow(d) < 3 || nlevels(d$group) < 2) return(NA)
  
  ## n per group
  n_by_group <- table(d$group)
  n_df <- data.frame(
    group = names(n_by_group),
    n     = as.integer(n_by_group),
    stringsAsFactors = FALSE
  )
  
  ## Fit model
  mm <- try(lm(y ~ group + center + plate + sex + age + vaccinated, data = d), silent = TRUE)
  if (inherits(mm, "try-error")) return(NA)
  
  ## emmeans
  emm_g <- try(emmeans(mm, ~ group, weights = "proportional"), silent = TRUE)
  if (inherits(emm_g, "try-error")) return(NA)
  
  # Group means
  sm_g <- as.data.frame(summary(emm_g))
  dd1 <- data.frame(
    group    = as.character(sm_g$group),
    log2mean = sm_g$emmean,
    log2se   = sm_g$SE,
    stringsAsFactors = FALSE
  )
  
  ## Global mean: equal-weight average of group means
  ng <- length(levels(d$group))
  w  <- rep(1/ng, ng)
  gm <- contrast(emm_g, list("global.mean" = w), adjust = "none")
  sm_gm <- as.data.frame(summary(gm))
  dd0 <- c(
    "global.mean" = sm_gm$estimate,
    "global.se"   = sm_gm$SE
  )
  
  ## Contrasts
  lv_present <- levels(emm_g@grid$group)
  
  # for each contrast definition restrict to lv_present order
  cont_list_present <- lapply(cont_list, function(w) {
    w_named <- setNames(as.numeric(w), levg)
    as.numeric(w_named[lv_present])
  })  
  
  applicable <- sapply(cont_list, function(w) {
    w_named <- setNames(as.numeric(w), levg)
    
    nz <- names(w_named)[abs(w_named) > 0]
    all(nz %in% lv_present)
  })
  
  cont_list_present <- cont_list_present[applicable]
  
  if (length(cont_list_present) == 0) {
    dd2 <- data.frame(
      log2fc = numeric(0), log2se = numeric(0),
      log2_lo = numeric(0), log2_hi = numeric(0),
      t = numeric(0), pv = numeric(0),
      stringsAsFactors = FALSE
    )
  } else {
    ctr <- try(contrast(emm_g, method = cont_list_present, adjust = "none"), silent = TRUE)
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
  }
  
  
  ## ycorr for visualization: subtract center + plate term contributions
  tm <- predict(mm, type = "terms")
  batch_terms <- intersect(colnames(tm), c("center", "plate","sex","age","vaccinated"))
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

gr_levels_all <- c("Mild","Severe","Critical")

dn <- t(sapply(li, function(o) {
  nn <- o$n
  out <- unlist(lapply(gr_levels_all, function(lv) {
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

gr_levels_all <- c("Mild","Severe","Critical")

lr_means <- lapply(li, function(o) o$sm1)

dm <- t(sapply(lr_means, function(o) {
  o <- as.data.frame(o)
  out <- unlist(lapply(gr_levels_all, function(lv) {
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

contrast_names <- names(cont_list)

lr_tests <- lapply(li, function(o) o$sm2)

dt_raw <- t(sapply(lr_tests, function(o) {
  o <- as.data.frame(o)
  
  # compute FC and CI on original scale
  if (nrow(o) > 0) {
    o$fc    <- 2^(o$log2fc)
    o$fc_lo <- 2^(o$log2_lo)
    o$fc_hi <- 2^(o$log2_hi)
  }
  
  out <- unlist(lapply(contrast_names, function(nm) {
    if (!nm %in% rownames(o)) {
      # contrast not available for this protein -> NA block
      rr <- rep(NA_real_, 9)
      names(rr) <- c("log2fc","log2se","log2_lo","log2_hi","t","pv","fc","fc_lo","fc_hi")
    } else {
      rr <- o[nm, c("log2fc","log2se","log2_lo","log2_hi","t","pv","fc","fc_lo","fc_hi")]
      rr <- as.numeric(rr)
      names(rr) <- c("log2fc","log2se","log2_lo","log2_hi","t","pv","fc","fc_lo","fc_hi")
    }
    setNames(rr, paste0(nm, ".", names(rr)))
  }), use.names = TRUE)
  
  out
}))

dt_raw <- as.data.frame(dt_raw, check.names = FALSE)
saveRDS(dt_raw, paste0(.res, "dtests_emm_all_raw.rds"))

###################################################################################
## Primary comparison: Severe-Mild with BH across proteins where it exists

rm(list=ls())
dt_raw <- readRDS(paste0(.res, "dtests_emm_all_raw.rds"))

pv_col <- "Severe-Mild.pv"
adj_col <- "Severe-Mild.adj.pv"

dt_primary <- dt_raw[, grep("^Severe-Mild\\.", colnames(dt_raw), value = TRUE), drop = FALSE]

# BH only over proteins with non-missing pv
pv <- dt_primary[[pv_col]]
adj <- rep(NA_real_, length(pv))
ok <- !is.na(pv)
adj[ok] <- p.adjust(pv[ok], method = "BH")
dt_primary[[adj_col]] <- adj

# order cols nicely
ord <- c("Severe-Mild.log2fc","Severe-Mild.log2se","Severe-Mild.log2_lo","Severe-Mild.log2_hi",
         "Severe-Mild.t","Severe-Mild.fc","Severe-Mild.fc_lo","Severe-Mild.fc_hi",
         "Severe-Mild.pv","Severe-Mild.adj.pv")
ord <- ord[ord %in% colnames(dt_primary)]
dt_primary <- dt_primary[, ord, drop = FALSE]

saveRDS(dt_primary, paste0(.res, "dtests_primary_emm.rds"))


###################################################################################
## Secondary comparisons: Critical-Mild and Critical-Severe
## BH is done separately per contrast across proteins where it exists
rm(list=ls())
dt_raw <- readRDS(paste0(.res, "dtests_emm_all_raw.rds"))

make_secondary <- function(dt_raw, cname) {
  cols <- grep(paste0("^", cname, "\\."), colnames(dt_raw), value = TRUE)
  if (length(cols) == 0) return(NULL)
  
  dt <- dt_raw[, cols, drop = FALSE]
  pv_col  <- paste0(cname, ".pv")
  adj_col <- paste0(cname, ".adj.pv")
  
  pv <- dt[[pv_col]]
  adj <- rep(NA_real_, length(pv))
  ok <- !is.na(pv)
  adj[ok] <- p.adjust(pv[ok], method = "BH")
  dt[[adj_col]] <- adj
  
  ord <- c(paste0(cname, c(".log2fc",".log2se",".log2_lo",".log2_hi",".t",
                           ".fc",".fc_lo",".fc_hi",".pv",".adj.pv")))
  ord <- ord[ord %in% colnames(dt)]
  dt <- dt[, ord, drop = FALSE]
  dt
}

dt_cm <- make_secondary(dt_raw, "Critical-Mild")
dt_cs <- make_secondary(dt_raw, "Critical-Severe")

dt_secondary <- cbind(
  if (!is.null(dt_cm)) dt_cm else NULL,
  if (!is.null(dt_cs)) dt_cs else NULL
)

saveRDS(dt_secondary, paste0(.res, "dtests_secondary_emm.rds"))


###################################################################################
## Extract global means

rm(list=ls())
d.list <- readRDS(paste0(.res, "dataForAnalysis.rds"))
li <- readRDS(paste0(.res, "lmod.group_emm.rds"))

dmg <- t(sapply(li, function(o) o$sm0))
saveRDS(dmg, paste0(.res, "dmeans.glob_emm.rds"))

###################################################################################
## Build dataForAnalysis with ycorr merged
rm(list=ls())
d.list <- readRDS(paste0(.res, "dataForAnalysis.rds"))
li <- readRDS(paste0(.res, "lmod.group_emm.rds"))


lr.corr <- lapply(names(li), function(prot) {
  conc <- data.frame(sample.id = names(li[[prot]]$ycorr),
                     ycorr     = as.numeric(li[[prot]]$ycorr),
                     stringsAsFactors = FALSE)
  
  d <- d.list[[prot]]
  d$sample.id <- as.character(d$sample.id)
  d$y <- log2(d$conc.dilCor)
  
  d$ycorr <- conc$ycorr[match(d$sample.id, conc$sample.id)]
  d
})
names(lr.corr) <- names(li)

saveRDS(lr.corr, paste0(.res, "dataForAnalysis.ycorr_emm.rds"))


###################################################################################
## Merge main results table (primary comparison: Severe-Mild)
rm(list=ls())
dm  <- readRDS(paste0(.res, "dmeans_emm.rds"))
dmg <- readRDS(paste0(.res, "dmeans.glob_emm.rds"))
dn  <- readRDS(paste0(.res, "dcounts_emm.rds"))
dt1 <- readRDS(paste0(.res, "dtests_primary_emm.rds"))

nd <- readRDS(paste0(.dat, "1_DIA/nd.fp.rds"))
df <- fData(nd)

dr <- data.frame(gene.symbol=c("LBP", "FGL1", "SAA2", "APOC1",
                               "HRG","PGLYRP2","SAA1",
                               "SERPINA3", "ORM1"),
                 stringsAsFactors=FALSE)
rownames(dr) <- dr[, 1]
rownames(dr)[rownames(dr)=='HRG'] <- "HPRG"

dr$protein.id <- df[match(dr$gene.symbol, df$gene.symbol), "protein.group"]
dr <- dr[, 2:1]
dr$gene.name <- df[match(dr$gene.symbol, df$gene.symbol), "gene.name"]
rownames(dr)[rownames(dr) %in% c("FGL1","SAA2","APOC1","SAA1","SERPINA3","ORM1")] <-
  c("FGL-1","SAA-2","APOC-1","SAA-1","SERPINA-3","ORM-1")

# merge blocks
dr_main <- cbind(dr, dmg[match(rownames(dr), rownames(dmg)), ])
dr_main <- cbind(dr_main, dm [match(rownames(dr), rownames(dm)),  ])
dr_main <- cbind(dr_main, dn [match(rownames(dr), rownames(dn)),  ])
dr_main <- cbind(dr_main, dt1[match(rownames(dr), rownames(dt1)), ])

saveRDS(dr_main, paste0(.res, "dcomp_primary_Severe-Mild.rds"))


###################################################################################
## Merge SECONDARY results table (Critical contrasts only)
rm(list=ls())
dm  <- readRDS(paste0(.res, "dmeans_emm.rds"))
dmg <- readRDS(paste0(.res, "dmeans.glob_emm.rds"))
dn  <- readRDS(paste0(.res, "dcounts_emm.rds"))
dt2 <- readRDS(paste0(.res, "dtests_secondary_emm.rds"))

nd <- readRDS(paste0(.dat, "1_DIA/nd.fp.rds"))
df <- fData(nd)

sec_prots <- intersect(rownames(dt2)[rowSums(!is.na(dt2)) > 0], rownames(dm))

dr2 <- data.frame(gene.symbol=c("SERPINA3","ORM1","MSN"), stringsAsFactors=FALSE)
rownames(dr2) <- dr2$gene.symbol
dr2$protein.id <- df[match(dr2$gene.symbol, df$gene.symbol), "protein.group"]
dr2$gene.name  <- df[match(dr2$gene.symbol, df$gene.symbol), "gene.name"]

rownames(dr2)[rownames(dr2) == "SERPINA3"] <- "SERPINA-3"
rownames(dr2)[rownames(dr2) == "ORM1"]     <- "ORM-1"

dr2 <- dr2[intersect(rownames(dr2), rownames(dt2)), , drop = FALSE]
dr2 <- dr2[, c("protein.id","gene.symbol","gene.name"), drop = FALSE]

dr_sec <- cbind(dr2, dmg[match(rownames(dr2), rownames(dmg)), ])
dr_sec <- cbind(dr_sec, dm [match(rownames(dr2), rownames(dm)),  ])
dr_sec <- cbind(dr_sec, dn [match(rownames(dr2), rownames(dn)),  ])
dr_sec <- cbind(dr_sec, dt2[match(rownames(dr2), rownames(dt2)), ])

saveRDS(dr_sec, paste0(.res, "dcomp_secondary_Critical_contrasts.rds"))


################################################################################################################
# Table for reporting - primary (Severe vs Mild)
rm(list=ls())

dr <- readRDS(paste0(.res, "dcomp_primary_Severe-Mild.rds"))

tab <- data.frame(
  Protein = rownames(dr),
  `UniProt ID`  = dr$protein.id,
  `Gene symbol` = dr$gene.symbol,
  `Gene name`   = dr$gene.name,
  `n (Mild / Severe)` = paste0(dr$Mild.n, " / ", dr$Severe.n),
  `Mild mean (SE)`   = sprintf("%.2f (%.2f)", dr$Mild.log2mean, dr$Mild.log2se),
  `Severe mean (SE)` = sprintf("%.2f (%.2f)", dr$Severe.log2mean, dr$Severe.log2se),
  `Δlog2 (95% CI)` = sprintf("%+.2f [%.2f, %.2f]",
                             dr$`Severe-Mild.log2fc`,
                             dr$`Severe-Mild.log2_lo`,
                             dr$`Severe-Mild.log2_hi`),
  `FC (95% CI)` = sprintf("%.2f [%.2f, %.2f]",
                          dr$`Severe-Mild.fc`,
                          dr$`Severe-Mild.fc_lo`,
                          dr$`Severe-Mild.fc_hi`),
  `p` = formatC(dr$`Severe-Mild.pv`, format = "e", digits = 2),
  `p (BH)` = formatC(dr$`Severe-Mild.adj.pv`, format = "e", digits = 2),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# For styling (numeric)
p_nom_num <- dr$`Severe-Mild.pv`
p_bh_num  <- dr$`Severe-Mild.adj.pv`

gt_tab_primary <-
  gt(tab) %>%
  tab_header(
    title = "ELISA validation (adjusted for plate and center)",
    subtitle = "Primary comparison: Severe − Mild (adjusted means on log2 scale; FC = 2^Δlog2; BH-FDR across proteins with this contrast)"
  ) %>%
  cols_align(
    align = "left",
    columns = c(Protein, `UniProt ID`, `Gene symbol`, `Gene name`)
  ) %>%
  cols_align(
    align = "right",
    columns = c(`n (Mild / Severe)`, `Mild mean (SE)`, `Severe mean (SE)`,
                `Δlog2 (95% CI)`, `FC (95% CI)`, `p`, `p (BH)`)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = which(p_nom_num < 0.05), columns = `p`)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = which(p_bh_num < 0.05), columns = `p (BH)`)
  ) %>%
  tab_footnote(
    footnote = "Linear model per protein on log2(conc.dilCor): y ~ group + center + plate. Marginal means via emmeans with weights='proportional'.",
    locations = cells_column_labels(columns = c(`Mild mean (SE)`, `Severe mean (SE)`))
  ) %>%
  tab_footnote(
    footnote = "Δlog2 is the adjusted difference (Severe − Mild). 95% CI on log2 scale; FC and CI are 2^Δlog2 and 2^CI.",
    locations = cells_column_labels(columns = c(`Δlog2 (95% CI)`, `FC (95% CI)`))
  ) %>%
  tab_options(
    table.font.size = px(12),
    heading.title.font.size = px(14),
    heading.subtitle.font.size = px(12),
    data_row.padding = px(4)
  )

gt_tab_primary

gtsave(gt_tab_primary, filename = paste0(.res, "ELISA_table_primary_adjSexAgeVacc.html"))

################################################################################################################
# Table for reporting - SUPPLEMENTARY (Critical contrasts)
rm(list=ls())

dr <- readRDS(paste0(.res, "dcomp_secondary_Critical_contrasts.rds"))

# helper to safely format numeric vectors with NAs
fmt_e <- function(x) ifelse(is.na(x), "", formatC(x, format = "e", digits = 2))
fmt_num <- function(x, dig = 2) ifelse(is.na(x), "", formatC(x, format = "f", digits = dig))
fmt_pm <- function(est, lo, hi, dig = 2, signed = FALSE) {
  ifelse(is.na(est), "",
         sprintf(paste0(if (signed) "%+." else "%.", dig, "f [%.", dig, "f, %.", dig, "f]"),
                 est, lo, hi))
}

# build columns
tab <- data.frame(
  Protein = rownames(dr),
  `UniProt ID`  = dr$protein.id,
  `Gene symbol` = dr$gene.symbol,
  `Gene name`   = dr$gene.name,
  `n (Mild / Severe / Critical)` = paste0(dr$Mild.n, " / ", dr$Severe.n, " / ", dr$Critical.n),
  
  # Critical - Mild
  `Δlog2 C−M (95% CI)` = fmt_pm(dr$`Critical-Mild.log2fc`, dr$`Critical-Mild.log2_lo`, dr$`Critical-Mild.log2_hi`, dig = 2, signed = TRUE),
  `FC C−M (95% CI)`    = fmt_pm(dr$`Critical-Mild.fc`,   dr$`Critical-Mild.fc_lo`,   dr$`Critical-Mild.fc_hi`,   dig = 2, signed = FALSE),
  `p C−M`              = fmt_e(dr$`Critical-Mild.pv`),
  `p C−M (BH)`         = fmt_e(dr$`Critical-Mild.adj.pv`),
  
  # Critical - Severe
  `Δlog2 C−S (95% CI)` = fmt_pm(dr$`Critical-Severe.log2fc`, dr$`Critical-Severe.log2_lo`, dr$`Critical-Severe.log2_hi`, dig = 2, signed = TRUE),
  `FC C−S (95% CI)`    = fmt_pm(dr$`Critical-Severe.fc`,   dr$`Critical-Severe.fc_lo`,   dr$`Critical-Severe.fc_hi`,   dig = 2, signed = FALSE),
  `p C−S`              = fmt_e(dr$`Critical-Severe.pv`),
  `p C−S (BH)`         = fmt_e(dr$`Critical-Severe.adj.pv`),
  
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# numeric p for styling (handle missing cols safely)
p_cm_nom <- dr$`Critical-Mild.pv`
p_cm_bh  <- dr$`Critical-Mild.adj.pv`
p_cs_nom <- dr$`Critical-Severe.pv`
p_cs_bh  <- dr$`Critical-Severe.adj.pv`

gt_tab_secondary <-
  gt(tab) %>%
  tab_header(
    title = "ELISA validation (adjusted for plate and center)",
    subtitle = "Supplementary comparisons: Critical vs Mild and Critical vs Severe (BH-FDR computed within each contrast family across proteins with that contrast)"
  ) %>%
  cols_align(
    align = "left",
    columns = c(Protein, `UniProt ID`, `Gene symbol`, `Gene name`)
  ) %>%
  cols_align(
    align = "right",
    columns = setdiff(colnames(tab), c("Protein", "UniProt ID", "Gene symbol", "Gene name"))
  ) %>%
  # bold nominal and BH p-values per contrast
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = which(!is.na(p_cm_nom) & p_cm_nom < 0.05), columns = `p C−M`)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = which(!is.na(p_cm_bh) & p_cm_bh < 0.05), columns = `p C−M (BH)`)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = which(!is.na(p_cs_nom) & p_cs_nom < 0.05), columns = `p C−S`)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = which(!is.na(p_cs_bh) & p_cs_bh < 0.05), columns = `p C−S (BH)`)
  ) %>%
  tab_footnote(
    footnote = "Linear model per protein on log2(conc.dilCor): y ~ group + center + plate. Marginal means via emmeans with weights='proportional'.",
    locations = cells_column_labels(columns = c(`Δlog2 C−M (95% CI)`, `Δlog2 C−S (95% CI)`))
  ) %>%
  tab_footnote(
    footnote = "Δlog2 contrasts are adjusted differences; FC and CIs are computed as 2^Δlog2 and 2^CI.",
    locations = cells_column_labels(columns = c(`FC C−M (95% CI)`, `FC C−S (95% CI)`))
  ) %>%
  tab_options(
    table.font.size = px(12),
    heading.title.font.size = px(14),
    heading.subtitle.font.size = px(12),
    data_row.padding = px(4)
  )

gt_tab_secondary

gtsave(gt_tab_secondary, filename = paste0(.res, "ELISA_table_secondary_adjSexAgeVacc.html"))



################################################################################################################
# Plots for proteins
################################################################################################################
rm(list=ls())

d.list <- readRDS(paste0(.res, "dataForAnalysis.rds"))
dt_raw      <- readRDS(paste0(.res, "dtests_emm_all_raw.rds"))
dt_primary  <- readRDS(paste0(.res, "dtests_primary_emm.rds"))
dt_secondary<- readRDS(paste0(.res, "dtests_secondary_emm.rds"))

# Merge p(adjusted) columns into dt_raw
add_cols_by_rownames <- function(base, add) {
  if (is.null(add) || ncol(add) == 0) return(base)
  ix <- match(rownames(base), rownames(add))
  cbind(base, add[ix, setdiff(colnames(add), colnames(base)), drop = FALSE])
}
dt_all <- dt_raw
dt_all <- add_cols_by_rownames(dt_all, dt_primary)
dt_all <- add_cols_by_rownames(dt_all, dt_secondary)

# ycorr list
ycorr <- readRDS(paste0(.res, "dataForAnalysis.ycorr_emm.rds"))

# Group levels + colors
grp_levels <- c("Mild", "Severe", "Critical")
cols <- c("Mild" = "orange3", "Severe" = "blueviolet", "Critical" = "red3")

#  UNITS 
unit_map <- c(
  "SERPINA-3" = "ug/ml",
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

# Helper: stars 
stars <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01,  "**",
                ifelse(p < 0.05,  "*", "")))
}

# Long data for points/violins (original scale)
df_plot <- do.call(rbind, lapply(names(ycorr), function(p) {
  d <- ycorr[[p]]
  d <- d[!is.na(d$ycorr) & !is.na(d$group), , drop = FALSE]
  
  # keep only known groups; drop others
  d <- d[d$group %in% grp_levels, , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  
  # drop unused factor levels in this protein
  d$group <- factor(as.character(d$group), levels = grp_levels)
  d <- d[!is.na(d$group), , drop = FALSE]
  
  d$protein <- p
  d$facet   <- facet_lab(p)
  d$y_plot  <- 2^(d$ycorr)  # back-transform to concentration
  d[, c("protein","facet","group","y_plot","sample.id")]
}))
df_plot <- as.data.frame(df_plot, stringsAsFactors = FALSE)
df_plot$group <- factor(df_plot$group, levels = grp_levels)

# facet order: alphabetical by protein name in the list
facet_order <- unique(vapply(unique(df_plot$protein), facet_lab, character(1)))
df_plot$facet <- factor(df_plot$facet, levels = facet_order)


# 2) Overlay: adjusted means + 95% CI
dm <- readRDS(paste0(.res, "dmeans_emm.rds"))

mean_df <- do.call(rbind, lapply(rownames(dm), function(p) {
  
  if (!p %in% unique(df_plot$protein)) return(NULL)
  
  # group-wise mean/se on log2 scale (may include NA for missing group)
  m1 <- c(Mild     = dm[p, "Mild.log2mean"],
          Severe   = dm[p, "Severe.log2mean"],
          Critical = dm[p, "Critical.log2mean"])
  
  s1 <- c(Mild     = dm[p, "Mild.log2se"],
          Severe   = dm[p, "Severe.log2se"],
          Critical = dm[p, "Critical.log2se"])
  
  # keep only available groups
  ok <- !is.na(m1) & !is.na(s1)
  if (!any(ok)) return(NULL)
  
  grp <- factor(names(m1)[ok], levels = grp_levels)
  
  lo_log2 <- m1[ok] - 1.96 * s1[ok]
  hi_log2 <- m1[ok] + 1.96 * s1[ok]
  
  data.frame(
    protein = p,
    facet = facet_lab(p),
    group = grp,
    mean  = 2^(m1[ok]),
    lo    = 2^(lo_log2),
    hi    = 2^(hi_log2),
    stringsAsFactors = FALSE
  )
}))
mean_df <- as.data.frame(mean_df, stringsAsFactors = FALSE)
mean_df$facet <- factor(mean_df$facet, levels = facet_order)
mean_df$group <- factor(mean_df$group, levels = grp_levels)


# 3) Brackets for multiple comparisons (Severe-Mild, Critical-Mild, Critical-Severe)
ymax <- tapply(df_plot$y_plot, df_plot$facet, max, na.rm = TRUE)

make_brackets <- function(contrast_name, x1, x2,
                          y_mult0_nom = 1.02, y_mult1_nom = 1.06, y_mult_lab_nom = 1.09,
                          y_mult0_bh  = 1.08, y_mult1_bh  = 1.12, y_mult_lab_bh  = 1.15,
                          label_prefix = "") {
  
  pv_col  <- paste0(contrast_name, ".pv")
  adj_col <- paste0(contrast_name, ".adj.pv")
  
  # if raw p not available -> return empty
  if (!(pv_col %in% colnames(dt_all))) {
    empty <- data.frame(facet=factor(character(0), levels=facet_order),
                        x1=numeric(0), x2=numeric(0), y0=numeric(0), y=numeric(0),
                        ylab=numeric(0), lab=character(0), stringsAsFactors = FALSE)
    return(list(nom = empty, bh = empty))
  }
  
  pv <- dt_all[[pv_col]]
  names(pv) <- rownames(dt_all)
  
  adj <- NULL
  if (adj_col %in% colnames(dt_all)) {
    adj <- dt_all[[adj_col]]
    names(adj) <- rownames(dt_all)
  }
  
  # BH significant
  bh_prots <- character(0)
  if (!is.null(adj)) {
    bh_prots <- names(adj)[!is.na(adj) & adj < 0.05]
  }
  
  # Nominal-only significant (raw p < 0.05 but NOT BH-significant)
  nom_prots <- names(pv)[!is.na(pv) & pv < 0.05]
  nom_prots <- setdiff(nom_prots, bh_prots)
  
  # helper to build df
  build_df <- function(prots, y0m, ym, ylabm, lab_fun) {
    if (!length(prots)) {
      return(data.frame(facet=factor(character(0), levels=facet_order),
                        x1=numeric(0), x2=numeric(0), y0=numeric(0), y=numeric(0),
                        ylab=numeric(0), lab=character(0), stringsAsFactors = FALSE))
    }
    facs <- vapply(prots, facet_lab, character(1))
    ytop <- as.numeric(ymax[facs])
    data.frame(
      facet = factor(facs, levels = facet_order),
      x1 = x1, x2 = x2,
      y0 = ytop * y0m,
      y  = ytop * ym,
      ylab = ytop * ylabm,
      lab = lab_fun(prots),
      stringsAsFactors = FALSE
    )
  }
  
  # labels
  lab_bh <- function(prots) {
    paste0(label_prefix, stars(adj[prots]), "  BH p=",
           formatC(adj[prots], format="e", digits=2))
  }
  lab_nom <- function(prots) {
    paste0(label_prefix, stars(pv[prots]), "  p=",
           formatC(pv[prots], format="e", digits=2))
  }
  
  bh_df  <- build_df(bh_prots,  y_mult0_bh,  y_mult1_bh,  y_mult_lab_bh,  lab_bh)
  nom_df <- build_df(nom_prots, y_mult0_nom, y_mult1_nom, y_mult_lab_nom, lab_nom)
  
  list(nom = nom_df, bh = bh_df)
}

# Build bracket sets per contrast 
br_SM <- make_brackets("Severe-Mild",     x1=1, x2=2, label_prefix = "")
br_CM <- make_brackets("Critical-Mild",   x1=1, x2=3, label_prefix = "C−M: ")
br_CS <- make_brackets("Critical-Severe", x1=2, x2=3, label_prefix = "C−S: ")

# convenience objects for plotting layers
nom_SM <- br_SM$nom; bh_SM <- br_SM$bh
nom_CM <- br_CM$nom; bh_CM <- br_CM$bh
nom_CS <- br_CS$nom; bh_CS <- br_CS$bh


# Caption
cap_raw <- paste0(
  "Each panel shows batch-corrected ELISA concentrations (points) by clinical group (Mild/Severe/Critical when available). ",
  "Batch correction was performed on the log2 scale by subtracting fitted center and plate contributions from a linear model; ",
  "concentrations are displayed after back-transforming (2^ycorr). ",
  "Adjusted group means (points) and 95% CIs (bars) are from models y ~ group + center + plate + sex + age + vaccinated using emmeans (weights='proportional'). ",
  "Brackets indicate nominal (p<0.05) and BH-FDR significance (within each contrast family across proteins where the contrast is available)."
)
cap <- paste(strwrap(cap_raw, width = 150), collapse = "\n")

#Plot
p <- ggplot(df_plot, aes(x = group, y = y_plot)) +
  
  geom_violin(
    aes(fill = group),
    alpha = 0.06, color = "grey35",
    width = 0.72, trim = FALSE, linewidth = 0.25
  ) +
  geom_jitter(
    aes(color = group),
    width = 0.10, height = 0,
    alpha = 0.75, size = 1.15
  ) +
  
  geom_errorbar(
    data = mean_df,
    aes(x = group, y = mean, ymin = lo, ymax = hi),
    inherit.aes = FALSE,
    position = position_nudge(x = 0.25),
    width = 0.05, linewidth = 0.65, color = "black"
  ) +
  geom_point(
    data = mean_df,
    aes(x = group, y = mean, fill = group),
    inherit.aes = FALSE,
    position = position_nudge(x = 0.25),
    shape = 21, size = 2.8, stroke = 0.45, color = "black"
  ) +
  
  # ---------- Bracket layers: nominal then BH (for each contrast) ----------
# Nominal SM
geom_segment(data = nom_SM, aes(x = x1, xend = x2, y = y0, yend = y0),
             inherit.aes = FALSE, linewidth = 0.7, color = "black") +
  geom_segment(data = nom_SM, aes(x = x1, xend = x1, y = y0, yend = y),
               inherit.aes = FALSE, linewidth = 0.7, color = "black") +
  geom_segment(data = nom_SM, aes(x = x2, xend = x2, y = y0, yend = y),
               inherit.aes = FALSE, linewidth = 0.7, color = "black") +
  geom_text(data = nom_SM, aes(x = (x1 + x2)/2, y = ylab, label = lab),
            inherit.aes = FALSE, size = 3.1, color = "black") +
  
  # BH SM
  geom_segment(data = bh_SM, aes(x = x1, xend = x2, y = y0, yend = y0),
               inherit.aes = FALSE, linewidth = 0.7, color = "black") +
  geom_segment(data = bh_SM, aes(x = x1, xend = x1, y = y0, yend = y),
               inherit.aes = FALSE, linewidth = 0.7, color = "black") +
  geom_segment(data = bh_SM, aes(x = x2, xend = x2, y = y0, yend = y),
               inherit.aes = FALSE, linewidth = 0.7, color = "black") +
  geom_text(data = bh_SM, aes(x = (x1 + x2)/2, y = ylab, label = lab),
            inherit.aes = FALSE, size = 3.1, color = "black") +
  
  # Nominal CM
  geom_segment(data = nom_CM, aes(x = x1, xend = x2, y = y0, yend = y0),
               inherit.aes = FALSE, linewidth = 0.7, color = "black") +
  geom_segment(data = nom_CM, aes(x = x1, xend = x1, y = y0, yend = y),
               inherit.aes = FALSE, linewidth = 0.7, color = "black") +
  geom_segment(data = nom_CM, aes(x = x2, xend = x2, y = y0, yend = y),
               inherit.aes = FALSE, linewidth = 0.7, color = "black") +
  geom_text(data = nom_CM, aes(x = (x1 + x2)/2, y = ylab, label = lab),
            inherit.aes = FALSE, size = 3.1, color = "black") +
  
  # BH CM
  geom_segment(data = bh_CM, aes(x = x1, xend = x2, y = y0, yend = y0),
               inherit.aes = FALSE, linewidth = 0.7, color = "black") +
  geom_segment(data = bh_CM, aes(x = x1, xend = x1, y = y0, yend = y),
               inherit.aes = FALSE, linewidth = 0.7, color = "black") +
  geom_segment(data = bh_CM, aes(x = x2, xend = x2, y = y0, yend = y),
               inherit.aes = FALSE, linewidth = 0.7, color = "black") +
  geom_text(data = bh_CM, aes(x = (x1 + x2)/2, y = ylab, label = lab),
            inherit.aes = FALSE, size = 3.1, color = "black") +
  
  # Nominal CS
  geom_segment(data = nom_CS, aes(x = x1, xend = x2, y = y0, yend = y0),
               inherit.aes = FALSE, linewidth = 0.7, color = "black") +
  geom_segment(data = nom_CS, aes(x = x1, xend = x1, y = y0, yend = y),
               inherit.aes = FALSE, linewidth = 0.7, color = "black") +
  geom_segment(data = nom_CS, aes(x = x2, xend = x2, y = y0, yend = y),
               inherit.aes = FALSE, linewidth = 0.7, color = "black") +
  geom_text(data = nom_CS, aes(x = (x1 + x2)/2, y = ylab, label = lab),
            inherit.aes = FALSE, size = 3.1, color = "black") +
  
  # BH CS
  geom_segment(data = bh_CS, aes(x = x1, xend = x2, y = y0, yend = y0),
               inherit.aes = FALSE, linewidth = 0.7, color = "black") +
  geom_segment(data = bh_CS, aes(x = x1, xend = x1, y = y0, yend = y),
               inherit.aes = FALSE, linewidth = 0.7, color = "black") +
  geom_segment(data = bh_CS, aes(x = x2, xend = x2, y = y0, yend = y),
               inherit.aes = FALSE, linewidth = 0.7, color = "black") +
  geom_text(data = bh_CS, aes(x = (x1 + x2)/2, y = ylab, label = lab),
            inherit.aes = FALSE, size = 3.1, color = "black") +
  
  facet_wrap(~ facet, scales = "free_y") +
  expand_limits(y = 0) +
  coord_cartesian(ylim = c(0, NA)) +
  scale_color_manual(values = cols, drop = FALSE) +
  scale_fill_manual(values = cols, drop = FALSE) +
  # labs(x = NULL, y = "Batch-corrected concentration (original scale)", caption = cap) +
  labs(x = NULL, y = "Batch-corrected concentration (original scale)") +
  
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey92", color = "black", linewidth = 0.4),
    strip.text = element_text(face = "bold", size = 10),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    axis.text.x = element_text(size = 9),
    # plot.caption = element_text(hjust = 0, size = 9, margin = ggplot2::margin(t = 8)),
    plot.margin = ggplot2::margin(10, 10, 10, 10)
  )

p
ggsave(
  filename = paste0(.res, "ELISA_groups_ycorr_originalScale_facets_adjSexAgeVacc.png"),
  plot = p,
  width = 1000, height = 600, units = "px", dpi = 96
)



################################################################################################################
# Plot in LOG2 scale (ycorr), y axis in original scale
################################################################################################################
rm(list=ls())
ycorr <- readRDS(paste0(.res, "dataForAnalysis.ycorr_emm.rds"))
dm    <- readRDS(paste0(.res, "dmeans_emm.rds"))

dt_raw       <- readRDS(paste0(.res, "dtests_emm_all_raw.rds"))
dt_primary   <- readRDS(paste0(.res, "dtests_primary_emm.rds"))
dt_secondary <- readRDS(paste0(.res, "dtests_secondary_emm.rds"))

# Merge adjusted p columns into dt_raw
add_cols_by_rownames <- function(base, add) {
  if (is.null(add) || ncol(add) == 0) return(base)
  ix <- match(rownames(base), rownames(add))
  new_cols <- setdiff(colnames(add), colnames(base))
  if (!length(new_cols)) return(base)
  cbind(base, add[ix, new_cols, drop = FALSE])
}

dt_all <- dt_raw
dt_all <- add_cols_by_rownames(dt_all, dt_primary)
dt_all <- add_cols_by_rownames(dt_all, dt_secondary)

# Settings
grp_levels <- c("Mild", "Severe", "Critical")
cols <- c("Mild"="orange3","Severe"="blueviolet","Critical"="red3")

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
fmt_p <- function(p) {
  p <- as.numeric(p)
  ifelse(is.na(p), NA_character_,
         ifelse(p < 1e-16, "<1e-16", formatC(p, format="e", digits=2)))
}

# y tick labels: convert log2 to original scale
fmt_num <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_character_, length(x))
  ok <- !is.na(x)
  x0 <- x[ok]
  dec <- ifelse(x0 < 10, 2,
                ifelse(x0 < 100, 1, 0))
  out[ok] <- format(round(x0, dec), big.mark = ",", trim = TRUE, scientific = FALSE)
  out
}

# 1) Data in long format
df_plot <- do.call(rbind, lapply(names(ycorr), function(p) {
  d <- ycorr[[p]]
  d <- d[!is.na(d$ycorr) & !is.na(d$group), , drop=FALSE]
  d <- d[d$group %in% grp_levels, , drop=FALSE]
  if (nrow(d) == 0) return(NULL)
  
  d$group   <- factor(as.character(d$group), levels = grp_levels)
  d <- d[!is.na(d$group), , drop=FALSE]
  
  d$protein <- p
  d$facet   <- facet_lab(p)
  d$y_log2  <- d$ycorr
  
  d[, c("protein","facet","group","y_log2","sample.id")]
}))
df_plot <- as.data.frame(df_plot, stringsAsFactors = FALSE)
df_plot$group <- factor(df_plot$group, levels = grp_levels)

facet_order <- unique(vapply(unique(df_plot$protein), facet_lab, character(1)))
df_plot$facet <- factor(df_plot$facet, levels = facet_order)


#A djusted means + 95% CI (on log2 scale)

mean_df <- do.call(rbind, lapply(rownames(dm), function(p) {
  
  if (!p %in% unique(df_plot$protein)) return(NULL)
  
  m1 <- c(Mild     = dm[p, "Mild.log2mean"],
          Severe   = dm[p, "Severe.log2mean"],
          Critical = dm[p, "Critical.log2mean"])
  
  s1 <- c(Mild     = dm[p, "Mild.log2se"],
          Severe   = dm[p, "Severe.log2se"],
          Critical = dm[p, "Critical.log2se"])
  
  ok <- !is.na(m1) & !is.na(s1)
  if (!any(ok)) return(NULL)
  
  grp <- factor(names(m1)[ok], levels = grp_levels)
  lo  <- m1[ok] - 1.96 * s1[ok]
  hi  <- m1[ok] + 1.96 * s1[ok]
  
  data.frame(
    protein = p,
    facet   = facet_lab(p),
    group   = grp,
    mean    = m1[ok],
    lo      = lo,
    hi      = hi,
    stringsAsFactors = FALSE
  )
}))
mean_df <- as.data.frame(mean_df, stringsAsFactors = FALSE)
mean_df$facet <- factor(mean_df$facet, levels = facet_order)
mean_df$group <- factor(mean_df$group, levels = grp_levels)

# Brackets for multiple comparisons
ymax <- tapply(df_plot$y_log2, df_plot$facet, max, na.rm = TRUE)
make_brackets_log2 <- function(contrast_name, x1, x2,
                               add0_nom = 0.05, add1_nom = 0.10, addlab_nom = 0.14,
                               add0_bh  = 0.13, add1_bh  = 0.18, addlab_bh  = 0.22,
                               label_prefix = "") {
  
  pv_col  <- paste0(contrast_name, ".pv")
  adj_col <- paste0(contrast_name, ".adj.pv")
  
  empty <- data.frame(facet=factor(character(0), levels=facet_order),
                      x1=numeric(0), x2=numeric(0), y0=numeric(0), y=numeric(0),
                      ylab=numeric(0), lab=character(0), stringsAsFactors=FALSE)
  
  if (!(pv_col %in% colnames(dt_all))) return(list(nom=empty, bh=empty))
  
  pv <- as.numeric(dt_all[[pv_col]])
  names(pv) <- rownames(dt_all)
  
  adj <- NULL
  if (adj_col %in% colnames(dt_all)) {
    adj <- as.numeric(dt_all[[adj_col]])
    names(adj) <- rownames(dt_all)
  }
  
  bh_prots <- character(0)
  if (!is.null(adj)) bh_prots <- names(adj)[!is.na(adj) & adj < 0.05]
  
  nom_prots <- names(pv)[!is.na(pv) & pv < 0.05]
  nom_prots <- setdiff(nom_prots, bh_prots)
  
  build_df <- function(prots, add0, add1, addlab, lab_vec) {
    if (!length(prots)) return(empty)
    facs <- vapply(prots, facet_lab, character(1))
    ytop <- as.numeric(ymax[facs])
    data.frame(
      facet=factor(facs, levels=facet_order),
      x1=x1, x2=x2,
      y0=ytop + add0,
      y =ytop + add1,
      ylab=ytop + addlab,
      lab=lab_vec[prots],
      stringsAsFactors=FALSE
    )
  }
  
  lab_bh  <- setNames(paste0(label_prefix, stars(adj), "  BH p=", fmt_p(adj)), names(adj))
  lab_nom <- setNames(paste0(label_prefix, stars(pv),  "  p=",     fmt_p(pv)),  names(pv))
  
  list(
    nom = build_df(nom_prots, add0_nom, add1_nom, addlab_nom, lab_nom),
    bh  = if (is.null(adj)) empty else build_df(bh_prots, add0_bh, add1_bh, addlab_bh, lab_bh)
  )
}
# x positions: Mild=1, Severe=2, Critical=3
br_SM <- make_brackets_log2("Severe-Mild",     x1=1, x2=2, label_prefix="")
br_CM <- make_brackets_log2("Critical-Mild",   x1=1, x2=3, label_prefix="C−M: ",
                            add0_nom=0.20, add1_nom=0.25, addlab_nom=0.29,
                            add0_bh =0.28, add1_bh =0.33, addlab_bh =0.37)
br_CS <- make_brackets_log2("Critical-Severe", x1=2, x2=3, label_prefix="C−S: ",
                            add0_nom=0.35, add1_nom=0.40, addlab_nom=0.44,
                            add0_bh =0.43, add1_bh =0.48, addlab_bh =0.52)

# Caption
cap_raw <- paste0(
  "Each panel shows batch-corrected ELISA values on the log2 scale (ycorr) by clinical group. ",
  "Batch correction was performed by subtracting fitted center and plate contributions from y ~ group + center + plate + sex + age + vaccinated. ",
  "Adjusted group means (points) and 95% CIs (bars) are from emmeans (weights='proportional'). ",
  "Y-axis is log2, but tick labels are shown on the original concentration scale (2^log2). ",
  "Brackets show BH significance when BH p<0.05; otherwise raw p is shown only if nominally significant."
)
cap <- paste(strwrap(cap_raw, width=150), collapse="\n")

# Plot
p <- ggplot(df_plot, aes(x=group, y=y_log2)) +
  geom_violin(aes(fill=group), alpha=0.06, color="grey35",
              width=0.72, trim=FALSE, linewidth=0.25) +
  geom_jitter(aes(color=group), width=0.10, height=0,
              alpha=0.75, size=1.10) +
  
  geom_errorbar(data=mean_df,
                aes(x=group, y=mean, ymin=lo, ymax=hi),
                inherit.aes=FALSE,
                position=position_nudge(x=0.25),
                width=0.05, linewidth=0.60, color="black") +
  geom_point(data=mean_df,
             aes(x=group, y=mean, fill=group),
             inherit.aes=FALSE,
             position=position_nudge(x=0.25),
             shape=21, size=2.6, stroke=0.45, color="black") +
  
  # --- brackets (nominal first, then BH) ---
  # Severe-Mild
  geom_segment(data=br_SM$nom, aes(x=x1, xend=x2, y=y0, yend=y0), inherit.aes=FALSE, linewidth=0.7) +
  geom_segment(data=br_SM$nom, aes(x=x1, xend=x1, y=y0, yend=y),  inherit.aes=FALSE, linewidth=0.7) +
  geom_segment(data=br_SM$nom, aes(x=x2, xend=x2, y=y0, yend=y),  inherit.aes=FALSE, linewidth=0.7) +
  geom_text(   data=br_SM$nom, aes(x=(x1+x2)/2, y=ylab, label=lab), inherit.aes=FALSE, size=3.0) +
  
  geom_segment(data=br_SM$bh,  aes(x=x1, xend=x2, y=y0, yend=y0), inherit.aes=FALSE, linewidth=0.7) +
  geom_segment(data=br_SM$bh,  aes(x=x1, xend=x1, y=y0, yend=y),  inherit.aes=FALSE, linewidth=0.7) +
  geom_segment(data=br_SM$bh,  aes(x=x2, xend=x2, y=y0, yend=y),  inherit.aes=FALSE, linewidth=0.7) +
  geom_text(   data=br_SM$bh,  aes(x=(x1+x2)/2, y=ylab, label=lab), inherit.aes=FALSE, size=3.0) +
  
  # Critical-Mild
  geom_segment(data=br_CM$nom, aes(x=x1, xend=x2, y=y0, yend=y0), inherit.aes=FALSE, linewidth=0.7) +
  geom_segment(data=br_CM$nom, aes(x=x1, xend=x1, y=y0, yend=y),  inherit.aes=FALSE, linewidth=0.7) +
  geom_segment(data=br_CM$nom, aes(x=x2, xend=x2, y=y0, yend=y),  inherit.aes=FALSE, linewidth=0.7) +
  geom_text(   data=br_CM$nom, aes(x=(x1+x2)/2, y=ylab, label=lab), inherit.aes=FALSE, size=3.0) +
  
  geom_segment(data=br_CM$bh,  aes(x=x1, xend=x2, y=y0, yend=y0), inherit.aes=FALSE, linewidth=0.7) +
  geom_segment(data=br_CM$bh,  aes(x=x1, xend=x1, y=y0, yend=y),  inherit.aes=FALSE, linewidth=0.7) +
  geom_segment(data=br_CM$bh,  aes(x=x2, xend=x2, y=y0, yend=y),  inherit.aes=FALSE, linewidth=0.7) +
  geom_text(   data=br_CM$bh,  aes(x=(x1+x2)/2, y=ylab, label=lab), inherit.aes=FALSE, size=3.0) +
  
  # Critical-Severe
  geom_segment(data=br_CS$nom, aes(x=x1, xend=x2, y=y0, yend=y0), inherit.aes=FALSE, linewidth=0.7) +
  geom_segment(data=br_CS$nom, aes(x=x1, xend=x1, y=y0, yend=y),  inherit.aes=FALSE, linewidth=0.7) +
  geom_segment(data=br_CS$nom, aes(x=x2, xend=x2, y=y0, yend=y),  inherit.aes=FALSE, linewidth=0.7) +
  geom_text(   data=br_CS$nom, aes(x=(x1+x2)/2, y=ylab, label=lab), inherit.aes=FALSE, size=3.0) +
  
  geom_segment(data=br_CS$bh,  aes(x=x1, xend=x2, y=y0, yend=y0), inherit.aes=FALSE, linewidth=0.7) +
  geom_segment(data=br_CS$bh,  aes(x=x1, xend=x1, y=y0, yend=y),  inherit.aes=FALSE, linewidth=0.7) +
  geom_segment(data=br_CS$bh,  aes(x=x2, xend=x2, y=y0, yend=y),  inherit.aes=FALSE, linewidth=0.7) +
  geom_text(   data=br_CS$bh,  aes(x=(x1+x2)/2, y=ylab, label=lab), inherit.aes=FALSE, size=3.0) +
  
  facet_wrap(~ facet, scales="free_y") +
  scale_color_manual(values=cols, drop=FALSE) +
  scale_fill_manual(values=cols, drop=FALSE) +
  
  # show ORIGINAL SCALE labels (no scientific)
  scale_y_continuous(labels = function(z) fmt_num(2^z),
                     expand = expansion(mult = c(0.02, 0.14))) +
  
  labs(x=NULL, y="Batch-corrected value (log2; labels on original scale)", caption=cap) +
  coord_cartesian(clip="off") +
  theme_classic(base_size=12) +
  theme(
    legend.position="none",
    strip.background=element_rect(fill="grey92", color="black", linewidth=0.4),
    strip.text=element_text(face="bold", size=10),
    panel.border=element_rect(fill=NA, color="black", linewidth=0.5),
    axis.line=element_line(color="black", linewidth=0.4),
    axis.ticks=element_line(color="black", linewidth=0.4),
    axis.text.x=element_text(size=9),
    plot.margin=ggplot2::margin(10, 10, 10, 10),
    plot.caption=element_text(hjust=0, size=9, margin=ggplot2::margin(t=8))
  )

p
ggsave(
  filename = paste0(.res, "ELISA_groups_ycorr_log2Scale_originalLabels_facets_adjSexAgeVacc.png"),
  plot = p,
  width = 1000, height = 600, units = "px", dpi = 96
)








################################################################################################################
################################################################################################################
################################################################################################################
