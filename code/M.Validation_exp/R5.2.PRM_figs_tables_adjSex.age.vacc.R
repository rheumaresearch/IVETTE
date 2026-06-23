#################################################################################################################
#################################################################################################################
### R.For_drafts/R5.2.PRM_figs_tables_adjSex.age.vacc.R
### Project: jcalvet_202110_marato
#################################################################################################################
#################################################################################################################
  # Libraries, paths and options

rm(list=ls())
.res <- paste0(.res0, "M.Validation_exp/R5.2.PRM_figs_tables_adjSex.age.vacc/")
dir.create(.res, F, T)

#################################################################################################################
  # Merge differential expression results for PRM

rm(list=ls())
d1 <- readRDS(paste0(.res0, "M.Validation_exp/M5.3.Diff_expression_PRM_1st.round_MaxLQF_no.DIA_adj.sex.age.vaccine/dcomp.rds"))
d2 <- readRDS(paste0(.res0, "M.Validation_exp/M12.3.Diff_expression_PRM_2nd.round_MaxLQF_no.DIA_adj.sex.age.vaccine/dcomp.rds"))
d3 <- readRDS(paste0(.res0, "M.Validation_exp/M20.3.Diff_expression_PRM_3rd.round_MaxLQF_no.DIA_adj.sex.age.vaccine/dcomp.rds"))
d4 <- readRDS(paste0(.res0, "M.Validation_exp/M25.2.Diff_expression_PRM_4th.round_MaxLQF_no.DIA_adj.sex.age.vaccine/dcomp.rds"))


d1 <- d1[d1$protein.id!='iRT-Kit_WR_fusion', ]
d2 <- d2[d2$protein.id!='iRT-Kit_WR_fusion', ]
d3 <- d3[d3$protein.id!='iRT-Kit_WR_fusion', ]
d4 <- d4[d4$protein.id!='iRT-Kit_WR_fusion', ]

d <- rbind(d1, d2, d3, d4)

### Redo multitesting considering all proteins
dr <- data.frame(d, check.names=F)

# npv <- colnames(d)[regexpr("\\.adj\\.pv$", colnames(d)) > 0]
# comps <- gsub("\\.adj\\.pv$", "", npv)
# npv <- gsub("\\.adj\\.pv$", ".pv", npv)
# 
# 
# for (o in npv)
# {
#     dr[, gsub("\\.pv", ".adj.pv", o)] <- NA
#     dr[, gsub("\\.pv", ".adj.pv", o)] <- p.adjust(dr[, o], "BH")
# }

saveRDS(dr, paste0(.res, "dcomp.rds"))

################################################################################################################
  # Results for validated proteins in Severe.vs.Mild

rm(list=ls())
da <- readRDS(paste0(.res, "dcomp.rds"))
protSel <- readRDS(paste0(.dat,"2_PRM/dproteins.sel.prm_20230710.rds"))

# Get proteins tested for this comparison according to DIA results.
prots <- protSel[protSel$`Severe-Mild`=="X","protein.group"]


nvs <- c("protein.id", "gene.symbol", "gene.name", "chrom",
         "Severe-Mild.fc", "Severe-Mild.pv", "Severe-Mild.adj.pv")

d <- da[match(prots, da$protein.id), nvs]

# Remove non quantified in PRM or failed (NA)
d <- d[!is.na(d$protein.id),]

# Get the significant results
d <- d[d$`Severe-Mild.adj.pv`<0.1 & abs(d$`Severe-Mild.fc`) > 1.20,]

# Order rows
protsOrds <- c("P01011", "P02763", "P0DJI8", "P18428", "Q08830", "P0DJI9", "Q96PD5", "P02654",
               "P04196", "P13796","P02750", "P02765", "P02766", "P00738", "P29622", "P00740")

d <- d[match(protsOrds,d$protein.id),]
d <- d[!is.na(d$protein.id),]

# d$"Critical-Severe.fc" <- formatC(d$"Critical-Severe.fc", format='f', dig=2)
d$"Severe-Mild.fc" <- formatC(d$"Severe-Mild.fc", format='f', dig=2)
# d$"Hospitalized-Not.hospitalized.fc" <- formatC(d$"Hospitalized-Not.hospitalized.fc", format='f', dig=2)

# d$"Critical-Severe.pv" <- base::format.pval(d$"Critical-Severe.pv")
d$"Severe-Mild.pv" <- base::format.pval(d$"Severe-Mild.pv")
# d$"Hospitalized-Not.hospitalized.pv" <- formatC(d$"Hospitalized-Not.hospitalized.pv", format='f', dig=2)

# d$"Critical-Severe.adj.pv" <- base::format.pval(d$"Critical-Severe.adj.pv")
d$"Severe-Mild.adj.pv" <- base::format.pval(d$"Severe-Mild.adj.pv")
# d$"Hospitalized-Not.hospitalized.adj.pv" <- formatC(d$"Hospitalized-Not.hospitalized.adj.pv", format='f', dig=2)

write.csv(d,paste0(.res,"PRM_valProt.univ.csv"),row.names = F)
saveRDS(d, paste0(.res, "dprm.val.prm.SevereMild.rds"))

################################################################################################################
  # Results from DIA Severe  vs Mild

rm(list=ls())
da <- readRDS(paste0(.res0, "G.Discovery_DIA_FP_MaxLFQ/G2.Diff_expression/dcomp.rds"))
d.prm.ms <- readRDS(paste0(.res, "dprm.val.prm.SevereMild.rds"))

# prots <- c("P26038", "P01011", "P02763", "P0DJI8", "P18428", "Q08830", "P0DJI9", "Q96PD5", "P02654",
#            "P04196", "P13796", "P02765", "P02766", "P00738", "P29622", "P00740")
prots <- d.prm.ms$protein.id
nvs <- c("protein.group",
         "Severe-Mild.fc", "Severe-Mild.pv", "Severe-Mild.adj.pv")

d <- da[match(prots, da$protein.group), nvs]

# d$"Critical-Severe.fc" <- formatC(d$"Critical-Severe.fc", format='f', dig=2)
d$"Severe-Mild.fc" <- formatC(d$"Severe-Mild.fc", format='f', dig=2)

# d$"Critical-Severe.pv" <- base::format.pval(d$"Critical-Severe.pv")
d$"Severe-Mild.pv" <- base::format.pval(d$"Severe-Mild.pv")

# d$"Critical-Severe.adj.pv" <- base::format.pval(d$"Critical-Severe.adj.pv")
d$"Severe-Mild.adj.pv" <- base::format.pval(d$"Severe-Mild.adj.pv")


saveRDS(d, paste0(.res, "dprm.val.dia.prm.SevereMild.rds"))

################################################################################################################
  # Merge results from PRM and DIA

rm(list=ls())
d <- readRDS(paste0(.res, "dprm.val.prm.SevereMild.rds"))
dd <- readRDS(paste0(.res, "dprm.val.dia.prm.SevereMild.rds"))

colnames(d)[-c(1:4)] <- paste0("PRM.", colnames(d)[-c(1:4)])
colnames(dd)[-c(1)] <- paste0("DIA.", colnames(dd)[-c(1)])

d <- cbind(d, dd[match(d$protein.id, dd$protein.group), -1])
for (o in colnames(d))
    d[, o] <- as.character(d[, o])
d <- cbind(rownames(d), d)
d <- rbind(colnames(d), d)


p <- openPage(paste0(.res,  "PRM_DIA_validated_SevereMild.html"))
hwrite(paste0("<br><br>Results from PRM and dia for validated proteins<br><br>"), p)
hwrite(d, page=p, center=F, row.names=F, col.names=F,
       col.style=c("text-align:left", rep("text-align:center", ncol(d)-1)),
       col.width=c("400px", rep("150px", ncol(d)-1)))
closePage(p)

saveRDS(d,paste0(.res,"dprm.val.dia.prm.SevereMild.merge.rds"))


################################################################################################################
# Results for validated proteins in Critical.vs.Severe

rm(list=ls())
da <- readRDS(paste0(.res, "dcomp.rds"))
protSel <- readRDS(paste0(.dat,"2_PRM/dproteins.sel.prm_20230710.rds"))

# Get proteins tested for this comparison according to DIA results.
prots <- protSel[protSel$`Critical-Severe`=="X","protein.group"]

nvs <- c("protein.id", "gene.symbol", "gene.name", "chrom",
         "Critical-Severe.fc", "Critical-Severe.pv", "Critical-Severe.adj.pv")

d <- da[match(prots, da$protein.id), nvs]

# Remove non quantified in PRM or failed (NA)
d <- d[!is.na(d$protein.id),]
# d <- d[order(d$`Critical-Severe.pv`),]

# Get the significant results
d <- d[d$`Critical-Severe.pv`<0.05 & abs(d$`Critical-Severe.fc`) > 1.20,]


d$"Critical-Severe.fc" <- formatC(d$"Critical-Severe.fc", format='f', dig=2)
# d$"Severe-Mild.fc" <- formatC(d$"Severe-Mild.fc", format='f', dig=2)
# d$"Hospitalized-Not.hospitalized.fc" <- formatC(d$"Hospitalized-Not.hospitalized.fc", format='f', dig=2)

d$"Critical-Severe.pv" <- base::format.pval(d$"Critical-Severe.pv")
# d$"Severe-Mild.pv" <- base::format.pval(d$"Severe-Mild.pv")
# d$"Hospitalized-Not.hospitalized.pv" <- formatC(d$"Hospitalized-Not.hospitalized.pv", format='f', dig=2)

d$"Critical-Severe.adj.pv" <- base::format.pval(d$"Critical-Severe.adj.pv")
# d$"Severe-Mild.adj.pv" <- base::format.pval(d$"Severe-Mild.adj.pv")
# d$"Hospitalized-Not.hospitalized.adj.pv" <- formatC(d$"Hospitalized-Not.hospitalized.adj.pv", format='f', dig=2)

write.csv(d,paste0(.res,"PRM_valProt_CritSev.csv"),row.names = F)
saveRDS(d, paste0(.res, "dprm.val.prm.CritSevere.rds"))

################################################################################################################
# Results from DIA Critical vs Severe

rm(list=ls())
da <- readRDS(paste0(.res0, "G.Discovery_DIA_FP_MaxLFQ/G2.Diff_expression/dcomp.rds"))
d.prm.cs <- readRDS(paste0(.res, "dprm.val.prm.CritSevere.rds"))

prots <- d.prm.cs$protein.id
nvs <- c("protein.group",
         "Critical-Severe.fc", "Critical-Severe.pv", "Critical-Severe.adj.pv")
d <- da[match(prots, da$protein.group), nvs]

d$"Critical-Severe.fc" <- formatC(d$"Critical-Severe.fc", format='f', dig=2)
# d$"Severe-Mild.fc" <- formatC(d$"Severe-Mild.fc", format='f', dig=2)

d$"Critical-Severe.pv" <- base::format.pval(d$"Critical-Severe.pv")
# d$"Severe-Mild.pv" <- base::format.pval(d$"Severe-Mild.pv")

d$"Critical-Severe.adj.pv" <- base::format.pval(d$"Critical-Severe.adj.pv")
# d$"Severe-Mild.adj.pv" <- base::format.pval(d$"Severe-Mild.adj.pv")


saveRDS(d, paste0(.res, "dprm.val.dia.prm.CritSevere.rds"))

################################################################################################################
# Merge results from PRM and DIA

rm(list=ls())
d <- readRDS(paste0(.res, "dprm.val.prm.CritSevere.rds"))
dd <- readRDS(paste0(.res, "dprm.val.dia.prm.CritSevere.rds"))

colnames(d)[-c(1:4)] <- paste0("PRM.", colnames(d)[-c(1:4)])
colnames(dd)[-c(1)] <- paste0("DIA.", colnames(dd)[-c(1)])

d <- cbind(d, dd[match(d$protein.id, dd$protein.group), -1])
for (o in colnames(d))
  d[, o] <- as.character(d[, o])
d <- cbind(rownames(d), d)
d <- rbind(colnames(d), d)


p <- openPage(paste0(.res,  "PRM_DIA_validated_CriticalSevere.html"))
hwrite(paste0("<br><br>Results from PRM and dia for validated proteins<br><br>"), p)
hwrite(d, page=p, center=F, row.names=F, col.names=F,
       col.style=c("text-align:left", rep("text-align:center", ncol(d)-1)),
       col.width=c("400px", rep("150px", ncol(d)-1)))
closePage(p)

saveRDS(d,paste0(.res,"dprm.val.dia.prm.CritSevere.merge.rds"))


################################################################################################################
# Results for validated proteins in Asymptomatic

rm(list=ls())
da <- readRDS(paste0(.res, "dcomp.rds"))
protSel <- readRDS(paste0(.dat,"2_PRM/dproteins.sel.prm_20230710.rds"))

# Get proteins tested for this comparison according to DIA results.
prots <- protSel[protSel$`Mild-Asymptomatic`=="X","protein.group"]


nvs <- c("protein.id", "gene.symbol", "gene.name", "chrom",
         "Mild-Asymptomatic.fc", "Mild-Asymptomatic.pv", "Mild-Asymptomatic.adj.pv")

d <- da[match(prots, da$protein.id), nvs]

# Remove non quantified in PRM or failed (NA)
d <- d[!is.na(d$protein.id),]
# d <- d[order(d$`Mild-Asymptomatic.pv`),]

# Get the significant results
d <- d[d$`Mild-Asymptomatic.pv`<0.1 & abs(d$`Mild-Asymptomatic.fc`) > 1.20,]


d$"Mild-Asymptomatic.fc" <- formatC(d$"Mild-Asymptomatic.fc", format='f', dig=2)
# d$"Severe-Mild.fc" <- formatC(d$"Severe-Mild.fc", format='f', dig=2)
# d$"Hospitalized-Not.hospitalized.fc" <- formatC(d$"Hospitalized-Not.hospitalized.fc", format='f', dig=2)

d$"Mild-Asymptomatic.pv" <- base::format.pval(d$"Mild-Asymptomatic.pv")
# d$"Severe-Mild.pv" <- base::format.pval(d$"Severe-Mild.pv")
# d$"Hospitalized-Not.hospitalized.pv" <- formatC(d$"Hospitalized-Not.hospitalized.pv", format='f', dig=2)

d$"Mild-Asymptomatic.adj.pv" <- base::format.pval(d$"Mild-Asymptomatic.adj.pv")
# d$"Severe-Mild.adj.pv" <- base::format.pval(d$"Severe-Mild.adj.pv")
# d$"Hospitalized-Not.hospitalized.adj.pv" <- formatC(d$"Hospitalized-Not.hospitalized.adj.pv", format='f', dig=2)

write.csv(d,paste0(.res,"PRM_valProt.adj_MildAsym.csv"),row.names = F)
saveRDS(d, paste0(.res, "dprm.val.prm.MildAsym.rds"))

################################################################################################################
# Results from DIA Mild vs Asymptomatic

rm(list=ls())
da <- readRDS(paste0(.res0, "G.Discovery_DIA_FP_MaxLFQ/G2.Diff_expression/dcomp.rds"))
d.prm.cs <- readRDS(paste0(.res, "dprm.val.prm.MildAsym.rds"))

prots <- d.prm.cs$protein.id
nvs <- c("protein.group",
         "Mild-Asymptomatic.fc", "Mild-Asymptomatic.pv", "Mild-Asymptomatic.adj.pv")
d <- da[match(prots, da$protein.group), nvs]

d$"Mild-Asymptomatic.fc" <- formatC(d$"Mild-Asymptomatic.fc", format='f', dig=2)
# d$"Severe-Mild.fc" <- formatC(d$"Severe-Mild.fc", format='f', dig=2)

d$"Mild-Asymptomatic.pv" <- base::format.pval(d$"Mild-Asymptomatic.pv")
# d$"Severe-Mild.pv" <- base::format.pval(d$"Severe-Mild.pv")

d$"Mild-Asymptomatic.adj.pv" <- base::format.pval(d$"Mild-Asymptomatic.adj.pv")
# d$"Severe-Mild.adj.pv" <- base::format.pval(d$"Severe-Mild.adj.pv")


saveRDS(d, paste0(.res, "dprm.val.dia.prm.MildAsym.rds"))

################################################################################################################
# Merge results from PRM and DIA

rm(list=ls())
d <- readRDS(paste0(.res, "dprm.val.prm.MildAsym.rds"))
dd <- readRDS(paste0(.res, "dprm.val.dia.prm.MildAsym.rds"))

colnames(d)[-c(1:4)] <- paste0("PRM.", colnames(d)[-c(1:4)])
colnames(dd)[-c(1)] <- paste0("DIA.", colnames(dd)[-c(1)])

d <- cbind(d, dd[match(d$protein.id, dd$protein.group), -1])
for (o in colnames(d))
  d[, o] <- as.character(d[, o])
d <- cbind(rownames(d), d)
d <- rbind(colnames(d), d)


p <- openPage(paste0(.res,  "PRM_DIA_validated_MildAsym.html"))
hwrite(paste0("<br><br>Results from PRM and dia for validated proteins<br><br>"), p)
hwrite(d, page=p, center=F, row.names=F, col.names=F,
       col.style=c("text-align:left", rep("text-align:center", ncol(d)-1)),
       col.width=c("400px", rep("150px", ncol(d)-1)))
closePage(p)

saveRDS(d,paste0(.res,"dprm.val.dia.prm.MildAsym.merge.rds"))



################################################################################################################
# ADD Univariable AUCs to dataset
rm(list=ls())
d.sm <- readRDS(paste0(.res,"dprm.val.dia.prm.SevereMild.merge.rds"))[-1,-1]
d.cs <- readRDS(paste0(.res,"dprm.val.dia.prm.CritSevere.merge.rds"))[-1,-1]
d.ma <- readRDS(paste0(.res,"dprm.val.dia.prm.MildAsym.merge.rds"))[-1,-1]
auc.dia <- readRDS(paste0(.res0,"P.Multiv_predictor/P5.Predict_univ/lauc.univ.dia.rds"))
auc.prm <- readRDS(paste0(.res0,"P.Multiv_predictor/P5.Predict_univ/lauc.univ.prm.rds"))

r1 <- t(sapply(auc.dia, function(l) sapply(l, function(o) o$rc["auc"])))
r2 <- t(sapply(auc.prm, function(l) sapply(l, function(o) o$rc["auc"])))

r1 <- cbind(rownames(r1), r1)
r2 <- cbind(rownames(r2), r2)


# Severe vs Mild
s1 <- match(rownames(r1),d.sm$gene.symbol)
s2 <- match(rownames(r2),d.sm$gene.symbol)
res.sm <- cbind(d.sm,
                "DIA.AUC"=r1[,"sev.mild.auc"][d.sm$gene.symbol],
                "PRM.AUC"=r2[,"sev.mild.auc"][d.sm$gene.symbol])

v <- c("protein.id","gene.symbol","gene.name","chrom",
       grep("DIA\\.",names(res.sm),value = T),
       grep("PRM\\.",names(res.sm),value = T))
res.sm <- res.sm[,v]

# Critical vs Severe
s1 <- match(rownames(r1),d.cs$gene.symbol)
s2 <- match(rownames(r2),d.cs$gene.symbol)
res.cs <- cbind(d.cs,
                "DIA.AUC"=r1[,"sev.mild.auc"][d.cs$gene.symbol],
                "PRM.AUC"=r2[,"sev.mild.auc"][d.cs$gene.symbol])

v <- c("protein.id","gene.symbol","gene.name","chrom",
       grep("DIA\\.",names(res.cs),value = T),
       grep("PRM\\.",names(res.cs),value = T))
res.cs <- res.cs[,v]

# Mild vs Asymptomatic
s1 <- match(rownames(r1),d.ma$gene.symbol)
s2 <- match(rownames(r2),d.ma$gene.symbol)
res.ma <- cbind(d.ma,
                "DIA.AUC"=r1[,"sev.mild.auc"][d.ma$gene.symbol],
                "PRM.AUC"=r2[,"sev.mild.auc"][d.ma$gene.symbol])

v <- c("protein.id","gene.symbol","gene.name","chrom",
       grep("DIA\\.",names(res.ma),value = T),
       grep("PRM\\.",names(res.ma),value = T))
res.ma <- res.ma[,v]


saveRDS(res.sm, paste0(.res,"dprm.val.dia.prm.SevereMild.merge.withAUC.rds"))
saveRDS(res.cs, paste0(.res,"dprm.val.dia.prm.CriticalSevere.merge.withAUC.rds"))
saveRDS(res.ma, paste0(.res,"dprm.val.dia.prm.MildAsymp.merge.withAUC.rds"))

################################################################################################################
# Formatting tables for draft - Severe vs Mild
rm(list=ls())
res.sm <- readRDS(paste0(.res,"dprm.val.dia.prm.SevereMild.merge.withAUC.rds"))

tab <- res.sm
names(tab) <- c(
  "UniProt ID",
  "Gene symbol",
  "Gene name",
  "Chromosome location",
  "DIA FC",
  "DIA P-value",
  "DIA Adjusted P-value",
  "DIA AUC [95% CI]",
  "PRM FC",
  "PRM P-value",
  "PRM Adjusted P-value",
  "PRM AUC [95% CI]"
)

# Format FC values: 2 digits
tab[["DIA FC"]] <- sprintf("%.2f", as.numeric(tab[["DIA FC"]]))
tab[["PRM FC"]] <- sprintf("%.2f", as.numeric(tab[["PRM FC"]]))

# Format p-values: 4 digits, <0.0001 if smaller
tab[["DIA P-value"]] <- ifelse(
  as.numeric(tab[["DIA P-value"]]) < 0.0001,
  "<0.0001",
  sprintf("%.4f", as.numeric(tab[["DIA P-value"]]))
)

tab[["DIA Adjusted P-value"]] <- ifelse(
  as.numeric(tab[["DIA Adjusted P-value"]]) < 0.0001,
  "<0.0001",
  sprintf("%.4f", as.numeric(tab[["DIA Adjusted P-value"]]))
)

tab[["PRM P-value"]] <- ifelse(
  as.numeric(tab[["PRM P-value"]]) < 0.0001,
  "<0.0001",
  sprintf("%.4f", as.numeric(tab[["PRM P-value"]]))
)

tab[["PRM Adjusted P-value"]] <- ifelse(
  as.numeric(tab[["PRM Adjusted P-value"]]) < 0.0001,
  "<0.0001",
  sprintf("%.4f", as.numeric(tab[["PRM Adjusted P-value"]]))
)

# Format AUC values: 3 digits
dia_auc_num <- regmatches(
  tab[["DIA AUC [95% CI]"]],
  gregexpr("[0-9]+\\.?[0-9]*", tab[["DIA AUC [95% CI]"]])
)

tab[["DIA AUC [95% CI]"]] <- vapply(
  dia_auc_num,
  function(x) {
    x <- as.numeric(x)
    sprintf("%.3f [%.3f, %.3f]", x[1], x[2], x[3])
  },
  character(1)
)

prm_auc_num <- regmatches(
  tab[["PRM AUC [95% CI]"]],
  gregexpr("[0-9]+\\.?[0-9]*", tab[["PRM AUC [95% CI]"]])
)

tab[["PRM AUC [95% CI]"]] <- vapply(
  prm_auc_num,
  function(x) {
    x <- as.numeric(x)
    sprintf("%.3f [%.3f, %.3f]", x[1], x[2], x[3])
  },
  character(1)
)

# Gt table
gt_tab_sm <-
  gt(tab) |>
  cols_label(
    `DIA FC` = "FC",
    `DIA P-value` = "P-value",
    `DIA Adjusted P-value` = "Adjusted P-value",
    `DIA AUC [95% CI]` = "AUC [95% CI]",
    `PRM FC` = "FC",
    `PRM P-value` = "P-value",
    `PRM Adjusted P-value` = "Adjusted P-value",
    `PRM AUC [95% CI]` = "AUC [95% CI]"
  ) |>
  tab_spanner(
    label = "Discovery DIA",
    columns = c(
      `DIA FC`,
      `DIA P-value`,
      `DIA Adjusted P-value`,
      `DIA AUC [95% CI]`
    )
  ) |>
  tab_spanner(
    label = "Validation PRM",
    columns = c(
      `PRM FC`,
      `PRM P-value`,
      `PRM Adjusted P-value`,
      `PRM AUC [95% CI]`
    )
  ) |>
  cols_align(
    align = "left",
    columns = c(
      `UniProt ID`,
      `Gene symbol`,
      `Gene name`,
      `Chromosome location`
    )
  ) |>
  cols_align(
    align = "center",
    columns = c(
      `DIA FC`,
      `DIA P-value`,
      `DIA Adjusted P-value`,
      `DIA AUC [95% CI]`,
      `PRM FC`,
      `PRM P-value`,
      `PRM Adjusted P-value`,
      `PRM AUC [95% CI]`
    )
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      rows = 1:9,
      columns = c(
        `UniProt ID`,
        `Gene symbol`,
        `Gene name`,
        `Chromosome location`
      )
    )
  ) |>
  tab_options(
    table.font.size = px(12),
    column_labels.font.weight = "bold",
    data_row.padding = px(4)
  )

gtsave(
  data = gt_tab_sm,
  filename = paste0(.res,"validation_table_severe_mild.html")
)
gtsave(
  data = gt_tab_sm,
  filename = paste0(.res,"validation_table_severe_mild.png"),
  vwidth = 3000,
  vheight = 1600,
  zoom = 2
)

gt_tab_sm




################################################################################################################
# Formatting tables for draft - Suplmentary: Crit.Sever - Mild,Asympt
rm(list=ls())
res.cs <- readRDS(paste0(.res,"dprm.val.dia.prm.CriticalSevere.merge.withAUC.rds"))
res.ma <- readRDS(paste0(.res,"dprm.val.dia.prm.MildAsymp.merge.withAUC.rds"))

tab1 <- res.cs
tab2 <- res.ma
names(tab1) <- c(
  "UniProt ID",
  "Gene symbol",
  "Gene name",
  "Chromosome location",
  "DIA FC",
  "DIA P-value",
  "DIA Adjusted P-value",
  "DIA AUC [95% CI]",
  "PRM FC",
  "PRM P-value",
  "PRM Adjusted P-value",
  "PRM AUC [95% CI]"
)
names(tab2) <- names(tab1)
tab1$Comparison <- "Critical vs Severe"
tab2$Comparison <- "Mild vs Asymptomatic"

tab <- rbind(tab1, tab2)

# Format FC values: 2 digits
tab[["DIA FC"]] <- sprintf("%.2f", as.numeric(tab[["DIA FC"]]))
tab[["PRM FC"]] <- sprintf("%.2f", as.numeric(tab[["PRM FC"]]))

# Format p-values: 4 digits, <0.0001 if smaller
tab[["DIA P-value"]] <- ifelse(
  as.numeric(tab[["DIA P-value"]]) < 0.0001,
  "<0.0001",
  sprintf("%.4f", as.numeric(tab[["DIA P-value"]]))
)

tab[["DIA Adjusted P-value"]] <- ifelse(
  as.numeric(tab[["DIA Adjusted P-value"]]) < 0.0001,
  "<0.0001",
  sprintf("%.4f", as.numeric(tab[["DIA Adjusted P-value"]]))
)

tab[["PRM P-value"]] <- ifelse(
  as.numeric(tab[["PRM P-value"]]) < 0.0001,
  "<0.0001",
  sprintf("%.4f", as.numeric(tab[["PRM P-value"]]))
)

tab[["PRM Adjusted P-value"]] <- ifelse(
  as.numeric(tab[["PRM Adjusted P-value"]]) < 0.0001,
  "<0.0001",
  sprintf("%.4f", as.numeric(tab[["PRM Adjusted P-value"]]))
)

tab <- tab[tab$Comparison == "Critical vs Severe",]
gt_tab_supp <-
  gt(tab, groupname_col = "Comparison") |>
  cols_label(
    `DIA FC` = "FC",
    `DIA P-value` = "P-value",
    `DIA Adjusted P-value` = "Adjusted P-value",
    `DIA AUC [95% CI]` = "AUC [95% CI]",
    `PRM FC` = "FC",
    `PRM P-value` = "P-value",
    `PRM Adjusted P-value` = "Adjusted P-value",
    `PRM AUC [95% CI]` = "AUC [95% CI]"
  ) |>
  tab_spanner(
    label = "Discovery DIA",
    columns = c(
      `DIA FC`,
      `DIA P-value`,
      `DIA Adjusted P-value`,
      `DIA AUC [95% CI]`
    )
  ) |>
  tab_spanner(
    label = "Validation PRM",
    columns = c(
      `PRM FC`,
      `PRM P-value`,
      `PRM Adjusted P-value`,
      `PRM AUC [95% CI]`
    )
  ) |>
  cols_align(
    align = "center",
    columns = everything()
  ) |>
  tab_options(
    table.font.size = px(12),
    column_labels.font.weight = "bold",
    data_row.padding = px(4)
  )

gtsave(
  data = gt_tab_supp,
  filename = paste0(.res,"supplementary_critSever_mildAsympt_results.html")
)
gtsave(
  data = gt_tab_supp,
  filename = paste0(.res,"supplementary_critSever_mildAsympt_results.png"),
  vwidth = 3000,
  vheight = 1600,
  zoom = 2
)

gt_tab_supp


################################################################################################################
# Results for validated proteins ordinal models
#  Validation results according to ordinal trends

rm(list=ls())
dd <- readRDS(paste0(.res0, "G.Discovery_DIA_FP_MaxLFQ/G4.Diff_expression_ordinal/dcomp.rds"))
dp <- readRDS(paste0(.res0, "M.Validation_exp/M26.3.Diff_expression_PRM_all_ordinal_adjSexAgeVacc/dcomp.rds"))
nd <- readRDS(paste0(.res0, "M.Validation_exp/M26.3.Diff_expression_PRM_all_ordinal_adjSexAgeVacc/nd.prm.all.rds"))

# M.Validation_exp/M26.3.Diff_expression_PRM_all_ordinal_adjSexAgeVacc/

prots <- intersect(rownames(dd), rownames(dp))

ld <- dd[prots, "model.aic"]
lp  <- dp[prots, "model.aic"]

ds <- data.frame(prot=prots, dia=ld, prm=lp)


nrow(ds)
nrow(ds[ds[, 2]==ds[, 3],])
nrow(ds[ds[, 2]==ds[, 3] & ds[, 2]=='ordinal',])

dss <- ds[ds[, 2]==ds[, 3] & ds[, 2]=='ordinal', ]

protv <- dss[, 1]

nvs <- c("protein.group","lineal.fc", "lineal.pv", "lineal.adj.pv")
dds <- dd[protv, nvs]

nvs <- c("protein.id", "lineal.fc", "lineal.pv", "lineal.adj.pv")
dps <- dp[protv, nvs]

colnames(dds)[-1] <- paste0("DIA.", colnames(dds)[-1])
colnames(dps)[-1] <- paste0("PRM.", colnames(dps)[-1])

dr <- cbind(dps, dds[match(dps$protein.id, dds$protein.group), -1])
dr <- cbind(dr, fData(nd)[match(dr$protein.id, fData(nd)$protein.id),
                          c("gene.symbol", "gene.name", "chrom")])
dr <- dr[, c(1, ncol(dr)-2, ncol(dr)-1, ncol(dr), 2:(ncol(dr)-3))]

dr <- dr[order(abs(dr$DIA.lineal.fc), decreasing=T), ]

dr <- dr[, c(1:4, 8:10, 5:7)]

dr <- dr[dr$PRM.lineal.pv < 0.05, ]

dr$"DIA.lineal.fc" <- formatC(dr$"DIA.lineal.fc", format='f', dig=2)
dr$"PRM.lineal.fc" <- formatC(dr$"PRM.lineal.fc", format='f', dig=2)

dr$"DIA.lineal.pv" <- base::format.pval(dr$"DIA.lineal.pv")
dr$"PRM.lineal.pv" <- base::format.pval(dr$"PRM.lineal.pv")

dr$"DIA.lineal.adj.pv" <- base::format.pval(dr$"DIA.lineal.adj.pv")
dr$"PRM.lineal.adj.pv" <- base::format.pval(dr$"PRM.lineal.adj.pv")


dr <- rbind(colnames(dr), dr)

p <- openPage(paste0(.res,  "PRM_DIA_validated_ordinal_adjSexAgeVacc.html"))
hwrite(paste0("<br><br>Results from PRM and dia for validated proteins with linear trend<br><br>"), p)
hwrite(dr, page=p, center=F, row.names=F, col.names=F,
       col.style=c("text-align:left", rep("text-align:center", ncol(dr)-1)),
       col.width=c("150px", rep("150px", ncol(dr)-1)))
closePage(p)


table(dd[rownames(dp), "model.aic"])

## Gt format

tab.lineal <- dr[-1, ]

tab.lineal <- tab.lineal[, c(
  "protein.id",
  "gene.symbol",
  "gene.name",
  "chrom",
  "DIA.lineal.fc",
  "DIA.lineal.pv",
  "DIA.lineal.adj.pv",
  "PRM.lineal.fc",
  "PRM.lineal.pv",
  "PRM.lineal.adj.pv"
)]

names(tab.lineal) <- c(
  "UniProt ID",
  "Gene symbol",
  "Gene name",
  "Chromosome location",
  "DIA FC",
  "DIA P-value",
  "DIA Adjusted P-value",
  "PRM FC",
  "PRM P-value",
  "PRM Adjusted P-value"
)

tab.lineal[["DIA FC"]] <- sprintf("%.2f",as.numeric(tab.lineal[["DIA FC"]]))
tab.lineal[["PRM FC"]] <- sprintf("%.2f",as.numeric(tab.lineal[["PRM FC"]]))

tab.lineal[["DIA P-value"]] <- ifelse(
  as.numeric(tab.lineal[["DIA P-value"]]) < 0.0001,
  "<0.0001",
  sprintf("%.4f", as.numeric(tab.lineal[["DIA P-value"]]))
)

tab.lineal[["DIA Adjusted P-value"]] <- ifelse(
  as.numeric(tab.lineal[["DIA Adjusted P-value"]]) < 0.0001,
  "<0.0001",
  sprintf("%.4f", as.numeric(tab.lineal[["DIA Adjusted P-value"]]))
)

tab.lineal[["PRM P-value"]] <- ifelse(
  as.numeric(tab.lineal[["PRM P-value"]]) < 0.0001,
  "<0.0001",
  sprintf("%.4f", as.numeric(tab.lineal[["PRM P-value"]]))
)

tab.lineal[["PRM Adjusted P-value"]] <- ifelse(
  as.numeric(tab.lineal[["PRM Adjusted P-value"]]) < 0.0001,
  "<0.0001",
  sprintf("%.4f", as.numeric(tab.lineal[["PRM Adjusted P-value"]]))
)

gt_tab_lineal <-
  gt(tab.lineal) |>
  cols_label(
    `DIA FC` = "FC",
    `DIA P-value` = "P-value",
    `DIA Adjusted P-value` = "Adjusted P-value",
    `PRM FC` = "FC",
    `PRM P-value` = "P-value",
    `PRM Adjusted P-value` = "Adjusted P-value"
  ) |>
  tab_spanner(
    label = "Discovery DIA",
    columns = c(
      `DIA FC`,
      `DIA P-value`,
      `DIA Adjusted P-value`
    )
  ) |>
  tab_spanner(
    label = "Validation PRM",
    columns = c(
      `PRM FC`,
      `PRM P-value`,
      `PRM Adjusted P-value`
    )
  ) |>
  cols_align(
    align = "center",
    columns = everything()
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      rows = 2,
      columns = c(
        `UniProt ID`,
        `Gene symbol`,
        `Gene name`,
        `Chromosome location`
      )
    )
  ) |>
  tab_options(
    table.font.size = px(12),
    column_labels.font.weight = "bold",
    data_row.padding = px(4)
  )

gt_tab_lineal

gtsave(
  data = gt_tab_lineal,
  filename = paste0(.res,"ordinal_severity_linear_results.html")
)
gtsave(
  data = gt_tab_lineal,
  filename = paste0(.res,"ordinal_severity_linear_results.png")
)




################################################################################################################
  # SERPINA3 - PRM

rm(list=ls())
nd <- readRDS(paste0(.res0, "O.Concordance_DIA_PRM_complete/O1.Merge_PRM_data/nd.prm.all.rds"))
dc <- readRDS(paste0(.res, "dcomp.rds"))
source("code/functions/make_transparent.R")

nprot <- "P01011"
nx <- "serpina3"
xlab <- "SERPINA3"


d <- pData(nd)
d[, nx] <- exprs(nd)[nprot, ]


cols <- c("green3", "deepskyblue3", "orange3", "blueviolet", "red3")
names(cols) <- levels(d$group)

d <- d[d$group!='Control', ]
d$group <- droplevels(d$group)

o <- d[, nx]
x <- as.numeric(d$group) + runif(length(o), -0.25, 0.25)
names(x) <- rownames(d)

d <- d[!is.na(d[, nx]), ]
d$group <- droplevels(d$group)
cols <- cols[levels(d$group)]


meds <- as.numeric(dc[nprot, paste0(levels(d$group), ".log2mean")])
names(meds) <- levels(d$group)
ees <- as.numeric(dc[nprot, paste0(levels(d$group), ".log2se")])
names(ees) <- levels(d$group)

pdf(paste0(.res, nx, ".pdf"), width=11, height=7)
par(mar=c(5, 10, 5, 3))
plot(x, o, pch=16, col=makeTransparent(cols[as.numeric(d$group)], alpha=0.4), xlim=c(0.5, nlevels(d$group)+1),
         ylab="MaxLQF intensity (log2 scale)\n\n", xlab='', sub='', main=xlab,
         cex.lab=1.2, cex.axis=1.3, cex.sub=1.3, cex.main=1.5, axes=F, cex=2.3)
points(1:length(meds) + 0.45, meds, pch=5, col=cols, cex=2.5)
for(n in 1:length(meds))
        arrows(x0=n + 0.45, x1=n + 0.45, y0=meds[n] - ees[n], y1=meds[n] + ees[n], angle=90, length=0.2, code=3,
               col=cols[n], lwd=4)
axis(1, at=c(-100, 1:nlevels(d$group)+0.25, nlevels(d$group)+10),
     lab=c("", paste("\n", gsub("-", "\n", levels(d$group))), ""), cex.axis=1.3)
axis(2, at=seq(1, 1000, 1), labels=2^seq(1, 1000, 1), las=2)
#legend(x='bottomright', inset=0.01, pch=16, col=cols, legend=levels(d$group), cex=1.2, pt.cex=2)
dev.off()


################################################################################################################
  # ORM1 - PRM

rm(list=ls())
nd <- readRDS(paste0(.res0, "O.Concordance_DIA_PRM_complete/O1.Merge_PRM_data/nd.prm.all.rds"))
dc <- readRDS(paste0(.res, "dcomp.rds"))
source("code/functions/make_transparent.R")

nprot <- "P02763"
nx <- "orm1"
xlab <- "ORM1"


d <- pData(nd)
d[, nx] <- exprs(nd)[nprot, ]

cols <- c("green3", "deepskyblue3", "orange3", "blueviolet", "red3")
names(cols) <- levels(d$group)

d <- d[d$group!='Control', ]
d$group <- droplevels(d$group)

o <- d[, nx]
x <- as.numeric(d$group) + runif(length(o), -0.25, 0.25)
names(x) <- rownames(d)

d <- d[!is.na(d[, nx]), ]
d$group <- droplevels(d$group)
cols <- cols[levels(d$group)]


meds <- as.numeric(dc[nprot, paste0(levels(d$group), ".log2mean")])
names(meds) <- levels(d$group)
ees <- as.numeric(dc[nprot, paste0(levels(d$group), ".log2se")])
names(ees) <- levels(d$group)

pdf(paste0(.res, nx, ".pdf"), width=11, height=7)
par(mar=c(5, 10, 5, 3))
plot(x, o, pch=16, col=makeTransparent(cols[as.numeric(d$group)], alpha=0.4), xlim=c(0.5, nlevels(d$group)+1),
         ylab="MaxLQF intensity (log2 scale)\n\n", xlab='', sub='', main=xlab,
         cex.lab=1.2, cex.axis=1.3, cex.sub=1.3, cex.main=1.5, axes=F, cex=2.3)
points(1:length(meds) + 0.45, meds, pch=5, col=cols, cex=2.5)
for(n in 1:length(meds))
        arrows(x0=n + 0.45, x1=n + 0.45, y0=meds[n] - ees[n], y1=meds[n] + ees[n], angle=90, length=0.2, code=3,
               col=cols[n], lwd=4)
axis(1, at=c(-100, 1:nlevels(d$group)+0.25, nlevels(d$group)+10),
     lab=c("", paste("\n", gsub("-", "\n", levels(d$group))), ""), cex.axis=1.3)
axis(2, at=seq(1, 1000, 1), labels=2^seq(1, 1000, 1), las=2)
#legend(x='bottomright', inset=0.01, pch=16, col=cols, legend=levels(d$group), cex=1.2, pt.cex=2)
dev.off()

################################################################################################################
  # TTR - PRM

rm(list=ls())
nd <- readRDS(paste0(.res0, "O.Concordance_DIA_PRM_complete/O1.Merge_PRM_data/nd.prm.all.rds"))
dc <- readRDS(paste0(.res, "dcomp.rds"))
source("code/functions/make_transparent.R")

nprot <- "P02766"
nx <- "ttr"
xlab <- "TTR"


d <- pData(nd)
d[, nx] <- exprs(nd)[nprot, ]

cols <- c("green3", "deepskyblue3", "orange3", "blueviolet", "red3")
names(cols) <- levels(d$group)

d <- d[d$group!='Control', ]
d$group <- droplevels(d$group)

o <- d[, nx]
x <- as.numeric(d$group) + runif(length(o), -0.25, 0.25)
names(x) <- rownames(d)


d <- d[!is.na(d[, nx]), ]
d$group <- droplevels(d$group)
cols <- cols[levels(d$group)]


meds <- as.numeric(dc[nprot, paste0(levels(d$group), ".log2mean")])
names(meds) <- levels(d$group)
ees <- as.numeric(dc[nprot, paste0(levels(d$group), ".log2se")])
names(ees) <- levels(d$group)

pdf(paste0(.res, nx, ".pdf"), width=11, height=7)
par(mar=c(5, 10, 5, 3))
plot(x, o, pch=16, col=makeTransparent(cols[as.numeric(d$group)], alpha=0.4), xlim=c(0.5, nlevels(d$group)+1),
         ylab="MaxLQF intensity (log2 scale)\n\n", xlab='', sub='', main=xlab,
         cex.lab=1.2, cex.axis=1.3, cex.sub=1.3, cex.main=1.5, axes=F, cex=2.3)
points(1:length(meds) + 0.45, meds, pch=5, col=cols, cex=2.5)
for(n in 1:length(meds))
        arrows(x0=n + 0.45, x1=n + 0.45, y0=meds[n] - ees[n], y1=meds[n] + ees[n], angle=90, length=0.2, code=3,
               col=cols[n], lwd=4)
axis(1, at=c(-100, 1:nlevels(d$group)+0.25, nlevels(d$group)+10),
     lab=c("", paste("\n", gsub("-", "\n", levels(d$group))), ""), cex.axis=1.3)
axis(2, at=seq(1, 1000, 1), labels=2^seq(1, 1000, 1), las=2)
#legend(x='bottomright', inset=0.01, pch=16, col=cols, legend=levels(d$group), cex=1.2, pt.cex=2)
dev.off()

################################################################################################################
  # TNXB - PRM

rm(list=ls())
nd <- readRDS(paste0(.res0, "O.Concordance_DIA_PRM_complete/O1.Merge_PRM_data/nd.prm.all.rds"))
dc <- readRDS(paste0(.res, "dcomp.rds"))
source("code/functions/make_transparent.R")

nprot <- "P22105"
nx <- "tnxb"
xlab <- "TNXB"


d <- pData(nd)
d[, nx] <- exprs(nd)[nprot, ]

cols <- c("green3", "deepskyblue3", "orange3", "blueviolet", "red3")
names(cols) <- levels(d$group)

d <- d[d$group!='Control', ]
d$group <- droplevels(d$group)

o <- d[, nx]
x <- as.numeric(d$group) + runif(length(o), -0.25, 0.25)
names(x) <- rownames(d)

d <- d[!is.na(d[, nx]), ]
d$group <- droplevels(d$group)
cols <- cols[levels(d$group)]


meds <- as.numeric(dc[nprot, paste0(levels(d$group), ".log2mean")])
names(meds) <- levels(d$group)
ees <- as.numeric(dc[nprot, paste0(levels(d$group), ".log2se")])
names(ees) <- levels(d$group)

pdf(paste0(.res, nx, ".pdf"), width=11, height=7)
par(mar=c(5, 10, 5, 3))
plot(x, o, pch=16, col=makeTransparent(cols[as.numeric(d$group)], alpha=0.4), xlim=c(0.5, nlevels(d$group)+1),
         ylab="MaxLQF intensity (log2 scale)\n\n", xlab='', sub='', main=xlab,
         cex.lab=1.2, cex.axis=1.3, cex.sub=1.3, cex.main=1.5, axes=F, cex=2.3)
points(1:length(meds) + 0.45, meds, pch=5, col=cols, cex=2.5)
for(n in 1:length(meds))
        arrows(x0=n + 0.45, x1=n + 0.45, y0=meds[n] - ees[n], y1=meds[n] + ees[n], angle=90, length=0.2, code=3,
               col=cols[n], lwd=4)
axis(1, at=c(-100, 1:nlevels(d$group)+0.25, nlevels(d$group)+10),
     lab=c("", paste("\n", gsub("-", "\n", levels(d$group))), ""), cex.axis=1.3)
axis(2, at=seq(1, 1000, 1), labels=2^seq(1, 1000, 1), las=2)
#legend(x='bottomright', inset=0.01, pch=16, col=cols, legend=levels(d$group), cex=1.2, pt.cex=2)
dev.off()

################################################################################################################
  # APOC1 - PRM

rm(list=ls())
nd <- readRDS(paste0(.res0, "O.Concordance_DIA_PRM_complete/O1.Merge_PRM_data/nd.prm.all.rds"))
dc <- readRDS(paste0(.res, "dcomp.rds"))
source("code/functions/make_transparent.R")

nprot <- "P02654"
nx <- "apoc1"
xlab <- "APOC1"


d <- pData(nd)
d[, nx] <- exprs(nd)[nprot, ]


cols <- c("green3", "deepskyblue3", "orange3", "blueviolet", "red3")
names(cols) <- levels(d$group)

d <- d[d$group!='Control', ]
d$group <- droplevels(d$group)

o <- d[, nx]
x <- as.numeric(d$group) + runif(length(o), -0.25, 0.25)
names(x) <- rownames(d)


d <- d[!is.na(d[, nx]), ]
d$group <- droplevels(d$group)
cols <- cols[levels(d$group)]


meds <- as.numeric(dc[nprot, paste0(levels(d$group), ".log2mean")])
names(meds) <- levels(d$group)
ees <- as.numeric(dc[nprot, paste0(levels(d$group), ".log2se")])
names(ees) <- levels(d$group)

pdf(paste0(.res, nx, ".pdf"), width=11, height=7)
par(mar=c(5, 10, 5, 3))
plot(x, o, pch=16, col=makeTransparent(cols[as.numeric(d$group)], alpha=0.4), xlim=c(0.5, nlevels(d$group)+1),
         ylab="MaxLQF intensity (log2 scale)\n\n", xlab='', sub='', main=xlab,
         cex.lab=1.2, cex.axis=1.3, cex.sub=1.3, cex.main=1.5, axes=F, cex=2.3)
points(1:length(meds) + 0.45, meds, pch=5, col=cols, cex=2.5)
for(n in 1:length(meds))
        arrows(x0=n + 0.45, x1=n + 0.45, y0=meds[n] - ees[n], y1=meds[n] + ees[n], angle=90, length=0.2, code=3,
               col=cols[n], lwd=4)
axis(1, at=c(-100, 1:nlevels(d$group)+0.25, nlevels(d$group)+10),
     lab=c("", paste("\n", gsub("-", "\n", levels(d$group))), ""), cex.axis=1.3)
axis(2, at=seq(1, 1000, 1), labels=2^seq(1, 1000, 1), las=2)
#legend(x='bottomright', inset=0.01, pch=16, col=cols, legend=levels(d$group), cex=1.2, pt.cex=2)
dev.off()

################################################################################################################
  # SAA1 - PRM

rm(list=ls())
nd <- readRDS(paste0(.res0, "O.Concordance_DIA_PRM_complete/O1.Merge_PRM_data/nd.prm.all.rds"))
dc <- readRDS(paste0(.res, "dcomp.rds"))
source("code/functions/make_transparent.R")

nprot <- "P0DJI8"
nx <- "saa1"
xlab <- "SAA1"


d <- pData(nd)
d[, nx] <- exprs(nd)[nprot, ]

cols <- c("green3", "deepskyblue3", "orange3", "blueviolet", "red3")
names(cols) <- levels(d$group)

d <- d[d$group!='Control', ]
d$group <- droplevels(d$group)

o <- d[, nx]
x <- as.numeric(d$group) + runif(length(o), -0.25, 0.25)
names(x) <- rownames(d)

d <- d[!is.na(d[, nx]), ]
d$group <- droplevels(d$group)
cols <- cols[levels(d$group)]


meds <- as.numeric(dc[nprot, paste0(levels(d$group), ".log2mean")])
names(meds) <- levels(d$group)
ees <- as.numeric(dc[nprot, paste0(levels(d$group), ".log2se")])
names(ees) <- levels(d$group)

pdf(paste0(.res, nx, ".pdf"), width=11, height=7)
par(mar=c(5, 10, 5, 3))
plot(x, o, pch=16, col=makeTransparent(cols[as.numeric(d$group)], alpha=0.4), xlim=c(0.5, nlevels(d$group)+1),
         ylab="MaxLQF intensity (log2 scale)\n\n", xlab='', sub='', main=xlab,
         cex.lab=1.2, cex.axis=1.3, cex.sub=1.3, cex.main=1.5, axes=F, cex=2.3)
points(1:length(meds) + 0.45, meds, pch=5, col=cols, cex=2.5)
for(n in 1:length(meds))
        arrows(x0=n + 0.45, x1=n + 0.45, y0=meds[n] - ees[n], y1=meds[n] + ees[n], angle=90, length=0.2, code=3,
               col=cols[n], lwd=4)
axis(1, at=c(-100, 1:nlevels(d$group)+0.25, nlevels(d$group)+10),
     lab=c("", paste("\n", gsub("-", "\n", levels(d$group))), ""), cex.axis=1.3)
axis(2, at=seq(1, 1000, 1), labels=2^seq(1, 1000, 1), las=2)
#legend(x='bottomright', inset=0.01, pch=16, col=cols, legend=levels(d$group), cex=1.2, pt.cex=2)
dev.off()

################################################################################################################
  # Specificity for PRM proteins -> controls

rm(list=ls())
nd <- readRDS(paste0(.res0, "O.Concordance_DIA_PRM_complete/O1.Merge_PRM_data/nd.prm.all.rds"))
dc <- readRDS(paste0(.res, "dcomp.rds"))
source("code/functions/make_transparent.R")

nd <- nd[regexpr("iRT-Kit_WR_fusion", featureNames(nd)) < 0, ]
                 

nd <- nd[, nd$group%in%c("Control", "Asymptomatic")]
nd$group <- droplevels(nd$group)

r <- do.call(c, mclapply(1:nrow(nd), function(j)
{
    x <- exprs(nd)[j, ]
    if (median(x[nd$group=='Asymptomatic'], na.rm=T) > median(x[nd$group=='Control'], na.rm=T))
    {
        th <- quantile(x[nd$group=='Control'], 0.90, na.rm=T)
        tab <- table(nd$group, factor(as.numeric(!is.na(x) & x > th), levels=0:1))
        r <- tab[2, 2]/sum(tab[2,  ])
    }else
    {
        th <- quantile(x[nd$group=='Control'], 0.10, na.rm=T)
        tab <- table(nd$group, factor(as.numeric(!is.na(x) & x < th), levels=0:1))
        r <- tab[2, 2]/sum(tab[2,  ])
    }
    r
}, mc.cores=12))
names(r) <- fData(nd)$gene.symbol
sort(r)
summary(r) # sensitivitis around 30% for correctly detect asymptomatics using a threshold
          # that ensures a specificity of 90% for controls.


################################################################################################################
################################################################################################################
################################################################################################################




