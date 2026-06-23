#################################################################################################################
#################################################################################################################
### G.Discovery_DIA_FP_MaxLFQ/G5.3.ROAST_GSA_summary.R
### Project: jcalvet_202110_marato
#################################################################################################################
#################################################################################################################
  # Libraries, paths and options

rm(list=ls())
.res <- paste0(.res0, "G.Discovery_DIA_FP_MaxLFQ/G5.3.ROAST_GSA_summary/")
dir.create(.res, F, T)

################################################################################################################
  # Differential expression results - including ordinal models

rm(list=ls())
d1 <- readRDS(paste0(.res0, "G.Discovery_DIA_FP_MaxLFQ/G2.Diff_expression/dcomp.rds"))
d2 <- readRDS(paste0(.res0, "G.Discovery_DIA_FP_MaxLFQ/G4.Diff_expression_ordinal/dcomp.rds"))

d <- cbind(d1[, c("protein.group", "protein.name", "gene.symbol", "model.aic",
                  colnames(d1)[regexpr("\\.fc$", colnames(d1)) > 0 |
                               regexpr("\\.pv$", colnames(d1)) > 0 |
                               regexpr("\\.adj.pv$", colnames(d1)) > 0])],
           d2[rownames(d1),
              c(colnames(d2)[regexpr("\\.fc$", colnames(d2)) > 0 |
                               regexpr("\\.pv$", colnames(d2)) > 0 |
                               regexpr("\\.adj.pv$", colnames(d2)) > 0])])

saveRDS(d, paste0(.res, "dcomp.all.rds"))

################################################################################################################
  # Functions differentially enriched in each comparison

rm(list=ls())
lg <- readRDS(paste0(.res0, "G.Discovery_DIA_FP_MaxLFQ/G5.2.ROAST_GSA_linear/lrgsa.maxmean.all.rds"))
dc <- readRDS(paste0(.res, "dcomp.all.rds"))


lgs <- sapply(lg, function(l)
             sapply(l, function(o) rownames(o$res[o$res[, "pval"] < 0.05, , drop=F]), simplify=F), simplify=F)

gs <- unique(as.character(unlist(unlist(lgs))))
comps <- names(lgs)
dr <- matrix(NA, ncol=length(comps), nrow=length(gs))
rownames(dr) <- gs
colnames(dr) <- comps
for (i in 1:nrow(dr))
{
    for (k in 1:ncol(dr))
    {
        dr[i, k] <- c("", "X")[rownames(dr)[i]%in%unlist(lgs[colnames(dr)[k]]) + 1]
    }
}

drc <- dr
drc <- cbind(genesets=rownames(drc), drc)
write.table(drc, paste0(.res, "dgenesets.sign.comp.csv"), sep=';', row.names=F, quote=F)

################################################################################################################
  # Genes in significant genesets

rm(list=ls())
lg <- readRDS(paste0(.res0, "G.Discovery_DIA_FP_MaxLFQ/G5.2.ROAST_GSA_linear/lrgsa.maxmean.all.rds"))
dc <- readRDS(paste0(.res, "dcomp.all.rds"))
nd <- readRDS(paste0(.res0, "P.Multiv_predictor/P1.Data_setup/nd.dia.corr.rds"))


### Proteins validated by PRM

ps <- c("P01011", "P02763", "P0DJI8",
        "P18428", "Q08830",
        "P0DJI9",
        "Q96PD5", "P02654", "P04196",
        "P13796", "P02765","P02766", "P00738", "P29622", "P00740",
        "P26038",
        "P22105", "Q6UXB8", "O75636", "P29401", "P15169", "P07357")
psg <- fData(nd)[match(ps, fData(nd)$protein.group), "gene.symbol"]


### Comparisons to summarize

comps <- unique(gsub("\\.pv", "", gsub("\\.adj", "", colnames(dc)[regexpr("\\.pv$", colnames(dc)) > 0])))
comps <- comps[comps!='Hospitalized-Not.hospitalized']


### Selected genesets - different thresholds for each collection

lgs <- sapply(lg, function(l)
    sapply(names(l), function(o, l)
    {
        if (o == "GOBP_v5") th <- 0.01 else th <-0.05
        rownames(l[[o]]$res[l[[o]]$res[, "pval"] < th, , drop=F])
    }, l, simplify=F), simplify=F)
r <- t(sapply(lgs, function(l) sapply(l, length)))
r

lgs2 <- sapply(1:length(lgs), function(j)
{
    l <- lgs[[j]]
    lr <- sapply(1:length(l), function(k, l)
    {
        if (length(l[[k]]) > 0)
        {
            paste0(names(l)[k], " - ", unlist(l[k]))
        }else
        {
            NULL
        }
    }, l, simplify=F)
    names(lr) <- names(l)
    lr[sapply(lr, length) > 0]
}, simplify=F)
names(lgs2) <- names(lgs)

lgs <- lgs2


gs <- unique(unlist(unlist(lgs)))

dr <- matrix(NA, ncol=length(comps), nrow=length(gs))
rownames(dr) <- gs
colnames(dr) <- comps
for (i in 1:nrow(dr))
{
    for (k in 1:ncol(dr))
    {
        dr[i, k] <- c("", "X")[rownames(dr)[i]%in%unlist(lgs[colnames(dr)[k]]) + 1]
    }
}


### Genes significant and involved in selected genesets

lgg <- sapply(lg[[1]], function(o) o$index)

ng <- unique(unlist(sapply(comps, function(o) dc[dc[, paste0(o, ".pv")] < 0.01, "gene.symbol"])))
ng <- ng[ng!='NA' & !is.na(ng)]
#ng <- ng[ng%in%unlist(unlist(lgg))]

lggs <- do.call(c, lgg)
names(lggs) <- gsub("v5\\.", "v5 - ", names(lggs))
lggs <- lggs[gs]
lggs <- sapply(lggs, function(o) o[o%in%ng])
lggs <- sapply(lggs, function(o) c(o[o%in%psg], " ||| ", o[!o%in%psg]))

dd <- do.call(rbind, strsplit(rownames(dr), split=' \\- '))
colnames(dd) <- c("Collection", "gset")

dr <- data.frame(dd, dr, stringsAsFactors=F)
dr$genes.in <- sapply(lggs[rownames(dr)], paste0, collapse=', ')
drc <- dr
write.table(drc, paste0(.res, "dgenesets.sign.comp.csv"), sep=';', row.names=F, quote=F)


dr <- dr[order(dr$Collection, dr$gset, substring(dr$genes.in, 1, 10)), ]
drc <- dr
write.table(drc, paste0(.res, "dgenesets.sign.comp_2.csv"), sep=';', row.names=F, quote=F)


### Genesets of genes validated by PRM

lgg <- sapply(lg[[1]], function(o) o$index)

ng <- psg

lggs <- do.call(c, lgg)
names(lggs) <- gsub("v5\\.", "v5 - ", names(lggs))
lggs <- lggs[gs]
lggs <- sapply(lggs, function(o) o[o%in%ng])
lggs <- sapply(lggs, function(o) c(o[o%in%psg]))
lggs <- lggs[sapply(lggs, length) > 0]

dr <- data.frame(gset=rep(names(lggs), sapply(lggs, length)),
                 genes=unlist(lggs), stringsAsFactors=F)
dr <- dr[order(dr$gset), ]

drc <- dr
write.table(drc, paste0(.res, "dgenesets.PRM.val.comp.csv"), sep=';', row.names=F, quote=F)

tab <- as.matrix(table(dr[, 1]))
tab <- tab[order(tab, decreasing=T), , drop=F]


dr <- dr[order(dr$Collection, dr$gset, substring(dr$genes.in, 1, 10)), ]
drc <- dr
write.table(drc, paste0(.res, "dgenesets.sign.comp_2.csv"), sep=';', row.names=F, quote=F)



################################################################################################################
################################################################################################################
################################################################################################################


