#################################################################################################################
#################################################################################################################
### G.Discovery_DIA_FP_MaxLFQ/G5.2.ROAST_GSA_linear.R
### Project: jcalvet_202110_marato
#################################################################################################################
#################################################################################################################
  # Libraries, paths and options

rm(list=ls())
.res <- paste0(.res0, "G.Discovery_DIA_FP_MaxLFQ/G5.2.ROAST_GSA_linear/")
dir.create(.res, F, T)

#################################################################################################################
  # Prepare data

rm(list=ls())
nd <- readRDS(paste0(.res0, "G.Discovery_DIA_FP_MaxLFQ/G5.ROAST_GSA/nd.eof.roastgsa.rds"))
source("code/functions/make.contrasts.R")


nd$groupn <- as.numeric(nd$group) - 1

form <- "y ~ groupn + quant.batch"
form <- gsub("y ", "", form)

cont.mat <- "groupn"
saveRDS(cont.mat, paste0(.res, "cont.mat.rds"))

saveRDS(nd, paste0(.res, "nd.eof.roastgsa.rds"))
saveRDS(form, paste0(.res, "form.rds"))

#################################################################################################################
  # ROAST-GSA

rm(list=ls())
nd <- readRDS(paste0(.res, "nd.eof.roastgsa.rds"))
form <- readRDS(paste0(.res, "form.rds"))
cont.mat <- readRDS(paste0(.res, "cont.mat.rds"))
source("code/functions/make.contrasts.R")


### Genesets

gspath <- paste0(.dat, "20240628_genesets_irb/")
gsets <- list.files(gspath)
gsets <- gsets[regexpr("\\.gmt", gsets) > 0]
gsets <- gsub("\\.gmt", "", gsets)


### ROAST-GSA

lr <- list()
lrr <- list()
for(g in gsets)
{
    print(g)
    lrr[[g]] <-
        roastgsa(exprs(nd), form = form, covar = pData(nd),
                 contrast = cont.mat,
                 gspath = gspath,
                 gsetsel = g, nrot = 5000, mccores = 12, set.statistic = "maxmean",
                 normalizeScores=TRUE)
}
lr[[1]] <- lrr
names(lr) <- "lineal"

saveRDS(lr, file=paste0(.res, "lrgsa.maxmean.rds"))

#################################################################################################################
  # Finding effective signature size

rm(list=ls())
lrg <- readRDS(file=paste0(.res, "lrgsa.maxmean.rds"))
nd <- readRDS(paste0(.res, "nd.eof.roastgsa.rds"))

lr <- list()
for(o in names(lrg))
{
    lgs <- lrg[[o]][[1]]
    lr[[o]] <- varrotrand(lgs, exprs(nd),
                          testedsizes = c(3:4, seq(5, 50, by=5), seq(60, 200, by=10)),
                          nrep = 200)
}

saveRDS(lr, paste0(.res, "lsd.rot.Rdata"))

#################################################################################################################
 # List of genes for genesets

rm(list=ls())
lrg <- readRDS(file=paste0(.res, "lrgsa.maxmean.rds"))
lc <- readRDS(paste0(.res0, "G.Discovery_DIA_FP_MaxLFQ/G4.Diff_expression_ordinal/lcomp.rds"))

ngs <- "gene.symbol"
rdir <- paste0(.res, "gene_lists")

dir.create(paste0(rdir), F)

for(o in names(lrg))
{
    d <- data.frame(lc[[o]], check.names=F)
    d <- d[!is.na(d$gene.symbol) & d$gene.symbol!='NA', ]
    d <- d[!is.na(d[, paste0(o, ".pv")]), ]
    nv2 <- paste0(o, c(".fc", ".pv", ".adj.pv"))
    rownames(d) <- d$gene.symbol
    lg <- lrg[[o]]
    dir.create(paste0(rdir, "/", o), F)
    for (g in names(lg))
    {
        res <- lg[[g]]
        ndir <- paste0(rdir, "/", o, "/", g, "/")
        dir.create(ndir, F)
        lr <- mclapply(1:min(200,nrow(res$res)), function(k)
        {
            path <- rownames(res$res)[k]
            gs <- res$index[[path]]
            ps <- rownames(d)[d[,ngs] %in%gs]
            drr <- d[rownames(d)%in%ps, ]
            modt <- res$stats[order(-res$stats)]
            rnk <- 1:length(modt)
            names(rnk) <- names(modt)
            index <- sapply(res$index, function(o) rnk[o], simplify=F)
            index2 <- unlist(index)
            names(index2) <- do.call(c, sapply(index, function(o) names(o)))
            index3 <- unique(index2)
            names(index3) <- unique(names(index2))
            rank.pos <- index3[match(drr[,ngs], names(index3))]
            drr <- data.frame(rank.pos = rank.pos,
                              rank.neg = length(modt) - rank.pos + 1, drr,
                              stringsAsFactors=F, check.names=F)
            drr <- drr[, regexpr("\\.log2se$", colnames(drr)) < 0]
            drr <- drr[order(drr$rank.pos), ]
            drr$"plot<br>(only if called significant)" <- "view"
            nv2 <- paste0(o, c(".fc", ".pv", ".adj.pv"))
            nv3 <- "plot<br>(only if called significant)"
            nv4 <- colnames(drr)[regexpr("\\.log2mean$", colnames(drr)) > 0]
            nv1 <- c("rank.pos", "rank.neg", "gene.symbol", "protein.group", "protein.name",
                     "global.mean", "global.se")
            nvs <- c(nv1, nv2, nv3, nv4)
            drr <- drr[, nvs]
            nvs <- colnames(drr)[regexpr("fc", colnames(drr)) > 0 |
                                 regexpr("log2mean", colnames(drr)) > 0 |
                                 regexpr("global", colnames(drr)) > 0]
            for (oo in nvs) drr[, oo] <- round(drr[, oo], 3)
            nvs <- colnames(drr)[regexpr("pv", colnames(drr)) > 0]
            for (oo in nvs) drr[, oo] <- base::format.pval((drr[, oo]))
            drr$rank.pos[!rownames(drr)%in%ps] <- ""
            drr$rank.neg[!rownames(drr)%in%ps] <- ""
            links <- vector('list',length=ncol(drr))
            names(links) <- colnames(drr)
            plots <- links
            links$protein.group <- c(NA,
                                     paste0("https://www.uniprot.org/uniprot/",
                                            drr$protein.group))
            links$gene.symbol<- c(NA,
                             paste("http://www.ncbi.nlm.nih.gov/gene/?term=",
                                   drr$gene.symbol, sep=""))
            links$"plot<br>(only if called significant)" <-
                plots$plot <- c(NA, paste0("../../../../", "G2.Diff_expression/Stripcharts_zoom/",
                                           drr$protein.group, ".png"))
            drr <- rbind(colnames(drr), drr)
            is.plot <- file.exists(paste0(ndir,  links$"plot<br>(only if called significant)"))
            drr$"plot<br>(only if called significant)"[!is.plot][-1] <- ""
            links$"plot<br>(only if called significant)"[!is.plot] <- NA
            p <- openPage(paste0(ndir,
                                 gsub("", "", gsub("[[:punct:]]", "_", path)),
                                 "_genes.html"))
            hwrite(paste0("<br>", path, "<br><br>"), p);
            hwrite(drr, page=p, center=F, row.names=F, col.names=F,
                   col.style=c("text-align:center", rep("text-align:center",
                                                        ncol(drr)-1)),
                   col.width=c(100, rep(125, ncol(drr)-1)),
                   col.links=links)
            closePage(p)
        }, mc.cores=12)
    }
}


#################################################################################################################
  # Plots

rm(list=ls())
lrg <- readRDS(file=paste0(.res, "lrgsa.maxmean.rds"))
lse <- readRDS(paste0(.res, "lsd.rot.Rdata"))
nd <- readRDS(paste0(.res, "nd.eof.roastgsa.rds"))
form <- readRDS(paste0(.res, "form.rds"))
cont.mat <- readRDS( paste0(.res, "cont.mat.rds"))
source("code/functions/write.html.mod.R")
load("code/functions/write.html.mod.Rdata")

gdir <- paste0(.res, "gene_lists/")
rdir <- paste0(.res, "rgsa_results/")
dir.create(paste0(rdir), F)

cols <- c("deepskyblue3", "orange3", "blueviolet", "red3")
gr <- pData(nd)$group

lr <- mclapply(names(lrg), function(o)
{
    lgs <- lrg[[o]]
    print(o)
    for (g in names(lgs))
    {
        print(g)
        lg <- lgs[[g]]
        fs <- paste0("../gene_lists/", o, "/", g, "/",
                     gsub("[[:punct:]]", "_", rownames(lg$res)), "_genes.html")
        htmlrgsa(lg,
                 htmlpath = rdir,
                 htmlname = paste0(o, "_", g, ".html"),
                 plotpath =paste0("plots_", o, "_", g, "/"),
                 plotstats = TRUE, plotgsea = TRUE, indheatmap = TRUE, ploteffsize = TRUE,
                 y = exprs(nd),
                 whplots = rownames(lg$res)[1:min(200, nrow(lg$res))],
                 tit = paste0("ROAST-GSA results for ", o, " comparison - ", g, "<br>maxmean statistic"),
                 margins = c(5,16),
                 geneDEhtmlfiles  = fs,
                 varrot =  lse[[o]],
                 intvar = "group",
                 adj.var = c("quant.batchquant.batch.2", "quant.batchquant.batch.2", "quant.batchquant.batch.4"),
                 mycol = c("white", cols),
                 sorttable=sorttable, dragtable=dragtable,
                 typeheatmap = "ggplot2")
    }
}, mc.cores=3)

################################################################################################################
  # Heatmaps by genesets collections

rm(list=ls())
lrg <- readRDS(file=paste0(.res, "lrgsa.maxmean.rds"))
nd <- readRDS(paste0(.res, "nd.eof.roastgsa.rds"))

cols <- c("deepskyblue3", "orange3", "blueviolet", "red3")

dir.create(paste0(.res, "heatmaps_gsets"), F)
lr <- mclapply(names(lrg), function(o)
{
    lgs <- lrg[[o]]
    print(o)
    for (g in names(lgs))
    {
        print(g)
        lg <- lgs[[g]]
        png(paste0(.res, "heatmaps_gsets/", o, "_", g, ".png"), width = 1000, height=1500)
        aux <- heatmaprgsa_hm(lg, exprs(nd),  intvar = "group",
                           adj.var = c("quant.batchquant.batch.2", "quant.batchquant.batch.2", "quant.batchquant.batch.4"),
                           whplot = rownames(lg$res)[1:min(100, nrow(lg$res))],
                           toplot = TRUE, pathwaylevel = TRUE, fdrkey = TRUE,
                           mycol = c("white", cols),
                           dendrogram =  "n", col=  bluered(100), trace='none',
                           margins=c(8,20), notecol='black', notecex=1, keysize=0.9, cexCol=1.5,
                           Rowv = FALSE, Colv = FALSE, las=2)
        dev.off()
    }
}, mc.cores=3)

#################################################################################################################
  # HTML summary - add previous results

rm(list=ls())
lg1 <- readRDS(paste0(.res, "lrgsa.maxmean.rds"))
lg2 <- readRDS(paste0(.res0, "G.Discovery_DIA_FP_MaxLFQ/G5.ROAST_GSA/lrgsa.maxmean.rds"))
source("code/functions/write.html.mod.R")
load("code/functions/write.html.mod.Rdata")

lg <- c(lg2, lg1)

saveRDS(lg, paste0(.res, "lrgsa.maxmean.all.rds"))


d <- data.frame(gsets=names(lg), stringsAsFactors=F)
for (o in names(lg[[1]]))
{
    d[, o] <- "view"
    d[, paste0(o, "<br>heatmap")] <- "view"
}
rownames(d) <- d[, 1]

links <- plots <- vector('list',length=ncol(d))
names(links) <- names(plots) <- colnames(d)
for (o in names(lg[[1]]))
{
    links[[o]] <- paste0("../G5.ROAST_GSA/rgsa_results/", rownames(d), "_", o, ".html")
    links[[o]][nrow(d)] <- paste0("rgsa_results/", rownames(d)[nrow(d)], "_", o, ".html")    
    links[[paste0(o, "<br>heatmap")]] <- paste0("../G5.ROAST_GSA/heatmaps_gsets/", rownames(d), "_", o, ".png")
    links[[paste0(o, "<br>heatmap")]][nrow(d)] <- paste0("heatmaps_gsets/", rownames(d)[nrow(d)], "_", o, ".png")
    plots[[paste0(o, "<br>heatmap")]] <- paste0("../G5.ROAST_GSA/heatmaps_gsets/", rownames(d), "_", o, ".png")
    plots[[paste0(o, "<br>heatmap")]][nrow(d)] <- paste0("heatmaps_gsets/", rownames(d)[nrow(d)], "_", o, ".png")
}

write.html.mod(d, file=paste0(.res, "/RGSA_maxmean_summary.html"),
               links=links,tiny.pic=plots,
                              title="<br>RGSA results maxmean statistics<br><br>",tiny.pic.size = 100)



################################################################################################################
# List of GO terms for REVIGO
rm(list=ls())
lg <- readRDS(paste0(.res, "lrgsa.maxmean.rds"))

lr <- lapply(names(lg),function(cont){
  l <- lg[[cont]]
  aux <- l[c("GOBP_v5","GOCC_v5","GOMF_v5")]
  goList <- lapply(names(aux), function(j){
    ll <- aux[[j]]
    goterms <- gsub("*(GO\\d+)\\s.*$","\\1",rownames(ll$res))
    goterms <- gsub("GO","GO:",goterms)
    go.dd <- data.frame(goterm = goterms,
                        database = j,
                        pvalue = ll$res[,"pval"],
                        contrast = cont
    )
    go.dd
  })
  res <- do.call(rbind,goList)
  res
})

goTermsRes <- do.call(rbind,lr)
head(goTermsRes)
write.table(goTermsRes,paste0(.res,"GOterms.bp.mf.cc.lineal.csv"),
            sep=",")

goTermsRes.signif <- goTermsRes[goTermsRes$pvalue<0.05,]
summary(goTermsRes.signif)
write.table(goTermsRes.signif,paste0(.res,"GOterms.bp.mf.cc.lineal.signif.csv"),
            sep=",")


################################################################################################################
################################################################################################################
################################################################################################################


