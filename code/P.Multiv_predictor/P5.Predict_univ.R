#################################################################################################################
#################################################################################################################
### P5.Predict_univ.R
### Project: jcalvet_202110_marato
#################################################################################################################
#################################################################################################################
  # Libraries, paths and options

rm(list=ls())
.res <- paste0(.res0, "P.Multiv_predictor/P5.Predict_univ/")
dir.create(.res, F, T)

#################################################################################################################
  # Prepare data

rm(list=ls())
ndd <- readRDS(paste0(.res0, "P.Multiv_predictor/P1.Data_setup/nd.dia.corr.rds"))
ndp <- readRDS(paste0(.res0, "P.Multiv_predictor/P1.Data_setup/nd.prm.all.corr.calib.rds"))

samps <- intersect(sampleNames(ndd),sampleNames(ndp))
ndp <- ndp[, !sampleNames(ndp)%in%samps]
ndp <- ndp[, ndp$group!='Control']


nprot <- c("P01011", "P02763", "P0DJI8", "P18428", "Q08830", "P0DJI9", "Q96PD5", "P02654", "P04196",
           "P13796", "P02765","P02766", "P00738", "P29622", "P00740",
           "P22105", "P26038", "Q6UXB8", "O75636", "P29401", "P15169", "P07357",
           "P26038")

ndd <- ndd[nprot, ]
ndp <- ndp[nprot, ]

saveRDS(ndd, paste0(.res, "nd.dia.rds"))
saveRDS(ndp, paste0(.res, "nd.prm.rds"))

################################################################################################################
  # AUC for DIA

rm(list=ls())
nd <- readRDS(paste0(.res, "nd.dia.rds"))
source("code/functions/myfuns.R")


dir.create(paste0(.res, "auc.univ_dia"), F)
dir.create(paste0(.res, "auc.univ_dia/mild_asymp/"), F)
dir.create(paste0(.res, "auc.univ_dia/sev_mild/"), F)
dir.create(paste0(.res, "auc.univ_dia/cr_sev/"), F)
RNGkind("L'Ecuyer-CMRG")
set.seed(284591)
lr <- sapply(1:nrow(nd), function(j)
{
    nds <- nd[, nd$group%in%c("Asymptomatic", "Mild")]
    nds$group <- droplevels(nds$group)
    x <- exprs(nds)[j, ]
    gr <- nds$group
    s <- !is.na(x) & !is.na(gr)
    x <- x[s]
    gr <- gr[s]
    th <- sort(x)[sum(gr==levels(gr)[1], na.rm=T)]
    probs <- x
    prs <- factor(as.numeric(probs > th), levels=0:1)
    r1 <- facc(gr, prs, probs, fn=paste0(.res, "auc.univ_dia/mild_asymp/",
                                         tolower(fData(nd)$gene.symbol[j]), ".png"), th=th)
    nds <- nd[, nd$group%in%c("Mild", "Severe")]
    nds$group <- droplevels(nds$group)
    x <- exprs(nds)[j, ]
    gr <- nds$group
    s <- !is.na(x) & !is.na(gr)
    x <- x[s]
    gr <- gr[s]
    th <- sort(x)[sum(gr==levels(gr)[1], na.rm=T)]
    probs <- x
    prs <- factor(as.numeric(probs > th), levels=0:1)
    r2 <- facc(gr, prs, probs, fn=paste0(.res, "auc.univ_dia/sev_mild/",
                                        tolower(fData(nd)$gene.symbol[j]), ".png"), th=th)
    nds <- nd[, nd$group%in%c("Severe", "Critical")]
    nds$group <- droplevels(nds$group)
    x <- exprs(nds)[j, ]
    gr <- nds$group
    s <- !is.na(x) & !is.na(gr)
    x <- x[s]
    gr <- gr[s]
    th <- sort(x)[sum(gr==levels(gr)[1], na.rm=T)]
    probs <- x
    prs <- factor(as.numeric(probs > th), levels=0:1)
    r3 <- facc(gr, prs, probs, fn=paste0(.res, "auc.univ_dia/cr_sev/",
                                        tolower(fData(nd)$gene.symbol[j]), ".png"), th=th)
    list(mild.asymp=r1, sev.mild=r2, cr.sev=r3)
}, simplify=F)

names(lr) <- fData(nd)$gene.symbol

saveRDS(lr, paste0(.res, "lauc.univ.dia.rds"))

################################################################################################################
  # AUC for PRM

rm(list=ls())
nd <- readRDS(paste0(.res, "nd.prm.rds"))
source("code/functions/myfuns.R")


dir.create(paste0(.res, "auc.univ_prm"), F)
dir.create(paste0(.res, "auc.univ_prm/mild_asymp/"), F)
dir.create(paste0(.res, "auc.univ_prm/sev_mild/"), F)
dir.create(paste0(.res, "auc.univ_prm/cr_sev/"), F)
RNGkind("L'Ecuyer-CMRG")
set.seed(284591)
lr <- sapply(1:nrow(nd), function(j)
{
    nds <- nd[, nd$group%in%c("Asymptomatic", "Mild")]
    nds$group <- droplevels(nds$group)
    x <- exprs(nds)[j, ]
    gr <- nds$group
    s <- !is.na(x) & !is.na(gr)
    x <- x[s]
    gr <- gr[s]
    th <- sort(x)[sum(gr==levels(gr)[1], na.rm=T)]
    probs <- x
    prs <- factor(as.numeric(probs > th), levels=0:1)
    r1 <- facc(gr, prs, probs, fn=paste0(.res, "auc.univ_prm/mild_asymp/",
                                         tolower(fData(nd)$gene.symbol[j]), ".png"), th=th)
    nds <- nd[, nd$group%in%c("Mild", "Severe")]
    nds$group <- droplevels(nds$group)
    x <- exprs(nds)[j, ]
    gr <- nds$group
    s <- !is.na(x) & !is.na(gr)
    x <- x[s]
    gr <- gr[s]
    th <- sort(x)[sum(gr==levels(gr)[1], na.rm=T)]
    probs <- x
    prs <- factor(as.numeric(probs > th), levels=0:1)
    r2 <- facc(gr, prs, probs, fn=paste0(.res, "auc.univ_prm/sev_mild/",
                                         tolower(fData(nd)$gene.symbol[j]), ".png"), th=th)
    nds <- nd[, nd$group%in%c("Severe", "Critical")]
    nds$group <- droplevels(nds$group)
    x <- exprs(nds)[j, ]
    gr <- nds$group
    s <- !is.na(x) & !is.na(gr)
    x <- x[s]
    gr <- gr[s]
    th <- sort(x)[sum(gr==levels(gr)[1], na.rm=T)]
    probs <- x
    prs <- factor(as.numeric(probs > th), levels=0:1)
    r3 <- facc(gr, prs, probs, fn=paste0(.res, "auc.univ_prm/cr_sev/",
                                         tolower(fData(nd)$gene.symbol[j]), ".png"), th=th)
    list(mild.asymp=r1, sev.mild=r2, cr.sev=r3)
}, simplify=F)

names(lr) <- fData(nd)$gene.symbol


saveRDS(lr, paste0(.res, "lauc.univ.prm.rds"))


### Exemple

source("code/functions/myfuns.R")
nds <- nd[, nd$group%in%c("Mild", "Severe")]
nds$group <- droplevels(nds$group)
j <- which(fData(nd)$gene.symbol=='APOC1')
x <- exprs(nds)[j, ]
gr <- nds$group
s <- !is.na(x) & !is.na(gr)
x <- x[s]
gr <- gr[s]
th <- sort(x)[sum(gr==levels(gr)[1], na.rm=T)]
probs <- x
prs <- factor(as.numeric(probs > th), levels=0:1)
r2 <- facc2(gr, prs, probs, fn=paste0(.res, "auc.univ_prm/sev_mild/",
                                     tolower(fData(nd)$gene.symbol[j]), ".pdf"), th=th)



################################################################################################################
  # Compare results between DIA and PRM

rm(list=ls())
lr1 <- readRDS(paste0(.res, "lauc.univ.dia.rds"))
lr2 <- readRDS(paste0(.res, "lauc.univ.prm.rds"))


r1 <- t(sapply(lr1, function(l) sapply(l, function(o) o$rc["auc"])))
r2 <- t(sapply(lr2, function(l) sapply(l, function(o) o$rc["auc"])))

r1 <- cbind(rownames(r1), r1)
r2 <- cbind(rownames(r2), r2)

write.table(r1, paste0(.res, "auc_dia.csv"), sep=';', row.names=F)
write.table(r2, paste0(.res, "auc_prm.csv"), sep=';', row.names=F)


r <- cbind(unique(r1), unique(r2[match(rownames(r1), rownames(r2)), -1]))
r <- rbind(colnames(r), r)
r[1, ] <- gsub("\\.", " - ", gsub("\\.auc", "", r[1, ]))
r[1, ] <- gsub("cr", "Critical", r[1, ])
r[1, ] <- gsub("sev", "Severe", r[1, ])
r[1, ] <- gsub("mild", "Mild", r[1, ])
r[1, ] <- gsub("asymp", "Asymptomatic", r[1, ])
r <- rbind(c("", "DIA", rep("", 2), "PRM", rep("", 2)), r)
r <- gsub(" \\[", "<br>\\[", r)

p <- openPage(paste0(.res, 'AUC_univ.html'))
hwrite(paste('<br>AUC univariate in DIA and PRM for validated <br><br>'), p)
hwrite(r, page=p, center=F, row.names=F, col.names=F,
       col.style=c("text-align:center", rep("text-align:center", ncol(r)-1)),
       col.width=c(100, rep(150, ncol(r)-1)))
closePage(p)


################################################################################################################
################################################################################################################
################################################################################################################


