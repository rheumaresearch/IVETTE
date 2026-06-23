#################################################################################################################
#################################################################################################################
### O.Concordance_DIA_PRM_complete/O4.Concordance_DIA_PRM.R
### Project: jcalvet_202110_marato
#################################################################################################################
#################################################################################################################
  # Libraries, paths and options

rm(list=ls())
.res <- paste0(.res0, "O.Concordance_DIA_PRM_complete/O4.Concordance_DIA_PRM/")
dir.create(.res, F, T)

#################################################################################################################
  # Expression sets for sample concordance - corrected by batch and overlapping samples and proteins

rm(list=ls())
nd10 <- readRDS(paste0(.dat, "1_DIA/nd.fp.rds"))
nd20 <- readRDS(paste0(.dat, "4_Concordance_DIA.PRM/nd.prm.all.rds"))

nd <- nd10
nv <- c("1")
nvr <- c("quant.batch")
form <- paste("~ ", paste(c(nv, nvr), collapse=' + '))
design <- model.matrix(as.formula(form), nd)
design[, -1] <- apply(design[, -1], 2, function(o) o - mean(o))
formr <- paste("~ ", paste(nvr, collapse=' + '))
designr <- model.matrix(as.formula(formr), nd)[, -1, drop=F]
designr <- apply(designr, 2, function(o) o - mean(o))
lmf <- lmFit(exprs(nd)[, rownames(design), drop=F], design)
betas <- coef(lmf)
exprs(nd) <- exprs(nd) - betas[, colnames(designr), drop=F]%*%t(designr)
nd1 <- nd

nd <- nd20
nv <- c("1")
nvr <- c("quant.batch")
form <- paste("~ ", paste(c(nv, nvr), collapse=' + '))
design <- model.matrix(as.formula(form), nd)
design[, -1] <- apply(design[, -1], 2, function(o) o - mean(o))
formr <- paste("~ ", paste(nvr, collapse=' + '))
designr <- model.matrix(as.formula(formr), nd)[, -1, drop=F]
designr <- apply(designr, 2, function(o) o - mean(o))
lmf <- lmFit(exprs(nd)[, rownames(design), drop=F], design)
betas <- coef(lmf)
exprs(nd) <- exprs(nd) - betas[, colnames(designr), drop=F]%*%t(designr)
nd2 <- nd

samps <- intersect(sampleNames(nd1), sampleNames(nd2))
prots <- intersect(featureNames(nd1), featureNames(nd2))

nd1 <- nd1[prots, samps]
nd2 <- nd2[prots, samps]


saveRDS(nd1, paste0(.res, "nd.dia.corr.rds"))
saveRDS(nd2, paste0(.res, "nd.prm.corr.rds"))

#################################################################################################################
  # Patients descriptives

rm(list=ls())
nd <- readRDS(paste0(.res0, "O.Concordance_DIA_PRM_complete/O3.Diff_exprs_hosp_DIA_sel/nd.log.filt.rds"))
source("code/functions/clin_descrip.R")

nd$group2 <- factor(as.numeric(nd$group%in%c("Severe", "Critical")),
                    levels=0:1, labels=c("Not.hospitalized", "Hospitalized"))

d <- pData(nd)

ny <- "group2"
nxs <- c("sex", "age", "aht", "dm", "obesity", 
         "cortis", "remdesivir", "toci")

nv <- nxs
names(nv) <- nv


### P-values

pvs <- sapply(nv, function(o)
{
    s <- rep(TRUE, nrow(d))
    y <- d[s, o]
    x <- droplevels(d[s, ny])
    if (is.numeric(y))
    {
        r <- kruskal.test(y ~ x)$p.value
    }else if (is.factor(y))
    {
        r <- fisher.test(table(y, x))$p.value        
    }else if (o=="treatment")
    {
        r <- NA
    }
    r
})


### All samples

s <- d$group%in%c("Asymptomatic", "Mild")
d$bilateral.pneumony <- factor(as.character(d$bilateral.pneumony), levels=c("No", "Yes"))
d$bilateral.pneumony[s] <- "No"
d$nimv.oaf[s] <- "No"
d$icu[s] <- "No"
d$exitus[s] <- "No"
d$cortis[s] <- "No"
d$remdesivir[s] <- "No"
d$toci[s] <- "No"


dir.create(paste(.res, "all/", sep=''), F)
r <- clin.descrip(nv, d, ndir=paste0(.res, "all/"), bc=bcs)
ra <- r


### Levels

lr <- list()
for (j in 1:nlevels(d$group2))
{
    s <- !is.na(d$group) & d[, ny]==levels(d[, ny])[j]
    dir.create(paste(.res, levels(d[, ny])[j], "/", sep=''), F)
    r <- clin.descrip(nv, d[s, ], ndir=paste0(.res, levels(d[, ny])[j], "/"), bc=bcs)
    lr[[levels(d[, ny])[j]]] <- r
}


### Merge

ng <- summary(d[, ny])
ng <- ng[names(ng)!="NA's"]
ngg <- paste0(ng, "<br>(", formatC(ng/sum(ng)*100, format='f', dig=1), "%)")
r <- cbind(ra[, -4], All=ra[, 4], do.call(cbind, sapply(lr, function(o) o[, 4], simplify=F)))
r <- cbind(r, "P-value"=rep("", nrow(r)))
r[match(names(pvs), r[, 1]), ncol(r)] <- base::format.pval(pvs)
r[is.na(r[, ncol(r)]) | gsub(" ", "", r[, ncol(r)])=='NA', ncol(r)] <- ""
colnames(r) <- c("Variab.", "N", "Groups", paste0("All<br>n=", nrow(d)),
                 paste0(levels(d[, ny]), "<br>n=", ngg), "P-value")
r <- rbind(colnames(r), rep("", ncol(r)), r)


### HTML summary

p <- openPage(paste0(.res,  "Descriptives_discovery.html"))
hwrite(paste0('<br><br>Descriptives - Univariate association with condition group<br><br>',
              "Measures: Median(Minimum, Maximum) for continouous variables<br>",
              "N(%) for categorical variables<br><br>"), p)
hwrite(r, page=p, center=F, row.names=F, col.names=F,
       col.style=c("text-align:center", rep("text-align:center", ncol(r)-1)),
       col.width=c("250px", "100px", rep('200px', ncol(r)-3)))
closePage(p)

#################################################################################################################
  # Number of proteins differentially expressed

rm(list=ls())
dc1 <- readRDS(paste0(.res0, "O.Concordance_DIA_PRM_complete/O3.Diff_exprs_hosp_DIA_sel/dcomp.rds"))
dc2 <- readRDS(paste0(.res0, "O.Concordance_DIA_PRM_complete/O2.Diff_exprs_hosp_PRM/dcomp.rds"))
nd1 <- readRDS(paste0(.res0, "O.Concordance_DIA_PRM_complete/O3.Diff_exprs_hosp_DIA_sel/nd.log.filt.rds"))
nd2 <- readRDS(paste0(.res0, "O.Concordance_DIA_PRM_complete/O2.Diff_exprs_hosp_PRM/nd.log.filt.rds"))
source("code/functions/icc.R")


### Number of proteins with no missings

prot1.miss <- apply(exprs(nd1), 1, function(o) sum(is.na(o)))
prot2.miss <- apply(exprs(nd2), 1, function(o) sum(is.na(o)))

prot1.nom <- prot1.miss[prot1.miss==0]
prot2.nom <- prot2.miss[prot2.miss==0]

nprot1.nom <- length(prot1.nom)
nprot2.nom <- length(prot2.nom)

pcprot1.nom <- formatC(nprot1.nom/nrow(nd1)*100, format='f', dig=1)
pcprot2.nom <- formatC(nprot2.nom/nrow(nd2)*100, format='f', dig=1)


### Differentially expressed - FC 1.20, FDR < 10%

s <- !is.na(dc1$"Hospitalized-Not.hospitalized.fc") &
     !is.na(dc1$"Hospitalized-Not.hospitalized.adj.pv") &
     abs(dc1$"Hospitalized-Not.hospitalized.fc") >= 1.20 &
     dc1$"Hospitalized-Not.hospitalized.adj.pv" < 0.10
prot1.sig <- dc1[s, "protein.group"]
nprot1.sig <- sum(s)
pcprot1.sig <- formatC(nprot1.sig/nrow(dc1)*100, format='f', dig=1)

s <- !is.na(dc1$"Hospitalized-Not.hospitalized.fc") &
     !is.na(dc1$"Hospitalized-Not.hospitalized.adj.pv") &
     dc1$"Hospitalized-Not.hospitalized.fc" >= 1.20 &
     dc1$"Hospitalized-Not.hospitalized.adj.pv" < 0.10
prot1.sig.p <- dc1[s, "protein.group"]
nprot1.sig.p <- sum(s)
pcprot1.sig.p <- formatC(nprot1.sig.p/nrow(dc1)*100, format='f', dig=1)

s <- !is.na(dc1$"Hospitalized-Not.hospitalized.fc") &
     !is.na(dc1$"Hospitalized-Not.hospitalized.adj.pv") &
     dc1$"Hospitalized-Not.hospitalized.fc" <= -1.20 &
     dc1$"Hospitalized-Not.hospitalized.adj.pv" < 0.10
prot1.sig.n <- dc1[s, "protein.group"]
nprot1.sig.n <- sum(s)
pcprot1.sig.n <- formatC(nprot1.sig.n/nrow(dc1)*100, format='f', dig=1)


s <- !is.na(dc2$"Hospitalized-Not.hospitalized.fc") &
     !is.na(dc2$"Hospitalized-Not.hospitalized.adj.pv") &
     abs(dc2$"Hospitalized-Not.hospitalized.fc") >= 1.20 &
     dc2$"Hospitalized-Not.hospitalized.adj.pv" < 0.10
prot2.sig <- dc2[s, "protein.id"]
nprot2.sig <- sum(s)
pcprot2.sig <- formatC(nprot2.sig/nrow(dc2)*100, format='f', dig=1)

s <- !is.na(dc2$"Hospitalized-Not.hospitalized.fc") &
     !is.na(dc2$"Hospitalized-Not.hospitalized.adj.pv") &
     dc2$"Hospitalized-Not.hospitalized.fc" >= 1.20 &
     dc2$"Hospitalized-Not.hospitalized.adj.pv" < 0.10
prot2.sig.p <- dc2[s, "protein.id"]
nprot2.sig.p <- sum(s)
pcprot2.sig.p <- formatC(nprot2.sig.p/nrow(dc2)*100, format='f', dig=1)

s <- !is.na(dc2$"Hospitalized-Not.hospitalized.fc") &
     !is.na(dc2$"Hospitalized-Not.hospitalized.adj.pv") &
     dc2$"Hospitalized-Not.hospitalized.fc" <= -1.20 &
     dc2$"Hospitalized-Not.hospitalized.adj.pv" < 0.10
prot2.sig.n <- dc2[s, "protein.id"]
nprot2.sig.n <- sum(s)
pcprot2.sig.n <- formatC(nprot2.sig.n/nrow(dc2)*100, format='f', dig=1)

lr <- list(dia.p=prot1.sig.p, dia.n=prot1.sig.n, dia.a=prot1.sig,
           prm.p=prot2.sig.p, prm.n=prot2.sig.n, prm.2=prot2.sig)
saveRDS(lr, paste0(.res, "lprot.sign.dia.prm.rds"))


### Differentially expressed - FC 1.20, pv < 0.05

s <- !is.na(dc1$"Hospitalized-Not.hospitalized.fc") &
     !is.na(dc1$"Hospitalized-Not.hospitalized.adj.pv") &
     abs(dc1$"Hospitalized-Not.hospitalized.fc") >= 1.20 &
     dc1$"Hospitalized-Not.hospitalized.pv" < 0.05
prot1.sig2 <- dc1[s, "protein.group"]
nprot1.sig2 <- sum(s)
pcprot1.sig2 <- formatC(nprot1.sig2/nrow(dc1)*100, format='f', dig=1)

s <- !is.na(dc1$"Hospitalized-Not.hospitalized.fc") &
     !is.na(dc1$"Hospitalized-Not.hospitalized.adj.pv") &
     dc1$"Hospitalized-Not.hospitalized.fc" >= 1.20 &
     dc1$"Hospitalized-Not.hospitalized.pv" < 0.05
prot1.sig.p2 <- dc1[s, "protein.group"]
nprot1.sig.p2 <- sum(s)
pcprot1.sig.p2 <- formatC(nprot1.sig.p2/nrow(dc1)*100, format='f', dig=1)

s <- !is.na(dc1$"Hospitalized-Not.hospitalized.fc") &
     !is.na(dc1$"Hospitalized-Not.hospitalized.adj.pv") &
     dc1$"Hospitalized-Not.hospitalized.fc" <= -1.20 &
     dc1$"Hospitalized-Not.hospitalized.pv" < 0.05
prot1.sig.n2 <- dc1[s, "protein.group"]
nprot1.sig.n2 <- sum(s)
pcprot1.sig.n2 <- formatC(nprot1.sig.n2/nrow(dc1)*100, format='f', dig=1)


s <- !is.na(dc2$"Hospitalized-Not.hospitalized.fc") &
     !is.na(dc2$"Hospitalized-Not.hospitalized.adj.pv") &
     abs(dc2$"Hospitalized-Not.hospitalized.fc") >= 1.20 &
     dc2$"Hospitalized-Not.hospitalized.pv" < 0.05
prot2.sig2 <- dc2[s, "protein.id"]
nprot2.sig2 <- sum(s)
pcprot2.sig2 <- formatC(nprot2.sig2/nrow(dc2)*100, format='f', dig=1)

s <- !is.na(dc2$"Hospitalized-Not.hospitalized.fc") &
     !is.na(dc2$"Hospitalized-Not.hospitalized.adj.pv") &
     dc2$"Hospitalized-Not.hospitalized.fc" >= 1.20 &
     dc2$"Hospitalized-Not.hospitalized.pv" < 0.05
prot2.sig.p2 <- dc2[s, "protein.id"]
nprot2.sig.p2 <- sum(s)
pcprot2.sig.p2 <- formatC(nprot2.sig.p2/nrow(dc2)*100, format='f', dig=1)

s <- !is.na(dc2$"Hospitalized-Not.hospitalized.fc") &
     !is.na(dc2$"Hospitalized-Not.hospitalized.adj.pv") &
     dc2$"Hospitalized-Not.hospitalized.fc" <= -1.20 &
     dc2$"Hospitalized-Not.hospitalized.pv" < 0.05
prot2.sig.n2 <- dc2[s, "protein.id"]
nprot2.sig.n2 <- sum(s)
pcprot2.sig.n2 <- formatC(nprot2.sig.n2/nrow(dc2)*100, format='f', dig=1)


### Intersections

nproti.nom <- length(intersect(names(prot1.nom), names(prot2.nom)))
pcproti.nom <- formatC(nproti.nom/nrow(dc1)*100, format='f', dig=1)


nproti.sig <- length(intersect(prot1.sig.p, prot2.sig.p)) +
               length(intersect(prot1.sig.n, prot2.sig.n))
pcproti.sig <- formatC(nproti.sig/nrow(dc1)*100, format='f', dig=1)

nproti.sig2 <- length(intersect(prot1.sig.p2, prot2.sig.p2)) +
               length(intersect(prot1.sig.n2, prot2.sig.n2))
pcproti.sig2 <- formatC(nproti.sig2/nrow(dc1)*100, format='f', dig=1)


### Create table


tab <- cbind("DIA"=
                 c(nrow(dc1),
                   paste0(nprot1.nom, "<br>(", pcprot1.nom, "%)"),
                   paste0(nprot1.sig, "<br>(", pcprot1.sig, "%)"),
                   paste0(nprot1.sig2, "<br>(", pcprot1.sig2, "%)")),
             "PRM"=
                 c(nrow(dc2),
                   paste0(nprot2.nom, "<br>(", pcprot2.nom, "%)"),
                   paste0(nprot2.sig, "<br>(", pcprot2.sig, "%)"),
                   paste0(nprot2.sig2, "<br>(", pcprot2.sig2, "%)")),
             "Intersection"=
                 c(nrow(dc1),
                   paste0(nproti.nom, "<br>(", pcproti.nom, "%)"),
                   paste0(nproti.sig, "<br>(", pcproti.sig, "%)"),
                   paste0(nproti.sig2, "<br>(", pcproti.sig2, "%)")))
rownames(tab) <- c("Total", "Proteins with no missings",
                   "Diferentially expressed - FC 1.20; 10%FDR",
                   "Diferentially expressed - FC 1.20; p-value < 0.05")


tabc <- rbind(colnames(tab), tab)
tabc <- cbind(rownames(tabc), tabc)


### Save as HTML

p <- openPage(paste0(.res, 'Summary_proteins_identified.html'))
hwrite(paste('<br>Summary of identified proteins across technologies<br><br>'), p)
hwrite(tabc, page=p, center=F, row.names=F, col.names=F,
       col.style=c("text-align:center", rep("text-align:center", ncol(tabc)-1)),
       col.width=c(320, rep(160, ncol(tabc)-1)))
closePage(p)

#################################################################################################################
  # FC correlations - DIA - DDA-PD

rm(list=ls())
dc1 <- readRDS(paste0(.res0, "O.Concordance_DIA_PRM_complete/O3.Diff_exprs_hosp_DIA_sel/dcomp.rds"))
dc2 <- readRDS(paste0(.res0, "O.Concordance_DIA_PRM_complete/O2.Diff_exprs_hosp_PRM/dcomp.rds"))
source("code/functions/icc.R")
source("code/functions/make_transparent.R")

comps <- colnames(dc1)[regexpr("\\.t$", colnames(dc1)) > 0]
comps <- gsub("\\.t$", "", comps)

.fpp <- function(x, y)
{
    alfa <- 0.05
    m <- lm(y ~ x)
    a <- m$coef[1]
    b <- m$coef[2]
    sm <- summary(m)$coef
    R2 <- summary(m)$r.squared
    ree <- summary(m)$sigma
    s <- !is.na(x) & !is.na(y)
    ccc.v <- ic.ccc(x=x[s], y=y[s], alfa=0.05)
    dif95 <- qt(1-alfa/2, df=sum(s)-2)*ree
    cr <- cor.test(x, y)
    pv.int <- sm[1,4]
    pv.slope <- pt(abs(sm[2, 1]-1)/sm[2,2], df=m$df, lower.tail=F)*2
    ac <- paste(formatC(m$coef[1], format='f', dig=3), " (",
                formatC(confint(m)[1, 1], format='f', dig=3), " to ",
                formatC(confint(m)[1, 2], format='f', dig=3), ")", sep='')
    bc <- paste(formatC(m$coef[2], format='f', dig=3), " (",
                formatC(confint(m)[2, 1], format='f', dig=3), " to ",
                formatC(confint(m)[2, 2], format='f', dig=3), ")", sep='')
    R2c <- formatC(R2, format='f', dig=3)
    reec <- formatC(ree, format='f', dig=3)
    ccc.c <- paste(formatC(ccc.v["ccc"], format='f', dig=3), " (",
                   formatC(ccc.v["llci"], format='f', dig=3), " to ",
                   formatC(ccc.v["upci"], format='f', dig=3), ")", sep='')
    cbc <- formatC(ccc.v["cb"], format='f', dig=3)
    uc <- formatC(ccc.v["u"], format='f', dig=3)
    vc <- formatC(ccc.v["v"], format='f', dig=3)
    crc <- paste(formatC(cr$estimate, format='f', dig=3), " (",
                 formatC(cr$conf.int[1], format='f', dig=3), " to ",
                 formatC(cr$conf.int[2], format='f', dig=3), ")", sep='')
    maxdifc <- paste(formatC(dif95, format='f', dig=2),
                     " (", formatC(dif95*100/diff(range(x, na.rm=T)),
                                   format='f', dig=2), "%)", sep='')
    c(
#                     paste("Int.[95%CI] =", ac),
#                     paste("Slope[95%CI] =", bc),
#                     paste("R2 =", R2c),
#                     paste("Resid.ee =", reec),
        paste("CCC =", ccc.c),
        paste("PC  = ", crc, sep=''),

#        paste("U =", uc),
#        paste("V = ", vc),
        paste("Cb =", cbc))
#        paste("Max.dif.95% =", maxdifc, sep=''))
}
fpp <- .fpp

comp <- comps[1]
d1 <- dc1[, c("protein.group", paste0(comp, ".log2fc"), paste0(comp, ".pv"), paste0(comp, ".adj.pv"))]
d2 <- dc2[, c("protein.id", paste0(comp, ".log2fc"), paste0(comp, ".pv"), paste0(comp, ".adj.pv"))]
colnames(d1) <- colnames(d2) <- c("prot", "log2fc", "pv", "adj.pv")
prots12 <- intersect(rownames(d1), rownames(d2))
pdf(paste0(.res, "log2FC_corr_", comp, "_DIA_PRM.pdf"), width=8, height=8)
prots <- prots12
x <- d1[prots, "log2fc"]
y <- d2[prots, "log2fc"]
ccc <- fpp(x, y)
m <- lm(y ~ x)
cr <- cor(x, y, method='spearman', use='pairwise.complete.obs')
xats <- 2^seq(-100, 100, 1)
plot(x, y, ylab="PRM", xlab="DIA", axes=F)
axis(1, at=log2(xats), labels=xats,
     lwd=1.5)
axis(2, at=log2(xats), labels=xats, lwd=1.5, las=2)
box(lwd=1.5)
legend(x='topleft', inset=0.01, bty='n', pch=NA,
       legend=c(ccc))
legend(x='bottomright', inset=0.01, bty='n', pch=NA,
       legend=c(paste0(length(prots), " proteins")))
abline(a=m$coef[1], b=m$coef[2], lty=3, col='red4', lwd=3)
abline(v=0, h=0, lty=2, col='gray30', lwd=1.5)
abline(a=0, b=1, lty=3, lwd=3, col='gray30')
dev.off()


pdf(paste0(.res, "log2FC_corr_", comp, "_DIA_PRM_zoom.pdf"), width=8, height=8)
par(mar=c(5, 5, 5, 5))
prots <- prots12
x <- d1[prots, "log2fc"]
y <- d2[prots, "log2fc"]
ccc <- fpp(x, y)
m <- lm(y ~ x)
cr <- cor(x, y, method='spearman', use='pairwise.complete.obs')
xats <- 2^seq(-100, 100, 1)
plot(x, y,
     ylab="PRM (FC)",
     xlab="DIA (FC)", ylim=c(-1.20, 3), xlim=c(-1.20, 3), axes=F,
     cex.axis=1.75, cex.lab=1.75,
     col=makeTransparent(rep("blue3", length(x)), 0.01), cex=2.5, pch=16)
axis(1, at=log2(xats), labels=xats, lwd=1.5, cex.axis=1.7)
axis(2, at=log2(xats), labels=xats, lwd=1.5, las=2, cex.axis=1.7)
box(lwd=1.5)
abline(a=m$coef[1], b=m$coef[2], lty=2, col='red2', lwd=2.3)
abline(a=0, b=1, lty=2, lwd=2.5, col='gray30')
abline(v=0, h=0, lty=3, col='gray30', lwd=2.5)
legend(x='topleft', inset=0.01, bty='n', pch=NA, legend=c(ccc), cex=1.5)
legend(x='bottomright', inset=0.01, pch=NA,
       bty='n', 
       lty=c(2, 2), col=rev(c("gray30", "red2")), seg.len=3.5,
       legend=rev(c("Perfect concordance", "Regression line")),
       cex=1.5,
       title=paste0(length(prots), " proteins"),
       lwd=4)
dev.off()

#################################################################################################################
  # Sample concordance - DIA vs PRM

rm(list=ls())
nd1 <- readRDS(paste0(.res, "nd.dia.corr.rds"))
nd2 <- readRDS(paste0(.res, "nd.prm.corr.rds"))
source("code/functions/icc.R")
source("code/functions/make_transparent.R")


nd1$group2 <- factor(as.numeric(nd1$group%in%c("Severe", "Critical")),
                    levels=0:1, labels=c("Not.hospitalized", "Hospitalized"))
nd2$group2 <- factor(as.numeric(nd2$group%in%c("Severe", "Critical")),
                    levels=0:1, labels=c("Not.hospitalized", "Hospitalized"))


fpp <- function(x, y)
{
    alfa <- 0.05
    m <- lm(y ~ x)
    a <- m$coef[1]
    b <- m$coef[2]
    sm <- summary(m)$coef
    R2 <- summary(m)$r.squared
    ree <- summary(m)$sigma
    s <- !is.na(x) & !is.na(y)
    ccc.v <- ic.ccc(x=x[s], y=y[s], alfa=0.05)
    dif95 <- qt(1-alfa/2, df=sum(s)-2)*ree
    cr <- cor.test(x, y)
    crr <- cor(cbind(x, y, as.numeric(nd1$group2)-1))
    pr <- partial.r(crr, c(1, 2), 3)
    cp <- corr.p(pr[1, 2, drop=F], n=length(x)-1, adjust='none')
    cp <- as.numeric(c(cp$ci[c(2, 1, 3)], cp$p[1, 1]))
    names(cp) <- c("cr", "l", "u", "p")
    pv.int <- sm[1,4]
    pv.slope <- pt(abs(sm[2, 1]-1)/sm[2,2], df=m$df, lower.tail=F)*2
    ac <- paste(formatC(m$coef[1], format='f', dig=3), " [",
                formatC(confint(m)[1, 1], format='f', dig=3), ", ",
                formatC(confint(m)[1, 2], format='f', dig=3), "]", sep='')
    bc <- paste(formatC(m$coef[2], format='f', dig=3), " [",
                formatC(confint(m)[2, 1], format='f', dig=3), ", ",
                formatC(confint(m)[2, 2], format='f', dig=3), "]", sep='')
    R2c <- formatC(R2, format='f', dig=3)
    reec <- formatC(ree, format='f', dig=3)
    ccc.c <- paste(formatC(ccc.v["ccc"], format='f', dig=3), " [",
                   formatC(ccc.v["llci"], format='f', dig=3), ", ",
                   formatC(ccc.v["upci"], format='f', dig=3), "]", sep='')
    cbc <- formatC(ccc.v["cb"], format='f', dig=3)
    uc <- formatC(ccc.v["u"], format='f', dig=3)
    vc <- formatC(ccc.v["v"], format='f', dig=3)
    crc <- paste(formatC(cr$estimate, format='f', dig=3), " [",
                 formatC(cr$conf.int[1], format='f', dig=3), ", ",
                 formatC(cr$conf.int[2], format='f', dig=3), "]", sep='')
    maxdifc <- paste(formatC(dif95, format='f', dig=2),
                     " (", formatC(dif95*100/diff(range(x, na.rm=T)),
                                   format='f', dig=2), "%)", sep='')
    r1 <- c(paste("CCC [95%CI] =", ccc.c),
           paste("Cb =", cbc),
           paste("U =", uc),
           paste("V = ", vc),
           paste("PC[95%CI] = ", crc, sep=''),
           paste("Max.dif.95% =", maxdifc, sep=''),
           paste("Int.[95%CI] =", ac),
           paste("Slope[95%CI] =", bc),
           paste("R2 =", R2c),
           paste("Resid.ee =", reec))
    r2 <- c(ccc.v["ccc"], ccc.v["cb"], ccc.v["u"], ccc.v["v"], cr=cr$estimate, cp=cp[1])
    list(r1, r2)
}


samps <- colnames(nd1)
prots <- intersect(featureNames(nd1), featureNames(nd2))

nd1 <- nd1[prots, samps]
nd2 <- nd2[prots, samps]

s1 <- apply(exprs(nd1), 1, function(o) sum(!is.na(o)) >= 10)
s2 <- apply(exprs(nd2), 1, function(o) sum(!is.na(o)) >= 10)
s <- s1 & s2

nd1 <- nd1[s, ]
nd2 <- nd2[s, ]

sampleNames(nd1) <- paste0(sampleNames(nd1), "-DIA")
sampleNames(nd2) <- paste0(sampleNames(nd2), "-PRM")


dir.create(paste0(.res, "ccc_dia_dda.pd"), F)
lr <- mclapply(1:nrow(nd1), function(i)
{
    x <- exprs(nd1)[i, ]
    y <- exprs(nd2)[i, ]
    s <- !is.na(x) & !is.na(y)
    x <- x[s]
    y <- y[s]
    if (sum(s) >= 10)
    {
        ylim <- xlim <- range(c(x, y))
        png(paste0(.res, "ccc_dia_dda.pd/", rownames(nd1)[i], ".png"), width=800, height=800)
        r <- fpp(x, y)
        ccc <- r[[1]]
        m <- lm(y ~ x)
        cr <- cor(x, y, method='spearman', use='pairwise.complete.obs')
        plot(x, y, ylab="DDA-PD", xlab="DIA", xlim=xlim, ylim=ylim)
        legend(x='bottomright', inset=0.01, bty='n', pch=NA,
               legend=c(ccc))
        abline(a=m$coef[1], b=m$coef[2], lty=3, col='red4', lwd=3)
        abline(v=0, h=0, lty=3, col='gray30')
        abline(a=0, b=1, lty=3, col='gray30')    
        dev.off()
        r[[2]]
    }else
    {
        NA
    }
}, mc.cores=12)
names(lr) <- rownames(nd1)

r <- do.call(rbind, lr)
r <- r[!is.na(r[, 1]), ]

saveRDS(r, paste0(.res, "res.conc.dia.prm.rds"))


### Version for report - density

pdf(paste0(.res, "ccc_dia_prm.pdf"), width=8, height=8)
par(mfrow=c(2, 2))
x <- r[, 1]
xt <- gtools::logit(x, -1, 1)
dens <- density(xt, n=10000, adjust = 3)
dens$xo <- gtools::inv.logit(dens$x, -1, 1)
dens$yo <- dens$y  / (1/(dens$xo*(1 - dens$xo)))
dens$yo[dens$yo < 0] <- 0
plot(dens$xo, dens$yo, main='CCC', type='l', lwd=4, col='red4',
     ylab='Density', xlab='')
x <- r[, 2]
xt <- gtools::logit(x, 0, 1)
dens <- density(xt, n=10000, adjust = 1)
dens$xo <- gtools::inv.logit(dens$x, 0, 1)
dens$yo <- dens$y  / (1/(dens$xo*(1 - dens$xo)))
dens$yo[dens$yo < 0] <- 0
plot(dens$xo, dens$yo, main='Cb', type='l', lwd=4, col='red4',
     ylab='Density', xlab='')
x <- r[, 5]
xt <- gtools::logit(x, -1, 1)
dens <- density(xt, n=10000, adjust = 1)
dens$xo <- gtools::inv.logit(dens$x, -1, 1)
dens$yo <- dens$y  / (1/(dens$xo*(1 - dens$xo)))
dens$yo[dens$yo < 0] <- 0
plot(dens$xo, dens$yo, main='PC', type='l', lwd=4, col='red4',
     ylab='Density', xlab='')
dev.off()


### Version for report - boxplots

rr <- r[, 1:5]
colnames(rr) <- c("CCC", "Cb", "U", "V", "PC")
rr <- rr[, c(1, 5, 2, 4, 3)]
rr[, "U"] <- 1.5*((r[, "u"] - (-8.1))/(12 - (-8.1)))  -0.5
rr[, "V"] <- 1.5*((r[, "v"] - 0.24)/(2.8 - 0.24)) - 0.5
lr <- as.list(as.data.frame(rr))
pdf(paste0(.res, "ccc_dia_prm.pdf"), width=11, height=8)
par(mar=c(5, 7, 5, 6))
stripchart(lr, vertical=T, method='jitter', jitter=0.20, at=c(0.9, 2, 3.1, 5.1, 7.1),
           pch=1,
           col=makeTransparent(rep("gray70", nrow(rr)), 0.1), xlim=c(0.4, 7.5), cex=1.5, axes=F)
bp <- boxplot(rr, pch='', at=c(0.9, 2, 3.1, 5.1, 7.1),
        col=NA,
        border="red4", lwd=3, add=T, names=NA, axes=F)
axis(1, at=c(-100, 100), labels=c("", ""), lwd=2.5)
axis(2, at=seq(-10, 10, 0.20), lwd=2.5, las=2, cex.axis=1.25)
axis(2, pos=6.25, at=1.5*(seq(-30, 30, 5) - -(8.1))/(12 - (-8.1)) - 0.5,
     labels=seq(-30, 30, 5), lwd=2.5, las=2, cex.axis=1.25)
axis(2, pos=4.25, at=1.5*(c(-3, seq(-2, 3, 0.5) - 0.24)/(2.8 - 0.24)) - 0.5,
     labels=c(-3, seq(-2, 3, 0.5)), lwd=2.5, las=2, cex.axis=1.25)
text(x=-0.7, y=0.25, label="Concordance values", srt=90, xpd=T, cex=1.75)
box(lwd=2.5)
text(x=c(0.9, 2, 3.1, 5.1, 7.1), y=min(rr) - diff(range(rr))*0.09, labels=colnames(rr),
     xpd=T, cex=1.8)
segments(x0=4.25, x1=6.25, y0=1.5*((1 - 0.24)/(2.8 - 0.24)) - 0.5, y1=1.5*((1 - 0.24)/(2.8 - 0.24)) - 0.5,
         lwd=2, lty=2, col='gray30')
segments(x0=6.25, x1=7.75, y0=1.5*((0 - (-8.1))/(12 - (-8.1))) - 0.5,
         y1=1.5*((0 - (-8.1))/(12 - (-8.1))) - 0.5,
         lwd=2, lty=2, col='gray30')
dev.off()


### Descriptives

rr <- t(apply(r, 2, function(o)
{
    med <- median(o)
    mad <- mad(o)
    qts <- quantile(o, c(0, 0.10, 0.25, 0.75, 0.90, 1))
    c(med, mad, qts)
}))
rownames(rr) <- c("CCC", "Cb", "U", "V", "PC", "PartCorr")
rr <- rr[c(1, 5, 2, 4, 3), ]
colnames(rr) <- c("Median", "MAD", "Minimum", "10% Perc.", "25% Perc.", "75% Perc.", "90% Perc.", "Maximum")
rrc <- formatC(rr, format='f', dig=3)

rrc[, 1] <- paste0(rrc[, 1], "<br>(", rrc[, 2], ")")
colnames(rrc)[1] <- "Median<br>(MAD)"
rrc <- rrc[, -2]


rrc <- rbind(colnames(rrc), rrc)
rrc <- cbind(rownames(rrc), rrc)

p <- openPage(paste0(.res, 'Summary_concordance.html'))
hwrite(paste('<br>Summary of concordance indexex across proteins<br><br>'), p)
hwrite(rrc, page=p, center=F, row.names=F, col.names=F,
       col.style=c("text-align:center", rep("text-align:center", ncol(rrc)-1)),
       col.width=c(75, rep(110, ncol(rrc)-1)))
closePage(p)


### Summary for thresholds

ths <- seq(0.1, 0.9, 0.05)
ths2 <- c(seq(0, 5, 0.5), seq(7, 12, 1))

rr <- sapply(1:ncol(r), function(j)
{
    o <- r[, j]
    if (colnames(r)[j] !='u')
    {
        sapply(ths, function(th, o)
        {
            n <- sum(o > th)
            pc <- n/length(o)*100
            paste0(n, "<br>(", formatC(pc, format='f', dig=1), "%)")
        }, o)
    }else
    {
        sapply(ths2, function(th, o)
        {
            n <- sum(abs(o) < th)
            pc <- n/length(o)*100
            paste0(n, "<br>(", formatC(pc, format='f', dig=1), "%)")
        }, o)
    }
})
colnames(rr) <- colnames(r)
colnames(rr) <- c("CCC", "Cb", "U", "V", "PC", "PartCorr")
rr <- rr[, c(1, 5, 2, 4, 3)]
rownames(rr) <- ths
rr <- cbind("Threshold"=rownames(rr), rr, "Threshold"=ths2)
rr <- rbind(colnames(rr), rr)

p <- openPage(paste0(.res, 'Corr_V_by_thresholds.html'))
hwrite(paste('<br>Partial correlation and V index by thresholds<br><br>'), p)
hwrite(rr, page=p, center=F, row.names=F, col.names=F,
       col.style=c("text-align:center", rep("text-align:center", ncol(rr)-1)),
       col.width=c(75, rep(110, ncol(rr)-2), 75))
closePage(p)

#################################################################################################################
  # Sample concordance - DIA vs DDA-PD - protein exemples

rm(list=ls())
nd1 <- readRDS(paste0(.res, "nd.dia.corr.rds"))
nd2 <- readRDS(paste0(.res, "nd.prm.corr.rds"))
conc <- readRDS(paste0(.res, "res.conc.dia.prm.rds"))
source("code/functions/icc.R")
source("code/functions/make_transparent.R")

nd1$group2 <- factor(as.numeric(nd1$group%in%c("Severe", "Critical")),
                    levels=0:1, labels=c("Not.hospitalized", "Hospitalized"))
nd2$group2 <- factor(as.numeric(nd2$group%in%c("Severe", "Critical")),
                    levels=0:1, labels=c("Not.hospitalized", "Hospitalized"))

samps <- colnames(nd1)
prots <- intersect(featureNames(nd1), featureNames(nd2))

nd1 <- nd1[prots, samps]
nd2 <- nd2[prots, samps]

s1 <- apply(exprs(nd1), 1, function(o) sum(!is.na(o)) >= 10)
s2 <- apply(exprs(nd2), 1, function(o) sum(!is.na(o)) >= 10)
s <- s1 & s2

med <- apply(conc[, c("cb", "cr.cor")], 2, median)
dst <- as.matrix(dist(rbind(med, conc), method='euclidean'))

#dsts <- dst[1, -1]
#dsts <- dsts[order(dsts, decreasing=F)]
#prot <- names(dsts)[6]

conc[conc[, "cr.cor"] > 0.85 & conc[, "cr.cor"] < 0.90,]
prot <- "P15169"


### Version 1 - before and after calibration

pdf(paste0(.res, prot, "_DIA_PRM.pdf"), width=9, height=9)
par(mar=c(6, 6, 5, 5))
x <- exprs(nd1)[prot, ]
y <- exprs(nd2)[prot, ]
ccc <- .fpp(x, y)
m <- lm(y ~ x)
cr <- cor(x, y, method='spearman', use='pairwise.complete.obs')
plot(x, y,
     ylab="",
     xlab="",
#     ylim=c(17, 25), xlim=c(18, 21),
#     ylim=c(14, 19), xlim=c(14, 19),     
     axes=F,
     cex.axis=1.75, cex.lab=1.75,
     col=makeTransparent(rep("blue3", length(x)), 0.01), cex=2.5, pch=16)
axis(1, lwd=2, cex.axis=1.4)
axis(2, lwd=2, las=2, cex.axis=1.4)
box(lwd=2)
legend(x='topleft', inset=0.005, bty='n', pch=NA, legend=c(ccc), cex=1.3)
legend(x='bottomright', inset=0.01, pch=NA,
       bty='n', 
       lty=c(2, 2), col=c("gray30", "red4"), seg.len=3,
       legend=c("Perfect concordance", "Regression line"),
       cex=1.5,
       title=paste0(length(prots), " proteins"),
       lwd=3)
abline(a=m$coef[1], b=m$coef[2], lty=2, col='red4', lwd=1.75)
abline(a=0, b=1, lty=2, lwd=1.75, col='gray30')
abline(v=0, h=0, lty=3, col='gray30', lwd=2.5)
text(x=14.5, y=20.5, label="DDA - MaxLQF (log2)", srt=90, xpd=T, cex=1.75)
text(x=20.5, y=14.5, label="DIA - MaxLQF (log2)", xpd=T, cex=1.75)
dev.off()


pdf(paste0(.res, prot, "_DIA_PRM_calib.pdf"), width=9, height=9)
par(mar=c(6, 6, 5, 5))
x <- exprs(nd1)[prot, ]
y <- exprs(nd2)[prot, ]
m <- lm(y ~ x)
y <- (y - m$coef[1])/m$coef[2]
ccc <- .fpp(x, y)
cr <- cor(x, y, method='spearman', use='pairwise.complete.obs')
plot(x, y,
     ylab="",
     xlab="",
#     ylim=c(17, 21), xlim=c(17, 21),
#     ylim=c(14, 19), xlim=c(14, 19),     
     axes=F,
     cex.axis=1.75, cex.lab=1.75,
     col=makeTransparent(rep("blue3", length(x)), 0.01), cex=2.5, pch=16)
axis(1, lwd=2, cex.axis=1.4)
axis(2, lwd=2, las=2, cex.axis=1.4)
box(lwd=2)
legend(x='topleft', inset=0.005, bty='n', pch=NA, legend=c(ccc), cex=1.3)
legend(x='bottomright', inset=0.01, pch=NA,
       bty='n', 
       lty=c(2, 2), col=c("gray30", "red4"), seg.len=3,
       legend=c("Perfect concordance", "Regression line"),
       cex=1.5,
       title=paste0(ncol(nd1), " samples"),
       lwd=3)
abline(a=m$coef[1], b=m$coef[2], lty=2, col='red4', lwd=1.75)
abline(a=0, b=1, lty=2, lwd=1.75, col='gray30')
abline(v=0, h=0, lty=3, col='gray30', lwd=2.5)
text(x=14.5, y=20.5, label="DDA - MaxLQF (log2)", srt=90, xpd=T, cex=1.75)
text(x=20.5, y=14.5, label="DIA - MaxLQF (log2)", xpd=T, cex=1.75)
dev.off()


### Version 2 - before and after calibration

pdf(paste0(.res, prot, "_DIA_PRM_2.pdf"), width=9.5, height=9)
par(mar=c(6, 7, 5, 5))
x <- exprs(nd1)[prot, ]
y <- exprs(nd2)[prot, ]
ccc <- .fpp(x, y)
m <- lm(y ~ x)
cr <- cor(x, y, method='spearman', use='pairwise.complete.obs')
plot(x, y,
     ylab="",
     xlab="",
#     ylim=c(17, 21), xlim=c(17, 21),
     ylim=c(18, 22), xlim=c(17.75, 21),
     axes=F,
     cex.axis=1.75, cex.lab=1.75,
     col=makeTransparent(rep("blue3", length(x)), 0.01), cex=2.5, pch=16)
axis(1, lwd=2, cex.axis=1.8, at=seq(18, 21, 1))
axis(2, lwd=2, las=2, cex.axis=1.8)
box(lwd=2.2)
legend(x='topleft', inset=-0.005, bty='n', pch=NA, legend=c(ccc), cex=1.6)
legend(x='bottomright', inset=0.01, pch=NA,
       bty='n', 
       lty=c(2, 2), col=rev(c("gray30", "red3")), seg.len=3.5,
       legend=rev(c("Perfect concordance", "Regression line")),
       cex=1.75,
       title=paste0(ncol(nd1), " samples"),
       lwd=4)
abline(a=m$coef[1], b=m$coef[2], lty=2, col='red3', lwd=4)
abline(a=0, b=1, lty=2, col='gray30', lwd=4)
abline(v=0, h=0, lty=3, col='gray30', lwd=2.5)
text(x=17.15, y=20, label="PRM - Calibrated MaxLQF (log2)", srt=90, xpd=T, cex=2)
text(x=19.5, y=17.3, label="DIA - MaxLQF (log2)", xpd=T, cex=2)
dev.off()


pdf(paste0(.res, prot, "_DIA_PRM_calib_2.pdf"), width=9.5, height=9)
par(mar=c(6, 7, 5, 5))
x <- exprs(nd1)[prot, ]
y <- exprs(nd2)[prot, ]
m <- lm(y ~ x)
y <- (y - m$coef[1])/m$coef[2]
ccc <- .fpp(x, y)
cr <- cor(x, y, method='spearman', use='pairwise.complete.obs')
plot(x, y,
     ylab="",
     xlab="",
#     ylim=c(17, 21), xlim=c(17, 21),
     ylim=c(18, 22), xlim=c(17.75, 21),
     axes=F,
     cex.axis=1.75, cex.lab=1.75,
     col=makeTransparent(rep("blue3", length(x)), 0.01), cex=2.5, pch=16)
axis(1, lwd=2, cex.axis=1.8, at=seq(15, 22, 1))
axis(2, lwd=2, las=2, cex.axis=1.8)
box(lwd=2.2)
legend(x='topleft', inset=-0.005, bty='n', pch=NA, legend=c(ccc), cex=1.6)
legend(x='bottomright', inset=0.01, pch=NA,
       bty='n', 
       lty=c(2, 2), col=rev(c("gray30", "red3")), seg.len=3.5,
       legend=rev(c("Perfect concordance", "Regression line")),
       cex=1.75,
       title=paste0(ncol(nd1), " samples"),
       lwd=4)
abline(a=m$coef[1], b=m$coef[2], lty=2, col='red3', lwd=4)
abline(a=0, b=1, lty=2, col='gray30', lwd=4)
abline(v=0, h=0, lty=3, col='gray30', lwd=2.5)
text(x=17.15, y=20, label="PRM - Calibrated MaxLQF (log2)", srt=90, xpd=T, cex=2)
text(x=19.5, y=17.3, label="DIA - MaxLQF (log2)", xpd=T, cex=2)
dev.off()

################################################################################################################
################################################################################################################
################################################################################################################


