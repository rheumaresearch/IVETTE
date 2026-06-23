#################################################################################################################
#################################################################################################################
### P.Multiv_predictor/P1.Data_setup.R
### Project: jcalvet_202110_marato
#################################################################################################################
#################################################################################################################
  # Libraries, paths and options

rm(list=ls())
.res <- paste0(.res0, "P.Multiv_predictor/P1.Data_setup/")
dir.create(.res, F, T)

#################################################################################################################
  # Parameters for PRM calibration to DIA

rm(list=ls())
nd1 <- readRDS(paste0(.res0, "O.Concordance_DIA_PRM_complete/O4.Concordance_DIA_PRM/nd.dia.corr.rds"))
nd2 <- readRDS(paste0(.res0, "O.Concordance_DIA_PRM_complete/O4.Concordance_DIA_PRM/nd.prm.corr.rds"))

samps <- intersect(sampleNames(nd1), sampleNames(nd2))
ps <- intersect(featureNames(nd1), featureNames(nd2))

nd1 <- nd1[ps, samps]
nd2 <- nd2[ps, samps]

bs <- t(sapply(featureNames(nd1), function(o)
{
    x <- exprs(nd1)[o, ]
    y <- exprs(nd2)[o, ]    
    m <- lm(y ~ x)
    c(a=as.numeric(m$coef[1]), b=as.numeric(m$coef[2]))
}))

saveRDS(bs, paste0(.res, "bs.calib.prm.rds"))

#################################################################################################################
  # Correct and calibrate PRM dataset to DIA

rm(list=ls())
nd <- readRDS(paste0(.res0, "O.Concordance_DIA_PRM_complete/O1.Merge_PRM_data/nd.prm.all.rds"))
bs <- readRDS(paste0(.res, "bs.calib.prm.rds"))

nd <- nd[regexpr("iRT-Kit_WR_fusion", featureNames(nd)) < 0, ]


### Correct by batch

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


### Calibrate

M <- t(sapply(rownames(bs), function(o)
{
    bss <- bs[o, ]
    y <- exprs(nd)[o, ]
    (y - bss[1])/bss[2]
}))


nd <- new('ExpressionSet', exprs=M, featureData=new('AnnotatedDataFrame', fData(nd)[rownames(M), ]),
          phenoData=new('AnnotatedDataFrame', pData(nd)),
          annotation="PRM-Protein-normarea-MaxLFQ-corrected-calibrated")


saveRDS(nd, paste0(.res, "nd.prm.all.corr.calib.rds"))

#################################################################################################################
  # Correct DIA dataset by technical effects and impute missings

rm(list=ls())
nd <- readRDS(paste0(.res0, "G.Discovery_DIA_FP_MaxLFQ/G2.Diff_expression/nd.log.filt.rds"))


### Correct by quantification batch

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


### Impute missings with eof

pc <- eof(t(exprs(nd)), recursive=T, scaled=FALSE, nu=ncol(nd)-1)
Mr <- t(eofRecon(pc, uncenter=TRUE, unscale=TRUE))
rownames(Mr) <- featureNames(nd)
colnames(Mr) <- sampleNames(nd)

ndr <- new('ExpressionSet', exprs=Mr, featureData=new('AnnotatedDataFrame', fData(nd)),
           phenoData=new('AnnotatedDataFrame', pData(nd)), annotation="DIA-FP-MaxLFQ-corr.eof")


saveRDS(nd, paste0(.res, "nd.dia.corr.rds"))
saveRDS(ndr, paste0(.res, "nd.dia.corr.eof.rds"))

#################################################################################################################
  # Train DIA -> validated in PRM

rm(list=ls())
ndd <- readRDS(paste0(.res, "nd.dia.corr.eof.rds"))
ndp <- readRDS(paste0(.res, "nd.prm.all.corr.calib.rds"))


samps <- intersect(sampleNames(ndd),sampleNames(ndp))
ndp <- ndp[, !sampleNames(ndp)%in%samps]


ps <- names(which(apply(exprs(ndp), 1, function(o) sum(is.na(o)) == 0)))
ndd <- ndd[ps, ]
ndp <- ndp[ps, ]


Mtr <- exprs(ndd)
ytr <- factor(as.numeric(ndd$group%in%c("Severe", "Critical")),
              levels=0:1, labels=c("Not-hospitalized", "Hospitalized"))
Mts <- exprs(ndp)
yts <- factor(as.numeric(ndp$group%in%c("Severe", "Critical")),
            levels=0:1, labels=c("Not-hospitalized", "Hospitalized"))


rf <- randomForest(x=t(Mtr), y=ytr,
                   ntree=10000,
                   mtry=floor(sqrt(nrow(Mtr))),
                   replace=FALSE, 
                   importance=TRUE,
                   do.trace=100)

prs <- predict(rf, t(Mts))
probs <- predict(rf, t(Mts), type='prob')

tab <- table(yts, prs )
acc <- sum(diag(tab))/sum(tab)
sens <- tab[2, 2]/sum(tab[2, ])
spec <- tab[1, 1]/sum(tab[1, ])
phi <- phi(tab, 3)
c(acc=acc, sens=sens, spec=spec, phi=phi)


roc.r <- roc(yts, probs[, 2])


auc <- as.numeric(roc.r$auc)
ci <- as.numeric(ci.auc(roc.r, method='bootstrap'))[c(1, 3)]
aucc <- paste0(formatC(auc, format='f', dig=3), " [",
               formatC(ci[1], format='f', dig=3), ", ",
               formatC(ci[2], format='f', dig=3), "]")
ra <- cbind(roc.r$thresholds, roc.r$sensitivities, roc.r$specificities)
r <- coords(roc.r, x="best", best.method='closest.topleft', transpose = TRUE)
th <- r[1]
sens <- r[3]
spec <- r[2]
tab <- table(yts, prs)
acc <- sum(diag(tab))/sum(tab)
accs <- cbind(roc.r$sens, roc.r$spec)
perc <- which.min(abs(sort(probs)-th))/length(probs)

png(paste0(.res, "pp.png"), width=800, height=800)
plot(roc.r, cex.axis=1.2, cex.lab=1.1, col='black')
legend(x='bottomright', inset=0.01, lty=0,
       legend=c(paste0("AUC = ", aucc),
                "",
                paste0("Threshold. = ", formatC(th, format='f', dig=3)),
                paste0("Sens. = ", formatC(sens, format='f', dig=3)),
                paste0("Spec. = ", formatC(spec, format='f', dig=3)),
                paste0("Acc. = ", formatC(acc, format='f', dig=3))),
       bty='n', cex=1)
abline(h=sens, v=spec, lty=3, col='gray40', lwd=2)
dev.off()


#######################

s <- ndd$group%in%c("Mild", "Severe")
Mtr <- exprs(ndd)[, s]
ytr <- droplevels(ndd$group[s])

s <- ndp$group%in%c("Mild", "Severe")
Mts <- exprs(ndp)[, s]
yts <- droplevels(ndp$group[s])


rf <- randomForest(x=t(Mtr), y=ytr,
                   ntree=10000,
                   mtry=floor(sqrt(nrow(Mtr))),
                   strata=ytr, 
                   sampsize = ceiling(.8*ncol(Mtr)),
                   replace=FALSE, 
                   importance=TRUE,
                   do.trace=100)

prs <- predict(rf, t(Mts))
probs <- predict(rf, t(Mts), type='prob')

tab <- table(yts, prs )
acc <- sum(diag(tab))/sum(tab)
sens <- tab[2, 2]/sum(tab[2, ])
spec <- tab[1, 1]/sum(tab[1, ])
phi <- phi(tab, 3)
c(acc=acc, sens=sens, spec=spec, phi=phi)


roc.r <- roc(yts, probs[, 2])


auc <- as.numeric(roc.r$auc)
ci <- as.numeric(ci.auc(roc.r, method='bootstrap'))[c(1, 3)]
aucc <- paste0(formatC(auc, format='f', dig=3), " [",
               formatC(ci[1], format='f', dig=3), ", ",
               formatC(ci[2], format='f', dig=3), "]")
ra <- cbind(roc.r$thresholds, roc.r$sensitivities, roc.r$specificities)
r <- coords(roc.r, x="best", best.method='closest.topleft', transpose = TRUE)
th <- r[1]
sens <- r[3]
spec <- r[2]
tab <- table(yts, prs)
acc <- sum(diag(tab))/sum(tab)
accs <- cbind(roc.r$sens, roc.r$spec)
perc <- which.min(abs(sort(probs)-th))/length(probs)

png(paste0(.res, "pp.png"), width=800, height=800)
plot(roc.r, cex.axis=1.2, cex.lab=1.1, col='black')
legend(x='bottomright', inset=0.01, lty=0,
       legend=c(paste0("AUC = ", aucc),
                "",
                paste0("Threshold. = ", formatC(th, format='f', dig=3)),
                paste0("Sens. = ", formatC(sens, format='f', dig=3)),
                paste0("Spec. = ", formatC(spec, format='f', dig=3)),
                paste0("Acc. = ", formatC(acc, format='f', dig=3))),
       bty='n', cex=1)
abline(h=sens, v=spec, lty=3, col='gray40', lwd=2)
dev.off()






################################################################################################################
################################################################################################################
################################################################################################################


