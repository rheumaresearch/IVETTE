#################################################################################################################
#################################################################################################################
### P2.Multiv_pred_hosp.R
### Project: jcalvet_202110_marato
#################################################################################################################
#################################################################################################################
  # Libraries, paths and options

rm(list=ls())
.res <- paste0(.res0, "P.Multiv_predictor/P2.Multiv_pred_hosp_2/")
dir.create(.res, F, T)

#################################################################################################################
  # Prepare data

rm(list=ls())
ndd <- readRDS(paste0(.res0, "P.Multiv_predictor/P1.Data_setup/nd.dia.corr.eof.rds"))
ndp <- readRDS(paste0(.res0, "P.Multiv_predictor/P1.Data_setup/nd.prm.all.corr.calib.rds"))


#sort(apply(exprs(ndp), 1, function(o) sum(is.na(o))))
#do.call(rbind, tapply(exprs(ndp)["P05164", ], ndp$group,
#                      function(o) table(factor(as.numeric(is.na(o)), levels=0:1, labels=c("Value", "Missing")))))


ps <- c(names(which(apply(exprs(ndp), 1, function(o) sum(is.na(o)) == 0))),
        "P26038", "P29401")
ndd <- ndd[ps, ]
ndp <- ndp[ps, ]

ndp <- ndp[, apply(exprs(ndp), 2, function(o) all(!is.na(o)))]


Mtr <- exprs(ndd)
ytr <- factor(as.numeric(ndd$group%in%c("Severe", "Critical")),
              levels=0:1, labels=c("Not-hospitalized", "Hospitalized"))
Mts <- exprs(ndp)
yts <- factor(as.numeric(ndp$group%in%c("Severe", "Critical")),
            levels=0:1, labels=c("Not-hospitalized", "Hospitalized"))


names(ytr) <- sampleNames(ndd)
names(yts) <- sampleNames(ndp)


saveRDS(Mtr, paste0(.res, "Mtr.rds"))
saveRDS(Mts, paste0(.res, "Mts.rds"))
saveRDS(ytr, paste0(.res, "ytr.rds"))
saveRDS(yts, paste0(.res, "yts.rds"))


### Folds for training cv

set.seed(295483)
nfolds <- 10
id1 <- names(ytr[ytr==levels(ytr)[1]])
id2 <- names(ytr[ytr==levels(ytr)[2]])
id1 <- sample(id1)
id2 <- sample(id2)
idall <- c(id1, id2)
fd1 <- rep(1:nfolds, ceiling(length(id1)/nfolds))[1:length(id1)]
names(fd1) <- id1
fd2 <- rep(1:nfolds, ceiling(length(id2)/nfolds))[1:length(id2)]
names(fd2) <- id2
lfolds <- sapply(1:nfolds, function(j) names(c(fd1[fd1!=j], fd2[fd2!=j])))
lfolds <- sapply(lfolds, function(o)
{
    id1 <- sample(o[o%in%names(ytr)[ytr==levels(ytr)[1]]])
    id2 <- sample(o[o%in%names(ytr)[ytr==levels(ytr)[2]]])
    fd1 <- rep(1:nfolds, ceiling(length(id1)/nfolds))[1:length(id1)]
    names(fd1) <- id1
    fd2 <- rep(1:nfolds, ceiling(length(id2)/nfolds))[1:length(id2)]
    names(fd2) <- id2
    sapply(1:nfolds, function(j) names(c(fd1[fd1==j], fd2[fd2==j])))
}, simplify=F)

saveRDS(lfolds, paste0(.res, "lfods.tr.dia.rds"))

#################################################################################################################
  # Overall prediction - RF

rm(list=ls())
Mtr <- readRDS(paste0(.res, "Mtr.rds"))
Mts <- readRDS(paste0(.res, "Mts.rds"))
ytr <- readRDS(paste0(.res, "ytr.rds"))
yts <- readRDS(paste0(.res, "yts.rds"))
ndd <- readRDS(paste0(.res0, "P.Multiv_predictor/P1.Data_setup/nd.dia.corr.rds"))
ndp <- readRDS(paste0(.res0, "P.Multiv_predictor/P1.Data_setup/nd.prm.all.corr.calib.rds"))
source("code/functions/myfuns.R")


### Add sex and age 
ndd <- ndd[, names(ytr)]
ndp <- ndp[, names(yts)]
nvs <- c("sex","age")

d1 <- pData(ndd)[, nvs]
d1 <- d1[apply(d1, 1, function(o) all(!is.na(o))), ]
s1 <- match(rownames(d1),rownames(t(Mtr)))
d.Mtr <- cbind(d1,t(Mtr)[rownames(d1),])
identical(rownames(d.Mtr),names(ytr))

d2 <- pData(ndp)[, nvs]
d2 <- d2[apply(d2, 1, function(o) all(!is.na(o))), ]  
s2 <- match(rownames(d2),rownames(t(Mts)))
d.Mts <- cbind(d2,t(Mts)[rownames(d2),])
identical(rownames(d.Mts),names(yts))
# d2$y <- yts[rownames(d2)]



### Training -> DIA

set.seed(674932)
rf1 <- randomForest(x=d.Mtr, y=ytr,
                   ntree=10000,
                   mtry=floor(sqrt(ncol(d.Mtr))),
                   replace=FALSE, 
                   importance=TRUE,
                   do.trace=100)

RNGkind("L'Ecuyer-CMRG")
set.seed(38291)
r1 <- facc(ytr, rf1$predicted, rf1$votes[, 2], fn=paste0(.res, "ROC_hosp_rf_train_sexAge.png"))


### Test -> PRM

set.seed(38291)
rf2 <- randomForest(x=d.Mts, y=yts,
                   ntree=10000,
                   mtry=floor(sqrt(ncol(d.Mtr))),
                   replace=FALSE, 
                   importance=TRUE,
                   do.trace=100)

r2 <- facc(yts, rf2$predicted, rf2$votes[, 2],
          fn=paste0(.res, "ROC_hosp_rf_test_sexAge_1.png"))


RNGkind("L'Ecuyer-CMRG")
set.seed(284591)
r3 <- facc(yts, predict(rf1, d.Mts), predict(rf1, d.Mts, type='prob')[, 2],
          fn=paste0(.res, "ROC_hosp_rf_test_sexAge_2.png"))




#################################################################################################################
  # Overall prediction - RF - only validated proteins by PRM

rm(list=ls())
Mtr <- readRDS(paste0(.res, "Mtr.rds"))
Mts <- readRDS(paste0(.res, "Mts.rds"))
ytr <- readRDS(paste0(.res, "ytr.rds"))
yts <- readRDS(paste0(.res, "yts.rds"))
ndd <- readRDS(paste0(.res0, "P.Multiv_predictor/P1.Data_setup/nd.dia.corr.rds"))
ndp <- readRDS(paste0(.res0, "P.Multiv_predictor/P1.Data_setup/nd.prm.all.corr.calib.rds"))
source("code/functions/myfuns.R")


### Select proteins

ps <- c("P01011", "P02763", "P0DJI8",
        "P18428", "Q08830",
#        "P0DJI9",
        "Q96PD5", "P02654", "P04196",
        "P13796", "P02765","P02766", "P00738",
        "P26038",
        "P22105", "Q6UXB8", "P29401", "P33151","P15169", "P07358")
Mtr <- Mtr[ps, ]
Mts <- Mts[ps, ]


### Add sex and age 
ndd <- ndd[, names(ytr)]
ndp <- ndp[, names(yts)]
nvs <- c("sex","age")

d1 <- pData(ndd)[, nvs]
d1 <- d1[apply(d1, 1, function(o) all(!is.na(o))), ]
s1 <- match(rownames(d1),rownames(t(Mtr)))
d.Mtr <- cbind(d1,t(Mtr)[rownames(d1),])
identical(rownames(d.Mtr),names(ytr))

d2 <- pData(ndp)[, nvs]
d2 <- d2[apply(d2, 1, function(o) all(!is.na(o))), ]  
s2 <- match(rownames(d2),rownames(t(Mts)))
d.Mts <- cbind(d2,t(Mts)[rownames(d2),])
identical(rownames(d.Mts),names(yts))
# d2$y <- yts[rownames(d2)]


### Training -> DIA

set.seed(674932)
rf1 <- randomForest(x=d.Mtr, y=ytr,
                   ntree=10000,
                   mtry=floor(sqrt(ncol(d.Mtr))),
                   replace=FALSE, 
                   importance=TRUE,
                   do.trace=100)

RNGkind("L'Ecuyer-CMRG")
set.seed(38291)
r1 <- facc(ytr, rf1$predicted, rf1$votes[, 2], fn=paste0(.res, "ROC_hosp_rf_train_sexAge_valid.png"))


### Test -> PRM

set.seed(38291)
rf2 <- randomForest(x=d.Mts, y=yts,
                   ntree=10000,
                   mtry=floor(sqrt(ncol(d.Mts))),
                   replace=FALSE, 
                   importance=TRUE,
                   do.trace=100)

r2 <- facc(yts, rf2$predicted, rf2$votes[, 2],
          fn=paste0(.res, "ROC_hosp_rf_test_sexAge_valid_1.png"))


RNGkind("L'Ecuyer-CMRG")
set.seed(284591)
r3 <- facc(yts, predict(rf1, d.Mts), predict(rf1, d.Mts, type='prob')[, 2],
          fn=paste0(.res, "ROC_hosp_rf_test_sexAge_valid_2.png"))





################################################################################################################
################################################################################################################
################################################################################################################


# #################################################################################################################
# # Lasso selection - any PRM protein
# 
# rm(list=ls())
# Mtr <- readRDS(paste0(.res, "Mtr.rds"))
# Mts <- readRDS(paste0(.res, "Mts.rds"))
# ytr <- readRDS(paste0(.res, "ytr.rds"))
# yts <- readRDS(paste0(.res, "yts.rds"))
# lfolds <- readRDS(paste0(.res, "lfods.tr.dia.rds"))
# ndp <- readRDS(paste0(.res0, "P.Multiv_predictor/P1.Data_setup/nd.prm.all.corr.calib.rds"))
# source("code/functions/myfuns.R")
# 
# 
# ### Training -> DIA
# 
# ids <- colnames(Mtr)
# prs <- do.call(rbind, sapply(1:length(lfolds), function(j)
# {
#   lf <- lfolds[[j]]
#   idtr <- unique(unlist(lf))
#   idts <- ids[!ids%in%idtr]
#   fid <- unlist(sapply(1:length(lf), function(k) rep(k, length(lf[[k]])), simplify=F))
#   mcv <- cv.glmnet(x=t(Mtr[, idtr]), y=ytr[idtr], family = "binomial", foldid = fid)
#   m <- glmnet(x=t(Mtr[, idtr]), y=ytr[idtr], family = "binomial", lambda=mcv$lambda.1se)
#   prs <- predict(m, newx = t(Mtr[, idts]), family='binomial', lambda=mcv$lambda.1se, type='response')
# }))
# prs <- prs[,1]
# 
# tab <- table(ytr[names(prs)], prs > 0.5)
# r1 <- facc(ytr[names(prs)], factor(as.numeric(prs > 0.5)), prs,
#            fn=paste0(.res, "ROC_hosp_lasso_train.png"))
# 
# 
# fid <- sapply(ids, function(o) which(sapply(lfolds, function(l, o) !o%in%unlist(l), o)))
# mcv <- cv.glmnet(x=t(Mtr), y=ytr, family = "binomial", foldid = fid)
# m <- glmnet(x=t(Mtr), y=ytr, family = "binomial", lambda=mcv$lambda.1se)
# 
# bs <- as.matrix(coef(m))
# bs <- bs[bs[, 1] > 0, ,drop= F]
# prots <- fData(ndp)[rownames(bs), c("protein.id", "gene.symbol", "gene.name")]
# prots <- cbind(prots, beta=bs[rownames(prots), ])
# 
# sink(paste0(.res, "lasso_selected.proteins_all.prm.txt"))
# print(prots)
# sink()
# 
# 
# ### Test - PRM
# 
# prs <- predict(m, newx = t(Mts), family='binomial', lambda=mcv$lambda.1se, type='response')
# 
# source("code/functions/myfuns.R")
# rts <- facc(yts[rownames(prs)], factor(as.numeric(prs[, 1] > 0.5)), prs[, 1],
#             fn=paste0(.res, "ROC_hosp_lasso_test.png"))
# 
# #################################################################################################################
# # Lasso selection - proteins selected for ELISA
# 
# rm(list=ls())
# Mtr <- readRDS(paste0(.res, "Mtr.rds"))
# Mts <- readRDS(paste0(.res, "Mts.rds"))
# ytr <- readRDS(paste0(.res, "ytr.rds"))
# yts <- readRDS(paste0(.res, "yts.rds"))
# lfolds <- readRDS(paste0(.res, "lfods.tr.dia.rds"))
# ndp <- readRDS(paste0(.res0, "P.Multiv_predictor/P1.Data_setup/nd.prm.all.corr.calib.rds"))
# source("code/functions/myfuns.R")
# 
# 
# ### Select proteins
# 
# ps <- c(
#   #        "P26038",
#   "P01011", "P02763", "P0DJI8", "P18428", "Q08830",
#   #        "P0DJI9",
#   "Q96PD5",
#   "P02654", "P04196")
# Mtr <- Mtr[ps, ]
# Mts <- Mts[ps, ]
# 
# 
# ### Training -> DIA
# 
# ids <- colnames(Mtr)
# prs <- do.call(rbind, sapply(1:length(lfolds), function(j)
# {
#   lf <- lfolds[[j]]
#   idtr <- unique(unlist(lf))
#   idts <- ids[!ids%in%idtr]
#   fid <- unlist(sapply(1:length(lf), function(k) rep(k, length(lf[[k]])), simplify=F))
#   mcv <- cv.glmnet(x=t(Mtr[, idtr]), y=ytr[idtr], family = "binomial", foldid = fid)
#   m <- glmnet(x=t(Mtr[, idtr]), y=ytr[idtr], family = "binomial", lambda=mcv$lambda.1se)
#   prs <- predict(m, newx = t(Mtr[, idts]), family='binomial', lambda=mcv$lambda.1se, type='response')
# }))[, 1]
# 
# 
# tab <- table(ytr[names(prs)], prs > 0.5)
# r1 <- facc(ytr[names(prs)], factor(as.numeric(prs > 0.5)), prs,
#            fn=paste0(.res, "ROC_hosp_lasso_train_2.png"))
# 
# 
# fid <- sapply(ids, function(o) which(sapply(lfolds, function(l, o) !o%in%unlist(l), o)))
# mcv <- cv.glmnet(x=t(Mtr), y=ytr, family = "binomial", foldid = fid)
# m <- glmnet(x=t(Mtr), y=ytr, family = "binomial", lambda=mcv$lambda.1se)
# 
# bs <- as.matrix(coef(m))
# bs <- bs[bs[, 1] > 0, ,drop= F]
# prots <- fData(ndp)[rownames(bs), c("protein.id", "gene.symbol", "gene.name")]
# prots <- cbind(prots, beta=bs[rownames(prots), ])
# 
# sink(paste0(.res, "lasso_selected.proteins_elisa.only.txt"))
# print(prots)
# sink()
# 
# 
# ### Test - PRM
# 
# prs <- predict(m, newx = t(Mts), family='binomial', lambda=mcv$lambda.1se, type='response')
# 
# r2 <- facc(yts[rownames(prs)], factor(as.numeric(prs[, 1] > 0.5)), prs[, 1],
#            fn=paste0(.res, "ROC_hosp_lasso_test_2.png"))



#################################################################################################################
# Standard logistic regression selection - proteins selected for ELISA

# rm(list=ls())
# Mtr <- readRDS(paste0(.res, "Mtr.rds"))
# Mts <- readRDS(paste0(.res, "Mts.rds"))
# ytr <- readRDS(paste0(.res, "ytr.rds"))
# yts <- readRDS(paste0(.res, "yts.rds"))
# lfolds <- readRDS(paste0(.res, "lfods.tr.dia.rds"))
# ndp <- readRDS(paste0(.res0, "P.Multiv_predictor/P1.Data_setup/nd.prm.all.corr.calib.rds"))
# source("code/functions/myfuns.R")
# 
# 
# ### Select proteins
# 
# ps <- c(
#   #        "P26038",
#   "P01011", "P02763", "P0DJI8", "P18428", "Q08830",
#   #        "P0DJI9",
#   "Q96PD5",
#   "P02654", "P04196")
# Mtr <- Mtr[ps, ]
# Mts <- Mts[ps, ]
# 
# 
# ### Training -> DIA
# 
# ids <- colnames(Mtr)
# prs <- do.call(c, sapply(1:length(lfolds), function(j)
# {
#   lf <- lfolds[[j]]
#   idtr <- unique(unlist(lf))
#   idts <- ids[!ids%in%idtr]
#   form <- paste0("ytr[idtr] ~ ", paste0(ps, collapse=' + '))
#   m <- glm(as.formula(form), data=data.frame(t(Mtr[, idtr, drop=F])), family='binomial')
#   prs <- predict(m, newdata = data.frame(t(Mtr[, idts, drop=F])), type='response')
# }))
# 
# 
# tab <- table(ytr[names(prs)], prs > 0.5)
# r1 <- facc(ytr[names(prs)], factor(as.numeric(prs > 0.5)), prs,
#            fn=paste0(.res, "ROC_hosp_reglog_train.png"))
# 
# 
# form <- paste0("ytr ~ ", paste0(ps, collapse=' + '))
# m <- glm(as.formula(form), data=data.frame(t(Mtr)), family='binomial')
# m0 <- glm(ytr ~ 1, data=data.frame(t(Mtr)), family='binomial')
# mst <- step(m0, scope=list(lower=m0, upper=m), k=0.7)
# mst
# 
# #mm <- glm(ytr ~ P01011 + P02654 + P04196, data=data.frame(t(Mtr)), family='binomial')
# 
# 
# ### Test - PRM
# 
# prs <- predict(mst, newdata = data.frame(t(Mts)), type='response')
# r2 <- facc(yts[names(prs)], factor(as.numeric(prs > 0.5)), prs,
#            fn=paste0(.res, "ROC_hosp_reglog_test.png"))



################################################################################################################
# Classic parameters - RF

rm(list=ls())
ndd <- readRDS(paste0(.res0, "P.Multiv_predictor/P1.Data_setup/nd.dia.corr.eof.rds"))
ndp <- readRDS(paste0(.res0, "P.Multiv_predictor/P1.Data_setup/nd.prm.all.corr.calib.rds"))
ytr <- readRDS(paste0(.res, "ytr.rds"))
yts <- readRDS(paste0(.res, "yts.rds"))
lfolds <- readRDS(paste0(.res, "lfods.tr.dia.rds"))
source("code/functions/myfuns.R")

ndd <- ndd[, names(ytr)]
ndp <- ndp[, names(yts)]

nvs <- c(#"dm", "obesity", "dlp", "aht",
  "lymphocytes", "neutrophils", "leukocytes", "ferritin", "pcr", "dd", "ldh",
  "center")
d1 <- pData(ndd)[, nvs]
d2 <- pData(ndp)[, nvs]
d1$set <- "DIA"
d2$set <- "PRM"
d <- rbind(d1, d2)
d$set <- factor(d$set)
y <- c(ytr, yts)
d$y <- y[rownames(d)]

d <- d[apply(d, 1, function(o) all(!is.na(o))), ]
d$center <- NULL
d$set <- NULL

summary(d)


### RF model

M <- d[, colnames(d)!="y"]
y <- d$y
names(y) <- rownames(M)

set.seed(72911)
rf1 <- randomForest(x=M, y=y,
                    ntree=10000,
                    mtry=floor(sqrt(ncol(M))),
                    strata=y,
                    sampsize=ceiling(rep(min(summary(y)), 2)*0.8),
                    replace=FALSE, 
                    importance=TRUE,
                    do.trace=100)


rf1$importance[order(rf1$importance[, 3], decreasing=T), ]

RNGkind("L'Ecuyer-CMRG")
set.seed(38291)
r1 <- facc(y, rf1$predicted, rf1$votes[, 2], fn=paste0(.res, "ROC_hosp_rf_parmclass.png"))


### RF model - add PRM proteins

ps <-  c("P01011", "P04196", "P18428")
for (o in ps)
{
  x1 <- exprs(ndd)[o, ]
  names(x1) <- sampleNames(ndd)
  x2 <- exprs(ndp)[o, ]
  names(x2) <- sampleNames(ndp)
  x <- c(x1, x2)
  d[, o] <- x[rownames(d)]
}

M <- d[, colnames(d)!="y"]
y <- d$y
names(y) <- rownames(M)

set.seed(72911)
rf2 <- randomForest(x=M, y=y,
                    ntree=10000,
                    mtry=floor(sqrt(ncol(M))),
                    strata=y,
                    sampsize=ceiling(rep(min(summary(y)), 2)*0.8),
                    replace=FALSE, 
                    importance=TRUE,
                    do.trace=100)


RNGkind("L'Ecuyer-CMRG")
set.seed(38922)
r2 <- facc(y, rf2$predicted, rf2$votes[, 2], fn=paste0(.res, "ROC_hosp_rf_parmclass_prots.png"))

rf1$importance[order(rf1$importance[, 3], decreasing=T), ]
rf2$importance[order(rf2$importance[, 3], decreasing=T), ]

################################################################################################################
# Classic parameters - logistic regression

# rm(list=ls())
# ndd <- readRDS(paste0(.res0, "P.Multiv_predictor/P1.Data_setup/nd.dia.corr.eof.rds"))
# ndp <- readRDS(paste0(.res0, "P.Multiv_predictor/P1.Data_setup/nd.prm.all.corr.calib.rds"))
# ytr <- readRDS(paste0(.res, "ytr.rds"))
# yts <- readRDS(paste0(.res, "yts.rds"))
# lfolds <- readRDS(paste0(.res, "lfods.tr.dia.rds"))
# source("code/functions/myfuns.R")
# 
# ndd <- ndd[, names(ytr)]
# ndp <- ndp[, names(yts)]
# 
# nvs <- c(#"dm", "obesity", "dlp", "aht",
#   "lymphocytes", "neutrophils", "leukocytes", "ferritin", "pcr", "dd", "ldh")
# d1 <- pData(ndd)[, nvs]
# d2 <- pData(ndp)[, nvs]
# d1$set <- "DIA"
# d2$set <- "PRM"
# d <- rbind(d1, d2)
# d$set <- factor(d$set)
# y <- c(ytr, yts)
# d$y <- y[rownames(d)]
# 
# d <- d[apply(d, 1, function(o) all(!is.na(o))), ]
# d$center <- NULL
# d$set <- NULL
# 
# summary(d)
# 
# 
# ### Folds
# 
# set.seed(729291)
# nfolds <- 10
# id1 <- rownames(d)[d$y==levels(d$y)[1]]
# id2 <- rownames(d)[d$y==levels(d$y)[2]]
# id1 <- sample(id1)
# id2 <- sample(id2)
# idall <- c(id1, id2)
# fd1 <- rep(1:nfolds, ceiling(length(id1)/nfolds))[1:length(id1)]
# names(fd1) <- id1
# fd2 <- rep(1:nfolds, ceiling(length(id2)/nfolds))[1:length(id2)]
# names(fd2) <- id2
# lfolds <- sapply(1:nfolds, function(j) names(c(fd1[fd1!=j], fd2[fd2!=j])))
# 
# 
# ### Logistic regression
# 
# 
# M <- d[, colnames(d)!="y"]
# y <- d$y
# names(y) <- rownames(M)
# 
# ids <- rownames(M)
# probs <- do.call(c, sapply(1:length(lfolds), function(j)
# {
#   idtr <- lfolds[[j]]
#   idts <- ids[!ids%in%idtr]
#   form <- paste0("y[idtr] ~ ", paste0(nvs, collapse=' + '))
#   m <- glm(as.formula(form), data=M[idtr, , drop=F], family='binomial')
#   prs <- predict(m, newdata = M[idts, , drop=F], type='response')
# }))
# 
# co <- (summary(y)/length(y))[2]
# prs <- factor(as.numeric(probs > co), levels=0:1, labels=levels(y))
# 
# RNGkind("L'Ecuyer-CMRG")
# set.seed(95022)
# r1 <- facc(y, prs, probs, fn=paste0(.res, "ROC_hosp_reglog_parmclass.png"), th=co)
# 
# 
# ### Add PRM proteins
# 
# ps <-  c("P01011", "P04196", "P18428")
# for (o in ps)
# {
#   x1 <- exprs(ndd)[o, ]
#   names(x1) <- sampleNames(ndd)
#   x2 <- exprs(ndp)[o, ]
#   names(x2) <- sampleNames(ndp)
#   x <- c(x1, x2)
#   d[, o] <- x[rownames(d)]
# }
# 
# M <- d[, colnames(d)!="y"]
# y <- d$y
# names(y) <- rownames(M)
# 
# ids <- rownames(M)
# probs <- do.call(c, sapply(1:length(lfolds), function(j)
# {
#   idtr <- lfolds[[j]]
#   idts <- ids[!ids%in%idtr]
#   form <- paste0("y[idtr] ~ ", paste0(c(nvs, ps), collapse=' + '))
#   m <- glm(as.formula(form), data=M[idtr, , drop=F], family='binomial')
#   prs <- predict(m, newdata = M[idts, , drop=F], type='response')
# }))
# 
# co <- (summary(y)/length(y))[2]
# prs <- factor(as.numeric(probs > co), levels=0:1, labels=levels(y))
# 
# RNGkind("L'Ecuyer-CMRG")
# set.seed(95022)
# r2 <- facc(y, prs, probs, fn=paste0(.res, "ROC_hosp_reglog_parmclass_prots.png"), th=co)
