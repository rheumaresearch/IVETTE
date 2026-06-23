#################################################################################################################
#################################################################################################################
### O.Concordance_DIA_PRM_complete/O3.Diff_exprs_hosp_DIA_sel.R
### Project: jcalvet_202110_marato
#################################################################################################################
#################################################################################################################
  # Libraries, paths and options

rm(list=ls())
.res <- paste0(.res0, "O.Concordance_DIA_PRM_complete/O3.Diff_exprs_hosp_DIA_sel/")
dir.create(.res, F, T)

#################################################################################################################
  # Select data

rm(list=ls())
nd <- readRDS(paste0(.dat, "1_DIA/nd.fp.rds"))
nd2 <- readRDS(paste0(.dat, "4_Concordance_DIA.PRM/nd.prm.all.rds"))

nd$group2 <- factor(as.numeric(nd$group%in%c("Severe", "Critical")),
                    levels=0:1, labels=c("Not.hospitalized", "Hospitalized"))

samps <- intersect(sampleNames(nd), sampleNames(nd2))
prots <- intersect(featureNames(nd), featureNames(nd2))

nd <- nd[prots, samps]

s <- sapply(1:nrow(nd), function(j)
{
    y <- exprs(nd)[j, ]
    nval <- tapply(y, nd$group2, function(o) sum(!is.na(o)))
    sum(nval >= 3) >= 2
})

nds <- nd[s, ]

saveRDS(nds, paste0(.res, "nd.log.filt.rds"))

#################################################################################################################
  # Contrast matrix

rm(list=ls())

conts <- c("Hospitalized-Not.hospitalized"="group2Hospitalized - group2Not.hospitalized = 0")


saveRDS(conts, paste0(.res, "conts.rds"))

#################################################################################################################
  # Comparisons

rm(list=ls())
conts <- readRDS( paste0(.res, "conts.rds"))
nd <- readRDS(paste0(.res, "nd.log.filt.rds"))


levg <- levels(nd$group2)

#i <- which(rownames(nd)=='P05556')
#i <- 150

form <- "y ~ -1 + group2 + quant.batch"
lr <- mclapply(1:nrow(nd), function(i)
{

    d <- cbind(pData(nd), y=exprs(nd)[i, ])


    ### Categorical models for all

    m <- try(lm(as.formula(form), data=d), silent=T)
    
    if (class(m)!='try-error' && all(!is.na(summary(m)$coef[, 2])))
    {

        bs <- coef(m)
        nbs <- names(coef(m))
        nva <- c("quant.batch", "pch1", "col.hemol")
        nva.in <- sapply(nva, function(o, form) any(regexpr(o, nbs) > 0), nbs)
        if (any(nva.in))
        {
            nvas <- names(which(nva.in))
            fmm0 <- paste0(" ~ ", paste0(nvas, collapse=' + '))
            dmm <- data.frame(model.matrix(as.formula(fmm0), data=d))[, -1, drop=F]
            dmm <- dmm[, colnames(dmm)%in%names(bs), drop=F]
            fmm1 <- fmm0
            dmm1 <- model.matrix(as.formula(fmm1), data=pData(nd))[, -1, drop=F]
            dmm1 <- dmm1[, colnames(dmm)]
            for (o in colnames(dmm1)) dmm[, o] <- dmm[, o] - mean(dmm1[, o])
            dmm2 <- as.matrix(dmm)
            fmm <- form
            fmm <- gsub("quant.batch", paste0(colnames(dmm)[regexpr("^quant.batch", colnames(dmm)) > 0], collapse=' + '), fmm)
            fmm <- gsub("col.hemol", paste0(colnames(dmm)[regexpr("^col.hemol", colnames(dmm)) > 0], collapse=' + '), fmm)            
            if (!"pch1"%in%nvas)
                fmm <- gsub(" \\+ pch1", "", fmm)
            if (!"col.hemol"%in%nvas)
                fmm <- gsub(" \\+ col.hemol", "", fmm)
            if (!"quant.batch"%in%nvas)
                fmm <- gsub(" \\+ quant.batch", "", fmm)
            dmm[, "y"] <- d[rownames(dmm), "y"]
            dmm[, "group2"] <- d[rownames(dmm), "group2"]
            dmm[, "patid"] <- d[rownames(dmm), "patid"]
            dmm[, "peptide"] <- d[rownames(dmm), "peptide"]
            mm <- lm(as.formula(fmm), data=dmm)
        }else
        {
            mm <- m
        }

        ### Adjusted values

        dmm2 <- as.matrix(dmm2)
        ycorr <- d$y - as.numeric(dmm2%*%as.matrix(coef(mm)[colnames(dmm2)]))
        names(ycorr) <- rownames(d)

        
        ### Group means
        
        d <<- d

        sm <- summary(mm)$coef
        dd1 <- cbind(sm[regexpr("^group2", rownames(sm)) > 0, c("Estimate", "Std. Error")])
        dd1 <- data.frame(group=gsub("^group2", "", rownames(dd1)), dd1, stringsAsFactors=F)
        colnames(dd1) <- c("group2", "log2mean", "log2se")

        contg <- paste0(paste0(rownames(dd1), "/", nrow(dd1) , collapse=' + '), " = 0")
        gl <- try(glht(mm, linfct = contg), silent=T)

        if (class(gl)!="try-error")
        {    
            sm0 <- summary(gl, test=adjusted("none"))
            dd0 <- c(sm0$test$coef, sm0$test$sigma)
            names(dd0) <- c("global.mean", "global.se")

            ### Contrasts

            levr <- levg[!paste0("group2", levg)%in%nbs]
            whs <- which(sapply(conts, function(o)
            {
                !any(sapply(levr, function(l, o) regexpr(l, o) > 0, o))
            }))
            contss <- conts[whs]
            gl <- try(glht(m, linfct = contss), silent=T)
            if (class(gl)!='try-error')
            {
                
                test <- summary(gl, test=adjusted("none"))
                ci <- as.data.frame(confint(test)$confint)
                pvs <- as.numeric(test$test$pvalues)
                se <- test$test$sigma
                dd2 <- cbind(ci[, 1], se, test$test$tstat, pvs)
                colnames(dd2) <- c("log2fc", "log2se", "t", "pv")

                ### AIC

                rm(test)
                rm(gl)
                rm(d)
                
                ### Results
                
                r <- list(ycorr=ycorr, sm1=dd1, sm2=dd2, sm0=dd0, form=form, m=m)
            }else
            {
                r <- NA
            }
        }else
        {
            r <- NA
        }
    }else
    {
        r <- NA
    }
    r
}, mc.cores=24)

table(sapply(lr,  class))

#which(sapply(lr,  class)=='try-error')[1]
#which(sapply(lr,  class)=='logical')[1]

names(lr) <- featureNames(nd)


lr <- lr[!is.na(lr)]

saveRDS(lr, paste0(.res, "lmod.group.rds"))

#################################################################################################################
  # Extract group means

rm(list=ls())
li <- readRDS(paste0(.res, "lmod.group.rds"))
nd <- readRDS(paste0(.res, "nd.log.filt.rds"))


lr <- sapply(li, function(o) o$sm1, simplify=F)
dr <- t(sapply(lr, function(o)
{
    pp <- do.call(c, sapply(levels(nd$group2), function(l, o)
    {
        o <- data.frame(o)
        ll <- paste0("group2", l)
        r <- as.numeric(c(o[ll, 2], o[ll, 3]))
        names(r) <- paste0(l, ".", colnames(o)[2:3])
        r
    }, o, simplify=F, USE.NAMES=F))
}))


saveRDS(dr, paste0(.res, "dmeans.rds"))

#################################################################################################################
  # Extract test results

rm(list=ls())
li <- readRDS(paste0(.res, "lmod.group.rds"))
nd <- readRDS(paste0(.res, "nd.log.filt.rds"))
conts <- readRDS( paste0(.res, "conts.rds"))

lr <- sapply(li, function(o) o$sm2, simplify=F)

dr <- t(sapply(lr, function(o)
{
    r <- do.call(c, sapply(1:length(conts), function(j, o)
    {
        o <- data.frame(o)
        o$fc <- 2^abs(o$log2fc)*sign(o$log2fc)
        o <- o[, c(1:3, 5, 4)]
        l <- conts[j]
        ll <- gsub(" \\= 0", "", l)
        r <- as.numeric(o[ll, ])
        names(r) <- paste0(names(conts)[j], ".", colnames(o))
        r
    }, o, simplify=F, USE.NAMES=F))
}))

### Multitesting

npv <- colnames(dr)[regexpr("\\.pv$", colnames(dr)) > 0]
comps <- gsub("\\.pv$", "", npv)

dr <- data.frame(dr, check.names=F)

for (o in npv)
{
    dr[, gsub("\\.pv", ".adj.pv", o)] <- p.adjust(dr[, o], "BH")
}

nv <- paste0(rep(comps, each=6), c(".log2fc", ".log2se", ".t", ".fc", ".pv", ".adj.pv"))
dr <- dr[, nv]

saveRDS(dr, paste0(.res, "dtests.rds"))

#################################################################################################################
  # Extract global means

rm(list=ls())
li <- readRDS(paste0(.res, "lmod.group.rds"))
nd <- readRDS(paste0(.res, "nd.log.filt.rds"))
conts <- readRDS( paste0(.res, "conts.rds"))

dr <- t(sapply(li, function(o) o$sm0))

saveRDS(dr, paste0(.res, "dmeans.glob.rds"))

#################################################################################################################
  #  Merge results

rm(list=ls())
dm <- readRDS(paste0(.res, "dmeans.rds"))
dt <- readRDS(paste0(.res, "dtests.rds"))
nd <- readRDS(paste0(.res, "nd.log.filt.rds"))


### All results

dr <- fData(nd)
rownames(dr) <- dr$protein.group
dr <- cbind(dr, dm[match(rownames(dr), rownames(dm)), ])
dr <- cbind(dr, dt[match(rownames(dr), rownames(dt)), ])


saveRDS(dr, paste0(.res, "dcomp.rds"))
write.table(dr, paste0(.res, "dgroup.results.csv"), sep=';', row.names=T)

################################################################################################################
################################################################################################################
################################################################################################################


