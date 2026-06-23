#################################################################################################################
#################################################################################################################
### G.Discovery_DIA_FP_MaxLFQ/G2.Diff_expression.R
### Project: jcalvet_202110_marato
#################################################################################################################
#################################################################################################################
  # Libraries, paths and options

rm(list=ls())
.res <- paste0(.res0, "G.Discovery_DIA_FP_MaxLFQ/G2.Diff_expression/")
dir.create(.res, F, T)

#################################################################################################################
  # Remove proteins with value in less than a minimum o samples per group  and center adjusting variables

rm(list=ls())
nd <- readRDS(paste0(.dat, "1_DIA/nd.fp.rds"))

s <- sapply(1:nrow(nd), function(j)
{
    y <- exprs(nd)[j, ]
    nval <- tapply(y, nd$group, function(o) sum(!is.na(o)))
    sum(nval >= 3) >= 2
})

nds <- nd[s, ]

saveRDS(nds, paste0(.res, "nd.log.filt.rds"))

#################################################################################################################
  # Contrast matrix

rm(list=ls())

conts <- c("Mild-Asymptomatic"="groupMild - groupAsymptomatic = 0",
           "Severe-Mild"="groupSevere - groupMild = 0",
           "Critical-Severe"="groupCritical - groupSevere = 0",
           "Hospitalized-Not.hospitalized"="(groupCritical + groupSevere)/2 - (groupMild + groupAsymptomatic)/2 = 0")


### Matrix for merged groups

conts.mean <- c("Hospitalized"="(groupCritical + groupSevere)/2  = 0",
                "Not.hospitalized"="(groupMild + groupAsymptomatic)/2 = 0")


saveRDS(conts, paste0(.res, "conts.rds"))
saveRDS(conts.mean, paste0(.res, "conts.means.rds"))

#################################################################################################################
  # Comparisons

rm(list=ls())
conts <- readRDS( paste0(.res, "conts.rds"))
conts.means <- readRDS( paste0(.res, "conts.means.rds"))
nd <- readRDS(paste0(.res, "nd.log.filt.rds"))

nd$groupn <- as.numeric(nd$group)
nd$group2 <- factor(as.numeric(nd$group%in%c("Severe", "Critical")), levels=0:1,
                    labels=c("Not-hospitalized", "Hospitalized"))


levg <- levels(nd$group)

#i <- which(rownames(nd)=='P05556')
#i <- 150

form <- "y ~ -1 + group + quant.batch"
saveRDS(form, paste0(.res, "form.rds"))

lr <- mclapply(1:nrow(nd), function(i)
{

    d <- cbind(pData(nd), y=exprs(nd)[i, ])

    ### Best model -> only for labelling
    
    aicg <- try(AIC(lm(y ~ group + quant.batch, data=d)), silent=T)
    aic2 <- try(AIC(lm(y ~ group2 + quant.batch, data=d)), silent=T)
    aicn <- try(AIC(lm(y ~ groupn + quant.batch, data=d)), silent=T)

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
            dmm[, "group"] <- d[rownames(dmm), "group"]
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
        dd1 <- cbind(sm[regexpr("^group", rownames(sm)) > 0, c("Estimate", "Std. Error")])
        dd1 <- data.frame(group=gsub("^group", "", rownames(dd1)), dd1, stringsAsFactors=F)
        colnames(dd1) <- c("group", "log2mean", "log2se")

        contg <- paste0(paste0(rownames(dd1), "/", nrow(dd1) , collapse=' + '), " = 0")
        gl <- try(glht(mm, linfct = contg), silent=T)

        gl2 <- try(glht(mm, linfct = conts.means), silent=T)
        
        if (class(gl2)!='try-error')
        {
            test <- summary(gl2, test=adjusted("none"))
            dd1b <- data.frame(group=names(conts.means), test$test$coef, test$test$sigma, stringsAsFactors=F)
        }else
        {
            dd1b <- rbind(rep(NA, 3), rep(NA, 3))
        }
        colnames(dd1b) <- c("group", "log2mean", "log2se")
        rownames(dd1b) <- names(conts.means)
       
        if (class(gl)!="try-error")
        {    
            sm0 <- summary(gl, test=adjusted("none"))
            dd0 <- c(sm0$test$coef, sm0$test$sigma)
            names(dd0) <- c("global.mean", "global.se")

            ### Contrasts

            levr <- levg[!paste0("group", levg)%in%nbs]
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
                
                r <- list(ycorr=ycorr, sm1=dd1, sm1b=dd1b, sm2=dd2, sm0=dd0, aic=aicg, aic2=aic2, aicn=aicn,
                          form=form, m=m)
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
    pp <- do.call(c, sapply(levels(nd$group), function(l, o)
    {
        o <- data.frame(o)
        ll <- paste0("group", l)
        r <- as.numeric(c(o[ll, 2], o[ll, 3]))
        names(r) <- paste0(l, ".", colnames(o)[2:3])
        r
    }, o, simplify=F, USE.NAMES=F))
}))

lr <- sapply(li, function(o) o$sm1b, simplify=F)
dr2 <- t(sapply(lr, function(o)
{
    pp <- do.call(c, sapply(c("Hospitalized", "Not.hospitalized"), function(l, o)
    {
        o <- data.frame(o)
        r <- as.numeric(c(o[l, 2], o[l, 3]))
        names(r) <- paste0(l, ".", colnames(o)[2:3])
        r
    }, o, simplify=F, USE.NAMES=F))
}))


dr <- cbind(dr, dr2[rownames(dr), ])

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
  # Model type according to AIC

rm(list=ls())
li <- readRDS(paste0(.res, "lmod.group.rds"))
nd <- readRDS(paste0(.res, "nd.log.filt.rds"))
conts <- readRDS( paste0(.res, "conts.rds"))

dr <- t(sapply(li, function(o) c(o$aic, o$aic2, o$aicn), simplify=T))
colnames(dr) <- c("categ", "2groups", "ordinal")
dr <- data.frame(dr, check.names=F)

dr$model.aic <- colnames(dr)[apply(dr, 1, which.min)]

saveRDS(dr, paste0(.res, "dmodel.aic.rds"))    

#################################################################################################################
  # Extract global means

rm(list=ls())
li <- readRDS(paste0(.res, "lmod.group.rds"))
nd <- readRDS(paste0(.res, "nd.log.filt.rds"))
conts <- readRDS( paste0(.res, "conts.rds"))

dr <- t(sapply(li, function(o) o$sm0))

saveRDS(dr, paste0(.res, "dmeans.glob.rds"))

#################################################################################################################
  # Matrix of corrected expression values at protein level

rm(list=ls())
li <- readRDS(paste0(.res, "lmod.group.rds"))
nd <- readRDS(paste0(.res, "nd.log.filt.rds"))


M <- t(sapply(li, function(o) o$ycorr, simplify=T))

length(unique(fData(nd)$protein.group))
df <- unique(fData(nd)[, c("protein.group", "protein.name", "gene.symbol")])
rownames(df) <- df$protein.group

ps <- intersect(rownames(df), rownames(M))
df <- df[ps, ]

M <- as.matrix(M[ps, sampleNames(nd)])
nda <- new('ExpressionSet', exprs=M, featureData=new('AnnotatedDataFrame', df),
                      phenoData=new('AnnotatedDataFrame', pData(nd)), annotation="PD-MaxLFQprotein_corrected")


saveRDS(nda, paste0(.res, "nd.log.prot.corr.rds"))

#################################################################################################################
  #  Merge results

rm(list=ls())
dm <- readRDS(paste0(.res, "dmeans.rds"))
dmg <- readRDS(paste0(.res, "dmeans.glob.rds"))
dt <- readRDS(paste0(.res, "dtests.rds"))
da <- readRDS(paste0(.res, "dmodel.aic.rds"))
nd <- readRDS(paste0(.res, "nd.log.prot.corr.rds"))


### All results

dr <- fData(nd)
rownames(dr) <- dr$protein.group
dr <- cbind(dr, dmg[match(rownames(dr), rownames(dmg)), ])
dr <- cbind(dr, dm[match(rownames(dr), rownames(dm)), ])
dr <- cbind(dr, dt[match(rownames(dr), rownames(dt)), ])
dr <- cbind(dr, "model.aic"=da[match(rownames(dr), rownames(da)), "model.aic"])


saveRDS(dr, paste0(.res, "dcomp.rds"))
write.table(dr, paste0(.res, "dgroup.results.csv"), sep=';', row.names=T)

#################################################################################################################
  # Create tables for report

rm(list=ls())
dc <- readRDS(paste0(.res, "dcomp.rds"))
conts <- readRDS( paste0(.res, "conts.rds"))

ncomp <- colnames(dc)[regexpr("\\.adj\\.pv$", colnames(dc)) > 0]
ncomp <- gsub("\\.adj\\.\\.pv$", "",ncomp)

vsel <- c("protein.group", "protein.name", "gene.symbol",         
          "global.mean", "global.se")                          

dc$plot <- ""
dc$"plot zoom" <- ""

lr <- sapply(1:length(conts), function(i)
      {
          o <- conts[i]
          if (i!=4)
          {
              ngr <- unique(gsub("^group", "", strsplit(o, split=' [\\-] ')[[1]]))
              ngr <- unique(gsub("^group", "", unlist(strsplit(ngr, split=' [\\+] '))))
              ngr <- rev(gsub(" \\= 0", "", ngr))
          }else
          {
              ngr <- c("Hospitalized", "Not.hospitalized")
          }
          nv1 <- paste0(rep(ngr, each=2), c(".log2mean", ".log2se"))          
          nv2 <- paste0(rep(names(conts[i]), each=3), c(".fc", ".pv", ".adj.pv"))
          ds <- dc[, c(vsel, nv2, "model.aic", "plot", "plot zoom", nv1)]
          ds <- ds[order(ds[, paste0(names(conts)[i], ".pv")]), ]
          ds
      }, simplify=F)

names(lr) <- names(conts)

saveRDS(lr, paste0(.res, "lcomp.rds"))

#################################################################################################################
  # Significant proteins by thresholds

rm(list=ls())
d <- readRDS(paste0(.res, "dcomp.rds"))

### Comparisons

nadjp <- colnames(d)[regexpr("\\.adj\\.pv", colnames(d))>0]
nfc <- colnames(d)[regexpr("\\.log2fc", colnames(d))>0]
comps <- gsub("\\.adj\\.pv", "", nadjp)


### Thresholds

th.pv <- c(0.01, 0.05, 0.10, 0.20, 0.25)
th.fc <- c(1, 1.20, 1.25, 1.5, 1.75, 2, 3, 5, 10)


### Compute number of significant genes by thresholds


lr <- sapply(comps, function(o1)
{
    r <- sapply(th.fc, function(o2, o1)
    {
        sapply(th.pv, function(o3, o2, o1)
        {
            ps <- rownames(d[!is.na(d[, paste0(o1, ".fc")]) & abs(d[, paste0(o1, ".fc")]) > o2 &
                             !is.na(d[, paste0(o1, ".adj.pv")]) & d[, paste0(o1, ".adj.pv")] < o3, ])
            ps.u <- rownames(d[!is.na(d[, paste0(o1, ".fc")]) & d[, paste0(o1, ".fc")] > o2 &
                               !is.na(d[, paste0(o1, ".adj.pv")]) & d[, paste0(o1, ".adj.pv")] < o3, ])
            ps.d <- rownames(d[!is.na(d[, paste0(o1, ".fc")]) & d[, paste0(o1, ".fc")] < -o2 &
                               !is.na(d[, paste0(o1, ".adj.pv")]) & d[, paste0(o1, ".adj.pv")] < o3, ])
            nps <- length(ps)
            nps.u <- length(ps.u)
            nps.d <- length(ps.d)
            paste0(nps, "<br>", nps.u, "<br>", nps.d)
        }, o2, o1)
    }, o1)
    colnames(r) <- paste("absFC > ", formatC(th.fc, format='f', dig=2), "<br>Up<br>Down")
    rownames(r) <- paste("FDR < ", formatC(th.pv, format='f', dig=3))
    r <- rbind(colnames(r), r)
    r <- cbind(rownames(r), r)
    r
}, simplify=F)

saveRDS(lr, paste0(.res, "ln.sign.rds"))


### Save as html for report

p <- openPage(paste0(.res, 'Nb_sign_by_thresholds.html'))
hwrite(paste('<br>Number of significant proteins by thresholds<br><br>'), p)
for (j in 1:length(lr))
{
    hwrite(paste('<br><br><br>', names(lr)[j],'<br><br>'), p)
    hwrite(lr[[j]], page=p, center=F, row.names=F, col.names=F,
           col.style=c("text-align:center", rep("text-align:center", ncol(lr[[j]])-1)),
           col.width=c(100, rep(125, ncol(lr[[j]])-1)))
}
closePage(p)

#################################################################################################################
  # Significant proteins by thresholds - raw p-value

rm(list=ls())
d <- readRDS(paste0(.res, "dcomp.rds"))

### Comparisons

nadjp <- colnames(d)[regexpr("\\.adj\\.pv", colnames(d))>0]
nfc <- colnames(d)[regexpr("\\.log2fc", colnames(d))>0]
comps <- gsub("\\.adj\\.pv", "", nadjp)


### Thresholds

th.pv <- c(0.0001, 0.001, 0.01, 0.05, 1)
th.fc <- c(1, 1.20, 1.25, 1.5, 1.75, 2, 3, 5, 10)


### Compute number of significant genes by thresholds


lr <- sapply(comps, function(o1)
{
    r <- sapply(th.fc, function(o2, o1)
    {
        sapply(th.pv, function(o3, o2, o1)
        {
            ps <- rownames(d[!is.na(d[, paste0(o1, ".fc")]) & abs(d[, paste0(o1, ".fc")]) > o2 &
                             !is.na(d[, paste0(o1, ".pv")]) & d[, paste0(o1, ".pv")] < o3, ])
            ps.u <- rownames(d[!is.na(d[, paste0(o1, ".fc")]) & d[, paste0(o1, ".fc")] > o2 &
                               !is.na(d[, paste0(o1, ".pv")]) & d[, paste0(o1, ".pv")] < o3, ])
            ps.d <- rownames(d[!is.na(d[, paste0(o1, ".fc")]) & d[, paste0(o1, ".fc")] < -o2 &
                               !is.na(d[, paste0(o1, ".pv")]) & d[, paste0(o1, ".pv")] < o3, ])
            nps <- length(ps)
            nps.u <- length(ps.u)
            nps.d <- length(ps.d)
            paste0(nps, "<br>", nps.u, "<br>", nps.d)
        }, o2, o1)
    }, o1)
    colnames(r) <- paste("absFC > ", formatC(th.fc, format='f', dig=2), "<br>Up<br>Down")
    rownames(r) <- paste("raw p-value < ", formatC(th.pv, format='f', dig=4))
    r <- rbind(colnames(r), r)
    r <- cbind(rownames(r), r)
    r
}, simplify=F)

saveRDS(lr, paste0(.res, "ln.sign.raw.rds"))


### Save as html for report

p <- openPage(paste0(.res, 'Nb_sign_by_thresholds_raw.pvalue.html'))
hwrite(paste('<br>Number of significant proteins by thresholds - raw p-value<br><br>'), p)
for (j in 1:length(lr))
{
    hwrite(paste('<br><br><br>', names(lr)[j],'<br><br>'), p)
    hwrite(lr[[j]], page=p, center=F, row.names=F, col.names=F,
           col.style=c("text-align:center", rep("text-align:center", ncol(lr[[j]])-1)),
           col.width=c(100, rep(125, ncol(lr[[j]])-1)))
}
closePage(p)

################################################################################################################
  # Stripcharts

rm(list=ls())
d <- readRDS(paste0(.res, "dcomp.rds"))
nd <- readRDS(paste0(.res, "nd.log.filt.rds"))
source("code/functions/make_transparent.R")


#system(paste0("rm -R ", .res, 'Stripcharts'))
pss <- unique(rownames(nd))
dir.create(paste0(.res, 'Stripcharts'), F)
cols <- c("deepskyblue3", "orange3", "blueviolet", "red3")
gr <- pData(nd)$group
levels(gr) <- gsub("\\:", "\\\n", levels(gr))
#pat <- pData(nd)$num.estudi
#patu <- unique(pat)

dir.create(paste0(.res, 'Stripcharts_zoom'), F)
pp <- mclapply(pss, function(ps)
{
    o <- exprs(nd)[ps, ]
    g <- ps
    g2 <- fData(nd)[ps, "gene.symbol"]          
    x <- as.numeric(nd$group) + runif(length(o), -0.25, 0.25)
    names(x) <- names(o)
    meds <- as.numeric(d[ps, paste0(levels(nd$group), ".log2mean")])
    names(meds) <- levels(nd$group)
    ees <- as.numeric(d[ps, paste0(levels(nd$group), ".log2se")])
    names(ees) <- levels(nd$group)
    png(paste0(.res, '/Stripcharts_zoom/', g, ".png"), width=900, height=600)
    par(mar=c(5, 10, 5, 3))
    plot(x, o, pch=16, col=makeTransparent(cols[as.numeric(nd$group)], 0.7), xlim=c(0.5, nlevels(nd$group)+1.5),
         ylab='Intensity\n\n', xlab='', sub='', main=paste0(g, " / ", g2, " - zoom"),
         cex.lab=1.3, cex.axis=1.3, cex.sub=1.3, cex.main=1.5, axes=F, cex=2.5)
    points(1:length(meds), meds, pch=5, col=cols, cex=4)
    for(n in 1:length(meds))
        arrows(x0=n, x1=n, y0=meds[n]-ees[n], y1=meds[n]+ees[n], angle=90, length=0.2, code=3,
               col=cols[n], lwd=4)
    axis(1, at=c(-100, 1:nlevels(nd$group), nlevels(nd$group)+10),
         lab=c("", paste("\n", gsub("-", "\n", levels(nd$group))), ""))
    axis(2, at=seq(1, 1000, 1), labels=2^seq(1, 1000, 1), las=2)
    legend(x='bottomright', inset=0.01, pch=16, col=cols, legend=levels(nd$group), cex=1.5, pt.cex=2)
    dev.off()
}, mc.cores=12)


ylim <- quantile(as.numeric(exprs(nd)), c(0.001, 0.999), na.rm=T)
ylim <- c(14, 26)
dir.create(paste0(.res, 'Stripcharts'), F)
pp <- mclapply(pss, function(ps)
{
    o <- exprs(nd)[ps, ]
    g <- ps
    g2 <- fData(nd)[ps, "gene.symbol"]          
    x <- as.numeric(nd$group) + runif(length(o), -0.25, 0.25)
    names(x) <- names(o)
    meds <- as.numeric(d[ps, paste0(levels(nd$group), ".log2mean")])
    names(meds) <- levels(nd$group)
    ees <- as.numeric(d[ps, paste0(levels(nd$group), ".log2se")])
    names(ees) <- levels(nd$group)
    png(paste0(.res, '/Stripcharts/', g, ".png"), width=900, height=600)
    par(mar=c(5, 10, 5, 3))
    plot(x, o, pch=16, col=makeTransparent(cols[as.numeric(nd$group)], 0.7), xlim=c(0.5, nlevels(nd$group)+1.5),
         ylab='Intensity\n\n', xlab='', sub='', main=paste0(g, " / ", g2, ""),
         cex.lab=1.3, cex.axis=1.3, cex.sub=1.3, cex.main=1.5, axes=F, cex=2.5, ylim=ylim)
    points(1:length(meds), meds, pch=5, col=cols, cex=4)
    for(n in 1:length(meds))
        arrows(x0=n, x1=n, y0=meds[n]-ees[n], y1=meds[n]+ees[n], angle=90, length=0.2, code=3,
               col=cols[n], lwd=4)
    axis(1, at=c(-100, 1:nlevels(nd$group), nlevels(nd$group)+10),
         lab=c("", paste("\n", gsub("-", "\n", levels(nd$group))), ""))
    axis(2, at=seq(1, 1000, 1), labels=2^seq(1, 1000, 1), las=2)
    legend(x='bottomright', inset=0.01, pch=16, col=cols, legend=levels(nd$group), cex=1.5, pt.cex=2)
    dev.off()
}, mc.cores=12)

################################################################################################################
  # QQ-plot pvals

rm(list=ls())
d <- readRDS(paste0(.res, "dcomp.rds"))
source("code/functions/utils.R")

npv <- colnames(d)[regexpr("\\.pv", colnames(d))>0]
npv <- unique(gsub("adj\\.", "", npv))
comps <- gsub("\\.pv", "", npv)

dir.create(paste(.res, "plot_pvals", sep=''), F)
for (j in 1:length(npv))
{
    png(paste0(.res, "plot_pvals/", comps[j], ".png"), width=1500, height=800);
    par(mfrow=c(1, 2));
    qqpval(d[, npv[j]], main=paste0('Class comparison p-values\n', comps[j]))
    hist(d[, npv[j]], main=paste0('Class comparison p-values\n', comps[j]))
    dev.off();
}

################################################################################################################
  # Create results tables

rm(list=ls())
lc <- readRDS(paste0(.res, "lcomp.rds"))
nd <- readRDS(paste0(.res, "nd.log.prot.corr.rds"))
source("code/functions/write.html.mod.R")


lr <- lc
pp <- mclapply(1:length(lr), function(j)
{
    drr <- lr[[j]][1:200, ]
    links <- vector('list',length=ncol(drr))
    names(links) <- colnames(drr)
    plots <- links;
    links$gene.symbol <- paste0("http://www.ncbi.nlm.nih.gov/gene/?term=",
                                drr$gene.symbol)
    links$protein.group <- paste0("https://www.uniprot.org/uniprot/",
                                  drr$protein.group)
    links$plot <- plots$plot <- paste0('Stripcharts/', rownames(drr), ".png")
    links$"plot zoom" <- plots$"plot zoom" <- paste0('Stripcharts_zoom/', rownames(drr), ".png")
    tit <- paste("Results for ",  names(lr)[j])
    write.html.mod(drr, file=paste(.res, "/Diff.expr_", names(lr)[j], ".html",  sep=''),
                   links=links,tiny.pic=plots,
                   title=tit,tiny.pic.size = 200)
}, mc.cores=12)

################################################################################################################
# Volcano plots for results

rm(list=ls())
lc <- readRDS(paste0(.res, "lcomp.rds"))
nd <- readRDS(paste0(.res, "nd.log.filt.rds"))
source("code/functions/write.html.mod.R")

# Common settings
thfc     <- 1.20
thpv_nom <- 0.05
thfdr_B  <- 0.10

col_fdr  <- "red3"       
col_pv   <- "tomato2"    
col_ns   <- "gray70"     

pt_size  <- 2            
xlims    <- c(-3, 3)     
ylims    <- c(0, 4.5)    


# Helper function

prep_data <- function(lc, comp, thfc, sig_col, sig_thr, sig_label) {
  d <- lc[[comp]]
  d <- d[, c("protein.group", "gene.symbol", paste0(comp, c(".fc", ".pv", ".adj.pv")))]
  colnames(d) <- gsub(paste0(comp, "\\."), "", colnames(d))
  d$log2fc <- log2(abs(d$fc)) * sign(d$fc)
  d <- d[!is.na(d$fc), ]
  d$logpv <- -log10(d$pv)
  d$group <- "Not significant"
  d$group[d$fc >  thfc & d[[sig_col]] < sig_thr] <- sig_label
  d$group[d$fc < -thfc & d[[sig_col]] < sig_thr] <- sig_label
  d$label <- ifelse(d$group != "Not significant", d$gene.symbol, "")
  d$label[d$label == "NA"] <- d$protein.group[d$label == "NA"]
  d
}

base_theme <- function() {
  theme_classic() +
    theme(legend.position   = "none",
          panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.8),
          axis.title        = element_text(size = 16),
          axis.text         = element_text(size = 12),
          plot.title        = element_text(size = 14, face = "bold"),
          axis.ticks.length = unit(0.3, "cm"))
}


# Mild vs Asymptomatic  (nominal p < 0.05)

dA <- prep_data(lc, "Mild-Asymptomatic", thfc,
                sig_col = "pv", sig_thr = thpv_nom,
                sig_label = "|FC| > 1.20, p < 0.05")

pA <- ggplot(dA) +
  geom_point(aes(log2fc, logpv, colour = group), size = pt_size) +
  scale_colour_manual(name = "",
                      values = c("|FC| > 1.20, p < 0.05" = col_pv,
                                 "Not significant"        = col_ns),
                      drop = FALSE) +
  geom_text_repel(aes(log2fc, logpv, label = label),
                  size = 2.7, col = "black", segment.color = "grey80", max.overlaps = 20) +
  geom_vline(xintercept = log2(c(1/thfc, thfc)), color = "gray40", linewidth = 0.8, linetype = 3) +
  geom_hline(yintercept = -log10(thpv_nom), color = "gray40", linewidth = 0.8, linetype = 3) +
  coord_cartesian(xlim = xlims, ylim = ylims) +
  labs(x = "\nlog2(FC)", y = "-log10(p-value)\n", title = "Mild vs. Asymptomatic") +
  base_theme()


# Severe vs Mild  (FDR < 10%)

dB <- prep_data(lc, "Severe-Mild", thfc,
                sig_col = "adj.pv", sig_thr = thfdr_B,
                sig_label = "|FC| > 1.20, FDR < 10%")

thpv_B <- max(dB$pv[dB$adj.pv < thfdr_B])   

pB <- ggplot(dB) +
  geom_point(aes(log2fc, logpv, colour = group), size = pt_size) +
  scale_colour_manual(name = "",
                      values = c("|FC| > 1.20, FDR < 10%" = col_fdr,
                                 "Not significant"         = col_ns),
                      drop = FALSE) +
  geom_text_repel(aes(log2fc, logpv, label = label),
                  size = 2.7, col = "black", segment.color = "grey80", max.overlaps = 15) +
  geom_vline(xintercept = log2(c(1/thfc, thfc)), color = "gray40", linewidth = 0.8, linetype = 3) +
  geom_hline(yintercept = -log10(thpv_nom), color = "gray40", linewidth = 0.8, linetype = 3) +
  geom_hline(yintercept = -log10(thpv_B),   color = "gray40", linewidth = 0.8, linetype = 2) +
  coord_cartesian(xlim = xlims, ylim = ylims) +
  labs(x = "\nlog2(FC)", y = "-log10(p-value)\n", title = "Severe vs. Mild") +
  base_theme()


# Critical vs Severe  (nominal p < 0.05)

dC <- prep_data(lc, "Critical-Severe", thfc,
                sig_col = "pv", sig_thr = thpv_nom,
                sig_label = "|FC| > 1.20, p < 0.05")

pC <- ggplot(dC) +
  geom_point(aes(log2fc, logpv, colour = group), size = pt_size) +
  scale_colour_manual(name = "",
                      values = c("|FC| > 1.20, p < 0.05" = col_pv,
                                 "Not significant"        = col_ns),
                      drop = FALSE) +
  geom_text_repel(aes(log2fc, logpv, label = label),
                  size = 2.7, col = "black", segment.color = "grey80", max.overlaps = 20) +
  geom_vline(xintercept = log2(c(1/thfc, thfc)), color = "gray40", linewidth = 0.8, linetype = 3) +
  geom_hline(yintercept = -log10(thpv_nom), color = "gray40", linewidth = 0.8, linetype = 3) +
  coord_cartesian(xlim = xlims, ylim = ylims) +
  labs(x = "\nlog2(FC)", y = "-log10(p-value)\n", title = "Critical vs. Severe") +
  base_theme()


# Shared legend
d_legend <- data.frame(
  x     = c(0, 0, 0),
  y     = c(0, 0, 0),
  group = c("|FC| > 1.20, FDR < 10%", "|FC| > 1.20, p < 0.05", "Not significant")
)
d_legend$group <- factor(d_legend$group,
                         levels = c("|FC| > 1.20, FDR < 10%", "|FC| > 1.20, p < 0.05", "Not significant"))

p_legend <- ggplot(d_legend) +
  geom_point(aes(x, y, colour = group), size = 4) +
  scale_colour_manual(name = "",
                      values = c("|FC| > 1.20, FDR < 10%" = col_fdr,
                                 "|FC| > 1.20, p < 0.05"  = col_pv,
                                 "Not significant"         = col_ns)) +
  theme_void() +
  theme(legend.position  = "bottom",
        legend.text      = element_text(size = 12),
        legend.key.size  = unit(0.5, "cm")) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 4)))

# Extract legend as a grob
legend_grob <- cowplot::get_legend(p_legend)


# Combined layout
plots_row <- (pA | pB | pC) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

combined <- cowplot::plot_grid(
  plots_row,
  legend_grob,
  ncol    = 1,
  rel_heights = c(10, 0.8)
)

ggsave(paste0(.res, "volcano_combined.pdf"), combined, width = 21, height = 8)
ggsave(paste0(.res, "volcano_combined.png"), combined, width = 21, height = 8, dpi = 300)



################################################################################################################
  # Create results tables for possible markers of differences in extraction

rm(list=ls())
lc <- readRDS(paste0(.res, "lcomp.rds"))
nd <- readRDS(paste0(.res, "nd.log.prot.corr.rds"))
source("code/functions/write.html.mod.R")


d <- fData(nd)
prots <- d[regexpr("epith", tolower(d$protein.name)) > 0 | regexpr("skin", tolower(d$protein.name)) > 0, 
           "protein.group"]
prots <- c(prots, "P68871", "P69891", "P02042", "P69905", "P02768")


### and remove due to not properly quantified?

lr <- lc
pp <- mclapply(1:length(lr), function(j)
{
    drr <- lr[[j]][prots, ]
    links <- vector('list',length=ncol(drr))
    names(links) <- colnames(drr)
    plots <- links;
    links$gene.symbol <- paste0("http://www.ncbi.nlm.nih.gov/gene/?term=",
                                drr$gene.symbol)
    links$protein.group <- paste0("https://www.uniprot.org/uniprot/",
                                  drr$protein.group)
    links$plot <- plots$plot <- paste0('Stripcharts/', rownames(drr), ".png")
    links$"plot zoom" <- plots$"plot zoom" <- paste0('Stripcharts_zoom/', rownames(drr), ".png")
    tit <- paste("Results for ",  names(lr)[j])
    write.html.mod(drr, file=paste(.res, "/Diff.expr_", names(lr)[j], "_extract_markers.html",  sep=''),
                   links=links,tiny.pic=plots,
                   title=tit,tiny.pic.size = 200)
}, mc.cores=12)

################################################################################################################
################################################################################################################
################################################################################################################


# ################################################################################################################
# # Volcano plot - Severe vs Mild
# 
# rm(list=ls())
# lc <- readRDS(paste0(.res, "lcomp.rds"))
# nd <- readRDS(paste0(.res, "nd.log.filt.rds"))
# source("code/functions/write.html.mod.R")
# 
# 
# comp <- "Severe-Mild"
# d <- lc[[comp]]
# d <- d[, c("protein.group", "gene.symbol", paste0(comp, c(".fc", ".pv", ".adj.pv")))]
# colnames(d) <- gsub(paste0(comp, "\\."), "", colnames(d))
# d$log2fc <- log2(abs(d$fc))*sign(d$fc)
# d <- d[!is.na(d$fc), ]
# 
# d <- d[order(abs(d$fc), decreasing=T), ]
# d$num <- 1:nrow(d)
# 
# 
# thfdr <- 0.10
# thpv <- max(d$pv[d$adj.pv < thfdr])
# thpv2 <- max(d$pv[d$adj.pv < 0.05])
# thfc <- 1.25
# d$logpv <- -log10(d$pv)
# 
# #cols <- brewer.pal(n = nlevels(d$group)-1 +4, name="Paired")[-c(3, 4, 5, 6)]
# #cols <- plasma(nlevels(d$group)-1)
# #cols[length(cols)] <- "yellow3"
# cols <- c("red3", "gray70")
# ps <- c(2, 5)
# d$num <- as.character(d$num)
# d$sign <- "No.sign"
# s <- d$fc > thfc & d$adj.pv < thfdr
# d$sign[s] <- "up"
# s <- d$fc < -thfc & d$adj.pv < thfdr
# d$sign[s] <- "down"
# d$sign <- factor(d$sign, levels=rev(c("No.sign", "down", "up")))
# mx <- ceiling(max(d$logpv))
# d$group <- factor(as.numeric(d$sign=="No.sign"), levels=0:1,
#                   labels=c(paste0("FC > 1.25, FDR < ", thfdr*100, "%"), "Not sign."))
# d$gene.symbol[d$group=='Not sign.'] <- ""
# s <- d$gene.symbol=="NA"
# d$gene.symbol[s] <- d$protein.group[s]
# d$num[d$sign=="No.sign"] <- ""
# d$protein.group[d$sign=="No.sign"] <- ""
# 
# log2fc <- d$log2fc
# logpv <- d$logpv
# adj.pv <- d$adj.pv
# gene.symbol <- d$gene.symbol
# plot <- ggplot(d) +
#   geom_point(aes(log2fc, logpv, colour=group), size=ps[as.numeric(d$adj.pv < thpv2)+1]*0.8) +
#   #    xlim(c(-2.8, 2.8)) +
#   #    ylim(c(0, mx)) +
#   scale_colour_manual(name="", values=cols) + 
#   geom_text_repel(aes(log2fc, logpv, label = gene.symbol),
#                   size=2.7, col="black", segment.color = "grey80", max.overlaps=15) +
#   geom_vline(xintercept=log2(c(1/thfc, thfc)), color='gray40', size=0.8, linetype=3) +
#   geom_hline(yintercept= -log10(thpv), color='gray40', size=0.8, linetype=3) +
#   ggtitle(paste0("")) +
#   theme(plot.title=element_text(size=18)) +
#   labs(x='\nlog2(FC)', y='-log10(p-value)\n') +
#   theme(plot.margin = unit(c(1, 1, 3, 1), "cm")) +
#   theme_classic() +
#   guides(color = guide_legend(override.aes = list(size=5))) + 
#   theme(legend.text=element_text(size=13)) + 
#   theme(panel.border = element_rect(colour = "black", fill=NA, size=0.8)) +
#   theme(axis.title=element_text(size=20)) +
#   theme(axis.ticks.length = unit(0.3, "cm")) +
#   theme(axis.text = element_text(size=15))
# pdf(paste0(.res, "volcano_", comp, ".pdf"), width=13, height=8)
# print(plot)
# dev.off()
# graphics.off()
# 
# ################################################################################################################
# # Volcano plot - Mild vs Asymptomatic
# 
# rm(list=ls())
# lc <- readRDS(paste0(.res, "lcomp.rds"))
# nd <- readRDS(paste0(.res, "nd.log.filt.rds"))
# source("code/functions/write.html.mod.R")
# 
# 
# comp <- "Mild-Asymptomatic"
# d <- lc[[comp]]
# d <- d[, c("protein.group", "gene.symbol", paste0(comp, c(".fc", ".pv", ".adj.pv")))]
# colnames(d) <- gsub(paste0(comp, "\\."), "", colnames(d))
# d$log2fc <- log2(abs(d$fc))*sign(d$fc)
# d <- d[!is.na(d$fc), ]
# 
# d <- d[order(abs(d$fc), decreasing=T), ]
# d$num <- 1:nrow(d)
# 
# 
# thpv <- 0.05
# thfdr <- max(d$adj.pv[d$pv < thpv])
# thfdr <- 0.25
# thpv <- max(d$pv[d$adj.pv < thfdr])
# #thpv2 <- max(d$pv[d$adj.pv < 0.05])
# thfc <- 1.20
# d$logpv <- -log10(d$pv)
# 
# #cols <- brewer.pal(n = nlevels(d$group)-1 +4, name="Paired")[-c(3, 4, 5, 6)]
# #cols <- plasma(nlevels(d$group)-1)
# #cols[length(cols)] <- "yellow3"
# cols <- c("red3", "gray70")
# ps <- c(2, 5)
# d$num <- as.character(d$num)
# d$sign <- "No.sign"
# s <- d$fc > thfc & d$adj.pv < thfdr
# d$sign[s] <- "up"
# s <- d$fc < -thfc & d$adj.pv < thfdr
# d$sign[s] <- "down"
# d$sign <- factor(d$sign, levels=rev(c("No.sign", "down", "up")))
# mx <- ceiling(max(d$logpv))
# d$group <- factor(as.numeric(d$sign=="No.sign"), levels=0:1,
#                   labels=c(paste0("FC > 1.25, FDR < ", thfdr*100, "%"), "Not sign."))
# d$gene.symbol[d$group=='Not sign.'] <- ""
# s <- d$gene.symbol=="NA"
# d$gene.symbol[s] <- d$protein.group[s]
# d$num[d$sign=="No.sign"] <- ""
# d$protein.group[d$sign=="No.sign"] <- ""
# 
# log2fc <- d$log2fc
# logpv <- d$logpv
# adj.pv <- d$adj.pv
# gene.symbol <- d$gene.symbol
# plot <- ggplot(d) +
#   geom_point(aes(log2fc, logpv, colour=group), size=ps[as.numeric(d$adj.pv < thpv)+1]*0.8) +
#   #    xlim(c(-2.8, 2.8)) +
#   #    ylim(c(0, mx)) +
#   scale_colour_manual(name="", values=cols) + 
#   geom_text_repel(aes(log2fc, logpv, label = gene.symbol),
#                   size=2.7, col="black", segment.color = "grey80", max.overlaps=20) +
#   geom_vline(xintercept=log2(c(1/thfc, thfc)), color='gray40', size=0.8, linetype=3) +
#   geom_hline(yintercept= -log10(thpv), color='gray40', size=0.8, linetype=3) +
#   ggtitle(paste0("")) +
#   theme(plot.title=element_text(size=18)) +
#   labs(x='\nlog2(FC)', y='-log10(p-value)\n') +
#   theme(plot.margin = unit(c(1, 1, 3, 1), "cm")) +
#   theme_classic() +
#   guides(color = guide_legend(override.aes = list(size=5))) + 
#   theme(legend.text=element_text(size=13)) + 
#   theme(panel.border = element_rect(colour = "black", fill=NA, size=0.8)) +
#   theme(axis.title=element_text(size=20)) +
#   theme(axis.ticks.length = unit(0.3, "cm")) +
#   theme(axis.text = element_text(size=15))
# pdf(paste0(.res, "volcano_", comp, ".pdf"), width=13, height=8)
# print(plot)
# dev.off()
# graphics.off()
# 
# ################################################################################################################
# # Volcano plot - Critic vs Severe
# 
# rm(list=ls())
# lc <- readRDS(paste0(.res, "lcomp.rds"))
# nd <- readRDS(paste0(.res, "nd.log.filt.rds"))
# source("code/functions/write.html.mod.R")
# 
# 
# comp <- "Critical-Severe"
# d <- lc[[comp]]
# d <- d[, c("protein.group", "gene.symbol", paste0(comp, c(".fc", ".pv", ".adj.pv")))]
# colnames(d) <- gsub(paste0(comp, "\\."), "", colnames(d))
# d$log2fc <- log2(abs(d$fc))*sign(d$fc)
# d <- d[!is.na(d$fc), ]
# 
# d <- d[order(abs(d$fc), decreasing=T), ]
# d$num <- 1:nrow(d)
# 
# 
# thpv <- 0.05
# #thfdr <- max(d$adj.pv[d$pv < thpv])
# #thfdr <- 0.25
# #thpv <- max(d$pv[d$pv < thpv])
# #thpv2 <- max(d$pv[d$adj.pv < 0.05])
# thfc <- 1.20
# d$logpv <- -log10(d$pv)
# 
# #cols <- brewer.pal(n = nlevels(d$group)-1 +4, name="Paired")[-c(3, 4, 5, 6)]
# #cols <- plasma(nlevels(d$group)-1)
# #cols[length(cols)] <- "yellow3"
# cols <- c("red3", "gray70")
# ps <- c(2, 5)
# d$num <- as.character(d$num)
# d$sign <- "No.sign"
# s <- d$fc > thfc & d$pv < thpv
# d$sign[s] <- "up"
# s <- d$fc < -thfc & d$pv < thpv
# d$sign[s] <- "down"
# d$sign <- factor(d$sign, levels=rev(c("No.sign", "down", "up")))
# mx <- ceiling(max(d$logpv))
# d$group <- factor(as.numeric(d$sign=="No.sign"), levels=0:1,
#                   labels=c(paste0("FC > 1.25, pv < ", thpv), "Not sign."))
# d$gene.symbol[d$group=='Not sign.'] <- ""
# s <- d$gene.symbol=="NA"
# d$gene.symbol[s] <- d$protein.group[s]
# d$num[d$sign=="No.sign"] <- ""
# d$protein.group[d$sign=="No.sign"] <- ""
# 
# log2fc <- d$log2fc
# logpv <- d$logpv
# adj.pv <- d$adj.pv
# gene.symbol <- d$gene.symbol
# plot <- ggplot(d) +
#   geom_point(aes(log2fc, logpv, colour=group), size=ps[as.numeric(d$adj.pv < thpv)+1]*0.8) +
#   #    xlim(c(-2.8, 2.8)) +
#   #    ylim(c(0, mx)) +
#   scale_colour_manual(name="", values=cols) + 
#   geom_text_repel(aes(log2fc, logpv, label = gene.symbol),
#                   size=2.7, col="black", segment.color = "grey80", max.overlaps=20) +
#   geom_vline(xintercept=log2(c(1/thfc, thfc)), color='gray40', size=0.8, linetype=3) +
#   geom_hline(yintercept= -log10(thpv), color='gray40', size=0.8, linetype=3) +
#   ggtitle(paste0("")) +
#   theme(plot.title=element_text(size=18)) +
#   labs(x='\nlog2(FC)', y='-log10(p-value)\n') +
#   theme(plot.margin = unit(c(1, 1, 3, 1), "cm")) +
#   theme_classic() +
#   guides(color = guide_legend(override.aes = list(size=5))) + 
#   theme(legend.text=element_text(size=13)) + 
#   theme(panel.border = element_rect(colour = "black", fill=NA, size=0.8)) +
#   theme(axis.title=element_text(size=20)) +
#   theme(axis.ticks.length = unit(0.3, "cm")) +
#   theme(axis.text = element_text(size=15))
# pdf(paste0(.res, "volcano_", comp, ".pdf"), width=13, height=8)
# print(plot)
# dev.off()
# graphics.off()
# 
# ################################################################################################################
# # Volcano plot - Hospitalized vs Not.hospitalized
# 
# rm(list=ls())
# lc <- readRDS(paste0(.res, "lcomp.rds"))
# nd <- readRDS(paste0(.res, "nd.log.filt.rds"))
# source("code/functions/write.html.mod.R")
# 
# 
# comp <- "Hospitalized-Not.hospitalized"
# d <- lc[[comp]]
# d <- d[, c("protein.group", "gene.symbol", paste0(comp, c(".fc", ".pv", ".adj.pv")))]
# colnames(d) <- gsub(paste0(comp, "\\."), "", colnames(d))
# d$log2fc <- log2(abs(d$fc))*sign(d$fc)
# d <- d[!is.na(d$fc), ]
# 
# d <- d[order(abs(d$fc), decreasing=T), ]
# d$num <- 1:nrow(d)
# 
# 
# thfdr <- 0.10
# thpv <- max(d$pv[d$adj.pv < thfdr])
# thpv2 <- max(d$pv[d$adj.pv < 0.05])
# thfc <- 1.25
# d$logpv <- -log10(d$pv)
# 
# #cols <- brewer.pal(n = nlevels(d$group)-1 +4, name="Paired")[-c(3, 4, 5, 6)]
# #cols <- plasma(nlevels(d$group)-1)
# #cols[length(cols)] <- "yellow3"
# cols <- c("red3", "gray70")
# ps <- c(2, 5)
# d$num <- as.character(d$num)
# d$sign <- "No.sign"
# s <- d$fc > thfc & d$adj.pv < thfdr
# d$sign[s] <- "up"
# s <- d$fc < -thfc & d$adj.pv < thfdr
# d$sign[s] <- "down"
# d$sign <- factor(d$sign, levels=rev(c("No.sign", "down", "up")))
# mx <- ceiling(max(d$logpv))
# d$group <- factor(as.numeric(d$sign=="No.sign"), levels=0:1,
#                   labels=c(paste0("FC > 1.25, FDR < ", thfdr*100, "%"), "Not sign."))
# d$gene.symbol[d$group=='Not sign.'] <- ""
# s <- d$gene.symbol=="NA"
# d$gene.symbol[s] <- d$protein.group[s]
# d$num[d$sign=="No.sign"] <- ""
# d$protein.group[d$sign=="No.sign"] <- ""
# 
# log2fc <- d$log2fc
# logpv <- d$logpv
# adj.pv <- d$adj.pv
# gene.symbol <- d$gene.symbol
# plot <- ggplot(d) +
#   geom_point(aes(log2fc, logpv, colour=group), size=ps[as.numeric(d$adj.pv < thpv2)+1]*0.5) +
#   #    xlim(c(-2.8, 2.8)) +
#   #    ylim(c(0, mx)) +
#   scale_colour_manual(name="", values=cols) + 
#   geom_text_repel(aes(log2fc, logpv, label = gene.symbol),
#                   size=2.7, col="black", segment.color = "grey80", max.overlaps=30) +
#   geom_vline(xintercept=log2(c(1/thfc, thfc)), color='gray40', size=0.8, linetype=3) +
#   geom_hline(yintercept= -log10(thpv), color='gray40', size=0.8, linetype=3) +
#   ggtitle(paste0("")) +
#   theme(plot.title=element_text(size=18)) +
#   labs(x='\nlog2(FC)', y='-log10(p-value)\n') +
#   theme(plot.margin = unit(c(1, 1, 3, 1), "cm")) +
#   theme_classic() +
#   guides(color = guide_legend(override.aes = list(size=5))) + 
#   theme(legend.text=element_text(size=13)) + 
#   theme(panel.border = element_rect(colour = "black", fill=NA, size=0.8)) +
#   theme(axis.title=element_text(size=20)) +
#   theme(axis.ticks.length = unit(0.3, "cm")) +
#   theme(axis.text = element_text(size=15))
# pdf(paste0(.res, "volcano_", comp, ".pdf"), width=13, height=8)
# print(plot)
# dev.off()
# graphics.off()
