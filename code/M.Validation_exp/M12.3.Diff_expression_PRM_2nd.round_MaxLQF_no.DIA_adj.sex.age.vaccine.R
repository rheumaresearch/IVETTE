#################################################################################################################
#################################################################################################################
### M12.3.Diff_expression_PRM_2nd.round_MaxLQF_no.DIA_adj.sex.age.vaccine.R
### Project: jcalvet_202110_marato
#################################################################################################################
#################################################################################################################
  # Libraries, paths and options

rm(list=ls())
.res <- paste0(.res0, "M.Validation_exp/M12.3.Diff_expression_PRM_2nd.round_MaxLQF_no.DIA_adj.sex.age.vaccine/")
dir.create(.res, F, T)

#################################################################################################################
  # Select data and apply filters

rm(list=ls())
nd <- readRDS(paste0(.dat, "2_PRM/nd.prot.no0.norm.area.maxlfq.round2.rds"))

nds <- nd[, nd$group!='Control']
nds$group <- droplevels(nds$group)

nds$total.area.log2 <- log2(nds$total.area)

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

form <- "y ~ -1 + group + center + sex + age + vaccinated + quant.batch"
lr <- mclapply(1:nrow(nd), function(i)
{

    d <- cbind(pData(nd), y=exprs(nd)[i, ])

    ### Best model -> only for labelling
    
    aicg <- try(AIC(lm(y ~ group + center + sex + vaccinated + age + quant.batch, data=d)), silent=T)
    aic2 <- try(AIC(lm(y ~ group2 + center + sex + vaccinated + age +  quant.batch, data=d)), silent=T)
    aicn <- try(AIC(lm(y ~ groupn + center + sex + vaccinated + age + quant.batch, data=d)), silent=T)

    ### Categorical models for all

    m <- try(lm(as.formula(form), data=d), silent=T)
    
    if (class(m)!='try-error' && all(!is.na(summary(m)$coef[, 2])))
    {

        bs <- coef(m)
        nbs <- names(coef(m))
        nva <- c("center", "sex", "age", "vaccinated", "obesity", "dm", "dlp", "aht", "quant.batch", "pch1", "col.hemol")
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
            fmm <- gsub("center", paste0(colnames(dmm)[regexpr("^center", colnames(dmm)) > 0], collapse=' + '), fmm)            
            fmm <- gsub("sex", paste0(colnames(dmm)[regexpr("^sex", colnames(dmm)) > 0], collapse=' + '), fmm)
            fmm <- gsub("vaccinated", paste0(colnames(dmm)[regexpr("^vaccinated", colnames(dmm)) > 0], collapse=' + '), fmm)
            fmm <- gsub("obesity", paste0(colnames(dmm)[regexpr("^obesity", colnames(dmm)) > 0], collapse=' + '), fmm)
            fmm <- gsub("dm", paste0(colnames(dmm)[regexpr("^dm", colnames(dmm)) > 0], collapse=' + '), fmm)
            fmm <- gsub("dlp", paste0(colnames(dmm)[regexpr("^dlp", colnames(dmm)) > 0], collapse=' + '), fmm)
            fmm <- gsub("aht", paste0(colnames(dmm)[regexpr("^aht", colnames(dmm)) > 0], collapse=' + '), fmm)
            if (!"quant.batch"%in%nvas)
                fmm <- gsub(" \\+ quant.batch", "", fmm)
            if (!"center"%in%nvas)
                fmm <- gsub(" \\+ center", "", fmm)
            if (!"sex"%in%nvas)
                fmm <- gsub(" \\+ sex", "", fmm)
            if (!"vaccinated"%in%nvas)
                fmm <- gsub(" \\+ vaccinated", "", fmm)
            if (!"obesity"%in%nvas)
                fmm <- gsub(" \\+ obesity", "", fmm)
            if (!"dm"%in%nvas)
                fmm <- gsub(" \\+ dm", "", fmm)
            if (!"dlp"%in%nvas)
                fmm <- gsub(" \\+ dlp", "", fmm)
            if (!"aht"%in%nvas)
                fmm <- gsub(" \\+ aht", "", fmm)
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

length(unique(fData(nd)$protein.id))
df <- unique(fData(nd)[, c("protein.id", "protein.name", "gene.symbol")])
rownames(df) <- df$protein.id

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
nd <- readRDS(paste0(.res, "nd.log.filt.rds"))


### All results

dr <- fData(nd)
rownames(dr) <- dr$protein.id
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
da <- readRDS(paste0(.dat, "2_PRM/dproteins.sel.prm.rds"))

colnames(da) <- paste0("disc.", colnames(da))

ncomp <- colnames(dc)[regexpr("\\.adj\\.pv$", colnames(dc)) > 0]
ncomp <- gsub("\\.adj\\.\\.pv$", "",ncomp)

vsel <- c("protein.id", "protein.name", "gene.symbol",         
          "global.mean", "global.se")                          

dc$plot <- ""
dc$"plot zoom" <- ""

conts2 <- conts
names(conts2)[names(conts2)=='Hospitalized-Not.hospitalized'] <- "Severe-Mild"

lr <- sapply(1:length(conts2), function(i)
      {
          o <- conts2[i]
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
          ds <- dc[, c(vsel, nv2, "model.aic", "plot", "plot zoom")]
          ds <- cbind(ds, da[rownames(ds), c(colnames(da)[c(6:10)], paste0("disc.", names(o), ".fc"))])
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
ylim <- c(9, 28)
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
nd <- readRDS(paste0(.res, "nd.log.filt.rds"))
source("code/functions/write.html.mod.R")


lr <- lc
pp <- mclapply(1:length(lr), function(j)
{
    drr <- lr[[j]]
    links <- vector('list',length=ncol(drr))
    names(links) <- colnames(drr)
    plots <- links;
    links$gene.symbol <- paste0("http://www.ncbi.nlm.nih.gov/gene/?term=",
                                drr$gene.symbol)
    links$protein.id <- paste0("https://www.uniprot.org/uniprot/",
                                  drr$protein.id)
    links$plot <- plots$plot <- paste0('Stripcharts/', rownames(drr), ".png")
    links$"plot zoom" <- plots$"plot zoom" <- paste0('Stripcharts_zoom/', rownames(drr), ".png")
    tit <- paste("Results for ",  names(lr)[j])
    write.html.mod(drr, file=paste(.res, "/Diff.expr_", names(lr)[j], ".html",  sep=''),
                   links=links,tiny.pic=plots,
                   title=tit,tiny.pic.size = 200)
}, mc.cores=12)


################################################################################################################
################################################################################################################
################################################################################################################

