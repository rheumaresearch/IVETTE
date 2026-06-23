#################################################################################################################
#################################################################################################################
### A.Read_data/A4.Descriptives_all.R
### Project: jcalvet_202110_marato
#################################################################################################################
#################################################################################################################
  # Libraries, paths and options

rm(list=ls())
.res <- paste0(.res0, "A.Read_data/A4.Descriptives_all/")
dir.create(.res, F, T)

#################################################################################################################
  # Descriptives - All

rm(list=ls())
d <- readRDS(paste0(.dat, "0_metadata/data.all.rds"))
source("code/functions/clin_descrip.R")

ny <- "group"
nxs <- c("disc", "prm", "elisa",
         "center", "sex", "age", 
         "previous.covid", "vaccinated", 
         "aht", "dm", "obesity", "dlp",
         "ferritin", "pcr", "dd", "ldh", "lymphocytes", "neutrophils", "leukocytes")
#         "cortis", "remdesivir", "toci",         
#         "bilateral.pneumony", "nimv.oaf", "icu", "exitus")
nv <- nxs
names(nv) <- c("Discovery - DIA", "PRM validation", "ELISA validation",
               "Center", "Sex", "Age at diagnosi", 
               "Previous SARS-CoV-2 infection", "Vaccinated against SARS-CoV-2",
               "Arterial Hypertension", "Diabetes Mellitus", "Obesity", "Dyslipidemia",
               "Ferritin", "C-Reactive Protein (CRP)", "D-Dimer (DD)", "Lactate Dehydrogenase (LDH)",               
               "Lymphocytes", "Neutrophils", "Leukocytes")
#               "Glucocorticoids", "Remdesivir", "Tocilizumab",
#               "Bilateral Pneumony", "Non-Invasive Mechanical Ventilation (NIMV) / High Flow Oxygen Therapy (HFO)",
#               "Intensive Care Unit (ICU)", "Exitus")


for (o in nv[!nv%in%c("center", "series")])
{
    if (is.factor(d[, o])) d[, o] <- relevel(d[, o], ref=levels(d[, o])[2])
}



### P-values

set.seed(513335)
pvs <- sapply(nv, function(o)
{
    s <- rep(TRUE, nrow(d))
    y <- d[s, o]
    x <- droplevels(d[s, "group"])
    s <- x != "Control"
    y <- y[s]
    x <- droplevels(x[s])
    if (is.numeric(y))
    {
        r <- kruskal.test(y ~ x)$p.value
    }else if (is.factor(y))
    {
        r <- try(fisher.test(table(y, x), simulate.p.value=TRUE, B=10000)$p.value, silent=T)        
        if (class(r)=='try-error') r <- NA
    }else if (o=="treatment")
    {
        r <- NA
    }
    r
})


### All samples

dir.create(paste(.res, "all_all/", sep=''), F)
r <- clin.descrip(nv, d, ndir=paste0(.res, "all_all/"), bc=bcs)
ra <- r


### Levels

lr <- list()
for (j in 1:nlevels(d$group))
{
    s <- !is.na(d$group) & d$group==levels(d$group)[j]
    dir.create(paste(.res, levels(d$group)[j], "_all/", sep=''), F)
    r <- clin.descrip(nv, d[s, ], ndir=paste0(.res, levels(d$group)[j], "_all/"), bc=bcs)
    lr[[levels(d$group)[j]]] <- r
}


### Merge

ng <- summary(d$group)
ng <- ng[names(ng)!="NA's"]
ngg <- paste0(ng, "<br>(", formatC(ng/sum(ng)*100, format='f', dig=1), "%)")
r <- cbind(ra[, -4], All=ra[, 4], do.call(cbind, sapply(lr, function(o) o[, 4], simplify=F)))
r <- cbind(r, "P-value"=rep("", nrow(r)))
r[match(names(pvs), r[, 1]), ncol(r)] <- base::format.pval(pvs)
r[is.na(r[, ncol(r)]) | gsub(" ", "", r[, ncol(r)])=='NA', ncol(r)] <- ""
colnames(r) <- c("Variab.", "N", "Groups", paste0("All<br>n=", nrow(d)),
                 paste0(levels(d$group), "<br>n=", ngg), "P-value")
r <- rbind(colnames(r), rep("", ncol(r)), r)



### HTML summary

r <- r[!r[, 3]%in%c("Male", "No"), ]
r[r[, 1]=="Sex", 1] <- paste0(r[r[, 1]=="Sex", c(1, 3)], collapse=' - ')
r[r[, 1]=="Center", 1] <- paste0(r[r[, 1]=="Center", c(1, 3)], collapse=' - ')
r <- r[, -3]
r[regexpr("NA", r) > 0] <- ""
r[regexpr("NaN", r) > 0] <- ""
r <- r[apply(r, 1, function(o) any(o!='')), ]

p <- openPage(paste0(.res,  "Descriptives_all.html"))
hwrite(paste0('<br><br>Descriptives - Univariate association with condition group<br><br>',
              "Measures: Median(Minimum, Maximum) for continouous variables<br>",
              "N(%) for categorical variables<br><br>"), p)
hwrite(r, page=p, center=F, row.names=F, col.names=F,
       col.style=c("text-align:center", rep("text-align:center", ncol(r)-1)),
       col.width=c("250px", "100px", rep('200px', ncol(r)-3)))
closePage(p)


#################################################################################################################
# Descriptives - All - No cohort sizes

rm(list=ls())
d <- readRDS(paste0(.dat, "0_metadata/data.all.rds"))
source("code/functions/clin_descrip.R")

ny <- "group"
nxs <- c("center", "sex", "age", 
         "previous.covid", "vaccinated", 
         "aht", "dm", "obesity", "dlp",
         "ferritin", "pcr", "dd", "ldh", "lymphocytes", "neutrophils", "leukocytes")
#         "cortis", "remdesivir", "toci",         
#         "bilateral.pneumony", "nimv.oaf", "icu", "exitus")
nv <- nxs
names(nv) <- c("Center", "Sex", "Age at diagnosis", 
               "Previous SARS-CoV-2 infection", "Vaccinated against SARS-CoV-2",
               "Arterial Hypertension", "Diabetes Mellitus", "Obesity", "Dyslipidemia",
               "Ferritin", "C-Reactive Protein (CRP)", "D-Dimer (DD)", "Lactate Dehydrogenase (LDH)",               
               "Lymphocytes", "Neutrophils", "Leukocytes")
#               "Glucocorticoids", "Remdesivir", "Tocilizumab",
#               "Bilateral Pneumony", "Non-Invasive Mechanical Ventilation (NIMV) / High Flow Oxygen Therapy (HFO)",
#               "Intensive Care Unit (ICU)", "Exitus")


for (o in nv[!nv%in%c("center", "series")])
{
  if (is.factor(d[, o])) d[, o] <- relevel(d[, o], ref=levels(d[, o])[2])
}



### P-values

set.seed(513335)
pvs <- sapply(nv, function(o)
{
  s <- rep(TRUE, nrow(d))
  y <- d[s, o]
  x <- droplevels(d[s, "group"])
  s <- x != "Control"
  y <- y[s]
  x <- droplevels(x[s])
  if (is.numeric(y))
  {
    r <- kruskal.test(y ~ x)$p.value
  }else if (is.factor(y))
  {
    r <- try(fisher.test(table(y, x), simulate.p.value=TRUE, B=10000)$p.value, silent=T)        
    if (class(r)=='try-error') r <- NA
  }else if (o=="treatment")
  {
    r <- NA
  }
  r
})


### All samples

dir.create(paste(.res, "all_all/", sep=''), F)
r <- clin.descrip(nv, d, ndir=paste0(.res, "all_all/"), bc=bcs)
ra <- r


### Levels

lr <- list()
for (j in 1:nlevels(d$group))
{
  s <- !is.na(d$group) & d$group==levels(d$group)[j]
  dir.create(paste(.res, levels(d$group)[j], "_all/", sep=''), F)
  r <- clin.descrip(nv, d[s, ], ndir=paste0(.res, levels(d$group)[j], "_all/"), bc=bcs)
  lr[[levels(d$group)[j]]] <- r
}


### Merge

ng <- summary(d$group)
ng <- ng[names(ng)!="NA's"]
ngg <- paste0(ng, "<br>(", formatC(ng/sum(ng)*100, format='f', dig=1), "%)")
r <- cbind(ra[, -4], All=ra[, 4], do.call(cbind, sapply(lr, function(o) o[, 4], simplify=F)))
r <- cbind(r, "P-value"=rep("", nrow(r)))
r[match(names(pvs), r[, 1]), ncol(r)] <- base::format.pval(pvs)
r[is.na(r[, ncol(r)]) | gsub(" ", "", r[, ncol(r)])=='NA', ncol(r)] <- ""
colnames(r) <- c("Variab.", "N", "Groups", paste0("All<br>n=", nrow(d)),
                 paste0(levels(d$group), "<br>n=", ngg), "P-value")
r <- rbind(colnames(r), rep("", ncol(r)), r)



### HTML summary

r <- r[!r[, 3]%in%c("Male", "No"), ]
r[r[, 1]=="Sex", 1] <- paste0(r[r[, 1]=="Sex", c(1, 3)], collapse=' - ')
r[grep("Center",r[,1])+1, 2] <- r[grep("Center",r[,1]), 2]
r[grep("Center",r[,1]), "P-value"] <- ""
r[grep("Center",r[,1])+1, 1] <- r[grep("Center",r[,1]), 1]
r[r[, 1]=="Center", 1] <- apply(r[r[, 1]=="Center", c(1, 3)],1,function(o){
  paste0(o, collapse=' - ')})
r <- r[, -3]
r[regexpr("NA", r) > 0] <- ""
r[regexpr("NaN", r) > 0] <- ""
r <- r[apply(r, 1, function(o) any(o!='')), ]

p <- openPage(paste0(.res,  "Descriptives_all_noCohorts.html"))
hwrite(paste0('<br><br>Descriptives - Univariate association with condition group<br><br>',
              "Measures: Median(Minimum, Maximum) for continouous variables<br>",
              "N(%) for categorical variables<br><br>"), p)
hwrite(r, page=p, center=F, row.names=F, col.names=F,
       col.style=c("text-align:center", rep("text-align:center", ncol(r)-1)),
       col.width=c("250px", "100px", rep('200px', ncol(r)-3)))
closePage(p)




#################################################################################################################
  # Descriptives of hospitalized patients

rm(list=ls())
d <- readRDS(paste0(.dat, "0_metadata/data.all.rds"))
source("code/functions/clin_descrip.R")


d <- d[d$group%in%c("Severe", "Critical"), ]
for (o in colnames(d))
    if (is.factor(d[, o])) d[, o] <- droplevels(d[, o])


ny <- "group"
nxs <- c("cortis", "remdesivir", "toci",         
         "bilateral.pneumony", "nimv.oaf", "icu", "exitus")
nv <- nxs
names(nv) <- c("Glucocorticoids", "Remdesivir", "Tocilizumab",
               "Bilateral Pneumony", "Non-Invasive Mechanical Ventilation (NIMV) / High Flow Oxygen Therapy (HFO)",
               "Intensive Care Unit (ICU)", "Exitus")


for (o in nv[!nv%in%c("center", "series")])
{
    if (is.factor(d[, o])) d[, o] <- relevel(d[, o], ref=levels(d[, o])[2])
}



### P-values

pvs <- sapply(nv, function(o)
{
    s <- rep(TRUE, nrow(d))
    y <- d[s, o]
    x <- droplevels(d[s, "group"])
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

dir.create(paste(.res, "all_hosp/", sep=''), F)
r <- clin.descrip(nv, d, ndir=paste0(.res, "all_hosp/"), bc=bcs)
ra <- r


### Levels

lr <- list()
for (j in 1:nlevels(d$group))
{
    s <- !is.na(d$group) & d$group==levels(d$group)[j]
    dir.create(paste(.res, levels(d$group)[j], "_hosp/", sep=''), F)
    r <- clin.descrip(nv, d[s, ], ndir=paste0(.res, levels(d$group)[j], "_hosp/"), bc=bcs)
    lr[[levels(d$group)[j]]] <- r
}


### Merge

ng <- summary(d$group)
ng <- ng[names(ng)!="NA's"]
ngg <- paste0(ng, "<br>(", formatC(ng/sum(ng)*100, format='f', dig=1), "%)")
r <- cbind(ra[, -4], All=ra[, 4], do.call(cbind, sapply(lr, function(o) o[, 4], simplify=F)))
r <- cbind(r, "P-value"=rep("", nrow(r)))
r[match(names(pvs), r[, 1]), ncol(r)] <- base::format.pval(pvs)
r[is.na(r[, ncol(r)]) | gsub(" ", "", r[, ncol(r)])=='NA', ncol(r)] <- ""
colnames(r) <- c("Variab.", "N", "Groups", paste0("All<br>n=", nrow(d)),
                 paste0(levels(d$group), "<br>n=", ngg), "P-value")
r <- rbind(colnames(r), rep("", ncol(r)), r)



### HTML summary

r <- r[!r[, 3]%in%c("Male", "No"), ]
r[r[, 1]=="Sex", 1] <- paste0(r[r[, 1]=="Sex", c(1, 3)], collapse=' - ')
r[r[, 1]=="Center", 1] <- paste0(r[r[, 1]=="Center", c(1, 3)], collapse=' - ')
r <- r[, -3]
r[regexpr("NA", r) > 0] <- ""
r[regexpr("NaN", r) > 0] <- ""
#r <- r[, -c(2:3)]
r <- r[apply(r, 1, function(o) any(o!='')), ]


p <- openPage(paste0(.res,  "Descriptives_hosp.html"))
hwrite(paste0('<br><br>Descriptives - Univariate association with condition group',
              'Hospitalized patients<br><br>',
              "Measures: Median(Minimum, Maximum) for continouous variables<br>",
              "N(%) for categorical variables<br><br>"), p)
hwrite(r, page=p, center=F, row.names=F, col.names=F,
       col.style=c("text-align:center", rep("text-align:center", ncol(r)-1)),
       col.width=c("350px", rep('150px', ncol(r)-1)))
closePage(p)


#################################################################################################################
  # Descriptives - all, blood determinations

rm(list=ls())
d <- readRDS(paste0(.dat, "0_metadata/data.all.rds"))
source("code/functions/clin_descrip.R")

ny <- "group"
nxs <- c("series", "ferritin", "pcr", "dd", "ldh", "lymphocytes", "neutrophils", "leukocytes")
nv <- nxs
names(nv) <- c("Series", "Ferritin", "C-Reactive Protein (CRP)", "D-Dimer (DD)", "Lactate Dehydrogenase (LDH)",               
               "Lymphocytes", "Neutrophils", "Leukocytes")


#s <- apply(d[, nv[-1]], 1, function(o) all(!is.na(o)))
#d <- d[s, ]

#s1 <- apply(d[, c("ferritin", "pcr", "dd", "ldh")], 1, function(o) any(!is.na(o)))
#s2 <- apply(d[, c("lymphocytes", "neutrophils", "leukocytes")], 1, function(o) any(!is.na(o)))
#d <- d[s1 & s2, ]

d <- d[d$group!='Control', ]

for (o in colnames(d))
    if (is.factor(d[, o])) d[, o] <- droplevels(d[, o])

r <- t(sapply(nv[-1], function(o)
{
    tapply(d[, o], d$group, function(o) round(sum(!is.na(o))/length(o)*100))
}))

r <- t(sapply(nv[-1], function(o)
{
    tapply(d[, o], d$group, function(o) sum(!is.na(o)))
}))


nv <- nv[-1]

### P-values

pvs <- sapply(nv, function(o)
{
    s <- rep(TRUE, nrow(d))
    y <- d[s, o]
    x <- droplevels(d[s, "group"])
    if (is.numeric(y))
    {
        r <- kruskal.test(y ~ x)$p.value
    }else if (is.factor(y))
    {
        r <- try(fisher.test(table(y, x), simulate.p.value=10000)$p.value, silent=T)
        if (class(r)=='try-error') r <- NA
    }else if (o=="treatment")
    {
        r <- NA
    }
    r
})


### All samples

dir.create(paste(.res, "all_parm/", sep=''), F)
r <- clin.descrip(nv, d, ndir=paste0(.res, "all_parm/"), bc=bcs)
ra <- r


### Levels

lr <- list()
for (j in 1:nlevels(d$group))
{
    s <- !is.na(d$group) & d$group==levels(d$group)[j]
    dir.create(paste(.res, levels(d$group)[j], "_parm/", sep=''), F)
    r <- clin.descrip(nv, d[s, ], ndir=paste0(.res, levels(d$group)[j], "_parm/"), bc=bcs)
    lr[[levels(d$group)[j]]] <- r
}


### Merge

ng <- summary(d$group)
ng <- ng[names(ng)!="NA's"]
ngg <- paste0(ng, "<br>(", formatC(ng/sum(ng)*100, format='f', dig=1), "%)")
r <- cbind(ra[, -4], All=ra[, 4], do.call(cbind, sapply(lr, function(o) o[, 4], simplify=F)))
r <- cbind(r, "P-value"=rep("", nrow(r)))
r[match(names(pvs), r[, 1]), ncol(r)] <- base::format.pval(pvs)
r[is.na(r[, ncol(r)]) | gsub(" ", "", r[, ncol(r)])=='NA', ncol(r)] <- ""
colnames(r) <- c("Variab.", "N", "Groups", paste0("All<br>n=", nrow(d)),
                 paste0(levels(d$group), "<br>n=", ngg), "P-value")
r <- rbind(colnames(r), rep("", ncol(r)), r)



### HTML summary

r <- r[!r[, 3]%in%c("Male", "No"), ]
r[r[, 1]=="Sex", 1] <- paste0(r[r[, 1]=="Sex", c(1, 3)], collapse=' - ')
r[r[, 1]=="Center", 1] <- paste0(r[r[, 1]=="Center", c(1, 3)], collapse=' - ')
r[r[, 1]=="Series", 1] <- paste0(r[r[, 1]=="Series", c(1, 3)], collapse=' - ')
r <- r[, -3]
r[regexpr("NA", r) > 0] <- ""
r[regexpr("NaN", r) > 0] <- ""
r <- r[apply(r, 1, function(o) any(o!='')), ]


ns <- t(sapply(nv, function(o)
{

    r <- tapply(d[, o], d$group, function(o) sum(!is.na(o)))
    r <- paste0(r, "<br>(", formatC(r/ng*100, format='f', dig=1), "%)")
    
}))

rr <- r[, 1:3]
for (j in 1:ncol(ns))
{
    rr <- cbind(rr, N=c("N obs. (%)", ns[, j]), r[, j+3])
}
rr <- cbind(rr, r[, (j+1+3):ncol(r)])

rr[-1, 2] <- paste0(as.numeric(rr[-1, 2]),
                    "<br>(", formatC(as.numeric(rr[-1, 2])/sum(ng)*100), "%)")


r <- rr

p <- openPage(paste0(.res,  "Descriptives_parm.html"))
hwrite(paste0('<br><br>Descriptives - Univariate association with condition group - blood determinations<br><br>',
              "Measures: Median(Minimum, Maximum) for continouous variables<br>",
              "N(%) for categorical variables<br><br>"), p)
hwrite(r, page=p, center=F, row.names=F, col.names=F,
       col.style=c("text-align:center", rep("text-align:center", ncol(r)-1)),
       col.width=c("250px", "100px", rep('200px', ncol(r)-3)))
closePage(p)



################################################################################################################
################################################################################################################
################################################################################################################




