#################################################################################################################
#################################################################################################################
### A.Read_data/A3.Descriptives_elisa.R
### Project: jcalvet_202110_marato
#################################################################################################################
#################################################################################################################
  # Libraries, paths and options

rm(list=ls())
.res <- paste0(.res0, "A.Read_data/A3.Descriptives_elisa/")
dir.create(.res, F, T)

#################################################################################################################
  # Descriptives for ELISA series - All
  # - No previous infection

rm(list=ls())
d <- readRDS(paste0(.dat, "0_metadata/data.all.rds"))
d <- d[d$elisa == "Yes",]
# nd <- readRDS(paste0(.dat, "1_DIA/nd.fp.rds"))
source("code/functions/clin_descrip.R")


#table(nd$study.id%in%d$study.id)
#table(d$study.id%in%nd$study.id)
#summary(d[d$study.id%in%ndd$patient.id, "group"])
#summary(d[d$study.id%in%nd$study.id, "group"])
#table(nd$study.id%in%ndd$patient.id)

# dd <- pData(nd)
# d <- d[!d$study.id%in%dd$patient.id, ]

ny <- "group"
nxs <- c("sex", "age",
         "center", "previous.covid", "vaccinated", 
         "aht", "dm", "obesity", "dlp",
         "ferritin", "pcr", "dd", "ldh", "lymphocytes", "neutrophils", "leukocytes")
#         "cortis", "remdesivir", "toci",         
#         "bilateral.pneumony", "nimv.oaf", "icu", "exitus")
nv <- nxs
names(nv) <- c("Sex", "Age at diagnosi",
               "Center", "Previous SARS-CoV-2 infection", "Vaccinated against SARS-CoV-2",
               "Arterial Hypertension", "Diabetes Mellitus", "Obesity", "Dyslipidemia",
               "Ferritin", "C-Reactive Protein (CRP)", "D-Dimer (DD)", "Lactate Dehydrogenase (LDH)",               
               "Lymphocytes", "Neutrophils", "Leukocytes")
#               "Glucocorticoids", "Remdesivir", "Tocilizumab",
#               "Bilateral Pneumony", "Non-Invasive Mechanical Ventilation (NIMV) / High Flow Oxygen Therapy (HFO)",
#               "Intensive Care Unit (ICU)", "Exitus")


d$center <- factor(d$center, levels=c("HUPT", "IRSICaixa"), labels=c("Parc Taulí University Hospital", "IRSICaixa"))
d$previous.covid <- factor(d$previous.covid, levels=c("No", "Yes"))
d$vaccinated <- factor(d$vaccinated, levels=c("No", "Yes"))
d$exitus <- factor(d$exitus, levels=c("No", "Yes"))

d$exitus[d$group!='Critical'] <- NA    ### one Mild exitus -> email from jcalvet in 2024-05-29

d$bilateral.pneumony <- factor(d$bilateral.pneumony, levels=c("No", "Yes"))

d$nimv.oaf <- factor(d$nimv.oaf, levels=c("No", "Yes"))
s <- !d$group%in%c("Critical")
d$nimv.oaf[s] <- NA
s <- d$group%in%c("Critical") & is.na(d$nimv.oaf)
d$nimv.oaf[s] <- "No"

d$icu <- factor(d$icu, levels=c("No", "Yes"))
s <- is.na(d$icu)
d$icu[s] <- "No"
s <- !d$group%in%c("Critical")
d$icu[s] <- NA

for (o in c("cortis", "remdesivir", "toci"))
{
    s <- is.na(d[, o]) & d$group%in%c("Severe", "Critical")
    d[s, o] <- "No"
}

#nxs <- c("lymphocytes", "neutrophils", "leukocytes", "ferritin", "pcr", "dd", "ldh")
#for (o in nx)
#    table(d$group, !is.na(d[, o]), exclude=NULL)
#for (o in nxs)
#    d[d$group=='Asymptomatic', o] <- NA


for (o in nv[!nv%in%c("center")])
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

p <- openPage(paste0(.res,  "Descriptives_elisa_all.html"))
hwrite(paste0('<br><br>Descriptives - Univariate association with condition group - ELISA<br><br>',
              "Measures: Median(Minimum, Maximum) for continouous variables<br>",
              "N(%) for categorical variables<br><br>"), p)
hwrite(r, page=p, center=F, row.names=F, col.names=F,
       col.style=c("text-align:center", rep("text-align:center", ncol(r)-1)),
       col.width=c("250px", "100px", rep('200px', ncol(r)-3)))
closePage(p)




################################################################################################################
################################################################################################################
################################################################################################################




