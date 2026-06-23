#################################################################################################################
#################################################################################################################
### R.For_drafts/R6.DIA_overall_view.R
### Project: jcalvet_202110_marato
#################################################################################################################
#################################################################################################################
  # Libraries, paths and options

rm(list=ls())
.res <- paste0(.res0, "G.Discovery_DIA_FP_MaxLFQ/R6.DIA_overall_view/")
dir.create(.res, F, T)

#################################################################################################################
### eof version of DIA experiment
rm(list=ls())
nd <- readRDS(paste0(.dat, "1_DIA/nd.fp.rds"))

pc <- eof(t(exprs(nd)), recursive=T, scaled=FALSE, nu=ncol(nd)-1)
Mr <- t(eofRecon(pc, uncenter=TRUE, unscale=TRUE))
rownames(Mr) <- featureNames(nd)
colnames(Mr) <- sampleNames(nd)

df <- fData(nd)
dp <- pData(nd)
ndr <- new('ExpressionSet', exprs=Mr, featureData=new('AnnotatedDataFrame', df),
           phenoData=new('AnnotatedDataFrame', dp), annotation="FP-MaxLFQ-eof")

saveRDS(ndr, paste0(.res, "nd.fp.eof.rds"))


#################################################################################################################
  # Merge differential expression results for PRM

rm(list=ls())
nd <- readRDS(paste0(.res, "nd.fp.eof.rds"))
# source("code/functions/assoc.pc.R")
source("code/functions/make_transparent.R")

pc <- prcomp(t(exprs(nd)))
pcvar <- pc$sdev^2/sum(pc$sdev^2)*100


nv <- c("1")
nvr <- c("quant.batch")
form <- paste("~ ", paste(c(nv, nvr), collapse=' + '))
design <- model.matrix(as.formula(form), nd)
formr <- paste("~ ", paste(nvr, collapse=' + '))
designr <- model.matrix(as.formula(formr), nd)[, -1, drop=F]
lmf <- lmFit(exprs(nd)[, rownames(design), drop=F], design)
betas <- coef(lmf)
exprs(nd) <- exprs(nd) - betas[, colnames(designr), drop=F]%*%t(designr)


### PC1 vs PC2

cols <- c("deepskyblue3", "orange3", "blueviolet", "red3")
pdf(paste0(.res, "pc_1_2_protein_corr.pdf"), width=8.5, height=8)
par(mar=c(5, 6, 5, 5))
plot(pc$x[, 1], pc$x[, 2], col=makeTransparent(cols[as.numeric(nd$group)], 0.7), pch=16, cex=2.25,
     xlim=c(-24, 24), ylim=c(-24, 24), 
     main='Principal components 1 and 2\nProtein level - corrected',
     xlab=paste0("PC1 (", formatC(pcvar[1], format='f', dig=1), "%)"),
     ylab=paste0("PC2 (", formatC(pcvar[2], format='f', dig=1), "%)"),
     axes=F, cex.lab=1.75)
axis(1, lwd=2.5, cex.axis=1.25)
axis(2, lwd=2.5, las=2, cex.axis=1.25)
abline(v=0, lty=1, lwd=4, col=makeTransparent('gray80', 0.1))
abline(h=0, lty=1, lwd=4, col=makeTransparent('gray80', 0.1))
box(lwd=2.5)
legend(x='topleft', inset=0.01, pch=16,
       col=rev(cols), pt.cex=2,
       legend=rev(levels(nd$group)), cex=1.25)
dev.off()




################################################################################################################
################################################################################################################
################################################################################################################




