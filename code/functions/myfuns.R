#################################################################################################################
#################################################################################################################
### myfuns.R: useful functions
### Author: Toni Berenguer
#################################################################################################################
#################################################################################################################

facc <- function(y, prs, probs, fn, th=0.5)
{

    tab <- table(y, prs)
    acc <- sum(diag(tab))/sum(tab)
    sens <- tab[2, 2]/sum(tab[2, ])
    spec <- tab[1, 1]/sum(tab[1, ])
    phi <- phi(tab, 3)
    roc.r <- roc(y, probs, quiet=T)
    auc <- as.numeric(roc.r$auc)
    accs <- c(acc=acc, sens=sens, spec=spec, phi=phi, auc=auc)
    th1 <- min(probs[prs==levels(prs)[2]])
    
    
    r <- coords(roc.r, x="best", best.method='closest.topleft', transpose = TRUE)
    th2 <- r[1]
    sens2 <- r[3]
    spec2 <- r[2]
    tab2 <- table(y, factor(as.numeric(probs > th), levels=0:1))    
    acc2 <- sum(diag(tab2))/sum(tab2)
    accs2 <- c(th=as.numeric(th2), acc=as.numeric(acc2), sens=as.numeric(sens2), spec=as.numeric(spec2))
    
    B <- 1000
    rb <- do.call(rbind, mclapply(1:B, function(j, y, prs, probs)
    {
      b1 <- sample(which(y==levels(y)[1]), replace=T)
      b2 <- sample(which(y==levels(y)[2]), replace=T)
      b <- c(b1, b2)
        tabb <- table(y[b], prs[b])
        accb <- sum(diag(tabb))/sum(tabb)
        sensb <- tabb[2, 2]/sum(tabb[2, ])
        specb <- tabb[1, 1]/sum(tabb[1, ])
        phib <- phi(tabb, 3)
        roc.rb <- roc(y[b], probs[b], quiet=T)
        aucb <- as.numeric(roc.rb$auc)
        c(acc=accb, sens=sensb, spec=specb, phi=phib, auc=aucb)
    }, y, prs, probs, mc.cores=8))

    cis <- t(apply(rb, 2, quantile, c(0.025, 0.975),na.rm=T))
    accs <- cbind(accs, cis[names(accs), ])

    rc <- apply(formatC(accs, format='f', dig=3), 1,
                function(o) paste0(o[1], " [", o[2], ", ", o[3], "]"))

    rc2 <- formatC(accs2, format='f', dig=3)

    
    png(fn, width=800, height=800)
    par(mar=c(7, 9, 4, 4))
    plot(1 - roc.r$specificities, roc.r$sensitivities, type='l',
         cex.axis=2, cex.lab=1.5, col='red4', lwd=4, cex.lab=2.5, axes=F,
         xlab='', sub=paste0("\nSpecificity"), ylab=paste0("Sensitivity\n"), cex.sub=2.5)
    axis(1, at=seq(0, 1, 0.1), labels=seq(1, 0, -0.1), lwd=2, cex.lab=1.5, cex.axis=1.8)
    axis(2, at=seq(0, 1, 0.1), labels=seq(0, 1, 0.1), lwd=2, cex.lab=1.5, cex.axis=1.8, las=2)
    box(lwd=2)
    legend(x='bottomright', inset=0.01, lty=0,
           legend=c(paste0("AUC = ", rc["auc"]),
                    "",
                    paste0("Threshold. = ", formatC(th, format='f', dig=3)),
                    paste0("Sens. = ", rc["sens"]),
                    paste0("Spec. = ", rc["spec"]),
                    paste0("Acc. = ", rc["acc"])),
           bty='n', cex=1.75)
    abline(h=sens, v=1-spec, lty=3, col='gray40', lwd=3)
    abline(a=0, b=1, lty=1, col='gray60', lwd=3)
    dev.off()

    list(accs=accs, rc=rc, accs.th=accs2)
}


facc2 <- function(y, prs, probs, fn, th=0.5)
{

    tab <- table(y, prs)
    acc <- sum(diag(tab))/sum(tab)
    sens <- tab[2, 2]/sum(tab[2, ])
    spec <- tab[1, 1]/sum(tab[1, ])
    phi <- phi(tab, 3)
    roc.r <- roc(y, probs, quiet=T)
    auc <- as.numeric(roc.r$auc)
    accs <- c(acc=acc, sens=sens, spec=spec, phi=phi, auc=auc)
    th1 <- min(probs[prs==levels(prs)[2]])
    

    r <- coords(roc.r, x="best", best.method='closest.topleft', transpose = TRUE)
    th2 <- r[1]
    sens2 <- r[3]
    spec2 <- r[2]
    tab2 <- table(y, factor(as.numeric(probs > th), levels=0:1))    
    acc2 <- sum(diag(tab2))/sum(tab2)
    accs2 <- c(th=as.numeric(th2), acc=as.numeric(acc2), sens=as.numeric(sens2), spec=as.numeric(spec2))
    
    B <- 1000
    rb <- do.call(rbind, mclapply(1:B, function(j, y, prs, probs)
    {
        b1 <- sample(which(y==levels(y)[1]), replace=T)
        b2 <- sample(which(y==levels(y)[2]), replace=T)
        b <- c(b1, b2)
        tabb <- table(y[b], prs[b])
        accb <- sum(diag(tabb))/sum(tabb)
        sensb <- tabb[2, 2]/sum(tabb[2, ])
        specb <- tabb[1, 1]/sum(tabb[1, ])
        phib <- phi(tabb, 3)
        roc.rb <- roc(y[b], probs[b], quiet=T)
        aucb <- as.numeric(roc.rb$auc)
        c(acc=accb, sens=sensb, spec=specb, phi=phib, auc=aucb)
    }, y, prs, probs, mc.cores=8))

    cis <- t(apply(rb, 2, quantile, c(0.025, 0.975)))
    accs <- cbind(accs, cis[names(accs), ])

    rc <- apply(formatC(accs, format='f', dig=3), 1,
                function(o) paste0(o[1], " [", o[2], ", ", o[3], "]"))

    rc2 <- formatC(accs2, format='f', dig=3)

    
    pdf(fn, width=9, height=8)
    par(mar=c(7, 9, 4, 4))
    plot(1 - roc.r$specificities, roc.r$sensitivities, type='l',
         cex.axis=2, cex.lab=1.5, col='red4', lwd=6, cex.lab=2, axes=F,
         xlab='', sub=paste0("\nSpecificity"), ylab=paste0("Sensitivity\n"), cex.sub=2)
    axis(1, at=seq(0, 1, 0.1), labels=seq(1, 0, -0.1), lwd=2, cex.lab=1.5, cex.axis=1.8)
    axis(2, at=seq(0, 1, 0.1), labels=seq(0, 1, 0.1), lwd=2, cex.lab=1.5, cex.axis=1.8, las=2)
    box(lwd=2)
    legend(x='bottomright', inset=0.01, lty=0,
           legend=c(paste0("AUC = ", rc["auc"])),
           bty='n', cex=1.75)
    abline(h=sens, v=1-spec, lty=3, col='gray40', lwd=3)
    abline(a=0, b=1, lty=1, col='gray60', lwd=3)
    dev.off()

    list(accs=accs, rc=rc, accs.th=accs2)
}


facc3 <- function(y, prs, probs, fn, th=0.5)
{

    tab <- table(y, prs)
    acc <- sum(diag(tab))/sum(tab)
    sens <- tab[2, 2]/sum(tab[2, ])
    spec <- tab[1, 1]/sum(tab[1, ])
    phi <- phi(tab, 3)
    roc.r <- roc(y, probs, quiet=T)
    auc <- as.numeric(roc.r$auc)
    accs <- c(acc=acc, sens=sens, spec=spec, phi=phi, auc=auc)
    th1 <- min(probs[prs==levels(prs)[2]])
    

    r <- coords(roc.r, x="best", best.method='closest.topleft', transpose = TRUE)
    th2 <- r[1]
    sens2 <- r[3]
    spec2 <- r[2]
    tab2 <- table(y, factor(as.numeric(probs > th), levels=0:1))    
    acc2 <- sum(diag(tab2))/sum(tab2)
    accs2 <- c(th=as.numeric(th2), acc=as.numeric(acc2), sens=as.numeric(sens2), spec=as.numeric(spec2))
    
    B <- 1000
    rb <- do.call(rbind, mclapply(1:B, function(j, y, prs, probs)
    {
        b1 <- sample(which(y==levels(y)[1]), replace=T)
        b2 <- sample(which(y==levels(y)[2]), replace=T)
        b <- c(b1, b2)
        tabb <- table(y[b], prs[b])
        accb <- sum(diag(tabb))/sum(tabb)
        sensb <- tabb[2, 2]/sum(tabb[2, ])
        specb <- tabb[1, 1]/sum(tabb[1, ])
        phib <- phi(tabb, 3)
        roc.rb <- roc(y[b], probs[b], quiet=T)
        aucb <- as.numeric(roc.rb$auc)
        c(acc=accb, sens=sensb, spec=specb, phi=phib, auc=aucb)
    }, y, prs, probs, mc.cores=8))

    cis <- t(apply(rb, 2, quantile, c(0.025, 0.975)))
    accs <- cbind(accs, cis[names(accs), ])

    rc <- apply(formatC(accs, format='f', dig=3), 1,
                function(o) paste0(o[1], " [", o[2], ", ", o[3], "]"))

    rc2 <- formatC(accs2, format='f', dig=3)

    
    pdf(fn, width=9, height=8)
    par(mar=c(7, 9, 4, 4))
    plot(1 - roc.r$specificities, roc.r$sensitivities, type='l',
         cex.axis=2, cex.lab=1.5, col='red4', lwd=6, cex.lab=2, axes=F,
         xlab='', sub=paste0("\nSpecificity"), ylab=paste0("Sensitivity\n"), cex.sub=2)
    axis(1, at=seq(0, 1, 0.1), labels=seq(1, 0, -0.1), lwd=2, cex.lab=1.5, cex.axis=1.8)
    axis(2, at=seq(0, 1, 0.1), labels=seq(0, 1, 0.1), lwd=2, cex.lab=1.5, cex.axis=1.8, las=2)
    box(lwd=2)
    legend(x='bottomright', inset=0.01, lty=0,
           legend=c(paste0("AUC = ", rc["auc"]),
                    "",
                    paste0("Threshold. = ", formatC(th, format='f', dig=3)),
                    paste0("Sens. = ", rc["sens"]),
                    paste0("Spec. = ", rc["spec"]),
                    paste0("Acc. = ", rc["acc"])),
           bty='n', cex=1.5)
    abline(h=sens, v=1-spec, lty=3, col='gray40', lwd=3)
    abline(a=0, b=1, lty=1, col='gray60', lwd=3)
    dev.off()

    list(accs=accs, rc=rc, accs.th=accs2)
}

#################################################################################################################
#################################################################################################################
#################################################################################################################
