#################################################################################################################
#################################################################################################################
### clin_descrip.R
###
### Functions for descriptive analysis of clinical variables
###
#################################################################################################################
#################################################################################################################


clin.descrip <- function(nv,d, ndir, bc=NULL)
{
	r <- sapply(1:length(nv), function(j, d, nv, bc)
                    clin.descrip.1(d[, nv[j]], names(nv)[j], ndir, nv[j], bc=bc),
                    d, nv, bc, simplify=F);
	r <- do.call(rbind, r);
        r <- r[, c(1, 3, 2, 4)]
	rownames(r) <- NULL;
	colnames(r) <- NULL;
	r;
}


clin.descrip.1 <- function(x, xlab, ndir, nx, bc=NULL)
{

    ns <- sum(!is.na(x))
    if (is.numeric(x) | is.integer(x))
    {
        r <- c(median(x, na.rm=T), quantile(x, c(0, 1), na.rm=T))
        r <- paste(formatC(r[1], format='f', dig=2), "<br>(",
                   formatC(r[2], format='f', dig=2), ", ", formatC(r[3], format='f', dig=2), ")", sep='')
        r <- c(xlab, "", ns, r);
        r <- t(r)
        r <- rbind(r, rep("", ncol(r)))
        if (sum(!is.na(x)) >= 3 )
        {
            png(paste(ndir, nx, ".png", sep=''), width=500, height=800)
            par(mfrow=c(2, 1))
            boxplot(x, main=paste("Boxplot\n", xlab, sep=''),
                    cex.main=1.2, pch='',  
                    col='lightblue', xlab="", ylab='', cex.lab=1.2)
            stripchart(x, vertical=T, method='jitter',
                       cex.main=1.2, 
                       xlab="", ylab="", cex.lab=1.2, add=T)
            s <- !is.na(x)
            dens <- density(x[s])
            plot(dens, main=paste("Density\n", xlab, sep=''), cex.main=1.2, cex.lab=1.2,
                 ylim=c(-max(dens$y)/100, max(dens$y)))
            abline(h=0)
            text(x=x[s], y=0, label='|', cex=1.2)
            dev.off()
        }
    }else if (is.factor(x))
    {
        fr <- table(x);
        tot <- sum(fr);
        pc <- fr/tot*100
        r <- paste(fr, "<br>(", formatC(as.numeric(pc), format='f', dig=1), "%)", sep='');
        r <- cbind(levels(x), c(ns, rep("", nlevels(x)-1)), r);
        r <- cbind(c(xlab, rep("", nrow(r)-1)), r)
        r <- rbind(r, rep("", ncol(r)))                
        png(paste(ndir, nx, ".png", sep=''), width=500, height=500)
        bp <- barplot(pc, ylim=c(0, 100), ylab='%', sub="", xlab=xlab, col='lightblue', 
                      main=paste("Barplot - ", xlab, sep=''),
                      cex.names=1.5, cex.main=1.5, cex.lab=1.2, plot=T)
        text(x=bp[, 1], y=pc + 3, label=fr, cex=1.5)
        dev.off()
    }
    r
}




#################################################################################################################
#################################################################################################################
#################################################################################################################

