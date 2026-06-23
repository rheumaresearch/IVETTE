#################################################################################################################
#################################################################################################################
### write.html.mod.R
### 
#################################################################################################################
#################################################################################################################

write.html.mod <- function (x, links, tiny.pic, tiny.pic.size = 100, title = "",
    file, digits = 3, col.align='center', cellpadding=10, sorting=F)
{
    stopifnot(class(x) == "data.frame")
    if (missing(links))
        links <- vector("list", ncol(x))
    if (missing(tiny.pic))
        tiny.pic <- vector("list", ncol(x))
    stopifnot(class(links) == "list")
    stopifnot(class(tiny.pic) == "list")
    stopifnot(length(links) == ncol(x))
    stopifnot(length(tiny.pic) == ncol(x))
    stopifnot(!missing(file))
    column.class <- unlist(lapply(x, class))
    for (j in 1:ncol(x)) {
        if (column.class[j] == "factor")
            x[, j] <- as.character(x[, j])
        if (column.class[j] == "numeric")
            x[, j] <- round(x[, j], digits = digits)
    }
    cat("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">\n",
        sep = "", file = file)
    cat("<html>\n", file = file, append = T)
    cat("<body>\n", file = file, append = T)
    cat(paste("<CAPTION ALIGN=\"top\"><right><B>", title, "</B><br></center></CAPTION><BR>\n"),
        sep = "", file = file, append = T)
    cat(paste("<TABLE border=1 cellpadding=", cellpadding,">\n", sep=''), file = file, append = T)
    cat("<TR>\n", file = file, append = T)
    for (j in 1:ncol(x)) {
        cat("<TH>", file = file, append = T)
        cat(colnames(x)[j], file = file, append = T)
        cat("</TH>\n", file = file, append = T)
    }
    cat("</TR>\n", file = file, append = T)
    for (i in 1:nrow(x)) {
        cat("<TR>\n", file = file, append = T)
        for (j in 1:ncol(x)) {
            cat(paste("<TD align=", col.align, ">", sep=''), file = file, append = T)
            if (is.null(links[[j]]) & is.null(tiny.pic[[j]])) {
                cat(x[i, j], file = file, append = T)
            }
            else if (is.null(links[[j]]) & !is.null(tiny.pic[[j]])) {
                cat(paste("<A HREF=\"", links[[j]][[i]], "\"><img src=\"",
                  tiny.pic[[j]][[i]], "\" height=\"", tiny.pic.size,
                  "\" width=\"", tiny.pic.size, "\" /></A>",
                  sep = ""), file = file, append = T)
            }
            else if (!is.null(links[[j]]) & is.null(tiny.pic[[j]])) {
                cat(paste("<A HREF=\"", links[[j]][[i]], "\">",
                  x[i, j], "</A>", sep = ""), file = file, append = T)
            }
            else if (!is.null(links[[j]]) & !is.null(tiny.pic[[j]])) {
                cat(paste("<A HREF=\"", links[[j]][[i]], "\"><img src=\"",
                  tiny.pic[[j]][[i]], "\" height=\"", tiny.pic.size,
                  "\" width=\"", tiny.pic.size, "\" /></A>",
                  sep = ""), file = file, append = T)
            }
            cat("</TD>\n", file = file, append = T)
        }
        cat("</TR>\n", file = file, append = T)
    }
    cat("</TABLE>\n", file = file, append = T)
    cat("</body>\n", file = file, append = T)
    cat("</html>\n", file = file, append = T)
    if (sorting)
        sortDragHtmlTable(filename = file)
}


################################################################################################################


sortDragHtmlTable <- function (filename) 
{
    lastSlashPos <- gregexpr(.Platform$file.sep, filename)[[1]][length(gregexpr(.Platform$file.sep, 
        filename)[[1]])]
    outputdir <- ifelse(lastSlashPos == -1, getwd(), substr(filename, 
        0, lastSlashPos))
    fileCode <- ""
#    data(sorttable)
#    writeLines(sorttable, con = paste(outputdir, "sorttable.js", 
#        sep = .Platform$file.sep))
#    data(dragtable)
#    writeLines(dragtable, con = paste(outputdir, "dragtable.js", 
#        sep = .Platform$file.sep))
    tmpTxt <- readLines(filename)
#    tmpTxt[1] <- paste("<script src=\"sorttable.js\"></script>\n", 
#        tmpTxt[1])
#    tmpTxt[1] <- paste("<script src=\"dragtable.js\"></script>\n", 
#        tmpTxt[1])
#    tmpTxt <- sub("TABLE", "TABLE class=\"draggable sortable\"", 
#        tmpTxt)
    writeLines(tmpTxt, con = filename)
}



################################################################################################################
################################################################################################################
################################################################################################################




