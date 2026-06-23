make_contrast_list <- function(conts, levg) {
  # IMPORTANT: parent must include base operators like + - /
  basis_env <- new.env(parent = baseenv())
  
  # basis vectors for each group level
  for (k in seq_along(levg)) {
    v <- rep(0, length(levg))
    v[k] <- 1
    # allow both "Mild" and "groupMild" symbols in contrast strings
    assign(levg[k], v, envir = basis_env)
    assign(paste0("group", levg[k]), v, envir = basis_env)
  }
  
  out <- vector("list", length(conts))
  names(out) <- names(conts)
  
  for (j in seq_along(conts)) {
    expr <- conts[j]
    # strip trailing "= 0" if present
    expr <- gsub("\\s*=\\s*0\\s*$", "", expr)
    # evaluate expression -> numeric vector length = #group levels
    out[[j]] <- eval(parse(text = expr), envir = basis_env)
  }
  
  out
}
