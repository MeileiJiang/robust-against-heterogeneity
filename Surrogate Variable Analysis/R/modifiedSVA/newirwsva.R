# Iterative Reweighted Surogate Variable Analysis

newirwsva = function (dat, mod, mod0 = NULL, n.sv, B = 5) 
{
  n <- ncol(dat)
  m <- nrow(dat)
  if (is.null(mod0)) {
    mod0 <- mod[, 1]
  }
  Id <- diag(n)
  resid <- dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*% 
                      t(mod))
  uu <- eigen(t(resid) %*% resid)
  vv <- uu$vectors
  ndf <- n - dim(mod)[2]
  pprob <- rep(1, m)
  one <- rep(1, n)
  Id <- diag(n)
  df1 <- dim(mod)[2] + n.sv
  df0 <- dim(mod0)[2] + n.sv
  rm(resid)
  cat(paste("Iteration (out of", B, "):"))
  for (i in 1:B) {
    mod.b <- cbind(mod, uu$vectors[, 1:n.sv])
    mod0.b <- cbind(mod0, uu$vectors[, 1:n.sv])
    ptmp <- f.pvalue(dat, mod.b, mod0.b)
    pprob.b <- (1 - edge.lfdr(ptmp))
    mod.gam <- cbind(mod0, uu$vectors[, 1:n.sv])
    mod0.gam <- cbind(mod0)
    ptmp <- f.pvalue(dat, mod.gam, mod0.gam)
    pprob.gam <- (1 - edge.lfdr(ptmp))
    pprob <- pprob.gam * (1 - pprob.b)
    dats <- dat * pprob.gam
    resids <- dats %*% (Id - mod %*% solve(t(mod) %*% mod) %*% 
                        t(mod))
    uu <- eigen(t(resids) %*% resids)
    cat(paste(i, " "))
  }
  sv = svd(resids)$v[, 1:n.sv]
  retval <- list(sv = sv, pprob.gam = pprob.gam, pprob.b = pprob.b, 
                 n.sv = n.sv)
  return(retval)
}
