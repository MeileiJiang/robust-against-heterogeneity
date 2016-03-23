nullirwsva = function (dat, n.sv, B = 5) 
{
  n <- ncol(dat)
  m <- nrow(dat)
  mod = as.matrix(rep(1, n), ncol = 1)
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

  rm(resid)
  cat(paste("Iteration (out of", B, "):"))
  for (i in 1:B) {
    mod.gam <- cbind(mod, uu$vectors[, 1:n.sv])
    mod0.gam <- cbind(mod)
    ptmp <- f.pvalue(dat, mod.gam, mod0.gam)
    pprob.gam <- (1 - edge.lfdr(ptmp))
    dats <- dat * pprob.gam
    uu <- eigen(t(dats) %*% dats)
    cat(paste(i, " "))
  }
  sv = svd(dats)$v[, 1:n.sv]
  retval <- list(sv = sv, pprob.gam = pprob.gam,  
                 n.sv = n.sv)
  return(retval)
}
