####### EM Algorithm for Finite Mixture of MTLMM #######
FMMTLMM.EM = function(Data, X, Z, g=g, init.para, cor.type = c("UNC", "DEC", "CAR1", "ARp"), true.cls=NULL, tol = 1e-6, max.iter=500, comm.nu=F, per=100)
{
  begin = proc.time()[1]
  p = ncol(X)
  q = ncol(Z)
  n =  length(unique(Data$Subject))
  sj = numeric(n)
  for(j in 1: n) sj[j] = nrow(Data[Data$Subject==j, ])
  nj = sj * r
  cumsum.nj = cumsum(nj)
  cor.type = cor.type[1]

  y = NULL
  if(r == 1){
   for(j in 1: n) y = c(y, Data$Var1[Data$Subject == j])
  } else{
   for(j in 1: n) y = c(y, Data$Var1[Data$Subject == j], Data$Var2[Data$Subject == j])
  } 
  N = length(y)
  TrZ = matrix(0, ncol=n*q, nrow=N)
  TrZ[1: cumsum.nj[1], 1: q] = Z[1: cumsum.nj[1], ]
  for(j in 2: n) TrZ[(cumsum.nj[j-1]+1): cumsum.nj[j], ((j-1)*q+1): (j*q)] = Z[(cumsum.nj[j-1]+1): cumsum.nj[j], ]

#initial values of parameter
  w = init.para$w
  Beta = matrix(init.para$Beta[,1:g], ncol=g)
  DD = init.para$DD
  Sigma = init.para$Sigma
  Phi = init.para$Phi
  ga = init.para$ga
  nu = init.para$nu

  Lam = TLam = TLam.inv = SigCi = TSigCi = TSigCi.inv = as.list(g)
  b = matrix(rnorm(n*q*g), n*q, g)
  wden = matrix(NA, n, g)
  tau = kappa = matrix(NA, n, g)
  cumsum.q = cumsum(rep(q, n))
  Delta = matrix(NA, n, g)
  vech.D = vech.posi(q)
  vech.Sig = vech.posi(r)

#observed log-likelihood
  for(i in 1: g){
    TLam[[i]] = TSigCi[[i]] = matrix(0, ncol=N, nrow=N)
    for(j in 1: n){
      if(j == 1)idx = 1: cumsum.nj[1]
      else idx = (cumsum.nj[j-1]+1): cumsum.nj[j]
      SigCi[[i]] = kronecker(Sigma[[i]], cor.fn(Phi[i], ga=ga[i], dim=sj[j], type=cor.type, Ti=Data$Time[Data$Subject == j]))
      Lam[[i]] = matrix(Z[idx, ], ncol=q) %*% DD[,,i] %*% t(matrix(Z[idx, ], ncol=q)) + SigCi[[i]]
      TSigCi[[i]][idx, idx] = SigCi[[i]]
      TLam[[i]][idx, idx] = Lam[[i]]
  }
    TLam.inv[[i]] = solve(TLam[[i]])
  }
  for(i in 1: g){for(j in 1: n){
    if(j == 1) idx = 1: cumsum.nj[1]
    else idx = (cumsum.nj[j-1]+1): cumsum.nj[j]
    if(length(idx) == 1) wden[j,i] = w[i] * dt((y[idx]-X[idx, ]%*%Beta[,i])/sqrt(TLam[[i]][idx, idx]), df = nu[i], log=F) / sqrt(TLam[[i]][idx, idx])
    else wden[j,i] = w[i] * dmvt(y[idx], X[idx, ]%*%Beta[,i], TLam[[i]][idx, idx], df = nu[i], log=F)
  }}
  indv.den = rowSums(wden)
  lnL = loglik.old = sum(log(indv.den))
  y.cent = y - X %*% Beta
  iter = 0
  cat(paste(rep("=", 50), sep = "", collapse = ""), "\n")
  cat("Finite mixture of multivariate T linear mixed model with ", g, "components and ", cor.type, " errors: ", "\n")
  cat("iter = ", iter, ",\t obs.loglik = ", loglik.old, sep = "", "\n")
  repeat
  {
#E-step
    iter = iter+1
    U = wden / indv.den
    Ni = colSums(U)
    w = Ni / n
    for(i in 1: g){
     Delta[1, i] = t(y.cent[1: cumsum.nj[1], i]) %*% TLam.inv[[i]][1: cumsum.nj[1], 1: cumsum.nj[1]] %*% y.cent[1: cumsum.nj[1], i]
     for(j in 2: n) Delta[j, i] = t(y.cent[(cumsum.nj[j-1]+1): cumsum.nj[j], i]) %*% TLam.inv[[i]][(cumsum.nj[j-1]+1): cumsum.nj[j], (cumsum.nj[j-1]+1): cumsum.nj[j]] %*% y.cent[(cumsum.nj[j-1]+1): cumsum.nj[j], i]
     tau[,i] = (nu[i]+nj)/(nu[i]+Delta[,i])
     kappa[,i] = digamma((nu[i]+nj)/2) - log((nu[i] + Delta[,i])/2)
    }
# CM-steps:
    for(i in 1: g){
# Beta
     Beta[,i] = solve(t(rep(U[,i]*tau[,i], nj) * X) %*% TLam.inv[[i]] %*% X) %*% (t(rep(U[,i]*tau[,i], nj) * X) %*% TLam.inv[[i]] %*% y)    # ECM
# DD
     TD = kronecker(diag(n), DD[,,i])
     b[,i] = TD %*% t(TrZ) %*% TLam.inv[[i]] %*% y.cent[,i]
     b.c = matrix(b[,i], ncol=n)
     TVb = solve(solve(TD) + t(TrZ) %*% solve(TSigCi[[i]]) %*% TrZ)
     b.c = matrix(b[,i], ncol=n)
     Sum.Vb = 0
     for(j in 1: n) Sum.Vb = Sum.Vb + U[j, i] * TVb[((j-1)*q+1):(j*q), ((j-1)*q+1):(j*q)]
     sum.B = (rep(U[,i]*tau[,i], each=q)*b.c) %*% t(b.c) + Sum.Vb
     DD[,,i] = sum.B / Ni[i]
    }
# Sigma
    y.cent = y - X %*% Beta
    e = y.cent - TrZ %*% b
    for(i in 1: g){
     ZVbZ = t(TrZ) %*% TLam.inv[[i]] %*% TrZ
     TD = kronecker(diag(n), DD[,,i])
     TD.inv = kronecker(diag(n), solve(DD[,,i]))
     e.c = NULL
     for(j in 2: n) e.c = c(e.c, rep(0, N), e[,i][(cumsum.nj[j-1]+1): cumsum.nj[j]])
     e.c = matrix(c(e[,i][1: cumsum.nj[1]], e.c), ncol=n)
     TPsi = rep(U[,i], nj) * (rep(tau[,i], nj) *e.c %*% t(e.c) + (TSigCi[[i]] - TSigCi[[i]] %*% TLam.inv[[i]] %*% TSigCi[[i]]))
     if(r==1)
     {
       Ce=0
       Ce=sum(solve(cor.fn(Phi[i], ga=ga[i], dim=sj[1], type=cor.type, Ti=Data$Time[Data$Subject == 1])) * TPsi[1: cumsum.nj[1], 1: cumsum.nj[1]])
       for(j in 2: n) Ce = Ce + sum(solve(cor.fn(Phi[i], ga=ga[i], dim=sj[1], type=cor.type, Ti=Data$Time[Data$Subject == j])) * TPsi[(cumsum.nj[j-1]+1): cumsum.nj[j], (cumsum.nj[j-1]+1): cumsum.nj[j]])
       Sigma[[i]] = Ce / sum(sj*U[,i])
     } else{
     for(v in 1: r)for(l in 1: r)
     {
       Ce=0
       Ce=sum(solve(cor.fn(Phi[i], ga=ga[i], dim=sj[1], type=cor.type, Ti=Data$Time[Data$Subject == 1]))*TPsi[1: cumsum.nj[1], 1: cumsum.nj[1]][((v-1)*sj[1]+1): (v*sj[1]), ((l-1)*sj[1]+1): (l*sj[1])])
       for(j in 2: n) Ce = Ce + sum(solve(cor.fn(Phi[i], ga=ga[i], dim=sj[j], type=cor.type, Ti=Data$Time[Data$Subject == j]))*TPsi[(cumsum.nj[j-1]+1): cumsum.nj[j], (cumsum.nj[j-1]+1): cumsum.nj[j]][((v-1)*sj[j]+1): (v*sj[j]), ((l-1)*sj[j]+1): (l*sj[j])])
       Sigma[[i]][v, l] = Ce / sum(sj*U[,i])
     }
     }
     if(cor.type == "CAR1"){             # ECME
     ga[i] = 1
     Phi[i] = optim(par = Phi[i], fn = phi.t.loglik, method = "L-BFGS-B", lower = 1e-6, upper = 1-1e-6,
                    gup=i, y=y,  w=w, Beta=Beta, DD=DD, Sigma=Sigma, Phi=Phi, ga=ga, nu=nu, X=X, Z=Z, cor.type=cor.type, sj=sj, cumsum.nj=cumsum.nj)$par
     }
     if(cor.type == "DEC"){            # ECME
     DEC.est = optim(par = c(Phi[i], ga[i]), fn = phi.t.loglik, method = "L-BFGS-B", lower = c(1e-6, 1e-6), upper = c(1-1e-6, Inf),
                    gup=i, y=y, w=w, Beta=Beta, DD=DD, Sigma=Sigma, Phi=Phi, ga=ga, nu=nu, X=X, Z=Z, cor.type=cor.type, sj=sj, cumsum.nj=cumsum.nj)$par
     Phi[i] = DEC.est[1]; ga[i] = DEC.est[2]
     }
     if(cor.type == "ARp"){         # ECM
     P = length(Phi[i])
     Phi[i] = optim(par = Phi[i], fn = phi.t.loglik, method = "L-BFGS-B", lower = rep(-1 + 1e-6, P), upper = rep(1 - 1e-6, P),
                    gup=i, y=y, w=w, Beta=Beta, DD=DD, Sigma=Sigma, Phi=Phi, ga=ga, nu=nu, X=X, Z=Z, cor.type=cor.type, sj=sj, cumsum.nj=cumsum.nj)$par
     }
     if(cor.type == "UNC") Phi = rep(1e-6, g)
     if(comm.nu == F) 
       if(g == 1){
       nu = optim(par = nu, fn = nu1.loglik, method = "L-BFGS-B", lower = 2, upper = 200, y=y, X=X, Beta=Beta, TLam=TLam, cumsum.nj=cumsum.nj)$par
       } else{
       nu[i] = optim(par = nu[i], fn = nu.loglik, method = "L-BFGS-B", lower = 2, upper = 200, gup=i, nu=nu, y=y, X=X, w=w, Beta=Beta, TLam=TLam, cumsum.nj=cumsum.nj)$par
    }}
    if(comm.nu == T) nu = rep(optim(par = nu[1], fn = Enu.loglik, method = "L-BFGS-B", lower = 2, upper = 200, y=y, X=X, w=w, Beta=Beta, TLam=TLam, cumsum.nj=cumsum.nj)$par, g) 
#observed log-likelihood
    for(i in 1: g){
      TLam[[i]] = TSigCi[[i]] = matrix(0, ncol=N, nrow=N)
      for(j in 1: n){
        if(j == 1)idx = 1: cumsum.nj[1]
        else idx = (cumsum.nj[j-1]+1): cumsum.nj[j]
        SigCi[[i]] = kronecker(Sigma[[i]], cor.fn(Phi[i], ga=ga[i], dim=sj[j], type=cor.type[1], Ti=Data$Time[Data$Subject == j]))
        Lam[[i]] = matrix(Z[idx, ], ncol=q) %*% DD[,,i] %*% t(matrix(Z[idx, ], ncol=q)) + SigCi[[i]]
        TSigCi[[i]][idx, idx] = SigCi[[i]]
        TLam[[i]][idx, idx] = Lam[[i]]
    }
    TLam.inv[[i]] = solve(TLam[[i]])
    }
    for(i in 1: g){for(j in 1: n){
      if(j == 1) idx = 1: cumsum.nj[1]
      else idx = (cumsum.nj[j-1]+1): cumsum.nj[j]
      if(length(idx) == 1) wden[j,i] = w[i] * dt((y[idx]-X[idx, ]%*%Beta[,i])/sqrt(TLam[[i]][idx, idx]), df = nu[i], log=F) / sqrt(TLam[[i]][idx, idx])
      else wden[j,i] = w[i] * dmvt(y[idx], X[idx, ]%*%Beta[,i], TLam[[i]][idx, idx], df = nu[i], log=F)
    }}
    indv.den = rowSums(wden)
    loglik.new = sum(log(indv.den))
    lnL = c(lnL, loglik.new)
    diff = loglik.new - loglik.old
    para = c(w, as.vector(Beta))
    for(i in 1: g) para = c(para, DD[,,i][vech.D])
    for(i in 1: g) para = c(para, Sigma[[i]][vech.Sig])
    para = c(para, as.vector(Phi), ga, nu)
#    if(iter%%per == 0) cat("iter = ", iter, ",\t obs.loglik = ", loglik.new, ",\t diff = ", diff, ",\t w = ", w, ",\t nu = ", round(nu, 3), ",\t Phi = ", round(Phi,3), ",\t ga = ", round(ga,3), sep = " ", "\n")
    if(diff < tol || iter >= max.iter) break
    loglik.old = loglik.new
  }
  end = proc.time()[1]
# Parameter estimation
  para.est = list(w = w, Beta = Beta, Sigma = Sigma, DD = DD, Phi = Phi, ga = ga, nu = nu, para = para)
  est.out = MI.fmmtlmm(para.est, cor.type=cor.type, comm.nu=comm.nu)

# Model selection
  if(cor.type == "UNC"){ m = g * (p + q*(q+1)/2 + r*(r+1)/2 + 1) + (g-1)
  } else if(cor.type == "DEC"){ m = g * (p + q*(q+1)/2 + r*(r+1)/2 + 1) + length(as.vector(Phi)) + length(ga) + (g-1)
  } else{ m = g * (p + q*(q+1)/2 + r*(r+1)/2 + 1) + length(as.vector(Phi)) + (g-1)}
  if(comm.nu == T) m = m - g + 1
  aic = 2 * m - 2 * loglik.new
  bic = m * log(n) - 2 * loglik.new
  icl = bic - 2*sum(U*log(U+1e-300))
  awe = icl + 3*m + m*log(n)
  model.inf = list(loglik = loglik.new, iter.lnL = lnL, aic = aic, bic = bic, icl = icl, awe = awe, m = m, CPU=end-begin)
# Estimation
  Yest = X %*% Beta + TrZ %*% b
  if(g == 1){
  Yfit1 = Yfit2 = Yest
  bb1 = bb2 = b
  tau.hat1 = tau.hat2 = tau 
  } else{
  Yfit1 = rowSums(apply(U, 2, rep, nj) * Yest)
  bb1 = matrix(rowSums(apply(U, 2, rep, each=q) * b), ncol=q, byrow=T)
  z = ifelse(t(apply(U, 1, order))==g, 1, 0)
  Yfit2 = rowSums(apply(z, 2, rep, nj) * Yest)
  bb2 = matrix(rowSums(apply(z, 2, rep, each=q) * b), ncol=q, byrow=T)
  tau.hat1 = rowSums(U * tau)
  tau.hat2 = rowSums(z * tau)
  }
  Yfit = list(Yfit1 = Yfit1, Yfit2 = Yfit2)
  b = list(b1 = bb1, b2 = bb2)
  tauhat = list(tau1 = tau.hat1, tau2 = tau.hat2)
# Cluster
  post.clus = matrix(apply(U, 1, order), nrow = g)[g,]
  if(length(true.cls) != 0 & g >= 2){
    CCR = 1 - classError(true.cls, post.clus)$errorRate
    ARI = adjustedRandIndex(true.cls, post.clus)
  } else{
      CCR = ARI = NULL
  } 
  pre.cls = list(wden = wden, u.hat = U, post.cls = post.clus, CCR = CCR, ARI = ARI)
  cat("loglik = ", loglik.new, ",\t AIC = ", aic, ",\t BIC = ", bic, ",\t ICL = ", icl, ",\t AWE = ", awe, ',\t CCR = ', CCR, ',\t ARI = ', ARI, sep = "", "\n")
  return(list(model.inf = model.inf, para.est = para.est, est.out = est.out, b = b, tau=tauhat, Yfit = Yfit, pre.cls = pre.cls))
}

phi.t.Qfn = function(phiga, ycent, u, Di, Sigma.i, nu.i, Z, cor.type, sj, cumsum.nj)
{
  if(cor.type == 'DEC'){ phi = phiga[1]; ga = phiga[2]}
  else{ phi = phiga; ga = 1}
  nj = sj * r
  Lam = t(t(Z[1:cumsum.nj[1], ])) %*% Di %*% t(Z[1:cumsum.nj[1], ])+ kronecker(Sigma.i, cor.fn(phi, ga=ga, dim=sj[1], type=cor.type, Ti=Data$Time[Data$Subject == 1]))
  Lam.inv = solve(Lam)
  Q1.phi.fn = u[1] * (log(det(Lam)) + (nu.i+nj[1])/2 * log(1+(t(ycent[1: cumsum.nj[1]])%*% Lam.inv %*% ycent[1: cumsum.nj[1]])/nu.i))
  for(j in 2: n){
    Lam = t(t(Z[(cumsum.nj[j-1]+1): cumsum.nj[j], ])) %*% Di %*% t(Z[(cumsum.nj[j-1]+1): cumsum.nj[j], ]) + kronecker(Sigma.i, cor.fn(phi, ga=ga, dim=sj[j], type=cor.type, Ti=Data$Time[Data$Subject == j]))
    Lam.inv = solve(Lam)
    Q1.phi.fn = Q1.phi.fn + u[j] * (log(det(Lam)) + (nu.i+sj[j])/2 * log(1+(t(ycent[(cumsum.nj[j-1]+1):cumsum.nj[j]])%*% Lam.inv %*% ycent[(cumsum.nj[j-1]+1): cumsum.nj[j]])/nu.i))
  }
  return(Q1.phi.fn)
}

phi.t.loglik = function(phiga, gup, y, w, Beta, DD, Sigma, Phi, ga, nu, X=X, Z=Z, cor.type, sj, cumsum.nj)
{
  if(cor.type == 'DEC'){ Phi[gup] = phiga[1]; ga[gup] = phiga[2]}
  else{ Phi[gup] = phiga; ga[gup] = 1}
  N = length(y)
  TLam = as.list(g)
  wden = matrix(NA, n, g)
  for(i in 1: g){
    TLam[[i]] = matrix(0, ncol=N, nrow=N)
    for(j in 1: n){
      if(j == 1)idx = 1: cumsum.nj[1]
      else idx = (cumsum.nj[j-1]+1): cumsum.nj[j]
      TLam[[i]][idx, idx] = matrix(Z[idx, ], ncol=q) %*% DD[,,i] %*% t(matrix(Z[idx, ], ncol=q)) + kronecker(Sigma[[i]], cor.fn(Phi[i], ga=ga[i], dim=sj[j], type=cor.type, Ti=Data$Time[Data$Subject == j]))
  }}
  for(i in 1: g){for(j in 1: n){
    if(j == 1) idx = 1: cumsum.nj[1]
    else idx = (cumsum.nj[j-1]+1): cumsum.nj[j]
    if(length(idx) == 1) wden[j,i] = w[i] * dt((y[idx]-X[idx, ]%*%Beta[,i])/sqrt(TLam[[i]][idx, idx]), df = nu[i], log=F) / sqrt(TLam[[i]][idx, idx])
    else wden[j,i] = w[i] * dmvt(y[idx], X[idx, ]%*%Beta[,i], TLam[[i]][idx, idx], df = nu[i], log=F)
  }}
  phi.loglik = sum(log(rowSums(wden)))
  return(-phi.loglik)
}

nu.Q.fn = function(nu.i, ycent, TLam.inv.i, u, nj, cumsum.nj)
{
 lndel = u[1] * (nu.i+nj[1])/2 * log(1+(t(ycent[1: cumsum.nj[1]])%*% TLam.inv.i[1:cumsum.nj[1], 1:cumsum.nj[1]] %*% ycent[1: cumsum.nj[1]])/nu.i)
 for(j in 2: n) lndel = lndel + u[j] * (nu.i+nj[j])/2 * log(1+(t(ycent[(cumsum.nj[j-1]+1):cumsum.nj[j]])%*% TLam.inv.i[(cumsum.nj[j-1]+1): cumsum.nj[j], (cumsum.nj[j-1]+1): cumsum.nj[j]] %*% ycent[(cumsum.nj[j-1]+1): cumsum.nj[j]])/nu.i)
 Q1.nu.fn = sum(u*(lgamma((nu.i+sj)/2) - lgamma(nu.i/2) - 0.5*nj*log(nu.i))) - lndel
 return(-Q1.nu.fn)
}

nu.loglik = function(nu.i, gup, nu, y, X, w, Beta, TLam, cumsum.nj)
{
  nu[gup] = nu.i
  wden = matrix(NA, n, g)
  for(i in 1: g){for(j in 1: n){
    if(j == 1) idx = 1: cumsum.nj[1]
    else idx = (cumsum.nj[j-1]+1): cumsum.nj[j]
    if(length(idx) == 1) wden[j,i] = w[i] * dt((y[idx]-X[idx, ]%*%Beta[,i])/sqrt(TLam[[i]][idx, idx]), df = nu[i], log=F) / sqrt(TLam[[i]][idx, idx])
    else wden[j,i] = w[i] * dmvt(y[idx], X[idx, ]%*%Beta[,i], TLam[[i]][idx, idx], df = nu[i], log=F)
  }}
  nu.loglik = sum(log(rowSums(wden)))
  return(-nu.loglik)
}

nu1.loglik = function(nu, y, X, Beta, TLam, cumsum.nj)
{
  i = 1
  wden = numeric(n)
  for(j in 1: n){
    if(j == 1) idx = 1: cumsum.nj[1]
    else idx = (cumsum.nj[j-1]+1): cumsum.nj[j]
    if(length(idx) == 1) wden[j] = dt((y[idx]-X[idx, ]%*%Beta[,i])/sqrt(TLam[[i]][idx, idx]), df = nu, log=F) / sqrt(TLam[[i]][idx, idx])
    else wden[j] = dmvt(y[idx], X[idx, ]%*%Beta[,i], TLam[[i]][idx, idx], df = nu, log=F)
  }
  nu.loglik = sum(log(wden))
  return(-nu.loglik)
}

Enu.loglik = function(nu, y, X, w, Beta, TLam, cumsum.nj)
{
  wden = matrix(NA, n, g)
  for(i in 1: g){for(j in 1: n){
    if(j == 1) idx = 1: cumsum.nj[1]
    else idx = (cumsum.nj[j-1]+1): cumsum.nj[j]
    if(length(idx) == 1) wden[j,i] = w[i] * dt((y[idx]-X[idx, ]%*%Beta[,i])/sqrt(TLam[[i]][idx, idx]), df = nu, log=F) / sqrt(TLam[[i]][idx, idx])
    else wden[j,i] = w[i] * dmvt(y[idx], X[idx, ]%*%Beta[,i], TLam[[i]][idx, idx], df = nu, log=F)
  }}
  nu.loglik = sum(log(rowSums(wden)))
  return(-nu.loglik)
}

MI.fmmtlmm = function(para.est, cor.type=c("UNC", "DEC", "CAR1", "ARp"), comm.nu=F)
{
   w = para.est$w
   Beta = para.est$Beta
   DD = para.est$DD
   Sigma = para.est$Sigma
   Phi = para.est$Phi
   ga = para.est$ga
   nu = para.est$nu

   g = ncol(Beta); p = nrow(Beta); q = nrow(DD[,,1]); r = nrow(Sigma[[1]]); P = length(Phi[1])
   n =  length(unique(Data$Subject))
   y = NULL
   for(j in 1: n) y = c(y, Data$Var1[Data$Subject == j], Data$Var2[Data$Subject == j])
   N = length(y)

   sj = numeric(n)
   for(j in 1: n) sj[j] = nrow(Data[Data$Subject==j, ])
   nj = sj * r
   cumsum.nj = cumsum(nj)
   cor.type = cor.type[1]

   TrZ = matrix(0, ncol=n*q, nrow=N)
   TrZ[1: cumsum.nj[1], 1: q] = Z[1: cumsum.nj[1], ]
   for(j in 2: n) TrZ[(cumsum.nj[j-1]+1): cumsum.nj[j], ((j-1)*q+1): (j*q)] = Z[(cumsum.nj[j-1]+1): cumsum.nj[j], ]
   d1 = q*(q+1)/2
   d2 = r*(r+1)/2
   if(cor.type == "DEC" | cor.type == "CAR1"){ m1 = d1 + d2 + 2
   } else m1 = d1 + d2 + P
   m = 1+ p + m1 + 1
   vech.D = vech.posi(q)
   vech.Sig = vech.posi(r)

   Delta = matrix(NA, n, g)
   tau = ka1 = matrix(NA, n, g)
   wden = matrix(NA, n, g)
   Lam = TLam = TLam.inv = SigCi = TSigCi = TSigCi.inv = as.list(g)
   for(i in 1: g){
     TLam[[i]] = TSigCi[[i]] = matrix(0, ncol=N, nrow=N)
       for(j in 1: n){
        if(j == 1) idx = 1: cumsum.nj[1]
        else idx = (cumsum.nj[j-1]+1): cumsum.nj[j]
        SigCi[[i]] = kronecker(Sigma[[i]], cor.fn(Phi[i], ga=ga[i], dim=sj[j], type=cor.type[1], Ti=Data$Time[Data$Subject == j]))
        Lam[[i]] = matrix(Z[idx, ], ncol=q) %*% DD[,,i] %*% t(matrix(Z[idx, ], ncol=q)) + SigCi[[i]]
        TSigCi[[i]][idx, idx] = SigCi[[i]]
        TLam[[i]][idx, idx] = Lam[[i]]
       }
     TLam.inv[[i]] = solve(TLam[[i]])
   }
   for(i in 1: g){for(j in 1: n){
     if(j == 1) idx = 1: cumsum.nj[1]
     else idx = (cumsum.nj[j-1]+1): cumsum.nj[j]
     if(length(idx) == 1) wden[j,i] = w[i] * dt((y[idx]-X[idx, ]%*%Beta[,i])/sqrt(TLam[[i]][idx, idx]), df = nu[i], log=F) / sqrt(TLam[[i]][idx, idx])
     else wden[j,i] = w[i] * dmvt(y[idx], X[idx, ]%*%Beta[,i], TLam[[i]][idx, idx], df = nu[i], log=F)
   }}
   indv.den = rowSums(wden)
   U = wden / indv.den
   ng = colSums(U)
   y.cent = y - X %*% Beta
   for(i in 1: g){
     Delta[1, i] = t(y.cent[1: cumsum.nj[1], i]) %*% TLam.inv[[i]][1: cumsum.nj[1], 1: cumsum.nj[1]] %*% y.cent[1: cumsum.nj[1], i]
     for(j in 2: n) Delta[j, i] = t(y.cent[(cumsum.nj[j-1]+1): cumsum.nj[j], i]) %*% TLam.inv[[i]][(cumsum.nj[j-1]+1): cumsum.nj[j], (cumsum.nj[j-1]+1): cumsum.nj[j]] %*% y.cent[(cumsum.nj[j-1]+1): cumsum.nj[j], i]
     tau[,i] = (nu[i]+nj)/(nu[i]+Delta[,i])
     ka1[,i] = Delta[,i] / nu[i]
   }

   Si = as.list(g)
   for(i in 1: g){
   Tdot.L = array(0, dim=c(N,N,m1)) 
   for(l in 1: d1)
   {
     TZDZ = matrix(0, ncol=N, nrow=N)
     dot.DD = matrix(0, q, q)
     dot.DD[matrix(vech.D[l, ], 1)] = dot.DD[matrix(rev(vech.D[l, ]), 1)] = 1
     TZDZ[1: cumsum.nj[1], 1: cumsum.nj[1]] = Z[1: cumsum.nj[1], ] %*% dot.DD %*% t(Z[1: cumsum.nj[1], ])
     for(j in 2: n)
     {
       rZj = Z[(cumsum.nj[j-1]+1): cumsum.nj[j], ]
       TZDZ[(cumsum.nj[j-1]+1): cumsum.nj[j], (cumsum.nj[j-1]+1): cumsum.nj[j]] = rZj %*% dot.DD %*% t(rZj)
     }
       Tdot.L[,,l] = TZDZ
   }
   for(l in 1: d2)
   {
     dot.Sig=matrix(0, r, r)
     dot.Sig[matrix(vech.Sig[l, ], 1)] = dot.Sig[matrix(rev(vech.Sig[l, ]), 1)] = 1
     Tdot.L[,,(d1+l)][1: cumsum.nj[1], 1: cumsum.nj[1]]=kronecker(dot.Sig, cor.fn(Phi[i], ga=ga[i], dim=sj[1], type=cor.type, Ti=Data$Time[Data$Subject == 1]))
     for(j in 2: n) Tdot.L[,, (d1+l)][(cumsum.nj[j-1]+1): cumsum.nj[j], (cumsum.nj[j-1]+1): cumsum.nj[j]]=kronecker(dot.Sig, cor.fn(Phi[i], ga=ga[i], dim=sj[j], type=cor.type, Ti=Data$Time[Data$Subject == j]))
   }
   if(cor.type == "DEC" | cor.type == "CAR1"){
     Tdot.L[,,(d1+d2+1)][1: cumsum.nj[1], 1: cumsum.nj[1]] = kronecker(Sigma[[i]], DEC.dot.phi(Phi[i], ga[i], Data$Time[Data$Subject == 1]))
     for(j in 2: n) Tdot.L[,,(d1+d2+1)][(cumsum.nj[j-1]+1): cumsum.nj[j], (cumsum.nj[j-1]+1): cumsum.nj[j]] = kronecker(Sigma[[i]], DEC.dot.phi(Phi[i], ga[i], Data$Time[Data$Subject == j]))
     Tdot.L[,,(d1+d2+2)][1: cumsum.nj[1], 1: cumsum.nj[1]] = kronecker(Sigma[[i]], DEC.dot.ga(Phi[i], ga[i], Data$Time[Data$Subject == 1]))
     for(j in 2: n) Tdot.L[,,(d1+d2+2)][(cumsum.nj[j-1]+1): cumsum.nj[j], (cumsum.nj[j-1]+1): cumsum.nj[j]] = kronecker(Sigma[[i]], DEC.dot.ga(Phi[i], ga[i], Data$Time[Data$Subject == j]))
   }
   else{
   for(l in 1: P){
     Tdot.L[,, (d1+d2+l)][1: cumsum.nj[1], 1: cumsum.nj[1]] = kronecker(Sigma[[i]], arp.Ci.dot(Phi[i], dim=sj[1], l))
     for(j in 2: n) Tdot.L[,, (d1+d2+l)][(cumsum.nj[j-1]+1): cumsum.nj[j], (cumsum.nj[j-1]+1): cumsum.nj[j]] = kronecker(Sigma[[i]], arp.Ci.dot(Phi[i], dim=sj[j], l))
   }}
   Linv.dotL = as.list(numeric(m1))
   for(s in 1: m1) Linv.dotL[[s]] = TLam.inv[[i]] %*% Tdot.L[,,s]

   Si[[i]] = matrix(NA, nrow=m, ncol=n)
   Si[[i]][1,] = U[,i]/w[i] - U[,g]/w[g]
   for(j in 1: n){
   if(j == 1) idx = 1: cumsum.nj[1]
   else idx = (cumsum.nj[j-1]+1): cumsum.nj[j]
   Si[[i]][2:(p+1), j] = U[j, i] * tau[j, i] * t(matrix(X[idx, ], ncol=p)) %*% TLam.inv[[i]][idx, idx] %*% y.cent[idx, i]
   for(s in 1: m1) Si[[i]][(1+p+s), j] = (U[j, i] * tau[j,i] * sum(diag(y.cent[idx, i] %*% t(y.cent[idx, i]) %*% Linv.dotL[[s]][idx, idx] %*% TLam.inv[[i]][idx, idx])) - U[j, i] * sum(diag(Linv.dotL[[s]][idx, idx])) )/2
   }
   if(comm.nu == T){
    Si[[i]][m, ] = U[,1] * (digamma((nu[1]+nj)/2) - digamma(nu[1]/2) - nj/nu[1] - log(1+ka1[,1]) + tau[,1]*ka1[,1]) / 2 
    for(k in 2: g) Si[[i]][m, ] = Si[[i]][m, ] + U[,k] * (digamma((nu[k]+nj)/2) - digamma(nu[k]/2) - nj/nu[k] - log(1+ka1[,k]) + tau[,k]*ka1[,k]) / 2 
   } else{
    Si[[i]][m, ] = U[,i] * (digamma((nu[i]+nj)/2) - digamma(nu[i]/2) - nj/nu[i] - log(1+ka1[,i]) + tau[,i]*ka1[,i]) / 2
   }}
# Meilijson (1989) formula
   EIs = NULL
   for(i in 1: g) EIs = rbind(EIs, Si[[i]])
   FI = EIs[,1] %*% t(EIs[,1])
   for(j in 2: n) FI = FI + EIs[,j] %*% t(EIs[,j])

   Est = NULL
   if(cor.type == 'UNC'){
    Det = ((1+p+d1+d2+1): (1+p+d1+d2+P))
    for(i in 2:g) Det = c(Det, ((1+p+d1+d2+1): (1+p+d1+d2+P))+(i-1)*m)
    if(comm.nu == T) Det = c(Det, m*(1:(g-1)))
    Det = c(Det, ((g-1)*m+1))
    Cov = try(solve(FI[-Det, -Det]), silent=F)
    if(comm.nu == T){
        for(i in 1: g) Est = c(Est, c(w[i], Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig]))
        Est = c(Est, nu[1])
    } else{
      for(i in 1: g) Est = c(Est, c(w[i], Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig], nu[i]))
    }
    Est = Est[-c((g-1)*m+1)]
    }
    if(cor.type == 'CAR1'){
    Det = (m-1)
    for(i in 2: g) Det = c(Det, (m-1)+(i-1)*m)
    if(comm.nu == T) Det = c(Det, m*(1:(g-1)))
    Det = c(Det, ((g-1)*m+1))
    Cov = try(solve(FI[-Det, -Det]), silent=F)
    if(comm.nu == T){
        for(i in 1: g) Est = c(Est, c(w[i], Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig], Phi[i]))
        Est = c(Est, nu[1])
    } else{
    for(i in 1: g) Est = c(Est, c(w[i], Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig], Phi[i], nu[i]))
    }
    Est = Est[-c((g-1)*m+1)]
    }
    if(cor.type == "DEC"){
    if(comm.nu == T){ Det = c(m*(1:(g-1)), ((g-1)*m+1))
    Cov = try(solve(FI[-Det, -Det]), silent=F)
    } else{
    Det = (g-1)*m+1
    Cov = try(solve(FI[-Det, -Det]), silent=F)
    }
    if(comm.nu == T){
        for(i in 1: g) Est = c(Est, c(w[i], Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig], Phi[i], ga[i]))
        Est = c(Est, nu[1])
    } else{
    for(i in 1: g) Est = c(Est, c(w[i], Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig], Phi[i], ga[i], nu[i]))
    }
    Est = Est[-Det]
    }
    if(cor.type == "ARp"){
    if(comm.nu == T){ Det = c(m*(1:(g-1)), ((g-1)*m+1))
    Cov = try(solve(FI[-Det, -Det]), silent=F)
    } else{
    Det = (g-1)*m+1
    Cov = try(solve(FI[-Det, -Det]), silent=F)
    }
    if(comm.nu == T){
        for(i in 1: g) Est = c(Est, c(w[i], Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig], Phi[i]))
        Est = c(Est, nu[1])
    } else{
     for(i in 1: g) Est = c(Est, c(w[i], Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig], Phi[i], nu[i]))
    }
    Est = Est[-Det]
    }
   if(class(Cov)[1] == "try-error"){
    Sd = NA
   } else Sd = sqrt(diag(Cov))
    out = rbind(Est, Sd)
    return(list(FI=FI, Cov=Cov, out=out))
}

# score vector
Sij.fmmtlmm = function(subj, wi, yi.cent, Delta.i, tau.i, TLami.inv, Linv.dotL, nu.i, X, m, m1, p, nj)
{
   cumsum.nj = cumsum(nj)
   if(subj == 1) idx = 1: cumsum.nj[1]
   else idx = (cumsum.nj[subj-1]+1): cumsum.nj[subj]
   sj = numeric(m)
   sj[1] = 1/wi
   sj[2:(p+1)] = tau.i[subj] * t(matrix(X[idx, ], ncol=p)) %*% TLami.inv[idx, idx] %*% yi.cent[idx]
   for(s in 1: m1) sj[(1+p+s)] = ( tau.i[subj]*(t(yi.cent[idx]) %*% Linv.dotL[[s]][idx, idx] %*% TLami.inv[idx, idx] %*% yi.cent[idx]) - sum(diag(Linv.dotL[[s]][idx, idx])) )/2
   sj[m] = (digamma((nu.i+nj[subj])/2) - digamma(nu.i/2) - nj[subj]/nu.i - log(1+Delta.i[subj]/nu.i) + tau.i[subj]*Delta.i[subj]/nu.i) / 2
   return(sj)
}
