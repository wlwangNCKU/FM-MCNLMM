####### EM Algorithm for Extended Finite Mixture of MTLMM #######
EFMMTLMM.EM = function(Data, X, Z, V, g=g, init.para, cor.type = c("UNC", "DEC", "CAR1", "ARp"), true.cls=NULL, tol = 1e-6, max.iter=500, comm.nu=F, per=100)
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
  Beta = init.para$Beta[,1:g]
  DD = init.para$DD
  Sigma = init.para$Sigma
  Phi = init.para$Phi
  ga = init.para$ga
  nu = init.para$nu
  psi = init.para$psi

  Lam = TLam = TLam.inv = SigCi = TSigCi = TSigCi.inv = as.list(g)
  b = matrix(rnorm(n*q*g), n*q, g)
  w = wden = matrix(NA, n, g)
  tau = kappa = matrix(NA, n, g)
  cumsum.q = cumsum(rep(q, n))
  Delta = matrix(NA, n, g)
  vech.D = vech.posi(q)
  vech.Sig = vech.posi(r)

#observed log-likelihood
  for(j in 1: n){
    eta = kronecker(diag(g-1), t(V[j,])) %*% psi
    w[j, ] = c(exp(eta), 1)/(1+sum(exp(eta[1:(g-1)])))
  }
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
    if(length(idx) == 1) wden[j,i] = w[j,i] * dt((y[idx]-X[idx, ]%*%Beta[,i])/sqrt(TLam[[i]][idx, idx]), df = nu[i], log=F) / sqrt(TLam[[i]][idx, idx])
    else wden[j,i] = w[j,i] * dmvt(y[idx], X[idx, ]%*%Beta[,i], TLam[[i]][idx, idx], df = nu[i], log=F)
  }}
  indv.den = rowSums(wden)
  lnL = loglik.old = sum(log(indv.den))
  y.cent = y - X %*% Beta
  iter = 0
  cat(paste(rep("=", 50), sep = "", collapse = ""), "\n")
  cat("Extended finite mixture of multivariate T linear mixed model with ", g, "components and ", cor.type, " errors: ", "\n")
  cat("iter = ", iter, ",\t obs.loglik = ", loglik.old, sep = "", "\n")
  repeat
  {
#E-step
    iter = iter+1
    U = wden / indv.den
    Ni = colSums(U)
    for(i in 1: g){
     Delta[1, i] = t(y.cent[1: cumsum.nj[1], i]) %*% TLam.inv[[i]][1: cumsum.nj[1], 1: cumsum.nj[1]] %*% y.cent[1: cumsum.nj[1], i]
     for(j in 2: n) Delta[j, i] = t(y.cent[(cumsum.nj[j-1]+1): cumsum.nj[j], i]) %*% TLam.inv[[i]][(cumsum.nj[j-1]+1): cumsum.nj[j], (cumsum.nj[j-1]+1): cumsum.nj[j]] %*% y.cent[(cumsum.nj[j-1]+1): cumsum.nj[j], i]
     tau[,i] = (nu[i]+nj)/(nu[i]+Delta[,i])
     kappa[,i] = digamma((nu[i]+nj)/2) - log((nu[i] + Delta[,i])/2)
    }

# CM-steps:
# psi
    S.psi = J.psi = 0
    for(j in 1: n){
       Kj = kronecker(diag(g-1), t(V[j,]))
       S.psi = S.psi + t(Kj) %*% (U[j,-g] - w[j,-g])
       J.psi = J.psi + t(Kj) %*% Kj
    }
    psi = psi + 4*solve(J.psi) %*% S.psi
    for(i in 1: g){
# Beta
     Beta[,i] = solve(t(rep(U[,i]*tau[,i], nj) * X) %*% TLam.inv[[i]] %*% X) %*% (t(rep(U[,i]*tau[,i], nj) * X) %*% TLam.inv[[i]] %*% y)    # ECM
#     Beta[,i] = solve(t(rep(U[,i]*tau[,i], nj)* X) %*% TSigCi.inv[[i]] %*% X) %*% (t(rep(U[,i]*tau[,i], nj) * X) %*% TSigCi.inv[[i]] %*% (y - TrZ %*%b[,i]))    # AECM
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
      Phi[i] = optim(par = Phi[i], fn = phi.t.Eloglik, method = "L-BFGS-B", lower = 1e-6, upper = 1-1e-6,
                    gup=i, y=y,  w=w, Beta=Beta, DD=DD, Sigma=Sigma, Phi=Phi, ga=ga, nu=nu, X=X, Z=Z, cor.type=cor.type, sj=sj, cumsum.nj=cumsum.nj)$par
     }
     if(cor.type == "DEC"){            # ECME
      DEC.est = optim(par = c(Phi[i], ga[i]), fn = phi.t.Eloglik, method = "L-BFGS-B", lower = c(1e-6, 1e-6), upper = c(1-1e-6, Inf),
                    gup=i, y=y, w=w, Beta=Beta, DD=DD, Sigma=Sigma, Phi=Phi, ga=ga, nu=nu, X=X, Z=Z, cor.type=cor.type, sj=sj, cumsum.nj=cumsum.nj)$par
      Phi[i] = DEC.est[1]; ga[i] = DEC.est[2]
     }
     if(cor.type == "ARp"){         # ECM
      P = length(Phi[i])
      Phi[i] = optim(par = Phi[i], fn = phi.t.Eloglik, method = "L-BFGS-B", lower = rep(-1 + 1e-6, P), upper = rep(1 - 1e-6, P),
                    gup=i, y=y, w=w, Beta=Beta, DD=DD, Sigma=Sigma, Phi=Phi, ga=ga, nu=nu, X=X, Z=Z, cor.type=cor.type, sj=sj, cumsum.nj=cumsum.nj)$par
     }
     if(cor.type == "UNC") Phi = rep(1e-6, g)
     if(comm.nu == F) #nu[i] = optim(par = nu[i], fn = nu.Q.fn, method = "L-BFGS-B", lower = 2, upper = 200, ycent=y.cent[,i], TLam.inv.i=TLam.inv[[i]], u=U[,i], nj=nj, cumsum.nj=cumsum.nj)$par
       nu[i] = optim(par = nu[i], fn = nu.Eloglik, method = "L-BFGS-B", lower = 2, upper = 200, gup=i, nu=nu, y=y, X=X, w=w, Beta=Beta, TLam=TLam, cumsum.nj=cumsum.nj)$par
    }
    if(comm.nu == T) nu = rep(optim(par = nu[1], fn = Enu.Eloglik, method = "L-BFGS-B", lower = 2, upper = 200, y=y, X=X, w=w, Beta=Beta, TLam=TLam, cumsum.nj=cumsum.nj)$par, g) 
#observed log-likelihood
    for(j in 1: n){
      eta = kronecker(diag(g-1), t(V[j,])) %*% psi
      w[j, ] = c(exp(eta), 1)/(1+sum(exp(eta[1:(g-1)])))
    }
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
      if(length(idx) == 1) wden[j,i] = w[j,i] * dt((y[idx]-X[idx, ]%*%Beta[,i])/sqrt(TLam[[i]][idx, idx]), df = nu[i], log=F) / sqrt(TLam[[i]][idx, idx])
      else wden[j,i] = w[j,i] * dmvt(y[idx], X[idx, ]%*%Beta[,i], TLam[[i]][idx, idx], df = nu[i], log=F)
    }}
    indv.den = rowSums(wden)
    loglik.new = sum(log(indv.den))
    lnL = c(lnL, loglik.new)
    diff = loglik.new - loglik.old
    para = c(as.vector(Beta))
    for(i in 1: g) para = c(para, DD[,,i][vech.D])
    for(i in 1: g) para = c(para, Sigma[[i]][vech.Sig])
    para = c(para, as.vector(Phi), ga, nu, psi)
#    if(iter%%per == 0) cat("iter = ", iter, ",\t obs.loglik = ", loglik.new, ",\t diff = ", diff, ",\t psi = ", psi, ",\t nu = ", round(nu, 3), ",\t Phi = ", round(Phi,3), ",\t ga = ", round(ga,3), sep = " ", "\n")
    if(diff < tol || iter >= max.iter) break
    loglik.old = loglik.new
  }
  end = proc.time()[1]
# Parameter estimation
  para.est = list(Beta = Beta, Sigma = Sigma, DD = DD, Phi = Phi, ga = ga, nu = nu, psi = psi, para = para)
  est.out = MI.efmmtlmm(para.est, cor.type=cor.type, comm.nu=comm.nu)

# Model selection
  v = ncol(V)*(g-1)
  if(cor.type == "UNC"){ m = g * (p + q*(q+1)/2 + r*(r+1)/2 + 1) + v
  } else if(cor.type == "DEC"){ m = g * (p + q*(q+1)/2 + r*(r+1)/2 + 1) + length(as.vector(Phi)) + length(ga) + v
  } else{ m = g * (p + q*(q+1)/2 + r*(r+1)/2 + 1) + length(as.vector(Phi)) + v } 
  if(comm.nu == T) m = m - g + 1
  aic = 2 * m - 2 * loglik.new
  bic = m * log(n) - 2 * loglik.new
  icl = bic - 2*sum(U*log(U+1e-300))
  awe = icl + 3*m + m*log(n)
  model.inf = list(loglik = loglik.new, iter.lnL = lnL, aic = aic, bic = bic, icl = icl, awe = awe, m = m, CPU=end-begin)
# Estimation
  Yest = X %*% Beta + TrZ %*% b
  Yfit1 = rowSums(apply(U, 2, rep, nj) * Yest)
  bb1 = matrix(rowSums(apply(U, 2, rep, each=q) * b), ncol=q, byrow=T)
  z = ifelse(t(apply(U, 1, order))==g, 1, 0)
  Yfit2 = rowSums(apply(z, 2, rep, nj) * Yest)
  Yfit = list(Yfit1 = Yfit1, Yfit2 = Yfit2)
  bb2 = matrix(rowSums(apply(z, 2, rep, each=q) * b), ncol=q, byrow=T)
  b = list(b1 = bb1, b2 = bb2)
  tau.hat1 = rowSums(U * tau)
  tau.hat2 = rowSums(z * tau)
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

phi.t.Eloglik = function(phiga, gup, y, w, Beta, DD, Sigma, Phi, ga, nu, X=X, Z=Z, cor.type, sj, cumsum.nj)
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
    if(length(idx) == 1) wden[j,i] = w[j,i] * dt((y[idx]-X[idx, ]%*%Beta[,i])/sqrt(TLam[[i]][idx, idx]), df = nu[i], log=F) / sqrt(TLam[[i]][idx, idx])
    else wden[j,i] = w[j,i] * dmvt(y[idx], X[idx, ]%*%Beta[,i], TLam[[i]][idx, idx], df = nu[i], log=F)
  }}
  phi.loglik = sum(log(rowSums(wden)))
  return(-phi.loglik)
}

nu.Eloglik = function(nu.i, gup, nu, y, X, w, Beta, TLam, cumsum.nj)
{
  nu[gup] = nu.i
  wden = matrix(NA, n, g)
  for(i in 1: g){for(j in 1: n){
    if(j == 1) idx = 1: cumsum.nj[1]
    else idx = (cumsum.nj[j-1]+1): cumsum.nj[j]
    if(length(idx) == 1) wden[j,i] = w[j,i] * dt((y[idx]-X[idx, ]%*%Beta[,i])/sqrt(TLam[[i]][idx, idx]), df = nu[i], log=F) / sqrt(TLam[[i]][idx, idx])
    else wden[j,i] = w[j,i] * dmvt(y[idx], X[idx, ]%*%Beta[,i], TLam[[i]][idx, idx], df = nu[i], log=F)
  }}
  nu.loglik = sum(log(rowSums(wden)))
  return(-nu.loglik)
}

Enu.Eloglik = function(nu, y, X, w, Beta, TLam, cumsum.nj)
{
  wden = matrix(NA, n, g)
  for(i in 1: g){for(j in 1: n){
    if(j == 1) idx = 1: cumsum.nj[1]
    else idx = (cumsum.nj[j-1]+1): cumsum.nj[j]
    if(length(idx) == 1) wden[j,i] = w[j,i] * dt((y[idx]-X[idx, ]%*%Beta[,i])/sqrt(TLam[[i]][idx, idx]), df = nu, log=F) / sqrt(TLam[[i]][idx, idx])
    else wden[j,i] = w[j,i] * dmvt(y[idx], X[idx, ]%*%Beta[,i], TLam[[i]][idx, idx], df = nu, log=F)
  }}
  nu.loglik = sum(log(rowSums(wden)))
  return(-nu.loglik)
}

MI.efmmtlmm = function(para.est, cor.type=c("UNC", "DEC", "CAR1", "ARp"), comm.nu=F)
{
   Beta = para.est$Beta
   DD = para.est$DD
   Sigma = para.est$Sigma
   Phi = para.est$Phi
   ga = para.est$ga
   nu = para.est$nu
   psi = para.est$psi

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
   m = p + m1 + 1
   vech.D = vech.posi(q)
   vech.Sig = vech.posi(r)

   Delta = matrix(NA, n, g)
   tau = ka1 = matrix(NA, n, g)
   w = wden = matrix(NA, n, g)
   for(j in 1: n){
     eta = kronecker(diag(g-1), t(V[j,])) %*% psi
     w[j, ] =  c(exp(eta), 1)/(1+sum(exp(eta[1:(g-1)])))
   }
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
     if(length(idx) == 1) wden[j,i] = w[j,i] * dt((y[idx]-X[idx, ]%*%Beta[,i])/sqrt(TLam[[i]][idx, idx]), df = nu[i], log=F) / sqrt(TLam[[i]][idx, idx])
     else wden[j,i] = w[j,i] * dmvt(y[idx], X[idx, ]%*%Beta[,i], TLam[[i]][idx, idx], df = nu[i], log=F)
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
   for(j in 1: n){
   if(j == 1) idx = 1: cumsum.nj[1]
   else idx = (cumsum.nj[j-1]+1): cumsum.nj[j]
   Si[[i]][1:p, j] = U[j, i] * tau[j, i] * t(matrix(X[idx, ], ncol=p)) %*% TLam.inv[[i]][idx, idx] %*% y.cent[idx, i]
   for(s in 1: m1) Si[[i]][(p+s), j] = (U[j, i] * tau[j,i] * sum(diag(y.cent[idx, i] %*% t(y.cent[idx, i]) %*% Linv.dotL[[s]][idx, idx] %*% TLam.inv[[i]][idx, idx])) - U[j, i] * sum(diag(Linv.dotL[[s]][idx, idx])) )/2
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
   v = ncol(V)*(g-1)
   S.psi = matrix(NA, nrow=v, ncol=n)
   for(j in 1: n){
     Kj = kronecker(diag(g-1), t(V[j,]))
     S.psi[,j] = t(Kj) %*% (U[j,-g] - w[j,-g])
   }
   EIs = rbind(EIs, S.psi)
   FI = EIs[,1] %*% t(EIs[,1])
   for(j in 2: n) FI = FI + EIs[,j] %*% t(EIs[,j])

   Est = NULL
   if(cor.type == 'UNC'){
    Det = ((p+d1+d2+1): (p+d1+d2+P))
    for(i in 2:g) Det = c(Det, ((p+d1+d2+1): (p+d1+d2+P))+(i-1)*m)
    if(comm.nu == T) Det = c(Det, m*(1:(g-1)))
    Cov = try(solve(FI[-Det, -Det]), silent=F)
    if(comm.nu == T){
        for(i in 1: g) Est = c(Est, c(Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig]))
        Est = c(Est, nu[1])
    } else{
      for(i in 1: g) Est = c(Est, c(Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig], nu[i]))
    }}
    if(cor.type == 'CAR1'){
    Det = (m-1)
    for(i in 2: g) Det = c(Det, (m-1)+(i-1)*m)
    if(comm.nu == T) Det = c(Det, m*(1:(g-1)))
    Cov = try(solve(FI[-Det, -Det]), silent=F)
    if(comm.nu == T){
        for(i in 1: g) Est = c(Est, c(Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig], Phi[i]))
        Est = c(Est, nu[1])
    } else{
    for(i in 1: g) Est = c(Est, c(Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig], Phi[i], nu[i]))
    }}
    if(cor.type == "DEC"){
    if(comm.nu == T){ Det = m*(1:(g-1))
    Cov = try(solve(FI[-Det, -Det]), silent=F)
    } else Cov = try(solve(FI), silent=F)
    if(comm.nu == T){
        for(i in 1: g) Est = c(Est, c(Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig], Phi[i], ga[i]))
        Est = c(Est, nu[1])
    } else{
    for(i in 1: g) Est = c(Est, c(Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig], Phi[i], ga[i], nu[i]))
    }}
    if(cor.type == "ARp"){
    if(comm.nu == T){ Det = m*(1:(g-1))
    Cov = try(solve(FI[-Det, -Det]), silent=F)
    } else Cov = try(solve(FI), silent=F)
    if(comm.nu == T){
        for(i in 1: g) Est = c(Est, c(Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig], Phi[i]))
        Est = c(Est, nu[1])
    } else{
     for(i in 1: g) Est = c(Est, c(Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig], Phi[i], nu[i]))
    }}
    if(class(Cov)[1] == "try-error"){
     Sd = NA
    } else Sd = sqrt(diag(Cov))
    out = rbind(c(Est, psi), Sd)
    return(list(FI=FI, Cov=Cov, out=out))
}
