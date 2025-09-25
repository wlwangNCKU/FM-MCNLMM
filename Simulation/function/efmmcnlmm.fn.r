# generate random sample from truncated multivariate normal independent distribution
rmcn = function(n, mu, Sigma, nu=NULL, rho=NULL)
{
  p = length(mu)
  if(p == 1){
   Sd = sqrt(Sigma)
  } else{
  Lam = diag(sqrt(diag(Sigma)))
  Lam.inv = diag(1/sqrt(diag(Sigma)))
  R = Lam.inv %*%Sigma%*% Lam.inv
  }
  Y = matrix(NA, nrow=n, ncol=p)
  a = runif(n)
  tau = ifelse(a<=nu, rho, 1)
  tau.sq = sqrt(tau)
  if(p==1){
    Y = mu + 1/tau.sq * Sd * rnorm(n, mean=0, sd=1)
  } else Y = t(mu + 1/tau.sq * Lam %*% t(rmvnorm(n, mean=rep(0,p), sigma=R)))
  return(Y)
}

# Generate EFM-MCNLMM Data:
refmmcnlmm = function(n, para, cor.type, sj=sj, rXi=rXi, rZi=rZi)
{
   Beta = para$Beta
   DD = para$DD
   Sigma = para$Sigma
   Phi = para$Phi
   ga = para$ga
   nu = para$nu
   rho = para$rho
   psi = para$psi
   g = ncol(Beta)
   
   w = matrix(NA, n, g)
   V.samp = matrix(1, ncol=2, nrow=n)
   gs = sample(1:3, n, replace=T)
   for(j in 1: n){
     V.samp[j, ] = c(1, rnorm(1, gs[j]^2, sd=0.5*gs[j]))
     eta = kronecker(diag(g-1), t(V.samp[j,])) %*% psi
     w[j, ] = c(exp(eta), 1)/(1+sum(exp(eta[1:(g-1)])))
   }
   cls = rep(0, n)
   for(j in 1:n) cls[j] = which(rmultinom(n=1, size=1, prob=w[j,])==1)
   Y = NULL
   for(j in 1: n) Y = rbind(Y, rmcn(1, mu=c(rXi %*% Beta[,cls[j]]), Sigma=(rZi %*% DD[,,cls[j]] %*% t(rZi)+ kronecker(Sigma[[cls[j]]], cor.fn(Phi[cls[j]], ga=ga[cls[j]], dim=sj, type=cor.type[1], Ti=1:sj))), nu=nu[cls[j]], rho=rho[cls[j]]))

   Ymat = V = gcls = NULL
   for(i in 1: g){
     Ymat = rbind(Ymat, Y[which(cls==i), ])
     V = rbind(V, V.samp[which(cls==i), ])
     gcls = c(gcls, cls[which(cls==i)])
   }
   Data = list(Y=Ymat, cls=gcls, V=V)
   return(Data)
}

####### EM Algorithm for Extended Finite Mixture of MCNLMM #######
EFMMCNLMM.EM = function(Data, X, Z, V, g=g, init.para, cor.type = c("UNC", "DEC", "CAR1", "ARp"), true.cls=NULL, tol = 1e-6, max.iter=500, per=100)
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
  rho = init.para$rho
  g = ncol(Beta)

  Lam = TLam = TLam.inv = SigCi = TSigCi = TSigCi.inv = as.list(g)
  b = matrix(rnorm(n*q*g), n*q, g)
  w = wden = wrho = w1 = matrix(NA, n, g)
  xi = tau = Uxi = matrix(NA, n, g)
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
    mu = X[idx, ]%*%Beta[,i]
    if(length(idx) == 1){
     wden[j,i] = w[j,i] * (nu[i]*dnorm(y[idx], mu, sd=sqrt(TLam[[i]][idx, idx]/rho[i])) + (1 - nu[i])*dnorm(y[idx], mu, sd=sqrt(TLam[[i]][idx, idx])))
    } else wden[j,i] = w[j,i] * (nu[i]*dmvnorm(y[idx], mu, TLam[[i]][idx, idx]/rho[i]) + (1 - nu[i])*dmvnorm(y[idx], mu, TLam[[i]][idx, idx]))
  }}
  indv.den = rowSums(wden)
  lnL = loglik.old = sum(log(indv.den))
  y.cent = y - X %*% Beta
  iter = 0
  cat(paste(rep("=", 50), sep = "", collapse = ""), "\n")
  cat("Extended finite mixture of multivariate CN linear mixed model with ", g, "components and ", cor.type, " errors: ", "\n")
  cat("iter = ", iter, ",\t obs.loglik = ", loglik.old, sep = "", "\n")
  repeat
  {
#E-step
    iter = iter+1
    U = wden / indv.den
    Ni = colSums(U)
    for(i in 1: g){
    for(j in 1: n){
    if(j == 1) idx = 1: cumsum.nj[1]
    else idx = (cumsum.nj[j-1]+1): cumsum.nj[j]
    mu = X[idx, ]%*%Beta[,i]
    if(length(idx) == 1){
     wrho[j,i] = nu[i] * dnorm(y[idx], mu, sd=sqrt(TLam[[i]][idx, idx]/rho[i]))
     w1[j,i] = (1 - nu[i]) * dnorm(y[idx], mu, sd=sqrt(TLam[[i]][idx, idx]))
    } else{
     wrho[j,i] = nu[i] * dmvnorm(y[idx], mu, TLam[[i]][idx, idx]/rho[i])
     w1[j,i] = (1 - nu[i]) * dmvnorm(y[idx], mu, TLam[[i]][idx, idx])
    }}
     xi[,i] = wrho[,i] / (wrho[,i] + w1[,i])
     tau[,i] = xi[,i]* (rho[i] - 1) + 1
     Uxi[,i] = U[,i]*xi[,i]
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
     Beta[,i] = solve(t(rep(U[,i]*tau[,i], nj) * X) %*% TLam.inv[[i]] %*% X) %*% (t(rep(U[,i]*tau[,i], nj) * X) %*% TLam.inv[[i]] %*% y)
# DD
     TD = kronecker(diag(n), DD[,,i])
     b[,i] = TD %*% t(TrZ) %*% TLam.inv[[i]] %*% y.cent[,i]
     b.c = matrix(b[,i], ncol=n)
     TVb = solve(solve(TD) + t(TrZ) %*% solve(TSigCi[[i]]) %*% TrZ)
     b.c = matrix(b[,i], ncol=n)
     Sum.Vb = 0
     for(j in 1: n) Sum.Vb = Sum.Vb + U[j, i] * tau[j, i] * (1+(1/rho[i]-1)*xi[j,i]) * TVb[((j-1)*q+1):(j*q), ((j-1)*q+1):(j*q)]
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
     nu[i] = sum(Uxi[,i]) / sum(U[,i])

     deo = Uxi[1,i] * tr(TLam.inv[[i]][1: cumsum.nj[1], 1: cumsum.nj[1]] %*% (y.cent[1: cumsum.nj[1], i] %*% t(y.cent[1: cumsum.nj[1], i])))
     for(j in 2: n){
      idx = (cumsum.nj[j-1]+1): cumsum.nj[j]
      deo = deo + Uxi[j,i] * tr(TLam.inv[[i]][idx, idx] %*% (y.cent[idx, i] %*% t(y.cent[idx, i])))
     }
     rho[i] = min(1, sum(nj * Uxi[,i])/deo)
     
     if(cor.type == "CAR1"){             # ECME
     ga[i] = 1
     Phi[i] = optim(par = Phi[i], fn = phi.cn.Eloglik, method = "L-BFGS-B", lower = 1e-6, upper = 1-1e-6,
                    gup=i, y=y,  w=w, Beta=Beta, DD=DD, Sigma=Sigma, nu=nu, rho=rho, Phi=Phi, ga=ga, X=X, Z=Z, cor.type=cor.type, sj=sj, cumsum.nj=cumsum.nj)$par
     }
     if(cor.type == "DEC"){            # ECME
     DEC.est = optim(par = c(Phi[i], ga[i]), fn = phi.cn.Eloglik, method = "L-BFGS-B", lower = c(1e-6, 1e-6), upper = c(1-1e-6, Inf),
                    gup=i, y=y, w=w, Beta=Beta, DD=DD, Sigma=Sigma, nu=nu, rho=rho, Phi=Phi, ga=ga, X=X, Z=Z, cor.type=cor.type, sj=sj, cumsum.nj=cumsum.nj)$par
     Phi[i] = DEC.est[1]; ga[i] = DEC.est[2]
     }
     if(cor.type == "ARp"){         # ECM
     P = length(Phi[i])
     Phi[i] = optim(par = Phi[i], fn = phi.cn.Eloglik, method = "L-BFGS-B", lower = rep(-1 + 1e-6, P), upper = rep(1 - 1e-6, P),
                    gup=i, y=y, w=w, Beta=Beta, DD=DD, Sigma=Sigma, nu=nu, rho=rho, Phi=Phi, ga=ga, X=X, Z=Z, cor.type=cor.type, sj=sj, cumsum.nj=cumsum.nj)$par
     }
     if(cor.type == "UNC") Phi = rep(1e-6, g)
    }
#observed log-likelihood
    for(j in 1: n){
      eta = kronecker(diag(g-1), t(V[j,])) %*% psi
      w[j, ] = c(exp(eta), 1)/(1+sum(exp(eta[1:(g-1)])))
    }
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
    mu = X[idx, ]%*%Beta[,i]
    if(length(idx) == 1){
     wden[j,i] = w[j,i] * (nu[i]*dnorm(y[idx], mu, sd=sqrt(TLam[[i]][idx, idx]/rho[i])) + (1 - nu[i])*dnorm(y[idx], mu, sd=sqrt(TLam[[i]][idx, idx])))
    } else wden[j,i] = w[j,i] * (nu[i]*dmvnorm(y[idx], mu, TLam[[i]][idx, idx]/rho[i]) + (1 - nu[i])*dmvnorm(y[idx], mu, TLam[[i]][idx, idx]))
    }}
    indv.den = rowSums(wden)
    loglik.new = sum(log(indv.den))
    lnL = c(lnL, loglik.new)
    diff = loglik.new - loglik.old
    para = c(as.vector(Beta))
    for(i in 1: g) para = c(para, DD[,,i][vech.D])
    for(i in 1: g) para = c(para, Sigma[[i]][vech.Sig])
    para = c(para, as.vector(Phi), ga, nu, rho, psi)
#    if(iter%%per == 0) cat("iter = ", iter, ",\t obs.loglik = ", loglik.new, ",\t diff = ", diff, ",\t psi = ", psi, ",\t nu = ", round(nu, 3), ",\t rho = ", round(rho, 3), ",\t Phi = ", round(Phi,3), ",\t ga = ", round(ga,3), sep = " ", "\n")
    if(diff < tol || iter >= max.iter) break
    loglik.old = loglik.new
  }
  end = proc.time()[1]
# Parameter estimation
  para.est = list(Beta = Beta, Sigma = Sigma, DD = DD, Phi = Phi, ga = ga, nu = nu, rho = rho, psi = psi, para = para)
  est.out = MI.efmmcnlmm(para.est, cor.type=cor.type)
# Model selection
  v = ncol(V)*(g-1)
  if(cor.type == "UNC"){ m = g * (p + q*(q+1)/2 + r*(r+1)/2 + 2) + v
  } else if(cor.type == "DEC"){ m = g * (p + q*(q+1)/2 + r*(r+1)/2 + 2) + length(as.vector(Phi)) + length(ga) + v
  } else{ m = g * (p + q*(q+1)/2 + r*(r+1)/2 + 2) + length(as.vector(Phi)) + v}
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
  xi.hat1 = rowSums(U * xi)
  xi.hat2 = rowSums(z * xi)
  xihat = list(xi1 = xi.hat1, xi2 = xi.hat2)
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
  return(list(model.inf = model.inf, para.est = para.est, est.out = est.out, b = b, xi=xihat, Yfit = Yfit, pre.cls = pre.cls))
}

phi.cn.Eloglik = function(phiga, gup, y, w, Beta, DD, Sigma, nu, rho, Phi, ga, X=X, Z=Z, cor.type, sj, cumsum.nj)
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
      mu = X[idx, ]%*%Beta[,i]
      if(length(idx) == 1){
       wden[j,i] = w[j,i] * (nu[i]*dnorm(y[idx], mu, sd=sqrt(TLam[[i]][idx, idx]/rho[i])) + (1 - nu[i])*dnorm(y[idx], mu, sd=sqrt(TLam[[i]][idx, idx])))
      } else wden[j,i] = w[j,i] * (nu[i]*dmvnorm(y[idx], mu, TLam[[i]][idx, idx]/rho[i]) + (1 - nu[i])*dmvnorm(y[idx], mu, TLam[[i]][idx, idx]))
  }}
  phi.loglik = sum(log(rowSums(wden)))
  return(-phi.loglik)
}

# information matrix
MI.efmmcnlmm = function(para.est, cor.type=c("UNC", "DEC", "CAR1", "ARp"))
{
   Beta = para.est$Beta
   DD = para.est$DD
   Sigma = para.est$Sigma
   Phi = para.est$Phi
   ga = para.est$ga
   nu = para.est$nu
   rho = para.est$rho
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
   m = p + m1 + 2
   vech.D = vech.posi(q)
   vech.Sig = vech.posi(r)

   w = wden = wrho = w1 = matrix(NA, n, g)
   xi = tau = Uxi = matrix(NA, n, g)
   for(j in 1: n){
     eta = kronecker(diag(g-1), t(V[j,])) %*% psi
     w[j, ] = c(exp(eta), 1)/(1+sum(exp(eta[1:(g-1)])))
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
    if(j == 1){ idx = 1: cumsum.nj[1]
    } else idx = (cumsum.nj[j-1]+1): cumsum.nj[j]
      mu = X[idx, ]%*%Beta[,i]
      if(length(idx) == 1){
      wrho[j,i] = nu[i] * dnorm(y[idx], mu, sd=sqrt(TLam[[i]][idx, idx]/rho[i]))
      w1[j,i] = (1 - nu[i]) * dnorm(y[idx], mu, sd=sqrt(TLam[[i]][idx, idx]))
      wden[j,i] = w[j,i] * (wrho[j,i] + w1[j,i])
      } else{ 
      wrho[j,i] = nu[i] * dmvnorm(y[idx], mu, TLam[[i]][idx, idx]/rho[i])
      w1[j,i] = (1 - nu[i]) * dmvnorm(y[idx], mu, TLam[[i]][idx, idx])
      wden[j,i] = w[j,i] * (wrho[j,i] + w1[j,i])
     }
   }}
   indv.den = rowSums(wden)
   U = wden / indv.den
   ng = colSums(U)
   y.cent = y - X %*% Beta
   for(i in 1: g){
     xi[,i] = wrho[,i] / (wrho[,i] + w1[,i])
     tau[,i] = xi[,i]* (rho[i] - 1) + 1
     Uxi[,i] = U[,i]*xi[,i]
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
   Si[[i]][m, j] = Uxi[j,i] * (nj[j]/rho[i] - t(y.cent[idx, i]) %*% TLam.inv[[i]][idx, idx] %*% y.cent[idx, i]) / 2
   }
   Si[[i]][(m-1), ] = (Uxi[,i] -U[,i]*nu[i])/(nu[i]*(1-nu[i]))
   }
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
    Cov = try(solve(FI[-Det, -Det]), silent=F)
    for(i in 1: g) Est = c(Est, c(Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig], nu[i], rho[i]))
   }
   if(cor.type == 'CAR1'){
    Det = (m-2)
    for(i in 2: g) Det = c(Det, (m-2)+(i-1)*m)
    Cov = try(solve(FI[-Det, -Det]), silent=F)
    for(i in 1: g) Est = c(Est, c(Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig], Phi[i], nu[i], rho[i]))
   }
   if(cor.type == "DEC"){
    Cov = try(solve(FI), silent=F)
    for(i in 1: g) Est = c(Est, c(Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig], Phi[i], ga[i], nu[i], rho[i]))
   }
   if(cor.type == "ARp"){
    Cov = try(solve(FI), silent=F)
    for(i in 1: g) Est = c(Est, c(Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig], Phi[i], nu[i], rho[i]))
   }
   if(class(Cov)[1] == "try-error"){
    Sd = NA
   } else Sd = sqrt(diag(Cov))
   out = rbind(c(Est, psi), Sd)
   return(list(FI=FI, Cov=Cov, out=out))
}
