####### EM Algorithm for Finite Mixture of MCNLMM #######
FMMCNLMM.EM = function(Data, y, X, Z, g=g, init.para, cor.type = c("UNC", "DEC", "CAR1", "ARp"), true.cls=NULL, tol = 1e-6, max.iter=500, per=10)
{
  begin = proc.time()[1]
#initial values of parameter
  w = init.para$w
  Beta = init.para$Beta[,1:g]
  DD = init.para$DD
  Sigma = init.para$Sigma
  Phi = init.para$Phi
  ga = init.para$ga
  nu = init.para$nu
  rho = init.para$rho

  r = ncol(Sigma[[1]])
  p = ncol(X)
  q = ncol(Z)
  n =  length(unique(Data$Subject))
  sj = numeric(n)
  for(j in 1: n) sj[j] = nrow(Data[Data$Subject==j, ])
  nj = sj * r
  cumsum.nj = cumsum(nj)
  cor.type = cor.type[1]
  N = length(y)
  TrZ = matrix(0, ncol=n*q, nrow=N)
  TrZ[1: cumsum.nj[1], 1: q] = Z[1: cumsum.nj[1], ]
  for(j in 2: n) TrZ[(cumsum.nj[j-1]+1): cumsum.nj[j], ((j-1)*q+1): (j*q)] = Z[(cumsum.nj[j-1]+1): cumsum.nj[j], ]

  Lam = TLam = TLam.inv = SigCi = TSigCi = TSigCi.inv = as.list(g)
  b = matrix(rnorm(n*q*g), n*q, g)
  wden = wrho = w1 = matrix(NA, n, g)
  xi = tau = Uxi = matrix(NA, n, g)
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
    mu = X[idx, ]%*%Beta[,i]
    if(length(idx) == 1){
     wden[j,i] = w[i] * (nu[i]*dnorm(y[idx], mu, sd=sqrt(TLam[[i]][idx, idx]/rho[i])) + (1 - nu[i])*dnorm(y[idx], mu, sd=sqrt(TLam[[i]][idx, idx])))
    } else wden[j,i] = w[i] * (nu[i]*dmvnorm(y[idx], mu, TLam[[i]][idx, idx]/rho[i]) + (1 - nu[i])*dmvnorm(y[idx], mu, TLam[[i]][idx, idx]))
  }}
  indv.den = rowSums(wden)
  lnL = loglik.old = sum(log(indv.den))
  y.cent = y - X %*% Beta
  iter = 0
  cat(paste(rep("=", 50), sep = "", collapse = ""), "\n")
  cat("Finite mixture of multivariate CN linear mixed model with ", g, "components and ", cor.type, " errors: ", "\n")
  cat("iter = ", iter, ",\t obs.loglik = ", loglik.old, sep = "", "\n")
  repeat
  {
#E-step
    iter = iter+1
    U = wden / indv.den
    Ni = colSums(U)
    w = Ni / n
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
     Phi[i] = optim(par = Phi[i], fn = phi.cn.loglik, method = "L-BFGS-B", lower = 1e-6, upper = 1-1e-6,
                    gup=i, y=y,  w=w, Beta=Beta, DD=DD, Sigma=Sigma, nu=nu, rho=rho, Phi=Phi, ga=ga, X=X, Z=Z, cor.type=cor.type, sj=sj, cumsum.nj=cumsum.nj)$par
     }
     if(cor.type == "DEC"){            # ECME
     DEC.est = optim(par = c(Phi[i], ga[i]), fn = phi.cn.loglik, method = "L-BFGS-B", lower = c(1e-6, 1e-6), upper = c(1-1e-6, Inf),
                    gup=i, y=y, w=w, Beta=Beta, DD=DD, Sigma=Sigma, nu=nu, rho=rho, Phi=Phi, ga=ga, X=X, Z=Z, cor.type=cor.type, sj=sj, cumsum.nj=cumsum.nj)$par
     Phi[i] = DEC.est[1]; ga[i] = DEC.est[2]
     }
     if(cor.type == "ARp"){         # ECM
     P = length(Phi[i])
     Phi[i] = optim(par = Phi[i], fn = phi.cn.loglik, method = "L-BFGS-B", lower = rep(-1 + 1e-6, P), upper = rep(1 - 1e-6, P),
                    gup=i, y=y, w=w, Beta=Beta, DD=DD, Sigma=Sigma, nu=nu, rho=rho, Phi=Phi, ga=ga, X=X, Z=Z, cor.type=cor.type, sj=sj, cumsum.nj=cumsum.nj)$par
     }
     if(cor.type == "UNC") Phi = rep(1e-6, g)
    }
#observed log-likelihood
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
     wden[j,i] = w[i] * (nu[i]*dnorm(y[idx], mu, sd=sqrt(TLam[[i]][idx, idx]/rho[i])) + (1 - nu[i])*dnorm(y[idx], mu, sd=sqrt(TLam[[i]][idx, idx])))
    } else wden[j,i] = w[i] * (nu[i]*dmvnorm(y[idx], mu, TLam[[i]][idx, idx]/rho[i]) + (1 - nu[i])*dmvnorm(y[idx], mu, TLam[[i]][idx, idx]))
    }}
    indv.den = rowSums(wden)
    loglik.new = sum(log(indv.den))
    lnL = c(lnL, loglik.new)
    diff = loglik.new - loglik.old
    para = c(w, as.vector(Beta))
    for(i in 1: g) para = c(para, DD[,,i][vech.D])
    for(i in 1: g) para = c(para, Sigma[[i]][vech.Sig])
    para = c(para, as.vector(Phi), ga, nu, rho)
    if(iter%%per == 0) cat("iter = ", iter, ",\t obs.loglik = ", loglik.new, ",\t diff = ", diff, ",\t w = ", w, ",\t nu = ", round(nu, 3), ",\t rho = ", round(rho, 3), ",\t Phi = ", round(Phi,3), ",\t ga = ", round(ga,3), sep = " ", "\n")
    if(diff < tol || iter >= max.iter) break
    loglik.old = loglik.new
  }
  end = proc.time()[1]
# Parameter estimation
  cat(rep("=", 20), "Mixture of Linear mixed modles with ", cor.type[1], " errors", rep("=", 20), sep = "", "\n")
  cat("It took", end - begin, "seconds.\n")
  cat("iter = ", iter, ",\t obs.loglik = ", loglik.new, sep = "", "\n")
  para.est = list(w = w, Beta = Beta, Sigma = Sigma, DD = DD, Phi = Phi, ga = ga, nu = nu, rho = rho, para = para)
  est.out = MI.fmmcnlmm(para.est, y, cor.type=cor.type)
# Model selection
  if(cor.type == "UNC"){ m = g * (p + q*(q+1)/2 + r*(r+1)/2 + 2) + (g-1)
  } else if(cor.type == "DEC"){ m = g * (p + q*(q+1)/2 + r*(r+1)/2 + 2) + length(as.vector(Phi)) + length(ga) + (g-1)
  } else{ m = g * (p + q*(q+1)/2 + r*(r+1)/2 + 2) + length(as.vector(Phi)) + (g-1)}
  aic = 2 * m - 2 * loglik.new
  bic = m * log(n) - 2 * loglik.new
  icl = bic - 2*sum(U*log(U+1e-300))
  awe = icl + 3*m + m*log(n)
  cat(paste(rep("=", 50), sep = "", collapse = ""), "\n")
  model.inf = list(loglik = loglik.new, iter.lnL = lnL, aic = aic, bic = bic, icl = icl, awe = awe, m = m)
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
  if(length(true.cls) != 0){
    CCR = 1 - classError(true.cls, post.clus)$errorRate
    ARI = adjustedRandIndex(true.cls, post.clus)
  } else{
      CCR = ARI = NULL
  }
  pre.cls = list(wden = wden, u.hat = U, post.cls = post.clus, CCR = CCR, ARI = ARI)
  return(list(model.inf = model.inf, para.est = para.est, est.out = est.out, b = b, xi = xihat, Yfit = Yfit, pre.cls = pre.cls))
}

phi.cn.loglik = function(phiga, gup, y, w, Beta, DD, Sigma, nu, rho, Phi, ga, X=X, Z=Z, cor.type, sj, cumsum.nj)
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
       wden[j,i] = w[i] * (nu[i]*dnorm(y[idx], mu, sd=sqrt(TLam[[i]][idx, idx]/rho[i])) + (1 - nu[i])*dnorm(y[idx], mu, sd=sqrt(TLam[[i]][idx, idx])))
      } else wden[j,i] = w[i] * (nu[i]*dmvnorm(y[idx], mu, TLam[[i]][idx, idx]/rho[i]) + (1 - nu[i])*dmvnorm(y[idx], mu, TLam[[i]][idx, idx]))
  }}
  phi.loglik = sum(log(rowSums(wden)))
  return(-phi.loglik)
}

# information matrix
MI.fmmcnlmm = function(para.est, y, cor.type=c("UNC", "DEC", "CAR1", "ARp"))
{
   w = para.est$w
   Beta = para.est$Beta
   DD = para.est$DD
   Sigma = para.est$Sigma
   Phi = para.est$Phi
   ga = para.est$ga
   nu = para.est$nu
   rho = para.est$rho

   g = ncol(Beta); p = nrow(Beta); q = nrow(DD[,,1]); r = nrow(Sigma[[1]]); P = length(Phi[1])
   n =  length(unique(Data$Subject))
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
   m = 1+ p + m1 + 2
   vech.D = vech.posi(q)
   vech.Sig = vech.posi(r)

   wden = wrho = w1 = matrix(NA, n, g)
   xi = tau = Uxi = matrix(NA, n, g)
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
      wden[j,i] = w[i] * (wrho[j,i] + w1[j,i])
      } else{ 
      wrho[j,i] = nu[i] * dmvnorm(y[idx], mu, TLam[[i]][idx, idx]/rho[i])
      w1[j,i] = (1 - nu[i]) * dmvnorm(y[idx], mu, TLam[[i]][idx, idx])
      wden[j,i] = w[i] * (wrho[j,i] + w1[j,i])
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
   Si[[i]][1,] = U[,i]/w[i] - U[,g]/w[g]
   for(j in 1: n){
   if(j == 1) idx = 1: cumsum.nj[1]
   else idx = (cumsum.nj[j-1]+1): cumsum.nj[j]
   Si[[i]][2:(p+1), j] = U[j, i] * tau[j, i] * t(matrix(X[idx, ], ncol=p)) %*% TLam.inv[[i]][idx, idx] %*% y.cent[idx, i]
   for(s in 1: m1) Si[[i]][(1+p+s), j] = (U[j, i] * tau[j,i] * sum(diag(y.cent[idx, i] %*% t(y.cent[idx, i]) %*% Linv.dotL[[s]][idx, idx] %*% TLam.inv[[i]][idx, idx])) - U[j, i] * sum(diag(Linv.dotL[[s]][idx, idx])) )/2
   Si[[i]][m, j] = Uxi[j,i] * (nj[j]/rho[i] - t(y.cent[idx, i]) %*% TLam.inv[[i]][idx, idx] %*% y.cent[idx, i]) / 2
   }
   Si[[i]][(m-1), ] = (Uxi[,i] -U[,i]*nu[i])/(nu[i]*(1-nu[i]))
   }
# Meilijson (1989) formula
   EIs = NULL
   for(i in 1: g) EIs = rbind(EIs, Si[[i]])
   FI = EIs[,1] %*% t(EIs[,1])
   for(j in 2: n) FI = FI + EIs[,j] %*% t(EIs[,j])

   Est = NULL
   if(cor.type == 'UNC'){
    Det = ((1+p+d1+d2+1): (1+p+d1+d2+P))
    for(i in 2:g) Det = c(Det, ((1+p+d1+d2+1): (1+p+d1+d2+P))+(i-1)*m)
    Det = c(Det, ((g-1)*m+1))
    Cov = try(solve(FI[-Det, -Det]), silent=F)
    for(i in 1: g) Est = c(Est, c(w[i], Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig], nu[i], rho[i]))
    Est = Est[-c((g-1)*m+1)]
   }
   if(cor.type == 'CAR1'){
    Det = (m-2)
    for(i in 2: g) Det = c(Det, (m-2)+(i-1)*m)
    Det = c(Det, ((g-1)*m+1))
    Cov = try(solve(FI[-Det, -Det]), silent=F)
    for(i in 1: g) Est = c(Est, c(w[i], Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig], Phi[i], nu[i], rho[i]))
    Est = Est[-c((g-1)*m+1)]
   }
   if(cor.type == "DEC"){
    Det = (g-1)*m+1
    Cov = try(solve(FI[-Det, -Det]), silent=F)
    for(i in 1: g) Est = c(Est, c(w[i], Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig], Phi[i], ga[i], nu[i], rho[i]))
    Est = Est[-Det]   
   }
   if(cor.type == "ARp"){
    Det = (g-1)*m+1
    Cov = try(solve(FI[-Det, -Det]), silent=F)
    for(i in 1: g) Est = c(Est, c(w[i], Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig], Phi[i], nu[i], rho[i]))
    Est = Est[-Det]
   }
   if(class(Cov)[1] == "try-error"){
    Sd = NA
   } else Sd = sqrt(diag(Cov))
   out = rbind(Est, Sd)
   return(list(FI=FI, Cov=Cov, out=out))
}
