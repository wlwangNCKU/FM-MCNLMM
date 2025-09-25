################################################################################
#
#   Filename: simExd.r
#   Purpose: conduct Simulation 2 and re-generate the intermediate results 
#   Input data files: Simulation/function/basicfn.r;
#                     Simulation/function/fmmlmm.fn.r;
#                     Simulation/function/fmmtlmm.fn.r;
#                     Simulation//function/fmmcnlmm.fn.r
#                     Simulation/function/efmmlmm.fn.r;
#                     Simulation/function/efmmtlmm.fn.r;
#                     Simulation//function/efmmcnlmm.fn.r
#   Output data files: Simulation/result/ESIM1/...;
#                      Simulation/result/ESIM2/...;
#                      Simulation/result/ESIM3/...;
#                      Simulation/result/ESIM4/...;
#                      Simulation/result/ESIM5/...;
#                      Simulation/result/ESIM6/...; 
#                      Simulation/result/ESIM7/...;
#                      Simulation/result/ESIM8/... 
#   Required R packages: nlme; mvtnorm; MASS
#
################################################################################

source(paste(WD.PATH, 'Simulation/function/basicfn.r', sep=''))
source(paste(WD.PATH, 'Simulation/function/fmmlmm.fn.r', sep=''))
source(paste(WD.PATH, 'Simulation/function/fmmtlmm.fn.r', sep=''))
source(paste(WD.PATH, 'Simulation/function/fmmcnlmm.fn.r', sep=''))
source(paste(WD.PATH, 'Simulation/function/efmmlmm.fn.r', sep=''))
source(paste(WD.PATH, 'Simulation/function/efmmtlmm.fn.r', sep=''))
source(paste(WD.PATH, 'Simulation/function/efmmcnlmm.fn.r', sep=''))

# basic setting
g=3; p=4; q=2; r=2; sj=7
vech.D = vech.posi(q)
vech.Sig = vech.posi(r)
Xi = matrix(c(rep(1, sj), 1: sj), nrow=sj)
Zi = matrix(Xi[, 1], nrow=sj)
rXi = kronecker(diag(r), Xi)
rZi = kronecker(diag(r), Zi)
X11 = rep(c(rep(1, sj), rep(0, sj)), n)
X12 = rep(c(1:sj, rep(0, sj)), n)
X21 = rep(c(rep(0, sj), rep(1, sj)), n)
X22 = rep(c(rep(0, sj), 1: sj), n)
X = cbind(X11, X12, X21, X22)
Z = X[, c(1, 3)]

# presummed parameters
Beta = matrix(c(1, 2, 2, -1, -1, 1, 0, -2, -2, 3, 1, -3), ncol=g)
DD = array(NA, dim=c(q,q,g))
for(i in 1: g)  DD[,,i] = matrix(c(2, 0.25, 0.25, 2), q, q, byrow=T)  # UNC
Sigma = as.list(g)
for(i in 1: g){
Sigma[[i]] = i*diag(c(1,2))
Sigma[[i]][1, 2] = Sigma[[i]][2, 1] = 0.5 * sqrt(prod(diag(Sigma[[i]][1:2, 1:2])))
}
Phi = rep(1e-6, g); ga = rep(1, g)  # UNC
para = list(Beta=Beta, DD=DD, Sigma=Sigma, nu=nu, rho=rho, Phi=Phi, ga=ga, psi=psi)

true.par = c(as.vector(Beta))
for(i in 1: g) true.par = c(true.par, DD[,,i][vech.D])
for(i in 1: g) true.par = c(true.par, Sigma[[i]][vech.Sig])
true.par = c(true.par, Phi, ga)
true.parComm = c(true.par, psi)
true.parC = c(true.par, nu, rho)

while(Rep <= total.rep)
{
cat(paste(c(rep('-', 15), rep(' ', 3), 'The ', Rep, ' time simulation: sample size = ', n, '; nu = ', nu[1], '; rho = ', rho[1], '; psi = ', psi, rep(' ', 3), rep('-', 15)), sep = ' ', collapse = ''), '\n')
# generate data 
genY = refmmcnlmm(n, para, cor.type='UNC', sj=sj, rXi=rXi, rZi=rZi)  
resp = cbind(as.vector(t(genY$Y[, 1: sj])),  as.vector(t(genY$Y[, (sj+1): (2*sj)])))
Time = rep((1: sj), n)
cls = rep(genY$cls, each=sj)
Subject = rep(1: n, each =sj)

sim.data = data.frame(cbind(Subject, Time, resp, cls))
colnames(sim.data) = c('Subject', 'Time', 'Var1', 'Var2', 'cls')
Data = groupedData(Var1+Var2 ~ Time | Subject, data = sim.data)
cls = genY$cls
V = genY$V

# fitting extended models
init.para = list(Beta=Beta+runif(1), DD=DD, Sigma=Sigma, Phi=Phi, ga=rep(1,g), nu=rep(30, g), psi=psi)
efit10 = try(EFMMLMM.EM(Data, X, Z, V, g=g, init.para, cor.type = "UNC", true.cls=cls, tol = 1e-3, max.iter=500), silent=F)
        if(class(efit10) == "try-error") next;      
efit20 = try(EFMMTLMM.EM(Data, X, Z, V, g=g, init.para, cor.type = "UNC", true.cls=cls, tol = 1e-3, max.iter=500), silent=F)
        if(class(efit20) == "try-error") next;    
init.para = list(Beta=Beta+runif(1), DD=DD, Sigma=Sigma, nu=rep(0.5, g), rho=rep(0.5, g), Phi=Phi, ga=rep(1, g), psi=psi)
efit30 = try(EFMMCNLMM.EM(Data, X, Z, V, g=g, init.para, cor.type = "UNC", true.cls=cls, tol = 1e-3, max.iter=500), silent=F)
        if(class(efit30) == "try-error") next;    

# fitting models
init.para = list(w=rep(1/g, g), Beta=Beta+runif(1), DD=DD, Sigma=Sigma, Phi=Phi, ga=rep(1,g), nu=rep(30, g))
fit10 = try(FMMLMM.EM(Data, X, Z, g=g, init.para, cor.type = "UNC", true.cls=cls, tol = 1e-3, max.iter=500), silent=F)
        if(class(fit10) == "try-error") next;    
fit20 = try(FMMTLMM.EM(Data, X, Z, g=g, init.para, cor.type = "UNC", true.cls=cls, tol = 1e-3, max.iter=500), silent=F)
        if(class(fit20) == "try-error") next;    
init.para = list(w=rep(1/g, g), Beta=Beta+runif(1), DD=DD, Sigma=Sigma, nu=rep(0.5, g), rho=rep(0.5, g), Phi=Phi, ga=rep(1, g))
fit30 = try(FMMCNLMM.EM(Data, X, Z, g=g, init.para, cor.type = "UNC", true.cls=cls, tol = 1e-3, max.iter=500), silent=F)
        if(class(fit30) == "try-error") next;    

# model selection
Tab1 = rbind(
c(efit10$model.inf$aic, efit20$model.inf$aic, efit30$model.inf$aic, fit10$model.inf$aic, fit20$model.inf$aic, fit30$model.inf$aic),
c(efit10$model.inf$bic, efit20$model.inf$bic, efit30$model.inf$bic, fit10$model.inf$bic, fit20$model.inf$bic, fit30$model.inf$bic),
c(efit10$model.inf$icl, efit20$model.inf$icl, efit30$model.inf$icl, fit10$model.inf$icl, fit20$model.inf$icl, fit30$model.inf$icl),
c(efit10$model.inf$awe, efit20$model.inf$awe, efit30$model.inf$awe, fit10$model.inf$awe, fit20$model.inf$awe, fit30$model.inf$awe))

Tab2 = cbind(
c(matrix(apply(Tab1[, 1:6], 1, order), nrow=6)[1,]),
c(matrix(apply(Tab1[, c(1:3,6)], 1, order), nrow=4)[1,]),
c(matrix(apply(Tab1[, c(1,3,6)], 1, order), nrow=3)[1,]),
c(matrix(apply(Tab1[, c(3,6)], 1, order), nrow=2)[1,]))

# classification 
Tab3 = rbind(
c(efit10$pre.cls$CCR, efit20$pre.cls$CCR, efit30$pre.cls$CCR, fit10$pre.cls$CCR, fit20$pre.cls$CCR, fit30$pre.cls$CCR),
c(efit10$pre.cls$ARI, efit20$pre.cls$ARI, efit30$pre.cls$ARI, fit10$pre.cls$ARI, fit20$pre.cls$ARI, fit30$pre.cls$ARI))

# parameter estimation
estEN = as.vector(t(efit10$est.out$out))
estET = as.vector(t(efit20$est.out$out))
estEC = as.vector(t(efit30$est.out$out))
 
# explore the numerical results
rname = c('AIC', 'BIC', 'ICL', 'AWE')
write.table(cbind(Rep, Tab1), paste(WD.PATH, 'Simulation/result/', pd, '/fit.txt', sep=""), append=T, row.names = rname, col.names = F)
write.table(cbind(Rep, Tab2), paste(WD.PATH, 'Simulation/result/', pd, '/selection.txt', sep=""), append=T, row.names = rname, col.names = F)
write.table(cbind(Rep, Tab3), paste(WD.PATH, 'Simulation/result/', pd, '/class.txt', sep=""), append=T, row.names = c('CCR', 'ARI'), col.names = F)
write(c(Rep, estEN), paste(WD.PATH, 'result/', pd, '/estEN.txt',sep=""), ncol=69, append=T)
write(c(Rep, estET), paste(WD.PATH, 'result/', pd, '/estET.txt',sep=""), ncol=75, append=T)
write(c(Rep, estEC), paste(WD.PATH, 'result/', pd, '/estEC.txt',sep=""), ncol=81, append=T)

Rep = Rep + 1
}
  