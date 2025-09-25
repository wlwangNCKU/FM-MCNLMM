################################################################################
#
#   Filename: sim.r
#   Purpose: conduct Simulation 1 and re-generate the intermediate results 
#   Input data files: Simulation/function/basicfn.r;
#                     Simulation/function/fmmlmm.fn.r;
#                     Simulation/function/fmmtlmm.fn.r;
#                     Simulation//function/fmmcnlmm.fn.r
#   Output data files: Simulation/result/SIM1/...;
#                      Simulation/result/SIM2/...;
#                      Simulation/result/SIM3/...;
#                      Simulation/result/SIM4/...;
#   Required R packages: nlme; mvtnorm; MASS
#
################################################################################

source(paste(WD.PATH, 'Simulation/function/basicfn.r', sep=''))
source(paste(WD.PATH, 'Simulation/function/fmmlmm.fn.r', sep=''))
source(paste(WD.PATH, 'Simulation/function/fmmtlmm.fn.r', sep=''))
source(paste(WD.PATH, 'Simulation/function/fmmcnlmm.fn.r', sep=''))

# basic setting
g=2; p=6; q=2; r=2
vech.D = vech.posi(q)
vech.Sig = vech.posi(r)
raX = read.table(paste(WD.PATH, 'Simulation/result/', x, '.txt', sep=""), header=T)

# presummed parameters
w = rep(1/g, g)
Beta = matrix(c(1, 2, 2, 2, -1.5, 1.5, -1, 1.5, 1, 0.5, -2, 0.5), ncol=g)
DD = array(NA, dim=c(q,q,g))
for(i in 1: g) DD[,,i] = matrix(c(2, 0.25, 0.25, 2), q, q, byrow=T)  # UNC
Sigma = as.list(g)
for(i in 1: g){
Sigma[[i]] = i * diag(c(1,2))
Sigma[[i]][1, 2] = Sigma[[i]][2, 1] = 0.5 * sqrt(prod(diag(Sigma[[i]][1:2, 1:2])))
}
Phi = rep(1e-6, g); ga = rep(1, g)  # UNC
para = list(w=w, Beta=Beta, DD=DD, Sigma=Sigma, nu=nu, rho=rho, Phi=Phi, ga=ga)
true.para = c(w, as.vector(Beta))
for(i in 1: g) true.para = c(true.para, DD[,,i][vech.D])
for(i in 1: g) true.para = c(true.para, Sigma[[i]][vech.Sig])
true.para = c(true.para, Phi, ga, nu, rho)

while(Rep <= total.rep)
{
cat(paste(c(rep('-', 15), rep(' ', 3), 'The ', Rep, ' time simulation: sample size = ', n, '; nu = ', nu, '; rho = ', rho, rep(' ', 3), rep('-', 15)), sep = '', collapse = ''), '\n')
# generate data 
genY = rfmmcnlmm(n, para, cor.type='UNC', raX=raX)
Data = genY$Data
X = genY$X; Z = genY$Z
cls = numeric(n)
for(j in 1: n) cls[j] = Data$cls[which(Data$Subjec == j)][1]

# model fitting
# g=1 #
init.para = list(w=1, Beta=Beta+runif(1), DD=DD, Sigma=Sigma, Phi=1e-6, ga=1, nu=30)
fitN1 = try(FMMLMM.EM(Data, X, Z, g=1, init.para, cor.type = "UNC", true.cls=cls, tol = 1e-3, max.iter=500), silent=F)
        if(class(fitN1) == "try-error") next;   
init.para = list(w=1, Beta=fitN1$para$Beta, DD=fitN1$para$DD, Sigma=fitN1$para$Sigma, Phi=1e-6, ga=1, nu=30)
fitT1 = try(FMMTLMM.EM(Data, X, Z, g=1, init.para, cor.type = "UNC", true.cls=cls, tol = 1e-3, max.iter=500), silent=F)
        if(class(fitT1) == "try-error") next;      
init.para = list(w=1, Beta=fitN1$para$Beta, DD=fitN1$para$DD, Sigma=fitN1$para$Sigma, nu=nu[1], rho=rho[1], Phi=1e-6, ga=1)
fitC1 = try(FMMCNLMM.EM(Data, X, Z, g=1, init.para, cor.type = "UNC", true.cls=cls, tol = 1e-3, max.iter=500), silent=F)
        if(class(fitC1) == "try-error") next;      

# g=2 #
init.para = list(w=w, Beta=Beta+runif(1), DD=DD, Sigma=Sigma, Phi=Phi, ga=rep(1,2), nu=rep(30, 2))
fitN2 = try(FMMLMM.EM(Data, X, Z, g=2, init.para, cor.type = "UNC", true.cls=cls, tol = 1e-3, max.iter=500), silent=F)
        if(class(fitN2) == "try-error") next;      
fitT2 = try(FMMTLMM.EM(Data, X, Z, g=2, init.para, cor.type = "UNC", true.cls=cls, tol = 1e-3, max.iter=500), silent=F)
        if(class(fitT2) == "try-error") next;      
init.para = list(w=w, Beta=Beta+runif(1), DD=DD, Sigma=Sigma, nu=nu, rho=rho, Phi=Phi, ga=rep(1, 2))
fitC2 = try(FMMCNLMM.EM(Data, X, Z, g=2, init.para, cor.type = "UNC", true.cls=cls, tol = 1e-3, max.iter=500), silent=F)
        if(class(fitC2) == "try-error") next;      

# g=3 #
Beta = matrix(c(1, 2, 2, 2, -1.5, 1.5, -1, 1.5, 1, 0.5, -2, 0.5, 0, 0, 0, 0, 0, 0), ncol=3)
DD = array(NA, dim=c(q,q,3))
for(i in 1: 3) DD[,,i] = i * matrix(c(2, 0.25, 0.25, 1), q, q, byrow=T)  # UNC
Sigma = as.list(3)
for(i in 1: 3){
Sigma[[i]] = i * diag(c(1,2))
Sigma[[i]][1, 2] = Sigma[[i]][2, 1] = 0.5 * sqrt(prod(diag(Sigma[[i]][1:2, 1:2])))
}
init.para = list(w=rep(1/3, 3), Beta=Beta+0.5*runif(1), DD=DD, Sigma=Sigma, Phi=Phi, ga=rep(1,3), nu=rep(30, 3))
fitN3 = try(FMMLMM.EM(Data, X, Z, g=3, init.para, cor.type = "UNC", true.cls=cls, tol = 1e-3, max.iter=500), silent=F)
        if(class(fitN3) == "try-error") next;      
fitT3 = try(FMMTLMM.EM(Data, X, Z, g=3, init.para, cor.type = "UNC", true.cls=cls, tol = 1e-3, max.iter=500), silent=F)
        if(class(fitT3) == "try-error") next;      
init.para = list(w=rep(1/3, 3), Beta=Beta+0.5*runif(1), DD=DD, Sigma=Sigma, nu=rep(nu[1], 3), rho=rep(rho[1], 3), Phi=Phi, ga=rep(1, 3))
fitC3 = try(FMMCNLMM.EM(Data, X, Z, g=3, init.para, cor.type = "UNC", true.cls=cls, tol = 1e-3, max.iter=500), silent=F)
        if(class(fitC3) == "try-error") next;      

# model selection
Tab1 = round(rbind(
c(fitN1$model.inf$aic, fitN2$model.inf$aic, fitN3$model.inf$aic,
  fitT1$model.inf$aic, fitT2$model.inf$aic, fitT3$model.inf$aic,
  fitC1$model.inf$aic, fitC2$model.inf$aic, fitC3$model.inf$aic),
c(fitN1$model.inf$bic, fitN2$model.inf$bic, fitN3$model.inf$bic,
  fitT1$model.inf$bic, fitT2$model.inf$bic, fitT3$model.inf$bic,
  fitC1$model.inf$bic, fitC2$model.inf$bic, fitC3$model.inf$bic),
c(fitN1$model.inf$icl, fitN2$model.inf$icl, fitN3$model.inf$icl,
  fitT1$model.inf$icl, fitT2$model.inf$icl, fitT3$model.inf$icl,
  fitC1$model.inf$icl, fitC2$model.inf$icl, fitC3$model.inf$icl),
c(fitN1$model.inf$awe, fitN2$model.inf$awe, fitN3$model.inf$awe,
  fitT1$model.inf$awe, fitT2$model.inf$awe, fitT3$model.inf$awe,
  fitC1$model.inf$awe, fitC2$model.inf$awe, fitC3$model.inf$awe)), 3)

Tab2 = cbind(
c(matrix(apply(Tab1[, 1:3], 1, order), nrow=3)[1,]),
c(matrix(apply(Tab1[, 4:6], 1, order), nrow=3)[1,]),
c(matrix(apply(Tab1[, 7:9], 1, order), nrow=3)[1,]),
c(matrix(apply(Tab1[, c(2,5,8)], 1, order), nrow=3)[1,]),
c(matrix(apply(Tab1[, c(2,8)], 1, order), nrow=2)[1,]))

# classification
Tab3 = rbind(
c(fitN2$pre.cls$CCR, fitT2$pre.cls$CCR, fitC2$pre.cls$CCR),
c(fitN2$pre.cls$ARI, fitT2$pre.cls$ARI, fitC2$pre.cls$ARI))

m = length(true.para) - 2*g
biasC = c(fitC2$para.est$para - true.para)

# explore the numerical results
rname = c('AIC', 'BIC', 'ICL', 'AWE')
write.table(cbind(Rep, Tab1), paste(WD.PATH, 'Simulation/result/', pd, '/fit.txt', sep=""), append=T, row.names = rname, col.names = F)
write.table(cbind(Rep, Tab2), paste(WD.PATH, 'Simulation/result/', pd, '/selection.txt', sep=""), append=T, row.names = rname, col.names = F)
write.table(cbind(Rep, Tab3), paste(WD.PATH, 'Simulation/result/', pd, '/class.txt', sep=""), append=T, row.names = c('CCR', 'ARI'), col.names = F)
write(c(Rep, biasC), paste(WD.PATH, 'Simulation/result/', pd, '/biasC.txt',sep=""), ncol=35, append=T)

Rep = Rep + 1
}
