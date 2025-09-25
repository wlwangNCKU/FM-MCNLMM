################################################################################
#
#   Filename: adni.run.r
#   Purpose: perform model fitting of the subset of ADNI data
#   Input data files: Application/data/ADNIMERGE.csv; 
#                     Application/function/basicfn.r;
#                     Application/function/fmmlmm.fn.r;
#                     Application/function/fmmcnlmm.fn.r;
#                     Application/function/efmmlmm.fn.r;
#                     Application/function/efmmcnlmm.fn.r
#   Output data files: Application/result/fitADNI.RData
#   Required R packages: mvtnorm; nlme; car 
#
################################################################################

library(nlme)
library(car)

source(paste(WD.PATH, 'Application/function/basicfn.r', sep=''))
source(paste(WD.PATH, 'Application/function/fmmlmm.fn.r', sep=''))
source(paste(WD.PATH, 'Application/function/fmmcnlmm.fn.r', sep=''))
source(paste(WD.PATH, 'Application/function/efmmlmm.fn.r', sep=''))
source(paste(WD.PATH, 'Application/function/efmmcnlmm.fn.r', sep=''))

###### import data ######
ADNI = read.csv(paste(WD.PATH, 'Application/data/ADNIMERGE.csv', sep=''), header = T)
Subj = sort(ADNI[which(ADNI$COLPROT=='ADNI1' & ADNI$VISCODE=='bl' & ADNI$PTRACCAT=='White' & ADNI$PTMARRY == 'Married'), 'RID'])
N = length(Subj)

# Preliminary
Subj1 = sort(unique(ADNI[which(ADNI$COLPROT=='ADNI1'), 'RID'])); N1 = length(Subj1); N1
Subj2 = sort(unique(ADNI[which(ADNI$COLPROT=='ADNIGO'), 'RID'])); N2 = length(Subj2); N2
Subj3 = sort(unique(ADNI[which(ADNI$COLPROT=='ADNI2'), 'RID'])); N3 = length(Subj3); N3
Subj4 = sort(unique(ADNI[which(ADNI$COLPROT=='ADNI3'), 'RID'])); N4 = length(Subj4); N4

###### Reorganized Data ###### 
ti = sort(unique(ADNI$Month))
len.ti = length(ti)
ADNI1 = NULL
for(i in 1: N)for(k in 1: len.ti){
 ind = which(ADNI$RID == Subj[i] & ADNI$Month == ti[k])
 if(length(ind) == 0) next;
 IDi = ADNI[ind, ][1, ]
 ADNI1 = rbind(ADNI1, cbind(Subject=i, IDi[order(IDi$M), ]))
}
dim(ADNI1)
length(unique(ADNI1$Subject))

ADNI1$Time = ADNI1$Month/12
ADNI1$Time2 = ADNI1$Time^2
ADNI1$Sex = as.numeric(ADNI1$PTGENDER=='Male')
ADNI1$sqADAS13 = sqrt(ADNI1$ADAS13)
ADNI1$sqMMSE = sqrt(ADNI1$MMSE)
ADNI1$log10MidTemp = log10(ADNI1$MidTemp)
ADNI1$IDX = recode(ADNI1$DX, " c('','CN')=1; 'MCI'=2; 'Dementia'=3 ")
ADNI1$DXb1 = recode(ADNI1$DX_bl, " 'CN'=1; 'LMCI'=2; 'AD'=3 ")
varNAME = c('RID', 'Month', 'Time', 'Time2', 'ADAS13', 'MidTemp', 'sqADAS13', 'log10MidTemp', 'DX', 'IDX', 'DX_bl', 'DXb1', 'AGE', 'PTEDUCAT', 'PTGENDER', 'Sex', 'PTMARRY', 'PTRACCAT')

ADNI10 = ADNI1[, varNAME]
dim(ADNI10)
ADNI10.narm = na.omit(ADNI10)
dim(ADNI10.narm)
colSums(apply(ADNI10, 2, is.na))
length(unique(ADNI10$RID))
length(unique(ADNI10.narm$RID))

selID = as.numeric(names(table(ADNI10.narm$RID))[which(table(ADNI10.narm$RID)>1)])
n = length(selID)
subADNI = NULL
for(j in 1: n) subADNI = rbind(subADNI, cbind(Subject=j, ADNI10.narm[which(ADNI10.narm$RID == selID[j]), ]))

r = 2
Data = groupedData(sqADAS13+log10MidTemp ~ Time|Subject, data=subADNI)
dim(Data)
cluster = as.numeric(levels(Data$Subject))[Data$Subject]
n = length(unique(cluster))

# Preliminary
sj= numeric(n)
for(j in 1: n) sj[j] = length(unique(Data$Time[Data$Subject == j]))
sum(sj)
cls = numeric(n)
for(j in 1: n) cls[j] = Data$IDX[Data$Subject==j][sj[j]]
print(table(cls))

# Response & Design matrix
y = NULL
for(j in 1: n) y = c(y, Data$sqADAS13[Data$Subject == j], Data$log10MidTemp[Data$Subject == j])

X = kronecker(diag(r), cbind(1, Data$Time, Data$Time2, Data$AGE, Data$PTEDUCAT, Data$Sex)[which(Data$Subject == 1), ])
for(j in 2: n) X = rbind(X, kronecker(diag(r), cbind(1, Data$Time, Data$Time2, Data$AGE, Data$PTEDUCAT, Data$Sex)[which(Data$Subject == j), ]))

V = NULL
for(j in 1: n) V = rbind(V, c(1, as.numeric(Data[Data$Subject==j, c('DXb1','Sex')][1,])))

# Initialization
g = 3
g1Data = Data[which(Data$IDX==1), ]
g2Data = Data[which(Data$IDX==2), ]
g3Data = Data[which(Data$IDX==3), ]

# Group 1 
fm11 = lme(sqADAS13 ~ Time + Time2 + AGE + PTEDUCAT + Sex, data=g1Data, random=~1, corAR1(), method='ML')
print(summary(fm11))
fm12 = lme(log10MidTemp ~ Time + Time2 + AGE + PTEDUCAT + Sex, data=g1Data, random=~1, corAR1(), method='ML')
print(summary(fm12))

# Group 2 
fm21 = lme(sqADAS13 ~ Time + Time2 + AGE + PTEDUCAT + Sex, data=g2Data, random=~1, corAR1(), method='ML')
print(summary(fm21))
fm22 = lme(log10MidTemp ~ Time + Time2 + AGE + PTEDUCAT + Sex, data=g2Data, random=~1, corAR1(), method='ML')
print(summary(fm22))

# Group 3 
fm31 = lme(sqADAS13 ~ Time + Time2 + AGE + PTEDUCAT + Sex, data=g3Data, random=~1, corAR1(), method='ML')
print(summary(fm31)) 
fm32 = lme(log10MidTemp ~ Time + Time2 + AGE + PTEDUCAT + Sex, data=g3Data, random=~1, corAR1(), method='ML')
print(summary(fm32))

Beta = matrix(0, nrow=ncol(X), ncol=3)
Beta[,1] = c(fm11$coefficients$fix, fm12$coefficients$fix)
Beta[,2] = c(fm21$coefficients$fix, fm22$coefficients$fix)
Beta[,3] = c(fm31$coefficients$fix, fm32$coefficients$fix)

Sigma = as.list(g)
Sigma[[1]] = cor(cbind(g1Data$sqADAS13, g1Data$log10MidTemp))
Sigma[[2]] = cor(cbind(g2Data$sqADAS13, g2Data$log10MidTemp))
Sigma[[3]] = cor(cbind(g3Data$sqADAS13, g3Data$log10MidTemp))
Phi = rep(1e-6, g)
ga = rep(1, g)
psi = rep(0, (g-1)*ncol(V))
nu = rep(0.5, g)
rho = rep(0.5, g)

Z = X[, seq(1,12,6)]
q = ncol(Z)
DD = array(NA, dim=c(q,q,g))
for(i in 1: g) DD[,,i] = 0.5*diag(q)

###### Fit Models ######
init.para = list(w = rep(1/g, g), Beta=Beta, DD=DD, Sigma=Sigma, Phi=Phi, ga=ga)
fit10 = FMMLMM.EM(Data, y, X, Z, g=g, init.para, cor.type = "UNC", true.cls=cls, tol = 1e-3, max.iter=500)   
fit11 = FMMLMM.EM(Data, y, X, Z, g=g, init.para, cor.type = "CAR1", true.cls=cls, tol = 1e-3, max.iter=500)
fit12 = FMMLMM.EM(Data, y, X, Z, g=g, init.para, cor.type = "DEC", true.cls=cls, tol = 1e-3, max.iter=500)  
save.image(paste(WD.PATH, 'Application/result/fitADNI.RData', sep=''))
             
init.para = list(w = fit10$para$w, Beta=fit10$para$Beta, DD=fit10$para$DD, Sigma=fit10$para$Sigma, Phi=Phi, ga=ga, nu=nu, rho=rho)
fit30 = FMMCNLMM.EM(Data, y, X, Z, g=g, init.para, cor.type = "UNC", true.cls=cls, tol = 1e-3, max.iter=500)   
fit31 = FMMCNLMM.EM(Data, y, X, Z, g=g, init.para, cor.type = "CAR1", true.cls=cls, tol = 1e-3, max.iter=500)
fit32 = FMMCNLMM.EM(Data, y, X, Z, g=g, init.para, cor.type = "DEC", true.cls=cls, tol = 1e-3, max.iter=500)  
save.image(paste(WD.PATH, 'Application/result/fitADNI.RData', sep=''))

###### Fit Extended Models ######
init.para = list(Beta=Beta, DD=DD, Sigma=Sigma, Phi=Phi, ga=ga, psi=psi)
efit10 = EFMMLMM.EM(Data, y, X, Z, V, g=g, init.para, cor.type = "UNC", true.cls=cls, tol = 1e-3, max.iter=500)   
efit11 = EFMMLMM.EM(Data, y, X, Z, V, g=g, init.para, cor.type = "CAR1", true.cls=cls, tol = 1e-3, max.iter=500)
efit12 = EFMMLMM.EM(Data, y, X, Z, V, g=g, init.para, cor.type = "DEC", true.cls=cls, tol = 1e-3, max.iter=500)  
save.image(paste(WD.PATH, 'Application/result/fitADNI.RData', sep=''))

init.para = list(Beta=efit10$para$Beta, DD=efit10$para$DD, Sigma=efit10$para$Sigma, Phi=Phi, ga=ga, nu=nu, rho=rho, psi=psi)
efit30 = EFMMCNLMM.EM(Data, y, X, Z, V, g=g, init.para, cor.type = "UNC", true.cls=cls, tol = 1e-3, max.iter=500)   
efit31 = EFMMCNLMM.EM(Data, y, X, Z, V, g=g, init.para, cor.type = "CAR1", true.cls=cls, tol = 1e-3, max.iter=500)
efit32 = EFMMCNLMM.EM(Data, y, X, Z, V, g=g, init.para, cor.type = "DEC", true.cls=cls, tol = 1e-3, max.iter=500)  
save.image(paste(WD.PATH, 'Application/result/fitADNI.RData', sep=''))
