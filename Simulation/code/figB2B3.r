################################################################################
#
#   Filename: figB2B3.r
#   Purpose: produce Supplementary Figures B.2 and B.3 showing the numerical 
#                    results for Simulation 2
#   Input data files: Simulation/result/ESIM1/estEN.txt; ./estET.txt; ./estEC.txt; 
#                     Simulation/result/ESIM2/estEN.txt; ./estET.txt; ./estEC.txt; 
#                     Simulation/result/ESIM3/estEN.txt; ./estET.txt; ./estEC.txt; 
#                     Simulation/result/ESIM4/estEN.txt; ./estET.txt; ./estEC.txt; 
#                     Simulation/result/ESIM5/estEN.txt; ./estET.txt; ./estEC.txt; 
#                     Simulation/result/ESIM6/estEN.txt; ./estET.txt; ./estEC.txt; 
#                     Simulation/result/ESIM7/estEN.txt; ./estET.txt; ./estEC.txt; 
#                     Simulation/result/ESIM8/estEN.txt; ./estET.txt; ./estEC.txt; 
#   Output data files: Simulation/result/figB2a1.eps; ./figB2a2.eps; ./figB2a3.eps;
#                      Simulation/result/figB3a1.eps; ./figB3a2.eps; ./figB3a3.eps;
#
################################################################################

vech.posi=function(dim) cbind(rep(1:dim, 1:dim), unlist(mapply(':', 1, 1:dim)))

# ESIM1
estEC1 = read.table(paste(WD.PATH, 'Simulation/result/ESIM1/estEC.txt', sep=""))[,1:41]
estEN1 = read.table(paste(WD.PATH, 'Simulation/result/ESIM1/estEN.txt', sep=""))[,1:35]
estET1 = read.table(paste(WD.PATH, 'Simulation/result/ESIM1/estET.txt', sep=""))[,1:38]

# ESIM2
estEC2 = read.table(paste(WD.PATH, 'Simulation/result/ESIM2/estEC.txt', sep=""))[,1:41]
estEN2 = read.table(paste(WD.PATH, 'Simulation/result/ESIM2/estEN.txt', sep=""))[,1:35]
estET2 = read.table(paste(WD.PATH, 'Simulation/result/ESIM2/estET.txt', sep=""))[,1:38]

# ESIM3
estEC3 = read.table(paste(WD.PATH, 'Simulation/result/ESIM3/estEC.txt', sep=""))[,1:41]
estEN3 = read.table(paste(WD.PATH, 'Simulation/result/ESIM3/estEN.txt', sep=""))[,1:35]
estET3 = read.table(paste(WD.PATH, 'Simulation/result/ESIM3/estET.txt', sep=""))[,1:38]

# ESIM4
estEC4 = read.table(paste(WD.PATH, 'Simulation/result/ESIM4/estEC.txt', sep=""))[,1:41]
estEN4 = read.table(paste(WD.PATH, 'Simulation/result/ESIM4/estEN.txt', sep=""))[,1:35]
estET4 = read.table(paste(WD.PATH, 'Simulation/result/ESIM4/estET.txt', sep=""))[,1:38]

# ESIM5
estEC5 = read.table(paste(WD.PATH, 'Simulation/result/ESIM5/estEC.txt', sep=""))[,1:41]
estEN5 = read.table(paste(WD.PATH, 'Simulation/result/ESIM5/estEN.txt', sep=""))[,1:35]
estET5 = read.table(paste(WD.PATH, 'Simulation/result/ESIM5/estET.txt', sep=""))[,1:38]

# ESIM6
estEC6 = read.table(paste(WD.PATH, 'Simulation/result/ESIM6/estEC.txt', sep=""))[,1:41]
estEN6 = read.table(paste(WD.PATH, 'Simulation/result/ESIM6/estEN.txt', sep=""))[,1:35]
estET6 = read.table(paste(WD.PATH, 'Simulation/result/ESIM6/estET.txt', sep=""))[,1:38]

# ESIM7
estEC7 = read.table(paste(WD.PATH, 'Simulation/result/ESIM7/estEC.txt', sep=""))[,1:41]
estEN7 = read.table(paste(WD.PATH, 'Simulation/result/ESIM7/estEN.txt', sep=""))[,1:35]
estET7 = read.table(paste(WD.PATH, 'Simulation/result/ESIM7/estET.txt', sep=""))[,1:38]

# ESIM8
estEC8 = read.table(paste(WD.PATH, 'Simulation/result/ESIM8/estEC.txt', sep=""))[,1:41]
estEN8 = read.table(paste(WD.PATH, 'Simulation/result/ESIM8/estEN.txt', sep=""))[,1:35]
estET8 = read.table(paste(WD.PATH, 'Simulation/result/ESIM8/estET.txt', sep=""))[,1:38]

# basic setting
g=3; p=4; q=2; r=2; sj=7
vech.D = vech.posi(q)
vech.Sig = vech.posi(r)

Beta = matrix(c(1, 2, 2, -1, -1, 1, 0, -2, -2, 3, 1, -3), ncol=g)
DD = array(NA, dim=c(q,q,g))
for(i in 1: g)  DD[,,i] = matrix(c(2, 0.25, 0.25, 2), q, q, byrow=T)
Sigma = as.list(g)
for(i in 1: g){
Sigma[[i]] = i*diag(c(1,2))
Sigma[[i]][1, 2] = Sigma[[i]][2, 1] = 0.5 * sqrt(prod(diag(Sigma[[i]][1:2, 1:2])))
}
nu = rep(0.25, 3)
rho = rep(0.2, 3)
psi = c(6, -1, 8, -1.5)

true.par = NULL
for(i in 1: g) true.par = c(true.par, Beta[,i], DD[,,i][vech.D], Sigma[[i]][vech.Sig])
true.parE = c(true.par, psi)
m = length(true.parE)
ETHETA = matrix(rep(true.parE, each=100), nrow=100)
true.par0 = c(true.par, rep(0, 4))
E0THETA = matrix(rep(true.par0, each=100), nrow=100)

rmN = c(12, 13, 24, 25, 36, 37)
rmT = c(13, 25, 37)

nameEC = c('Rep', c('be11', 'be12', 'be13', 'be14', 'd111', 'd121', 'd122', 'sig111', 'sig121', 'sig122', 'nu1', 'rho1',
          'be21', 'be22', 'be23', 'be24', 'd211', 'd221', 'd222', 'sig211', 'sig221', 'sig222', 'nu2', 'rho2',
          'be31', 'be32', 'be33', 'be34', 'd311', 'd321', 'd322', 'sig311', 'sig321', 'sig322', 'nu3', 'rho3',
          'psi1', 'psi2', 'psi3', 'psi4'))

colnames(estEN1) = colnames(estEN2) = colnames(estEN3) = colnames(estEN4) = 
colnames(estEN5) = colnames(estEN6) = colnames(estEN7) = colnames(estEN8) = nameEC[-rmN]

colnames(estET1) = colnames(estET2) = colnames(estET3) = colnames(estET4) = 
colnames(estET5) = colnames(estET6) = colnames(estET7) = colnames(estET8) = nameEC[-rmT]

colnames(estEC1) = colnames(estEC2) = colnames(estEC3) = colnames(estEC4) = 
colnames(estEC5) = colnames(estEC6) = colnames(estEC7) = colnames(estEC8) = nameEC

# EFM-MLMM
ENbias1 = (estEN1[,-1]-ETHETA)
ENbias2 = (estEN2[,-1]-ETHETA)
ENbias3 = (estEN3[,-1]-ETHETA)
ENbias4 = (estEN4[,-1]-ETHETA)
ENbias5 = (estEN5[,-1]-E0THETA)
ENbias6 = (estEN6[,-1]-E0THETA)
ENbias7 = (estEN7[,-1]-E0THETA)
ENbias8 = (estEN8[,-1]-E0THETA)

# EFM-MTLMM
ETbias1 = (estET1[, -c(1,12,23,34)]-ETHETA)
ETbias2 = (estET2[, -c(1,12,23,34)]-ETHETA)
ETbias3 = (estET3[, -c(1,12,23,34)]-ETHETA)
ETbias4 = (estET4[, -c(1,12,23,34)]-ETHETA)
ETbias5 = (estET5[, -c(1,12,23,34)]-E0THETA)
ETbias6 = (estET6[, -c(1,12,23,34)]-E0THETA)
ETbias7 = (estET7[, -c(1,12,23,34)]-E0THETA)
ETbias8 = (estET8[, -c(1,12,23,34)]-E0THETA)

# EFM-MCNLMM
ECbias1 = (estEC1[, -c(1, rmN)]-ETHETA)
ECbias2 = (estEC2[, -c(1, rmN)]-ETHETA)
ECbias3 = (estEC3[, -c(1, rmN)]-ETHETA)
ECbias4 = (estEC4[, -c(1, rmN)]-ETHETA)
ECbias5 = (estEC5[, -c(1, rmN)]-E0THETA)
ECbias6 = (estEC6[, -c(1, rmN)]-E0THETA)
ECbias7 = (estEC7[, -c(1, rmN)]-E0THETA)
ECbias8 = (estEC8[, -c(1, rmN)]-E0THETA)

par.nameN = c(expression(beta[11]), expression(beta[12]), expression(beta[13]), expression(beta[14]), 
             expression(d[111]), expression(d[121]), expression(d[122]), expression(sigma[111]), expression(sigma[121]), expression(sigma[122]),
             expression(beta[21]), expression(beta[22]), expression(beta[23]), expression(beta[24]), 
             expression(d[211]), expression(d[221]), expression(d[222]), expression(sigma[211]), expression(sigma[221]), expression(sigma[222]), 
             expression(beta[31]), expression(beta[32]), expression(beta[33]), expression(beta[34]), 
             expression(d[311]), expression(d[321]), expression(d[322]), expression(sigma[311]), expression(sigma[321]), expression(sigma[322]), 
             expression(psi[1]), expression(psi[2]), expression(psi[3]), expression(psi[4]))

## Case 1 ## 
bd11 = c(-3, -0.3, -3, -0.4, -2, 
         -5, -2, -1, -1, -2, 
         -2, -0.2, -2, -0.4, -2, 
         -5, -2, -1, -1, -1,  
         -2.5, -0.4, -2.5, -0.5, -2,  
         -5, -2, -1.5, -2, -3,
         -5, -2.5, -5, -3)
bd12 = c(4, 0.4, 4, 1, 5, 
         5, 5, 2.5, 2, 5, 
         2.5, 0.4, 2.5, 0.5, 15, 
         10, 17, 6, 6, 12, 
         3, 0.3, 2.5, 0.5, 8,
         4, 10, 7, 6, 15, 
         12, 0.5, 17, 1)

#win.graph(width=20, height=45)
postscript(paste(WD.PATH, 'Simulation/result/FigB2a1.eps',sep=""), width=25, height=5.3)
par(mfrow=c(2, 5), mar=c(2,2.5,2.5,0.5))
for(i in 1: 10){
Size = c(30, 75, 150, 300)
BS1 = rbind(
cbind(Size=Size[1], EMN=ENbias1[,i], EMT=ETbias1[,i], EMCN=ECbias1[,i]),
cbind(Size=Size[2], EMN=ENbias2[,i], EMT=ETbias2[,i], EMCN=ECbias2[,i]),
cbind(Size=Size[3], EMN=ENbias3[,i], EMT=ETbias3[,i], EMCN=ECbias3[,i]),
cbind(Size=Size[4], EMN=ENbias4[,i], EMT=ETbias4[,i], EMCN=ECbias4[,i]))               
boxplot(EMN~Size, data=BS1, col = "white", border="black", at = 1:4-0.25, boxwex = 0.2, ylim=c(bd11[i], bd12[i]),
        xlab='Sample Size', ylab='Bias', las=1, xaxt='n', main=par.nameN[i], cex.axis=1.2, cex.main=1.75)
#mtext('Bias', 2, line=2, cex=0.5, font=1)
abline(h=0, lwd=1.5, col='gray50')
boxplot(EMN~Size, data=BS1, col = "white", border="black", at = 1:4-0.25, boxwex = 0.2, add = T, yaxt='n', xaxt='n')
boxplot(EMT~Size, data=BS1, col = "lightblue", border="blue3", at = 1:4, boxwex = 0.2, add = T, yaxt='n', xaxt='n')
boxplot(EMCN~Size, data=BS1, col = "lightpink", border="red4", at = 1:4+0.25, boxwex = 0.2, add=T, yaxt='n', xaxt='n')
axis(1, 1:4, labels=Size, cex.axis=1.2)
}
dev.off()

postscript(paste(WD.PATH, 'Simulation/result/FigB2a2.eps',sep=""), width=25, height=5.3)
par(mfrow=c(2, 5), mar=c(2,2.5,2.5,0.5))
for(i in 11: 20){
Size = c(30, 75, 150, 300)
BS1 = rbind(
cbind(Size=Size[1], EMN=ENbias1[,i], EMT=ETbias1[,i], EMCN=ECbias1[,i]),
cbind(Size=Size[2], EMN=ENbias2[,i], EMT=ETbias2[,i], EMCN=ECbias2[,i]),
cbind(Size=Size[3], EMN=ENbias3[,i], EMT=ETbias3[,i], EMCN=ECbias3[,i]),
cbind(Size=Size[4], EMN=ENbias4[,i], EMT=ETbias4[,i], EMCN=ECbias4[,i]))               
boxplot(EMN~Size, data=BS1, col = "white", border="black", at = 1:4-0.25, boxwex = 0.2, ylim=c(bd11[i], bd12[i]),
        xlab='Sample Size', ylab='Bias', las=1, xaxt='n', main=par.nameN[i], cex.axis=1.2, cex.main=1.75)
abline(h=0, lwd=1.5, col='gray50')
boxplot(EMN~Size, data=BS1, col = "white", border="black", at = 1:4-0.25, boxwex = 0.2, add = T, yaxt='n', xaxt='n')
boxplot(EMT~Size, data=BS1, col = "lightblue", border="blue3", at = 1:4, boxwex = 0.2, add = T, yaxt='n', xaxt='n')
boxplot(EMCN~Size, data=BS1, col = "lightpink", border="red4", at = 1:4+0.25, boxwex = 0.2, add=T, yaxt='n', xaxt='n')
axis(1, 1:4, labels=Size, cex.axis=1.2)
}
dev.off()

postscript(paste(WD.PATH, 'Simulation/result/FigB2a3.eps',sep=""), width=25, height=14)
par(mfrow=c(3, 5), mar=c(2,2.5,2.5,0.5))
for(i in 21: 34){
Size = c(30, 75, 150, 300)
BS1 = rbind(
cbind(Size=Size[1], EMN=ENbias1[,i], EMT=ETbias1[,i], EMCN=ECbias1[,i]),
cbind(Size=Size[2], EMN=ENbias2[,i], EMT=ETbias2[,i], EMCN=ECbias2[,i]),
cbind(Size=Size[3], EMN=ENbias3[,i], EMT=ETbias3[,i], EMCN=ECbias3[,i]),
cbind(Size=Size[4], EMN=ENbias4[,i], EMT=ETbias4[,i], EMCN=ECbias4[,i]))               
boxplot(EMN~Size, data=BS1, col = "white", border="black", at = 1:4-0.25, boxwex = 0.2, ylim=c(bd11[i], bd12[i]),
        xlab='Sample Size', ylab='Bias', las=1, xaxt='n', main=par.nameN[i], cex.axis=1.2, cex.main=1.75)    
abline(h=0, lwd=1.5, col='gray50')
boxplot(EMN~Size, data=BS1, col = "white", border="black", at = 1:4-0.25, boxwex = 0.2, add = T, yaxt='n', xaxt='n')
boxplot(EMT~Size, data=BS1, col = "lightblue", border="blue3", at = 1:4, boxwex = 0.2, add = T, yaxt='n', xaxt='n')
boxplot(EMCN~Size, data=BS1, col = "lightpink", border="red4", at = 1:4+0.25, boxwex = 0.2, add=T, yaxt='n', xaxt='n')
axis(1, 1:4, labels=Size, cex.axis=1.2)
}
plot(0:1, 0:1, type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
legend('center', c('EFM-MLMM', 'EFM-MtLMM', 'EFM-MCNLMM'), fill=c('white', 'lightblue', 'lightpink'), border=c('black', 'blue3', 'red4'), bty='n', cex=1.75)
dev.off()

## Case 2 ##
bd21 = c(-2, -0.3, -3, -0.4, -2, 
         -2, -2, -1, -1, -2, 
         -2, -0.25, -2, -0.3, -2, 
         -5.5, -2, -1, -1.5, -2,  
         -3, -0.35, -3, -0.6, -2,  
         -5, -2, -1.5, -3, -2.5,
         -4, -0.5, -1.5, -0.4)
bd22 = c(3, 0.3, 3, 0.5, 5, 
         3, 5, 2, 2, 3, 
         3, 1, 2.5, 0.8, 20, 
         15, 25, 6, 5.5, 12, 
         2.5, 0.35, 3, 0.5, 10,
         6, 10, 8, 6, 15, 
         2, 0.5, 2, 0.3)

postscript(paste(WD.PATH, 'Simulation/result/FigB3a1.eps',sep=""), width=25, height=5.3)
par(mfrow=c(2, 5), mar=c(2,2.5,2.5,0.5))
for(i in 1: 10){
Size = c(30, 75, 150, 300)
BS2 = rbind(
cbind(Size=Size[1], EMN=ENbias5[,i], EMT=ETbias5[,i], EMCN=ECbias5[,i]),
cbind(Size=Size[2], EMN=ENbias6[,i], EMT=ETbias6[,i], EMCN=ECbias6[,i]),
cbind(Size=Size[3], EMN=ENbias7[,i], EMT=ETbias7[,i], EMCN=ECbias7[,i]),
cbind(Size=Size[4], EMN=ENbias8[,i], EMT=ETbias8[,i], EMCN=ECbias8[,i]))               
boxplot(EMN~Size, data=BS2, col = "white", border="black", at = 1:4-0.25, boxwex = 0.2, ylim=c(bd21[i], bd22[i]),
        xlab='Sample Size', ylab='Bias', las=1, xaxt='n', main=par.nameN[i], cex.axis=1.2, cex.main=1.75)
abline(h=0, lwd=1.5, col='gray50')
boxplot(EMN~Size, data=BS2, col = "white", border="black", at = 1:4-0.25, boxwex = 0.2, add = T, yaxt='n', xaxt='n')
boxplot(EMT~Size, data=BS2, col = "lightblue", border="blue3", at = 1:4, boxwex = 0.2, add = T, yaxt='n', xaxt='n')
boxplot(EMCN~Size, data=BS2, col = "lightpink", border="red4", at = 1:4+0.25, boxwex = 0.2, add=T, yaxt='n', xaxt='n')
axis(1, 1:4, labels=Size, cex.axis=1.2)
}
dev.off()

postscript(paste(WD.PATH, 'Simulation/result/FigB3a2.eps',sep=""), width=25, height=5.3)
par(mfrow=c(2, 5), mar=c(2,2.5,2.5,0.5))
for(i in 11: 20){
Size = c(30, 75, 150, 300)
BS2 = rbind(
cbind(Size=Size[1], EMN=ENbias5[,i], EMT=ETbias5[,i], EMCN=ECbias5[,i]),
cbind(Size=Size[2], EMN=ENbias6[,i], EMT=ETbias6[,i], EMCN=ECbias6[,i]),
cbind(Size=Size[3], EMN=ENbias7[,i], EMT=ETbias7[,i], EMCN=ECbias7[,i]),
cbind(Size=Size[4], EMN=ENbias8[,i], EMT=ETbias8[,i], EMCN=ECbias8[,i]))               
boxplot(EMN~Size, data=BS2, col = "white", border="black", at = 1:4-0.25, boxwex = 0.2, ylim=c(bd21[i], bd22[i]),
        xlab='Sample Size', ylab='Bias', las=1, xaxt='n', main=par.nameN[i], cex.axis=1.2, cex.main=1.75)
abline(h=0, lwd=1.5, col='gray50')
boxplot(EMN~Size, data=BS2, col = "white", border="black", at = 1:4-0.25, boxwex = 0.2, add = T, yaxt='n', xaxt='n')
boxplot(EMT~Size, data=BS2, col = "lightblue", border="blue3", at = 1:4, boxwex = 0.2, add = T, yaxt='n', xaxt='n')
boxplot(EMCN~Size, data=BS2, col = "lightpink", border="red4", at = 1:4+0.25, boxwex = 0.2, add=T, yaxt='n', xaxt='n')
axis(1, 1:4, labels=Size, cex.axis=1.2)
}
dev.off()

postscript(paste(WD.PATH, 'Simulation/result/FigB3a3.eps',sep=""), width=25, height=14)
par(mfrow=c(3, 5), mar=c(2,2.5,2.5,0.5))
for(i in 21: 34){
Size = c(30, 75, 150, 300)
BS2 = rbind(
cbind(Size=Size[1], EMN=ENbias5[,i], EMT=ETbias5[,i], EMCN=ECbias5[,i]),
cbind(Size=Size[2], EMN=ENbias6[,i], EMT=ETbias6[,i], EMCN=ECbias6[,i]),
cbind(Size=Size[3], EMN=ENbias7[,i], EMT=ETbias7[,i], EMCN=ECbias7[,i]),
cbind(Size=Size[4], EMN=ENbias8[,i], EMT=ETbias8[,i], EMCN=ECbias8[,i]))               
boxplot(EMN~Size, data=BS2, col = "white", border="black", at = 1:4-0.25, boxwex = 0.2, ylim=c(bd21[i], bd22[i]),
        xlab='Sample Size', ylab='Bias', las=1, xaxt='n', main=par.nameN[i], cex.axis=1.2, cex.main=1.75)
abline(h=0, lwd=1.5, col='gray50')
boxplot(EMN~Size, data=BS2, col = "white", border="black", at = 1:4-0.25, boxwex = 0.2, add = T, yaxt='n', xaxt='n')
boxplot(EMT~Size, data=BS2, col = "lightblue", border="blue3", at = 1:4, boxwex = 0.2, add = T, yaxt='n', xaxt='n')
boxplot(EMCN~Size, data=BS2, col = "lightpink", border="red4", at = 1:4+0.25, boxwex = 0.2, add=T, yaxt='n', xaxt='n')
axis(1, 1:4, labels=Size, cex.axis=1.2)
}
plot(0:1, 0:1, type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
legend('center', c('EFM-MLMM', 'EFM-MtLMM', 'EFM-MCNLMM'), fill=c('white', 'lightblue', 'lightpink'), border=c('black', 'blue3', 'red4'), bty='n', cex=1.75)
dev.off()
 