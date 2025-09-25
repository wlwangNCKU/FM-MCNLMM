################################################################################
#
#   Filename: figB1.r
#   Purpose: produce Supplementary Figure B.1 showing the numerical results 
#                    for Simulation 1
#   Input data files: Simulation/result/SIM1/baisC.txt;
#                     Simulation/result/SIM2/biasC.txt; 
#                     Simulation/result/SIM3/biasC.txt; 
#                     Simulation/result/SIM4/biasC.txt; 
#   Output data files: Simulation/result/figB1.eps;
#
################################################################################

# SIM1
biasC1 = read.table(paste(WD.PATH, 'Simulation/result/SIM1/biasC.txt', sep=""))
# SIM2
biasC2 = read.table(paste(WD.PATH, 'Simulation/result/SIM2/biasC.txt', sep=""))
# SIM3
biasC3 = read.table(paste(WD.PATH, 'Simulation/result/SIM3/biasC.txt', sep=""))
# SIM4
biasC4 = read.table(paste(WD.PATH, 'Simulation/result/SIM4/biasC.txt', sep=""))

idx = c(2, 4:9, 16:18, 22:24, 32, 34,
        3, 10:15, 19:21, 25:27, 33, 35)
bC1 = biasC1[, idx]
bC2 = biasC2[, idx]
bC3 = biasC3[, idx]
bC4 = biasC4[, idx]

nameC = c('w1', 'be11', 'be12', 'be13', 'be14', 'be15', 'be16', 'd111', 'd121', 'd122', 'sig111', 'sig121', 'sig122', 'nu1', 'rho1',
          'w2', 'be21', 'be22', 'be23', 'be24', 'be25', 'be26', 'd211', 'd221', 'd222', 'sig211', 'sig221', 'sig222', 'nu2', 'rho2')
colnames(bC1) = colnames(bC2) = colnames(bC3) = colnames(bC4) = nameC
par.name = c(expression(w[1]), expression(beta[11]), expression(beta[12]), expression(beta[13]), 
             expression(beta[14]), expression(beta[15]), expression(beta[16]), 
             expression(d[111]), expression(d[121]), expression(d[122]), expression(sigma[111]), expression(sigma[121]), 
             expression(sigma[122]), expression(nu[1]), expression(rho[1]),
             expression(w[2]), expression(beta[21]), expression(beta[22]), expression(beta[23]), 
             expression(beta[24]), expression(beta[25]), expression(beta[26]),
             expression(d[211]), expression(d[221]), expression(d[222]), expression(sigma[211]), expression(sigma[221]), 
             expression(sigma[222]), expression(nu[2]), expression(rho[2]))
m = length(nameC)

postscript(paste(WD.PATH, 'Simulation/result/figB1.eps',sep=""), width=7, height=12)
bd1 = c(-0.5, -5, -0.2, -3, -3, -0.4, -3, -2, -1, -2, -1, -1, -2, -0.5, -0.2,
        -0.3, -2, -0.2, -2, -2, -0.2, -1.5, -2, -2.5, -2, -1.2, -1, -2.5, -0.3, -0.15)
bd2 = c(0.3, 5, 0.2, 3, 3, 0.4, 3, 2, 1, 2, 2, 2, 4, 1, 0.8,
        0.5, 3, 0.6, 2, 2, 0.6, 1.5, 8.5, 6.5, 6, 1, 1.2, 2, 0.4, 0.3)
        
par(mfrow=c(6, 5), mar=c(2,2.5,2,0.5))
for(i in 1: m){
  All = cbind(bC1[,i], bC2[,i], bC3[,i], bC4[,i])
  boxplot(All, axes=F, pch=16, lty=1, col=0, main=par.name[i], ylim=c(bd1[i], bd2[i]))
  abline(h=0, lty=1, col='red4', lwd=1.5) 
  boxplot(All, pch=16, lty=1, col=0, names=c(20, 50, 100, 200), ylim=c(bd1[i], bd2[i]), ylab='', las=1, cex=0.25, cex.main=1.2, cex.axis=0.75, add=T)
}
dev.off()
