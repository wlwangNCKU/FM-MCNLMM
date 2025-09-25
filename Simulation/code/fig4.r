################################################################################
#
#   Filename: fig4.r
#   Purpose: produce Figure 4 showing the numerical results for Simulation 2
#   Input data files: Simulation/result/ESIM1/class.txt;
#                     Simulation/result/ESIM2/class.txt; 
#                     Simulation/result/ESIM3/class.txt; 
#                     Simulation/result/ESIM4/class.txt; 
#                     Simulation/result/ESIM5/class.txt; 
#                     Simulation/result/ESIM6/class.txt 
#                     Simulation/result/ESIM7/class.txt; 
#                     Simulation/result/ESIM8/class.txt 
#   Output data files: Simulation/result/fig4.eps
#
################################################################################

Tab31 = read.table(paste(WD.PATH, 'Simulation/result/ESIM1/class.txt', sep=""))
Tab32 = read.table(paste(WD.PATH, 'Simulation/result/ESIM2/class.txt', sep=""))
Tab33 = read.table(paste(WD.PATH, 'Simulation/result/ESIM3/class.txt', sep=""))
Tab34 = read.table(paste(WD.PATH, 'Simulation/result/ESIM4/class.txt', sep=""))
Tab35 = read.table(paste(WD.PATH, 'Simulation/result/ESIM5/class.txt', sep=""))
Tab36 = read.table(paste(WD.PATH, 'Simulation/result/ESIM6/class.txt', sep=""))
Tab37 = read.table(paste(WD.PATH, 'Simulation/result/ESIM7/class.txt', sep=""))
Tab38 = read.table(paste(WD.PATH, 'Simulation/result/ESIM8/class.txt', sep=""))

colnames(Tab31) = colnames(Tab32) = colnames(Tab33) = colnames(Tab34) =
colnames(Tab35) = colnames(Tab36) = colnames(Tab37) = colnames(Tab38) =
 c('criteria', 'Rep', 'EMN', 'EMT', 'EMCN', 'MN', 'MT', 'MCN')

library(vioplot)

Size = c(30, 75, 150, 300)
CCR1 = rbind(
cbind(Size=Size[1], Tab31[which(Tab31$criteria == 'CCR'), c(5,8)]),
cbind(Size=Size[2], Tab32[which(Tab32$criteria == 'CCR'), c(5,8)]),
cbind(Size=Size[3], Tab33[which(Tab33$criteria == 'CCR'), c(5,8)]),
cbind(Size=Size[4], Tab34[which(Tab34$criteria == 'CCR'), c(5,8)]))

CCR0 = rbind(
cbind(Size=Size[1], Tab35[which(Tab35$criteria == 'CCR'), c(5,8)]),
cbind(Size=Size[2], Tab36[which(Tab36$criteria == 'CCR'), c(5,8)]),
cbind(Size=Size[3], Tab37[which(Tab37$criteria == 'CCR'), c(5,8)]),
cbind(Size=Size[4], Tab38[which(Tab38$criteria == 'CCR'), c(5,8)]))

postscript(paste(WD.PATH, 'Simulation/result/fig4.eps', sep=''), width=10, height=5)
par(mfrow=c(1,2), mar=c(4,4,3,0.5))
ylim1 = c(0.8, 1)
vioplot(MCN~Size, data=CCR1, ylim=ylim1, col = "lightyellow", rectCol="green4", lineCol="green3", border="green3",
        plotCentre = "line", side = "left", xlab='Sample Size', ylab='CCR', las=1, xaxt='n',
        main=expression(paste('(a) covariate-dependent data: ', psi, '=(6, -1, 8, -1.5)', sep='')))
abline(v=1:4, lty=2, lwd=0.2, col='gray10')
abline(h=seq(0.8, 1, 0.05), lty=2, lwd=0.2, col='gray10')
vioplot(MCN~Size, data=CCR1, ylim=ylim1, col = "lightyellow", rectCol="green4", lineCol="green3", border="green3", plotCentre = "line", side = "left", las=1, add=T)
vioplot(EMCN~Size, data=CCR1, ylim=ylim1, col = "lightpink", rectCol="red", lineCol="red4", border="red4", plotCentre = "line", side = "right", las=1, add=T)
axis(1, 1:4, labels=paste('n=', Size, sep=''), cex.axis=1.2)
legend("bottomright", c('EFM-MCNLMM','  FM-MCNLMM'), fill=c('lightpink','lightyellow'), border=c('red4','green3'), cex=1, bty='n')

ylim1 = c(0.8, 1)
vioplot(MCN~Size, data=CCR0, ylim=ylim1, col = "lightyellow", rectCol="green4", lineCol="green3", border="green3",
        plotCentre = "line", side = "left", xlab='Sample Size', ylab='CCR', las=1, xaxt='n',
        main=expression(paste('(b) covariate-independent data: ', psi, '=(0, 0, 0, 0)', sep='')))
abline(v=1:4, lty=2, lwd=0.2, col='gray10')
abline(h=seq(0.8, 1, 0.05), lty=2, lwd=0.2, col='gray10')
vioplot(MCN~Size, data=CCR0, ylim=ylim1, col = "lightyellow", rectCol="green4", lineCol="green3", border="green3", plotCentre = "line", side = "left", las=1, add=T)
vioplot(EMCN~Size, data=CCR0, ylim=ylim1, col = "lightpink", rectCol="red", lineCol="red4", border="red4", plotCentre = "line", side = "right", las=1, add=T)
axis(1, 1:4, labels=paste('n=', Size, sep=''), cex.axis=1.2)
legend("bottomright", c('EFM-MCNLMM','  FM-MCNLMM'), fill=c('lightpink','lightyellow'), border=c('red4','green3'), cex=1, bty='n')
dev.off()
