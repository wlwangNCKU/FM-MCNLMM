################################################################################
#
#   Filename: fig3.r
#   Purpose: produce Figure 3 showing the numerical results for Simulation 1
#   Input data files: Simulation/result/SIM1/class.txt;
#                     Simulation/result/SIM2/class.txt; 
#                     Simulation/result/SIM3/class.txt; 
#                     Simulation/result/SIM4/class.txt; 
#   Output data files: Simulation/result/fig3.eps
#
################################################################################

# SIM1
Tab31 = read.table(paste(WD.PATH, 'Simulation/result/SIM1/class.txt', sep=""))
# SIM2
Tab32 = read.table(paste(WD.PATH, 'Simulation/result/SIM2/class.txt', sep=""))
# SIM3
Tab33 = read.table(paste(WD.PATH, 'Simulation/result/SIM3/class.txt', sep=""))
# SIM4
Tab34 = read.table(paste(WD.PATH, 'Simulation/result/SIM4/class.txt', sep=""))

name3 = c('criteria', 'Rep', 'FM-MLMM', 'FM-MtLMM', 'FM-MCNLMM')

colnames(Tab31) =  colnames(Tab32) = colnames(Tab33) = colnames(Tab34) = name3
setcol = c('gray10', 'blue', 'red4')
fullcol = c('gray100', 'cadetblue1', 'lightpink')

postscript(paste(WD.PATH, 'Simulation/result/fig3.eps',sep=""), width=12, height=3)
par(mfrow=c(1,4), mar=c(3,3,2,0.5), font.main=1)
boxplot(Tab31[which(Tab31 == 'CCR'), -c(1:2)], ylim=c(0.5,1), las=1, col=fullcol, border=setcol,
        main='n = 20', yaxt='n', pch=16, cex.main=2, cex.axis=0.9)
axis(2, seq(0.5, 1, 0.1), labels=c('0.50', '0.60', '0.70', '0.80', '0.90', '1.00'), las=1)
boxplot(Tab32[which(Tab32 == 'CCR'), -c(1:2)], ylim=c(0.5,1), las=1, col=fullcol, border=setcol,
        main='n = 50', yaxt='n', pch=16, cex.main=2, cex.axis=0.9)
axis(2, seq(0.5, 1, 0.1), labels=c('0.50', '0.60', '0.70', '0.80', '0.90', '1.00'), las=1)
boxplot(Tab33[which(Tab33 == 'CCR'), -c(1:2)], ylim=c(0.5,1), las=1, col=fullcol, border=setcol,
        main='n = 100', yaxt='n', pch=16, cex.main=2, cex.axis=0.9)
axis(2, seq(0.5, 1, 0.1), labels=c('0.50', '0.60', '0.70', '0.80', '0.90', '1.00'), las=1)
boxplot(Tab34[which(Tab34 == 'CCR'), -c(1:2)], ylim=c(0.5,1), las=1, col=fullcol, border=setcol,
        main='n = 200', yaxt='n', pch=16, cex.main=2, cex.axis=0.9)
axis(2, seq(0.5, 1, 0.1), labels=c('0.50', '0.60', '0.70', '0.80', '0.90', '1.00'), las=1)
dev.off()
