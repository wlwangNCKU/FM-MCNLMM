################################################################################
#
#   Filename: figB4.r
#   Purpose: produce Figure B.4 (fitted mean curves under the best model)
#            presented in Supplementary Appendix B.
#   Input data files: Application/result/fitADNI.RData
#   Output data files: Application/result/figB4.eps
#
################################################################################

load(paste(WD.PATH, 'Application/result/fitADNI.RData', sep=''))

# Estimates for the EFM-MCNLMM
est = efit32$para
pre.cls = efit32$pre.cls$post.cls
g1 = which(pre.cls == 1); n1 = length(g1)
g2 = which(pre.cls == 2); n2 = length(g2)
g3 = which(pre.cls == 3); n3 = length(g3)

# Create design matrices for each group
G1 = G2 = G3 = NULL
for(j in 1: n1) G1 = rbind(G1, Data[which(Data$Subject == g1[j]), ]) 
for(j in 1: n2) G2 = rbind(G2, Data[which(Data$Subject == g2[j]), ]) 
for(j in 1: n3) G3 = rbind(G3, Data[which(Data$Subject == g3[j]), ]) 
G10 = G1[which(G1$Sex == 0), ]; G11 = G1[which(G1$Sex == 1), ]
G20 = G2[which(G2$Sex == 0), ]; G21 = G2[which(G2$Sex == 1), ]
G30 = G3[which(G3$Sex == 0), ]; G31 = G3[which(G3$Sex == 1), ]

subj10 = unique(G10$Subject); n10 = length(subj10); subj11 = unique(G11$Subject); n11 = length(subj11)
subj20 = unique(G20$Subject); n20 = length(subj20); subj21 = unique(G21$Subject); n21 = length(subj21)
subj30 = unique(G30$Subject); n30 = length(subj30); subj31 = unique(G31$Subject); n31 = length(subj31)

age10 = edu10 = numeric(n10)
for(j in 1: n10){
 age10[j] = G10[which(G10$Subject == subj10[j]), 'AGE'][1]
 edu10[j] = G10[which(G10$Subject == subj10[j]), 'PTEDUCAT'][1]
}
age11 = edu11 = numeric(n11)
for(j in 1: n11){
 age11[j] = G11[which(G11$Subject == subj11[j]), 'AGE'][1]
 edu11[j] = G11[which(G11$Subject == subj11[j]), 'PTEDUCAT'][1]
}

age20 = edu20 = numeric(n20)
for(j in 1: n20){
 age20[j] = G20[which(G20$Subject == subj20[j]), 'AGE'][1]
 edu20[j] = G20[which(G20$Subject == subj20[j]), 'PTEDUCAT'][1]
}
age21 = edu21 = numeric(n21)
for(j in 1: n21){
 age21[j] = G21[which(G21$Subject == subj21[j]), 'AGE'][1]
 edu21[j] = G21[which(G21$Subject == subj21[j]), 'PTEDUCAT'][1]
}

age30 = edu30 = numeric(n30)
for(j in 1: n30){
 age30[j] = G30[which(G30$Subject == subj30[j]), 'AGE'][1]
 edu30[j] = G30[which(G30$Subject == subj30[j]), 'PTEDUCAT'][1]
}
age31 = edu31 = numeric(n31)
for(j in 1: n31){
 age31[j] = G31[which(G31$Subject == subj31[j]), 'AGE'][1]
 edu31[j] = G31[which(G31$Subject == subj31[j]), 'PTEDUCAT'][1]
}

r = 2
Ti = seq(0, 15, 0.1)
# Group 1, Female
X10 = kronecker(diag(r), cbind(1, Ti, Ti^2, mean(age10), mean(edu10), 0))  
b10 = colMeans(efit32$b$b1[subj10,])
# Group 1, Male
X11 = kronecker(diag(r), cbind(1, Ti, Ti^2, mean(age11), mean(edu11), 1))
b11 = colMeans(efit32$b$b1[subj11,])

# Group 2, Female
X20 = kronecker(diag(r), cbind(1, Ti, Ti^2, mean(age20), mean(edu20), 0))  
b20 = colMeans(efit32$b$b1[subj20,])
# Group 2, Male
X21 = kronecker(diag(r), cbind(1, Ti, Ti^2, mean(age21), mean(edu21), 1))
b21 = colMeans(efit32$b$b1[subj21,])

# Group 3, Female
X30 = kronecker(diag(r), cbind(1, Ti, Ti^2, mean(age30), mean(edu30), 0))  
b30 = colMeans(efit32$b$b1[subj30,])
# Group 3, Male
X31 = kronecker(diag(r), cbind(1, Ti, Ti^2, mean(age31), mean(edu31), 1))
b31 = colMeans(efit32$b$b1[subj31,])
Z0 = X10[, c(1,7)]

# Fitted responses for three predicted groups 
y10 = X10 %*% est$Beta[,1] + Z0 %*% b10
y11 = X11 %*% est$Beta[,1] + Z0 %*% b11

y20 = X20 %*% est$Beta[,2] + Z0 %*% b20
y21 = X21 %*% est$Beta[,2] + Z0 %*% b21

y30 = X30 %*% est$Beta[,3] + Z0 %*% b30
y31 = X31 %*% est$Beta[,3] + Z0 %*% b31

WD.PATH = paste(getwd(),"/FM-MCNLMM/",sep="")
postscript(paste(WD.PATH, 'Application/result/figB4.eps', sep=''), width=25, height=6)
par(mfrow=c(1, 2), mar=c(4,4.5,0.5, 0.5), cex.lab=1.2)
idx = 1: length(Ti)
sel = seq(1, 150, 15)
plot(Ti, y10[idx], ylim = c(min(Data$sqADAS13), max(Data$sqADAS13)), type='n', xlab='', ylab=expression(ADAS13^0.5), xaxt='n', las=1)
abline(h=seq(0, 9, 1), lty=1, col='gray', lwd=0.2)
abline(v=seq(0, 15, 1), lty=1, col='gray', lwd=0.2)
axis(1, seq(0,14,2), seq(0, 180, 24), cex.lab=0.5)
mtext('Time (Months)', 1, line=3, cex=1.2, font=1)
lines(Ti, y10[idx], lty=1, col='blue', lwd=1.5)
points(Ti[sel], y10[sel], pch=1, cex=1.2, col='blue')
lines(Ti, y20[idx], lty=2, col='green4', lwd=2)
points(Ti[sel], y20[sel], pch=1, cex=1.2, col='green4')
lines(Ti, y30[idx], lty=3, col='red3', lwd=3)
points(Ti[sel], y30[sel], pch=1, cex=1.2, col='red4')

lines(Ti, y11[idx], lty=1, col='blue', lwd=1.5)
points(Ti[sel], y11[sel], pch=16, cex=1.2, col='blue')
lines(Ti, y21[idx], lty=2, col='green4', lwd=2)
points(Ti[sel], y21[sel], pch=16, cex=1.2, col='green4')
lines(Ti, y31[idx], lty=3, col='red3', lwd=2)
points(Ti[sel], y31[sel], pch=16, cex=1.2, col='red4')
legend('bottom', c('Female in Group 1', 'Female in Group 2', 'Female in Group 3', 'Male in Group 1', 'Male in Group 2', 'Male in Group 3'), 
       lty=rep(c(1,2,3), 2), col=rep(c('blue','green4','red3'), 2), pch=rep(c(1,16), each=3), cex=0.9, lwd=2, ncol=2, bty='n')

plot(Ti, y10[-idx], ylim = c(min(Data$log10MidTemp), max(Data$log10MidTemp)), type='n', xlab='', ylab=expression(log[10](MidTemp)), xaxt='n', las=1)
abline(h=seq(3.9, 4.4, 0.05), lty=1, col='gray', lwd=0.2)
abline(v=seq(0, 15, 1), lty=1, col='gray', lwd=0.2)
axis(1, seq(0,14,2), seq(0, 180, 24), cex.lab=0.5)
mtext('Time (Months)', 1, line=3, cex=1.2, font=1)
lines(Ti, y10[-idx], lty=1, col='blue', lwd=1.5)
points(Ti[sel], y10[-idx][sel], pch=1, cex=1.2, col='blue')
lines(Ti, y20[-idx], lty=2, col='green4', lwd=2)
points(Ti[sel], y20[-idx][sel], pch=1, cex=1.2, col='green4')
lines(Ti, y30[-idx], lty=3, col='red3', lwd=3)
points(Ti[sel], y30[-idx][sel], pch=1, cex=1.2, col='red4')

lines(Ti, y11[-idx], lty=1, col='blue', lwd=1.5)
points(Ti[sel], y11[-idx][sel], pch=16, cex=1.2, col='blue')
lines(Ti, y21[-idx], lty=2, col='green4', lwd=2)
points(Ti[sel], y21[-idx][sel], pch=16, cex=1.2, col='green4')
lines(Ti, y31[-idx], lty=3, col='red3', lwd=2)
points(Ti[sel], y31[-idx][sel], pch=16, cex=1.2, col='red4')
legend('top', c('Female in Group 1', 'Female in Group 2', 'Female in Group 3', 'Male in Group 1', 'Male in Group 2', 'Male in Group 3'), 
       lty=rep(c(1,2,3), 2), col=rep(c('blue','green4','red3'), 2), pch=rep(c(1,16), each=3), lwd=2, cex=0.9, ncol=2, bty='n')
dev.off()
