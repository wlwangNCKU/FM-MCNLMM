################################################################################
#
#   Filename: fig6.r
#   Purpose: produce Figure 6 (summary plots for the estimated fixed and random 
#            effects under the best model and its multivariate normal analogy) 
#            presented in Section 6
#   Input data files: Application/result/fitADNI.RData
#   Output data files: Application/result/fig6.eps
#
################################################################################

load(paste(WD.PATH, 'Application/result/fitADNI.RData', sep=''))

# estimates of random effects for the EFM-MCNLMM
pre.cls = efit32$pre.cls$post.cls
pn = table(pre.cls)
bhat1 = efit32$b$b1[which(pre.cls==1), ] + matrix(rep(efit32$para$Beta[c(1,7),1], pn[1]), ncol=2, byrow=T)
bhat2 = efit32$b$b1[which(pre.cls==2), ] + matrix(rep(efit32$para$Beta[c(1,7),2], pn[2]), ncol=2, byrow=T)
bhat3 = efit32$b$b1[which(pre.cls==3), ] + matrix(rep(efit32$para$Beta[c(1,7),3], pn[3]), ncol=2, byrow=T)
b10 = c(bhat1[,1], bhat2[,1], bhat3[,1])
b20 = c(bhat1[,2], bhat2[,2], bhat3[,2])

# estimate bivariate CN density
bivCN = function(x1, x2, mu=rep(0,2), sigma=diag(2), nu, rho)
{
  p = 2
  sig21 = sigma[2,1]; sig11 = sigma[1,1]; sig22 = sigma[2,2]
  rho12 = sig21 / sqrt(sig11*sig22)
  del = ( (x1-mu[1])^2/sig11 + (x2-mu[2])^2/sig22 - 2*rho12*(x1-mu[1])*(x2-mu[2])/sqrt(sig11*sig22) ) / (1-rho12^2)
  den = 1/((2*pi)^(p/2)*sqrt(det(sigma))) * (nu*rho^(p/2)*exp(-rho*del/2) + (1-nu)*exp(-del/2))
  return(den)
}

est = efit32$para
b110 = seq(min(b10), max(b10), length=100)
b210 = seq(min(b20), max(b20), length=50)
den11 = outer(b110, b210, FUN = bivCN, mu=est$Beta[c(1,7),1], sigma=est$DD[1:2,1:2,1], nu=est$nu[1], rho=est$rho[1])
den12 = outer(b110, b210, FUN = bivCN, mu=est$Beta[c(1,7),2], sigma=est$DD[1:2,1:2,2], nu=est$nu[2], rho=est$rho[2])
den13 = outer(b110, b210, FUN = bivCN, mu=est$Beta[c(1,7),3], sigma=est$DD[1:2,1:2,3], nu=est$nu[3], rho=est$rho[3])

dprob = c(0.8, 0.9, 0.95, 0.975, 0.99)

# mis-classification subject
which(cls != pre.cls)
1 - length(which(cls != pre.cls))/n

# identify bad points
outID = which(efit32$xi$xi1 > 0.99)

betab = NULL
for(j in 1: n) betab = rbind(betab, efit32$b$b1[j,] + efit32$para$Beta[c(1,7), pre.cls[j]])

# estimates of random effects for the EFM-MLMM
Npre.cls = efit12$pre.cls$post.cls
Npn = table(pre.cls)
bhat1n = efit12$b$b1[which(pre.cls==1), ] + matrix(rep(efit12$para$Beta[c(1,7),1], Npn[1]), ncol=2, byrow=T)
bhat2n = efit12$b$b1[which(pre.cls==2), ] + matrix(rep(efit12$para$Beta[c(1,7),2], Npn[2]), ncol=2, byrow=T)
bhat3n = efit12$b$b1[which(pre.cls==3), ] + matrix(rep(efit12$para$Beta[c(1,7),3], Npn[3]), ncol=2, byrow=T)
b10n = c(bhat1n[,1], bhat2n[,1], bhat3n[,1])
b20n = c(bhat1n[,2], bhat2n[,2], bhat3n[,2])

# estimate bivariate CN density
bivN = function(x1, x2, mu=rep(0,2), sigma=diag(2))
{
  p = 2
  sig21 = sigma[2,1]; sig11 = sigma[1,1]; sig22 = sigma[2,2]
  rho12 = sig21 / sqrt(sig11*sig22)
  del = ( (x1-mu[1])^2/sig11 + (x2-mu[2])^2/sig22 - 2*rho12*(x1-mu[1])*(x2-mu[2])/sqrt(sig11*sig22) ) / (1-rho12^2)
  den = 1/((2*pi)^(p/2)*sqrt(det(sigma))) * exp(-del/2)
  return(den)
}

Nest = efit12$para
b110n = seq(min(b10n), max(b10n), length=100)
b210n = seq(min(b20n), max(b20n), length=50)
den11n = outer(b110n, b210n, FUN = bivN, mu=Nest$Beta[c(1,7),1], sigma=Nest$DD[1:2,1:2,1])
den12n = outer(b110n, b210n, FUN = bivN, mu=Nest$Beta[c(1,7),2], sigma=Nest$DD[1:2,1:2,2])
den13n = outer(b110n, b210n, FUN = bivN, mu=Nest$Beta[c(1,7),3], sigma=Nest$DD[1:2,1:2,3])

# mis-classification subject
which(cls != Npre.cls)
1 - length(which(cls != Npre.cls))/n

WD.PATH = paste(getwd(),"/FM-MCNLMM/",sep="")
postscript(paste(WD.PATH, 'Application/result/fig6.eps', sep=''), width=30, height=6)
layout(matrix(c(1,3,2,1,0,4, 5,7,6,5,0,8), 3, 4), widths=c(6.2,3,6,3), heights=c(1.2,3,7))
setcol = c('blue', 'green4', 'red3')

par(mar=c(0, 5, 0.5, 0.5), cex.lab=0.9)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext('(a) EFM-MCNLMM', 1, line=-1.5, cex=1.5, font=2)

par(mar=c(4, 5, 0, 0), cex.lab=0.9)
plot(bhat1[, 1:2], type='n', xlab=expression(beta[i10]+b[ij1]), ylab=expression(beta[i20]+b[ij2]),
     font.lab=2, xlim=c(min(b10), max(b10)), ylim=c(min(b20), max(b20)), cex.lab=1.5, las=1)
contour(b110, b210, den11, levels=quantile(c(den11), prob=dprob), col='cyan', lty=1, lwd=0.75, drawlabels = F, add=T)
contour(b110, b210, den12, levels=quantile(c(den12), prob=dprob), col='green', lty=2, lwd=0.75, drawlabels = F, add=T)
contour(b110, b210, den13, levels=quantile(c(den13), prob=dprob), col='pink', lty=3, lwd=2, drawlabels = F, add=T)
points(bhat1[, 1:2], col='blue', pch=1, cex=0.7, lwd=0.75)
points(bhat2[, 1:2], col='green4', pch=2, cex=0.7, lwd=0.75)
points(bhat3[, 1:2], col='red3', pch=3, cex=0.7, lwd=0.75)
text(betab[outID, 1], betab[outID, 2]+0.015, outID, col=1, cex=0.75)
legend(2, 4.15, c('Group 1', 'Group 2', 'Group 3'), fill=c('white', 'lightyellow', 'lightpink'), border=c('blue', 'green4', 'red3'), bty='n', cex=1)
legend(1.35, 4.15, c('', '', ''), lty=1:3, col=c('cyan','green','pink'), bty='n', lwd=c(0.75, 0.75, 2))
legend(1, 4.15, c('', '', ''), pch=c(1:3), col=c('blue','green4','red3'), bty='n', cex=1)

par(mar=c(0.05, 5, 0, 0.5), cex.lab=0.9)
blist10 = list(b1 = bhat1[,1], b2 = bhat2[,1], b3 = bhat3[,1])
h10 = lapply(blist10, hist, breaks = seq(min(b10), max(b10), length=40), plot=F)
t10 = rbind(h10[[1]]$density, h10[[2]]$density, h10[[3]]$density)
rownames(t10)=names(h10)
colnames(t10)=h10[[1]]$mids
barplot(t10[1,], ylim=c(0,max(t10)), col=0, border='blue', xaxt='n',  yaxt='n', las=1, space=0, xlab='', ylab='')
barplot(t10[2,], ylim=c(0,max(t10)), col='lightyellow', border='green4', xaxt='n', yaxt='n', las=1, space=0, add=T)
barplot(t10[3,], ylim=c(0,max(t10)), col='lightpink', border='red3', xaxt='n', yaxt='n', las=1, space=0, add=T)

par(mar=c(4, 0.05, 0, 0.5), cex.lab=0.9)
blist20 = list(b1 = bhat1[,2], b2 = bhat2[,2], b3 = bhat3[,2])
h20 = lapply(blist20, hist, breaks = seq(min(b20), max(b20), length=40), plot=F)
t20 = rbind(h20[[1]]$density, h20[[2]]$density, h20[[3]]$density)
rownames(t20)=names(h20)
colnames(t20)=h20[[1]]$mids
barplot(t20[2,], horiz =T, xlim=c(0,max(t20)), col='lightyellow', border='green4', xaxt='n', yaxt='n', las=1, space=0, xlab='', ylab='')
barplot(t20[3,], horiz =T, xlim=c(0,max(t20)), col='lightpink', border='red3', xaxt='n', yaxt='n', las=1, space=0, add=T)
barplot(t20[1,], horiz =T, xlim=c(0,max(t20)), col=0, border='blue', xaxt='n',  yaxt='n', las=1, space=0, add=T)

par(mar=c(0, 4.5, 0.5, 0.5), cex.lab=0.9)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext('(b) EFM-MLMM', 1, line=-1.5, cex=1.5, font=2)

par(mar=c(4, 5, 0, 0), cex.lab=0.9)
plot(bhat1n[, 1:2], type='n', xlab=expression(beta[i10]+b[ij1]), ylab=expression(beta[i20]+b[ij2]),
     font.lab=2, xlim=c(min(b10n), max(b10n)), ylim=c(min(b20n), max(b20n)), cex.lab=1.5, las=1)
contour(b110n, b210n, den11n, levels=quantile(c(den11n), prob=dprob), col='cyan', lty=1, lwd=0.75, drawlabels = F, add=T)
contour(b110n, b210n, den12n, levels=quantile(c(den12n), prob=dprob), col='green', lty=2, lwd=0.75, drawlabels = F, add=T)
contour(b110n, b210n, den13n, levels=quantile(c(den13n), prob=dprob), col='pink', lty=3, lwd=2, drawlabels = F, add=T)
points(bhat1n[, 1:2], col='blue', pch=1, cex=0.7, lwd=0.75)
points(bhat2n[, 1:2], col='green4', pch=2, cex=0.7, lwd=0.75)
points(bhat3n[, 1:2], col='red3', pch=3, cex=0.7, lwd=0.75)
legend(2.7, 4.23, c('Group 1', 'Group 2', 'Group 3'), fill=c('white', 'lightyellow', 'lightpink'), border=c('blue', 'green4', 'red3'), bty='n', cex=1)
legend(2.05, 4.23, c('', '', ''), lty=1:3, col=c('cyan','green','pink'), bty='n', lwd=c(0.75, 0.75, 2))
legend(1.7, 4.23, c('', '', ''), pch=c(1:3), col=c('blue','green4','red3'), bty='n', cex=1)

par(mar=c(0.05, 5, 0, 0.5), cex.lab=0.9)
blist10n = list(b1 = bhat1n[,1], b2 = bhat2n[,1], b3 = bhat3n[,1])
h10n = lapply(blist10n, hist, breaks = seq(min(b10n), max(b10n), length=40), plot=F)
t10n = rbind(h10n[[1]]$density, h10n[[2]]$density, h10n[[3]]$density)
rownames(t10n)=names(h10n)
colnames(t10n)=h10n[[1]]$mids
barplot(t10n[1,], ylim=c(0,max(t10n)), col=0, border='blue', xaxt='n',  yaxt='n', las=1, space=0, xlab='', ylab='')
barplot(t10n[2,], ylim=c(0,max(t10n)), col='lightyellow', border='green4', xaxt='n', yaxt='n', las=1, space=0, add=T)
barplot(t10n[3,], ylim=c(0,max(t10n)), col='lightpink', border='red3', xaxt='n', yaxt='n', las=1, space=0, add=T)

par(mar=c(4, 0.05, 0, 0.5), cex.lab=0.9)
blist20n = list(b1 = bhat1n[,2], b2 = bhat2n[,2], b3 = bhat3n[,2])
h20n = lapply(blist20n, hist, breaks = seq(min(b20n), max(b20n), length=40), plot=F)
t20n = rbind(h20n[[1]]$density, h20n[[2]]$density, h20n[[3]]$density)
rownames(t20n)=names(h20n)
colnames(t20n)=h20n[[1]]$mids
barplot(t20n[2,], horiz =T, xlim=c(0,max(t20n)), col='lightyellow', border='green4', xaxt='n', yaxt='n', las=1, space=0, xlab='', ylab='')
barplot(t20n[3,], horiz =T, xlim=c(0,max(t20n)), col='lightpink', border='red3', xaxt='n', yaxt='n', las=1, space=0, add=T)
barplot(t20n[1,], horiz =T, xlim=c(0,max(t20n)), col=0, border='blue', xaxt='n',  yaxt='n', las=1, space=0, add=T)
dev.off()
