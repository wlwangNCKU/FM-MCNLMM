################################################################################
#
#   Filename: fig1.r
#   Purpose: produce Figure 1 (trajectory plots for ADAS^0.5 and log10(MidTemp) 
#            scores) presented in Section 2
#   Input data files: Application/result/fitADNI.RData
#   Output data files: Application/result/fig1.eps
#
################################################################################

load(paste(WD.PATH, 'Application/result/fitADNI.RData', sep=''))

# identify bad points
outID = which(efit32$xi$xi1 > 0.99)

sex = numeric(n)
for(j in 1: n) sex[j] = Data$Sex[which(Data$Subject==j)][1]
sex0 = which(sex==0); n0 = length(sex0)
sex1 = which(sex==1); n1 = length(sex1)

outsex = sex[outID]
outsex1 = outID[outsex==1]; no1 = length(outsex1)
outsex0 = outID[outsex==0]; no0 = length(outsex0)

cls.col = c('skyblue', 'green', 'pink')
setcol = c('blue', 'green4', 'red3')

WD.PATH = paste(getwd(),"/FM-MCNLMM/",sep="")
postscript(paste(WD.PATH, 'Application/result/fig1.eps', sep=''), width=25, height=12)
layout(matrix(c(1:6), 2, 3, byrow=T), widths=c(0.9, 6, 6), heights=c(6, 6.3))
par(mar=c(0, 0.5, 2, 0), cex.lab=1.5)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext(expression(ADAS13^0.5), 2, line=-2, cex=1.2, font=2)

par(mar=c(0, 0, 2, 0))
plot(Data$Month[Data$Subject==1], Data$sqADAS13[Data$Subject==1], las=1, type='n', xlab='', ylab=expression(ADAS13^0.5), xaxt='n',
     main='Female', xlim=c(0, max(Data$Month)+10), ylim=c(min(Data$sqADAS13), max(Data$sqADAS13)), cex.main=2)
axis(1, seq(0, 180, 24), labels=F)
for(j in 1: n0){
  points(Data$Month[Data$Subject==sex0[j]], Data$sqADAS13[Data$Subject==sex0[j]], pch=cls[sex0[j]], col=cls.col[cls[sex0[j]]], cex=1.25)
  lines(Data$Month[Data$Subject==sex0[j]], Data$sqADAS13[Data$Subject==sex0[j]], lty=cls[sex0[j]], col=cls.col[cls[sex0[j]]], lwd=1.2)
}
legend("bottomright", c('CN', 'MCI', 'AD'), pch=1:3, lty=1:3, col=cls.col, cex=2, lwd=1.5, bty='n', ncol=3)
for(j in 1: no0){
  points(Data$Month[Data$Subject==outsex0[j]], Data$sqADAS13[Data$Subject==outsex0[j]], pch=cls[outsex0[j]], col=setcol[cls[outsex0[j]]], cex=1.5)
  lines(Data$Month[Data$Subject==outsex0[j]], Data$sqADAS13[Data$Subject==outsex0[j]], lty=cls[outsex0[j]], col=setcol[cls[outsex0[j]]], lwd=2)
  text(rev(Data$Month[Data$Subject==outsex0[j]])[1]+5, rev(Data$sqADAS13[Data$Subject==outsex0[j]])[1], outsex0[j], cex=1.1, col=setcol[cls[outsex0[j]]])
}

par(mar=c(0, 0, 2, 0.5))
plot(Data$Month[Data$Subject==1], Data$sqADAS13[Data$Subject==1], las=1, type='n', xlab='', ylab='', xaxt='n', yaxt='n',
     main='Male', xlim=c(0, max(Data$Month)+10), ylim=c(min(Data$sqADAS13), max(Data$sqADAS13)), cex.main=2)
axis(1, seq(0, 180, 24), labels=F)
axis(2, seq(0, 8, 2), labels=F)
for(j in 1: n1){
  points(Data$Month[Data$Subject==sex1[j]], Data$sqADAS13[Data$Subject==sex1[j]], pch=cls[sex1[j]], col=cls.col[cls[sex1[j]]], cex=1.25)
  lines(Data$Month[Data$Subject==sex1[j]], Data$sqADAS13[Data$Subject==sex1[j]], lty=cls[sex1[j]], col=cls.col[cls[sex1[j]]], lwd=1.2)
}
legend("bottomright", c('CN', 'MCI', 'AD'), pch=1:3, lty=1:3, col=cls.col, cex=2, lwd=1.5, bty='n', ncol=3)
for(j in 1: no1){
  points(Data$Month[Data$Subject==outsex1[j]], Data$sqADAS13[Data$Subject==outsex1[j]], pch=cls[outsex1[j]], col=setcol[cls[outsex1[j]]], cex=1.5)
  lines(Data$Month[Data$Subject==outsex1[j]], Data$sqADAS13[Data$Subject==outsex1[j]], lty=cls[outsex1[j]], col=setcol[cls[outsex1[j]]], lwd=2)
  text(rev(Data$Month[Data$Subject==outsex1[j]])[1]+5, rev(Data$sqADAS13[Data$Subject==outsex1[j]])[1], outsex1[j], cex=1.1, col=setcol[cls[outsex1[j]]])
}

par(mar=c(4, 0.5, 0.5, 0.5), cex.lab=1.2)
plot(c(0, 1), c(0, 3), type='n', xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
mtext(expression(log[10](MidTemp)), 2, line=-2, cex=1.2, font=2)

par(mar=c(4, 0, 0, 0))
plot(Data$Month[Data$Subject==1], Data$log10MidTemp[Data$Subject==1], las=1, type='n',  xlab='', ylab=expression(log[10](MidTemp)), xaxt='n',
     main='', xlim=c(0, max(Data$Month)+10), ylim=c(min(Data$log10MidTemp), max(Data$log10MidTemp)))
axis(1, seq(0, 180, 24), cex.lab=0.5)
mtext('Time (Months)', 1, line=3, cex=1.2, font=1)
for(j in 1: n0){
  points(Data$Month[Data$Subject==sex0[j]], Data$log10MidTemp[Data$Subject==sex0[j]], pch=cls[sex0[j]], col=cls.col[cls[sex0[j]]], cex=1.25)
  lines(Data$Month[Data$Subject==sex0[j]], Data$log10MidTemp[Data$Subject==sex0[j]], lty=cls[sex0[j]], col=cls.col[cls[sex0[j]]], lwd=1.2)
}
legend("bottomright", c('CN', 'MCI', 'AD'), pch=1:3, lty=1:3, col=cls.col, cex=2, lwd=1.5, bty='n', ncol=3)
for(j in 1: no0){
  points(Data$Month[Data$Subject==outsex0[j]], Data$log10MidTemp[Data$Subject==outsex0[j]], pch=cls[outsex0[j]], col=setcol[cls[outsex0[j]]], cex=1.5)
  lines(Data$Month[Data$Subject==outsex0[j]], Data$log10MidTemp[Data$Subject==outsex0[j]], lty=cls[outsex0[j]], col=setcol[cls[outsex0[j]]], lwd=2)
  text(rev(Data$Month[Data$Subject==outsex0[j]])[1]+5, rev(Data$log10MidTemp[Data$Subject==outsex0[j]])[1], outsex0[j], cex=1.1, col=setcol[cls[outsex0[j]]])
}

par(mar=c(4, 0, 0, 0.5))
plot(Data$Month[Data$Subject==1], Data$log10MidTemp[Data$Subject==1], las=1, type='n',  xlab='', ylab='', xaxt='n', yaxt='n',
     main='', xlim=c(0, max(Data$Month)+10), ylim=c(min(Data$log10MidTemp), max(Data$log10MidTemp)))
axis(1, seq(0, 180, 24), cex.lab=0.5)
mtext('Time (Months)', 1, line=3, cex=1.2, font=1)
axis(2, seq(4, 4.4, 0.1), labels=F)
for(j in 1: n1){
  points(Data$Month[Data$Subject==sex1[j]], Data$log10MidTemp[Data$Subject==sex1[j]], pch=cls[sex1[j]], col=cls.col[cls[sex1[j]]], cex=1.25)
  lines(Data$Month[Data$Subject==sex1[j]], Data$log10MidTemp[Data$Subject==sex1[j]], lty=cls[sex1[j]], col=cls.col[cls[sex1[j]]], lwd=1.2)
}
legend("bottomright", c('CN', 'MCI', 'AD'), pch=1:3, lty=1:3, col=cls.col, cex=2, lwd=1.5, bty='n', ncol=3)
for(j in 1: no1){
  points(Data$Month[Data$Subject==outsex1[j]], Data$log10MidTemp[Data$Subject==outsex1[j]], pch=cls[outsex1[j]], col=setcol[cls[outsex1[j]]], cex=1.5)
  lines(Data$Month[Data$Subject==outsex1[j]], Data$log10MidTemp[Data$Subject==outsex1[j]], lty=cls[outsex1[j]], col=setcol[cls[outsex1[j]]], lwd=2)
  text(rev(Data$Month[Data$Subject==outsex1[j]])[1]+5, rev(Data$log10MidTemp[Data$Subject==outsex1[j]])[1], outsex1[j], cex=1.1, col=setcol[cls[outsex1[j]]])
}
dev.off()
