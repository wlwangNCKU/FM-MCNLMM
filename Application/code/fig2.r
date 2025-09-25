################################################################################
#
#   Filename: fig2.r
#   Purpose: produce Figure 2 (summary plots for empirical Bayes estimates for
#            random effects and residuals) presented in Section 2
#   Input data files: Application/data/ADNIMERGE.csv
#   Output data files: Application/result/fig2.eps
#   Required R packages: nlme; car; ellipse 
#
################################################################################

library(nlme)
library(car)
library(ellipse)

###### import data ######
ADNI = read.csv(paste(WD.PATH, 'Application/data/ADNIMERGE.csv', sep=''), header = T)
Subj = sort(ADNI[which(ADNI$COLPROT=='ADNI1' & ADNI$VISCODE=='bl' & ADNI$PTRACCAT=='White' & ADNI$PTMARRY == 'Married'), 'RID'])
N = length(Subj)

# Reorganized Data: Remove the row with common visited month of a given subject
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
table(cls)
cls.col = c('blue', 'green4', 'red3')

### sqADAS13 ###
fm1 = lme(sqADAS13 ~ Time + Time2 + AGE + PTEDUCAT + Sex, data=Data, random=~1, corAR1(), method='ML')
fix11 = fm1$residuals[which(cls==1),1]; res11 = fm1$residuals[which(cls==1),2]
fix12 = fm1$residuals[which(cls==2),1]; res12 = fm1$residuals[which(cls==2),2]
fix13 = fm1$residuals[which(cls==3),1]; res13 = fm1$residuals[which(cls==3),2]

b11 = as.numeric(fm1$coefficients$random$Subject[which(cls==1),1])
b12 = as.numeric(fm1$coefficients$random$Subject[which(cls==2),1])
b13 = as.numeric(fm1$coefficients$random$Subject[which(cls==3),1])

### log10MidTemp ###
fm2 = lme(log10MidTemp ~ Time + Time2 + AGE + PTEDUCAT + Sex, data=Data, random=~1, corAR1(), method='ML')
fix21 = fm2$residuals[which(cls==1),1]; res21 = fm2$residuals[which(cls==1),2]
fix22 = fm2$residuals[which(cls==2),1]; res22 = fm2$residuals[which(cls==2),2]
fix23 = fm2$residuals[which(cls==3),1]; res23 = fm2$residuals[which(cls==3),2]

b21 = as.numeric(fm2$coefficients$random$Subject[which(cls==1),1])
b22 = as.numeric(fm2$coefficients$random$Subject[which(cls==2),1])
b23 = as.numeric(fm2$coefficients$random$Subject[which(cls==3),1])

WD.PATH = paste(getwd(),"/FM-MCNLMM/",sep="")
postscript(paste(WD.PATH, 'Application/result/fig2.eps', sep=''), width=30, height=10)
layout(matrix(c(4,4,3,3,2,2,1, rep(c(rep(8,3),rep(12,3),5), 3), rep(c(rep(9,3),rep(13,3),6),3), 
              rep(c(rep(10,3),rep(14,3),7),3), rep(c(rep(11,3),rep(15,3),0),3)), 7, 13), widths=c(5,rep(0.75,12)), heights=c(rep(1,6),5.5))
# scatter plot for random intercepts
par(mar=c(3.5, 4, 0, 0), cex.lab=0.9)
plot(b11, b21, pch=1, col='blue', cex=0.9, lwd=0.75, las=1, xlab='', ylab='', 
     xlim=c(min(c(b11, b12, b13)), max(c(b11, b12, b13))), ylim=c(min(c(b21, b22, b23)), max(c(b21, b22, b23))))
mtext(expression(b[j1]), 1, line=2.5, cex=0.9, font=1)
mtext(expression(b[j2]), 2, line=2.5, cex=0.9, font=1)
points(b12, b22, pch=2, col='green4', cex=0.9, lwd=0.75)
points(b13, b23, pch=3, col='red3', cex=0.9, lwd=0.75)
abline(v = c(mean(b11), mean(b12), mean(b13)), lty=c(1,2,3), col=c('cyan', 'green4', 'lightpink')) 
abline(h = c(mean(b21), mean(b22), mean(b23)), lty=c(1,2,3), col=c('cyan', 'green4', 'lightpink')) 
lines(ellipse(cor(b11, b21), scale=c(sd(b11), sd(b21)), centre = colMeans(cbind(b11, b21)), level=0.95), lty=1, col='blue', type='l', lwd=2)
lines(ellipse(cor(b12, b22), scale=c(sd(b12), sd(b22)), centre = colMeans(cbind(b12, b22)), level=0.95), lty=2, col='green4', type='l', lwd=2)
lines(ellipse(cor(b13, b23), scale=c(sd(b13), sd(b23)), centre = colMeans(cbind(b13, b23)), level=0.95), lty=3, col='red3', type='l', lwd=2)
legend("bottomleft", c('CN', 'MCI', 'AD'), pch=c(1:3), lty=c(1:3), lwd=2, col=c('blue','green4','red3'), bty='n', cex=2)

par(mar=c(0, 4, 0.2, 0), cex.lab=0.9)
bhat1 = list(b1 = b11, b2 = b12, b3 = b13)
h1 = lapply(bhat1, hist, breaks = seq(min(fm1$coefficients$random$Subject[,1]), max(fm1$coefficients$random$Subject[,1]), length=30), plot=F)
t1 = rbind(h1[[1]]$density, h1[[2]]$density, h1[[3]]$density)
rownames(t1) = names(h1)
colnames(t1) = h1[[1]]$mids
barplot(t1[1,], ylim=c(0,max(t1)), col=0, border='blue', xaxt='n', yaxt='n', las=1, space=0)
legend("topleft", 'CN', fill='white', border='blue', bty='n', cex=2)
barplot(t1[2,], ylim=c(0,max(t1)), col='lightyellow', border='green4', xaxt='n', yaxt='n', las=1, space=0, xlab='', ylab='')
legend("topleft", 'MCI', fill='lightyellow', border='green4', bty='n', cex=2)
barplot(t1[3,], ylim=c(0,max(t1)), col='lightpink', border='red3', xaxt='n', yaxt='n', las=1, space=0)
legend("topleft", 'AD', fill='lightpink', border='red3', bty='n', cex=2)

par(mar=c(3.5, 0, 0, 0.2), cex.lab=0.9)
bhat2 = list(b1 = b21, b2 = b22, b3 = b23)
h2 = lapply(bhat2, hist, breaks = seq(min(fm2$coefficients$random$Subject[,1]), max(fm2$coefficients$random$Subject[,1]), length=30), plot=F)
t2 = rbind(h2[[1]]$density, h2[[2]]$density, h2[[3]]$density)
rownames(t2) = names(h2)
colnames(t2) = h2[[1]]$mids
barplot(t2[1,], horiz =T, xlim=c(0,max(t2)), col=0, border='blue', xaxt='n', yaxt='n', las=1, space=0)
legend("bottomright", 'CN', fill='white', border='blue', bty='n', cex=2)
barplot(t2[2,], horiz =T, xlim=c(0,max(t2)), col='lightyellow', border='green4', xaxt='n',  yaxt='n', las=1, space=0, xlab='', ylab='')
legend("bottomright", 'MCI', fill='lightyellow', border='green4', bty='n', cex=2)
barplot(t2[3,], horiz =T, xlim=c(0,max(t2)), col='lightpink', border='red3', xaxt='n',  yaxt='n', las=1, space=0)
legend("bottomright", 'AD', fill='lightpink', border='red3', bty='n', cex=2)

par(mar=c(3, 4, 1.5, 0.5), cex.lab=0.9)
boxplot(list(CN=b11, MCI=b12, AD=b13), col=c('white','lightyellow','lightpink'), border=c('blue', 'green4', 'red3'), las=1, main=expression(hat(b)[j1]))
boxplot(list(CN=b21, MCI=b22, AD=b23), col=c('white','lightyellow','lightpink'), border=c('blue', 'green4', 'red3'), las=1, main=expression(hat(b)[j2]))
boxplot(list(CN=res11, MCI=res12, AD=res13), col=c('white','lightyellow','lightpink'), border=c('blue', 'green4', 'red3'), las=1, main=expression(hat(e)[j1t]))
boxplot(list(CN=res21, MCI=res22, AD=res23), col=c('white','lightyellow','lightpink'), border=c('blue', 'green4', 'red3'), las=1, main=expression(hat(e)[j1t]))

b1.hat = fm1$coefficients$random$Subject
qqnorm(b1.hat, main=expression(paste('Q-Q plot for ', hat(b)[j1], sep='')), pch=16, cex=0.3, las=1)
qqline(b1.hat)
b2.hat = fm2$coefficients$random$Subject
qqnorm(b2.hat, main=expression(paste('Q-Q plot for ', hat(b)[j2], sep='')), pch=16, cex=0.3, las=1)
qqline(b2.hat)
e1.hat = fm1$residuals
qqnorm(e1.hat, main=expression(paste('Q-Q plot for ', hat(e)[j1t], sep='')), pch=16, cex=0.3, las=1)
qqline(e1.hat)
e2.hat = fm2$residuals
qqnorm(e2.hat, main=expression(paste('Q-Q plot for ', hat(e)[j2t], sep='')), pch=16, cex=0.3, las=1)
qqline(e2.hat)
dev.off()
