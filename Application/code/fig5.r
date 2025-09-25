################################################################################
#
#   Filename: fig5.r
#   Purpose: produce Figure 5 (Confusion matrix plots for classifcation results) 
#            presented in Section 5
#   Input data files: Application/result/fitADNI.RData
#   Output data files: Application/result/fig5.eps
#   Required R packages: ggplot2
#
################################################################################

library(ggplot2)
load(paste(WD.PATH, 'Application/result/fitADNI.RData', sep=''))

# classification table
EN1 = table(factor(efit12$pre.cls$post.cls[which(cls==1)], level=c(1,2,3)))
EN2 = table(factor(efit12$pre.cls$post.cls[which(cls==2)], level=c(1,2,3)))
EN3 = table(factor(efit12$pre.cls$post.cls[which(cls==3)], level=c(1,2,3)))
EMNcls = rbind(EN1, EN2, EN3)
colnames(EMNcls)=1:3

CN1 = table(factor(fit32$pre.cls$post.cls[which(cls==1)], level=c(1,2,3)))
CN2 = table(factor(fit32$pre.cls$post.cls[which(cls==2)], level=c(1,2,3)))
CN3 = table(factor(fit32$pre.cls$post.cls[which(cls==3)], level=c(1,2,3)))
MCNcls = rbind(CN1, CN2, CN3)
colnames(MCNcls)=1:3

ECN1 = table(factor(efit32$pre.cls$post.cls[which(cls==1)], level=c(1,2,3)))
ECN2 = table(factor(efit32$pre.cls$post.cls[which(cls==2)], level=c(1,2,3)))
ECN3 = table(factor(efit32$pre.cls$post.cls[which(cls==3)], level=c(1,2,3)))
EMCNcls = rbind(ECN1, ECN2, ECN3)
colnames(EMCNcls)=1:3

Tab3.2 = m = cbind(MCNcls, EMNcls, EMCNcls)
rownames(Tab3.2) = c('CN', 'MCI', 'Dementia')
print(Tab3.2)

#--- Data ---
m1 <- m[,1:3]
m2 <- m[,4:6]
m3 <- m[,7:9]
p1 <- round(cbind(prop.table(m1), prop.table(m2), prop.table(m3)), 3)
p2 <- round(cbind(prop.table(m1, margin = 2), prop.table(m2, margin = 2), prop.table(m3, margin = 2)), 3)
p3 <- round(cbind(prop.table(m1, margin = 1), prop.table(m2, margin = 1), prop.table(m3, margin = 1)), 3)
bold <- cbind(diag(diag(m[,1:3])), diag(diag(m[,4:6])), diag(diag(m[,7:9])))
bold2 <- cbind(diag(diag(p1[,1:3])), diag(diag(p1[,4:6])), diag(diag(p1[,7:9])))
bold[bold == 0] <- NA
bold2[bold2 == 0] <- NA

df <- data.frame(x = rep(c("CN","MCI","Dementia"),9),
                 y = rep(c(1,1,1,2,2,2,3,3,3),3),
                 z = c(rep("FM-MCNLMM (CCR=57.8%)",9),rep("EFM-MLMM (CCR=61.3%)",9),rep("EFM-MCNLMM (CCR=64.8%)",9)),
                 n = c(m),
                 p1 = c(p1),
                 p2 = c(p2),
                 p3 = c(p3),
                 bold = c(bold),
                 bold2 = c(bold2))

df$x <- factor(df$x, levels = c("Dementia","MCI","CN"))

WD.PATH = paste(getwd(),"/FM-MCNLMM/",sep="")
postscript(paste(WD.PATH, 'Application/result/fig5.eps', sep=""), width = 40, height = 4.5)
ggplot(df, aes(y, x, fill = p1)) +
  geom_tile() +
  geom_text(aes(label = n), size = 5) +
  geom_text(aes(label = paste0("(",scales::percent(p1),")")), size = 5, vjust = 2) +
  geom_text(aes(label = scales::percent(p2)), size = 3, vjust = 6) +
  geom_text(aes(label = scales::percent(p3)), size = 3, vjust = 5, angle = 90) +
  geom_text(aes(label = bold), size = 5, fontface = "bold", colour = "#EC0000") +
  geom_text(aes(label = scales::percent(bold2)), size = 5, vjust = 2, fontface = "bold", colour = "#EC0000") +
  scale_fill_gradient(low = "#FFFFFF", high = "#87CEEB") +
  xlab("Classification") + ylab("AD-Diagnosis") + 
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        legend.title = element_blank(), 
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 12, face = "bold")) +
  facet_wrap(~ z) 
dev.off()
