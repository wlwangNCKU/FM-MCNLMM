rm(list=ls())
WD.PATH = paste(getwd(),"/FM-MCNLMM/",sep="")

# Example: Tables 1 and 2 # 
source(paste(WD.PATH, 'Application/code/adni.run.r', sep=''))
source(paste(WD.PATH, 'Application/code/Tables1-2.r', sep=''))

# Example: Figure 1 # 
source(paste(WD.PATH, 'Application/code/fig1.r', sep=''))

# Example: Figure 2 # 
source(paste(WD.PATH, 'Application/code/fig2.r', sep=''))

# Example: Figure 5 # 
source(paste(WD.PATH, 'Application/code/fig5.r', sep=''))

# Example: Figure 6 # 
source(paste(WD.PATH, 'Application/code/fig6.r', sep=''))

# Example: Supplementary Figure B.4 # 
source(paste(WD.PATH, 'Application/code/figB4.r', sep=''))

# Simulation 1 #
## SIM1 ##
Rep = 1; total.rep = 100; pd = 'SIM1' 
n = 20; x = 'X1'; nu = rep(0.25, 2); rho = rep(0.2, 2)
source(paste(WD.PATH, 'Simulation/code/sim.r', sep=''))

## SIM2 ## 
Rep = 1; total.rep = 100; pd = 'SIM2'       
n = 50; x = 'X2'; nu = rep(0.25, 2); rho = rep(0.2, 2)
source(paste(WD.PATH, 'Simulation/code/sim.r', sep=''))

## SIM3 ## 
Rep = 1; total.rep = 100; pd = 'SIM3'
n = 100; x = 'X3'; nu = rep(0.25, 2); rho = rep(0.2, 2)
source(paste(WD.PATH, 'Simulation/code/sim.r', sep=''))

## SIM4 ##
Rep = 1;  total.rep = 100; pd = 'SIM4' 
n = 200; x = 'X4'; nu = rep(0.25, 2); rho = rep(0.2, 2)
source(paste(WD.PATH, 'Simulation/code/sim.r', sep=''))

# Simulation 1: Figure 3 # 
source(paste(WD.PATH, 'Simulation/code/fig3.r', sep=''))

# Simulation 1: Supplementary Figure B.1 # 
source(paste(WD.PATH, 'Simulation/code/figB1.r', sep=''))

# Simulation 1: Supplementary Table B.1 # 
source(paste(WD.PATH, 'Simulation/code/TableB1.r', sep=''))

# Fix covariate: generate just one time #
n = 20 
#n = 50 
#n = 100 
#n = 200 
x = round(rnorm(n), 3)
#x = sample(0:1, n, replace=T)
tj = 10
aX = NULL
for(j in 1: n){ 
  Xj = cbind(rep(j, tj), rep(1, tj), 1: tj, rep(x[j], tj))
  aX = rbind(aX, Xj)
}
Unbal = sample(1: (n*tj), (n*tj*0.2), replace=F)
raX = aX[-Unbal, ]
write.table(raX, paste(WD.PATH, 'Simulation/result/X1.txt', sep=""), row.names = F, col.names = c('Subj', 'x0', 'x1', 'x2'))
#write.table(raX, paste(WD.PATH, 'Simulation/result/X2.txt', sep=""), row.names = F, col.names = c('Subj', 'x0', 'x1', 'x2'))
#write.table(raX, paste(WD.PATH, 'Simulation/result/X3.txt', sep=""), row.names = F, col.names = c('Subj', 'x0', 'x1', 'x2'))
#write.table(raX, paste(WD.PATH, 'Simulation/result/X4.txt', sep=""), row.names = F, col.names = c('Subj', 'x0', 'x1', 'x2'))

# Simulation 2 #
## Case 1 ## 
## ESIM1 ## 
Rep = 1; total.rep = 100
pd = 'ESIM1'; n = 30; nu = rep(0.25, 3); rho = rep(0.2, 3)
psi = c(6, -1, 8, -1.5)
source(paste(WD.PATH, 'Simulation/code/simExd.r', sep=''))

## ESIM2 ##
Rep = 1; total.rep = 100
pd = 'ESIM2'; n = 75; nu = rep(0.25, 3); rho = rep(0.2, 3)
psi = c(6, -1, 8, -1.5)
source(paste(WD.PATH, 'Simulation/code/simExd.r', sep=''))

## ESIM3 ## 
Rep = 1; total.rep = 100
pd = 'ESIM3'; n = 150; nu = rep(0.25, 3); rho = rep(0.2, 3)
psi = c(6, -1, 8, -1.5)
source(paste(WD.PATH, 'Simulation/code/simExd.r', sep=''))

## ESIM4 ## 
Rep = 1; total.rep = 100
pd = 'ESIM4'; n = 300; nu = rep(0.25, 3); rho = rep(0.2, 3)
psi = c(6, -1, 8, -1.5)
source(paste(WD.PATH, 'Simulation/code/simExd.r', sep=''))

## Case 2 ## 
## ESIM5 ## 
Rep = 1; total.rep = 100
pd = 'ESIM5'; n = 30; nu = rep(0.25, 3); rho = rep(0.2, 3)
psi = rep(0, 4)
source(paste(WD.PATH, 'Simulation/code/simExd.r', sep=''))

## ESIM6 ## 
Rep = 1; total.rep = 100
pd = 'ESIM6'; n = 75; nu = rep(0.25, 3); rho = rep(0.2, 3)
psi = rep(0, 4)
source(paste(WD.PATH, 'Simulation/code/simExd.r', sep=''))

## ESIM7 ## 
Rep = 1; total.rep = 100
pd = 'ESIM7'; n = 150; nu = rep(0.25, 3); rho = rep(0.2, 3)
psi = rep(0, 4)
source(paste(WD.PATH, 'Simulation/code/simExd.r', sep=''))

## ESIM8 ## 
Rep = 1; total.rep = 100
pd = 'ESIM8'; n = 300; nu = rep(0.25, 3); rho = rep(0.2, 3)
psi = rep(0, 4)
source(paste(WD.PATH, 'Simulation/code/simExd.r', sep=''))

# Simultion 2: Figure 4 # 
source(paste(WD.PATH, 'Simulation/code/fig4.r', sep=''))

# Simulation 2: Supplementary Figure B.2 & Supplementary Figure B.3 # 
source(paste(WD.PATH, 'Simulation/code/figB2B3.r', sep=''))

# Simulation 2: Supplementary Table B.2 # 
source(paste(WD.PATH, 'Simulation/code/TableB2.r', sep=''))
