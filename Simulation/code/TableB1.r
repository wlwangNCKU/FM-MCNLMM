################################################################################
#
#   Filename: TableB1.r
#   Purpose: produce Supplementary Table B.1 in Section 4.1 for Simulation 1
#   Input data files: Simulation/result/SIM1/fit.txt;
#                     Simulation/result/SIM1/selection.txt;
#                     Simulation/result/SIM2/fit.txt;
#                     Simulation/result/SIM2/selection.txt; 
#                     Simulation/result/SIM3/fit.txt;
#                     Simulation/result/SIM3/selection.txt; 
#                     Simulation/result/SIM4/fit.txt;
#                     Simulation/result/SIM4/selection.txt; 
#   Output data files: (print the results in R Console)
#
################################################################################

# SIM1
PATH1 = paste(WD.PATH, 'Simulation/result/SIM1/', sep='')
Tab11 = read.table(paste(PATH1, 'fit.txt', sep=""))
Tab21 = read.table(paste(PATH1, 'selection.txt', sep=""))

# SIM2
PATH2 = paste(WD.PATH, 'Simulation/result/SIM2/', sep='')
Tab12 = read.table(paste(PATH2, 'fit.txt', sep=""))
Tab22 = read.table(paste(PATH2, 'selection.txt', sep=""))

# SIM3
PATH3 = paste(WD.PATH, 'Simulation/result/SIM3/', sep='')
Tab13 = read.table(paste(PATH3, 'fit.txt', sep=""))
Tab23 = read.table(paste(PATH3, 'selection.txt', sep=""))

# SIM4
PATH4 = paste(WD.PATH, 'Simulation/result/SIM4/', sep='')
Tab14 = read.table(paste(PATH4, 'fit.txt', sep=""))
Tab24 = read.table(paste(PATH4, 'selection.txt', sep=""))

name1 = c('criteria', 'Rep', 'MN1', 'MN2', 'MN3', 'MT1', 'MT2', 'MT3', 'MCN1', 'MCN2', 'MCN3')
name2 = c('criteria', 'Rep', 'NG', 'TG', 'CNG', 'Distr', 'NorCN')

colnames(Tab11) = colnames(Tab12) = colnames(Tab13) = colnames(Tab14) = name1
colnames(Tab21) = colnames(Tab22) = colnames(Tab23) = colnames(Tab24) = name2

Exper1Tab1 = function(Tab1, Tab2)
{
summTab = rbind(colMeans(Tab1[which(Tab1$criteria == 'BIC'), -c(1:2)]),
                colMeans(Tab1[which(Tab1$criteria == 'ICL'), -c(1:2)]))
rownames(summTab) = c('BIC', 'ICL')

no.bic = apply(Tab2[which(Tab2$criteria == 'BIC'), -c(1:2)], 2, table)
no.icl = apply(Tab2[which(Tab2$criteria == 'ICL'), -c(1:2)], 2, table)
return(list(SummTab = round(summTab, 2), listTab =list(bic=no.bic, icl=no.icl)))
}

# Report the numerical results
# SIM1
cat('n = 20; (nu, rho)=(0.25, 0.2):', '\n')
print(Exper1Tab1(Tab11, Tab21))
# SIM2
cat('n = 50; (nu, rho)=(0.25, 0.2):', '\n')
print(Exper1Tab1(Tab12, Tab22))
# SIM3
cat('n = 100; (nu, rho)=(0.25, 0.2):', '\n')
print(Exper1Tab1(Tab13, Tab23))
# SIM4
cat('n = 200; (nu, rho)=(0, 1):', '\n')
print(Exper1Tab1(Tab14, Tab24))
