################################################################################
#
#   Filename: TableB2.r
#   Purpose: produce Supplementary Table B.2 in Section 4.2 for Simulation 2
#   Input data files: Simulation/result/ESIM1/fit.txt;
#                     Simulation/result/ESIM1/selection.txt;
#                     Simulation/result/ESIM2/fit.txt;
#                     Simulation/result/ESIM2/selection.txt; 
#                     Simulation/result/ESIM3/fit.txt;
#                     Simulation/result/ESIM3/selection.txt; 
#                     Simulation/result/ESIM4/fit.txt;
#                     Simulation/result/ESIM4/selection.txt; 
#                     Simulation/result/ESIM5/fit.txt;
#                     Simulation/result/ESIM5/selection.txt; 
#                     Simulation/result/ESIM6/fit.txt;
#                     Simulation/result/ESIM6/selection.txt; 
#                     Simulation/result/ESIM7/fit.txt;
#                     Simulation/result/ESIM7/selection.txt; 
#                     Simulation/result/ESIM8/fit.txt;
#                     Simulation/result/ESIM8/selection.txt
#   Output data files: (print the results in R Console)
#
################################################################################

Exper2Tab1 = function(PATH)
{
Tab1 = read.table(paste(PATH, 'fit.txt', sep=""))
Tab2 = read.table(paste(PATH, 'selection.txt', sep=""))[, c(1,6)]
colnames(Tab1) = c('criteria', 'Rep', 'EMN', 'EMT', 'EMCN', 'MN', 'MT', 'MCN')
colnames(Tab2) = c('criteria', 'Model')

Ave0 = rbind(colMeans(Tab1[which(Tab1$criteria == 'BIC'), -c(1:2)]),
             colMeans(Tab1[which(Tab1$criteria == 'ICL'), -c(1:2)]))

Tab0 = data.frame(criteria = rep(c('AIC', 'BIC', 'ICL', 'AWE'), 100),
                  Model = matrix(apply(Tab1[, -c(1:2)], 1, order), nrow=6)[1,])

ans0 = list(no.bic = table(Tab0[which(Tab0$criteria == 'BIC'), 2]),
            no.icl = table(Tab0[which(Tab0$criteria == 'ICL'), 2]))

rownames(Ave0) = c('BIC', 'ICL')
return(list(Average = round(Ave0, 2), Frequency = ans0))
}

# ESIM1
cat('n = 30; covariate-dependent:', '\n')
PATH1 = paste(WD.PATH, 'Simulation/result/ESIM1/', sep='')
print(Exper2Tab1(PATH1))
# ESIM2
cat('n = 75; covariate-dependent:', '\n')
PATH2 = paste(WD.PATH, 'Simulation/result/ESIM2/', sep='')
print(Exper2Tab1(PATH2))
# ESIM3
cat('n = 150; covariate-dependent:', '\n')
PATH3 = paste(WD.PATH, 'Simulation/result/ESIM3/', sep='')
print(Exper2Tab1(PATH3))
# ESIM4
cat('n = 300; covariate-dependent:', '\n')
PATH4 = paste(WD.PATH, 'Simulation/result/ESIM4/', sep='')
print(Exper2Tab1(PATH4))
# ESIM5
cat('n = 30; covariate-independent:', '\n')
PATH5 = paste(WD.PATH, 'Simulation/result/ESIM5/', sep='')
print(Exper2Tab1(PATH5))
# ESIM6
cat('n = 75; covariate-independent:', '\n')
PATH6 = paste(WD.PATH, 'Simulation/result/ESIM6/', sep='')
print(Exper2Tab1(PATH6))
# ESIM7
cat('n = 150; covariate-independent:', '\n')
PATH7 = paste(WD.PATH, 'Simulation/result/ESIM7/', sep='')
print(Exper2Tab1(PATH7))
# ESIM8
cat('n = 300; covariate-independent:', '\n')
PATH8 = paste(WD.PATH, 'Simulation/result/ESIM8/', sep='')
print(Exper2Tab1(PATH8))
