################################################################################
#
#   Filename: Tables1-2.r
#   Purpose: produce Table 1 (model selection criteria) and Table 2 (parameter 
#                    estimates and standard errors) presented in Section 5
#   Input data files: Application/result/fitADNI.RData
#   Output data files: Application/result/Table1.csv;
#                      Application/result/Table2a.csv;
#                      Application/result/Table2b.csv
#
################################################################################

load(paste(WD.PATH, 'Application/result/fitADNI.RData', sep=''))

# model selection criteria 
RI = rbind(
c(fit10$model.inf$m, fit11$model.inf$m, fit12$model.inf$m, fit30$model.inf$m, fit31$model.inf$m, fit32$model.inf$m, 
  efit10$model.inf$m, efit11$model.inf$m, efit12$model.inf$m, efit30$model.inf$m, efit31$model.inf$m, efit32$model.inf$m),
c(fit10$model.inf$loglik, fit11$model.inf$loglik, fit12$model.inf$loglik, fit30$model.inf$loglik, fit31$model.inf$loglik, fit32$model.inf$loglik, 
  efit10$model.inf$loglik, efit11$model.inf$loglik, efit12$model.inf$loglik, efit30$model.inf$loglik, efit31$model.inf$loglik, efit32$model.inf$loglik),
c(fit10$model.inf$bic, fit11$model.inf$bic, fit12$model.inf$bic, fit30$model.inf$bic, fit31$model.inf$bic, fit32$model.inf$bic, 
  efit10$model.inf$bic, efit11$model.inf$bic, efit12$model.inf$bic, efit30$model.inf$bic, efit31$model.inf$bic, efit32$model.inf$bic),
c(fit10$model.inf$icl, fit11$model.inf$icl, fit12$model.inf$icl, fit30$model.inf$icl, fit31$model.inf$icl, fit32$model.inf$icl, 
  efit10$model.inf$icl, efit11$model.inf$icl, efit12$model.inf$icl, efit30$model.inf$icl, efit31$model.inf$icl, efit32$model.inf$icl)) 

rowName = c('m', 'loglik', 'BIC', 'ICL')
colName = c('MN-UNC', 'MN-AR1', 'MN-DEC', 'MCN-UNC', 'MCN-AR1', 'MCN-DEC', 
            'EMN-UNC', 'EMN-AR1', 'EMN-DEC', 'EMCN-UNC', 'EMCN-AR1', 'EMCN-DEC') 
rownames(RI) = rowName
colnames(RI) = colName 
print(round(RI, 2))
WD.PATH = paste(getwd(),"/FM-MCNLMM/",sep="")
write.csv(round(RI, 2), paste(WD.PATH, 'Application/result/Table1.csv', sep=""), row.names = T)

# parameter estimates based on the best fitted model
Tab.est = cbind(t(efit32$est$out[, 1:22]), t(efit32$est$out[, 23:44]), t(efit32$est$out[, 45:66]))
Tab.psi = t(efit32$est$out[, 67:72]) 

colnames(Tab.est) = rep(c('Est', 'SD'),3)
rownames(Tab.est) = c(paste('beta.i', c(10,11,12,13,14,15, 20,21,22,23,24,25), sep=''), 
                      paste('D.i', c(11,21,22), sep=''), 
                      paste('Sigma.i', c(11,12,22), sep=''), 'phi.i', 'ga.i', 'nu.i', 'rho.i')
colnames(Tab.psi) = c('Est', 'SD')
rownames(Tab.psi) = paste('psi', c(10,11,12, 20,21,22), sep='')
print(round(Tab.est, 4))
print(round(Tab.psi, 4))

write.csv(round(Tab.est, 4), paste(WD.PATH, 'Application/result/Table2a.csv', sep=""), row.names = T)
write.csv(round(Tab.psi, 4), paste(WD.PATH, 'Application/result/Table2b.csv', sep=""), row.names = T)
 