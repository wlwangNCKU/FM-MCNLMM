# FM-MCNLMM
Supplement: Data and Code for "Grouped Multi-trajectory Modeling Using Finite Mixtures of Multivariate Contaminated Normal Linear Mixed Model"

#################################################################################################################################################
Source R codes and data for the manuscript: 
"Grouped Multi-trajectory Modeling Using Finite Mixtures of Multivariate Contaminated Normal Linear Mixed Model",
by Tsung-I Lin and Wan-Lun Wang*
#################################################################################################################################################

# Author responsible for the code #
For questions, comments or remarks about the code please contact responsible author, Wan-Lun Wang (wangwl@gs.ncku.edu.tw).

# Configurations #
The code was written/evaluated in R with the following software versions:
R version 4.1.1 (2021-08-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Traditional)_Taiwan.950  LC_CTYPE=Chinese (Traditional)_Taiwan.950   
[3] LC_MONETARY=Chinese (Traditional)_Taiwan.950 LC_NUMERIC=C                                
[5] LC_TIME=Chinese (Traditional)_Taiwan.950    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] ggplot2_3.4.3 ellipse_0.5.0 car_3.1-2     carData_3.0-5 nlme_3.1-152 

loaded via a namespace (and not attached):
 [1] magrittr_2.0.1   tidyselect_1.2.0 munsell_0.5.0    colorspace_2.1-0 lattice_0.20-44  R6_2.5.1        
 [7] rlang_1.1.0      fansi_0.5.0      dplyr_1.1.2      tools_4.1.1      grid_4.1.1       gtable_0.3.3    
[13] utf8_1.2.2       cli_3.6.1        withr_2.5.0      abind_1.4-5      tibble_3.2.1     lifecycle_1.0.3 
[19] farver_2.1.1     vctrs_0.6.1      glue_1.6.2       labeling_0.4.2   compiler_4.1.1   pillar_1.9.0    
[25] generics_0.1.2   scales_1.2.1     pkgconfig_2.0.3 

# Descriptions of the codes # 
Please extract the file "Data-and-Code-FM-MCNLMM.zip" to the "current working directory" of the R package.
The getwd() function shall determine an absolute pathname of the "current working directory".

Before running the codes **adni.run.r**, **fig2.r**, **fig4.r**, and **fig5.r**, one needs to install the following R packages:

    install.packages("car") Version 3.1-2
    install.packages("ellipse") Version 0.5.0
    install.packages("ggplot2") Version 3.4.3
    install.packages("mvtnorm")  Version 1.1-3
    install.packages("nlme")  Version 3.1-152
    install.packages("vioplot") Version 0.5.0
 
R codes for the implementation of the proposed methodology are provided.

# Data and Code for Application #
## Subfolder: ./Application/code ##
`./Application/code`
       contains main scripts of 
        - (1) **adni.run.r**: main script for reproducting the fitting results of 12 candidate models to the subset of ADNI data;
        - (2) **fig1.r**: main script for reproducting Figure 1 (run 'adni.run.r' first or load 'fitADNI.RData' directly, and then run 'fig1.r');
	    - (3) **fig2.r** main script for reproducting supplementary Figure 2 (read 'ADNIMERGE.csv' first, and then run 'fig2.r');
	    - (4) **fig5.r** main script for reproducting Figure 5 (load 'fitADNI.RData' first, and then run 'fig5.r');
	    - (5) **fig6.r**: main script for reproducting Figure 6 (load 'fitADNI.RData' first, and then run 'fig6.r');
	    - (6) **figB4.r**: main script for reproducting Figure B.4 (load 'fitADNI.RData' first, and then run 'figB4.r');
        - (7) **Tables1-2.r**: main script for reproducting Table 1 and Table 2 presented in Section 5.

## Subfolder: ./Application/data ##
`./Application/data`
      contains
      - **ADNIMERGE.csv**: the dataset from the Alzheimer's Disease Neuroimaging Initiative (ADNI) study.

## Subfolder: ./Application/function ##
`./Application/function`
 	 contains the program (function) of
     - (1) **basicfn.r**: main script collecting some basic functions for computing autocorrelation function, mis-classification rate, etc;
	 - (2) **fmmlmm.fn.r**: main script for running the AECM algorithm for fitting the FM-MLMM;
	 - (3) **fmmcnlmm.fn.r**: main script for running the AECM algorithm for fitting the FM-MCNLMM;
	 - (4) **efmmlmm.fn.r**: main script for running the AECM algorithm for fitting the EFM-MLMM;
	 - (5) **efmmcnlmm.fn.r**: main script for running the AECM algorithm for fitting the EFM-MCNLMM.

## Subfolder: ./Application/result ##
`./Application/result`
        contains
	- (1) **fitADNI.RData**: analysis results for the ADNI data based on 12 candidate models;
	- (2) **fig1.eps**: (Figure 1) trajectory plots for ADAS^0.5 and log10(MidTemp) scores for the cognitively normal (CN), mild cognitive impairment (MCI), and Alzheimers disease (AD) patients;
    - (3) **fig2.eps**: (Figure 2) scatter plots, the 95% confidence ellipses, histogram and boxplots of empirical Bayes estimates for random effects and residuals obtained by the fitted linear mixed-effects models;
	- (4) **fig5.eps**: (Figure 4) confusion matrix plots for classifcation results for the ADNI data based on three candidate models;
	- (5) **fig6.eps**: (Figure 5) scatter plots and summary histograms for the estimated fixed and random effects superimposed on a set of contour lines of the bivariate contaminated normal densities of each group;
	- (6) **figB4.eps**: (Figure B.4) fitted mean curves of the two outcomes across groups based on the 2-component EFM-MCNLMM with DEC errors;
	- (7) **Table1.csv**: (Table 1) summary of model selection criteria for the 12 candidate models;
	- (8) **Table2a.csv** & **Table2b.csv**: (Table 2) summary of parameter estimates along with their standard errors under the best fitted model.

# Data and Code for Simulation #
## Subfolder: ./Simulation/code ##
`./Simulation/code`
      contains
	- (1) **sim.r**: main script for re-generate simulation results for Simulation 1 with two (nu, rho) scenarios under various sample sizes;
	- (2) **fig3.r**: main script for reproducing Figure 3 that shows the results of comparision of CCR scores for Simulation 1;
	- (3) **figB1.r**: main script for reproducing Supplementary Figure B.1 that show the results of biases obtained by the fitted FM-MCNLMM for Simulation 1;
	- (4) **TableB1.r**: main script for reproducing Supplementary Table B.1 that prints the results of model selection criteria (BIC and ICL) for Simulation 1;
 	- (5) **simExd.r**: main script for re-generate simulation results for Simulation 2 with two psi scenarios under various sample sizes;
	- (6) **fig4.r**: main script for reproducing Figure 4 that shows the results of comparision of CCR scores for Simulation 2;
	- (7) **figB2B3.r**: main script for reproducing Supplementary Figure B.1 that show the results of biases obtained by the fitted FM-MCNLMM for Simulation 1;
	- (8) **TableB2.r**: main script for reproducing Supplementary Table B.2 that prints the results of model selection criteria (BIC and ICL) for Simulation 2.

## Subfolder: ./Simulation/function ##
`./Simulation/function`
 	 contains the program (function) of
     - (1) **basicfn.r**: main script for collecting some basic functions for computing autocorrelation function, mis-classification rate, etc;
	 - (2) **fmmlmm.fn.r**: main script for running the AECM algorithm for fitting the FM-MLMM;
	 - (3) **fmmtlmm.fn.r**: main script for running the AECM algorithm for fitting the FM-MtLMM;
	 - (4) **fmmcnlmm.fn.r**: main script for running the AECM algorithm for fitting the FM-MCNLMM;
	 - (5) **efmmlmm.fn.r**: main script for running the AECM algorithm for fitting the EFM-MLMM;
	 - (6) **efmmtlmm.fn.r**: main script for running the AECM algorithm for fitting the EFM-MtLMM;
	 - (7) **efmmcnlmm.fn.r**: main script for running the AECM algorithm for fitting the EFM-MCNLMM.

## Subfolder: ./Simulation/result ##
`./Simulation/result`
        contains
	- (1) 4 subsubfolders: `./SIM1`, `./SIM2`, `./SIM3`, and `./SIM4`,  storing the **class.txt**, **biasC.txt**, **estC.txt**, **fit.txt**, and **selection.txt** text files for Simulation 1;
	- (2) 8 subsubfolders: `./ESIM1`, `./ESIM2`, `./ESIM3`, `./ESIM4`, `./ESIM5`, `./ESIM6`, `./ESIM7`, and `./ESIM8`,  storing the **class.txt**, **estEN.txt**, **estET.txt**, **estEC.txt**, **fit.txt**, and **selection.txt** text files for Simulation 2;
	- (3) **fig3.eps**: (Figure 3) boxplots for the correct classification rates (CCR) obtained by fitting the FM-MLMM, FM-MtLMM and FM-MCNLMM;
	- (4) **fig4.eps**: (Figure 3) split violin plots for CCR obtained by fitting the 3-component EFM-MCNLMM and FM-MCNLMM;
	- (5) **figB1.eps**: (Supplementary Figure B.1) boxplots for biases of parameter estimates across sample sizes;
	- (6) **figB2a1.eps**, **figB2a2.eps**, **figB2a3.eps**: (Supplementary Figure B.2) boxplots for biases of parameter estimates obtained by the EFM-MLMM, EFM-MtLMM and EFM-MCNLMM to the data with covaraite-dependent weights across sample sizes;
	- (7) **figB3a1.eps**, **figB3a2.eps**, **figB3a3.eps**: (Supplementary Figure B.3) boxplots for biases of parameter estimates obtained by the EFM-MLMM, EFM-MtLMM and EFM-MCNLMM to the data with covaraite-independent weights across sample sizes;
	- (8) **X1.txt**, **X2.txt**, **X3.txt**, **X4.txt**: design matrices of fixed effects for n=20, 50, 100, and 200 subjects pre-specified for Simulation 1.

## Additional Remark ##
 - Note (1): One can directly run each "source(.)" described in **master.r** file in the seperate R session to obtain the results.
 - Note (2): The fitting results of the considered models for the ADNI dataset obtained by running **adni.run.r** have been stored in `.Application/result/fitADNI.Rdata`.
 - Note (3): To draw Figures 1, 5, 6, and B.4 in this paper, please load the **fitADNI.RData** file in subfolder `./Application/result/`, and then run the **fig1.r**, **fig5.r** and **fig6.r** scripts in subfolder `./Application/code/`. 
 - Note (4): Since 'adni.run.r' takes a long time to run, to reproduce Tables 1 and 2 in Section 5, please load the 'fitADNI.RData' file in subfolder "./Application/result/", and then run **Tables1-2.r** script in subfolder `./Application/code/`.
 - Note (5): Since 'sim.r' takes a long time to run, to reproduce numerical results in Section 4.1, we record these intermediately numerical results so that one can use the R codes **fig3.r**, **figB1.r**, and **TableB1.r** to obtain the final results based on files stored in `./Simulation/result/SIM1`, `./Simulation/result/SIM2`, `./Simulation/result/SIM3`, and `./Simulation/result/SIM4` subfolders.
 - Note (6): Since 'simExd.r' takes a long time to run, to reproduce numerical results in Section 4.2, we record these intermediately numerical results so that one can use the R codes **fig4.r**,  **figB2B3.r**, and **TableB2.r** to obtain the final results based on files stored in `./Simulation/result/ESIM1`, `./Simulation/result/ESIM2`, `./Simulation/result/ESIM3`, `./Simulation/result/ESIM4`, `./Simulation/result/ESIM5`, `./Simulation/result/ESIM6`, `./Simulation/result/ESIM7`, and `./Simulation/result/ESIM8` subfolders.
