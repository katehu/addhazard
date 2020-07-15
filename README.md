# addhazard

addhazard 1.2.0 updates:
1. allow either weight or selection probability to be the input argument
2. add row names for coefs in the summary table
3. break ties for users when ties = 'break'
4. generate wts variable in the example nwts dataset
5. add nwts2ph dataset, a hypothetical two-phase sampling based on nwts

addhazard 1.1.0 updates:

1. add a two-phase sampling example dataset, nwts.2ph.rda, to the data folder
2. add R file nwts.2ph.generator.R to generate the nwts.2ph example dataset
3. update documentations
4. add rownames to the summary stat where there is only one covariate
5. clean the typo in ah.2ph() on updating wts.2ph
6. add importFrom stats model.weight to ah.R file
7. add ah.fit.R but currently surppress the use of this function. 
   It will be activated in the later version. 
