#######################################################
# Author: Kate HU
# Date: July 27th, 2016
# Task: Generate nwts.2ph example dataset
#######################################################


nwts.2ph.generate <- function(nwtsco){
	# Create a hypothetical two-phase sampling (stratified sampling) dataset.
	# In this dataset, only people who have relapse and some of the controls have
	# their samples sent to the central laboratory for histology determination
	nwts.2ph <- nwtsco

	# Institutional histology is examined in the local hospital thus
	# it is reasonable to assume it is measured for all the samples.
	# Central histology is more expensive to obtain since the samples have to be
	# sent to the central laboratory and the work requires experienced lab scientists
	# thus  it is reasonable to assume only some samples were tested for central
	# histology. Based on the outcomes and institutional histology, we determine who
	# will be selected for central histology measurement, which follows two-phase sampling*
	# framework.

	# create strata based on institutional histology and disease status
	nwts.2ph$strt <- 1 + nwts.2ph$instit ### add a stratum containing all (relapsed) cases
	nwts.2ph$strt[nwts.2ph$relaps == 1] <- 3

	# phase II samples are obtained by Bernoulli Sampling
	# we oversample unfavorable institutional histology (instit = 1) and cases
	# selection probability Pi = 0.5 for instit = 0, Pi = 0.9 for instit = 1 and relaps = 1
	# assign phase II subsampling (Bernoulli) probabilities
	nwts.2ph$Pi <- 0.5 * (nwts.2ph$strt == 1) + 0.9 * (nwts.2ph$strt == 2) + 0.9 * (nwts.2ph$strt == 3)

	# generate phase II sampling indicators
	N <- dim(nwts.2ph)[1]
	nwts.2ph$in.ph2 <- rbinom(N, 1, nwts$Pi)

	# update the central histology variable. Those who were not selected
	# into phase II subsamples are assumed not have central histology measurement
	nwts.2ph$histol[nwts.2ph$in.ph2 == 0] <- NA
	return(nwts.2ph)
}

# generate the data 
# nwts.2ph <- nwts.2ph.generate(nwtsco)
# save the data
# save(nwts.2ph, file='nwts.2ph.rdata')
