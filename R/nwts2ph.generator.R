#######################################################
# Author: Kate HU
# Date: March 13th, 2017
# Task: Generate nwts2ph example dataset
#######################################################
#' This file genreate the example dataset nwts2ph
#' importFrom("stats", "rbinom")
#' @param data the full cohort data
#' @param seed the random seed we use for generating this dataset
#' @importFrom stats rbinom
#' @export nwts2ph.generate
nwts2ph.generate <- function(data, seed = 20){
	# Create a hypothetical two-phase sampling (stratified sampling) dataset.
	# In this dataset, only people who have relapse and some of the controls have
	# their samples sent to the central laboratory for histology determination
	set.seed(seed)
  nwts2ph <- data

	# Institutional histology is examined in the local hospital thus
	# it is reasonable to assume it is measured for all the samples.
	# Central histology is more expensive to obtain since the samples have to be
	# sent to the central laboratory and the work requires experienced lab scientists
	# thus  it is reasonable to assume only some samples were tested for central
	# histology. Based on the outcomes and institutional histology, we determine who
	# will be selected for central histology measurement, which follows two-phase sampling*
	# framework.

	# create strata based on institutional histology and disease status
	nwts2ph$strt <- 1 + nwts2ph$instit ### add a stratum containing all (relapsed) cases
	nwts2ph$strt[nwts2ph$relaps == 1] <- 3

	# phase II samples are obtained by Bernoulli Sampling
	# we oversample unfavorable institutional histology (instit = 1) and cases
	# selection probability Pi = 0.5 for instit = 0, Pi = 0.9 for instit = 1 and relaps = 1
	# assign phase II subsampling (Bernoulli) probabilities
	nwts2ph$Pi <- 0.5 * (nwts2ph$strt == 1) +
	              0.9 * (nwts2ph$strt == 2) +
	              0.9 * (nwts2ph$strt == 3)

	# generate phase II sampling indicators
	N <- dim(nwts2ph)[1]
	nwts2ph$in.ph2 <- rbinom(N, 1, nwts2ph$Pi)

	# update the central histology variable. Those who were not selected
	# into phase II subsamples are assumed not have central histology measurement
	nwts2ph$histol[nwts2ph$in.ph2 == 0] <- NA
	return(nwts2ph)
}

# generate the data
# nwts2ph <- nwts2ph.generate(nwtsco)
# save the data
#  devtools::use_data(nwts2ph)
