#!/usr/bin/Rscript

# This script estimates the phylogenetic heritability of the weighted 
# median ΔΔG value separately for ancestral- and archaeal-type ADKs.

library(MCMCglmm)

#####################
# F U N C T I O N S #
#####################

# This function checks if the two MCMCglmm chains converged on statistically 
# indistinguishable posterior distributions, based on the Potential 
# Scale Reduction Factor diagnostic (Gelman & Rubin, Stat. Sci., 1992).
check_convergence <- function(fit_1a, fit_1b)
{
	
	# Estimate the PSRF for each element of the fixed effects.
	# 
	# Values >= 1.1 indicate lack of convergence.
	# In that case, raise an error.
	gelman_Sol <- gelman.diag(mcmc.list(fit_1a$Sol, fit_1b$Sol), multivariate = FALSE)
	if ( length(gelman_Sol$psrf[
				!is.na(gelman_Sol$psrf[,1]) &
				gelman_Sol$psrf[,1] >= 1.1,
			]
		) > 0 
	)
	{
		stop("PROBLEM! The two runs have not converged!")
	}
	
	# Same thing for the random effects and the residual variance.
	gelman_VCV <- gelman.diag(mcmc.list(fit_1a$VCV, fit_1b$VCV), multivariate = FALSE)		
	if ( length(gelman_VCV$psrf[
				!is.na(gelman_VCV$psrf[,1]) &
				gelman_VCV$psrf[,1] >= 1.1,
			]
		) > 0 
	)
	{
		stop("PROBLEM! The two runs have not converged!")
	}
}

# This function checks that each parameter was adequately sampled from 
# the posterior, by calculating the Effective Sample Size.
check_ESS <- function(fit_1a, fit_1b)
{
	
	# Calculate the Effective Sample Size for each element of the 
	# fixed effects.
	merged_d_Sol <- rbind(fit_1a$Sol, fit_1b$Sol)
	ess_vals_Sol <- effectiveSize(merged_d_Sol)
	ess_vals_Sol <- ess_vals_Sol[ess_vals_Sol >= 1]
	
	# If there is a parameter with an Effective Sample Size below 
	# 200, raise an error.
	if ( length(ess_vals_Sol[ess_vals_Sol < 200]) > 0 )
	{
		stop("PROBLEM! The combined ESS is not big enough!")
	}

	# Same thing for the random effects and the residual variance.
	merged_d_VCV <- rbind(fit_1a$VCV, fit_1b$VCV)
	ess_vals_VCV <- effectiveSize(merged_d_VCV)
	ess_vals_VCV <- ess_vals_VCV[ess_vals_VCV >= 1]
	
	if ( length(ess_vals_VCV[ess_vals_VCV < 200]) > 0 )
	{
		stop("PROBLEM! The combined ESS is not big enough!")
	}	
}

####################
# M A I N  C O D E #
####################

# Read the time-calibrated phylogeny.
tree <- read.nexus('calibrated_species_tree.phy')

# Read the interspecific dataset of weighted median ΔΔG values.
dataset <- read.csv('ddG_vs_Temp_interspecific.csv', stringsAsFactors = FALSE)
row.names(dataset) <- dataset$Species

##################
# Ancestral ADKs #
##################

# Filter the dataset and phylogeny to only include ancestral ADKs.
ancestral_dat <- dataset[dataset$Extra_beta_strands == "NO",]
ancestral_tree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% row.names(ancestral_dat))])

# Fit 2 chains of the intercept-only model with MCMCglmm.
set.seed(1)
fit_anc_a <- MCMCglmm(ddG_w_med ~ 1, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = ancestral_dat,
    nitt = 10000000,
    burnin = 1000000,
    thin = 1000,
    ginverse=list(Species=inverseA(ancestral_tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = ancestral_dat$ddG_SEV
)

set.seed(2)
fit_anc_b <- MCMCglmm(ddG_w_med ~ 1, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = ancestral_dat,
    nitt = 10000000,
    burnin = 1000000,
    thin = 1000,
    ginverse=list(Species=inverseA(ancestral_tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = ancestral_dat$ddG_SEV
)

# Ensure that the two chains converged and that the Effective Sample 
# Size is high enough.
check_ESS(fit_anc_a, fit_anc_b)
check_convergence(fit_anc_a, fit_anc_b)

merged_chains_ancestral <- list(
	Sol = rbind(fit_anc_a$Sol, fit_anc_b$Sol),
	VCV = rbind(fit_anc_a$VCV, fit_anc_b$VCV)
)

# Calculate the phylogenetic heritability of ΔΔG for ancestral-type ADKs.
herit_ddG_anc <- merged_chains_ancestral$VCV[,'Species']/(merged_chains_ancestral$VCV[,'Species'] + merged_chains_ancestral$VCV[,'units'])

#################
# Archaeal ADKs #
#################

# Filter the dataset and phylogeny to only include archaeal ADKs.
archaeal_dat <- dataset[dataset$Extra_beta_strands == "YES",]
archaeal_tree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% row.names(archaeal_dat))])

# Fit 2 chains of the intercept-only model with MCMCglmm.
set.seed(1)
fit_arc_a <- MCMCglmm(ddG_w_med ~ 1, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = archaeal_dat,
    nitt = 10000000,
    burnin = 1000000,
    thin = 1000,
    ginverse=list(Species=inverseA(archaeal_tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = archaeal_dat$ddG_SEV
)

set.seed(2)
fit_arc_b <- MCMCglmm(ddG_w_med ~ 1, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = archaeal_dat,
    nitt = 10000000,
    burnin = 1000000,
    thin = 1000,
    ginverse=list(Species=inverseA(archaeal_tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = archaeal_dat$ddG_SEV
)

# Ensure that the two chains converged and that the Effective Sample 
# Size is high enough.
check_ESS(fit_arc_a, fit_arc_b)
check_convergence(fit_arc_a, fit_arc_b)

merged_chains_archaeal <- list(
	Sol = rbind(fit_arc_a$Sol, fit_arc_b$Sol),
	VCV = rbind(fit_arc_a$VCV, fit_arc_b$VCV)
)

# Calculate the phylogenetic heritability of ΔΔG for archaeal-type ADKs.
herit_ddG_arc <- merged_chains_archaeal$VCV[,'Species']/(merged_chains_archaeal$VCV[,'Species'] + merged_chains_archaeal$VCV[,'units'])
