#!/usr/bin/Rscript

# This script fits phylogenetic and non-phylogenetic regressions of 
# the weighted median ΔΔG value against temperature, LID length, ADK 
# type, or only an intercept, using MCMCglmm.
#
# 2 MCMC chains are run per model. The script checks that the chains have 
# converged on equivalent posterior distributions and that the 
# parameters are sufficiently sampled.
#
# The script also the calculates the R-squared value of the best 
# fitting model according to DIC.

library(MCMCglmm)

#####################
# F U N C T I O N S #
#####################

# This function returns the variance explained by the fixed effects in 
# an MCMCglmm.
calc_varF <- function(model)
{
	return(var(as.vector(apply(model$Sol,2,mean) %*% t(model$X))))
}

# This function calculates the R-squared value of an MCMCglmm without 
# random effects.
calc_r_squared_NO_RANDOM_EFFECTS <- function(model, varF)
{
	varResid <- mean(model$VCV[,'units'])
	return(varF/(varF + varResid))
}

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

#############################
# F I T T E D   M O D E L S #
#############################

# 1. Intercept-only

set.seed(1)
fit_1a <- MCMCglmm(ddG_w_med ~ 1, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(2)
fit_1b <- MCMCglmm(ddG_w_med ~ 1, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_1a, fit_1b)
check_convergence(fit_1a, fit_1b)

# 2. ~ Temperature

set.seed(3)
fit_2a <- MCMCglmm(ddG_w_med ~ Temperature, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(4)
fit_2b <- MCMCglmm(ddG_w_med ~ Temperature, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_2a, fit_2b)
check_convergence(fit_2a, fit_2b)

# 3. ~ LID_length

set.seed(5)
fit_3a <- MCMCglmm(ddG_w_med ~ LID_length, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(6)
fit_3b <- MCMCglmm(ddG_w_med ~ LID_length, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_3a, fit_3b)
check_convergence(fit_3a, fit_3b)

# 4. ~ Extra_beta_strands

set.seed(7)
fit_4a <- MCMCglmm(ddG_w_med ~ Extra_beta_strands, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(8)
fit_4b <- MCMCglmm(ddG_w_med ~ Extra_beta_strands, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_4a, fit_4b)
check_convergence(fit_4a, fit_4b)

# 5. ~ Temperature + LID_length

set.seed(9)
fit_5a <- MCMCglmm(ddG_w_med ~ Temperature + LID_length, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(10)
fit_5b <- MCMCglmm(ddG_w_med ~ Temperature + LID_length, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_5a, fit_5b)
check_convergence(fit_5a, fit_5b)

# 6. ~ Temperature + Extra_beta_strands

set.seed(11)
fit_6a <- MCMCglmm(ddG_w_med ~ Temperature + Extra_beta_strands, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(12)
fit_6b <- MCMCglmm(ddG_w_med ~ Temperature + Extra_beta_strands, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_6a, fit_6b)
check_convergence(fit_6a, fit_6b)

# 7. ~ LID_length + Extra_beta_strands

set.seed(13)
fit_7a <- MCMCglmm(ddG_w_med ~ LID_length + Extra_beta_strands, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(14)
fit_7b <- MCMCglmm(ddG_w_med ~ LID_length + Extra_beta_strands, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_7a, fit_7b)
check_convergence(fit_7a, fit_7b)

# 8. ~ Temperature + LID_length + Extra_beta_strands

set.seed(15)
fit_8a <- MCMCglmm(ddG_w_med ~ Temperature + LID_length + Extra_beta_strands, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(16)
fit_8b <- MCMCglmm(ddG_w_med ~ Temperature + LID_length + Extra_beta_strands, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_8a, fit_8b)
check_convergence(fit_8a, fit_8b)

# 9. ~ Temperature * LID_length + Extra_beta_strands

set.seed(17)
fit_9a <- MCMCglmm(ddG_w_med ~ Temperature * LID_length + Extra_beta_strands, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(18)
fit_9b <- MCMCglmm(ddG_w_med ~ Temperature * LID_length + Extra_beta_strands, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_9a, fit_9b)
check_convergence(fit_9a, fit_9b)

# 10. ~ Temperature * LID_length

set.seed(19)
fit_10a <- MCMCglmm(ddG_w_med ~ Temperature * LID_length, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(20)
fit_10b <- MCMCglmm(ddG_w_med ~ Temperature * LID_length, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_10a, fit_10b)
check_convergence(fit_10a, fit_10b)

# 11. ~ Temperature * Extra_beta_strands + LID_length

set.seed(21)
fit_11a <- MCMCglmm(ddG_w_med ~ Temperature * Extra_beta_strands + LID_length, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(22)
fit_11b <- MCMCglmm(ddG_w_med ~ Temperature * Extra_beta_strands + LID_length, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_11a, fit_11b)
check_convergence(fit_11a, fit_11b)

# 12. ~ Temperature * Extra_beta_strands

set.seed(23)
fit_12a <- MCMCglmm(ddG_w_med ~ Temperature * Extra_beta_strands, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(24)
fit_12b <- MCMCglmm(ddG_w_med ~ Temperature * Extra_beta_strands, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_12a, fit_12b)
check_convergence(fit_12a, fit_12b)

# 13. ~ Temperature * LID_length * Extra_beta_strands

set.seed(25)
fit_13a <- MCMCglmm(ddG_w_med ~ Temperature * LID_length * Extra_beta_strands, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(26)
fit_13b <- MCMCglmm(ddG_w_med ~ Temperature * LID_length * Extra_beta_strands, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_13a, fit_13b)
check_convergence(fit_13a, fit_13b)

# 14. ~ LID_length * Extra_beta_strands + Temperature

set.seed(27)
fit_14a <- MCMCglmm(ddG_w_med ~ LID_length * Extra_beta_strands + Temperature, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(28)
fit_14b <- MCMCglmm(ddG_w_med ~ LID_length * Extra_beta_strands + Temperature, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_14a, fit_14b)
check_convergence(fit_14a, fit_14b)

# 15. ~ LID_length * Extra_beta_strands

set.seed(29)
fit_15a <- MCMCglmm(ddG_w_med ~ LID_length * Extra_beta_strands, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(30)
fit_15b <- MCMCglmm(ddG_w_med ~ LID_length * Extra_beta_strands, family = "gaussian",
    prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
    random=~Species,
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    ginverse=list(Species=inverseA(tree, nodes = 'ALL', scale = TRUE)$Ainv),
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_15a, fit_15b)
check_convergence(fit_15a, fit_15b)

# WITHOUT PHYLOGENETIC CORRECTION

# 1. Intercept-only

set.seed(31)
fit_NON_PHYLO_1a <- MCMCglmm(ddG_w_med ~ 1, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(32)
fit_NON_PHYLO_1b <- MCMCglmm(ddG_w_med ~ 1, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_NON_PHYLO_1a, fit_NON_PHYLO_1b)
check_convergence(fit_NON_PHYLO_1a, fit_NON_PHYLO_1b)

# 2. ~ Temperature

set.seed(33)
fit_NON_PHYLO_2a <- MCMCglmm(ddG_w_med ~ Temperature, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(34)
fit_NON_PHYLO_2b <- MCMCglmm(ddG_w_med ~ Temperature, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_NON_PHYLO_2a, fit_NON_PHYLO_2b)
check_convergence(fit_NON_PHYLO_2a, fit_NON_PHYLO_2b)

# 3. ~ LID_length

set.seed(35)
fit_NON_PHYLO_3a <- MCMCglmm(ddG_w_med ~ LID_length, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(36)
fit_NON_PHYLO_3b <- MCMCglmm(ddG_w_med ~ LID_length, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_NON_PHYLO_3a, fit_NON_PHYLO_3b)
check_convergence(fit_NON_PHYLO_3a, fit_NON_PHYLO_3b)

# 4. ~ Extra_beta_strands

set.seed(37)
fit_NON_PHYLO_4a <- MCMCglmm(ddG_w_med ~ Extra_beta_strands, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(38)
fit_NON_PHYLO_4b <- MCMCglmm(ddG_w_med ~ Extra_beta_strands, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_NON_PHYLO_4a, fit_NON_PHYLO_4b)
check_convergence(fit_NON_PHYLO_4a, fit_NON_PHYLO_4b)

# 5. ~ Temperature + LID_length

set.seed(39)
fit_NON_PHYLO_5a <- MCMCglmm(ddG_w_med ~ Temperature + LID_length, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(40)
fit_NON_PHYLO_5b <- MCMCglmm(ddG_w_med ~ Temperature + LID_length, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_NON_PHYLO_5a, fit_NON_PHYLO_5b)
check_convergence(fit_NON_PHYLO_5a, fit_NON_PHYLO_5b)

# 6. ~ Temperature + Extra_beta_strands

set.seed(41)
fit_NON_PHYLO_6a <- MCMCglmm(ddG_w_med ~ Temperature + Extra_beta_strands, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(42)
fit_NON_PHYLO_6b <- MCMCglmm(ddG_w_med ~ Temperature + Extra_beta_strands, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_NON_PHYLO_6a, fit_NON_PHYLO_6b)
check_convergence(fit_NON_PHYLO_6a, fit_NON_PHYLO_6b)

# 7. ~ LID_length + Extra_beta_strands

set.seed(43)
fit_NON_PHYLO_7a <- MCMCglmm(ddG_w_med ~ LID_length + Extra_beta_strands, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(44)
fit_NON_PHYLO_7b <- MCMCglmm(ddG_w_med ~ LID_length + Extra_beta_strands, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_NON_PHYLO_7a, fit_NON_PHYLO_7b)
check_convergence(fit_NON_PHYLO_7a, fit_NON_PHYLO_7b)

# 8. ~ Temperature + LID_length + Extra_beta_strands

set.seed(45)
fit_NON_PHYLO_8a <- MCMCglmm(ddG_w_med ~ Temperature + LID_length + Extra_beta_strands, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(46)
fit_NON_PHYLO_8b <- MCMCglmm(ddG_w_med ~ Temperature + LID_length + Extra_beta_strands, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_NON_PHYLO_8a, fit_NON_PHYLO_8b)
check_convergence(fit_NON_PHYLO_8a, fit_NON_PHYLO_8b)

# 9. ~ Temperature * LID_length + Extra_beta_strands

set.seed(47)
fit_NON_PHYLO_9a <- MCMCglmm(ddG_w_med ~ Temperature * LID_length + Extra_beta_strands, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(48)
fit_NON_PHYLO_9b <- MCMCglmm(ddG_w_med ~ Temperature * LID_length + Extra_beta_strands, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_NON_PHYLO_9a, fit_NON_PHYLO_9b)
check_convergence(fit_NON_PHYLO_9a, fit_NON_PHYLO_9b)

# 10. ~ Temperature * LID_length

set.seed(49)
fit_NON_PHYLO_10a <- MCMCglmm(ddG_w_med ~ Temperature * LID_length, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(50)
fit_NON_PHYLO_10b <- MCMCglmm(ddG_w_med ~ Temperature * LID_length, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_NON_PHYLO_10a, fit_NON_PHYLO_10b)
check_convergence(fit_NON_PHYLO_10a, fit_NON_PHYLO_10b)

# 11. ~ Temperature * Extra_beta_strands + LID_length

set.seed(51)
fit_NON_PHYLO_11a <- MCMCglmm(ddG_w_med ~ Temperature * Extra_beta_strands + LID_length, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(52)
fit_NON_PHYLO_11b <- MCMCglmm(ddG_w_med ~ Temperature * Extra_beta_strands + LID_length, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_NON_PHYLO_11a, fit_NON_PHYLO_11b)
check_convergence(fit_NON_PHYLO_11a, fit_NON_PHYLO_11b)

# 12. ~ Temperature * Extra_beta_strands

set.seed(53)
fit_NON_PHYLO_12a <- MCMCglmm(ddG_w_med ~ Temperature * Extra_beta_strands, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(54)
fit_NON_PHYLO_12b <- MCMCglmm(ddG_w_med ~ Temperature * Extra_beta_strands, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_NON_PHYLO_12a, fit_NON_PHYLO_12b)
check_convergence(fit_NON_PHYLO_12a, fit_NON_PHYLO_12b)

# 13. ~ Temperature * LID_length * Extra_beta_strands

set.seed(55)
fit_NON_PHYLO_13a <- MCMCglmm(ddG_w_med ~ Temperature * LID_length * Extra_beta_strands, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(56)
fit_NON_PHYLO_13b <- MCMCglmm(ddG_w_med ~ Temperature * LID_length * Extra_beta_strands, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_NON_PHYLO_13a, fit_NON_PHYLO_13b)
check_convergence(fit_NON_PHYLO_13a, fit_NON_PHYLO_13b)

# 14. ~ LID_length * Extra_beta_strands + Temperature

set.seed(57)
fit_NON_PHYLO_14a <- MCMCglmm(ddG_w_med ~ LID_length * Extra_beta_strands + Temperature, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(58)
fit_NON_PHYLO_14b <- MCMCglmm(ddG_w_med ~ LID_length * Extra_beta_strands + Temperature, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_NON_PHYLO_14a, fit_NON_PHYLO_14b)
check_convergence(fit_NON_PHYLO_14a, fit_NON_PHYLO_14b)

# 15. ~ LID_length * Extra_beta_strands

set.seed(59)
fit_NON_PHYLO_15a <- MCMCglmm(ddG_w_med ~ LID_length * Extra_beta_strands, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

set.seed(60)
fit_NON_PHYLO_15b <- MCMCglmm(ddG_w_med ~ LID_length * Extra_beta_strands, family = "gaussian",
    prior=list(R=list(V=diag(1),nu=1.002)),
    data = dataset,
    nitt = 3000000,
    burnin = 300000,
    thin = 1000,
    verbose = TRUE,
    mev = dataset$ddG_SEV
)

check_ESS(fit_NON_PHYLO_15a, fit_NON_PHYLO_15b)
check_convergence(fit_NON_PHYLO_15a, fit_NON_PHYLO_15b)

# Merge the estimates of the two chains belonging to the best-fitting 
# model (according to DIC, excluding models whose coefficients have 
# 95% HPD intervals that include zero).
merged_chains <- list(
	Sol = rbind(fit_NON_PHYLO_6a$Sol, fit_NON_PHYLO_6b$Sol),
	VCV = rbind(fit_NON_PHYLO_6a$VCV, fit_NON_PHYLO_6b$VCV)
)

# Calculate the R-squared value.
r_squared <- calc_r_squared_NO_RANDOM_EFFECTS(merged_chains, mean(c(calc_varF(fit_NON_PHYLO_6a), calc_varF(fit_NON_PHYLO_6b))))
