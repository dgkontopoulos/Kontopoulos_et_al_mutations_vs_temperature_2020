#!/usr/bin/Rscript

# This script fits phylogenetic and non-phylogenetic regressions of 
# per-generation substitution rate against temperature, cell volume, or 
# only an intercept, using MCMCglmm.
#
# 2 MCMC chains are run per model. The script checks that the chains have 
# converged on equivalent posterior distributions and that the 
# parameters are sufficiently sampled.
#
# The quality of the fitted models is evaluated according to the 
# Deviance Information Criterion (DIC). Marginal and conditional 
# R-squared values are also calculated.

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

# This function calculates the marginal and conditional R-squared values 
# following Nakagawa & Schielzeth, Methods Ecol. Evol., 2013.
calc_r_squared <- function(model, varF)
{
	
	# Get the variance explained by random effects.
	varRandom <- mean(model$VCV[,'Species'])
	
	# Get the residual (unexplained) variance.
	varResid <- mean(model$VCV[,'units'])
	
	# Calculate the marginal and conditional R-squared values. 
	#
	# Marginal R-squared: fraction of variance explained by fixed effects only.
	# Conditional R-squared: fraction of variance explained by fixed and random effects.
	r_sq_marginal <- varF/(varF + varRandom + varResid)
	r_sq_conditional <- (varF + varRandom)/(varF + varRandom + varResid)
	
	return(list(marginal = r_sq_marginal, conditional = r_sq_conditional))
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

# Read the dataset of per-generation substitution rates, temperature, 
# and maximum/minimum cell volume estimates.
dataset <- read.csv('ln_K_t_gen_vs_Temp_and_Cell_Volume.csv', stringsAsFactors = FALSE)
row.names(dataset) <- dataset$Species

sim_volumes <- list()
resulting_fits <- list()

# Fit the models 20 times.
for ( i in 1:20 )
{
	cat("Now at ", i, " ...\n", sep = '')
	
	temp_dat <- dataset
	
	# Remove the 1 species for which there are no cell volume estimates.
	temp_tree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% row.names(temp_dat))])
	
	temp_dat$Volume <- rep(NA, nrow(temp_dat))
	
	# For each species, sample a cell volume value from a uniform distribution 
	# with bounds equal to the minimum and maximum cell volume estimates.
	set.seed(i)
	for ( j in 1:nrow(temp_dat) )
	{
		temp_dat$Volume[j] <- runif(1, temp_dat$MinVolume[j], temp_dat$MaxVolume[j])
	}
	
	# Fit the various models using the sampled cell volume values (when 
	# applicable).	
	set.seed(1)	
	fit_1a <- MCMCglmm(ln_K_t_gen ~ log(Volume) + log(Temperature), family = "gaussian",
		prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
		random=~Species,
		data = temp_dat,
		nitt = 1000000,
		burnin = 100000,
		thin = 1000,
		ginverse=list(Species=inverseA(temp_tree, nodes = 'ALL', scale = TRUE)$Ainv),
		verbose = TRUE
	)
	set.seed(2)
	fit_1b <- MCMCglmm(ln_K_t_gen ~ log(Volume) + log(Temperature), family = "gaussian",
		prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
		random=~Species,
		data = temp_dat,
		nitt = 1000000,
		burnin = 100000,
		thin = 1000,
		ginverse=list(Species=inverseA(temp_tree, nodes = 'ALL', scale = TRUE)$Ainv),
		verbose = TRUE
	)
	check_convergence(fit_1a, fit_1b)
	check_ESS(fit_1a, fit_1b)
	
	set.seed(3)
	fit_2a <- MCMCglmm(ln_K_t_gen ~ log(Volume) + log(Temperature), family = "gaussian",
		prior=list(R=list(V=diag(1),nu=1.002)),
		data = temp_dat,
		nitt = 1000000,
		burnin = 100000,
		thin = 1000,
		verbose = TRUE
	)
	set.seed(4)
	fit_2b <- MCMCglmm(ln_K_t_gen ~ log(Volume) + log(Temperature), family = "gaussian",
		prior=list(R=list(V=diag(1),nu=1.002)),
		data = temp_dat,
		nitt = 1000000,
		burnin = 100000,
		thin = 1000,
		verbose = TRUE
	)
	check_convergence(fit_2a, fit_2b)
	check_ESS(fit_2a, fit_2b)
	
	set.seed(5)	
	fit_3a <- MCMCglmm(ln_K_t_gen ~ log(Temperature), family = "gaussian",
		prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
		random=~Species,
		data = temp_dat,
		nitt = 1000000,
		burnin = 100000,
		thin = 1000,
		ginverse=list(Species=inverseA(temp_tree, nodes = 'ALL', scale = TRUE)$Ainv),
		verbose = TRUE
	)
	set.seed(6)
	fit_3b <- MCMCglmm(ln_K_t_gen ~ log(Temperature), family = "gaussian",
		prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
		random=~Species,
		data = temp_dat,
		nitt = 1000000,
		burnin = 100000,
		thin = 1000,
		ginverse=list(Species=inverseA(temp_tree, nodes = 'ALL', scale = TRUE)$Ainv),
		verbose = TRUE
	)
	check_convergence(fit_3a, fit_3b)
	check_ESS(fit_3a, fit_3b)
	
	set.seed(7)
	fit_4a <- MCMCglmm(ln_K_t_gen ~ log(Temperature), family = "gaussian",
		prior=list(R=list(V=diag(1),nu=1.002)),
		data = temp_dat,
		nitt = 1000000,
		burnin = 100000,
		thin = 1000,
		verbose = TRUE
	)
	set.seed(8)
	fit_4b <- MCMCglmm(ln_K_t_gen ~ log(Temperature), family = "gaussian",
		prior=list(R=list(V=diag(1),nu=1.002)),
		data = temp_dat,
		nitt = 1000000,
		burnin = 100000,
		thin = 1000,
		verbose = TRUE
	)
	check_convergence(fit_4a, fit_4b)
	check_ESS(fit_4a, fit_4b)
	
	set.seed(9)	
	fit_5a <- MCMCglmm(ln_K_t_gen ~ log(Volume), family = "gaussian",
		prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
		random=~Species,
		data = temp_dat,
		nitt = 1000000,
		burnin = 100000,
		thin = 1000,
		ginverse=list(Species=inverseA(temp_tree, nodes = 'ALL', scale = TRUE)$Ainv),
		verbose = TRUE
	)
	set.seed(10)
	fit_5b <- MCMCglmm(ln_K_t_gen ~ log(Volume), family = "gaussian",
		prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
		random=~Species,
		data = temp_dat,
		nitt = 1000000,
		burnin = 100000,
		thin = 1000,
		ginverse=list(Species=inverseA(temp_tree, nodes = 'ALL', scale = TRUE)$Ainv),
		verbose = TRUE
	)
	check_convergence(fit_5a, fit_5b)
	check_ESS(fit_5a, fit_5b)
	
	set.seed(11)
	fit_6a <- MCMCglmm(ln_K_t_gen ~ log(Volume), family = "gaussian",
		prior=list(R=list(V=diag(1),nu=1.002)),
		data = temp_dat,
		nitt = 1000000,
		burnin = 100000,
		thin = 1000,
		verbose = TRUE
	)
	set.seed(12)
	fit_6b <- MCMCglmm(ln_K_t_gen ~ log(Volume), family = "gaussian",
		prior=list(R=list(V=diag(1),nu=1.002)),
		data = temp_dat,
		nitt = 1000000,
		burnin = 100000,
		thin = 1000,
		verbose = TRUE
	)
	check_convergence(fit_6a, fit_6b)
	check_ESS(fit_6a, fit_6b)
	
	set.seed(13)	
	fit_7a <- MCMCglmm(ln_K_t_gen ~ 1, family = "gaussian",
		prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
		random=~Species,
		data = temp_dat,
		nitt = 1000000,
		burnin = 100000,
		thin = 1000,
		ginverse=list(Species=inverseA(temp_tree, nodes = 'ALL', scale = TRUE)$Ainv),
		verbose = TRUE
	)
	set.seed(14)
	fit_7b <- MCMCglmm(ln_K_t_gen ~ 1, family = "gaussian",
		prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
		random=~Species,
		data = temp_dat,
		nitt = 1000000,
		burnin = 100000,
		thin = 1000,
		ginverse=list(Species=inverseA(temp_tree, nodes = 'ALL', scale = TRUE)$Ainv),
		verbose = TRUE
	)
	check_convergence(fit_7a, fit_7b)
	check_ESS(fit_7a, fit_7b)
	
	set.seed(15)
	fit_8a <- MCMCglmm(ln_K_t_gen ~ 1, family = "gaussian",
		prior=list(R=list(V=diag(1),nu=1.002)),
		data = temp_dat,
		nitt = 1000000,
		burnin = 100000,
		thin = 1000,
		verbose = TRUE
	)
	set.seed(16)
	fit_8b <- MCMCglmm(ln_K_t_gen ~ 1, family = "gaussian",
		prior=list(R=list(V=diag(1),nu=1.002)),
		data = temp_dat,
		nitt = 1000000,
		burnin = 100000,
		thin = 1000,
		verbose = TRUE
	)
	check_convergence(fit_8a, fit_8b)
	check_ESS(fit_8a, fit_8b)
	
	set.seed(17)	
	fit_9a <- MCMCglmm(ln_K_t_gen ~ log(Volume) + poly(Temperature, 2, raw = TRUE), family = "gaussian",
		prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
		random=~Species,
		data = temp_dat,
		nitt = 1000000,
		burnin = 100000,
		thin = 1000,
		ginverse=list(Species=inverseA(temp_tree, nodes = 'ALL', scale = TRUE)$Ainv),
		verbose = TRUE
	)
	set.seed(18)
	fit_9b <- MCMCglmm(ln_K_t_gen ~ log(Volume) + poly(Temperature, 2, raw = TRUE), family = "gaussian",
		prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
		random=~Species,
		data = temp_dat,
		nitt = 1000000,
		burnin = 100000,
		thin = 1000,
		ginverse=list(Species=inverseA(temp_tree, nodes = 'ALL', scale = TRUE)$Ainv),
		verbose = TRUE
	)
	check_convergence(fit_9a, fit_9b)
	check_ESS(fit_9a, fit_9b)
	
	set.seed(19)
	fit_10a <- MCMCglmm(ln_K_t_gen ~ log(Volume) + poly(Temperature, 2, raw = TRUE), family = "gaussian",
		prior=list(R=list(V=diag(1),nu=1.002)),
		data = temp_dat,
		nitt = 1000000,
		burnin = 100000,
		thin = 1000,
		verbose = TRUE
	)
	set.seed(20)
	fit_10b <- MCMCglmm(ln_K_t_gen ~ log(Volume) + poly(Temperature, 2, raw = TRUE), family = "gaussian",
		prior=list(R=list(V=diag(1),nu=1.002)),
		data = temp_dat,
		nitt = 1000000,
		burnin = 100000,
		thin = 1000,
		verbose = TRUE
	)
	check_convergence(fit_10a, fit_10b)
	check_ESS(fit_10a, fit_10b)
	
	set.seed(21)	
	fit_11a <- MCMCglmm(ln_K_t_gen ~ poly(Temperature, 2, raw = TRUE), family = "gaussian",
		prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
		random=~Species,
		data = temp_dat,
		nitt = 1000000,
		burnin = 100000,
		thin = 1000,
		ginverse=list(Species=inverseA(temp_tree, nodes = 'ALL', scale = TRUE)$Ainv),
		verbose = TRUE
	)
	set.seed(22)
	fit_11b <- MCMCglmm(ln_K_t_gen ~ poly(Temperature, 2, raw = TRUE), family = "gaussian",
		prior=list(G=list(G1=list(V=diag(1),nu=1.002)), R=list(V=diag(1),nu=1.002)),
		random=~Species,
		data = temp_dat,
		nitt = 1000000,
		burnin = 100000,
		thin = 1000,
		ginverse=list(Species=inverseA(temp_tree, nodes = 'ALL', scale = TRUE)$Ainv),
		verbose = TRUE
	)
	check_convergence(fit_11a, fit_11b)
	check_ESS(fit_11a, fit_11b)
	
	set.seed(23)
	fit_12a <- MCMCglmm(ln_K_t_gen ~ poly(Temperature, 2, raw = TRUE), family = "gaussian",
		prior=list(R=list(V=diag(1),nu=1.002)),
		data = temp_dat,
		nitt = 1000000,
		burnin = 100000,
		thin = 1000,
		verbose = TRUE
	)
	set.seed(24)
	fit_12b <- MCMCglmm(ln_K_t_gen ~ poly(Temperature, 2, raw = TRUE), family = "gaussian",
		prior=list(R=list(V=diag(1),nu=1.002)),
		data = temp_dat,
		nitt = 1000000,
		burnin = 100000,
		thin = 1000,
		verbose = TRUE
	)
	check_convergence(fit_12a, fit_12b)
	check_ESS(fit_12a, fit_12b)
	
	# Save the sampled cell volume values and the model fits.
	sim_volumes[[i]] <- temp_dat$Volume
	resulting_fits[[i]] <- list(
		fit_1a,		fit_1b,		fit_2a,		fit_2b,		fit_3a,		fit_3b, 
		fit_4a,		fit_4b,		fit_5a,		fit_5b,		fit_6a,		fit_6b,
		fit_7a,		fit_7b,		fit_8a,		fit_8b,		fit_9a,		fit_9b,
		fit_10a,	fit_10b, 	fit_11a,	fit_11b,	fit_12a,	fit_12b
	)
}

# Create a table to store the DIC values of fitted models.
dic_table <- data.frame(
	ID = rep(NA, 20), 
	DIC_1 = rep(NA, 20), DIC_2 = rep(NA, 20), DIC_3 = rep(NA, 20), 
	DIC_4 = rep(NA, 20), DIC_5 = rep(NA, 20), DIC_6 = rep(NA, 20), 
	DIC_7 = rep(NA, 20), DIC_8 = rep(NA, 20), DIC_9 = rep(NA, 20), 
	DIC_10 = rep(NA, 20), DIC_11 = rep(NA, 20), DIC_12 = rep(NA, 20)
)

# Populate the DIC table.
for ( i in 1:20 )
{
	dic_table[i,1] <- i
	dic_table[i,2] <- mean(
		c(
			resulting_fits[[i]][[1]]$DIC, resulting_fits[[i]][[2]]$DIC
		)
	)
	dic_table[i,3] <- mean(
		c(

			resulting_fits[[i]][[3]]$DIC, resulting_fits[[i]][[4]]$DIC
		)
	)
	dic_table[i,4] <- mean(
		c(
			resulting_fits[[i]][[5]]$DIC, resulting_fits[[i]][[6]]$DIC
		)
	)
	dic_table[i,5] <- mean(
		c(
			resulting_fits[[i]][[7]]$DIC, resulting_fits[[i]][[8]]$DIC
		)
	)
	dic_table[i,6] <- mean(
		c(
			resulting_fits[[i]][[9]]$DIC, resulting_fits[[i]][[10]]$DIC
		)
	)
	dic_table[i,7] <- mean(
		c(
			resulting_fits[[i]][[11]]$DIC, resulting_fits[[i]][[12]]$DIC
		)
	)
	dic_table[i,8] <- mean(
		c(
			resulting_fits[[i]][[13]]$DIC, resulting_fits[[i]][[14]]$DIC
		)
	)
	dic_table[i,9] <- mean(
		c(
			resulting_fits[[i]][[15]]$DIC, resulting_fits[[i]][[16]]$DIC
		)
	)
	dic_table[i,10] <- mean(
		c(
			resulting_fits[[i]][[17]]$DIC, resulting_fits[[i]][[18]]$DIC
		)
	)
	dic_table[i,11] <- mean(
		c(
			resulting_fits[[i]][[19]]$DIC, resulting_fits[[i]][[20]]$DIC
		)
	)
	dic_table[i,12] <- mean(
		c(
			resulting_fits[[i]][[21]]$DIC, resulting_fits[[i]][[22]]$DIC
		)
	)
	dic_table[i,13] <- mean(
		c(
			resulting_fits[[i]][[23]]$DIC, resulting_fits[[i]][[24]]$DIC
		)
	)
}

# The best-fitting model is the #11, i.e., with a fixed effect of 
# temperature as a second-order polynomial and with a phylogenetic 
# random effect on the intercept.

# Merge the estimates of the two chains belonging to the 11th model.
joined_best_fit <- list(
	Sol = rbind(fit_11a$Sol, fit_11b$Sol),
	VCV = rbind(fit_11a$VCV, fit_11b$VCV)
)

# Calculate the marginal and conditional R-squared values.
r_squared_vals <- calc_r_squared(
	joined_best_fit, mean(c(calc_varF(fit_11a), calc_varF(fit_11b)))
)
