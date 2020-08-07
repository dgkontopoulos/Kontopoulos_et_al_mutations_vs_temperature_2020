#!/usr/bin/Rscript

# This script uses weighted least squares regression to estimate the 
# relationship between median ΔΔG and temperature for the ADKs of 
# 5 extant species and for that of the Last Universal Common Ancestor.
#
# AICc is used to identify the model that best represents each 
# relationship.

library(MuMIn)

#####################
# F U N C T I O N S #
#####################

# This function fits six different models to the relationship between 
# median ΔΔG and temperature. The best model is determined on the 
# basic of AICc.
fit_models <- function(current_dataset)
{
	lm_fit_empty <- lm(
		ddG_med ~ 1, data = current_dataset, weights = nframes
	)
	
	lm_fit_lin <- lm(
		ddG_med ~ Temperature, data = current_dataset, weights = nframes
	)
		
	lm_fit_log <- lm(
		ddG_med ~ log(Temperature), data = current_dataset, weights = nframes
	)
		
	lm_fit_exp <- lm(
		ddG_med ~ exp(Temperature), data = current_dataset, weights = nframes
	)
	
	lm_fit_sqrt <- lm(
		ddG_med ~ sqrt(Temperature), data = current_dataset, weights = nframes
	)
	
	lm_fit_squared <- lm(
		ddG_med ~ I(Temperature^2), data = current_dataset, weights = nframes
	)
	
	models <- list(
		lm_fit_empty,   lm_fit_lin,  lm_fit_log, 
		lm_fit_exp,     lm_fit_sqrt, lm_fit_squared
	)
	
	AICc_vals <- rep(NA, 6) 
	for ( i in 1:6 )
	{
		AICc_vals[i] <- AICc(models[[i]])
	}
	
	id_best <- which(AICc_vals == min(AICc_vals))
	
	return(models[[id_best]])
}

####################
# M A I N  C O D E #
####################

# Read the intraspecific dataset of ΔΔG vs temperature.
dataset <- read.csv("ddG_vs_Temp_intraspecific.csv", stringsAsFactors = FALSE)

# For each species, fit the models and return the most appropriate one.
best_fit_per_species <- list()

for ( current_species in unique(dataset$Species) )
{
	current_dataset <- dataset[dataset$Species == current_species,]
	
	best_fit_per_species[[current_species]] <- fit_models(current_dataset)	
}
