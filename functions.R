

get_ntax <- function(t, lambda, cutoff=500) {
	if(2*exp(lambda*t) < cutoff) {
		return(ape::Ntip(geiger::sim.bdtree(b=lambda, d=0, stop="time", t=t)))
		#return(ape::Ntip(TreeSim::sim.bd.age(age=t, numbsim=1, lambda=lambda, mu=0, frac = 1)[[1]]))
	} else {
		max_species <- min(round(max(1e5, 100*2*exp(lambda*t))),1e7)
		probabilities <- get_m_s_prob_nt_all_i(max_species_to_consider=max_species,netdiv=lambda, ef=0, time=t)
		probabilities <- probabilities/sum(probabilities) #normalize to handle cases where max_species is too small
		return(sim_from_m_s_prob_nt_all_i(probabilities))
	}
}

get_cv <- function(x) {
	return(sd(x)/mean(x))
}

get_r <- function(x, y) {
	return(cor(x, y))
}

get_empirical_ratio <- function(x1, x2) {
	return(mean(x1)/mean(x2))	
}

get_expected_ratio <- function(x1, x2) {
	return((mean(x1)/mean(x2)) * (1 + get_cv(x2) * (get_cv(x2) - get_r(x1, x2) * get_cv(x1))))	
}

get_m_s_alpha <- function(netdiv, ef, time) {
	return(ef*get_m_s_beta(netdiv, ef, time))
}	

get_m_s_beta <- function(netdiv, ef, time) {
	return((exp(netdiv*time)-1)/(exp(netdiv*time) - ef))
}



# Magallon and Sanderson 2001 equation 1b
get_m_s_prob_nt_eq_i <- function(final_species_count, netdiv, ef, time, starting_species_count=2) {
	j <- 0
	prob <- 0
	alpha <- get_m_s_alpha(netdiv, ef, time)
	beta <- get_m_s_beta(netdiv, ef, time)
	new_prob <- alpha^starting_species_count
	if(final_species_count>0) {
		for (j in 0:min(starting_species_count, final_species_count)) {
			#new_prob <- chooseMpfr(starting_species_count, j) * chooseMpfr(final_species_count+starting_species_count-j-1, 1) * (alpha^(starting_species_count-j)) * (beta^(final_species_count-j)) * ((1 - alpha - beta)^j)
			new_prob <- choose(starting_species_count, j) * choose(final_species_count+starting_species_count-j-1, 1) * (alpha^(starting_species_count-j)) * (beta^(final_species_count-j)) * ((1 - alpha - beta)^j)

			prob <- prob + new_prob
		}
	} else {
		prob <- new_prob
	}
	return(prob)
}

get_m_s_prob_nt_all_i <- function(max_species_to_consider, netdiv, ef, time, starting_species_count=2) {
	all_counts <- sequence(max_species_to_consider)
	names(all_counts) <- all_counts
	probs <- sapply(all_counts, get_m_s_prob_nt_eq_i, netdiv=netdiv, ef=ef, time=time, starting_species_count=starting_species_count)	
	return(probs)
}

sim_from_m_s_prob_nt_all_i <- function(probabilities) {
	winner_index <- sample(1:length(probabilities), size=1, prob=probabilities)
	return(as.numeric(names(probabilities)[winner_index]))
}

r_from_m_s_eq_7 <- function(ntax, ef, time) {
	return(
		(1/time)*
		(
			log(
				0.5*ntax*(1-ef^2) +
				2*ef + 
				0.5*(1-ef)*sqrt(ntax*(ntax*(ef^2)-8*ef+2*ntax*ef+ntax))
			) - log(2)
		)
	)
}

impute_missing <- function(yule_plot_data_all, lambda=0.1, nsims_original=25000, min_time=0, max_time=50) {
	yule_plot_data_all$imputed <- FALSE
	while(nrow(yule_plot_data_all) < nsims_original) {
		focal_time <- runif(min=min_time, max=max_time, n=1)
		prob_2 <- get_m_s_prob_nt_eq_i(final_species_count=2, netdiv=lambda, ef=0, time=focal_time)
		if(rbinom(1, 1, prob=as.numeric(prob_2))==1) {
			yule_plot_data_all <- rbind(yule_plot_data_all, data.frame(params=0, lower.CI=NA, upper.CI=NA, ages=focal_time, ntips=2, imputed=TRUE))
		}
	}
	return(yule_plot_data_all)
}

get_trees_hd <- function() {
	trees <- read.csv("data/henaodiaz_etal_trees-master/henao_diaz_etal_tree-samplesize.csv")
	trees$Clade.age <- as.numeric(trees$Clade.age)
	trees$Clade.size.randomized <- sample(trees$Clade.size)
	return(trees)
}

get_trees_tall_df <- function(trees) {
	trees_tall <- pivot_longer(trees, cols=c("Clade.size", "Clade.size.randomized"), names_to="Clade.size_type", values_to="Clade.size")
	trees_tall$Clade.age <- as.numeric(trees_tall$Clade.age)
	trees_tall$mom_rate <- (log(trees_tall$Clade.size) - log(2))/trees_tall$Clade.age
	trees_tall$log_mom_rate <- log(trees_tall$mom_rate)
	trees_tall$log_clade_age <- log(trees_tall$Clade.age)
	return(trees_tall)
}

get_expected_rates_bd_raw <- function() {
	netdiv_vector <- c(0.1)
	ef_vector <- c(0,0.5, .95)
	time_vector <- seq(from=1, to=10, by=.2)
	ntax_vector <- seq(from=0, to=1000, by=1)
	results <- data.frame()
	for (netdiv in netdiv_vector) {
		for(ef in ef_vector) {
			for (time in time_vector) {
				for (ntax in ntax_vector) {
					probN <- get_m_s_prob_nt_eq_i(final_species_count=ntax, netdiv=netdiv, ef=ef, time=time, starting_species_count=2) 
					local_results <- data.frame(netdiv=netdiv, ef=ef, time=time, ntax=ntax, log_probN=as.numeric(log(probN)), pure_birth_r=(log(ntax) - log(2))/time, bd_r=r_from_m_s_eq_7(ntax, ef, time))
					local_results$pure_birth_r_weighted <- 0
					try({local_results$pure_birth_r_weighted <- probN * local_results$pure_birth_r})
					local_results$bd_r_weighted <- 0
					try({local_results$bd_r_weighted <- probN * local_results$bd_r})
					results <- rbind(results, local_results)
				}
				save(results, file="results.Rsave")

			}
		}
	}
	return(results)
}

get_expected_rates_bd_summed <- function(results, minN=0) {
	
	netdiv_vector <- sort(unique(results$netdiv))
	ef_vector <- sort(unique(results$ef))
	time_vector <- sort(unique(results$time))
	results <- subset(results, ntax>=minN)
	
	expected_results <- data.frame()
	for (f_netdiv in netdiv_vector) {
		for(f_ef in ef_vector) {
			for (f_time in time_vector) {
				focal_results <- subset(results, netdiv==f_netdiv & ef==f_ef & time==f_time)
				focal_results <- focal_results[is.finite(focal_results$bd_r_weighted),]
				foca_results <- focal_results[is.finite(focal_results$pure_birth_r_weighted),]
				#focal_results <- focal_results[is.finite(Rmpfr::asNumeric(focal_results$bd_r_weighted)),]
				#focal_results <- focal_results[is.finite(Rmpfr::asNumeric(focal_results$pure_birth_r_weighted)),]
				#expected_results <- rbind(expected_results, data.frame(netdiv=f_netdiv, ef=f_ef, time=f_time, expected_pure_birth_r=Rmpfr::asNumeric(sum(focal_results$pure_birth_r_weighted)), expected_bd_r=Rmpfr::asNumeric(sum(focal_results$bd_r_weighted))))
				expected_results <- rbind(expected_results, data.frame(netdiv=f_netdiv, ef=f_ef, time=f_time, expected_pure_birth_r=sum(focal_results$pure_birth_r_weighted), expected_bd_r=sum(focal_results$bd_r_weighted)))

			}
		}
	}

	expected_results$appropriate_r <- expected_results$expected_bd_r
	expected_results$appropriate_r[expected_results$ef==0] <- expected_results$expected_pure_birth_r[expected_results$ef==0]	
	return(expected_results)
}

get_expected_results_no_correlation <- function() {
	time_vector <- seq(from=1, to=10, by=.2)
	ntax_vector <- seq(from=0, to=1000, by=1)
	results <- expand.grid(time=time_vector, ntax=ntax_vector)
	results$expected_pure_birth_r <- (log(results$ntax) - log(2))/results$time
	results$expected_bd_ef_0.95 <- r_from_m_s_eq_7(results$ntax, ef=0.95, results$time)
	return(results)
}

get_hd_fossils <- function() {
	fossils <- read.csv("https://raw.githubusercontent.com/mwpennell/macro_sadler/master/output/summary_paleo_results.csv")	
	return(fossils)
}

get_poisson_rate <- function() {
	r_vector <- c(0.1, 0.5)
	epsilon_0_vector <- c(0, 0.1, 0.5)
	t_vector <- seq(0.1, 10, 0.1)
	k_vector <- seq(0, 100, 1)
	df <- expand.grid(r=r_vector, epsilon_0=epsilon_0_vector, t=t_vector, k=k_vector)
	df$prob <- dpois(df$k, lambda=df$r*df$t)
	df$rate <- df$k/df$t
	df$rate_times_prob <- df$rate*df$prob
	return(df)
}

get_poisson_rate_summed <- function(df) {
	df <- df |> dplyr::group_by(r, epsilon_0, t) |> dplyr::summarize(expected_rate=sum(rate_times_prob))
	return(df)
}

do_time_error <- function(error_sd=0.1) {
	true_times <- runif(10000, 0, 10)
	observed_times <- rnorm(10000, true_times, sd=error_sd*true_times)
	true_times <- true_times[observed_times>0]
	observed_times <- observed_times[observed_times>0]
	true_x <- true_times*0.1
	recon_rate <- true_x/observed_times
	return(data.frame(true_times=true_times, observed_times=observed_times, true_x=true_x, recon_rate=recon_rate))
}


convert_poisson_to_negbinom_params <- function(lambda, variance_increase) {
	# variance = lambda*(1+lambda/r)
	# desired variance = variance_increase + lambda = lambda*(1+lambda/r)
	# variance_increase/lambda + 1 = 1 + lambda/r
	# variance_increase/lambda = lambda/r
	# r = lambda^2 / variance_increase
	r = lambda^2 / variance_increase
	p = r/(r+lambda)
	return(c(r=r, p=p))
}

compute_poisson_negbinom_example <- function() {
	true_rate <- 0.1
	time_vector <- seq(1, 10, 0.1)	
	variance_increase_vector <- c(0.00001, 0.5, 3)
	k_vector <- seq(0, 1000, 1)
	df <- expand.grid(true_rate=true_rate, time=time_vector, k=k_vector, variance_increase=variance_increase_vector)
	df$lambda <- df$true_rate*df$time
	df$variance <- df$lambda+df$variance_increase
	negbinom_params <- apply(df, 1, function(x) convert_poisson_to_negbinom_params(x["lambda"], x["variance_increase"]))
	df$r <- negbinom_params[1,]
	df$p <- negbinom_params[2,]
	df$prob_true <- dnbinom(df$k, size=df$r, prob=df$p)
	df$prob_pois <- dpois(df$k, lambda=df$lambda)
	df$rate <- df$k/df$time
	df$rate_times_prob_true <- df$rate*df$prob_true
	return(df)

}

get_neg_binom_example_summed <- function(df) {
	df <- df |> dplyr::group_by(true_rate, time, variance_increase) |> dplyr::summarize(recon_rate=sum(rate_times_prob_true))
	return(df)
}

compute_poisson_plus_poisson_example <- function() {
	r_vector <- c(0.1, 0.5)
	lambda0_vector <- c(0, 1, 5)
	t_vector <- seq(0.1, 50, 0.2)
	k_t_vector <- seq(0, 200, 1)
	k_0_vector <- seq(0, 200, 1)
	df <- expand.grid(r=r_vector, lambda0=lambda0_vector, t=t_vector, k_t=k_t_vector, k_0=k_0_vector)
	df$prob_k_t <- dpois(df$k_t, lambda=df$r*df$t)
	df$prob_k_0 <- dpois(df$k_0, lambda=df$lambda0)
	df$estimated_rate_raw <- (df$k_t+df$k_0)/df$t
	df$estimated_rate_times_prob_k_t_k_0 <- df$estimated_rate_raw*df$prob_k_t*df$prob_k_0
	return(df)
}

compute_constant_plus_time_error_example <- function() {
	lambda0_vector <- c(0, 1, 5)
	t_var_vector <- c(0, 0.1, 0.5)
	step_size <- 0.01
	t_quantile_vector <- seq(from=step_size, to=1-step_size, by=step_size)
	t_vector <- seq(0.1, 50, 0.2)
	k_0_vector <- seq(0, 200, 1)
	df <- expand.grid( lambda0=lambda0_vector, t=t_vector, k_0=k_0_vector, t_var=t_var_vector, t_quantile=t_quantile_vector)
	df$prob_k_0 <- dpois(df$k_0, lambda=df$lambda0)
	df$t_location <- log(df$t^2 / sqrt(df$t_var + df$t^2))
	df$t_shape <- sqrt(log(1 + (df$t_var / df$t^2)))
	df$t_from_quantile <- qlnorm(df$t_quantile, meanlog=df$t_location, sdlog=df$t_shape)
	df$estimated_rate_raw_true <- (df$k_0)/df$t
	df$estimated_rate_quantile <- (df$k_0)/df$t_from_quantile
	return(df)
}

compute_poisson_time_plus_time_error_example <- function() {
	max_time <- 20
	r_vector <- c(0, 0.1, 0.5)
	lambda0_vector <- c(0, 1, 5)
	t_var_vector <- c(0, 0.1, 0.5)
	step_size <- 0.01
	t_quantile_vector <- seq(from=step_size/2, to=1-step_size/2, by=step_size)
	t_vector <- exp(seq(from=log(0.1), to=log(max_time), length.out=50)) # do even spacing.
	#k_0_vector <- seq(0, 200, 1)
	#k_t_vector <- seq(0, 200, 1)
	k_0_vector <- seq(0, 100, 1)
	k_t_vector <- seq(0, 100, 1)
	df <- expand.grid(r=r_vector, lambda0=lambda0_vector, t=t_vector, k_0=k_0_vector, k_t=k_t_vector, t_var=t_var_vector)
	df$prob_k_t <- dpois(df$k_t, lambda=df$r*df$t)
	df$prob_k_0 <- dpois(df$k_0, lambda=df$lambda0)
	df$t_location <- log(df$t^2 / sqrt(df$t_var + df$t^2))
	df$t_shape <- sqrt(log(1 + (df$t_var / df$t^2)))
	df$estimated_rate_raw <- (df$k_t+df$k_0)/df$t
	df$estimated_rate_times_prob_k_t_k_0 <- df$estimated_rate_raw*df$prob_k_t*df$prob_k_0
	df$estimated_rate_quantile_prob_k_t_k_0_prob_quantile <- NA
	df$quantile_rate_numerator <- (df$k_t+df$k_0)*(df$prob_k_t*df$prob_k_0)
	pboptions(type="txt")
	df$estimated_rate_quantile_prob_k_t_k_0_prob_quantile <- mcmapply(FUN=function(t_quantile, t_location, t_shape, quantile_rate_numerator, number_of_quantiles) {
		t_from_quantile <- log(qlnorm(t_quantile, meanlog=t_location, sdlog=t_shape))
		local_rates <- (quantile_rate_numerator / t_from_quantile)/number_of_quantiles
		local_rates
	}, t_quantile=t_quantile_vector, t_location=df$t_location, t_shape=df$t_shape, quantile_rate_numerator=df$quantile_rate_numerator, number_of_quantiles=length(t_quantile_vector), mc.cores=6, mc.preschedule=TRUE, mc.cleanup=TRUE)
	return(df)
}

compute_poisson_with_simple_time_errors <- function() {
	max_time <- 20
	r_vector <- c(0, 0.1, 0.5)
	lambda0_vector <- c(0, 1, 5)
	t_var_vector <- c(0, 0.1, 0.5)
	t_vector <- exp(seq(from=log(0.1), to=log(max_time), length.out=50)) # do even spacing.
	k_0_vector <- seq(0, 100, 1)
	k_t_vector <- seq(0, 100, 1)
	df <- expand.grid(r=r_vector, lambda0=lambda0_vector, t=t_vector, k_0=k_0_vector, k_t=k_t_vector, t_var=t_var_vector, bias_direction=c("over", "under"))
	df$prob_k_t <- dpois(df$k_t, lambda=df$r*df$t)
	df$prob_k_0 <- dpois(df$k_0, lambda=df$lambda0)
	df$numerator_with_prob <- (df$k_t+df$k_0)*(df$prob_k_t*df$prob_k_0)
	df$t_estimated <- df$t
	over_ones <- df$bias_direction=="over"
	under_ones <- df$bias_direction=="under"
	df$t_estimated[over_ones] <- df$t_estimated[over_ones] * (1+df$t_var[over_ones])
	df$t_estimated[under_ones] <- df$t_estimated[under_ones] * (1-df$t_var[under_ones])
	return(df)
}

get_poisson_with_simple_time_errors_summed <- function(df) {
	df <- df |> dplyr::group_by(r, lambda0, t, t_var, bias_direction) |> dplyr::summarize(estimated_rate=sum(numerator_with_prob)/t_estimated)
	return(df)
}


# TODO Poisson, but rather than error in T, do it with T true, T 50% too high, T 50% too low

# Have variance in T be like 10% of true, or have this be bigger as we go deeper in tree, or have bias such that older ones are too old but never too young



get_constant_plus_time_error_summed <- function(df) {
	df <- df |> dplyr::group_by(lambda0, t, t_var) |> dplyr::summarize(estimated_rate_noise=mean(estimated_rate_quantile), estimated_rate_true=mean(estimated_rate_raw_true))
	return(df)
}

get_poisson_time_plus_time_error_summed <- function(df) {
	df <- df |> dplyr::group_by(lambda0, t, t_var, r) |> dplyr::summarize(estimated_rate=sum(estimated_rate_quantile_prob_k_t_k_0_prob_quantile))
	return(df)
}

falling_rock <- function() {
	t_vector <- seq(0, 60, 0.1)
	g <- 9.8 #m/s^2, just doing speed, not velocity
	x <- 0.5*g*t_vector^2	
	speed <- x/t_vector
	x_noisy <- rnorm(length(x), x, 0.2)
	t_noisy <- rnorm(length(t_vector), t_vector, 0.2)
	speed_noisy <- x_noisy/t_noisy
	
}

get_poisson_plus_poisson_rate_summed <- function(df) {
	df <- df |> dplyr::group_by(r, lambda0, t) |> dplyr::summarize(estimated_rate=sum(estimated_rate_times_prob_k_t_k_0))
	df$var_from_time <- df$t*df$r
	df$var_from_lambda0 <- df$lambda0
	df$var_total <- df$var_from_time + df$var_from_lambda0
	df$var_fraction_from_mserr <- df$var_from_lambda0/df$var_total
	return(df)
}

get_substitution_data <- function() {
	subs <- read.csv("data/digitized_rolland_etal_2023_fig1e.csv")
	colnames(subs) <- c("log_time_MY", "log_rate_sites_per_MY")
	subs$time_MY <- exp(subs$log_time_MY)
	subs$rate_sites_per_MY <- exp(subs$log_rate_sites_per_MY)
	subs$numerator <- subs$rate_sites_per_MY*subs$time_MY	
	subs$rate <- abs(subs$rate_sites_per_MY)
	subs$my <- subs$time_MY
	subs$weights <- 1
	return(subs)
}

get_gingerich_data <- function() {
	gingerich <- read.csv("data/digitized_Gingerich_1983_Fig1b.csv")
	colnames(gingerich) <- c("log_time_MY", "log_rate_darwins")
	gingerich$time_MY <- exp(gingerich$log_time_MY)
	gingerich$rate_darwins <- exp(gingerich$log_rate_darwins)
	gingerich$numerator <- gingerich$rate_darwins*gingerich$time_MY
	gingerich$my <- gingerich$time_MY
	gingerich$rate <- abs(gingerich$rate_darwins)
	gingerich$weights <- 1
	return(gingerich)
}

get_uyeda_data <- function(apply_bounds=TRUE) {
	raw <- read.csv("data/doi_10.5061_dryad.7d580__v1/Dryad7.csv")
	cleaned <- subset(raw, BodySizeCorrelated=="1")
	if(apply_bounds) { # we're going to use the bounds of Rolland et al. (2023)
		cleaned <- subset(cleaned, Darwins < 19000)
		cleaned <- subset(cleaned, Years < 1e8)
		cleaned <- subset(cleaned, log(abs(Darwins)) > -7)
	}
	cleaned$numerator <- abs(cleaned$Darwins)*(cleaned$Years*1e-6)
	cleaned$my <- (1e-6) * cleaned$Years
	cleaned$rate <- abs(cleaned$Darwins)
	cleaned$weights <- 1
	# raw2 <- read.csv("../blunderbus/data/doi_10.5061_dryad.7d580__v1/Phylogeniesbynode.csv")
	# cleaned2 <- data.frame(numerator=abs(raw2$d), my=raw2$Interval.my, rate=abs(raw2$d)/raw2$Log10Years, weights=1)
	# cleaned <- dplyr::bind_rows(cleaned, cleaned2)
	# cleaned <- subset(cleaned, my>0 & rate>0)
	return(cleaned)
}

get_diversification_data <- function() {
	#div_data <- read.csv("../henaodiaz_etal_trees-master/henao_diaz_etal_tree-samplesize.csv")
	div_data <-  read.csv("data/01-bd_df_plot.csv")
	#div_data$numerator <- log(div_data$Clade.size)-log(2)
	#div_data$numerator <- div_data$net.div * div_data$tot_bl # rate divided by total branch length
	div_data$numerator <- div_data$lambda * div_data$tot_bl # rate divided by total branch length
	div_data$my <- div_data$Clade.age # but this is what we actually plot it against
	div_data$total_time <- div_data$tot_bl
	div_data$rate <- div_data$net.div
	div_data$weights <- 1
	return(div_data)
}

aggregate_all_data <- function() {
	all_data <- data.frame()
	
	substitution_data <- get_substitution_data()
	all_data <- dplyr::bind_rows(all_data, data.frame(citation="Ho et al. (2005)", type="DNA substitution", my=substitution_data$my, numerator=substitution_data$numerator, rate=substitution_data$rate, xvals=substitution_data$my, rateunits="Substitutions/MY", xrateunits="Time (MY)", weights=substitution_data$weights))
	
	gingerich_data <- get_gingerich_data()
	all_data <- dplyr::bind_rows(all_data, data.frame(citation="Gingerich (1983)", type="Morphology", my=gingerich_data$my, numerator=gingerich_data$numerator, rate=gingerich_data$rate, xvals=gingerich_data$my, rateunits="Darwins", xrateunits="Time (MY)",weights=gingerich_data$weights))
	
	uyeda_data <- get_uyeda_data()
	all_data <- dplyr::bind_rows(all_data, data.frame(citation="Uyeda et al. (2011)",type="Body size", my=uyeda_data$my, numerator=uyeda_data$numerator, rate=uyeda_data$rate, xvals=uyeda_data$my, rateunits="Darwins", xrateunits="Time (MY)", weights=uyeda_data$weights))
	
	diversification_data <- get_diversification_data()
	all_data <- dplyr::bind_rows(all_data, data.frame(citation="Henao Diaz et al. (2019)", type="Diversification", my=diversification_data$my, numerator=diversification_data$numerator, total_time=diversification_data$total_time, rate=diversification_data$rate, xvals=diversification_data$my, rateunits="Species/MY", xrateunits="Time (MY)",weights=diversification_data$weights))
	
	#all_data <- dplyr::bind_rows(all_data, data.frame(citation="Henao Diaz et al. (2019) vs total time", type="Diversification", my=diversification_data$total_time, numerator=diversification_data$numerator, total_time=diversification_data$total_time, rate=diversification_data$rate, xvals=diversification_data$total_time, rateunits="Species/MY", xrateunits="Time (MY)",weights=diversification_data$weights))

	
	
	
	modern_extinction <- get_modern_extinction_no_gray(get_modern_extinction())
	all_data <- dplyr::bind_rows(all_data, data.frame(citation="Barnosky et al. (2011)", type="Modern extinction", my=modern_extinction$my, numerator=modern_extinction$numerator, rate=modern_extinction$rate, xvals=modern_extinction$my, rateunits="Extinctions/MY", xrateunits="Time (MY)", weights=modern_extinction$weights))
	
	data("yule_sim", package="hyperr8")
	yule_sim$type="Simulation"
	yule_sim$my=yule_sim$time
	yule_sim$numerator=yule_sim$rate*yule_sim$total_brlen
	yule_sim$xvals=yule_sim$my
	yule_sim$rateunits="Species/MY"
	yule_sim$xrateunits="Time (MY)"
	yule_sim$weights=1
	yule_sim$citation <- "Pure birth simulation"
	yule_sim$total_time <- yule_sim$total_brlen
	all_data <- dplyr::bind_rows(all_data, as.data.frame(yule_sim))
	
	
	# yule_sim_all_vs_total_time <- yule_sim
	# yule_sim_all_vs_total_time$my <- yule_sim_all_vs_total_time$total_time
	# yule_sim_all_vs_total_time$xvals <- yule_sim_all_vs_total_time$my
	# yule_sim_all_vs_total_time$citation <- "Pure birth simulation (total time)"
	# yule_sim_all_vs_total_time$time <- yule_sim_all_vs_total_time$total_time
	# all_data <- dplyr::bind_rows(all_data, as.data.frame(yule_sim_all_vs_total_time))
	
	# yule_sim_no_zeros <- yule_sim
	# yule_sim_no_zeros <- subset(yule_sim_no_zeros, rate>0)
	# yule_sim_no_zeros$citation <- "Pure birth simulation (no zeros)"
	# all_data <- dplyr::bind_rows(all_data, as.data.frame(yule_sim_no_zeros))
	
	all_data$denominator <- all_data$my
	all_data$denominator[!is.na(all_data$total_time)] <- all_data$total_time[!is.na(all_data$total_time)]
	all_data$time <- all_data$my
	all_data$datum_id <- sequence(nrow(all_data))
	return(all_data)
}

get_funny_yule <- function() {
	all_data <- data.frame()
	data("yule_sim", package="hyperr8")
	yule_sim$type="Simulation"
	yule_sim$my=yule_sim$time
	yule_sim$numerator=yule_sim$rate*yule_sim$total_brlen
	yule_sim$xvals=yule_sim$my
	yule_sim$rateunits="Species/MY"
	yule_sim$xrateunits="Time (MY)"
	yule_sim$weights=1
	yule_sim$citation <- "Pure birth simulation"
	yule_sim$total_time <- yule_sim$total_brlen
	
	yule_sim_all_vs_total_time <- yule_sim
	yule_sim_all_vs_total_time$my <- yule_sim_all_vs_total_time$total_time
	yule_sim_all_vs_total_time$xvals <- yule_sim_all_vs_total_time$my
	yule_sim_all_vs_total_time$citation <- "Pure birth simulation (total time)"
	yule_sim_all_vs_total_time$time <- yule_sim_all_vs_total_time$total_time
	all_data <- dplyr::bind_rows(all_data, as.data.frame(yule_sim_all_vs_total_time))
	
	yule_sim_no_zeros <- yule_sim
	yule_sim_no_zeros <- subset(yule_sim_no_zeros, rate>0)
	yule_sim_no_zeros$citation <- "Pure birth simulation (no zeros)"
	all_data <- dplyr::bind_rows(all_data, as.data.frame(yule_sim_no_zeros))
	return(all_data)
}

merge_yule_funny_and_regular <- function(hyperr8_analysis, hyperr8_analysis_yule_funny) {
	yule_regular <- subset(hyperr8_analysis, dataset=="Pure birth simulation")
	return(rbind(yule_regular, hyperr8_analysis_yule_funny))	
}

randomize_within_dataset <- function(all_data) {
	result <- data.frame()
	datasets <- unique(all_data$citation)
	for (dataset in datasets) {
		print(dataset)
		focal_data <- subset(all_data, citation==dataset)
		focal_data$rate <- sample(focal_data$numerator)/focal_data$denominator
		result <- rbind(result, focal_data)
	}
	return(result)
}


plot_individual_dataset <- function(df, dolog=FALSE, nreps=10) {
	alpha_value <- 0.1
	if(nrow(df)>1000) {
		alpha_value <- 0.01
	}
	if(nrow(df)<200) {
		alpha_value <- 0.5
	}
g <- ggplot(df, aes(x=xvals, y=rate)) + geom_point(alpha=alpha_value, shape=20) + xlab(unique(df$xrateunits)) + ylab(unique(df$rateunits)) + theme_bw() + scale_y_continuous(labels = scales::comma) #+ scale_x_continuous(labels = KJHSXDFGHJNJBVGCFXDZSAAZsxdczsA`SZ2) 

	if(dolog) {
		linear_model <- lm(log10(rate) ~ 1 + log10(my), data=df)	
		g <- g + scale_y_log10(labels = scales::comma) + scale_x_log10(labels = scales::comma) + geom_smooth(method="lm", span=1, col="red") #+ geom_abline(intercept=linear_model$coefficients[1], slope=linear_model$coefficients[2], colour="red")
	} else {
		g <- g + ylim(0, 1.3*max(df$rate))+ geom_smooth(colour='red', se=TRUE)
	}
	for (i in sequence(nreps)) {
		df_sampled  <- df
		denominator <- df_sampled$my
		if(!is.na(df_sampled$total_time[1])) {
			denominator <- df_sampled$total_time
		}
		#df_sampled$rate <- sample(df_sampled$numerator)/df_sampled$my
		df_sampled$rate <- sample(df_sampled$numerator)/denominator

		if(!dolog) {
			g <- g + geom_smooth(alpha=alpha_value, data=df_sampled, se=FALSE, colour="blue", linewidth=0.2)
		}
		else {
			linear_model <- lm(log10(rate) ~ 1 + log10(my), data=df_sampled)	
			g <- g + geom_abline(intercept=linear_model$coefficients[1], slope=linear_model$coefficients[2], colour="blue")
		}
	}
	if(dolog) {
		df$y <- log10(df$rate)
		df$x <- log10(df$xvals)
		linear_model <- lm(y ~ 1 + offset(-1*x), data=df)
		g <- g + geom_abline(intercept=linear_model$coefficients[1], slope=-1, colour="green", linetype="dashed") + coord_fixed(ylim=range(c(df$rate, df$rate_deviation)), clip='on')
	}
	return(g)
}

compute_lm_table <- function(all_data) {
	results <-data.frame()
	for (focal_citation in unique(all_data$citation)) {
		focal_data <- subset(all_data, citation==focal_citation)
		focal_data$log_rate <- log(focal_data$rate)
		focal_data$log_my <- log(focal_data$my)
		linear_model_fixed_slope <- lm(log_rate ~ 1 + offset(-1*log_my), data=focal_data)
		linear_model_floating_slope <- lm(log_rate ~ 1 + log_my, data=focal_data)
		AICs <- c(AIC(linear_model_fixed_slope), AIC(linear_model_floating_slope))
		names(AICs) <- c("fixed_slope", "floating_slope")
		deltaAICs <- AICs - min(AICs)
		AICweights <- exp(-0.5*deltaAICs)/sum(exp(-0.5*deltaAICs))
		results <- rbind(results, data.frame(
			citation=focal_citation, 
			floating_slope_slope=paste0(round(linear_model_floating_slope$coefficients[2], 2), " (", round(confint(linear_model_floating_slope)[2,1], 2), ", ", round(confint(linear_model_floating_slope)[2,2], 2), ")"),
			floating_slope_intercept=paste0(round(linear_model_floating_slope$coefficients[1], 2), " (", round(confint(linear_model_floating_slope)[1,1], 2), ", ", round(confint(linear_model_floating_slope)[1,2], 2), ")"),
			fixed_slope_intercept=paste0(round(linear_model_fixed_slope$coefficients[1], 2), " (", round(confint(linear_model_fixed_slope)[1,1], 2), ", ", round(confint(linear_model_fixed_slope)[1,2], 2), ")"),
			floating_slope_deltaAIC=deltaAICs["floating_slope"],
			fixed_slope_deltaAIC=deltaAICs["fixed_slope"],
			floating_slope_AICweight=AICweights["floating_slope"],
			fixed_slope_AICweight=AICweights["fixed_slope"]
		))
	}
	rownames(results) <- NULL
	return(results)
}

format_lm_table <- function(lm_table) {
	lm_table_pretty <- lm_table[,c('citation', 	'floating_slope_slope', 'floating_slope_intercept')]
	colnames(lm_table_pretty) <- c('Citation', 'Slope', 'Intercept')
	return(lm_table_pretty)
}

subtract_linear_fixed_slope_model <- function(all_data) {
	results <-data.frame()
	for (focal_citation in unique(all_data$citation)) {
		focal_data <- subset(all_data, citation==focal_citation)
		focal_data$log_rate <- log(focal_data$rate)
		focal_data$log_my <- log(focal_data$my)
		linear_model_fixed_slope <- lm(log_rate ~ 1 + offset(-1*log_my), data=focal_data)
		focal_data$log_rate_prediction <- predict(linear_model_fixed_slope, newdata=focal_data)
		focal_data$log_rate_deviation <- focal_data$log_rate - focal_data$log_rate_prediction
		focal_data$rate_deviation <- exp(focal_data$log_rate_deviation)
		results <- rbind(results, focal_data)
	}
	rownames(results) <- NULL
	return(results)
}

plot_subtracted_dataset <- function(df, dolog=FALSE) {
	alpha_value <- 0.1
	if(nrow(df)>1000) {
		alpha_value <- 0.01
	}
	if(nrow(df)<200) {
		alpha_value <- 0.5
	}
	g <- ggplot(df, aes(x=xvals, y=rate_deviation)) + geom_point(alpha=alpha_value, shape=20) + xlab(unique(df$xrateunits)) + ylab(unique(df$rateunits)) + theme_bw() + scale_y_continuous(labels = scales::comma) + scale_x_continuous(labels = scales::comma) 
	
	if(dolog) {
		linear_model <- lm(log10(rate_deviation) ~ 1 + log10(my), data=df)	
		g <- g + scale_y_log10(labels = scales::comma) + scale_x_log10(labels = scales::comma) + geom_smooth(method="lm", span=1, col="red") + coord_fixed(ylim=range(c(df$rate, df$rate_deviation)), clip='on') #+ geom_abline(intercept=linear_model$coefficients[1], slope=linear_model$coefficients[2], colour="red")
	} else {
		g <- g + ylim(0, 1.3*max(df$rate_deviation))+ geom_smooth(colour='red', se=TRUE)
	}
	return(g)
}

get_bm_model <- function() {
	sigma_0_vector <- c(0, 0.1, 0.5)
	sigma_t_vector <- c(0.1, 0.5)
	time_vector <- seq(0.1, 10, 0.1)
	df <- expand.grid(sigma_0=sigma_0_vector, sigma_t=sigma_t_vector, time=time_vector)	
	df$rate <- df$sigma_0/df$time + df$sigma_t/sqrt(df$time)
	return(df)
}

# From https://doi.org/10.1038/nature09678
get_modern_extinction <- function() {
	gray_dots <- read.csv("data/modernextinction_gray_dots.csv")	
	gray_and_orange_dots <- read.csv("data/modernextinction_gray_and_orange_dots.csv")
	gray_and_orange_and_red_dots <- read.csv("data/modernextinction_gray_and_orange_and_red_dots.csv")
	gray_and_orange_and_red_dots$type <- "General"
	gray_and_orange_and_red_dots$type[(nrow(gray_dots)+1) : nrow(gray_and_orange_dots)] <- "Post 2010"
	gray_and_orange_and_red_dots$type[(nrow(gray_and_orange_dots)+1) : nrow(gray_and_orange_and_red_dots)] <- "Pleistocene extinction"
	colnames(gray_and_orange_and_red_dots) <- c("years", "extinctions_per_my", "type")
	gray_and_orange_and_red_dots$my <- (10^-6) * gray_and_orange_and_red_dots$years
	gray_and_orange_and_red_dots$weights <- 1
	gray_and_orange_and_red_dots$rate <- gray_and_orange_and_red_dots$extinctions_per_my
	gray_and_orange_and_red_dots$numerator <- gray_and_orange_and_red_dots$my*gray_and_orange_and_red_dots$extinctions_per_my
	return(gray_and_orange_and_red_dots)
}

get_modern_extinction_no_gray <- function(x) {
	return(subset(x, type!="General"))	
}

simulate_counts <- function(chosen_model="flat", rate_multiplier=1) {
	max_time=10
	times <- sort(runif(1000000, 0, max_time))
	events <- rep(NA, length(times))
	if(chosen_model=="flat") {
		events <- rbinom(n=length(times), size=1, prob=0.05*rate_multiplier)
	}	
	if(chosen_model=="two_rates") {
		events <- c(rbinom(n=length(times)/2, size=1, prob=0.01*rate_multiplier), rbinom(n=length(times)/2, size=1, prob=0.09*rate_multiplier))
	}
	if(chosen_model=="smooth_increase") {
		events <- rbinom(n=length(times), size=1, prob=0.01*rate_multiplier + 0.09*rate_multiplier*times/max_time)
	}
	return(data.frame(times=times, events=events))	
}

prune_to_observed <- function(simulated_counts) {
	return(subset(simulated_counts, events==1))	
}

do_sliding_window_analysis <- function(simulated_counts, time_width=0.5, time_step_width=0.02, event_width=20) {
	results <- data.frame()

	# sliding window analysis
	time_start=0
	time_end=time_width+time_start
	while(time_end<max(simulated_counts$times)) {
		focal <- subset(simulated_counts, times>=time_start & times<time_end)
		focal_event_count <- sum(focal$events)
		results <- rbind(results, data.frame(time_start=time_start, time_end=time_end, time_mid=(time_start+time_end)/2, event_count=focal_event_count, rate=focal_event_count/time_width, method="sliding_time_window"))
		time_start <- time_start + time_step_width
		time_end <- time_end + time_step_width	
	}
	
	# time_to_present analysis
	time_end=max(simulated_counts$times)
	while(time_end>min(simulated_counts$times)) {
		focal <- subset(simulated_counts, times<=time_end)
		if(nrow(focal)>1) {
			time_start <- min(focal$times)
			total_time <- max(focal$times) - min(focal$times)

			focal_event_count <- sum(focal$events)
			results <- rbind(results, data.frame(time_start=time_start, time_end=time_end, time_mid=(time_start+time_end)/2, event_count=focal_event_count, rate=focal_event_count/total_time, method="time_to_present"))
		}
		time_end <- time_end - time_step_width
	}
	
	# step window analysis

	time_start=0
	time_end=time_width+time_start
	while(time_end<max(simulated_counts$times)) {
		focal <- subset(simulated_counts, times>=time_start & times<time_end)
		focal_event_count <- sum(focal$events)
		results <- rbind(results, data.frame(time_start=time_start, time_end=time_end, time_mid=(time_start+time_end)/2, event_count=focal_event_count, rate=focal_event_count/time_width, method="stepping_time_window"))
		time_start <- time_start + time_width
		time_end <- time_end + time_width	
	}
	
	#sliding event window analysis
	
	# row_end <- event_width
	# while(row_end <= nrow(simulated_counts)) {
	# 	focal <- simulated_counts[(1+row_end-event_width):row_end,]
	# 	focal_event_count <- sum(focal$events)
	# 	total_time <- max(focal$times) - min(focal$times)
	# 	results <- rbind(results, data.frame(time_start=min(focal$times), time_end=max(focal$times), time_mid=(min(focal$times)+max(focal$times))/2, event_count=focal_event_count, rate=focal_event_count/total_time, method="sliding_event_window"))
	# 	row_end <- row_end + 1
	# }
	
	# stepping event window analysis
	
	# row_end <- event_width
	# while(row_end <= nrow(simulated_counts)) {
	# 	focal <- simulated_counts[(1+row_end-event_width):row_end,]
	# 	focal_event_count <- sum(focal$events)
	# 	total_time <- max(focal$times) - min(focal$times)
	# 	results <- rbind(results, data.frame(time_start=min(focal$times), time_end=max(focal$times), time_mid=(min(focal$times)+max(focal$times))/2, event_count=focal_event_count, rate=focal_event_count/total_time, method="stepping_event_window"))
	# 	row_end <- row_end + event_width
	# }
	return(results)
}

do_modern_extinction_analysis <- function(year_span=519, nspecies=30873, E_per_MSY=1) {
	return(rpois(n=year_span, lambda=nspecies*E_per_MSY/1e6))
}

get_simulated_modern_extinction_rates <- function(year_min=10, year_max=519, nreps=100, nspecies=30873) {
	results <- data.frame()
	for (year_span in seq(from=year_min, to=year_max, by=1)) {
		for (rep in seq(from=1, to=nreps, by=1)) {
			focal <- do_modern_extinction_analysis(year_span=year_span, nspecies=nspecies)
			results <- rbind(results, data.frame(year_span=year_span, rep=rep, rate=sum(focal)*1e6/(year_span*nspecies)))
		}
	}
	return(results)
}

function_f1 <- function(par, focal_data, dolog=TRUE) {
	varepsilon_0_plus_k <- par[1]
	#log(\hat{r}(t)) = log(\varepsilon_0 + k) - log(\hat{t})
	if(dolog) {
		return(log(varepsilon_0_plus_k) - focal_data$log_time)
	} else {
		return(varepsilon_0_plus_k/exp(focal_data$log_time))
	}
}


	

function_f2_general <- function(par, focal_data, dolog=TRUE) {
	varepsilon_0 <- par[1]
	k <- par[2]
	if(dolog) {
		return(log(varepsilon_0 + k/exp(focal_data$log_time)) - focal_data$log_time)
	} else {
		return((varepsilon_0 + k/exp(focal_data$log_time))/exp(focal_data$log_time))
	}
}

function_f3_general <- function(par, focal_data, dolog=TRUE) {
	varepsilon_0 <- par[1]
	k <- par[2]
	if(dolog) {
		return(log(varepsilon_0 + k*exp(focal_data$log_time)) - focal_data$log_time)
	} else {
		return((varepsilon_0 + k*exp(focal_data$log_time))/exp(focal_data$log_time))
	}
}

function_f4_general <- function(par, focal_data, dolog=TRUE) {
	varepsilon_0 <- par[1]
	k <- par[2]
	a <- par[3]
	if(dolog) {
		return(log(varepsilon_0 + k*(focal_data$log_time)^a) - focal_data$log_time)
	} else {
		return((varepsilon_0 + k*(focal_data$log_time)^a)/exp(focal_data$log_time))
	}
}


optimize_rate_model<- function(focal_data, function_name, nparams, lb=-Inf, ub=Inf) {
	par=rep(1, nparams)
	if(any(is.finite(ub))) {
		par[is.finite(ub)] <- ub[is.finite(ub)]
	}
	model_distance <- function(par, focal_data) {
		predictions <- function_name(par, focal_data)
		difference <- sum((focal_data$log_rate - predictions)^2)
		if(!is.finite(difference)) {
			difference <- 1e10
		}
		neglnL <- 0.5*nrow(focal_data)*log(difference) #yes, see lnL at https://en.wikipedia.org/wiki/Akaike_information_criterion#Comparison_with_least_squares, which is -0.5*n*log(RSS), so we get rid of the negative sign
		return(neglnL)
	}
	#return(optim(par=par, fn=model_distance, df=df, lower=lb, upper=ub, method="L-BFGS-B"))
	result <- nloptr::nloptr(x0=par, eval_f=model_distance,  lb=lb, ub=ub, opts=list(algorithm="NLOPT_LN_SBPLX"), focal_data=focal_data)
	
	# starting with lower param values, since they're often small
	par2 <- c(0.1, 0.0001, 1)[1:nparams]
	if(any(is.finite(ub))) {
		par2[is.finite(ub)] <- ub[is.finite(ub)]
	}
	
	result2 <- nloptr::nloptr(x0=par2, eval_f=model_distance,  lb=lb, ub=ub, opts=list(algorithm="NLOPT_LN_SBPLX"), focal_data=focal_data)
	
	if(result2$objective < result$objective) {
		result <- result2
	}
	names(result$solution) <- c("e", "k", "a")[1:length(result$solution)]
	dentist_result <- dent_walk(par=result$solution, fn=model_distance, best_neglnL=result$objective, lower_bound=lb, upper_bound=ub, print_freq=1e6, focal_data=focal_data)
	result$dentist_result <- dentist_result
	return(result)
}




optimization_over_all_data <- function(all_data) {
	datasets <- unique(all_data$citation)
	all_data$time <- all_data$my
	all_data$log_rate <- log(all_data$rate)
	all_data$log_time <- log(all_data$time)
	results <- list()
	for(dataset in datasets) {
		focal_data <- subset(all_data, citation==dataset)
		lb <- -Inf
		ub <- Inf
		names(lb) <- c("e")
		names(ub) <- c("e")
		local_result <- optimize_rate_model(focal_data, function_f1, nparams=1, lb=lb, ub=ub)
		local_result$model <- "f1"
		local_result <- summarize_model(local_result, focal_data, function_f1)
		results[[length(results)+1]] <- local_result
		names(results)[length(results)] <- paste0(dataset, "_", local_result$model)

		#f2 and f3
		
		param_names <- c("e", "k", "a")
		for(select_count in sequence(length(param_names))) {
			print(paste0("select_count=", select_count))
			combinations <- combn(param_names, select_count)
			for (selected_index in sequence(ncol(combinations))) {
				selected_params <- combinations[,selected_index]
				lb=c(0, rep(-Inf, length(param_names)-1))
				ub=rep(Inf, length(param_names))
				names(lb) <- param_names
				names(ub) <- param_names
				tiny <- 0
				lb[!(param_names %in% selected_params)] <- 1-tiny
				ub[!(param_names %in% selected_params)] <- 1+tiny
				if(!("a" %in% selected_params)) {
					local_result <- optimize_rate_model(focal_data, function_f2_general, nparams=2, lb=lb, ub=ub)
					local_result$model <- paste0("f2_", paste0(selected_params, collapse=""))
					local_result <- summarize_model(local_result, focal_data, function_f2_general)
					results[[length(results)+1]] <- local_result
					names(results)[length(results)] <- paste0(dataset, "_", local_result$model)
					
					local_result <- optimize_rate_model(focal_data, function_f3_general, nparams=2, lb=lb, ub=ub)
					local_result$model <- paste0("f3_", paste0(selected_params, collapse=""))
					local_result <- summarize_model(local_result, focal_data, function_f3_general)
					results[[length(results)+1]] <- local_result
					names(results)[length(results)] <- paste0(dataset, "_", local_result$model)
				} else {
					local_result <- optimize_rate_model(focal_data, function_f4_general, nparams=3, lb=lb, ub=ub)
					local_result$model <- paste0("f4_", paste0(selected_params, collapse=""))
					local_result <- summarize_model(local_result, focal_data, function_f4_general)
					results[[length(results)+1]] <- local_result
					names(results)[length(results)] <- paste0(dataset, "_", local_result$model)	
				}	
			}	
		}
	}	
	return(results)
}

summarize_model <- function(local_result, focal_data, function_name) {
	local_result$n <- nrow(focal_data)
	local_result$AIC <- nrow(focal_data)*local_result$objective + 2*length(local_result$solution)
	local_result$numerator <- focal_data$numerator
	local_result$denominator <- focal_data$denominator
	local_result$total_time <- focal_data$total_time
	local_result$my <- focal_data$my
	solution <- local_result$solution
	solution_nomserr <- solution
	solution_nomserr[1] <- 0
	local_result$predicted_log_rate <- function_name(local_result$solution, focal_data)
	local_result$predicted_nonlog_rate <- function_name(local_result$solution, focal_data, dolog=FALSE)
	local_result$empirical_log_rate <- focal_data$log_rate
	local_result$empirical_nonlog_rate <- exp(focal_data$log_rate)
	local_result$predicted_log_rate_no_mserr <- rep(NA, length(local_result$empirical_log_rate))
	local_result$predicted_nonlog_rate_no_mserr <- rep(NA, length(local_result$empirical_log_rate))
	try({ local_result$predicted_log_rate_no_mserr <- function_name(solution_nomserr, focal_data)})
	try({ local_result$predicted_nonlog_rate_no_mserr <- function_name(solution_nomserr, focal_data, dolog=FALSE)})
	local_result$error_only_log_rate <- local_result$predicted_log_rate - local_result$predicted_log_rate_no_mser
	local_result$error_only_nonlog_rate <- local_result$predicted_nonlog_rate - local_result$predicted_nonlog_rate_no_mser
	parameters_no_epsilon <- local_result$par
	#parameters_no_epsilon[1] <- 0
	#local_result$predicted_log_rate_no_mserr <- function_name(parameters_no_epsilon, df)
	local_result$log_time <- focal_data$log_time
	return(local_result)
}

summarize_all_fitted_models <- function(minimization_approach_result) {
	results <- data.frame()
	for (focal_model in names(minimization_approach_result)) {
		data_name <- strsplit(focal_model, "_f")[[1]][1]
		focal_model_suffix <- strsplit(focal_model, "_f")[[1]][2]
		focal_result <- minimization_approach_result[[focal_model]]
		params <- rep(NA,3)
		solution <- focal_result$solution
		params[1:length(solution)] <- solution
		names(params) <- c("e", "k", "a")
		if(ncol(focal_result$dentist_result$all_ranges)<3) { # pad to handle only getting 1 or 2 params
			focal_result$dentist_result$all_ranges <- cbind(focal_result$dentist_result$all_ranges, rep(NA, nrow(focal_result$dentist_result$all_ranges)))
			focal_result$dentist_result$all_ranges <- cbind(focal_result$dentist_result$all_ranges, rep(NA, nrow(focal_result$dentist_result$all_ranges)))
		}
		focal_df <- data.frame(dataset=data_name, model=focal_model_suffix, n=focal_result$n, AIC=focal_result$AIC, nparams=length(focal_result$par), param_e=params['e'], param_k=params['k'], param_a=params['a'], param_e_lower = focal_result$dentist_result$all_ranges['lower.CI', 1], param_e_upper =  focal_result$dentist_result$all_ranges['upper.CI', 1], param_k_lower = focal_result$dentist_result$all_ranges['lower.CI', 2], param_k_upper =  focal_result$dentist_result$all_ranges['upper.CI', 2], param_a_lower = focal_result$dentist_result$all_ranges['lower.CI', 3], param_a_upper =  focal_result$dentist_result$all_ranges['upper.CI', 3], predicted_log_rate=focal_result$predicted_log_rate, empirical_log_rate=focal_result$empirical_log_rate, predicted_log_rate_no_mserr=focal_result$predicted_log_rate_no_mserr, error_only_log_rate=focal_result$error_only_log_rate, log_time=focal_result$log_time, numerator=focal_result$numerator, total_time=focal_result$total_time, my=focal_result$my, denominator=focal_result$denominator)
		focal_df_tall <- focal_df |> tidyr::pivot_longer(cols=c("predicted_log_rate", "empirical_log_rate", "predicted_log_rate_no_mserr", "error_only_log_rate"), names_to="rate_type", values_to="log_rate")
		#focal_df_tall <- focal_df_tall |> tidyr::pivot_longer(cols=c("predicted_nonlog_rate", "empirical_nonlog_rate", "predicted_nonlog_rate_no_mserr", "error_only_nonlog_rate"), names_to="rate_type", values_to="nonlog_rate")
		results <- rbind(results, focal_df_tall)
	}
	results$deltaAIC <- NA
	for (focal_dataset in unique(results$dataset)) {
		focal_rows <- which(results$dataset==focal_dataset)
		results$deltaAIC[focal_rows] <- results$AIC[focal_rows] - min(results$AIC[focal_rows])
	}
	return(results)	
}

summarize_all_models <- function(minimization_approach_result_summarized) {
	models <- minimization_approach_result_summarized |> dplyr::select(-rate_type) |> dplyr::select(-log_rate) |> dplyr::select(-log_time) |> dplyr::select(-nparams) |> dplyr::distinct()
	models <- models[order(models$dataset, models$deltaAIC),]
	return(models)
}

do_randomization_approach_idea <- function(all_data, nreps=5) {
	final_result <- data.frame()
	for(rep_index in sequence(nreps)) {
		local_result <- summarize_all_fitted_models(optimization_over_all_data(randomize_within_dataset(all_data)))
		local_result$rep <- rep_index
		final_result <- rbind(final_result, local_result)	
	}
	return(final_result)
}

merge_original_and_random <- function(minimization_models_summarized, randomized_data_models_summarized) {
	minimization_models_summarized$rep <- "Original"
	randomized_data_models_summarized$rep <- paste0("Rep ", randomized_data_models_summarized$rep)
	merged <- dplyr::bind_rows(minimization_models_summarized, randomized_data_models_summarized)
	return(merged)
}

merge_best_original_and_random <- function(minimization_models_best, randomized_data_models_best) {
	minimization_models_best$rep <- "Original"
	randomized_data_models_best$rep <- paste0("Rep ", randomized_data_models_best$rep)
	merged <- dplyr::bind_rows(minimization_models_best, randomized_data_models_best)
	return(merged)
}

get_unique_compared_to_original <- function(hyperr8_analysis) {
	hyperr8_analysis_best <- hyperr8_analysis |> dplyr::filter(model=="hmb") |> dplyr::group_by(dataset, model, rep, n, nfreeparams, param_h, param_b, param_m, param_h_lower, param_h_upper, param_m_lower, param_m_upper, param_b_lower, param_b_upper) |> dplyr::summarize(hyperbolic_component_proportion=mean(hyperbolic_component_proportion), linear_component_proportion=mean(linear_component_proportion), constant_component_proportion=mean(constant_component_proportion)) |>  dplyr::distinct(dataset, model, n, nfreeparams, param_h, param_m, param_b, param_h_lower, param_h_upper, param_m_lower, param_m_upper, param_b_lower, param_b_upper, rep, hyperbolic_component_proportion, linear_component_proportion, constant_component_proportion) |> dplyr::ungroup()
	return(hyperr8_analysis_best)
}

compute_percentiles <- function(hyperr8_unique_compare) {
	datasets <- unique(hyperr8_unique_compare$dataset)
	results <- data.frame()
	for (focal_dataset in datasets) {
		focal_data <- as.data.frame(subset(hyperr8_unique_compare, dataset==focal_dataset))
		original <- subset(focal_data, rep=="Original")
		reps <- subset(focal_data, rep!="Original")
		focal_columns <- colnames(hyperr8_unique_compare)
		focal_columns <- focal_columns[!focal_columns %in% c("dataset", "rep", "n", "nfreeparams", "model")]
		for (focal_column_name in focal_columns) {
			reps_value <- reps[,focal_column_name]
			original_value <- original[,focal_column_name]
			local_df <- data.frame(dataset=focal_dataset, parameter=focal_column_name, original_percentile=ecdf(reps_value)(original_value), original=original_value, rep_mean=mean(reps_value), rep_sd=sd(reps_value), rep_median=median(reps_value), rep_min=min(reps_value), rep_max=max(reps_value))
			results <- rbind(results, local_df)
			
		}
	}	
	return(results)
}

compute_per_datum_percentiles <- function(hyperr8_analysis) {
	print('starting')
	flush.console()
	hyperr8_analysis$empirical_percentile <- NA
	for (row_index in sequence(nrow(hyperr8_analysis))) {
		if(hyperr8_analysis$rep[row_index]=="Original") {
			randomized <- subset(hyperr8_analysis, datum_id==hyperr8_analysis$datum_id[row_index] & rep!="Original" & model==hyperr8_analysis$model[row_index])$empirical_rate
			original_value <- hyperr8_analysis$empirical_rate[row_index]
			hyperr8_analysis$empirical_percentile[row_index] <- ecdf(randomized)(original_value)
			cat(paste0("\r    ", row_index))
			flush.console()
		}
	}
	return(results)
}

compute_per_datum_percentiles_faster <- function(hyperr8_analysis) {
	reps <- subset(hyperr8_analysis, rep!="Original")
	original <- subset(hyperr8_analysis, rep=="Original")
	ecdf_reps <- reps|> group_by(datum_id, model) |> summarise(ecdfFun = list(ecdf(empirical_rate))) |> ungroup()
	results <- merge(original, ecdf_reps, by=c("datum_id", "model"))
	results$empirical_percentile <- results$ecdfFun |> purrr::map2_dbl(results$empirical_rate, ~ .x(.y))
	results$ecdfFun <- NULL
	results <- dplyr::bind_rows(results, reps)
	return(results)
}

summarize_per_datum_percentiles <- function(hyperr8_datum_percentiles) {
	best <- subset(hyperr8_datum_percentiles, rep=="Original" & deltaAIC==0)	
	results <- data.frame()
	datasets <- unique(best$dataset)
	for (focal_dataset in datasets) {
		focal_df <- subset(best, dataset==focal_dataset)	
		quantile_results <- quantile(focal_df$empirical_percentile, seq(from=0, to=1, length.out=21))
		results <- rbind(results, data.frame(dataset=focal_dataset, quantile=names(quantile_results), value=quantile_results))
	}
	results$quantile <- 0.01*as.numeric(gsub('%', '', results$quantile))
	rownames(results) <- NULL
	return(results)
}

filter_extrema <- function(all_data, min_percentile=0, max_percentile=1, min_age=0, max_age=Inf) {
	result <- data.frame()
	datasets <- unique(all_data$citation)
	for (focal_dataset in datasets) {
		focal_data <- subset(all_data, citation==focal_dataset)
		focal_data$cdf_times <- ecdf(focal_data$time)(focal_data$time)
		focal_data <- subset(focal_data, time>=min_age & time<=max_age & cdf_times>=min_percentile & cdf_times<=max_percentile)
		result <- rbind(result, focal_data)
	}
	return(result)
}

create_gganimate_plot <- function(hyperr8_analysis) {
	#short <- subset(hyperr8_analysis, n<100)
	short <- hyperr8_analysis
	
	short$rep_number <- as.numeric(gsub("Rep ", "", short$rep))
	
	short$rep_number[is.na(short$rep_number)] <- 0
	reps <- subset(short, rep_number>0)

	short <- subset(short, rep_number<=5)
	original <- subset(short, rep_number==0)
	
	potential_purge <- unique(subset(short, dataset=="Pure birth simulation")$datum_id)
	potential_purge <- sample(potential_purge, length(potential_purge)-5000, replace=FALSE)
	short <- subset(short, !(datum_id %in% potential_purge))
	
	short$grouptime <- as.factor(short$time)
	bounds <- reps |> dplyr::group_by(dataset, datum_id) |> dplyr::summarize(lower=quantile(empirical_rate, 0.025), upper=quantile(empirical_rate, 0.975), time=mean(time), hyperbolic_component_proportion=mean(hyperbolic_component_proportion)) |> dplyr::ungroup()
	
	g <- ggplot(short, aes(x=time, y=empirical_rate, colour=hyperbolic_component_proportion)) +  geom_point(data=short, aes(group=datum_id)) + scale_colour_gradient(low="blue", high="red") + facet_wrap(~dataset, scales="free") + theme_bw() + theme(legend.position="none") + xlab("h") + ylab("m") + scale_x_continuous(trans="log") +  labs(title = 'Rep: {round(frame_time)}', x = 'Time', y = 'Rate') + scale_y_continuous(trans="log") + transition_time(rep_number) +  enter_fade() + exit_fade() + ease_aes('cubic-in-out')  +  geom_smooth(data=bounds, aes(x=time,y=upper), se=FALSE, colour='black', lty="dashed", lwd=0.8) + geom_smooth(data=bounds, aes(x=time,y=lower), se=FALSE, colour='black', lty="dashed", lwd=0.8) 
	

	
	animate(g, height = 6, width = 12, units = "in", res = 150)	
  	anim_save(file="~/Downloads/hyperr8_animation.gif")
	system("open ~/Downloads/hyperr8_animation.gif")
}

subtract_medians <- function(hyperr8_analysis) {
	hyperr8_analysis$rep_number <- as.numeric(gsub("Rep ", "", hyperr8_analysis$rep))
	hyperr8_analysis$rep_number[is.na(hyperr8_analysis$rep_number)] <- 0
	original <- subset(hyperr8_analysis, rep_number==0)
	reps <- subset(hyperr8_analysis, rep_number>0)
	reps_aggregated <- reps |> group_by(datum_id) |> summarize(median_randomized_rate=median(empirical_rate), mean_randomized_rate=mean(empirical_rate)) |> ungroup()	
	original <- merge(original, reps_aggregated, by="datum_id")
	original$empirical_rate_minus_median_simulated <- original$empirical_rate - original$median_randomized_rate
	original$empirical_rate_minus_mean_simulated <- original$empirical_rate - original$mean_randomized_rate

	return(original)
}

plot_medians <- function(hyperr8_analysis_subtracted_medians) {
	g <- ggplot(hyperr8_analysis_subtracted_medians, aes(x=time, y=empirical_rate)) + geom_point(alpha=0.02) + scale_x_continuous(trans="log") + facet_wrap(~dataset, scales="free") + geom_smooth()
	return(g)
}

compute_ks_result <- function(hyperr8_analysis) {
	datasets <- unique(hyperr8_analysis$dataset)
	results <- data.frame()
	for (focal_dataset in datasets) {
		print(focal_dataset)
		focal_data <- subset(hyperr8_analysis, dataset==focal_dataset)
		original <- subset(focal_data, rep=="Original" & deltaAIC==0)
		filtered_data <- subset(focal_data, model==original$model[1])
		instances <- unique(filtered_data$rep)
		for (focal_instance in instances) {
			focal_rep <- subset(filtered_data, rep==focal_instance)
			other_rep <- subset(filtered_data, rep!=focal_instance)
			ks_result <- ks.test(focal_rep$empirical_rate, other_rep$empirical_rate)
			results <- rbind(results, data.frame(dataset=focal_dataset, focal_rep=focal_instance, D=ks_result$statistic, p=ks_result$p.value))
		}
	}
	rownames(results) <- NULL
	return(results)
}

# Thanks to Naomi O'Meara for the idea. This is also known as R^2
compute_coefficient_of_determination <- function(hyperr8_analysis) {
	results <- hyperr8_analysis |> group_by(dataset, rep, model, deltaAIC, param_h, param_m, param_b, param_h_lower, param_h_upper, param_m_lower, param_m_upper, param_b_lower, param_b_upper, n) |> summarize(r2=compute_coefficient_of_determination_with_lm(empirical_rate, predicted_rate)) |> ungroup()
	return(results)	
}


compute_coefficient_of_determination_with_lm <- function(observed, predicted) {
	log_observed <- log(observed)
	log_predicted <- log(predicted)
	if(any(observed==0)) {
		log_observed <- log1p(observed)
		log_predicted <- log1p(predicted)
	}
	return(summary(lm(log_observed~log_predicted))$r.squared)
}	

merge_r2_tables <- function(r2_results_random_params_summarized, r2_results_pretty) {
	#r2_results_random_params_summarized$dataset_model_rep <- paste0(r2_results_random_params_summarized$dataset, "_", r2_results_random_params_summarized$model, "_", r2_results_random_params_summarized$rep)
	
	#r2_results_pretty$dataset_model_rep <- paste0(r2_results_pretty$dataset, "_", r2_results_pretty$model, "_", r2_results_pretty$rep)
	merged <- dplyr::left_join(r2_results_random_params_summarized, r2_results_pretty, by=c("dataset", "model"))
	merged <- dplyr::select(merged, -r2.y) # due to rounding, these don't match
	colnames(merged) <- gsub("\\.x", "", colnames(merged))
	return(merged)
}

prettily_summarize_coefficient_of_determination <- function(r2_results, AIC_threshold=10) {
	results <- data.frame()
	r2_results <- r2_results |> dplyr::arrange(desc(n))
	datasets <- unique(r2_results$dataset)
	models <- unique(r2_results$model)
	for (focal_dataset in datasets) {
		for (focal_model in models) {
			focal_data <- subset(r2_results, dataset==focal_dataset & model==focal_model)
			original <- subset(focal_data, rep=="Original")	
			other_reps <- subset(focal_data, rep!="Original")
			final <- data.frame(
				dataset=focal_dataset,
				model=focal_model,
				deltaAIC=round(original$deltaAIC[1],2),
				r2=original$r2[1],
				r2_percentile_relative_to_randomized = NA,

				h = paste0(signif(original$param_h[1], 2), " (", signif(original$param_h_lower[1], 2), ", ", signif(original$param_h_upper[1], 2), ")"),
				h_percentile_relative_to_randomized = ecdf(other_reps$param_h)(original$param_h[1]),

				m = paste0(signif(original$param_m[1], 2), " (", signif(original$param_m_lower[1], 2), ", ", signif(original$param_m_upper[1], 2), ")"),
				m_percentile_relative_to_randomized = ecdf(other_reps$param_m)(original$param_m[1]),

				b = paste0(signif(original$param_b[1], 2), " (", signif(original$param_b_lower[1], 2), ", ", signif(original$param_b_upper[1], 2), ")"),
				b_percentile_relative_to_randomized = ecdf(other_reps$param_b)(original$param_b[1]),
				n=original$n[1]
			)
			try({final$r2_percentile_relative_to_randomized <- ecdf(other_reps$r2)(original$r2[1])}, silent=TRUE)
			results <- rbind(results, final)
		}	
	}
	#results <- subset(results, (model %in% c("hmb", "h00")) | deltaAIC<AIC_threshold)
	results <- results |> dplyr::arrange(desc(n), deltaAIC) |> dplyr::select(-n)
	return(results)
}

r2_from_prediction_from_randomizations <- function(hyperr8_analysis) {
	results <- data.frame()
	original <- subset(hyperr8_analysis, rep=="Original")
	datasets <- unique(original$dataset)
	models <- unique(original$model)
	for (focal_dataset in datasets) {
		print(focal_dataset)
		for (focal_model in models) {
			print(focal_model)
			focal_original <- subset(original, dataset==focal_dataset & model==focal_model)
			focal_all<- subset(hyperr8_analysis, dataset==focal_dataset & model==focal_model)
			for (focal_rep in unique(focal_all$rep)) {
				focal_rep_data <- subset(focal_all, rep==focal_rep)
				h <- focal_rep_data$param_h[1]
				m <- focal_rep_data$param_m[1]
				b <- focal_rep_data$param_b[1]
				focal_original$predicted_rate <- h/focal_original$time + m*focal_original$time + b
				r2 <- compute_coefficient_of_determination_with_lm(focal_original$empirical_rate, focal_original$predicted_rate)
				results <- rbind(results, data.frame(dataset=focal_dataset, model=focal_model, rep=focal_rep, r2=r2))
			}
		}
	}
	return(results)
}

summarize_r2_from_prediction_from_randomizations <- function(r2_results_random_params) {
	focal_models <- unique(r2_results_random_params$model)
	results <- data.frame()
	datasets <- unique(r2_results_random_params$dataset)
	for (focal_dataset in datasets) {
		for (focal_model in focal_models) {
			focal_data <- subset(r2_results_random_params, dataset==focal_dataset & model==focal_model)
			original <- subset(focal_data, rep=="Original")
			other_reps <- subset(focal_data, rep!="Original")
			final <- data.frame(
				dataset=focal_dataset,
				model=focal_model,
				r2=original$r2[1],
				r2_randomized_mean=NA,
				r2_randomized_sd=NA,
				r2_randomized_min=NA,
				r2_randomized_max=NA,
				r2_mean_amount_worse_with_randomized = NA,
				r2_min_amount_worse_with_randomized = NA,
				r2_max_amount_worse_with_randomized = NA
			)
			try({final$r2_randomized_mean <- mean(other_reps$r2, na.rm=TRUE)}, silent=TRUE)
			try({final$r2_randomized_sd <- sd(other_reps$r2, na.rm=TRUE)}, silent=TRUE)
			try({final$r2_randomized_min <- min(other_reps$r2, na.rm=TRUE)}, silent=TRUE)
			try({final$r2_randomized_max <- max(other_reps$r2, na.rm=TRUE)}, silent=TRUE)
			try({final$r2_mean_amount_worse_with_randomized <- mean( original$r2[1] - other_reps$r2, na.rm=TRUE)}, silent=TRUE)
			try({final$r2_min_amount_worse_with_randomized <- min( original$r2[1] - other_reps$r2, na.rm=TRUE)}, silent=TRUE)
			try({final$r2_max_amount_worse_with_randomized <- max( original$r2[1] - other_reps$r2, na.rm=TRUE)}, silent=TRUE)
			results <- rbind(results, final)
		}
	}
	return(results)
}

save_file_in_chunks <- function(hyperr8_analysis) {
	datasets <- unique(hyperr8_analysis$dataset)
	for (focal_dataset in datasets) {
		focal_data <- subset(hyperr8_analysis, dataset==focal_dataset)
		cleaned_name <- gsub("[[:space:].()]", "_", focal_dataset)
		write.csv(focal_data, file=gzfile(paste0("outputs/", cleaned_name, ".csv.gz")))
	}
}

