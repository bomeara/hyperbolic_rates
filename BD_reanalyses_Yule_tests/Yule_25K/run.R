library(TreeSim)
library(ape)
library(parallel)
library(ggplot2)
library(gridExtra)
library(dentist)
library(nloptr)
library(magrittr)
library(tidyverse)
library(viridis)
library(devtools)
library(latex2exp)

# load the results
source("Utils.R")


# Setting up the tree simulation
timed_sim.bd.age <- function(age, lambda, mu, time_limit=100) {
    setTimeLimit(time_limit)
    tree <- try(sim.bd.age(
        age = age,
        numbsim = 1,
        lambda = lambda,
        mu = mu,
        mrca = TRUE, complete=FALSE)[[1]])
    return(tree)
}

get_root_age <- function(phy){
    root_age <- max(branching.times(phy))
    return(root_age)
}

plot_nt <- function(ef, turn){
    pars <- quickConvert(c(ef, turn), c("ef", "turn"))
    ages <- seq(0, 25000, 0.1)
    n_t <- sapply(ages, function(x) get_n_t(2, pars[1], pars[2], x))
    plot(ages, log10(n_t), type = "l", xlab = "Time", ylab = "10^Number of tips", main = paste("ef =", ef, "turn =", turn))
}

GetML <- function(phy){
    opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = 100000, "ftol_rel" = .Machine$double.eps^.50)
    obj <- optimize(GetLikelihoodYule, interval=c(0,2),  tree=phy, rho=1, condition.type="survival", verbose=FALSE, tol=0.00000001)
    return(obj)
}

GetCI <- function(phy){
	if(Ntip(phy)==2){
		obj <- list(para=0, CI=c(0,0))
	}else{
		mle <- GetML(phy)
		sink("/dev/null")
		ll <- dent_walk(par=mle$minimum, fn=GetLikelihoodYule, best_neglnL=mle$objective, tree=phy, rho=1, condition.type="survival")
     	   	sink()
           	obj <- list(para=mle$minimum, CI=ll$all_ranges[c(2,3),1])
	}
	return(obj)
}


RunSims <- function(nrep){

	ef <- c(0)
	turn <- c(0.1)
	ages <- exp(c(runif(1000, log(1), log(50))))

	lhs_table <- expand.grid(ef, turn, ages)
	pars_table <- apply(lhs_table, 1, function(x) quickConvert(c(x[1], x[2]),
                                                        c("ef", "turn")))

	# Merge all the tables together and add age
	pars_table <- cbind(t(pars_table), age = lhs_table[,3])

	# A function to get the number of tips at a given time
	expected_taxa <- apply(pars_table, 1, function(x) get_n_t(2, x[1], x[2], x[6]))
	range(expected_taxa)
	# combine expected taaxa with the table
	pars_table <- cbind(pars_table, E_taxa = expected_taxa)

	# Accept only those with more than 4 taxa and fewer than 1000
	# Accepted_pars_table <- pars_table[pars_table[,7] > 4 & pars_table[,7] < 10000,]

	age_list <- sapply(ages, as.list)

	trees <- mclapply(age_list, function(x) timed_sim.bd.age(x, pars_table[1,1],
                              pars_table[1,2], time_limit = 100), mc.cores = 2)

	# get the mle parameters for each tree -- trees have to be >2 otherwise the rate is zero.
	finished_trees <- trees[unlist(lapply(trees, function(x) class(x) != "try-error"))]
	reconstructed_trees <- lapply(finished_trees, drop.extinct)
	#reconstructed_trees <- reconstructed_trees[unlist(lapply(reconstructed_trees,
                                       #function(x) length(x$tip.label) > 2))]

	fits <- lapply(reconstructed_trees, GetCI)

	params <- do.call(rbind, lapply(fits, function(x) x$para))
	ci <- do.call(rbind, lapply(fits, function(x) x$CI))

	# extract the ages of each tree
	ages <- unlist(lapply(reconstructed_trees, function(x) max(branching.times(x))))
	ntips <- unlist(lapply(reconstructed_trees, function(x) length(x$tip.label)))
	totl <- unlist(lapply(reconstructed_trees, function(x) sum(x$edge.length)))

	# organize the plot data
	yule_plot_data <- data.frame(params, ci, ages, ntips, totl)

	save(reconstructed_trees, yule_plot_data, file="12-yule_trees_bias_25K_uncen.Rsave")
}

RunSims()



