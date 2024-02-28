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


######################################################################################################################################
######################################################################################################################################
### Examining ascertainment bias and general estimator bias for Yule birth rate
######################################################################################################################################
######################################################################################################################################

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

ProbN_lt <- function(n, lambda, time){
  tmp <- (n - 1) * exp(-2*lambda*time) * ((1 - exp(-lambda*time))^(n - 2))
  return(tmp)
}

ProbN_G3 <- function(n, lambda, time){
  tmp <- c()
  for(index in 2:n){
	tmp <- rbind(tmp, ProbN_lt(index, lambda, time))
  }
  return (1 - colSums(tmp))
}

yuleMLE <- function(phy){
  mle <- (Ntip(phy)-2) / sum(phy$edge.length)
  return(mle)
}

# These columns represent lambda, mu, and mrca time
ef <- c(0)
turn <- c(0.1)
ages <- seq(from=1, to=50, length.out=1000)

lhs_table <- expand.grid(ef, turn, ages)

# Plot n_t as a function of time given ef and turn
# Convert the table into lambda and mu
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
							  pars_table[1,2], time_limit = 100), mc.cores = 4)
#save(trees, pars_table, file="yule_trees_bias.Rsave")

fits <- lapply(trees, yule)

# extract the parameters
params <- do.call(rbind, lapply(fits, function(x) x$lambda))

# extract the ages of each tree
ages <- unlist(lapply(trees, function(x) max(branching.times(x))))
ntips <- unlist(lapply(trees, function(x) length(x$tip.label)))
params.unbiased <- params * ((ntips)/(ntips-1))
params.unbiased[is.na(params.unbiased)] <- 0


# organize the plot data
yule_plot_data <- data.frame(params, params.unbiased, ages, ntips)
yule_plot_data_ab <- yule_plot_data[yule_plot_data$params>0,]
save(yule_plot_data, yule_plot_data_ab, file="07-yule_age_data_correction.Rsave")



######################################################################################################################################
######################################################################################################################################
### Plot code for fig. S8
######################################################################################################################################
######################################################################################################################################


######################################################################################################################################
### 1 x 1 panel -- Three different ways of plotting the line.
######################################################################################################################################
load("07-yule_age_data_correction.Rsave")
p_empirical_yule <- ggplot(yule_plot_data, aes(x = ages, y = params)) +
	geom_point(size=.75, alpha=.3) +
	theme_bw() +
	theme(panel.grid = element_blank()) +
	coord_cartesian(ylim = c(0, .4)) +
	geom_hline(yintercept=0.1, lty=3) +
	geom_smooth(aes(x = ages, y = params.unbiased), colour="#F0F921FF") +
	geom_smooth(aes(x = ages, y = params), colour="#9C179EFF") +
	geom_smooth(data=yule_plot_data_ab, aes(x = ages, y = params.unbiased), colour="#0D0887FF") +
	annotate(geom = 'text', label = c(expression(MLE[censored+biased]), expression(MLE[uncensored+biased]), expression(MLE[uncensored+unbiased])), x = c(30,30,30), y = c(.4,.36,.32), hjust = 0, vjust = 1, size=3) +
	annotate(geom="segment", x=c(29,29,29), y=c(.4,.36,.32), xend=c(27,27,27), yend=c(.4,.36,.32), colour=c("#0D0887FF", "#9C179EFF", "#F0F921FF")) +
	labs(y=expression(lambda~(Myr^-1)), x="Clade Age (Myr)")

out <- gridExtra::grid.arrange(p_empirical_yule, NULL, NULL, NULL, NULL, nrow = 3, ncol=2)
ggsave(file = "FigsS8.pdf", plot = out, width = 7.24, height = 7.24)

