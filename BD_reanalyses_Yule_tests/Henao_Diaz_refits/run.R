
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
### Making sure likelihood code is consistent with other implementations
######################################################################################################################################
######################################################################################################################################

library(hisse)
phy <- read.tree("whales_Steemanetal2009.tre")
opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = 100000, "ftol_rel" = .Machine$double.eps^.50)
obj.stad <- nloptr(x0=log(c(0.5,.5)), eval_f=GetLikelihood, ub=log(c(500,1.2)), lb=c(-100,-100), opts=opts, tree=phy, condition.type="survival", verbose=FALSE, log.convert=TRUE, par.class=c("turn", "ef"), rho=1)

obj.misse <- MiSSE(phy, f=1, turnover=c(1), eps=c(1), sann=FALSE, starting.vals=c(0.5, 0.5, 1))

comparison <- identical(round(-obj.stad$objective,4), round(obj.misse$loglik,4))
comparison

#Should be true -- it is 8/29/23


######################################################################################################################################
######################################################################################################################################
### Looping over the same set of trees from Henao Diaz et al (2019) getting MLE and 95% CI -- BIRTH-DEATH MODEL
######################################################################################################################################
######################################################################################################################################

# # organize meta data so that we can know the rho for each tree
all_files <- dir("henaodiaz_etal_trees-master/", full.names = TRUE)
tree_files <- all_files[-c(1,2)]
meta_files <- all_files[c(1,2)]
sample_size <- read.csv(meta_files[1])
meta_data <- read.csv(meta_files[2])
sample_size <- sample_size[match(meta_data$taxon.names, sample_size$Clade.name),]
index_tree_file <- unlist(sapply(meta_data$study.code, function(x) grep(paste0("_", x), tree_files)))
sf_tree_data <- data.frame(sample_size[,-5], study_code = meta_data$study.code, tree_file = tree_files[index_tree_file])
tree_row <- sf_tree_data[2,]


# # based on stadler 2013
get_mle_pars <- function(tree_row, ip=c(0.5, 0.5)){
   phy <- read.tree(tree_row[6])
   sf <- length(phy$tip.label)/as.numeric(tree_row[4])
   opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = 100000, "ftol_rel" = .Machine$double.eps^.50)
   obj <- nloptr(x0=log(ip), eval_f=GetLikelihood, ub=log(c(500,.99)), lb=c(-100,-100), opts=opts, tree=phy, condition.type="survival", verbose=FALSE, log.convert=TRUE, par.class=c("turn", "ef"), rho=sf)
   return(obj)
}


# A function to convert between parameters
quickConvert <- function(par, par.class){
	base_p <- c(lambda=NA, mu=NA, net.div=NA, turn=NA, ef=NA)
	base_p[match(par.class, names(base_p))] <- par
	p <- convertBetweenPars(base_p)
	names(p) <- c("lambda", "mu", "net.div", "turn", "ef")
	return(p)
}

quickCIconvert <- function(ci_row){
	tmp <- rbind(
		quickConvert(ci_row[c(1,3)], par.class = c("turn", "ef")),
		quickConvert(ci_row[c(1,4)], par.class = c("turn", "ef")),
		quickConvert(ci_row[c(2,3)], par.class = c("turn", "ef")),
		quickConvert(ci_row[c(2,4)], par.class = c("turn", "ef"))
		)
	min_est <- apply(tmp, 2, function(x) x[which.min(x)])
	names(min_est) <- paste0(names(min_est), "_lower")
	max_est <- apply(tmp, 2, function(x) x[which.max(x)])
	names(max_est) <- paste0(names(max_est), "_upper")
	
	return(c(min_est, max_est))
}

# # fit birthdeath models to the trees
# # fits <- lapply(trees, function(x) birthdeath)
fits <- apply(sf_tree_data, 1, get_mle_pars)

# # extract the parameters
params <- do.call(rbind, lapply(fits, function(x) c(x$objective, quickConvert(exp(x$solution), c("turn", "ef")))))
colnames(params)[1] <- "neglnL"
df_total <- data.frame(sf_tree_data, params)

# # get the CI
# # ci <- do.call(rbind, lapply(fits, function(x) x$CI[2,]))
# # tree_row <- unlist(c(df_total[1,]))
dent_ci <- function(tree_row){
	phy <- read.tree(tree_row[6])
	sf <- length(phy$tip.label)/as.numeric(tree_row[4])
	par <- setNames(as.numeric(tree_row[8:9]), c("turn", "ef"))
	out <- dent_walk(par, GetLikelihood, best_neglnL = as.numeric(tree_row[7]), tree=phy, nsteps = 2000, condition.type="survival", verbose=FALSE, log.convert=FALSE, par.class=c("turn", "ef"), rho=sf, sd = c(2,1))
	return(out)
}
ci <- apply(df_total, 1, dent_ci)

# # extract the ages of each tree
df_ci <- do.call(rbind, lapply(ci, function(x) unlist(c(x$all_ranges[2:3,]))))
colnames(df_ci) <- c("turn_lower", "turn_upper", "ef_lower", "ef_upper")
df_ci <- t(apply(df_ci, 1, quickCIconvert))
df_plot <- data.frame(df_total, df_ci)
write.csv(df_plot, file = "01-bd_df_plot.csv", row.names = FALSE)


######################################################################################################################################
######################################################################################################################################
### Looping over the same set of trees from Henao Diaz et al (2019) getting MLE and 95% CI -- YULE MODEL
######################################################################################################################################
######################################################################################################################################

# # written by Jeremy Danger Beaulieu

get_mle_pars_yule <- function(phy, sf){
	opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = 100000, "ftol_rel" = .Machine$double.eps^.50)
	obj <- optimize(GetLikelihoodYule, interval=c(0,5),  tree=phy, rho=sf, condition.type="survival", verbose=FALSE, tol=0.00000001)
	return(obj)
}


dent_ci_yule <- function(tree_row){
	phy <- read.tree(tree_row[6])
	sf <- length(phy$tip.label)/as.numeric(tree_row[4])
	mle <- get_mle_pars_yule(phy=phy, sf=sf)
	ll <- dent_walk(par=mle$minimum, fn=GetLikelihoodYule, best_neglnL=mle$objective, tree=phy, rho=sf, condition.type="survival")
	obj <- list(para=mle$minimum, neglnL=mle$objective, CI=ll$all_ranges[c(2,3),1])
}

fits <- apply(sf_tree_data, 1, dent_ci_yule)
lambda <- do.call(rbind, lapply(fits, function(x) x$para))
ci <- do.call(rbind, lapply(fits, function(x) x$CI))
neglnL <- do.call(rbind, lapply(fits, function(x) x$neglnL))
df_plot <- data.frame(sf_tree_data, neglnL=neglnL, lambda, ci)
write.csv(df_plot, file = "02-yule_df_plot.csv", row.names = FALSE)
