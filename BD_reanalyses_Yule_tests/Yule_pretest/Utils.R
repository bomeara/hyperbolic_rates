

p0 <- function(t,l,m,rho){
    1- rho*(l-m)/(rho*l+(l*(1-rho)-m)*exp(-(l-m)*t))
}


p1 <- function(t,l,m,rho){
    rho*(l-m)^2 * exp(-(l-m)*t)/(rho*l+(l*(1-rho)-m)*exp(-(l-m)*t))^2
}


qhelp <- function(t,l,m,rho){
    rho*l*(1-exp(-(l-m)*t))/(rho*l+(l*(1-rho)-m)*exp(-(l-m)*t))
}

stadler.pn <- function(birth.rate, death.rate, time, n){
    rho=1
    l <- birth.rate
    m <- death.rate
    t <- time
    part1 <- p1(t,l,m,rho)
    part2 <- qhelp(t,l,m,rho)
    pn.t <- part1 * part2^(n-1)
    return(pn.t)
}


#Based on Magallon and Sanderson Eq. 11
AscertainBias1 <- function(k, l, m, t, rho) {
    net.diver.rate <- l - m
    #Magallon and Sanderson 2001 -- Eq. 2a:
    exprt <- exp(net.diver.rate * t)
    beta <- (exprt - 1) / (exprt - (m/l))
    #Magallon and Sanderson 2001 -- Eq. 2b:
    alpha <- (m/l) * beta
    #Magallon and Sanderson 2001 -- Eq. 10a:
    probNgeNtax <- beta^(k-1)
    survival.prob <- (1 - alpha)^2
    probNgeNtax <- (beta^(k-2))*(k*(1 - alpha - beta + alpha*beta) + alpha + 2*beta-1)/(1 - alpha + 2*alpha*beta)
    combined.prob <- probNgeNtax * survival.prob
    return(combined.prob)
}


#Based on Jeremy messing around with Stadler. An approximation because we do not go to infinity (but the differences are trivial).
AscertainBias2 <- function(k, l, m, t, rho, end=10000000, get.partials=FALSE) {
    part1 <- p1(t,l,m,rho)
    part2 <- qhelp(t,l,m,rho)
    res.res <- numeric(10000000)
    for(n.index in k:end){
        res.res[n.index] <- (n.index-1)*(part1^2)*(part2^(n.index-2))
    }
    if(get.partials == TRUE){
        return(res.res)
    }else{
        prob.nplus <- sum(res.res)
        return(prob.nplus)
    }
}


#conditional probability based on survival only
survival.conditional.p <- function(time, turn, eps){
    l <- turn/(1+eps)
    m <- (turn * eps) / (1 + eps)
    conditional.p <- (1-p0(time,l,m,rho=1))^2
    return(conditional.p)
}


#conditional probability based on survival and n taxa
survival.n.conditional.p <- function(time, turn, eps, n){
    l <- turn/(1+eps)
    m <- (turn * eps) / (1 + eps)
    part1 <- p1(time,l,m,rho=1)
    part2 <- qhelp(time,l,m,rho=1)
    conditional.p <- (n-1)*(part1^2)*(part2)^(n-2)
    return(conditional.p)
}


GetLikelihood <- function(par, tree, par.class=c("turn", "ef"), rho=1, condition.type="survival", verbose=FALSE, log.convert=TRUE){
  base_p <- c(lambda=NA, mu=NA, net.div=NA, turn=NA, ef=NA)
  base_p[match(par.class, names(base_p))] <- par
    if(log.convert){
      base_p <- exp(base_p)
    }
    p <- convertBetweenPars(base_p)
    # return(p)
    n <- Ntip(tree)
    l <- p[1]
    m <- p[2]
    x <- getx(tree)
    x <- sort(x,decreasing=TRUE)
    t <- x

    lik <- 0
    for (i in 2:length(t)){
        lik <- lik+log(l*p1(t[i],l,m,rho))
    }
    if(condition.type == "survival"){
        lik.unconditioned <- lik + 2 * log(p1(t[1],l,m,rho))
        condition.p <- 2 * log(1-p0(t[1],l,m,rho))
    }
    if(condition.type == "exactlyn"){
        lik.unconditioned <- lik + 2*log(p1(t[1],l,m,rho))
        part1 <- p1(t[1],l,m,rho)
        part2 <- qhelp(t[1],l,m,rho)
        condition.p <- log((n-1)*(part1^2)*(part2)^(n-2))
    }
    if(condition.type == "beauOmeara"){
        lik.unconditioned <- lik + 2*log(p1(t[1],l,m,rho))
        condition.p <- log(AscertainBias2(n,l,m,t[1],rho=1))
    }
    if(condition.type == "magsan"){
        lik.unconditioned <- lik + 2*log(p1(t[1],l,m,rho))
        condition.p <- log(AscertainBias1(n,l,m,t[1],rho=1))
    }
    if(condition.type == "exactlyn.stem0"){
        lik.unconditioned <- lik + log(l*p1(t[1],l,m,rho))
        condition.p <- log(n * (p1(t[1],l,m,rho) / (1-p0(t[1],l,m,rho))))
    }
    if(condition.type == "exactlyn.stem1"){
        lik.unconditioned <- lik + log(l*p1(t[i],l,m,rho))
        condition.p <- log((qhelp(t=t[1], l=l, m=m, rho=rho))^(1-n))
    }
    if(condition.type == "none") {
      lik.unconditioned <- lik + log(l*p1(t[i],l,m,rho))
      condition.p <- log(1)
    }
    if(verbose==TRUE){
        print(paste("unconditional prob", lik.unconditioned))
        print(paste("conditional prob", condition.p))
    }

    if(condition.type == "exactlyn.stem0" | condition.type == "exactlyn.stem1" | condition.type == "exactlyn" | condition.type == "magsan" | condition.type == "morethann"){
        #tmp <- c(l+m, m/l, lik.unconditioned, log(condition.p))
        lik.conditioned.both.check <- lik.unconditioned - condition.p
        names(lik.conditioned.both.check) <- "llik"
        return(-lik.conditioned.both.check)
    }else{
        lik.conditioned.both.check <- lik.unconditioned - condition.p
        names(lik.conditioned.both.check) <- "llik"
        return(-lik.conditioned.both.check)
    }
}


getCI <- function(ntaxa, turn, frac){
  ip <- convertBetweenPars(c(NA, NA, NA, p))
  phy <- trees(ip[1:2], "bd", max.taxa=ntaxa)[[1]]
  # lik <- make.bd(phy)
  # fit_true <- find.mle(lik, c(.1, .05), method="optim", lower=0)
  opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = 100000, "ftol_rel" = .Machine$double.eps^.50)
  res <- list()
  condition.types <- c("none", "survival", "exactlyn", "magsan")
  for(i in 1:length(condition.types)){
    condition.type <- condition.types[i]
    obj <- nloptr(x0=log(ip[1:2]), eval_f=GetLikelihoodBD, ub=log(c(100,.99)), lb=c(-6,-6), opts=opts, tree=phy, condition.type=condition.type, verbose=FALSE, log.convert=TRUE)
    results <- setNames(exp(obj$solution), c("turn", "frac"))
    dented_results <- dent_walk(par = results, GetLikelihoodBD, best_neglnL = obj$objective, delta = 2, nsteps = 1000, tree=phy, condition.type=condition.type, verbose=FALSE, log.convert=FALSE)
    res[[i]] <- list(opt = obj, phy = phy, dent = dented_results)
  }
  names(res) <- condition.types
  return(res)
}

getCI_multi <- function(ntaxa, turn, frac, n_replicates, n_cores){
  multi_res <- mclapply(seq(n_replicates), function(x) getCI(ntaxa, turn, frac), mc.cores = n_cores)
  return(multi_res)
}

good_melt <- function(dent_res){
  tmp <- dent_res$all_ranges
  tmp <- cbind(measure = rownames(tmp), tmp)
  out <- melt(tmp)
  return(out)
}

getMesTable <- function(focal_rep){
  tmp <- do.call(rbind, lapply(focal_rep, function(x) good_melt(x$dent)))
  condition.type <- gsub("\\..*", "", rownames(tmp))
  rownames(tmp) <- NULL
  out <- as.data.frame(cbind(condition.type=condition.type, tmp))
  return(out)
}

# probabiliity of simulating n species
p_n <- function(l, m, C, n){
  a <- l^(n-1)
  b <- n*C
  c <- 1-exp(-(l-m)*C)
  d <- l-(m*exp(-(l-m)*C))
  p <- (a/b)*((c/d)^n)
  return(p)
}

log_pn <- function(par, l, m, n){
  return(-log(p_n(l, m, par, n)))
}

convert2Lambda <- function(pars){
  if(is.na(pars[1])){
    focal_pars <- sample(which(!is.na(pars)), size = 2, replace = FALSE)
    if(2 %in% focal_pars & 3 %in% focal_pars){
      # mu and div
      lambda <- pars[2] + pars[3]
    }
    if(2 %in% focal_pars & 4 %in% focal_pars){
      # mu and turn
      lambda <- pars[4] - pars[2]
    }
    if(2 %in% focal_pars & 5 %in% focal_pars){
      # mu and ef
      lambda <- pars[2]/pars[5]
    }
    if(3 %in% focal_pars & 4 %in% focal_pars){
      # div and turn
      lambda <- (pars[3]+pars[4])/2
    }
    if(3 %in% focal_pars & 5 %in% focal_pars){
      # div and ef
      lambda <- pars[3]/(1-pars[5])
    }
    if(4 %in% focal_pars & 5 %in% focal_pars){
      # turn and ef
      lambda <- pars[4]/(1+pars[5])
    }
  }else{
    lambda <- pars[1]
  }
  return(lambda)
}

convert2Mu <- function(pars){
  if(is.na(pars[2])){
    focal_pars <- sample(which(!is.na(pars)), size = 2, replace = FALSE)
    if(1 %in% focal_pars & 3 %in% focal_pars){
      # lambda and div
      mu <- pars[1] - pars[3]
    }
    if(1 %in% focal_pars & 4 %in% focal_pars){
      # lambda and turn
      mu <- pars[4] - pars[1]
    }
    if(1 %in% focal_pars & 5 %in% focal_pars){
      # lambda and ef
      mu <- pars[1]*pars[5]
    }
    if(3 %in% focal_pars & 4 %in% focal_pars){
      # div and turn
      mu <- (pars[4] - pars[3])/2
    }
    if(3 %in% focal_pars & 5 %in% focal_pars){
      # div and ef
      mu <- (pars[3]*pars[5])/(1-pars[5])
    }
    if(4 %in% focal_pars & 5 %in% focal_pars){
      # turn and ef
      mu <- (pars[4]*pars[5])/(1+pars[5])
    }
  }else{
    mu <- pars[2]
  }
  return(mu)
  
}

convertBetweenPars <- function(pars){
  # pars <- c("lambda", "mu", "net.div", "turn", "ef")
  if(length(which(!is.na(pars))) >= 3){
    warning("More than 2 paramaters are specified. Randomly choosing 2 for the calculations.")
  }
  if(is.na(pars[1])){
    lambda <- convert2Lambda(pars)
  }else{
    lambda <- pars[1]
  }
  if(is.na(pars[2])){
    mu <- convert2Mu(pars)
  }else{
    mu <- pars[2]
  }
  net.div <- lambda - mu
  turn <- lambda + mu
  ef <- mu/lambda
  out <- c(lambda=lambda, mu=mu, net.div=net.div, turn=turn, ef=ef)
  if(!setequal(round(out[which(!is.na(pars))], 5), round(pars[which(!is.na(pars))], 5))){
    stop("An error occured because the calculated output doesn't match the input. Please check that your input parameters can be combined in a way that is possible.")
  }
  return(out)
}

# pars <- c(0.05, 0.05, NA, NA, NA)
# convertBetweenPars(pars)

# phy <- phy_a
# node <- mrca 

getPathToNode <- function(phy, tip, node){
  nTip <- length(phy$tip.label)
  if(node == "root"){
    node <- nTip + 1
  }
  if(is.character(tip)){
    tip <- which(phy$tip.label == tip)
  }
  path <- 0
  count <- 1
  while (tip != node) {
    tip.ind <- which(phy$edge[, 2] == tip)
    path <- c(path, tip.ind)
    count <- count + 1
    tip <- phy$edge[tip.ind, 1]
  }
  path <- path[-1]
  return(path)
}

# slice the extra branches
trimOvershotEdges <- function(sliced_phy, time_slice){
  node_ages <- branching.times.age(sliced_phy)
  node_tip_ages <- sliced_phy$edge
  node_tip_ages[,1] <- node_ages[match(sliced_phy$edge[,1], names(branching.times.age(sliced_phy)))]
  node_tip_ages[,2] <- node_tip_ages[,1] - sliced_phy$edge.length
  overshot_index <- node_tip_ages[,2] <= time_slice
  sliced_phy$edge.length[overshot_index] <- node_tip_ages[overshot_index,1] - time_slice
  return(sliced_phy)
}


SingleSlice <- function(phy, time_slice, print = TRUE, return.trim = FALSE){
  sliced_phy <- phy
  node_to_collapse <- getNodeToCollapse(phy, time_slice)
  while(!is.na(node_to_collapse)){
    sliced_phy <- dropClade(sliced_phy, node_to_collapse)
    node_to_collapse <- getNodeToCollapse(sliced_phy, time_slice)
  }
  if(return.trim){
    sliced_phy <- trimOvershotEdges(sliced_phy, time_slice)
  }
  if(print){
    par(mfrow=c(1,2))
    plot(phy, show.tip.label = FALSE)
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    abline(v = max(branching.times.age(phy)) - time_slice)
    plot(sliced_phy, show.tip.label = FALSE, x.lim = lastPP$x.lim)
    abline(v = max(branching.times.age(phy)) - time_slice)
  }
  return(sliced_phy)
}

dropClade <- function(phy, node_to_collapse){
  desc_to_drop <- phy$tip.label[get.descendants(node_to_collapse, phy, tips.only = TRUE)]
  if(length(desc_to_drop) > 1){
    new_phy <- drop.tip(phy, desc_to_drop[-1], trim.internal = TRUE)
    return(new_phy)
  }else{
    return(phy)
  }
}

getNodeToCollapse <- function(phy, time_slice){
  node_to_collapse <- as.numeric(names(branching.times.age(phy))[which(branching.times.age(phy) < time_slice)])[1]
  return(node_to_collapse)
}

several_sample_sizes <- function(n, lambda, mu, focal_size, return.phy=FALSE, nreplicates=1, ncores=1){
  out <- mclapply(seq(nreplicates), function(x) single_sample_size(n, lambda, mu, focal_size, return.phy=return.phy), mc.cores = ncores)
  return(out)
}


# convert pars alternative
quickConvert <- function(par, par.class){
    base_p <- c(lambda=NA, mu=NA, net.div=NA, turn=NA, ef=NA)
    base_p[match(par.class, names(base_p))] <- par
    p <- convertBetweenPars(base_p)
    names(p) <- c("lambda", "mu", "net.div", "turn", "ef")
    return(p)
}

# fool proof way to get extant tip
get_ntip <- function(phy){
  node_ages <- branching.times.age(phy)
  node_tip_ages <- phy$edge
  node_tip_ages[,1] <- node_ages[match(phy$edge[,1], names(branching.times.age(phy)))]
  node_tip_ages[,2] <- node_tip_ages[,1] - phy$edge.length
  n_tip <- length(which(round(node_tip_ages[,2], 5) == 0))
  return(n_tip)
}

# fool proof way to get branchng times (fool proof here is non-ultrametric treess)
branching.times.age <- function(phy){
  phy <- ladderize(phy)
  n <- length(phy$tip.label)
  e1 <- phy$edge[, 1]
  e2 <- phy$edge[, 2]
  EL <- phy$edge.length
  N <- length(e1)
  xx <- numeric(phy$Nnode)
  interns <- which(e2 > n)
  for (i in interns){
    xx[e2[i] - n] <- xx[e1[i] - n] + EL[i]
  } 
  # kinda concerned this isn't the correct depth under other situations...
  all_tip_lengths <- lapply(phy$tip.label, function(x) sum(phy$edge.length[getPathToNode(phy, tip = x, "root")]))
  depth <- max(unlist(all_tip_lengths))
  xx <- depth - xx
  names(xx) <- if (is.null(phy$node.label)){
    (n + 1):(n + phy$Nnode)
  }else{
    phy$node.label
  } 
  return(xx)
}


# gets some summary statistics from a time phylogeny with extinction
get_time_tree_summ_stats <- function(phy, return.rate=TRUE){
  obs_speciations <- phy$Nnode - 1
  n_tip <- get_ntip(phy)
  obs_extinctions <- length(phy$tip.label) - n_tip
  obs_frac <- obs_extinctions/obs_speciations
  obs_turn <- obs_speciations + obs_extinctions
  total_time <- sum(phy$edge.length)
  if(return.rate){
    lam_mu <- c(lambda = obs_speciations/total_time, 
    mu = obs_extinctions/total_time)
    pars <- quickConvert(lam_mu, c("lambda", "mu"))
    out <- c(pars, 
             mrca = max(branching.times.age(phy)),
             n = n_tip)
             #net.div = (obs_speciations/total_time) - (obs_extinctions/total_time),
             #turn = obs_turn/total_time,
             #ef = obs_frac) 
  }else{
    out <- c(speciations = obs_speciations, 
             extinctions = obs_extinctions, 
             rel_extinct = obs_frac, 
             total_event = obs_turn)
  }
  return(out)
}

compute_total_time <- function(phy) {
  return(sum(phy$edge.length)) 
}

# get the expected number of species following a birth-death process with lambda and mu as inputs
get_n_t <- function(n0, lambda, mu, t){
  numer <- n0 * exp((lambda - mu) * t)
  alpha_t <- get_alpha_t(lambda, mu, t)
  denom <- 1 - (alpha_t^n0)
  n_t <- numer/denom
  return(n_t)
}

# calculate alpha_t
get_alpha_t <- function(lambda, mu, t){
  fraction <- mu/lambda
  beta_t <- (exp((lambda - mu) * t) - 1)/(exp((lambda - mu) * t) - fraction)
  alpha_t <- fraction * beta_t
  return(alpha_t)
}

# test particular aspects of a tree
test_tree <- function(tree) {
    # test that it is a phylogeny
    test_1 <- is.phylo(tree)
    # test that it has more than 3 tips
    test_2 <- length(tree$tip.label) > 3
    # has it passed all tests?
    answer <- test_1 & test_2
    return(answer)
}

# create_phy_set <- function(n, lambda, mu, time_props = c(.25, .50, .75)){
#   phy_set <- list()
#   phy_set[[1]] <- ladderize(sim.bd.taxa(n,numbsim,lambda,mu,complete = TRUE)[[1]])
#   phy_set[[3]] <- ladderize(drop.extinct(phy_set[[1]]))
#   mrca <- getMRCA(phy_set[[1]], phy_set[[3]]$tip.label)
#   phy_set[[2]] <- ladderize(keep.tip(phy_set[[1]], phy_set[[1]]$tip.label[get.descendants(mrca, phy_set[[1]], tips.only = TRUE)]))
#   for(i in 1:length(time_props)){
#     node_ages <- branching.times.age(phy_set[[1]])
#     time_slice <- max(node_ages) * time_props[i]
#     phy_set[[3+i]] <- SingleSlice(phy = phy_set[[1]], time_slice = time_slice, FALSE, TRUE)
#   }
#   names(phy_set) <- c("complete", "crown", "extant", paste0("complete_", time_props))
#   return(phy_set)
# }

create_phy_set <- function(phy){
  phy_set <- list()
  phy_set[[1]] <- ladderize(phy)
  phy_set[[3]] <- ladderize(drop.extinct(phy_set[[1]]))
  mrca <- getMRCA(phy_set[[1]], phy_set[[3]]$tip.label)
  phy_set[[2]] <- ladderize(keep.tip(phy_set[[1]], phy_set[[1]]$tip.label[get.descendants(mrca, phy_set[[1]], tips.only = TRUE)]))
  names(phy_set) <- c("complete", "crown", "extant")
  return(phy_set)
}


plot_phy_set <- function(phy_list){
  phy_a <- phy_list[[1]]
  phy_b <- phy_list[[2]]
  phy_c <- phy_list[[3]]
  mrca <- getMRCA(phy_a, phy_c$tip.label)
  par(mfrow=c(1,3))
  edge_color_a <- rep("black", dim(phy_a$edge)[1])
  edge_color_b <- rep("black", dim(phy_b$edge)[1])
  edge_color_c <- rep("black", dim(phy_c$edge)[1])
  edge_index_a <- unique(unlist(lapply(phy_c$tip.label, function(x) getPathToNode(phy_a, x, mrca))))
  edge_index_b <- unique(unlist(lapply(phy_c$tip.label, function(x) getPathToNode(phy_b, x, "root"))))
  edge_color_a[edge_index_a] <- "red"
  edge_color_b[edge_index_b] <- "red"
  edge_color_c <- "red"
  plot(phy_a, show.tip.label = FALSE, edge.color = edge_color_a); axisPhylo()
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  plot(phy_b, show.tip.label = FALSE, edge.color = edge_color_b); axisPhylo()
  plot(phy_c, show.tip.label = FALSE, edge.color = edge_color_c); axisPhylo()
}


# a function for organizing data from a web of science search
getMetaGearDat <- function(d, type="WoS"){
  if(type=="WoS"){
    out <- data.frame(AUTHORS = d$Authors,
                      YEAR = d$Publication.Year,
                      TITLE = d$Article.Title,
                      JOURNAL = d$Journal.Abbreviation,
                      VOLUME = d$Volume,
                      LPAGES = d$Start.Page,
                      UPAGES = d$End.Page,
                      DOI = d$DOI,
                      ABSTRACT = d$Abstract)
  }
  return(out)
}

# a function for finding duplicated DOIs
findDuplicateStudies <- function(doi_list){
  expand.grid(1:length(doi_list), 1:length(doi_list))
  
  duplicated_studies <- Reduce(intersect, doi_list)
  return(duplicated_studies)
}


GetLikelihoodYule <- function(par, tree, rho=1, condition.type="survival", verbose=FALSE){

    n <- Ntip(tree)
    l <- par
    m <- 0
    x <- getx(tree)
    x <- sort(x,decreasing=TRUE)
    t <- x

    lik <- 0
    for (i in 2:length(t)){
        lik <- lik+log(l*p1(t[i],l,m,rho))
    }
    if(condition.type == "survival"){
        lik.unconditioned <- lik + 2 * log(p1(t[1],l,m,rho))
        condition.p <- 2 * log(1-p0(t[1],l,m,rho))
    }
    if(condition.type == "exactlyn"){
        lik.unconditioned <- lik + 2*log(p1(t[1],l,m,rho))
        part1 <- p1(t[1],l,m,rho)
        part2 <- qhelp(t[1],l,m,rho)
        condition.p <- log((n-1)*(part1^2)*(part2)^(n-2))
    }
    if(condition.type == "beauOmeara"){
        lik.unconditioned <- lik + 2*log(p1(t[1],l,m,rho))
        condition.p <- log(AscertainBias2(n,l,m,t[1],rho=1))
    }
    if(condition.type == "magsan"){
        lik.unconditioned <- lik + 2*log(p1(t[1],l,m,rho))
        condition.p <- log(AscertainBias1(n,l,m,t[1],rho=1))
    }
    if(condition.type == "exactlyn.stem0"){
        lik.unconditioned <- lik + log(l*p1(t[1],l,m,rho))
        condition.p <- log(n * (p1(t[1],l,m,rho) / (1-p0(t[1],l,m,rho))))
    }
    if(condition.type == "exactlyn.stem1"){
        lik.unconditioned <- lik + log(l*p1(t[i],l,m,rho))
        condition.p <- log((qhelp(t=t[1], l=l, m=m, rho=rho))^(1-n))
    }
    if(condition.type == "none") {
      lik.unconditioned <- lik + log(l*p1(t[i],l,m,rho))
      condition.p <- log(1)
    }
    if(verbose==TRUE){
        print(paste("unconditional prob", lik.unconditioned))
        print(paste("conditional prob", condition.p))
    }

    if(condition.type == "exactlyn.stem0" | condition.type == "exactlyn.stem1" | condition.type == "exactlyn" | condition.type == "magsan" | condition.type == "morethann"){
        #tmp <- c(l+m, m/l, lik.unconditioned, log(condition.p))
        lik.conditioned.both.check <- lik.unconditioned - condition.p
        names(lik.conditioned.both.check) <- "llik"
        return(-lik.conditioned.both.check)
    }else{
        lik.conditioned.both.check <- lik.unconditioned - condition.p
        names(lik.conditioned.both.check) <- "llik"
        return(-lik.conditioned.both.check)
    }
}


GetLikelihoodEiffel <- function(par, tree, rho=1, condition.type="survival", verbose=FALSE){

    n <- Ntip(tree)
    turn <- par
    ef <- 0.99
    
    l <- turn / (1 + ef)
    m <- (turn * ef) / (1 + ef)
    
    x <- getx(tree)
    x <- sort(x,decreasing=TRUE)
    t <- x

    lik <- 0
    for (i in 2:length(t)){
        lik <- lik+log(l*p1(t[i],l,m,rho))
    }
    if(condition.type == "survival"){
        lik.unconditioned <- lik + 2 * log(p1(t[1],l,m,rho))
        condition.p <- 2 * log(1-p0(t[1],l,m,rho))
    }
    if(condition.type == "exactlyn"){
        lik.unconditioned <- lik + 2*log(p1(t[1],l,m,rho))
        part1 <- p1(t[1],l,m,rho)
        part2 <- qhelp(t[1],l,m,rho)
        condition.p <- log((n-1)*(part1^2)*(part2)^(n-2))
    }
    if(condition.type == "beauOmeara"){
        lik.unconditioned <- lik + 2*log(p1(t[1],l,m,rho))
        condition.p <- log(AscertainBias2(n,l,m,t[1],rho=1))
    }
    if(condition.type == "magsan"){
        lik.unconditioned <- lik + 2*log(p1(t[1],l,m,rho))
        condition.p <- log(AscertainBias1(n,l,m,t[1],rho=1))
    }
    if(condition.type == "exactlyn.stem0"){
        lik.unconditioned <- lik + log(l*p1(t[1],l,m,rho))
        condition.p <- log(n * (p1(t[1],l,m,rho) / (1-p0(t[1],l,m,rho))))
    }
    if(condition.type == "exactlyn.stem1"){
        lik.unconditioned <- lik + log(l*p1(t[i],l,m,rho))
        condition.p <- log((qhelp(t=t[1], l=l, m=m, rho=rho))^(1-n))
    }
    if(condition.type == "none") {
      lik.unconditioned <- lik + log(l*p1(t[i],l,m,rho))
      condition.p <- log(1)
    }
    if(verbose==TRUE){
        print(paste("unconditional prob", lik.unconditioned))
        print(paste("conditional prob", condition.p))
    }

    if(condition.type == "exactlyn.stem0" | condition.type == "exactlyn.stem1" | condition.type == "exactlyn" | condition.type == "magsan" | condition.type == "morethann"){
        #tmp <- c(l+m, m/l, lik.unconditioned, log(condition.p))
        lik.conditioned.both.check <- lik.unconditioned - condition.p
        names(lik.conditioned.both.check) <- "llik"
        return(-lik.conditioned.both.check)
    }else{
        lik.conditioned.both.check <- lik.unconditioned - condition.p
        names(lik.conditioned.both.check) <- "llik"
        return(-lik.conditioned.both.check)
    }
}




GetLikelihoodBD <- function(par, tree, rho=1, condition.type="survival", verbose=FALSE, log.convert=TRUE){

    if(log.convert){
      p.new <- exp(par)
    }else{
        p.new <- par
    }
    
    n <- Ntip(tree)
    l <- p.new[1] / (1 - p.new[2])
    m <- (p.new[1]*p.new[2]) / (1 - p.new[2])
    x <- getx(tree)
    x <- sort(x,decreasing=TRUE)
    t <- x

    lik <- 0
    for (i in 2:length(t)){
        lik <- lik+log(l*p1(t[i],l,m,rho))
    }
    if(condition.type == "survival"){
        lik.unconditioned <- lik + 2 * log(p1(t[1],l,m,rho))
        condition.p <- 2 * log(1-p0(t[1],l,m,rho))
    }
    if(condition.type == "exactlyn"){
        lik.unconditioned <- lik + 2*log(p1(t[1],l,m,rho))
        part1 <- p1(t[1],l,m,rho)
        part2 <- qhelp(t[1],l,m,rho)
        condition.p <- log((n-1)*(part1^2)*(part2)^(n-2))
    }
    if(condition.type == "beauOmeara"){
        lik.unconditioned <- lik + 2*log(p1(t[1],l,m,rho))
        condition.p <- log(AscertainBias2(n,l,m,t[1],rho=1))
    }
    if(condition.type == "magsan"){
        lik.unconditioned <- lik + 2*log(p1(t[1],l,m,rho))
        condition.p <- log(AscertainBias1(n,l,m,t[1],rho=1))
    }
    if(condition.type == "exactlyn.stem0"){
        lik.unconditioned <- lik + log(l*p1(t[1],l,m,rho))
        condition.p <- log(n * (p1(t[1],l,m,rho) / (1-p0(t[1],l,m,rho))))
    }
    if(condition.type == "exactlyn.stem1"){
        lik.unconditioned <- lik + log(l*p1(t[i],l,m,rho))
        condition.p <- log((qhelp(t=t[1], l=l, m=m, rho=rho))^(1-n))
    }
    if(condition.type == "none") {
      lik.unconditioned <- lik + log(l*p1(t[i],l,m,rho))
      condition.p <- log(1)
    }
    if(verbose==TRUE){
        print(paste("unconditional prob", lik.unconditioned))
        print(paste("conditional prob", condition.p))
    }

    if(condition.type == "exactlyn.stem0" | condition.type == "exactlyn.stem1" | condition.type == "exactlyn" | condition.type == "magsan" | condition.type == "morethann"){
        #tmp <- c(l+m, m/l, lik.unconditioned, log(condition.p))
        lik.conditioned.both.check <- lik.unconditioned - condition.p
        names(lik.conditioned.both.check) <- "llik"
        return(-lik.conditioned.both.check)
    }else{
        lik.conditioned.both.check <- lik.unconditioned - condition.p
        names(lik.conditioned.both.check) <- "llik"
        return(-lik.conditioned.both.check)
    }
}


#opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = 100000, "ftol_rel" = .Machine$double.eps^.50)
#obj <- nloptr(x0=log(c(.1,.5)), eval_f=GetLikelihoodBD, ub=log(c(100,.99)), lb=c(-6,-6), opts=opts, tree=trees[[1]], rho=1, condition.type="survival", verbose=FALSE, log.convert=TRUE)
