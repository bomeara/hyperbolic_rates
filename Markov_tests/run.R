

######################################################################################################################################
######################################################################################################################################
### Loading Packages
######################################################################################################################################
######################################################################################################################################

library(corHMM)
library(TreeSim)
library(TreeTools)
library(ggplot2)


######################################################################################################################################
######################################################################################################################################
### Basic Sim Function
######################################################################################################################################
######################################################################################################################################

DoRun <- function(nreps=1000, phy, markov.rate=0.1, error=0.1){
	res <- c()
	ages <- exp(runif(nreps, log(5), log(20)))
	Q <- rbind(c(-markov.rate, markov.rate), c(markov.rate, -markov.rate))
	for(age.index in 1:length(ages)){
		phy$edge.length <- phy$edge.length*ages[age.index]/max(branching.times(phy))
		traits <- corHMM:::simMarkov(phy, Q=Q, root.freqs=c(1,0))
		phy$node.label <- NULL
		phydat <- data.frame(taxon=phy$tip.label, Reg=traits$TipStates)
		if(length(table(phydat[,2]))>1){
			pp.noerror <- corHMM(phy, phydat, model="ER", rate.cat=1)
		}else{
			pp.noerror <- NULL
			pp.noerror$solution <- as.numeric(4)
		}
		
		error.absolute <- round(dim(phydat)[1]*error)
		taxa.sample <- sample(1:dim(phydat)[1], error.absolute)
		phydat.new <- phydat
		for(index in 1:length(taxa.sample)){
			state <- phydat[taxa.sample[index],2]
			if(state == 1){
				state <- 2
			}else{
				state <- 1
			}
			phydat.new[taxa.sample[index],2] <- state
		}
		pp.error <- corHMM(phy, phydat.new, model="ER", rate.cat=1)
#		if(pp.noerror$solution[2]>90){
#			save(phy, phydat.new, phydat, file=paste("Bad.age", age.index, ".Rsave", sep=""))
#		}
#		if(pp.error$solution[2]>90){
#			save(phy, phydat.new, phydat, file=paste("Bad.age", age.index, ".Rsave", sep=""))
#		}
		res <- rbind(res, c(ages[age.index], pp.noerror$solution[2], pp.error$solution[2], sum(phy$edge.length)))
	}
	return(res)
}


######################################################################################################################################
######################################################################################################################################
### Perform Sim on Pectinate tree and a Balanced tree
######################################################################################################################################
######################################################################################################################################

phy.pect <- PectinateTree(64)
phy.pect <- compute.brlen(phy.pect)
phy.pect$edge.length <- phy.pect$edge.length*5/max(branching.times(phy.pect))
sum(phy.pect$edge.length) * .01
phy.pect$edge.length <- phy.pect$edge.length*20/max(branching.times(phy.pect))
sum(phy.pect$edge.length) * .01
pp.pect <- DoRun(nreps=1000, phy=phy.pect, markov.rate=0.01, error=0.1)

phy.balan <- BalancedTree(64)
phy.balan <- compute.brlen(phy.balan)
phy.balan$edge.length <- phy.balan$edge.length*5/max(branching.times(phy.balan))
sum(phy.balan$edge.length) * .05
phy.balan$edge.length <- phy.balan$edge.length*20/max(branching.times(phy.balan))
sum(phy.balan$edge.length) * .05
pp.balan<- DoRun(nreps=1000, phy=phy.balan, markov.rate=0.05, error=0.1)

save(pp.pect, pp.balan, file="corHMMResults_samerate.Rsave")

######################################################################################################################################
######################################################################################################################################
### Plotting Results using Clade Age
######################################################################################################################################
######################################################################################################################################

load("corHMMResults.Rsave")

dat <- data.frame(age=pp.pect[,1], r_noerr=pp.pect[,2], r_err=pp.pect[,3])
#dat[is.na(dat[,2]),2] <- 0
dat <- dat[-which(dat[,3]>20),]

p_simulated.pect <- ggplot(dat, aes(x = age, y = r_noerr)) +
	geom_point(size=.75, alpha=.1) +
	theme_bw() +
	theme(panel.grid = element_blank(), legend.position = "none") +
	theme(panel.border = element_blank(), axis.line.x = element_line(size = 0.25, linetype = "solid", colour = "black"), axis.line.y = element_line(size = 0.25, linetype = "solid", colour = "black"), axis.ticks=element_line(size=0.25)) +
	#geom_errorbar(aes(ymin = lower.CI, ymax = upper.CI), size=0.2) +
	coord_cartesian(ylim = c(0.001, .2)) +
	scale_y_log10() +
	scale_x_log10() +
	geom_smooth() +
	geom_hline(yintercept=0.01, lty=3, col="red") +
	labs(y=expression(transitions~Myr^-1), x="", tag="A")


p_simulatederror.pect <- ggplot(dat, aes(x = age, y = r_err)) +
	geom_point(size=.75, alpha=.1) +
	theme_bw() +
	theme(panel.grid = element_blank(), legend.position = "none") +
	theme(panel.border = element_blank(), axis.line.x = element_line(size = 0.25, linetype = "solid", colour = "black"), axis.line.y = element_line(size = 0.25, linetype = "solid", colour = "black"), axis.ticks=element_line(size=0.25)) +
	#geom_errorbar(aes(ymin = lower.CI, ymax = upper.CI), size=0.2) +
coord_cartesian(ylim = c(0.001, .2)) +
	scale_y_log10() +
	scale_x_log10() +
	geom_smooth() +
	geom_hline(yintercept=0.01, lty=3, col="red") +
	labs(y="", x="",tag="B")


dat <- data.frame(age=pp.balan[,1], r_noerr=pp.balan[,2], r_err=pp.balan[,3])
#dat[is.na(dat[,2]),2] <- 0
#dat <- dat[-which(dat[,3]>90),]

p_simulated.balan <- ggplot(dat, aes(x = age, y = r_noerr)) +
	geom_point(size=.75, alpha=.1) +
	theme_bw() +
	theme(panel.grid = element_blank(), legend.position = "none") +
	theme(panel.border = element_blank(), axis.line.x = element_line(size = 0.25, linetype = "solid", colour = "black"), axis.line.y = element_line(size = 0.25, linetype = "solid", colour = "black"), axis.ticks=element_line(size=0.25)) +
	#geom_errorbar(aes(ymin = lower.CI, ymax = upper.CI), size=0.2) +
	#coord_cartesian(ylim = c(0, 1)) +
	coord_cartesian(ylim = c(0.01, 1)) +
	scale_y_log10() +
	scale_x_log10() +
	geom_smooth() +
	geom_hline(yintercept=0.05, lty=3, col="red") +
	labs(y=expression(transitions~Myr^-1), x="Clade Age (Myr)", tag="C")


p_simulatederror.balan <- ggplot(dat, aes(x = age, y = r_err)) +
	geom_point(size=.75, alpha=.1) +
	theme_bw() +
	theme(panel.grid = element_blank(), legend.position = "none") +
	theme(panel.border = element_blank(), axis.line.x = element_line(size = 0.25, linetype = "solid", colour = "black"), axis.line.y = element_line(size = 0.25, linetype = "solid", colour = "black"), axis.ticks=element_line(size=0.25)) +
	#geom_errorbar(aes(ymin = lower.CI, ymax = upper.CI), size=0.2) +
	#coord_cartesian(ylim = c(0, 1)) +
	coord_cartesian(ylim = c(0.01, 1)) +
	scale_y_log10() +
	scale_x_log10() +
	geom_smooth() +
	geom_hline(yintercept=0.05, lty=3, col="red") +
	labs(y="", x="Clade Age (Myr)", tag="D")

out <- gridExtra::grid.arrange(p_simulated.pect, p_simulatederror.pect, p_simulated.balan, p_simulatederror.balan, NULL, NULL, nrow = 3, ncol=2)
ggsave(file = "FigS9.pdf", plot = out, width = 7.24, height = 7.24)


######################################################################################################################################
######################################################################################################################################
### Plotting Results using Total Time
######################################################################################################################################
######################################################################################################################################

dat <- data.frame(age=pp.pect[,4], r_noerr=pp.pect[,2], r_err=pp.pect[,3])
#dat[is.na(dat[,2]),2] <- 0
dat <- dat[-which(dat[,3]>20),]

p_simulated.pect <- ggplot(dat, aes(x = age, y = r_noerr)) +
	geom_point(size=.75, alpha=.1) +
	theme_bw() +
	theme(panel.grid = element_blank(), legend.position = "none") +
	theme(panel.border = element_blank(), axis.line.x = element_line(size = 0.25, linetype = "solid", colour = "black"), axis.line.y = element_line(size = 0.25, linetype = "solid", colour = "black"), axis.ticks=element_line(size=0.25)) +
	#geom_errorbar(aes(ymin = lower.CI, ymax = upper.CI), size=0.2) +
	coord_cartesian(ylim = c(0.001, .3)) +
	#coord_cartesian(xlim = c(0, 700)) +
	scale_y_log10() +
	scale_x_log10() +
	geom_smooth() +
	geom_hline(yintercept=0.01, lty=3, col="red") +
	labs(y=expression(transitions~Myr^-1), x="", tag="A")


p_simulatederror.pect <- ggplot(dat, aes(x = age, y = r_err)) +
	geom_point(size=.75, alpha=.1) +
	theme_bw() +
	theme(panel.grid = element_blank(), legend.position = "none") +
	theme(panel.border = element_blank(), axis.line.x = element_line(size = 0.25, linetype = "solid", colour = "black"), axis.line.y = element_line(size = 0.25, linetype = "solid", colour = "black"), axis.ticks=element_line(size=0.25)) +
	#geom_errorbar(aes(ymin = lower.CI, ymax = upper.CI), size=0.2) +
	coord_cartesian(ylim = c(0.001, .3)) +
	#coord_cartesian(xlim = c(0, 700)) +
	scale_y_log10() +
	scale_x_log10() +
	geom_smooth() +
	geom_hline(yintercept=0.01, lty=3, col="red") +
	labs(y="", x="",tag="B")


dat <- data.frame(age=pp.balan[,4], r_noerr=pp.balan[,2], r_err=pp.balan[,3])
#dat[is.na(dat[,2]),2] <- 0
#dat <- dat[-which(dat[,3]>90),]

p_simulated.balan <- ggplot(dat, aes(x = age, y = r_noerr)) +
	geom_point(size=.75, alpha=.1) +
	theme_bw() +
	theme(panel.grid = element_blank(), legend.position = "none") +
	theme(panel.border = element_blank(), axis.line.x = element_line(size = 0.25, linetype = "solid", colour = "black"), axis.line.y = element_line(size = 0.25, linetype = "solid", colour = "black"), axis.ticks=element_line(size=0.25)) +
	#geom_errorbar(aes(ymin = lower.CI, ymax = upper.CI), size=0.2) +
	coord_cartesian(ylim = c(0.01, 1)) +
	#coord_cartesian(xlim = c(0, 700)) +
	scale_y_log10() +
	scale_x_log10() +
	geom_smooth() +
	geom_hline(yintercept=0.05, lty=3, col="red") +
	labs(y=expression(transitions~Myr^-1), x="Total time (Myr)", tag="C")


p_simulatederror.balan <- ggplot(dat, aes(x = age, y = r_err)) +
	geom_point(size=.75, alpha=.1) +
	theme_bw() +
	theme(panel.grid = element_blank(), legend.position = "none") +
	theme(panel.border = element_blank(), axis.line.x = element_line(size = 0.25, linetype = "solid", colour = "black"), axis.line.y = element_line(size = 0.25, linetype = "solid", colour = "black"), axis.ticks=element_line(size=0.25)) +
	#geom_errorbar(aes(ymin = lower.CI, ymax = upper.CI), size=0.2) +
	coord_cartesian(ylim = c(0.01, 1)) +
	#coord_cartesian(xlim = c(0, 700)) +
	scale_y_log10() +
	scale_x_log10() +
	geom_smooth() +
	geom_hline(yintercept=0.05, lty=3, col="red") +
	labs(y="", x="Total time (Myr)", tag="D")

out <- gridExtra::grid.arrange(p_simulated.pect, p_simulatederror.pect, p_simulated.balan, p_simulatederror.balan, NULL, NULL, nrow = 3, ncol=2)
ggsave(file = "FigS10.pdf", plot = out, width = 7.24, height = 7.24)
