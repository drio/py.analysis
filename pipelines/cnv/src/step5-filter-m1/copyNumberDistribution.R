# USAGE: copy.number.distributions.R to:
# Calculate genome-wide copy number distribution
# Calculate copy number distribution in control regions
# Calculate copy number distribution in non-control regions
# Calculate copy number distribution in non-control regions but only on autosomes
# Gain/loss cut-offs: using only control regions, normalize copy number (using an Stdev excluding windows with the 1% most extreme copy number) and define cut-offs as +-3 Stdev of the normalized distribution



##############################################################################
## Data preparation
##############################################################################

samples <- c("XCATS")

c.cols <- colors()[grep("orange", colors())][1:length(c)]
my.cols <- c.cols

## Import read depth and copy number per window
cn <- list()
cn.dist <- list()
cn.control.dist <- list()
cn.non.control.dist <- list()
cn.non.control.autosomes.dist <- list()

for (i in 1:length(samples)) {
	print(samples[i])
	## Import
  path.import1 <- paste("output.cw_norm.bed", sep = "")
  path.import2 <- paste("output.copynumber.bed", sep = "")
	x <- read.delim(path.import1)
	y <- read.delim(path.import2)
	names(x)[c(1, 4)] <- c("CHROM", "GC")
	names(y)[c(1, 4)] <- c("CHROM", "GC")
	# Merge
	z <- merge(x, y, merge = name(x)[1:4], sort = F)
	## Calculate copy number distribution in non-control regions but only on autosomes
	cn.non.control.autosomes.dist[[i]] <- data.frame(table(cut(z[z$CHROM != "chrUn" & z$CHROM != "chrX", ]$COPYNUMBER,
		breaks = c(seq(0, 100, 0.1), 1e3, 1e4, 1e5), right = T, include.lowest = T)))
	row.names(cn.non.control.autosomes.dist[[i]]) <- c(seq(0, 100, 0.1), 1e3, 1e4)

	cn[[i]] <- z
}


# Save data
save.image("copy.number.distribution.RData")


# estimate the maximum copy number for control/non-control in each dataset
k1 <- c()
k2 <- c()
for (i in 1:length(cn)) {
	k1 <- c(k1, max(cn[[i]][cn[[i]]$IS_CONTROL ==  "Y" , ]$COPYNUMBER))
	k2 <- c(k2, max(cn[[i]][cn[[i]]$IS_CONTROL ==  "N" , ]$COPYNUMBER))
}
max(k1)
max(k2)

##############################################################################
# Gain/loss cut-offs
##############################################################################

# Using only control regions, normalize copy number (using an Stdev excluding windows with the 1% most extreme copy number) and define cut-offs as +-3 Stdev of the normalized distribution
cutoffs <- data.frame(matrix(ncol = 7, nrow = 0))
names(cutoffs) <- c("sample", "mean", "stDev_99", "wins.excl", "p.wins.excl", "gain_cutoff", "loss_cutoff")
for (i in 1:length(samples)) {
	print(samples[i])
	# Subset control regions
	x <- cn[[i]]
	control <- x[x$IS_CONTROL == "Y", ]
	# Normalization with parameters
	# mean
	mu <- mean(control$COPYNUMBER)
	# Standard deviation (excluding windows with the 1% most extreme copy number)
	cn.q99 <- quantile(control$COPYNUMBER, probs = c(0.99))[[1]]
	stDev <- sd(control[control$COPYNUMBER < cn.q99, ]$COPYNUMBER)
	wins.excl <- length(control[control$COPYNUMBER > cn.q99, ]$COPYNUMBER)
	p.wins.excl <- (wins.excl * 100) / length(control$COPYNUMBER)
	control$COPYNUMBER_NORM <- (control$COPYNUMBER - mu) / stDev
	#print(c(mu, stDev))
 	# Cut-offs
	gain.cutoff <- (3 * stDev) + mu
	loss.cutoff <- (-3 * stDev) + mu
	y <- data.frame(samples[i], mu, stDev, wins.excl, p.wins.excl, gain.cutoff, loss.cutoff)
	names(y) <- names(cutoffs)
	cutoffs <- rbind(cutoffs, y)
}

cutoffs$species <- c(1:1)

data4plot <- cutoffs
data4plot <- data4plot[1:1, ]

# Export
write.table(data4plot, "sample.specific.cutoffs.txt", sep = "\t", quote = F, row.names = F)
