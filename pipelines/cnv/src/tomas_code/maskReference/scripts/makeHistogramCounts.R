## Import & edit
counts <- read.delim("counts.txt", h = F)
names(counts) <- c("chrom", "start", "end", "seq", "counts")


## Save & load
#save.image("makeHistogramCounts.RData")
#load("makeHistogramCounts.RData")


## How many kmers?
length(counts[,1])


## Cumulative distribution
counts.table <- data.frame(table(counts$counts))
names(counts.table) <- c("counts", "num.kmers")
counts.table <- transform(counts.table, freq = counts.table$num.kmers / sum(counts.table$num.kmers))
counts.table$cum.freq <- rep(0, length(counts.table[, 1]))
for (i in 1:length(counts.table[, 1])) {
	counts.table[i, ]$cum.freq <- sum(counts.table[1:i, ]$freq)
}


## Plot
#x11()
png("cumulativeDistribution.png", h = 700, w = 700, res = 100)
plot(counts.table$counts, counts.table$cum.freq, xlim = c(0, 100),
	ylab = "Cumulative percentage",
	xlab = "# placements for a given K-mer",
	main = "FelCat5 (K-mer size = 36, step = 5)")
abline(v = 19, col = "blue", lwd = 2)
legend("topright", "20 placements (89% cumulative percentage)", lty = 1, lwd = 2, col = "blue", bty = "n")
dev.off()


## Select those kmers with at least 20 placements
write.table(counts[counts$counts >= 20, ], "kmerCountsWithMrsFast_more20placements.txt", sep = "\t", quote = F, row.names = F, col.names = F)












## Histogram
x11()
hist(counts$counts, breaks = 1e5, xlim = c(0, 100))


png("figures/histograms.png", h = 960, w = 960, res = 100)
#x11(w = 12, h = 10)
par(mfrow = c(2, 1))
plot(counts.table$counts, counts.table$num.kmers,
	ylab = "# K-mers",
	main = "All K-mers",
	xlab = "# placement for a given K-mer")
abline(h = mean(counts.table$num.kmers), col = "blue", lwd = 2)
legend("topright", "Threshold: mean", lty = 1, lwd = 2, col = "blue", bty = "n")
plot(counts.table$counts, counts.table$num.kmers,
	ylab = "# K-mers",
	main = "K-mers with <50 counts",
	xlab = "# placements for a given K-mer",
	xlim = c(0, 50))
abline(h = mean(counts.table$num.kmers), col = "blue", lwd = 2)
dev.off()



