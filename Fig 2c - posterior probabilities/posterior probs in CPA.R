library(ggplot2)
library(ggthemes)
library(binom)
library(MLEcell)
library(clue)
library(pheatmap)
rm(list = ls())


#### load CPA data: ---------------
load("../data/Cell Pellet Array annotation and raw counts.RData")
load("../data/CPA16_RNAseq.RData")
counts = t(as.matrix(raw))
rm(raw)
badprobes = read.csv("../data/genes with efficiency 8-fold below average in old CPA panel.csv")[,2]
counts = counts[, !is.element(colnames(counts), badprobes)]
boxplot(log(annot$raw_totalCount)~annot$cell_line,las=2)

abline(h = 150, col = 2)
# remove the cell lines with failed FOVs:
failed.cell.lines = names(which(by(annot$raw_totalCount, annot$cell_line, median) < 150))
remove = is.element(annot$cell_line, failed.cell.lines)
annot = annot[!remove, ]
counts = counts[!remove, ]


#### supervised cell typing -----------------------

# get reference:
meanprofiles = MLEcell:::Estep(counts = counts, clust = annot$cell_line, neg = annot$negmean)

# without cohorts 
sup = insitutypeML(counts = counts,
                   neg = annot$negmean,
                   reference_profiles = meanprofiles, #CPA16_RNAseq,
                   cohort = NULL)

#### compare err rate vs. posterior prob:

# get relevant stats:
annot$celltype = sup$clust
annot$prob = sup$prob
annot$prob20 = cut(sup$prob, breaks = c(0,0.2,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.975,0.99, 0.999, 1))
annot$err = 1 * (annot$celltype != annot$cell_line)

# summary stats:
bins = levels(annot$prob20)
#bins.l = as.numeric(substr(bins,2,unlist(gregexpr(",", bins))-1))
#bins.u = as.numeric(substr(bins,unlist(gregexpr(",", bins))+1, nchar(bins)-1))

means = lci = uci = nominal = binsize = c()
for (bin in levels(annot$prob20)) {
    use = (annot$prob20 == bin) 
    tmp = t(binom.confint(x = sum(1 - annot$err[use]), n = sum(use), method = "wilson")[, c("mean", "lower", "upper")])
    means[bin] = tmp[1]
    lci[bin] = tmp[2]
    uci[bin] = tmp[3]
    nominal[bin] = mean(annot$prob[use])
    binsize[bin] = sum(use)
}

textlabs = rep("", length(binsize))
textlabs[1:4] = " cells"

svg("results/posterior probs vs. error rates.svg", width = 10, height = 6)
layout(matrix(1:2, 2), heights = c(10, 1))
par(mar = c(4,4.5,.1,.1))
plot(nominal, means, log = "xy",
     pch = 16, cex = 1.5, col = "darkblue",
     xlim = c(0.5,1.01), ylab = "Prop. correct", xlab = "", ylim = c(0.4,1.1),
     cex.lab = 1.5, xaxt = "n")
axis(1, at = nominal, las = 2, lab = c(names(nominal)[1:10], "", ""), cex = 0.5)
axis(1, at = c(0.9999973, 1.02), las = 2, lab = names(nominal)[11:12], cex = 0.5, tick = F)
for (i in 1:length(bins)) {
  lines(rep(nominal[i], 2), c(lci[i], uci[i]),
        col = alpha("darkblue", 0.75))
  text(nominal[i], uci[i]*1.02, paste0(binsize[i], textlabs[i]))
}
abline(0,1, col = "grey50")
par(mar = c(0,0,0,0))
frame()
legend("top", legend = "                Posterior probability", bty = "n", cex = 1.5)
dev.off()



# stats in manuscript:
mean(sup$prob == 1)
binom.confint(sum(annot$err[sup$prob == 1]), sum(sup$prob==1), method = "wilson")
