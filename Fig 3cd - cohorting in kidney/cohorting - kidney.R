library(RColorBrewer)
library(MLEcell)
library(binom)
library(vioplot)
library(scales)
library(rlist)
rm(list = ls())

source("getcellneighbors.R")

#### load data: ---------------
load("../data/lupus nephritis SP11_1139 data.RData")
refprofiles = readRDS("../data/kidney HCA reference matrix.RDs")


#### cohort based on immunofluorescence data ----------------------

## cohort based on the if data:
set.seed(0)
ifdata = as.matrix(annot[, c("Area","AspectRatio","Mean.CD298","Mean.PanCK","Mean.CD45","Mean.CD20","Mean.DAPI")])
ifcohort = fastCohorting(mat = ifdata)

## cohort based on spatial context:
set.seed(0)

# get neighbor cell IDs:
nmat = getcellneighbors(x = xy[, 1], y = xy[, 2], cell = annot$cell_ID, 
                        tissue = annot$tissue, n.neighbors = 50)
saveRDS(nmat, file = "results/nmat.RDS")
# get mean neighborhood expression (use original SCH code):
neighborhoodmean = counts * NA
for (id in rownames(neighborhoodmean)) {
  neighborhoodmean[id, ] = colMeans(counts[nmat[id, ], ])
}
saveRDS(neighborhoodmean, file = "results/neighborhoodmean.RDS")

# reduce to 10 PCs
neighborhoodPCs = prcomp(neighborhoodmean)$x[,1:10]
saveRDS(neighborhoodPCs, file = "results/neighborhoodPCs.RDS")

set.seed(0)
spatcohort = fastCohorting(mat = neighborhoodPCs)

## cohort based on IF and spatial context:
set.seed(0)
bothcohort = fastCohorting(mat = cbind(ifdata, neighborhoodPCs))


## random/uninformative cohorts:
randomcohort = sample(ifcohort, length(ifcohort), replace = FALSE)


#### supervised cell typing -----------------------


# without cohorts 
sup.nocohorts = insitutypeML(counts = counts,
                             neg = annot$negmean,
                             reference_profiles = refprofiles,
                             cohort = NULL)

# with IF cohorts
sup.ifcohorts = insitutypeML(counts = counts,
                             neg = annot$negmean,
                             reference_profiles = refprofiles,
                             cohort = ifcohort)

# with spatial context cohorts
sup.spatcohorts = insitutypeML(counts = counts,
                             neg = annot$negmean,
                             reference_profiles = refprofiles,
                             cohort = spatcohort)

# with BOTH IF and spatial cohorts
sup.bothcohorts = insitutypeML(counts = counts,
                             neg = annot$negmean,
                             reference_profiles = refprofiles,
                             cohort = bothcohort)

#with uninformative cohorts
sup.randomcohorts = insitutypeML(counts = counts,
                             neg = annot$negmean,
                             reference_profiles = refprofiles,
                             cohort = randomcohort)

save(sup.nocohorts, sup.ifcohorts, sup.spatcohorts, sup.bothcohorts, sup.randomcohorts, file = "results/supervised results.RData")


#### summarize results -----------------------------------


# rates of differences per expression level:

subsets = list(all = rep(TRUE, nrow(counts)),
               verylow = (rowSums(counts) <= 50),
               low = (rowSums(counts) > 50) & (rowSums(counts) <= 100),
               med = (rowSums(counts) > 100) & (rowSums(counts) <= 200),
               high = (rowSums(counts) > 200))

res = c()
for (name in names(subsets)) {
  use = subsets[[name]]

  ci.if = t(t(binom.confint(sum(sup.nocohorts$clust[use] != sup.ifcohorts$clust[use]), length(sup.nocohorts$clust[use]), methods = "wilson")[, c("mean", "lower", "upper")]))
  ci.spat = t(t(binom.confint(sum(sup.nocohorts$clust[use] != sup.spatcohorts$clust[use]), length(sup.nocohorts$clust[use]), methods = "wilson")[, c("mean", "lower", "upper")]))
  ci.both = t(t(binom.confint(sum(sup.nocohorts$clust[use] != sup.bothcohorts$clust[use]), length(sup.nocohorts$clust[use]), methods = "wilson")[, c("mean", "lower", "upper")]))
  ci.rand = t(t(binom.confint(sum(sup.nocohorts$clust[use] != sup.randomcohorts$clust[use]), length(sup.nocohorts$clust[use]), methods = "wilson")[, c("mean", "lower", "upper")]))

  res = rbind(res, c(ci.if, ci.spat, ci.both, ci.rand))
}
colnames(res) = c("mean IF", "lci.if", "uci.if", 
                  "mean spat", "lci.spat", "uci.spat", 
                  "mean both", "lci.both", "uci.both", 
                  "mean rand", "lci.rand", "uci.rand")
rownames(res) = names(subsets)
print(res)
print(res[, grepl("mean", colnames(res))])
write.csv(res, "results/mean and 95 conf intervals for changes from no cohort.csv")
res = res*100
# same cells aren't changed, right?
changed.by.if = (sup.nocohorts$clust != sup.ifcohorts$clust)
changed.by.spat = (sup.nocohorts$clust != sup.spatcohorts$clust)
table(changed.by.spat, changed.by.if)
table(changed.by.spat)
table(changed.by.if)


# summary barplot of change rates:
svg("results/barplot of changes from no cohorts.svg", width = 7, height = 4.5)
par(mar = c(5,6,2,1))
bpcols = c("orange", "dodgerblue3", "darkorchid1", "grey50")
bp = barplot(t(res[, grepl("mean", colnames(res))]), beside = T,
        ylab = "Percent of cell type assignments changed\nfrom gene expression-only cell typing", cex.lab = 1.2,
        names = c("all cells\n ", "0-50\ntranscripts", "51-100\ntranscripts", "101-200\ntranscripts", ">200\ntranscripts"),
        col = bpcols, ylim = c(0, max(res)), cex.names = 1)
legend("topright", pch = 15, col = c(NA, bpcols), cex = 1.2,
       legend = c("Alternative data used:", "Image data", "Spatial context", "Image + spatial", "Random cohorts"))
for (i in 1:5) {
  for (j in 1:4) {
    lci = res[i, which(grepl("lci", colnames(res)))[j]]
    uci = res[i, which(grepl("uci", colnames(res)))[j]]
    xpos = bp[j,i]
    lines(rep(xpos,2), c(lci, uci))
  }
}
dev.off()

#### what characterizes the changed cells? -------------------

# histogram of posterior probs:
cols2 = alpha(c("grey30", "orange"), 0.5)
svg("results/hist of posterior probs.svg", width = 7, height = 4.5)
hist(sup.nocohorts$prob[sup.nocohorts$clust == sup.bothcohorts$clust], 
     xlim = c(0,1), freq = F, breaks = 50, col = cols2[1], main = "", cex.lab = 1.2,
     xlab = "Cell type posterior probability from gene expression only")
hist(sup.nocohorts$prob[sup.nocohorts$clust != sup.bothcohorts$clust], 
     col = cols2[2], freq = F, breaks = 50, add = T)
legend("topleft", fill = rev(cols2), legend = c("Cell type changed by alternative data types", "Unchanged"),
       cex = 1.2)
dev.off()

# summary stats:
mean(sup.nocohorts$prob[sup.nocohorts$clust != sup.bothcohorts$clust])   
mean(sup.nocohorts$prob[sup.nocohorts$clust == sup.bothcohorts$clust])
mean(sup.nocohorts$prob[sup.nocohorts$clust != sup.bothcohorts$clust] > 0.8)   
mean(sup.nocohorts$prob[sup.nocohorts$clust == sup.bothcohorts$clust] > 0.8)

# posterior probs in changed:
mod = summary(lm(sup.nocohorts$prob[sup.nocohorts$clust != sup.bothcohorts$clust] ~ 1))
crit = qt(0.975, df = mod$df[2])
mod$coef[1,1] + c(0, -crit, crit)*mod$coef[1,2]

# posterior probs in unchanged:
mod = summary(lm(sup.nocohorts$prob[sup.nocohorts$clust == sup.bothcohorts$clust] ~ 1))
crit = qt(0.975, df = mod$df[2])
mod$coef[1,1] + c(0, -crit, crit)*mod$coef[1,2]


#### list the cells that changed ------------------------

changed = names(sup.nocohorts$clust)[sup.nocohorts$clust != sup.bothcohorts$clust]
cbind(sup.nocohorts$clust[changed], sup.bothcohorts$clust[changed], sup.nocohorts$prob[changed])[
  1:100,]

cbind(sup.nocohorts$clust[changed], sup.bothcohorts$clust[changed], sup.nocohorts$prob[changed])[
  sup.bothcohorts$clust[changed] == "Glomerular.endothelium",]  #sup.nocohorts$clust[changed] == "Podocyte",]

look = c( "c_5_8_3315","c_5_9_147", "c_5_8_2149", "c_5_8_2602") # "c_5_8_2863",
plot(xy, col = "grey70", cex = 0.1, pch = 16, asp = 1, 
     xlim = c(17,18), ylim = c(-25.5,-23.95))
points(xy[look, ], col = 1:4, pch = 16, cex = 2)

