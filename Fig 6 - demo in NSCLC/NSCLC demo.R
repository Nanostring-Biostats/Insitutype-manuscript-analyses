library(pheatmap)
library(RColorBrewer)
library(MLEcell)
library(binom)
library(vioplot)
library(scales)
library(rlist)
rm(list = ls())

source("getcellneighbors.R")

#### load data: ---------------

# load fundamental data:
annot <- readRDS("../data/NSCLC/cell annotations.rds")
annot$tissue2 = annot$tissue
annot$tissue2[grepl("Lung5", annot$tissue)] = "Lung5"
xy <- readRDS("../data/NSCLC/cell positions.rds")
raw <- readRDS("../data/NSCLC/raw_counts.rds")
annot$negmean <- readRDS("../data/NSCLC/negmean.rds")
um <- readRDS("../data/NSCLC/umap projection.rds")
um[annot$tissue == "Lung13", 1] = um[annot$tissue == "Lung13", 1] + 1
um[annot$tissue == "Lung12", 2] = um[annot$tissue == "Lung12", 2] - 1
counts <- t(as.matrix(raw))
rm(raw)

# load immune profiles:
ioprofiles <- readRDS("../data/NSCLC/ioprofiles.RDS")

# colors:
iocolors = readRDS(file = "../data/NSCLC/iocolors.RDS")


#### cohort based on image data ----------------------

## extract cells' image data:
set.seed(0)
ifdata = as.matrix(annot[, c("Area","AspectRatio","Mean.CD298","Mean.G","Mean.R","Mean.Y","Mean.DAPI")])
cohort = fastCohorting(mat = ifdata)
save(cohort, file = "results/cohort.RData")


#### semi-supervised cell typing -----------------------
# choose cluster number:
nclust = chooseClusterNumber(counts = counts,
                             neg = annot$negmean, 
                             bg = annot$background,
                             n_clusts = 6:20,
                             fixed_profiles = ioprofiles) # returns 11
plot(nclust$n_clusts, nclust$bic)
plot(nclust$n_clusts, nclust$aic)

# perform clustering:
res = insitutype(counts = counts,
                 neg = annot$negmean,
                 n_clusts = nclust$best_clust_number,
                 reference_profiles = ioprofiles,
                 cohort = cohort)
save(res, file = "results/NSCLC semi-reservised cell typing results.RData")

# set colors:
cols = iocolors[is.element(names(iocolors), unique(res$clust))]
cols[setdiff(unique(res$clust), names(cols))] = c(brewer.pal(11, "Set3"))

pheatmap(table(res$clust, annot$tissue))

# umap:
plot(um, cex = 0.1, pch = 16, col = cols[res$clust])
legend("bottomright", pch = 16, col = cols, legend = names(cols), cex = 0.5)

# zoomed umaps:
use = res$clust == "d"
plot(um, cex = 0.1, pch = 16, col = alpha(cols[res$clust], 0.3))
points(um[use, ], pch = 16, col = cols[res$clust[use]], cex = 0.3)

# flightpath:
fp0 = MLEcell::flightpath_layout(logliks = res$logliks, profiles = res$profiles)
plot(fp0$cell, pch = 16, cex = 0.1, col = cols[res$clust], xlab = "", ylab = "", xaxt = "n", yaxt = "n")
text(fp0$clustpos, rownames(fp0$clustpos), cex = 0.5)


### how to refine clusters:
# perform renaming:
res1 = refineClusters(logliks = res$logliks,
                      merges = c(#"macrophage" = "macrophage SPP1 pos",
                                 "k" = "macrophage SPP1 low",
                                 #"plasmablast" = "plasmablast1",
                                 "e" = "plasmablast IgG low",
                                 "i" = "fibroblast DUSP5 high",
                                 "a" = "PanCK+ 1",
                                 "b" = "PanCK+ 2",
                                 "c" = "PanCK+ 3",
                                 "f" = "PanCK+ 4",
                                 "g" = "PanCK+ 5",
                                 "j" = "PanCK+ 6",
                                 "h" = "healthy epithelial",
                                 "d" = "non-specific stroma"), 
                      counts = counts, neg = annot$negmean)
save(res1, file = "results/NSCLC renamed semi-sup results.RData")


# now delete cluster d:
res2 = refineClusters(logliks = res1$logliks,
                      merges = NULL, 
                      to_delete = "non-specific stroma", 
                      counts = counts, neg = annot$negmean)
save(res2, file = "results/NSCLC refined semi-sup results.RData")

#### plots ----------------------------------

# colors:
cols = iocolors[is.element(names(iocolors), unique(res2$clust))]
cols[setdiff(unique(res2$clust), names(cols))] = c(brewer.pal(11, "Set3"), brewer.pal(8, "Set2"))[1:length(setdiff(unique(res2$clust), names(cols)))]
cols["plasmablast IgG low"] = "dodgerblue3"
cols["macrophage SPP1 low"] = "chartreuse3"
cols["NK"] = "grey30"
cols["NK"] = "grey30"
cols["B-cell"] = "blue"
cols["non-specific stroma"] = "grey70"

# flightpath:
fp = MLEcell::flightpath_layout(logliks = res1$logliks, profiles = res1$profiles)
png("results/flightpath - renamed but not deleted.png", width = 5, height = 5, units = "in", res = 400)
par(mar = c(0,0,0,0))
plot(fp$cell, pch = 16, cex = 0.1, col = cols[res1$clust], 
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", xlim = c(-12.1,10.35))
text(fp$clustpos, 
     paste0(rownames(fp$clustpos), "\n", c("(new)", "(ref)")[
       1 + is.element(rownames(fp$clustpos), colnames(ioprofiles))]), 
       cex = 0.75)
dev.off()

# xy plots:

png("results/cell type xy.png", width = 14, height = 14, units = "in", res = 300)
plot(xy, pch = 16, cex = 0.5, col = cols[res2$clust])
abline(v = c(2,4,6,8,10,12,14), col = "grey")
abline(h = -c(2,4,6,8,10,12,14), col = "grey")
dev.off()

png("results/cell type xy - TLS.png", width = 3, height = 3, units = "in", res = 300)
xlim = c(0.5, 1.5)
ylim = c(-8.5, -7.5)
par(mar = c(0,0,0,0))
plot(xy, pch = 16, cex = 0.3, col = cols[res2$clust],
     xlim = xlim, ylim = ylim,
     asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
#vals = c(115196, 116065, 119178, 119306, 135670, 135701, 135853, 136065, 137161)
vals = c(114432, 114876, 115231, 115322, 115775, 116264, 116649, 116872, 119079,
         119281, 135552, 135748, 135802, 136066, 137161)
polygon(sweep(xy[vals[chull(xy[vals, ])], ], 2, c(.093,-.07), "+"))

lines(xlim[2] + c(-0.6, -0.1),
      rep(ylim[1] + 0.1 * diff(range(ylim)), 2))
text(xlim[2] - 0.3, ylim[1] + 0.1 * diff(range(ylim)) - 0.1, "0.5 mm")
dev.off()


png("results/cell type xy - neutrophil pocket.png", width = 3, height = 3, units = "in", res = 300)
xlim = c(2, 3.5)
ylim = c(-11.35, -10.1)
par(mar = c(0,0,0,0))
plot(xy, pch = 16, cex = 0.2, col = cols[res2$clust],
     xlim = xlim, ylim = ylim,
     asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
#vals = c(455124, 455356, 457639, 461559, 471204, 471502,
#         471528, 472195, 474791, 475702, 478302, 478457)
#polygon(sweep(xy[vals[chull(xy[vals, ])], ], 2, c(.093,-.07), "+"))

lines(xlim[2] + c(-0.6, -0.1),
      rep(ylim[1] + 0.1 * diff(range(ylim)), 2))
text(xlim[2] - 0.3, ylim[1] + 0.1 * diff(range(ylim)) - 0.1, "0.5 mm")
dev.off()




png("results/cell type xy - invading ts and macs.png", width = 3, height = 3, units = "in", res = 300)
xlim = c(13.7, 14.4)
ylim = c(-1, -0.1)
par(mar = c(0,0,0,0))
plot(xy, pch = 16, cex = 0.4, col = cols[res2$clust],
     xlim = xlim, ylim = ylim,
     asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
vals = c(386590, 386739, 386924, 387003, 387049, 387142, 387157, 388536)
polygon(sweep(xy[vals[chull(xy[vals, ])], ], 2, c(.093,-.04), "+"))
lines(xlim[2] + c(-0.6, -0.1),
      rep(ylim[1] + 0.1 * diff(range(ylim)), 2))
text(xlim[2] - 0.3, ylim[1] + 0.1 * diff(range(ylim)) - 0.1, "0.5 mm")
dev.off()


png("results/cell type xy - plasmablasts.png", width = 3, height = 3, units = "in", res = 300)
xlim = c(8.25,9.1)
ylim = c(-12, -11.15)
par(mar = c(0,0,0,0))
plot(xy, pch = 16, cex = 0.4, col = cols[res2$clust],
     xlim = xlim, ylim = ylim,
     asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
vals = c(386590, 386739, 386924, 387003, 387049, 387142, 387157, 388536)
polygon(sweep(xy[vals[chull(xy[vals, ])], ], 2, c(.093,-.04), "+"))
lines(xlim[2] + c(-0.6, -0.1),
      rep(ylim[1] + 0.1 * diff(range(ylim)), 2))
text(xlim[2] - 0.3, ylim[1] + 0.1 * diff(range(ylim)) - 0.1, "0.5 mm")
dev.off()


#### plot the two macrophage populations ---------------------------------------

png("results/cell type xy - macrophages - Lung9 - square.png", width = 6.5, height = 6.5, units = "in", res = 300)
xlim = c(0.5, 3.5)
ylim = c(-13, -10)
par(mar = c(0,0,0,0))
plot(xy, pch = 16, cex = 0.3, col = alpha(cols[res2$clust], 0.6),
     xlim = xlim, ylim = ylim,
     asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
points(xy[res2$clust == "macrophage", ], pch = 16, col = cols["macrophage"], cex = 0.5)
points(xy[res2$clust == "macrophage SPP1 low", ], pch = 16, col = cols["macrophage SPP1 low"], cex = 0.5)

lines(xlim[2] + c(-0.6, -0.1),
      rep(ylim[1] + 0.1 * diff(range(ylim)), 2))
text(xlim[2] - 0.3, ylim[1] + 0.1 * diff(range(ylim)) - 0.1, "0.5 mm")
dev.off()

png("results/cell type xy - macrophages - Lung9.png", width = 6.25, height = 3, units = "in", res = 300)
xlim = c(0.15, 3.8)
ylim = -12.5 + c(0, 1.7)
par(mar = c(0,0,0,0))
plot(xy, pch = 16, cex = 0.3, col = alpha(cols[res2$clust], 0.6),
     xlim = xlim, ylim = ylim,
     asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
points(xy[res2$clust == "macrophage", ], pch = 16, col = cols["macrophage"], cex = 0.5)
points(xy[res2$clust == "macrophage SPP1 low", ], pch = 16, col = cols["macrophage SPP1 low"], cex = 0.5)

lines(xlim[2] + c(-0.6, -0.1),
      rep(ylim[1] + 0.1 * diff(range(ylim)), 2))
text(xlim[2] - 0.3, ylim[1] + 0.1 * diff(range(ylim)) - 0.1, "0.5 mm")
dev.off()



png("results/cell type xy - macrophages - Lung6.png", width = 6.25, height = 3, units = "in", res = 300)
ylim = c(14, 15.5)-0.05
xlim = c(-3.5, -0.5)
par(mar = c(0,0,0,0))
plot(xy[, 2:1], pch = 16, cex = 0.3, col = alpha(cols[res2$clust], 0.6),
     xlim = xlim, ylim = ylim,
     asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
points(xy[res2$clust == "macrophage", 2:1], pch = 16, col = cols["macrophage"], cex = 0.5)
points(xy[res2$clust == "macrophage SPP1 low", 2:1], pch = 16, col = cols["macrophage SPP1 low"], cex = 0.5)

lines(xlim[2] + c(-0.6, -0.1),
      rep(ylim[1] + 0.1 * diff(range(ylim)), 2))
text(xlim[2] - 0.3, ylim[1] + 0.1 * diff(range(ylim)) - 0.1, "0.5 mm")
dev.off()


# summary stats:
tab = table(res2$clust, annot$tissue)[c("macrophage", "macrophage SPP1 low"), ]
round(sweep(tab, 2, colSums(tab), "/"), 3)

#### UMAP and flightpath -----------------------------------

# umap:
png("results/umap.png", width = 6.5, height = 5, units = "in", res = 400)
par(mar = c(0,0,0,0))
plot(um, pch = 16, cex = 0.1, col = cols[res2$clust], 
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", asp = 1)
for (tiss in unique(annot$tissue2)) {
  use = annot$tissue2 == tiss
  if (is.element(tiss, c("Lung5", "Lung6", "Lung9"))) {
    ypos = 8
  } else {
    ypos = -13
  }
  
  text(median(um[use,1]), ypos, tiss)
}
abline(h = -12, col = "grey30")
abline(v = 22, col = "grey30")
abline(v = 7.5, col = "grey30")
dev.off()


# confidence:
svg("results/posterior probs vioplot.svg",
    height = 3.5, width = 6)
par(mar = c(10,4,.5,.5))
vioplot(res2$prob ~ res2$clust, las = 2, ylim = c(0.4,1),
        col = cols[levels(as.factor(res2$clust))],
        xlab = "", ylab = "Posterior probability")
dev.off()

barplot(by(res2$prob > 0.95, res2$clust, mean), las = 2)
names(which(by(res2$prob > 0.9, res2$clust, mean) > 0.8))
names(which(by(res2$prob > 0.9, res2$clust, mean) <= 0.8))




#### plot results for cell pairs in question for supplement: --------------------------------
cpairs = list("fibroblast" = c("fibroblast", "fibroblast DUSP5 high"),
              "plasmablast" = c("plasmablast IgG low", "plasmablast"),
              "macrophage" = c("macrophage SPP1 low", "macrophage"))
for (name in names(cpairs)) {
  
  tempcols.legend = c("grey70", "orange", "dodgerblue2")
  names(tempcols.legend) = c("cancer", cpairs[[name]])
  
  svg(paste0("results/nsclc pairs supplement - ", name, " - legend.svg"), width = 3, height = 1.5)
  par(mar = c(0,0,0,0))
  frame()
  legend("center", col = tempcols.legend, legend = names(tempcols.legend), pch = 16)
  dev.off()
  
  tempcols = c("orange", "dodgerblue2")
  names(tempcols) = c(cpairs[[name]])
  tempcols[unique(res2$clust)[grepl("PanCK", unique(res2$clust))]] = "grey60"
  
  png(paste0("results/nsclc pairs supplement - ", name, "scatter-log.png"), width = 5, height = 5, units = "in", res = 400)
  par(mar = c(5,4,0.5,0.5))
  plot(res2$profiles[, cpairs[[name]]], log = "xy", col = 0)
  text(res2$profiles[, cpairs[[name]]], rownames(res2$profiles), cex = 0.6)
  dev.off()
  
  png(paste0("results/nsclc pairs supplement - ", name, "scatter-linear.png"), width = 5, height = 5, units = "in", res = 400)
  par(mar = c(5,4,0.5,0.5))
  plot(res2$profiles[, cpairs[[name]]], col = 0)
  text(res2$profiles[, cpairs[[name]]], rownames(res2$profiles), cex = 0.6)
  dev.off()
  
  png(paste0("results/nsclc pairs supplement - ", name, " umap.png"), 
      width = 6.5, height = 5, units = "in", res = 400)
  par(mar = c(0,0,0,0))
    show = is.element(res2$clust, cpairs[[name]])
  plot(um[show, ], pch = 16, cex = 0.1, col = tempcols[res2$clust[show]], 
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", asp = 1)
  for (tiss in unique(annot$tissue2)) {
    use = annot$tissue2 == tiss
    if (is.element(tiss, c("Lung5", "Lung6", "Lung9"))) {
      ypos = 8
    } else {
      ypos = -13
    }
    
    text(median(um[use,1]), ypos, tiss)
  }
  abline(h = -12, col = "grey30")
  abline(v = 22, col = "grey30")
  abline(v = 7.5, col = "grey30")
  legend(23, -20, pch = 16, col = tempcols, legend = names(tempcols), bty = "n")
  dev.off()
  
  
  png(paste0("results/nsclc pairs supplement - ", name, " xy.png"), 
      width = 7, height = 7, units = "in", res = 400)
  par(mar = c(0,0,0,0))
  show = is.element(res2$clust, cpairs[[name]]) | grepl("PanCK", res2$clust)
  plot(xy[show, ], pch = 16, cex = 0.1, col = tempcols[res2$clust[show]], 
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", asp = 1)
  for (tiss in unique(annot$tissue2)) {
    use = annot$tissue2 == tiss
    text(median(xy[use,1]), max(xy[use, 2]), tiss)
  }
  dev.off()
}


#### plots to zoom in on vignettes -----------------------------


# big profiles matrix:
pdf("results/big profiles heatmap.pdf", height = 30, width = 6)
pheatmap(sweep(res$profiles, 1, pmax(apply(res$profiles, 1, max), 1), "/"), 
         col = colorRampPalette(c("white", "darkblue"))(100),
         fontsize_row = 2)
dev.off()

# subclusters of a tumor:
show = is.element(res$clust, c("g", "a"))
plot(xy[show, ], col = cols[res$clust[show]], pch = 16, cex = 0.2,
     xlim = c(10,16), ylim = c(-4,0))
legend("topright", pch=16, col = cols[c("a", "g")], legend = c("a", "g"))
plot(1e-3+res$profiles[, c("g", "a")], log = "xy", col = 0)
text(1e-3+res$profiles[, c("g", "a")], rownames(res$profiles), cex = 0.5)
plot(1e-3+res$profiles[, c("g", "a")], col = 0)
text(1e-3+res$profiles[, c("g", "a")], rownames(res$profiles), cex = 0.5)

# look in space:
plot(xy[show, ], col = c("blue", "orange")[as.numeric(as.factor(res$clust[show]))], 
     pch = 16, cex = 0.2, xlim = c(11,15), ylim = c(-4,0))
# color by panck stain:
plot(xy[show, ], col = colorRampPalette(c("darkviolet", "grey90", "forestgreen"))(101)[
  1 + pmin(round(100*annot$Mean.G / quantile(annot$Mean.G, 0.95)),100)][show], 
  pch = 16, cex = 0.2, xlim = c(11,15), ylim = c(-4,0))


#### analyze panck+ clusters to assess error rate: ---------------------------

# get all panck+ clusters:
boxplot(annot$Mean.G ~ res$clust, las = 2)
boxplot(log2(1+annot$Mean.G) ~ res2$clust, las = 2)
panckclusts = names(which(by(annot$Mean.G, res2$clust, median) > 2^11))
panckclusts

tab = table(res2$clust, annot$tissue2)[panckclusts, ]
write.csv(sweep(tab, 1, rowSums(tab), "/"), file = "results/panck table data.csv")
svg("results/panck table heatmap.svg", width = 3.5, height = 4)
pheatmap(sweep(tab, 1, rowSums(tab), "/")[
  c("healthy epithelial", paste0("PanCK+ ", c(4,1,5,2,6,3))), paste0("Lung", c(5,6,9,12,13))], 
         display_numbers = T, legend = F,
         cluster_rows = F, cluster_cols = F,
         col = colorRampPalette(c("white", "cornflowerblue"))(100),
         main = "")#"Proportion of PanCK+\nclusters in each tissue")
dev.off()


tab200 = table(res2$clust[annot$totalcounts > 200], annot$tissue2[annot$totalcounts > 200])[panckclusts, ]
round(sweep(tab200, 1, rowSums(tab200), "/"), 3)

# error rate for each cluster: 
subsets = list(all = rep(TRUE, nrow(counts)),
               verylow = (rowSums(counts) <= 50),
               low = (rowSums(counts) > 50) & (rowSums(counts) <= 100),
               med = (rowSums(counts) > 100) & (rowSums(counts) <= 200),
               high = (rowSums(counts) > 200))
matches = c("PanCK+ 1" = "Lung6",
            "PanCK+ 2" = "Lung9",
            "PanCK+ 3" = "Lung13",
            "PanCK+ 4" = "Lung5",
            "PanCK+ 5" = "Lung6",
            "PanCK+ 6" = "Lung12")

errmat = lci = uci = matrix(NA, length(matches), length(subsets),
                            dimnames = list(names(matches), names(subsets)))
for (name in names(subsets)) {
  tab = table(res2$clust[subsets[[name]]], annot$tissue2[subsets[[name]]])[names(matches), matches]
  n.right = diag(tab)
  n = rowSums(tab[, !duplicated(colnames(tab))])
  errmat[, name] = 1 - (n.right / n)
  for (i in 1:nrow(tab)) {
    ci = t(binom.confint(x = n[i] - n.right[i], n = n[i], method = "wilson")[, c("lower", "upper")])
    lci[i, name] = ci[1]
    uci[i, name] = ci[2]
  }
  tab = tab[, !duplicated(colnames(tab))]
  pheatmap(sweep(tab, 1, rowSums(tab), "/"), main = name, display_numbers = T,
           cluster_rows = F, cluster_cols = F,
           col = colorRampPalette(c("white", "cornflowerblue"))(100))
}

colnames(errmat) = c("all cells",  "0-50\ntranscripts",  "51-100\ntranscripts", 
                     "101-200\ntranscripts", ">200\ntranscripts")
errcols = cols[rownames(errmat)]
svg("results/panck err barplot.svg", width = 7, height = 4)
par(mar = c(4.5,5,.5,.2))
bp = barplot(errmat, beside = T, col = errcols, cex.lab = 1.2,
             ylab = "Prop. assigned to wrong tissues")
for (i in 1:nrow(errmat)) {
  for (j in 1:ncol(errmat)) {
    lines(rep(bp[i,j], 2), c(lci[i,j], uci[i,j]))
  }
}
legend("topleft", fill = errcols, legend = names(errcols), cex = 0.85)
dev.off()


#### plot results ---------------------------------

# umap: 
plot(um, pch = 16, cex = 0.2, col = cols[res$clust], asp = 1)

# flightpath:
fp = MLEcell::flightpath_layout(logliks = res$logliks, profiles = res$profiles)
png("results/flightpath.png", width = 8, height = 8, units ="in", res = 400)
plot(fp$cell, pch = 16, cex = 0.3, col = cols[res$clust], xlab = "", ylab = "", xaxt = "n", yaxt = "n")
text(fp$clustpos, rownames(fp$clustpos), cex = 0.5)
dev.off()
