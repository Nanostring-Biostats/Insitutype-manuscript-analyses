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
cols = readRDS("../data/kidney cols.RDS")
refprofiles = readRDS("../data/HCA kidney profiles plus nstg immune profiles.RDS")


#### cohort based on image and location data ----------------------

## extract cells' image data:
set.seed(0)
ifdata = as.matrix(annot[, c("Area","AspectRatio","Mean.CD298","Mean.PanCK","Mean.CD45","Mean.CD20","Mean.DAPI")])

## extract cells' spatial context:
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

## cohort based on IF and spatial context:
set.seed(0)
bothcohort = fastCohorting(mat = cbind(ifdata, neighborhoodPCs))



#### supervised cell typing -----------------------
sup = insitutypeML(counts = counts,
                   neg = annot$negmean,
                   reference_profiles = refprofiles,
                   cohort = bothcohort)
save(sup, file = "results/kidney supervised cell typing results.RData")



#### plot results ---------------------------------

# umap: 
plot(um, pch = 16, cex = 0.2, col = cols[sup$clust], asp = 1)
for (cell in unique(sup$clust)) {
  iscell = sup$clust == cell
  text(median(um[iscell, 1]), median(um[iscell, 2]), cell, cex = 0.5)
}

# flightpath:
fp = flightpath_layout(logliks = sup$logliks, profiles = sup$profiles)
rownames(fp$clustpos)[rownames(fp$clustpos) == "MNP.b.non.classical.monocyte.derived"] = "Non.classical.monocyte"
rownames(fp$clustpos)[rownames(fp$clustpos) == "MNP.a.classical.monocyte.derived"] = "Classical.monocyte"
rownames(fp$clustpos)[rownames(fp$clustpos) == "Descending.vasa.recta.endothelium"] = "Desc.vasa.recta.endo"
rownames(fp$clustpos)[rownames(fp$clustpos) == "Ascending.vasa.recta.endothelium"] = "Asc.vasa.recta.endo"
rownames(fp$clustpos)[rownames(fp$clustpos) == "Peritubular.capillary.endothelium.1"] = "Peritub.capillary.endo.1"
rownames(fp$clustpos)[rownames(fp$clustpos) == "Peritubular.capillary.endothelium.2"] = "Peritub.capillary.endo.2"
rownames(fp$clustpos)[rownames(fp$clustpos) == "Thick.ascending.limb.of.Loop.of.Henle"] = "Loop.of.henle"
rownames(fp$clustpos)[rownames(fp$clustpos) == "Distinct.proximal.tubule.1"] = "Distinct.prox.tub.1"
rownames(fp$clustpos)[rownames(fp$clustpos) == "Distinct.proximal.tubule.2"] = "Distinct.prox.tub.2"
rownames(fp$clustpos)[rownames(fp$clustpos) == "Proliferating.Proximal.Tubule"] = "Prolif.prox.tub"


png("results/flightpath.png", width = 8, height = 7, units ="in", res = 400)
par(mar = c(0,0,0,0))
plot(fp$cell, pch = 16, cex = 0.3, col = cols[sup$clust], xlab = "", ylab = "", xaxt = "n", yaxt = "n",
     xlim = c(-12, 11.5))
text(fp$clustpos, rownames(fp$clustpos), cex = 1)
dev.off()

# xy:
png("results/cell type xy.png", width = 8, height = 14, units ="in", res = 400)
plot(xy, pch = 16, col = cols[sup$clust], cex = 0.2, 
     asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n") 
legend("topright", pch = 16, col = cols, legend = names(cols), cex = 0.7)
dev.off()

# revised xy:
plot(xy)
abline(8, -2)
abline(h = -21)
# rotate:
theta =  2*pi * 0.17
xy2 = xy
xy2[, 1] = xy[, 1] * cos(theta) - xy[, 2] * sign(theta)
xy2[, 2] = xy[, 1] * sin(theta) + xy[, 2] * cos(theta)
plot(xy2[sample(1:nrow(xy), 1000), ], asp = 1)
# condense:
inds1 = (xy2[,2] > 2.5) & (xy2[, 1] < 29)
inds2 = (xy2[,2] > 2.5) & (xy2[, 1] > 29)
xy2[inds1, 1] = xy2[inds1, 1] + 4
xy2[inds1, 2] = xy2[inds1, 2] -0.65
xy2[inds2, 1] = xy2[inds2, 1] + -2
xy2[inds2, 2] = xy2[inds2, 2] + 0.5
plot(xy2[sample(1:nrow(xy), 1000), ], asp = 1)

png("results/cell type xy2.png", width = 8, height = 4.25, units ="in", res = 400)
par(mar = c(0,0,0,0))
plot(xy2, pch = 16, col = cols[sup$clust], cex = 0.2, 
     asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n") 
lines(c(31, 32), rep(1, 2))
text(31.5, 0.9, "1 mm")
dev.off()



#### subtype myofibroblasts -------------
ismyo = sup$clust == "Myofibroblast"
set.seed(0)
myoclust = insitutype(counts = counts[ismyo, ],
                      neg = annot$negmean[ismyo],
                      n_clusts = 3)
save(myoclust, file = "results/myo subclusters.RData")

clusts = sup$clust
clusts[names(myoclust$clust)] = paste0("Myofibroblast_", myoclust$clust)
clusts[clusts == "Distinct.proximal.tubule.1"] = "Parietal.epithelium"



# plot myo subclusters in space:
tempclust = replace(sup$clust, TRUE, "other")
tempclust[names(myoclust$clust)] = paste0("Myofibroblast_", myoclust$clust)
tempcols = c("other" = alpha("grey70", 0.5),
             "Myofibroblast_a" = "red",
             "Myofibroblast_b" = "blue",
             "Myofibroblast_c" = "darkorchid3")#"Myofibroblast_d" = "forestgreen")

png("results/myofibroblast subcluster xy2.png", width = 2.5, height = .8, units ="in", res = 400)
isother = tempclust == "other"
par(mar = c(0,0,0,0))
plot(xy2[isother, ], pch = 16, col = tempcols[tempclust[isother]], cex = 0.2, 
     asp = 1, 
     xlab = "", ylab = "", xaxt = "n", yaxt = "n",
     xlim = c(26.5, 29), ylim = c(1.3, 2.1)) 
points(xy2[!isother, ], pch = 16, col = tempcols[tempclust[!isother]], cex = 0.2)
legend("topleft", pch = 16, col = tempcols, legend = names(tempcols), cex = 0.35)
lines(c(28, 29), rep(1.2, 2))
text(31.5, 0.9, "1 mm")
dev.off()


png("results/myofibroblast profiles pairs plot.png", res = 400, units = "in", width =7, height = 7)
pairs(myoclust$profiles, lab = paste0("Myofibroblast_", letters[1:3]),
      pch = 16, col = "darkblue", cex = 0.5)
dev.off()


#### proximal tubule subtypes --------------------------------

tempclust = replace(sup$clust, TRUE, "other")
tempclust[grepl("ubule", sup$clust)] = sup$clust[grepl("ubule", sup$clust)]
tempcols = c("other" = alpha("grey70", 0.5),
             "Distinct.proximal.tubule.1" = "dodgerblue4",
             "Distinct.proximal.tubule.2" = "blue",
             "Proximal.tubule" = "cornflowerblue",
             "Proliferating.Proximal.Tubule" = "red",
             "Connecting.tubule" = "orange")

png("results/proximal tubule subcluster xy.png", 2.5, height = .8, units ="in", res = 400)
isother = tempclust == "other"
par(mar = c(0,0,0,0))
plot(xy2[isother, ], pch = 16, col = tempcols[tempclust[isother]], cex = 0.2, 
     asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n",
     xlim = c(26.5, 29), ylim = c(1.3, 2.1)) 
points(xy2[!isother, ], pch = 16, col = tempcols[tempclust[!isother]], cex = 0.2)
legend("topleft", pch = 16, col = tempcols, 
       #legend = names(tempcols), 
       legend = c("other", "Distinct PT 1", "Distinct PT 2", "PT", "Prolif PT", "Connecting tub."), 
       cex = 0.35)
lines(c(28, 29), rep(1.2, 2))
text(31.5, 0.9, "1 mm")
dev.off()


png("results/prox tubules profiles pairs plot.png", res = 400, units = "in", width =7, height = 7)
par(mar = c(0,0,0,0))
pairs(refprofiles[is.element(rownames(refprofiles), colnames(counts)), setdiff(names(tempcols), "other")],
      pch = 16, col = "darkblue", cex = 0.5,
      labels = c("Distinct PT 1", "Distinct PT 2", "PT", "Prolif PT", "Connecting tub."))
dev.off()
