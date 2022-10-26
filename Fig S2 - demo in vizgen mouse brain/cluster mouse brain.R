rm(list = ls())
library(MLEcell)
library(pheatmap)
library(scales)
library(tictoc)
library(viridis)
library(RColorBrewer)

# load data:
load("../../SMI investigations/vizgen showcase/complete vizgene showcase data - counts and metadata.RData")

# cluster:
set.seed(0)
tic()
out = insitutype(counts = raw[, !grepl("Blank", colnames(raw))],
                 neg = rowMeans(raw[, grepl("Blank", colnames(raw))]),
                 n_clusts = 40)
toc() # 13771.301 sec elapsed
save(out, file = "results/clustering results.RData")


## plots
png("results/mouse brain xy plots.png", height = 8, width = 8, units = "in", res = 300)
par(mar = c(0,0,0,0))
par(mfrow = c(7, 6))
for (name in unique(out$clust)) {
  use = (out$clust == name) & (out$prob == 1)
  plot(annot$center_x[use], annot$center_y[use], pch = 16, col = "darkorchid4", cex = 0.05,
       asp = 1, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  text(median(range(annot$center_x[use])), max(annot$center_y[use]), name)
}
dev.off()

fp = MLEcell::flightpath_layout(logliks = out$logliks)
save(fp, file = "results/flightpath.RData")

png('results/flightpath2.png', width = 5, height = 5, units = "in", res = 300)
par(mar = c(0,0,0,0))
plot(fp$cellpos, pch = 16, cex = 0.2, col = alpha("darkorchid2", 0.1),# col = "darkorchid2", # 
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
text(fp$clustpos, rownames(fp$clustpos), cex = 0.6)#, col = "black")
dev.off()


mat = sweep(out$profiles, 1, pmax(apply(out$profiles,1,max), 0.05), "/")
png('results/profiles.png', width = 5, height = 5, units = "in", res = 300)
pheatmap(mat, col = viridis_pal(option = "B")(100), legend = F,
         show_rownames = F)
dev.off()
