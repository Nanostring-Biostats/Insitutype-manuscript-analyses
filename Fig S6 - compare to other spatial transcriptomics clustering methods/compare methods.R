library(scales)
library(pheatmap)
library(RColorBrewer)

# plot runtimes:
times = c(479.675, 25772.881) / 60 / 60
names(times) = c("Insitutype", "BASS")

svg("results/runtimes barplot.svg", width = 2.5, height = 4)
par(mar = c(4,5,0.1,0.1))
bp = barplot(times, ylim = c(0, max(times) + 0.5),
             col = alpha("dodgerblue4", 0.5),
             ylab = "Hours to cluster 97,809 cells",
             cex.lab = 1.2)
text(bp, times + 0.2, round(times, 2))
dev.off()

annot <- readRDS("../data/NSCLC/cell annotations.rds")
annot$negmean <- readRDS("../data/NSCLC/negmean.rds")
xy <- readRDS("../data/NSCLC/cell positions.rds")

# just take 1 slide:
keep = annot$tissue == "Lung5-5"
annot = annot[keep, ]
xy = xy[keep, ]


load("results/BASS cluster assignments.RData") # loads "bassclusts"
load("results/insitutype results.RData") # loads "res"
clusts = res$clust

# load semi-supervised results to help name clusters:
load("~/hades/pdanaher/Insitutype manuscript/demo in NSCLC/results/NSCLC refined semi-sup results.RData")
# (loads "res2")
semi = res2$clust[names(clusts)]
lymphoidcells = c(unique(semi)[grepl("T ", unique(semi))],
                  "Treg", "B-cell", "NK", "pDC")
semi[is.element(semi, lymphoidcells)] = "lymphoid"

pheatmap(table(semi, bassclusts), display_numbers = T,
         col = colorRampPalette(c("white", "cornflowerblue"))(100),
         breaks = seq(0, 2000, length.out = 101))
pheatmap(table(semi, clusts), display_numbers = T,
         col = colorRampPalette(c("white", "cornflowerblue"))(100),
         breaks = seq(0, 2000, length.out = 101))

# agreements between clustering methods:
svg("results/agreement heatmap.svg", width = 4, height = 4)
pheatmap(table(clusts, bassclusts), display_numbers = T,
         col = colorRampPalette(c("white", "cornflowerblue"))(100),
         breaks = seq(0, 2000, length.out = 101))
dev.off()


# align the clusts that agree:
bassclusts[bassclusts == "4"] = "l"
bassclusts[bassclusts == "9"] = "a"
bassclusts[bassclusts == "5"] = "b"
bassclusts[bassclusts == "3"] = "d"
bassclusts[bassclusts == "12"] = "h"
bassclusts[bassclusts == "7"] = "j"
bassclusts[bassclusts == "10"] = "e"
bassclusts[bassclusts == "8"] = "i"
bassclusts[bassclusts == "6"] = "g"


cols = brewer.pal(12, "Set3")
names(cols) = unique(clusts)
bcols = cols[is.element(names(cols), unique(bassclusts))]
bcols[setdiff(unique(bassclusts), names(cols))] = setdiff(cols, bcols)

png("results/bass vs insitutype.png", width = 8, height = 7, units = "in", res = 400)
par(mfrow = c(2, 2))
par(mar = c(0.3,0.3,0.3,0.3))
# plot insitutype results:
plot(um, pch = 16, cex = 0.2, col = cols[clusts],
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
for (id in names(table(semi)[table(semi) > 900])) {
  text(median(um[semi == id, 1]), median(um[semi == id, 2]), id)
}
plot(xy, pch = 16, cex = 0.2, col = cols[clusts],
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", asp = 1)
# plot BASS results:
plot(um, pch = 16, cex = 0.2, col = bcols[bassclusts],
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
for (id in names(table(semi)[table(semi) > 900])) {
  text(median(um[semi == id, 1]), median(um[semi == id, 2]), id)
}
plot(xy, pch = 16, cex = 0.2, col = bcols[bassclusts],
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", asp = 1)
dev.off()


