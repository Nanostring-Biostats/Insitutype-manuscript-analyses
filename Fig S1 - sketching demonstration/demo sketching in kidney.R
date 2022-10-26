
library(MLEcell)
library(RColorBrewer)
library(pheatmap)
library(scales)
library(uwot)

source("geoSketch.R")

#### load data ---------------------------
load("../data/lupus nephritis SP11_1139 data.RData")
refprofiles = readRDS("../data/kidney HCA reference matrix.RDs")
clust = readRDS("../data/lupus nephritis cell types.RDS")
cols = readRDS("../data/kidney cols.RDS")


#### run sketching -------------------------
set.seed(0)
sdat = prepDataForSketching(counts)
res = geoSketch(X = sdat, 
                N = 2000,
                alpha=0.1,
                max_iter=200,
                returnBins=FALSE,
                minCellsPerBin = 1,
                seed=NULL)

#### umap ---------------------------------
um = uwot::umap(sqrt(sweep(counts, 1, rowSums(counts), "/")))
rownames(um) = rownames(counts)
png("sketching umap.png", units = "in", res = 300, width = 6, height = 6)
par(mar = c(0,0,0,0))
plot(um, pch = 16, col = alpha("grey60", 0.5), cex = 0.2, asp = 1,
     xaxt = "n", yaxt = "n", xlab = "", ylab = "")
points(um[res, ], pch = 16, cex = 0.25, col = alpha("gold", 0.7))
legend("bottomleft", pch = 16, col = "gold", legend = "sketching subsample")
dev.off()

#### summarize Gini of granular cell types: ---------------
n.all = table(clust) 
n.sub = table(clust[res])
n.sub = n.sub[names(n.all)]
n.sub[is.na(n.sub)] = 0
names(n.sub) = names(n.all)
prop.all = n.all / sum(n.all)
prop.sub = n.sub / sum(n.sub)

DescTools::Gini(n.all)
DescTools::Gini(n.sub, na.rm = T)


lci.all = uci.all = lci.sub = uci.sub = c()
for (i in 1:length(n.all)) {
  ci = binom.confint(x = n.all[i], n = sum(n.all), methods = "wilson")
  lci.all[i] = ci$lower[1]
  uci.all[i] = ci$upper[1]
  ci = binom.confint(x = n.sub[i], n = sum(n.sub), methods = "wilson")
  lci.sub[i] = ci$lower[1]
  uci.sub[i] = ci$upper[1]
}

svg("barplot of granular cell type proportions.svg", width = 11, height = 6.5)
par(mar = c(20,5,1,1))
o = order(prop.all, decreasing = T)
bp = barplot(rbind(prop.all, prop.sub)[, o], beside = T, las = 2,
             ylab = "Proportion", cex.lab = 1.25,
             col = c("grey60", "gold"))
for (i in 1:length(o)) {
  lines(rep(bp[1, i], 2), c(lci.all[o[i]], uci.all[o[i]]))
  lines(rep(bp[2, i], 2), c(lci.sub[o[i]], uci.sub[o[i]]))
}
dev.off()



#### summarize Gini of coarse cell types: ---------------

# define coarse cell types:
clustmap = list("proximal tubule" = unique(clust[grepl("roxim", clust)]),
                "lymphoid" = c("NK.cell", "B.cell", "Plasmacytoid.dendritic.cell",
                            "", unique(clust[grepl("T.cell", clust)])),
                "myeloid" = c(unique(clust[grepl("MNP", clust)]),
                          "Neutrophil"),
                "fibroblast" = "Fibroblast",
                "myofibroblast" = "Myofibroblast",
                "endothelial" = unique(clust[grepl("ndoth", clust)]),
                "epithelial" = unique(clust[grepl("pithel", clust)]),
                "podocyte" = "Podocyte",
                "intercalated" =  unique(clust[grepl("intercal", clust)]),
                "other tubule" = c("Thick.ascending.limb.of.Loop.of.Henle", "Connecting.tubule", "Principal.cell"))

clustbig = clust
for (cl in names(clustmap)) {
  clustbig[is.element(clust, clustmap[[cl]])] = cl
}

# summarize Gini of coarse cell types:
n.all = table(clustbig)
n.sub = table(clustbig[res])
prop.all = n.all / sum(n.all)
prop.sub = n.sub / sum(n.sub)

DescTools::Gini(n.all)
DescTools::Gini(n.sub, na.rm = T)


lci.all = uci.all = lci.sub = uci.sub = c()
for (i in 1:length(n.all)) {
  ci = binom.confint(x = n.all[i], n = sum(n.all), methods = "wilson")
  lci.all[i] = ci$lower[1]
  uci.all[i] = ci$upper[1]
  ci = binom.confint(x = n.sub[i], n = sum(n.sub), methods = "wilson")
  lci.sub[i] = ci$lower[1]
  uci.sub[i] = ci$upper[1]
}

svg("barplot of coarse cell type proportions.svg", width = 11, height = 5)
par(mar = c(10,5,1,1))
o = order(prop.all, decreasing = T)
bp = barplot(rbind(prop.all, prop.sub)[, o], beside = T, las = 2,
             ylab = "Proportion", cex.lab = 1.25,
             col = c("grey60", "gold"))
for (i in 1:length(o)) {
  lines(rep(bp[1, i], 2), c(lci.all[o[i]], uci.all[o[i]]))
  lines(rep(bp[2, i], 2), c(lci.sub[o[i]], uci.sub[o[i]]))
}
dev.off()




















tab.all = table(clustbig) / nrow(counts)
tab.sub = table(clustbig[res])
tab.sub = tab.sub[match(names(tab.all), names(tab.sub))] / length(res)
names(tab.sub) = names(tab.all)
tab.sub[is.na(tab.sub)] = 0




bp = barplot(rbind(tab.all, tab.sub)[, order(tab.all)], beside = T, las = 2)
for ()



svg("panel a - abundance in fine cell types.svg")
nudge = tab.all * 0
nudge["myofibroblast"] = -0.02
nudge["endothelial"] = 0.01
nudge["other tubule"] = 0.01
plot(as.vector(tab.all), as.vector(tab.sub), col = 1, pch = 16, asp = 1,
     xlim = range(c(tab.all, tab.sub)), ylim = range(c(tab.all, tab.sub)),
     xlab = "Proportion of all cells", ylab = "Proportion of subsample")
text(as.vector(tab.all), as.vector(tab.sub) + nudge, names(tab.all), cex = 1)
abline(0,1)
dev.off()




tab.all = table(clustbig)
tab.sub = table(clustbig[res])
tab.sub = tab.sub[match(names(tab.all), names(tab.sub))]
names(tab.sub) = names(tab.all)
tab.sub[is.na(tab.sub)] = 0

binom.confint(x = tab.all["lymph"], n = sum(tab.all), methods = "wilson")[-1] * 100
binom.confint(x = tab.sub["lymph"], n = sum(tab.sub), methods = "wilson")[-1] * 100

