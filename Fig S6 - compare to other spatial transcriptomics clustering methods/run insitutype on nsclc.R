if (FALSE) {
  devtools::install_github("https://github.com/Nanostring-Biostats/InSituType")
}

rm(list = ls())
library(InSituType)
library(tictoc)
source("getcellneighbors.R")

#### load data -------------------------------

annot <- readRDS("../data/NSCLC/cell annotations.rds")
annot$negmean <- readRDS("../data/NSCLC/negmean.rds")
xy <- readRDS("../data/NSCLC/cell positions.rds")
raw <- readRDS("../data/NSCLC/raw_counts.rds")

# just take 1 slide:
keep = annot$tissue == "Lung5-5"
raw = raw[, keep]
annot = annot[keep, ]
xy = xy[keep, ]

#### cohort based on immunofluorescence and space ------------------------


## extract cells' image data:
set.seed(0)
ifdata = as.matrix(annot[, c("Area","AspectRatio","Mean.CD298","Mean.G","Mean.R","Mean.Y","Mean.DAPI")])

## extract spatial contexts:

# get neighbor cell IDs:
nmat = getcellneighbors(x = xy[, 1], y = xy[, 2], cell = annot$cell_ID, 
                        tissue = annot$tissue, n.neighbors = 50)
saveRDS(nmat, file = "results/nmat.RDS")
# get mean neighborhood expression (use original SCH code):
neighborhoodmean = as.matrix(t(raw) * NA)
for (id in rownames(neighborhoodmean)) {
  neighborhoodmean[id, ] = Matrix::rowMeans(raw[, nmat[id, ]])
}
saveRDS(neighborhoodmean, file = "results/neighborhoodmean.RDS")
# reduce to 10 PCs
neighborhoodPCs = prcomp(neighborhoodmean)$x[,1:10]
saveRDS(neighborhoodPCs, file = "results/neighborhoodPCs.RDS")


cohort = fastCohorting(mat = cbind(ifdata, neighborhoodPCs))
save(cohort, file = "results/cohort.RData")


#### run insitutype fully unsupervised --------------------------

tic()
res = InSituType::insitutype(counts = Matrix::t(raw), 
                 neg = annot$negmean,
                 cohort = cohort,
                 n_clusts = 12)
toc()
save(res, file = "results/insitutype results.RData")  #479.675 sec elapsed

