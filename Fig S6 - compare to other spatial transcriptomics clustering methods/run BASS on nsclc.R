rm(list = ls())

if (FALSE) {
  if(!require(devtools))
  install.packages(devtools)
  devtools::install_github("zhengli09/BASS")
}

library(BASS)
library(tictoc)

#### load data -------------------------------

annot <- readRDS("../data/NSCLC/cell annotations.rds")
annot$negmean <- readRDS("../data/NSCLC/negmean.rds")
xy <- readRDS("../data/NSCLC/cell positions.rds")
raw <- readRDS("../data/NSCLC/raw_counts.rds")
um <- readRDS("../data/NSCLC/umap projection.rds")

# just take 1 slide:
keep = annot$tissue == "Lung5-5"
raw = raw[, keep]
annot = annot[keep, ]
xy = xy[keep, ]
um = um[keep, ]

# create BASS object:
BASS = createBASSObject(X = list(as.matrix(raw)),
                        xy = list(xy),
                        C = 12,
                        R = 6,
                        beta_method = "fix",
                        beta = 1)
rm(raw)
tic()
# preprocess:
BASS <- BASS.preprocess(BASS)
# cluster:
BASS <- BASS.run(BASS)
BASS <- BASS.postprocess(BASS,
                         adjustLS = FALSE)
toc() # 13108.514 sec elapsed
str(BASS@results)


# extract clustering results:
bassclusts = BASS@results$c[[1]]
save(bassclusts, file = "results/BASS cluster assignments.RData")

