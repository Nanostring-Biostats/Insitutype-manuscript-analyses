# evaluate the likelihood model in a cell pellet array:

# package installation:
#library(devtools)
#install("../SMI-cell-typing")


library(MLEcell)
library(RColorBrewer)
library(pheatmap)
library(scales)

## load the CPA data:
load("../data/Cell Pellet Array annotation and raw counts.RData")
counts = t(as.matrix(raw))
rm(raw)


# estimate each cell's background tendency:
annot$bg = NA
for (celltype in unique(annot$cell_line)) {
  use = annot$cell_line == celltype
  bgmod <- stats::lm(annot$neg[use] ~ annot$raw_totalCount[use] - 1)
  annot$bg[use] <- bgmod$fitted
}

# get mean profile from each cell pellet:
profiles = MLEcell:::Estep(counts = counts,
                           clust = annot$cell_line,
                           neg = annot$negmean)


# function to get yhat matrix (adapted from MLEcell:::lldist function)
getyhat <- function(x, mat, bg = 0.01, size = 10, digits = 4) {
  # convert to matrix form if only a vector was input:
  if (is.vector(mat)) {
    mat <- matrix(mat, nrow = 1)
  }
  # calc scaling factor to put y on the scale of x:
  if ( is.vector( bg ) )
  {
    bgsub <- pmax( sweep( mat , 1 , bg , "-" ) , 0 )
  }
  else
  {
    bgsub <- pmax( mat - bg , 0 )
  }
  s <- rowSums( bgsub ) / sum( x )
  # override it if s is negative:
  s[s <= 0] = Matrix::rowSums(mat[s <= 0, , drop = FALSE]) / sum(x)
  
  # expected counts:
  if ( is.vector( bg ) )
  {
    yhat <- sweep( s %*% t( x ) , 1 , bg , "+" )
  }
  else
  {
    yhat <- s %*% t(x) + bg
  }
  return(yhat)
}

# get yhat matrix
yhat = counts * NA
for (celltype in unique(annot$cell_line)) {
  print(celltype)
  use = annot$cell_line == celltype
  temp = getyhat(mat = counts[use, ], x = profiles[, celltype], bg = annot$bg[use])  
  yhat[use, ] = temp 
}


#### estimate theta (size) for each gene * cell line ---------------------------


# function to estimate size:
estsize <- function(y, mu) {
  ll <- function(size, y, mu) {
    -sum(dnbinom(y, mu = mu, size = size, log = T))
  }
  op = optim(par = 1, fn = ll, y = y, mu = mu, method = "Brent", lower = 0.01, upper = 1000)$par
  return(op)
}


thetamat = matrix(NA, ncol(counts), ncol(profiles),
                  dimnames = list(colnames(counts), colnames(profiles)))
for (gene in colnames(counts)) {
  for (celltype in colnames(profiles)) {
    tempy = counts[annot$cell_line == celltype, gene]
    tempyhat = yhat[annot$cell_line == celltype, gene]
    thetamat[gene, celltype] = estsize(y = tempy, mu = tempyhat)
  }
}

plot(profiles + 1e-3, thetamat, log = "xy") # negligible trend
pheatmap(log2(thetamat))
# conclusions:
# some genes are poisson-like, but they're inconsistent (cell type specific)
# some genes are highly over-dispersed (low thetas) across cell types

# estimate theta:
#theta.est = round(median(thetamat), 1) # 1.5
theta.est = 10

#### goodness-of-fit plots ---------------------------------------------

yhattarget = c(0.1, 0.5, 1, 2, 5, 10)
indslist = list()
for (i in 1:length(yhattarget)) {
  #indslist[[i]] = which((yhat > 0.99 * yhattarget[i]) & (yhat <= 1.01 * yhattarget[i])
  indslist[[i]] = order(abs(yhat - yhattarget[i]))[1:2000]
}

# densities:
png("NB quantiles.png", width = 8, height = 3, units = "in", res = 600)
graphics::layout(mat = t(matrix(c(1:12, rep(13,6)), 6)), heights = c(2,2,.5), widths = c(5, rep(4,5)))
par(mar = c(2,4.5,2,0.5))
for (i in 1:length(indslist)) {

  if (i > 1) {
    par(mar = c(2,2,2,0.5))
  }
  inds = indslist[[i]]
  tempy = counts[inds]
  rand = rnbinom(20000, mu = mean(yhat[inds]), size = theta.est)
  breaks = -0.5+seq(0, max(c(rand, tempy))+1)
  hist(rand, col = alpha("cornflowerblue", 0.75), freq = F, breaks = breaks, border = NA,
       ylab = "Frequency", cex.lab = 1.3, 
       main = paste0("E(y) = ", yhattarget[i]), xlab = "",yaxt = "n")
  hist(tempy, col = alpha("orange", 0.5), freq = F, breaks = breaks, border = NA, add = T)
  if (i == 1) {
    legend("topright", fill = c(alpha("cornflowerblue", 0.75), alpha("orange", 0.5)), legend = c("Theoretical", "Observed"),
           bty = "n", cex = 0.8)
  }
}

# qq plots
par(mar = c(2,4.5,2,0.5))
for (i in 1:length(indslist)) {
  
  if (i > 1) {
    par(mar = c(2,2,2,0.5))
  }
  inds = indslist[[i]]
  tempy = counts[inds]
  rand = rnbinom(length(tempy), mu = mean(yhat[inds]), size = theta.est)
  plot(jitter(rand[order(rand)]), jitter(tempy[order(tempy)]), cex = 0.1,
       col = alpha("grey20", 0.75),
       ylab = "Observed quantiles", cex.lab = 1.3,
       main = paste0("E(y) = ", yhattarget[i]))
  abline(0,1)
}
par(mar = c(0,0,0,0))
frame()
legend("center", legend = "Theoretical quantiles", bty = "n", cex = 1.3)

dev.off()


# observation: high genes are less overdispersed than the model suggests



