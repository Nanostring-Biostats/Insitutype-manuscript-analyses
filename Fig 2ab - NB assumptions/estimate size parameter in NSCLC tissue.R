
library(MLEcell)
library(RColorBrewer)
library(pheatmap)


#### load the nsclc data: ----------------------------


# estimate background:

#### clustering: ------------------------------

cohorts = fastCohorting(mat = as.matrix(annot[, c("Area", "AspectRatio", "Mean.CD298", "Mean.PanCK", "Mean.CD45", "Mean.CD3", "Mean.DAPI")]), 
                        gaussian_transform = T)
str(cohorts)

# run insitutype
out = insitutype(counts = counts, 
                 neg = neg,
                 n_clusts = 20,
                 size = 10,
                 cohort = cohorts)
save(out, cohorts, file = "supervised cell typing results.RData")


#### estimate yhat: -----------------------------

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
for (celltype in unique(out$clust)) {
  print(celltype)
  use = out$clust == celltype
  temp = getyhat(mat = counts[use,], x = out$profiles[,celltype], bg = annot$bg[use])  
  yhat[use, ] = temp 
}

