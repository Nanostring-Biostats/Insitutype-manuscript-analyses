



#' Function to identify a cell's neighbors
#' 
#' Given cells' xy coordinates and metric of distance, identifies the closest cells.
#'  Does this for all cells in a dataset. 
#' @param x x coordinates of cells
#' @param y y coordinates of cells
#' @param cell Vector if cell IDs
#' @param n.neighbors Identify the closest n.neighbors cells
#' @param radius Identify all cells within this radius. Required if n.neighbors isn't provided. 
#' @return Two outputs:
#'  neighborlist: A list in which each element is a cell's neighbors 
#'  neighborvec: A vector in which each element is a concatenated vector of a cell's neighbors
getcellneighbors <- function(x, y, cell, tissue = NULL, n.neighbors = 10) {
  if (is.null(tissue)) {
    tissue = rep("a", length(x))
  }
  
  ## get distances between all cells:
  # for each cell, get closest k neighbors:
  getonecellneighbors = function(cid, x, y, cell, n.neighbors) {
    # calc distance to other cells
    xtemp = x[cell == cid]
    ytemp = y[cell == cid]
    tempdists = (x - xtemp)^2 + (y - ytemp)^2
    # get nearest bunch:
    nearest = as.character(cell[tempdists < quantile(tempdists, n.neighbors * 1.2 / length(tempdists))])
    nearest = setdiff(nearest, cid)
    # return them in order:
    nearest = nearest[order(tempdists[match(nearest, cell)])[seq_len(n.neighbors)]]
    return(nearest)
  }
  
  gettissueneighbors = function(tissuename) {
    
    use = tissue == tissuename
    tempx = x[use]
    tempy = y[use]
    tempcell = cell[use]
    
    tempneighbors = sapply(tempcell, getonecellneighbors, tempx, tempy, tempcell, n.neighbors)
    return(t(tempneighbors))
  }
  neighbors = rlist::list.rbind(lapply(unique(tissue), gettissueneighbors))
  return(neighbors[cell, ])
}


# for debugging and testing:
#x=  rnorm(2000)
#y = rnorm(2000)
#cell = paste0("c",1:2000)
#tissue = rep(c("a", "b"), each = 1000)
#x[tissue == "b"] = x[tissue == "b"] + 100 
#n.neighbors = 10
#neighbormat = getcellneighbors(x,y,cell,tissue,n.neighbors = 5)
