## Landscape_design.R
#
# Contains functions to generate dispersal matrix for various landscapes structures (circle, grid, Erdos-Renyi, RGG and OCN)
# for graphs that are not regular, the functions xxx_adj(...) returns an adjacency matrix (indices are 0 or 1) while xxx_disp_matrix(...) returns a dispersal matrix.

# libraries ##############################################################################################################

library(igraph)
library(bipartite)                
library(OCNet)


## regular graphs ########################################################################################################

# functions to generate regular graphs (i.e., graphs where all patches have the same number of connections): circular networks and grids (with 4 or 8 neighbours)


  disp_matrix_circle = function(n_patch){
    # generates a circular landscape with n_patch
    a = matrix(rep(0, n_patch**2), nrow = n_patch)
    influx = 1
    if (n_patch > 2){
      influx = 0.5
    }
    a[1,n_patch] = influx
    a[n_patch,1] = influx
    for (i in 1:n_patch){
      a[i, i] = -1
      if (i != 1){
        a[i, i-1] = influx
      }
      if(i != n_patch){
        a[i, i+1] = influx
      }
    }
    return(a)
  }

  ## Declaring some utility functions useful to generate grids:
    
    position2rowcol = function(position, size = 4){
      row = ceiling(position/size)
      column = (position-1)%%size + 1
      return(c(row, column))
    }
    
    grid.distance = function(pos1, pos2, size = 4, mod="nodes",border = T){
      pos1.rowcol = position2rowcol(pos1, size = size)
      pos2.rowcol = position2rowcol(pos2, size = size)
      
      dist.vector = abs(pos1.rowcol - pos2.rowcol)
      if (border == F & dist.vector[1] == size - 1){
        dist.vector[1] = 1
      }
      if (border == F & dist.vector[2] == size - 1){
        dist.vector[2] = 1
      }
      return(sum(dist.vector))
    }
    


  # Function to generate a regular graph on a grid with 4 neighbours
    disp_matrix_grid_4 = function(size = 4, border = F){
      # Generates a square matrix of size*size patches arranged on a grid
      # the "border" parameter controls if the landscape has borders (T) or is arranged on a torus (F)
      npatches = size*size
      neigh = matrix(rep(0, npatches**2), nrow = npatches)
      for (k in 1:npatches){
        for (l in 1:npatches){
          if (grid.distance(k, l, size = size, border = border) == 1){
            neigh[k, l] = 1
            neigh[l, k] = 1
          }
        }
      }
      for (k in 1:npatches){
        neigh[,k] = neigh[,k]/sum(neigh[,k])
      }
      for (k in 1:npatches){
        neigh[k,k] = -1
      }
      return(neigh)
    }
    
  # Function to generate a regular graph on a grid with 8 neighbours    
    disp_matrix_grid_8 = function(size = 4, border = F){
      npatches = size*size
      neigh = matrix(rep(0, npatches**2), nrow = npatches)
      for (k in 1:npatches){
        for (l in 1:npatches){
          if (grid.distance(k, l, size = size, border = border) %in% c(1,2)){
            neigh[k, l] = 1
            neigh[l, k] = 1
          }
        }
      }
      save=neigh
      
      for (i in 1:nrow(neigh)){
        for (j in 1:nrow(neigh)){
          if (j==i+2 | i==j+2 | j==i-2 | i==j-2 |
              j==i+2*size | j==i-2*size | i==j+2*size | i==j-2*size){
            neigh[i,j]=0
          }
        }
      }
      
      for (k in 1:npatches){
        neigh[,k] = neigh[,k]/sum(neigh[,k])
      }
      for (k in 1:npatches){
        neigh[k,k] = -1
      }
      
      return(neigh)
    }


## Erdos-Renyi graphs #########################################################################
# functions to generate erdos renyi (random) graphs, with a probability of connection "p"
    erdos_renyi_adj = function(n_patch, p){
      # returns an adjacency matrix (0, 1)
      output = matrix(rep(0, n_patch**2), nrow = n_patch)
      for (k in 1:(n_patch-1)){
        for(l in (k+1):n_patch){
          if(runif(1, 0, 1) < p){
            output[k,l] = 1
            output[l,k] = 1
          }  
        }
      }
      return(output)
    }
    
    erdos_renyi_disp_matrix = function(n_patch, p){
      # returns a connected dispersal matrix: each column sums to 1 (e.g., outgoing individuals are equally spread between links)
      adj_matrix = erdos_renyi_adj(n_patch, p)
      while (is_connected(graph_from_adjacency_matrix(adj_matrix)) != T){
        adj_matrix = erdos_renyi_adj(n_patch, p)
      }
      for (k in 1:n_patch){
        adj_matrix[,k] = adj_matrix[,k]/sum(adj_matrix[,k])
      }
      return(adj_matrix)
    }

    
## Random geometric graphs #########################################################################################################
# functions to generate random geometric graphs with a target connectivity c
    
    
    RGG_adj = function(n_patch, c, steps = 20){
      # n_patch: number of patches
      # c: target connectivity (number of links per patch)
      # steps: number of prospected threshold distance to reach this connectivity (more steps = more precise)
      points = data.frame(x = runif(n_patch, 0, 1), y = runif(n_patch, 0, 1))
      
      dist_matrix = matrix(0, nrow = n_patch, ncol = n_patch)
      for (k in 1:n_patch){
        for(l in 1:n_patch){
          dist = sqrt( (points$x[k]-points$x[l])**2 + (points$y[k]-points$y[l])**2 )
          dist_matrix[k,l] = dist
          dist_matrix[l,k] = dist
        }
        dist_matrix[k, k] = 1.2
      }
      
      # getting the minimal distance to connect all patches at least once
      list_min = apply(dist_matrix, 2, min)
      min_threshold = max(list_min)
      
      list_threshold = seq(min_threshold, sqrt(2), length.out = steps)
      list_connectivity = rep(0, steps)
      
      for (k in 1:steps){
        dist_threshold = list_threshold[k]
        adj = dist_matrix <= dist_threshold
        for(l in 1:n_patch){
          adj[l,l] = 0
        }
        
        list_connectivity[k] = sum(adj)/(2*n_patch)
        
      }
      
      thr_index = which(abs(list_connectivity-c) == min(abs(list_connectivity-c)))[1]
      thr = list_threshold[thr_index]
      
      output = dist_matrix <= thr
      for(l in 1:n_patch){
        output[l,l] = 0
      }
      return(output)
    }
    
    
    RGG_disp_matrix = function(n_patch, c, steps = 100){
      adj_matrix = RGG_adj(n_patch, c, steps)
      while (is_connected(graph_from_adjacency_matrix(adj_matrix)) != T){
        adj_matrix = RGG_adj(n_patch, c, steps)
      }
      
      for (k in 1:(n_patch)){
        adj_matrix[,k] = adj_matrix[,k]/sum(adj_matrix[,k])
      }
      return(adj_matrix)
    }
    

## Optimal Chanel Networks (OCNs) ############################################################################################
# functions to generate OCNs

OCN_adj = function(n_patch){
  OCN <- create_OCN(20, 20, outletPos = 1, cellsize = 500)
  OCN <- landscape_OCN(OCN, slope0 = 0.01)
  
  thr <- find_area_threshold_OCN(OCN)
  
  indThr <- which(abs(thr$nNodesAG - n_patch) == min(abs(thr$nNodesAG - n_patch)))
  indThr <- max(indThr) # pick the last ind_thr that satisfies the condition above
  thrA <- thr$thrValues[indThr] # corresponding threshold area
  OCN <- aggregate_OCN(OCN, thrA = thrA)
  
  output = as.matrix(OCN$AG$W)
  
  while(nrow(output) != n_patch){
    OCN <- create_OCN(20, 20, outletPos = 1, cellsize = 500)
    OCN <- landscape_OCN(OCN, slope0 = 0.01)
    
    thr <- find_area_threshold_OCN(OCN)
    indThr <- which(abs(thr$nNodesAG - n_patch) == min(abs(thr$nNodesAG - n_patch)))
    indThr <- max(indThr) # pick the last ind_thr that satisfies the condition above
    thrA <- thr$thrValues[indThr] # corresponding threshold area
    OCN <- aggregate_OCN(OCN, thrA = thrA)
    
    output = as.matrix(OCN$AG$W)
    print(nrow(output))
  }
  
  return(output + t(output))
}

OCN_disp_matrix = function(n_patch){
  adj_matrix = OCN_adj(n_patch)
  
  for (k in 1:(n_patch)){
    adj_matrix[,k] = adj_matrix[,k]/sum(adj_matrix[,k])
  }
  diag(adj_matrix)=-1
  
  return(adj_matrix)
}

