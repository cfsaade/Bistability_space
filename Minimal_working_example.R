## Minimal_working_example.R


## sources importation ############################################################
rm(list = ls())

library(igraph)
library(bipartite)                
library(OCNet)
library(deSolve)

source("Dynamics_functions.R")
source("Landscape_design.R")


## Generating landscapes and plotting them #######################################
## we generate each type of landscape from the main text, with 100 patches each

# Igraph options for a nicer plot
igraph_options( edge.arrow.size = 0, vertex.label = NA, vertex.size = 5)


# generating an plotting a regular graph
regular_graph = disp_matrix_grid_4(10, F)
plot(graph_from_adjacency_matrix(regular_graph>0)) #note that igraph is not well suited to plot a grid on a torus...


# generating an erdos_renyi graph
ER_graph = erdos_renyi_disp_matrix(100, 0.04)
plot(graph_from_adjacency_matrix(ER_graph>0))


# generating an RGG
RGG = RGG_disp_matrix(100, 4)
plot(graph_from_adjacency_matrix(RGG>0))


# generating an OCN (this may take a while)
OCN = OCN_disp_matrix(100)
plot(graph_from_adjacency_matrix(OCN>0), layout = layout_as_tree)

## Running the metapopulation model ###############################################
## let us run the metapopulation model in the RGG context, from random initial conditions

# declaring demographic parameters that give rise to bistability:
par = list(r=2, K = 50, B = 56, A = 25)


# taking rando initial biomasses:
n = runif(100, 0, 30)


# integrating "metapop_dynamics" over 300 time steps
out = ode(n,
          1:100,
          metapop_dynamics,
          model = NoyMeir,
          par = par,  # some parameters that give bistability, play around with this part
          disp_matrix = RGG,
          mu = 0.1) ### dispersal rate, play around with this part too


# taking a look at the dynamics:
plot(x = c(), y = c(), xlim = c(0, 100), ylim = c(0, 30), xlab = "time", ylab = "local density")
for (k in 2:101){
  points(out[,1], out[,k], type = "l")
}
