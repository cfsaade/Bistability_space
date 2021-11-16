## Dynamics_functions.R
#
# Contains functions to simulate the dynamics of a single patch ("NoyMeir") or a whole metapopulation ("metapop_dynamics")


## NoyMeir ###################################################################
##
## Returns the derivative of biomass in a single patch or a vector of local densities (n) for a given set of parameters (p)

NoyMeir = function(t, n, par){
  # n = vector of local densities plant density
  # par = vector of parameters
  
  r = par[[1]]    # growth rate
  K = par[[2]]    # carrying capacity
  B = par[[3]]    # comsumption intensity
  A = par[[4]]    # half consumption biomass
  
  dn = r*n*(1-n/K) - B*n/(A+n)
  
  return(list(dn))
}


## metapop_dynamics ###########################################################
##
## returns the dynamics (derivative at a given time) of an n_patch system
## integrate with the function "ode" from package "deSolve"

metapop_dynamics = function(t, n, model, par, disp_matrix, mu){
  # n = plant density (vector of n_patch elements)
  # model = the function for local dynamics (here always "NoyMeir" but other functions could be used)
  # parameters for the given model
  # disp matrix = a matrix giving the links between patches
  # mu = dispersal rate
  # t = dummy parameter necessary for deSolve
  
  dn = model(t, n, par)[[1]]
  
  dn = dn + disp_matrix%*%(mu*n)
  return(list(c(dn)))
}

