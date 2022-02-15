# Noé J Nava
# Research Agricultural Economist
# USDA Economic Research Service
# noe.nava@usda.gov
# February 14, 2022
# https://noejn2.github.io/

rm(list = ls())
options(scien = 9999)
library(tidyverse)

# Script to approximate elasticities for QUAIDS demand system
# following the work of Dong, Gould, Kaiser (2004)
# With four demos.
# The script relies on a simulaiton technique.

# For reference, rely on main paper (Nava and Dong, 2022)
# Section 3.3 Evaluation of predicted shares and demand elasticities

# Set-up -----
# Estimates from GAUSS estimations are storage in params
load("output/TabC1_params.RData")

# Mapping latent shares to observed shares - function
# This is equation (22) in Nava and Dong (2022)
obs.frm.lat.shr <- function(lat) {
  
  ncols <- length(lat)
  n <- length(lat[,1])
  
  for(c in 1:ncols) {
    
    lat[,c] <- ifelse(lat[,c] < 0, 0, lat[,c]) 
    
  }
  
  denominator <- rowSums(lat)
  
  obs <- as.data.frame(matrix(nrow = n))
  for(c in 1:ncols) {
    
    obsc <- lat[c]/denominator
    
    obsc <- as.data.frame(obsc)
    
    obs[,c] <- obsc[1]
    
    names(obs)[c] <- paste0("s", c, sep = "")
    
  }
  
  return(obs)
  
}

# Name of the parameters to loop over them
param <- c("alpha1",
           "alpha2",
           "alpha3",
           "lamda1",
           "lamda2",
           "lamda3",
           "theta11",
           "theta21",
           "theta31",
           "theta41",
           "theta12",
           "theta22",
           "theta32",
           "theta42",
           "theta13",
           "theta23",
           "theta33",
           "theta43",
           "beta1",
           "beta2",
           "beta3",
           "gamma11",
           "gamma21",
           "gamma22",
           "gamma31",
           "gamma32",
           "gamma33")

# Name of the variables to loop over them
varmu <- c("mulnp1", # ssb
           "mulnp2", # juice
           "mulnp3", # milk
           "mulnp4", # other
           "mulnw",  # expenditure 
           "muage",  # household age
           "musiz",  # household size
           "muedu",  # household head education
           "mufem")  # indicator if household is female

# For our simulation, we select seed, repetitions and a delta value
set.seed(0196) # last four digits of a 200 mxn bill (Souza, 2019)
reps <- 1000000
delta <- 1e-5

# Using main data for our analysis --- simulation is based on mean variables
ssb_data_2018 <- read_csv(file = 'data/ssb_dataset_2018.csv')

#prices --- 
mulnp1 <- mean(ssb_data_2018$lnp_ssb)
mulnp2 <- mean(ssb_data_2018$lnp_jui)
mulnp3 <- mean(ssb_data_2018$lnp_mil)
mulnp4 <- mean(ssb_data_2018$lnp_oth)

#expenditure
mulnw <- mean(ssb_data_2018$lnw, na.rm = T)

#demos
muage <- mean(log(ssb_data_2018$age))
musiz <- mean(1/ssb_data_2018$size)
muedu <- mean(ssb_data_2018$educ)
fem <- ifelse(ssb_data_2018$sex == 2, 1, 0)
mufem <- mean(fem)

# Define E(X,b) as the expected share:

# E(X,b) -- expected share without disturbance ----
### constant part
lnpindex <- 0 + alpha1*mulnp1 + alpha2*mulnp2 + alpha3*mulnp3 + alpha4*mulnp4
for(i in as.character(1:4)) {
  for(j in as.character(1:4)) {
    
    gammaij <- get(paste0("gamma", i, j))
    mulnpi  <- get(paste0("mulnp", i))
    mulnpj  <- get(paste0("mulnp", j))
    lnpindex <- lnpindex + 0.5*gammaij*mulnpi*mulnpj
    
  }
}

bofp <- 0
for(i in as.character(1:4)) {
  
  mulnpi <- get(paste0("mulnp", i))
  betai  <- get(paste0("beta", i))
  bofp <- bofp + mulnpi*betai
  
}
bofp <- exp(bofp)

s1 <- alpha1 + gamma11*mulnp1 + gamma12*mulnp2 + gamma13*mulnp3 + gamma14*mulnp4 + theta11*muage + theta21*musiz + theta31*muedu + theta41*mufem + beta1*(mulnw - lnpindex - (theta11*muage + theta21*musiz + theta31*muedu + theta41*mufem)) + (lamda1/bofp)*(mulnw - lnpindex - (theta11*muage + theta21*musiz + theta31*muedu + theta41*mufem))^2

s2 <- alpha2 + gamma21*mulnp1 + gamma22*mulnp2 + gamma23*mulnp3 + gamma24*mulnp4 + theta12*muage + theta22*musiz + theta32*muedu + theta42*mufem + beta2*(mulnw - lnpindex - (theta12*muage + theta22*musiz + theta32*muedu + theta42*mufem)) + (lamda2/bofp)*(mulnw - lnpindex - (theta12*muage + theta22*musiz + theta32*muedu + theta42*mufem))^2

s3 <- alpha3 + gamma31*mulnp1 + gamma32*mulnp2 + gamma33*mulnp3 + gamma34*mulnp4 + theta13*muage + theta23*musiz + theta33*muedu + theta43*mufem + beta3*(mulnw - lnpindex - (theta13*muage + theta23*musiz + theta33*muedu + theta43*mufem)) + (lamda3/bofp)*(mulnw - lnpindex - (theta13*muage + theta23*musiz + theta33*muedu + theta43*mufem))^2

s4 <- alpha4 + gamma41*mulnp1 + gamma42*mulnp2 + gamma43*mulnp3 + gamma44*mulnp4 + theta14*muage + theta24*musiz + theta34*muedu + theta44*mufem + beta4*(mulnw - lnpindex -(theta14*muage + theta24*musiz + theta34*muedu + theta44*mufem)) + (lamda4/bofp)*(mulnw - lnpindex - (theta14*muage + theta24*musiz + theta34*muedu + theta44*mufem))^2

# variable part
# Notice here that we rely on equation (12) in our paper
sigma <- fullsigma
sim_df <- mvtnorm::rmvnorm(n = reps, 
                           mean = rep(0,3), 
                           sigma = sigma, 
                           method = "chol")

sim_df <- as.data.frame(sim_df)
names(sim_df) <- c("e1", "e2", "e3")
# recovering fourth error term by mvtnorm assumption 
sim_df$e4 <- -(sim_df$e1 + sim_df$e2 + sim_df$e3)

# These are equation (21) in our paper 
sim_df$s1_lat <- s1 + sim_df$e1
sim_df$s2_lat <- s2 + sim_df$e2
sim_df$s3_lat <- s3 + sim_df$e3
sim_df$s4_lat <- s4 + sim_df$e4

obs_shares <- obs.frm.lat.shr(sim_df[5:8])
# Unaltered Expected shares (These are equation (23) in our paper)
E_s1 <- mean(obs_shares$s1)
E_s2 <- mean(obs_shares$s2)
E_s3 <- mean(obs_shares$s3)
E_s4 <- mean(obs_shares$s4)

# These are the repositories for our etas
# ETA repos
etas <- matrix(0, nrow = 4, ncol = length(varmu))
etas[1,1] <- -1
etas[2,2] <- -1
etas[3,3] <- -1
etas[4,4] <- -1
etas[,5]  <- 1

# ETA derivative wrt all of the parameters repos
etas1 <- matrix(0,nrow = length(varmu), ncol = length(param))
etas1[1,] <- rep(-1, length(param))
etas1[5,] <- rep(1, length(param))

etas2 <- matrix(0, nrow = length(varmu), ncol = length(param))
etas2[2,] <- rep(-1, length(param))
etas2[5,] <- rep(1, length(param))

etas3 <- matrix(0, nrow = length(varmu), ncol = length(param))
etas3[3,] <- rep(-1, length(param))
etas3[5,] <- rep(1, length(param))

etas4 <- matrix(0, nrow = length(varmu), ncol = length(param))
etas4[4,] <- rep(-1, length(param))
etas4[5,] <- rep(1, length(param))

# Elasticity estimation and S.E. -----
for (p in 1:length(varmu)) {

    # E(X + delta, b) ----
    # Disturbance is on the variable (adding delta)
    if(p %in% 1:6) {
      tru_varmu <- get(varmu[p]) # disturbed variable
      assign(varmu[p], log(exp(get(varmu[p])) + delta))
      }
  
    if(p %in% 7) {
      tru_varmu <- get(varmu[p]) # disturbed variable
      assign(varmu[p], ((get(varmu[p]))^-1 + delta)^-1)
      }
  
    if(p %in% 8:9) {
      tru_varmu <- get(varmu[p]) # disturbed variable
      assign(varmu[p], get(varmu[p])+ delta)
      }
  
    lnpindex <- 0 + alpha1*mulnp1 + alpha2*mulnp2 + alpha3*mulnp3 + alpha4*mulnp4
    for(i in as.character(1:4)) {
      for(j in as.character(1:4)) {
        
        gammaij <- get(paste0("gamma", i, j))
        mulnpi  <- get(paste0("mulnp", i))
        mulnpj  <- get(paste0("mulnp", j))
        lnpindex <- lnpindex + 0.5*gammaij*mulnpi*mulnpj
        
      }
    }
    
    bofp <- 0
    for(i in as.character(1:4)) {
      
      mulnpi <- get(paste0("mulnp", i))
      betai  <- get(paste0("beta", i))
      bofp <- bofp + mulnpi*betai
      
    }
    bofp <- exp(bofp)
    
    s1_sim <- alpha1 + gamma11*mulnp1 + gamma12*mulnp2 + gamma13*mulnp3 + gamma14*mulnp4 + theta11*muage + theta21*musiz + theta31*muedu + theta41*mufem + beta1*(mulnw - lnpindex - (theta11*muage + theta21*musiz + theta31*muedu + theta41*mufem)) + (lamda1/bofp)*(mulnw - lnpindex - (theta11*muage + theta21*musiz + theta31*muedu + theta41*mufem))^2 + sim_df$e1
    
    s2_sim <- alpha2 + gamma21*mulnp1 + gamma22*mulnp2 + gamma23*mulnp3 + gamma24*mulnp4 + theta12*muage + theta22*musiz + theta32*muedu + theta42*mufem + beta2*(mulnw - lnpindex - (theta12*muage + theta22*musiz + theta32*muedu + theta42*mufem)) + (lamda2/bofp)*(mulnw - lnpindex - (theta12*muage + theta22*musiz + theta32*muedu + theta42*mufem))^2 + sim_df$e2
    
    s3_sim <- alpha3 + gamma31*mulnp1 + gamma32*mulnp2 + gamma33*mulnp3 + gamma34*mulnp4 + theta13*muage + theta23*musiz + theta33*muedu + theta43*mufem + beta3*(mulnw - lnpindex - (theta13*muage + theta23*musiz + theta33*muedu + theta43*mufem)) + (lamda3/bofp)*(mulnw - lnpindex - (theta13*muage + theta23*musiz + theta33*muedu + theta43*mufem))^2 + sim_df$e3
    
    s4_sim <- alpha4 + gamma41*mulnp1 + gamma42*mulnp2 + gamma43*mulnp3 + gamma44*mulnp4 + theta14*muage + theta24*musiz + theta34*muedu + theta44*mufem + beta4*(mulnw - lnpindex -(theta14*muage + theta24*musiz + theta34*muedu + theta44*mufem)) + (lamda4/bofp)*(mulnw - lnpindex -(theta14*muage + theta24*musiz + theta34*muedu + theta44*mufem))^2 + sim_df$e4
    
    # Recalculating the expected share with the new shares
    
    sims <- cbind.data.frame(s1_sim, s2_sim, s3_sim, s4_sim)
    
    delta_obs_shares <- obs.frm.lat.shr(sims[1:4])
  
    # Altered Expected shares
    E_s1_x <- mean(delta_obs_shares$s1)
    E_s2_x <- mean(delta_obs_shares$s2)
    E_s3_x <- mean(delta_obs_shares$s3)
    E_s4_x <- mean(delta_obs_shares$s4)
    
    d_E_s1 <- E_s1_x - E_s1
    d_E_s2 <- E_s2_x - E_s2
    d_E_s3 <- E_s3_x - E_s3
    d_E_s4 <- E_s4_x - E_s4
  
    ###################################
    ####### Elasticity estimation #####
    ###################################
    
    # The following are Equation (24) and (25) For demos, you can think of either
    # (24) and (25) without the last term (the mu_j or the 1)
    # Notice that they are evaluated at the value of the paramter. That is,
    # if log(price) is transformed into price or 1/househodl_size is 
    # transformed into household size.
    if(p %in% 1:6) {
      etas[1,p] <- etas[1,p] + (d_E_s1/delta)*(exp(get(varmu[p])) + .5*delta)/(E_s1 + .5*d_E_s1)
      etas[2,p] <- etas[2,p] + (d_E_s2/delta)*(exp(get(varmu[p])) + .5*delta)/(E_s2 + .5*d_E_s2)
      etas[3,p] <- etas[3,p] + (d_E_s3/delta)*(exp(get(varmu[p])) + .5*delta)/(E_s3 + .5*d_E_s3)
      etas[4,p] <- etas[4,p] + (d_E_s4/delta)*(exp(get(varmu[p])) + .5*delta)/(E_s4 + .5*d_E_s4)
      }
    
    if(p %in% 7) {
      etas[1,p] <- etas[1,p] + (d_E_s1/delta)*(get(varmu[p])^-1 + .5*delta)/(E_s1 + .5*d_E_s1)
      etas[2,p] <- etas[2,p] + (d_E_s2/delta)*(get(varmu[p])^-1 + .5*delta)/(E_s2 + .5*d_E_s2)
      etas[3,p] <- etas[3,p] + (d_E_s3/delta)*(get(varmu[p])^-1 + .5*delta)/(E_s3 + .5*d_E_s3)
      etas[4,p] <- etas[4,p] + (d_E_s4/delta)*(get(varmu[p])^-1 + .5*delta)/(E_s4 + .5*d_E_s4)
      }
    
    if(p %in% 8:9) {
      etas[1,p] <- etas[1,p] + (d_E_s1/delta)*(get(varmu[p]) + .5*delta)/(E_s1 + .5*d_E_s1)
      etas[2,p] <- etas[2,p] + (d_E_s2/delta)*(get(varmu[p]) + .5*delta)/(E_s2 + .5*d_E_s2)
      etas[3,p] <- etas[3,p] + (d_E_s3/delta)*(get(varmu[p]) + .5*delta)/(E_s3 + .5*d_E_s3)
      etas[4,p] <- etas[4,p] + (d_E_s4/delta)*(get(varmu[p]) + .5*delta)/(E_s4 + .5*d_E_s4)
    }
    
    # Re-assingning initial value to the variable
    assign(varmu[p], tru_varmu) # baseline variable
  
    # Disturbance is on the parameter (adding delta)
    for(k in 1:length(param)) {
      # E(X + delta, b + delta) ----
      if(p %in% 1:6) {
        tru_varmu <- get(varmu[p]) # disturbed variable
        assign(varmu[p], log(exp(get(varmu[p])) + delta))
      }
      
      if(p %in% 7) {
        tru_varmu <- get(varmu[p]) # disturbed variable
        assign(varmu[p], ((get(varmu[p]))^-1 + delta)^-1)
      }
      
      if(p %in% 8:9) {
        tru_varmu <- get(varmu[p]) # disturbed variable
        assign(varmu[p], get(varmu[p])+ delta)
      }
      
      # Making the change in coefficient only 
      tru_param <- get(param[k]) # parameter
      assign(param[k], tru_param + delta)
      
      lnpindex <- 0 + alpha1*mulnp1 + alpha2*mulnp2 + alpha3*mulnp3 + alpha4*mulnp4
      for(i in as.character(1:4)) {
        for(j in as.character(1:4)) {
          
          gammaij <- get(paste0("gamma", i, j))
          mulnpi  <- get(paste0("mulnp", i))
          mulnpj  <- get(paste0("mulnp", j))
          lnpindex <- lnpindex + 0.5*gammaij*mulnpi*mulnpj
          
        }
      }
      
      bofp <- 0
      for(i in as.character(1:4)) {
        
        mulnpi <- get(paste0("mulnp", i))
        betai  <- get(paste0("beta", i))
        bofp <- bofp + mulnpi*betai
        
      }
      bofp <- exp(bofp)
      
      s1_sim <- alpha1 + gamma11*mulnp1 + gamma12*mulnp2 + gamma13*mulnp3 + gamma14*mulnp4 + theta11*muage + theta21*musiz + theta31*muedu + theta41*mufem + beta1*(mulnw - lnpindex - (theta11*muage + theta21*musiz + theta31*muedu + theta41*mufem)) + (lamda1/bofp)*(mulnw - lnpindex - (theta11*muage + theta21*musiz + theta31*muedu + theta41*mufem))^2 + sim_df$e1
      
      s2_sim <- alpha2 + gamma21*mulnp1 + gamma22*mulnp2 + gamma23*mulnp3 + gamma24*mulnp4 + theta12*muage + theta22*musiz + theta32*muedu + theta42*mufem + beta2*(mulnw - lnpindex - (theta12*muage + theta22*musiz + theta32*muedu + theta42*mufem)) + (lamda2/bofp)*(mulnw - lnpindex - (theta12*muage + theta22*musiz + theta32*muedu + theta42*mufem))^2 + sim_df$e2
      
      s3_sim <- alpha3 + gamma31*mulnp1 + gamma32*mulnp2 + gamma33*mulnp3 + gamma34*mulnp4 + theta13*muage + theta23*musiz + theta33*muedu + theta43*mufem + beta3*(mulnw - lnpindex - (theta13*muage + theta23*musiz + theta33*muedu + theta43*mufem)) + (lamda3/bofp)*(mulnw - lnpindex - (theta13*muage + theta23*musiz + theta33*muedu + theta43*mufem))^2 + sim_df$e3
      
      s4_sim <- alpha4 + gamma41*mulnp1 + gamma42*mulnp2 + gamma43*mulnp3 + gamma44*mulnp4 + theta14*muage + theta24*musiz + theta34*muedu + theta44*mufem + beta4*(mulnw - lnpindex -(theta14*muage + theta24*musiz + theta34*muedu + theta44*mufem)) + (lamda4/bofp)*(mulnw - lnpindex -(theta14*muage + theta24*musiz + theta34*muedu + theta44*mufem))^2 + sim_df$e4     
      # Recalculating the expected share with the new shares
      
      sims <- cbind.data.frame(s1_sim, s2_sim, s3_sim, s4_sim)
      
      delta_obs_shares <- obs.frm.lat.shr(sims[1:4])
      
      # Altered Expected shares
      E_s1_xb <- mean(delta_obs_shares$s1)
      E_s2_xb <- mean(delta_obs_shares$s2)
      E_s3_xb <- mean(delta_obs_shares$s3)
      E_s4_xb <- mean(delta_obs_shares$s4)
      
      assign(varmu[p], tru_varmu) # baseline variable
      
      # E(X, b + delta) ----
  
      lnpindex <- 0 + alpha1*mulnp1 + alpha2*mulnp2 + alpha3*mulnp3 + alpha4*mulnp4
      for(i in as.character(1:4)) {
        for(j in as.character(1:4)) {
          
          gammaij <- get(paste0("gamma", i, j))
          mulnpi  <- get(paste0("mulnp", i))
          mulnpj  <- get(paste0("mulnp", j))
          lnpindex <- lnpindex + 0.5*gammaij*mulnpi*mulnpj
          
        }
      }
      
      bofp <- 0
      for(i in as.character(1:4)) {
        
        mulnpi <- get(paste0("mulnp", i))
        betai  <- get(paste0("beta", i))
        bofp <- bofp + mulnpi*betai
        
      }
      bofp <- exp(bofp)
      
      s1_sim <- alpha1 + gamma11*mulnp1 + gamma12*mulnp2 + gamma13*mulnp3 + gamma14*mulnp4 + theta11*muage + theta21*musiz + theta31*muedu + theta41*mufem + beta1*(mulnw - lnpindex - (theta11*muage + theta21*musiz + theta31*muedu + theta41*mufem)) + (lamda1/bofp)*(mulnw - lnpindex - (theta11*muage + theta21*musiz + theta31*muedu + theta41*mufem))^2 + sim_df$e1
      
      s2_sim <- alpha2 + gamma21*mulnp1 + gamma22*mulnp2 + gamma23*mulnp3 + gamma24*mulnp4 + theta12*muage + theta22*musiz + theta32*muedu + theta42*mufem + beta2*(mulnw - lnpindex - (theta12*muage + theta22*musiz + theta32*muedu + theta42*mufem)) + (lamda2/bofp)*(mulnw - lnpindex - (theta12*muage + theta22*musiz + theta32*muedu + theta42*mufem))^2 + sim_df$e2
      
      s3_sim <- alpha3 + gamma31*mulnp1 + gamma32*mulnp2 + gamma33*mulnp3 + gamma34*mulnp4 + theta13*muage + theta23*musiz + theta33*muedu + theta43*mufem + beta3*(mulnw - lnpindex - (theta13*muage + theta23*musiz + theta33*muedu + theta43*mufem)) + (lamda3/bofp)*(mulnw - lnpindex - (theta13*muage + theta23*musiz + theta33*muedu + theta43*mufem))^2 + sim_df$e3
      
      s4_sim <- alpha4 + gamma41*mulnp1 + gamma42*mulnp2 + gamma43*mulnp3 + gamma44*mulnp4 + theta14*muage + theta24*musiz + theta34*muedu + theta44*mufem + beta4*(mulnw - lnpindex - (theta14*muage + theta24*musiz + theta34*muedu + theta44*mufem)) + (lamda4/bofp)*(mulnw - lnpindex -(theta14*muage + theta24*musiz + theta34*muedu + theta44*mufem))^2 + sim_df$e4
      
      # Recalculating the expected share with the new shares
      
      sims <- cbind.data.frame(s1_sim, s2_sim, s3_sim, s4_sim)
      
      delta_obs_shares <- obs.frm.lat.shr(sims[1:4])
      
      # Altered Expected shares
      E_s1_b <- mean(delta_obs_shares$s1)
      E_s2_b <- mean(delta_obs_shares$s2)
      E_s3_b <- mean(delta_obs_shares$s3)
      E_s4_b <- mean(delta_obs_shares$s4)
      
      d_E_s1_SE <- E_s1_xb - E_s1_b
      d_E_s2_SE <- E_s2_xb - E_s2_b
      d_E_s3_SE <- E_s3_xb - E_s3_b
      d_E_s4_SE <- E_s4_xb - E_s4_b
      
      ###################################
      ####### Elasticity estimation ##### --- S.E.
      ###################################
      
      # Here it is the new eta with the disturbed parameter value
      if(p %in% 1:6) {
        etas1[p,k] <- etas1[p,k] + (d_E_s1_SE/delta)*(exp(get(varmu[p])) + .5*delta)/(E_s1_b + .5*d_E_s1_SE)
        etas2[p,k] <- etas2[p,k] + (d_E_s2_SE/delta)*(exp(get(varmu[p])) + .5*delta)/(E_s2_b + .5*d_E_s2_SE)
        etas3[p,k] <- etas3[p,k] + (d_E_s3_SE/delta)*(exp(get(varmu[p])) + .5*delta)/(E_s3_b + .5*d_E_s3_SE)
        etas4[p,k] <- etas4[p,k] + (d_E_s4_SE/delta)*(exp(get(varmu[p])) + .5*delta)/(E_s4_b + .5*d_E_s4_SE)
      }
      
      if(p %in% 7) {
        etas1[p,k] <- etas1[p,k] + (d_E_s1_SE/delta)*(get(varmu[p])^-1 + .5*delta)/(E_s1_b + .5*d_E_s1_SE)
        etas2[p,k] <- etas2[p,k] + (d_E_s2_SE/delta)*(get(varmu[p])^-1 + .5*delta)/(E_s2_b + .5*d_E_s2_SE)
        etas3[p,k] <- etas3[p,k] + (d_E_s3_SE/delta)*(get(varmu[p])^-1 + .5*delta)/(E_s3_b + .5*d_E_s3_SE)
        etas4[p,k] <- etas4[p,k] + (d_E_s4_SE/delta)*(get(varmu[p])^-1 + .5*delta)/(E_s4_b + .5*d_E_s4_SE)
      }
      
      if(p %in% 8:9) {
        etas1[p,k] <- etas1[p,k] + (d_E_s1_SE/delta)*(get(varmu[p]) + .5*delta)/(E_s1_b + .5*d_E_s1_SE)
        etas2[p,k] <- etas2[p,k] + (d_E_s2_SE/delta)*(get(varmu[p]) + .5*delta)/(E_s2_b + .5*d_E_s2_SE)
        etas3[p,k] <- etas3[p,k] + (d_E_s3_SE/delta)*(get(varmu[p]) + .5*delta)/(E_s3_b + .5*d_E_s3_SE)
        etas4[p,k] <- etas4[p,k] + (d_E_s4_SE/delta)*(get(varmu[p]) + .5*delta)/(E_s4_b + .5*d_E_s4_SE)
      }
      
      # Here, we approximate the derivate of the elasticity wrt the parameter
      etas1[p,k] <- ((etas1[p,k] - etas[1,p])/delta)
      etas2[p,k] <- ((etas2[p,k] - etas[2,p])/delta)
      etas3[p,k] <- ((etas3[p,k] - etas[3,p])/delta)
      etas4[p,k] <- ((etas4[p,k] - etas[4,p])/delta)
      
      # Re-assigning the initial value of the parameter
      assign(param[k], tru_param) #baseline parameter
      
    }
}

# VC4 comes from the GAUSS and it is the variance-covariance matrix
cov <- readxl::read_xlsx(path = 'output/TabC1_vcov.xlsx', col_names = FALSE)
cov <- as.matrix(cov[1:27,1:27])

etas1a <- as.matrix(etas1)
etas2 <- as.matrix(etas2)
etas3 <- as.matrix(etas3)
etas4 <- as.matrix(etas4)

# We use the variances obtained from the delta method to recover S.E.s
SE_wrt1 <- sqrt(diag(etas1 %*% cov %*% t(etas1)))#/sqrt(length(ssb_data_2018$folioviv))
SE_wrt2 <- sqrt(diag(etas2 %*% cov %*% t(etas2)))#/sqrt(length(ssb_data_2018$folioviv))
SE_wrt3 <- sqrt(diag(etas3 %*% cov %*% t(etas3)))#/sqrt(length(ssb_data_2018$folioviv))
SE_wrt4 <- sqrt(diag(etas4 %*% cov %*% t(etas4)))#/sqrt(length(ssb_data_2018$folioviv))

############################
##### Printing Results #####
############################
output <- matrix(NA, nrow = 2*4, ncol = 9)
colnames(output) <- varmu
vars <- c("SSB", "s.e.", "Juice", "s.e.", "Milk", "s.e.", "Other", "s.e.")
rownames(output) <- vars
for(o in 1:4) {
  place <- seq(from = 1, to = 7, by = 2)[o]
  output[place,] <- etas[o,]
  SE_name <- paste0("SE_wrt", o)
  output[(place + 1),] <- get(SE_name)
}
output

saveRDS(output, file = 'output/etas.rds')
output <- as.data.frame(output)
output <- cbind.data.frame(vars, output)
readr::write_csv(output, file = 'output/etas.csv')

#############################
####### Hicksian Etas #######
#############################

# The hicksian elasticities recovered using the slutsky equation
# We also recovered their variances by 
# Var(eta_ij + eta_i*S_j) = var(eta_ij) + var(eta_i)*S_j^2
# we assume that cov(eta_i, S_j) = 0

etas_m <- etas[,1:4]
etas_b <- etas[,5]
E_shr <- c(E_s1, E_s2, E_s3, E_s4)

etas_b_mat <- matrix(NA, nrow = 4, ncol = 4)
etas_b_mat[,1] <- etas_b
etas_b_mat[,2] <- etas_b
etas_b_mat[,3] <- etas_b
etas_b_mat[,4] <- etas_b

E_shr_mat <- matrix(NA, nrow = 4, ncol = 4)
E_shr_mat[1,] <- E_shr
E_shr_mat[2,] <- E_shr
E_shr_mat[3,] <- E_shr
E_shr_mat[4,] <- E_shr

second_mat <- matrixcalc::hadamard.prod(etas_b_mat, E_shr_mat)

etas_h <- etas_m + second_mat

SE <- rbind(SE_wrt1*sqrt(length(ssb_data_2018$folioviv)), 
            SE_wrt2*sqrt(length(ssb_data_2018$folioviv)),
            SE_wrt3*sqrt(length(ssb_data_2018$folioviv)),
            SE_wrt4*sqrt(length(ssb_data_2018$folioviv)))

etas_m_SE <- SE[,1:4]
etas_b_SE <- SE[,5]

etas_b_SE_mat <- matrix(NA, nrow = 4, ncol = 4)
etas_b_SE_mat[,1] <- etas_b_SE
etas_b_SE_mat[,2] <- etas_b_SE
etas_b_SE_mat[,3] <- etas_b_SE
etas_b_SE_mat[,4] <- etas_b_SE

second_mat_SE <- matrixcalc::hadamard.prod(etas_b_SE_mat, E_shr_mat^2)

etas_h_SE <- (etas_m_SE + second_mat_SE)/sqrt(length(ssb_data_2018$folioviv))

### Printing results for hicksian etas
output_hic <- matrix(NA, nrow = 2*4, ncol = 4)
colnames(output_hic) <- varmu[1:4]
vars <- c("SSB", "s.e.", "Juice", "s.e.", "Milk", "s.e.", "Other", "s.e.")
rownames(output_hic) <- vars
for(o in 1:4) {
  place <- seq(from = 1, to = 7, by = 2)[o]
  output_hic[place,] <- etas_h[o,]
  output_hic[(place + 1),] <- etas_h_SE[o,]
  
}
output_hic

saveRDS(output_hic, file = 'output/etas_hic.rds')
output_hic <- as.data.frame(output_hic)
output_hic <- cbind.data.frame(output_hic, output)
readr::write_csv(output_hic, file = 'output/etas_hic.csv')
