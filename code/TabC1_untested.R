# Noé J Nava
# Research Agricultural Economist
# USDA Economic Research Service
# noe.nava@usda.gov
# February 14, 2022
# https://noejn2.github.io/
rm(list = ls())

# Notes:
# Consult citation in journal (JAAEA): 
# Nava and Dong (2022) The impact of taxing sugary soft-beverages: A censored 
# QUAI demand system approach. Journal of the Agricultural and Applied Economics Association.

# Script estimates the parameters of a censored QUAI demand system by parametrizing 
# a truncated log-likelihood. There are different cases as explained in the manuscript
# based on the 2^M - 1 = 15 purchasing regimes, where M = 4 is the number of goods
# in the basket of goods. These purchasing regimes can be classified into four
# cases depending on the number of zero-purchases. For instance, the case in which
# a household buys all the four goods represents a likelihood function of the 
# regular QUAIDS. For the case in which the household buys fewer than four, then
# we have three different parametrizations for each of these households.

# Therefore, the log-likelihood function, quaids.loglike(), loops over every single
# household and then calculates the log-likelihood contribution. This is not the 
# case for whenever the household buys all of the goods. This loop, however, is
# problematic since it requires a lot of computational power.

# There are some optimization options that you can play around. You can play around 
# with the step. It could be the case that in a random subset of data, the step is large 
# enough that the diagonal elements of the vcov matrix become negative. In addition, be # aware of not using zeros for the initial values since this practice will cause
# problems. As for the diagonal elements of the sigma, also do not use zeros.

# Another recommendation: Simplify the code a bit further if possible.
# Contact authors if succesfully done.

library(mvtnorm)
library(tidyverse)
library(rootSolve)
library(Matrix)
library(numDeriv)

# Optimization options:
c0         <- 0.001        # Convergence criteria
maxIter    <- 1000         # Maximum numnber of iterations to theck
step       <- 1            # Sufficiently small (usually 0.5 works fine)
reap       <- 1            # Report convergence diagnostics every `reap` times
samplesize <- .01          # Sample size. Use 1 for all the dataset
printRes   <- TRUE         # Print final results
saveRes    <- TRUE         # Saving results in outdir
outdir     <- 'R log-like studies/output/'
genesis    <- Sys.time()

# Using the dataset
data <- read_csv(file = 'data/ssb_dataset_2018.csv')
data <- na.omit(data)
# Renaming for simplicity
data <- data %>%
  rename(s1 = s_ssb,
         s2 = s_jui,
         s3 = s_mil,
         s4 = s_oth,
         lnp1 = lnp_ssb,
         lnp2 = lnp_jui,
         lnp3 = lnp_mil,
         lnp4 = lnp_oth) %>%
  mutate(age  = log(data$age),
         size = (1/data$size),
         sex  = ifelse(sex == 2, 1, 0)) %>%
  select(s1, s2, s3, s4, lnp1, lnp2, lnp3, lnp4, lnw, age, size, sex, educ)

# Focus on some of the purchasing regimes for testing
data <- data[sample(nrow(data), nrow(data)*samplesize),]

# Vectorizing for simplicity
s1 <- data$s1
s2 <- data$s2
s3 <- data$s3
s4 <- data$s4
lnp1 <- data$lnp1
lnp2 <- data$lnp2
lnp3 <- data$lnp3
lnp4 <- data$lnp4
lnw <- data$lnw
age <- data$age
size <- data$size
sex <- data$sex
educ <- data$educ

# Log likelihood function  
quaids.loglike <- function(param, s1, s2, s3, lnp1, lnp2, lnp3, lnp4, lnw, age, size, educ, sex){
  # Defining the parameters with their parameter restrictions
  # alphas
  a1 <- param[1]
  a2 <- param[2]
  a3 <- param[3]
  a4 <- (1 - a1 - a2 - a3)
  
  # betas
  b1 <- param[4]
  b2 <- param[5]
  b3 <- param[6]
  b4 <- (-b1 - b2 - b3)
  
  # gammas
  g11 <- param[7]  
  g12 <- param[8]
  g13 <- param[9]
  g14 <- (-g11 - g12 - g13)
  
  g21 <- g12
  g22 <- param[10]
  g23 <- param[11]
  g24 <- (-g21 - g22 - g23)
  
  g31 <- g13
  g32 <- g23
  g33 <- param[12]
  g34 <- (-g31 - g32 - g33)
  
  g41 <- g14
  g42 <- g24
  g43 <- g34
  g44 <- (-g41 - g42 - g43)
  
  # lambdas
  l1 <- param[13]
  l2 <- param[14]
  l3 <- param[15]
  l4 <- -(l1 + l2 + l3)
  
  # Thetas: Age
  t11 <- param[16]
  t12 <- param[17]
  t13 <- param[18]
  t14 <- (-t11 - t12 - t13)
  
  # Thetas: size
  t21 <- param[19]
  t22 <- param[20]
  t23 <- param[21]
  t24 <- (-t21 - t22 - t23)
  
  # Thetas: Education
  t31 <- param[22]
  t32 <- param[23]
  t33 <- param[24]
  t34 <- (-t31 - t32 - t33)
  
  # Thetas: sex
  t41 <- param[25]
  t42 <- param[26]
  t43 <- param[27]
  t44 <- (-t41 - t42 - t43)
  
  # Sigmas
  sig11 <- param[28]
  sig12 <- param[29]
  sig13 <- param[30]
  sig14 <- -(sig11 + sig12 + sig13)
  
  sig21 <- sig12
  sig22 <- param[31]
  sig23 <- param[32]
  sig24 <- -(sig21 + sig22 + sig23)
  
  sig31 <- sig13
  sig32 <- sig23
  sig33 <- param[33]
  sig34 <- -(sig31 + sig32 + sig33)
  
  sig41 <- sig14
  sig42 <- sig24
  sig43 <- sig34
  sig44 <- -(sig41 + sig42 + sig43)
  
  #### Bulding the parametric shares
  # Lnpindex (Ln a(p) where a0 = 0)
  lnpindex <- 0 + a1*lnp1 + a2*lnp2 + a3*lnp3 + a4*lnp4
  for(i in as.character(1:4)) {
    for(j in as.character(1:4)) {
      gij <- get(paste0("g", i, j))
      lnpi  <- get(paste0("lnp", i))
      lnpj  <- get(paste0("lnp", j))
      lnpindex <- lnpindex + 0.5*gij*lnpi*lnpj
    }
  }
  # b(p) price index
  bofp <- 0
  for(i in as.character(1:4)) {
    lnpi <- get(paste0("lnp", i))
    bi  <- get(paste0("b", i))
    bofp <- bofp + lnpi*bi
  }
  bofp <- exp(bofp)
  # System of Equations depicted in Eq. 8
  u1 <- a1 + g11*lnp1 + g12*lnp2 + g13*lnp3 + g14*lnp4 + t11*age + t21*size + t31*educ + t41*sex + b1*(lnw - lnpindex - (t11*age + t21*size + t31*educ + t41*sex)) + (l1/bofp)*(lnw - lnpindex - (t11*age + t21*size + t31*educ + t41*sex))^2
  
  u2 <- a2 + g21*lnp1 + g22*lnp2 + g23*lnp3 + g24*lnp4 + t12*age + t22*size + t32*educ + t42*sex + b2*(lnw - lnpindex - (t12*age + t22*size + t32*educ + t42*sex)) + (l2/bofp)*(lnw - lnpindex - (t12*age + t22*size + t32*educ + t42*sex))^2
  
  u3 <- a3 + g31*lnp1 + g32*lnp2 + g33*lnp3 + g34*lnp4 + t13*age + t23*size + t33*educ + t43*sex + b3*(lnw - lnpindex - (t13*age + t23*size + t33*educ + t43*sex)) + (l3/bofp)*(lnw - lnpindex - t13*age + t23*size + t33*educ + t43*sex)^2
  
  u4 <- a4 + g41*lnp1 + g42*lnp2 + g43*lnp3 + g44*lnp4 + t14*age + t24*size + t34*educ + t44*sex + b4*(lnw - lnpindex - (t14*age + t24*size + t34*educ + t44*sex)) + (l4/bofp)*(lnw - lnpindex - (t14*age + t24*size + t34*educ + t44*sex))^2
  
  U <- cbind(u1, u2, u3, u4)
  S <- cbind(s1, s2, s3, s4)
  
  # The variance matrix depicted in Eq. 10
  fullsigma <- c(sig11, sig12, sig13, sig14, sig21, sig22, sig23, sig24, sig31, sig32, sig33, sig34, sig41, sig42, sig43, sig44)
  fullsigma <- matrix(fullsigma, 4, 4)
  
  # Dummy for micro-regimes and index re-arragement
  d1 <- ifelse(s1, 1, 0)
  d2 <- ifelse(s2, 1, 0)
  d3 <- ifelse(s3, 1, 0)
  d4 <- ifelse(s4, 1, 0)
  d <- cbind(d1, d2, d3, d4)
  nu_bght <- rowSums(d) # Vector with number of purchased goods per household 
  
  # Log-likelihood values repository
  lf <- matrix(rep(0, length(s1)), ncol = 1, nrow = length(s1))

  #### Regime where all goods are purchased ###
  index <- nu_bght == 4 # Index for all households that buy all goods
  # Eq. 10
  sigma <- fullsigma[1:3,1:3]
  # Epsilons = Obs. share - param share (depicted in Eq. 18)
  e <- cbind(s1[index] - u1[index], s2[index] - u2[index], s3[index] - u3[index]) 
  ll <- dmvnorm(e, sigma = sigma, log = TRUE)
  
  # Depositing log-likelihood values in their positions
  lf[index] <- ll
  
  for(i in 1:length(s1)) {
    
    bght_index <- which(d[i,] == 1)    # indexing if good is bought
    zero_index <- which(d[i,] == 0)    # indexing if good is not bought
    index <- c(bght_index, zero_index) # index for re-arragement
    
    if(nu_bght[i] < 4) {
      # First, we re-arrange the our vectors and matrices
      sigma <- fullsigma[index, index] # sorted full sigma
      sigma <- sigma[1:3,1:3]          # sorted small sigma
      S_a   <- S[i,index]              # sorted shares
      U_a   <- U[i,index]              # sorted predicted shares
      U_a   <- matrix(U_a, ncol = 1, nrow = 4)
      
      a <- S_a[1:nu_bght[i]]/S_a[1]    # vector
      a <- matrix(a, ncol = 1, nrow = length(a)) 
      
      if(nu_bght[i] == 3) { ### Regime when only three goods are purchased ###
        
        AA1   <- diag(array(a))
        omega <- AA1 %*% solve(sigma) %*% t(AA1)
        s11   <- omega[1:nu_bght[i], 1:nu_bght[i]]
        s11   <- matrix(s11, nrow = length(1:nu_bght[i]), ncol = length(1:nu_bght[i]))
        
        U_bar <- U_a[1]
        II    <- matrix(rep(1,nu_bght[i]), ncol = 1, nrow = 3)
        JJ    <- U_a[1:nu_bght[i]]/(a * U_a[1])
        JJ    <- matrix(JJ, ncol = 1, nrow = 3)
        
        omega11 <- t(II) %*% s11 %*% II
        omega10 <- t(II) %*% s11 %*% JJ
        omega00 <- t(JJ) %*% s11 %*% JJ
        
        U_sta     <- (solve(omega11) %*% omega10)*U_bar
        omg_final <- solve(omega11)
        
        part1 <- t(U_bar) %*% omega00 %*% U_bar - t(U_sta) %*% omega11 %*% U_sta
        part2 <- exp(-0.5*part1)*((2*pi)^(0.5*(1 - nu_bght[i])))*(((det(sigma))^(-0.5))/((det(omg_final))^(-0.5)))
        
        D_2 <- diag((diag(omg_final)^0.5), 4 - nu_bght[i])
        RR <- cov2cor(omg_final)
        
        AA  <- 1/S_a[1]
        CC  <- AA %*% D_2
        PP0 <- 1
        AA[1] <- -AA[1]
        PP <- PP0 + AA*U_sta
        diagD <- rep(0, 4 - nu_bght[i])
        for(j in 1:(4 - nu_bght[i])) {
          diagD[j] <- (CC[j] %*% RR %*% t(CC[j]))^-0.5
        }
        diagD <- as.data.frame(diagD)
        DD <- diag(diagD, ncol = 4 - nu_bght[i], nrow = 4 - nu_bght[i])
        BB <- DD %*% PP
        R_c <- DD %*% CC %*% RR %*% t(CC) %*% DD
        
        # Depositing log-likelihood values in their positions
        HAJP <- log(part2) + log(pmvnorm(lower = c(-Inf,-Inf),
                                         upper = c(as.numeric(-BB), Inf),
                                         sigma = diag(2))[1])
        lf[i] <- HAJP
        
      }else{ ### Regime when less than three goods are purchased ###
        
        ones <- matrix(rep(1,4 - nu_bght[i] - 1), nrow = 4 - nu_bght[i] - 1, ncol = 1)
        AA1 <- rbind(a, ones)
        AA1 <- diag(AA1[1:length(AA1),])
        
        omega <- AA1 %*% solve(sigma) %*% t(AA1)
        
        s11 <- omega[1:nu_bght[i],1:nu_bght[i]]
        s11 <- matrix(s11,nrow = nu_bght[i], ncol = nu_bght[i])
        
        s10 <- omega[1:nu_bght[i],(nu_bght[i] +1):(4-1)]
        s10 <- matrix(s10, nrow = nu_bght[i], 
                      ncol = length((nu_bght[i] +1):(4-1)))
        
        s00 <- omega[(nu_bght[i]+1):(4-1),(nu_bght[i]+1):(4-1)]
        s00 <- matrix(s00, nrow = length((nu_bght[i] +1):(4-1)), 
                      ncol = length((nu_bght[i] +1):(4-1)))
        
        Ubar <- U_a[c(1, (nu_bght[i]+1):(4-1))]
        Ubar <- matrix(Ubar, nrow = length(c(1, (nu_bght[i]+1):(4-1))), ncol = 1)
        II   <- matrix(rep(1,nu_bght[i]), ncol = 1, nrow = nu_bght[i])
        JJ   <- U_a[1:nu_bght[i]] / (a %*% U_a[1])
        
        omega11 <- rbind(cbind(t(II) %*% s11 %*% II, t(II) %*% s10), cbind(t(s10) %*% II, s00))
        omega10 <- rbind(cbind(t(II) %*% s11 %*% JJ, t(II) %*% s10), cbind(t(s10) %*% JJ, s00))
        omega00 <- rbind(cbind(t(JJ) %*% s11 %*% JJ, t(JJ) %*% s10), cbind(t(s10) %*% JJ, s00))
        
        U_sta <- solve(omega11) %*% omega10 %*% Ubar
        
        omg_final <- solve(omega11)
        
        part1 <- t(Ubar) %*% omega00 %*% Ubar - t(U_sta) %*% omega11 %*% U_sta
        part2 <- exp(-0.5*part1)*((2*pi)^(0.5*(1 - nu_bght[i])))*(((det(sigma))^(-0.5))/((det(omg_final))^(-0.5)))
        
        D_2 <- diag((diag(omg_final)^0.5), nrow = 4 - nu_bght[i], ncol = 4 - nu_bght[i])
        RR <- cov2cor(omg_final)
        
        AA     <- diag(-1*rep(1,4 - nu_bght[i]), ncol = 4 - nu_bght[i], nrow = 4 - nu_bght[i])
        ones   <- matrix(rep(1, 4 - nu_bght[i] - 1), ncol = 4 - nu_bght[i] - 1, nrow = 1)
        #ones   <- t(ones)
        AA[1,] <- cbind((1/S_a[1]), ones)
        CC     <- AA %*% D_2
        zeros  <- matrix(rep(0, 4 - nu_bght[i] - 1), ncol = 1, nrow = 4 - nu_bght[i] - 1)
        PP0    <- rbind(1, zeros)
        AA[1,] <- -AA[1,]
        PP     <- PP0 + (AA %*% U_sta) 
        
        diagD <- rep(0, 4 - nu_bght[i])
        diagD <- matrix(diagD, nrow = 1, ncol = 4 - nu_bght[i])
        for(j in 1:(4 - nu_bght[i])) {
          row <- as.matrix(CC[j,])
          row <- t(row)
          diagD[j] <- (row %*% RR %*% t(row))^-0.5
        }
        
        DD  <- diag(array(diagD), ncol = length(diagD), nrow = length(diagD))
        BB  <- DD %*% PP
        R_c <- DD %*% CC %*% RR %*% t(CC) %*% DD
        
        # Depositing log-likelihood values in their positions
        HAJP <- log(part2) + log(pmvnorm(lower = rep(-Inf, length(-as.numeric(BB))), 
                                         upper = as.numeric(-BB), 
                                         corr = R_c)[1])
        lf[i] <- HAJP
      }
    }
  }
  # Return nx1 individual log-likelihood values
  return(lf)
}

# The initial guess (not using all zeros in non-constant values)
b0 <- rep(0, 3)                           # Alphas
b0 <- c(b0, rep(0.003, 3))                # Betas
b0 <- c(b0, 0.01, 0, 0.01, 0, 0, 0.01)    # Gammas
b0 <- c(b0, rep(0.001, 3))                # Lambdas
b0 <- c(b0, rep(0.02, 12))                # Thetas
b0 <- c(b0, 
        1, # Remember to keep a one in diag in vcov matrix
        0, 
        0, 
        1, # Remember to keep a one in diag in vcov matrix
        0, 
        1) # Remember to keep a one in diag in vcov matrix

# Gradient descent algorithm
for(iter in 1:maxIter) {
  # Calculate the gradient (direction of maximum)
  der <- jacobian(f = quaids.loglike,
                  x = b0,
                  s1 = s1,
                  s2 = s2,
                  s3 = s3,
                  lnp1 = lnp1,
                  lnp2 = lnp2,
                  lnp3 = lnp3,
                  lnp4 = lnp4,
                  lnw = lnw,                          
                  age = age,
                  size = size,
                  educ = educ,
                  sex = sex,
                  method = "simple")
  cb <- solve(t(der)%*%der)%*%colSums(der)
  c1 <- abs(sum(der))/length(s1)
  
  # Final reporting
  if(is.na(c1)) {c1 <- 0}
  if(c1 < c0) {
    cat("\n")
    cat("Criteria reached")
    cat("\n")
    cat("Iteration: ", iter)
    cat("\n")
    cat("Likelihood value: ", li1)
    cat("\n")
    cat("Sum of Abs. Grad/HH: ", c1)
    cat("\n")

    break
  }else{
      if(iter == maxIter) {
        cat("\n")
        cat("Maximum number of iterations reached")
        cat("\n")
        cat("Likelihood value: ", li1)
        cat("\n")
        cat("Sum of Abs. Grad/HH: ", c1)
        cat("\n")
        
        break
      }
    }
  # Taking the steps towards the direction marked by the gradient
  stp <- step
  li1 <- 0
  li2 <- 1
  b <- b0
  while(li1 < li2) {
    stp <- stp/2
    li1 <- quaids.loglike(b + stp*cb,                   
                          s1 = s1,
                          s2 = s2,
                          s3 = s3,
                          lnp1 = lnp1,
                          lnp2 = lnp2,
                          lnp3 = lnp3,
                          lnp4 = lnp4,
                          lnw = lnw,
                          age = age,
                          size = size,
                          educ = educ,
                          sex = sex)
    
    li2 <- quaids.loglike(b + .5*stp*cb,                   
                          s1 = s1,
                          s2 = s2,
                          s3 = s3,
                          lnp1 = lnp1,
                          lnp2 = lnp2,
                          lnp3 = lnp3,
                          lnp4 = lnp4,
                          lnw = lnw,
                          age = age,
                          size = size,
                          educ = educ,
                          sex = sex)
    
    li1 <- sum(li1) # Log-likelihood total contribution
    li2 <- sum(li2)
  }
  vcov <- solve(t(der) %*% der) # To be used to calculate SE
  # Updating and reporting
  b0 <- b + cb*stp
  if((iter/reap) %% 1 == 0) {# print when reminder is zero
    cat("\n")
    cat("Iteration: ", iter)
    cat("\n")
    cat("Likelihood value: ", li1)
    cat("\n")
    cat("Sum of Abs. Grad/HH: ", c1)
    cat("\n")
    cat("Estimated parameters:")
    cat("\n")
    cat(b0)
    cat("\n")
  }
}

# Computing the Standard Errors
SE <- sqrt(diag(vcov))

# Saving the results
if(saveRes){
  saveRDS(b, paste0(outdir, 'param_output.rds'))
  saveRDS(SE, paste0(outdir, 'SE_output.rds'))
}
### Printing the results
b <- format(b, digits = 3)
SE <- format(SE, digits = 3)

if(printRes) {
  cat("\n")
  cat("Paramater            (SE)")
  cat("\n")
  cat("a1:   ",  b[1], paste0("(",SE[1],")"))
  cat("\n")
  cat("a2:   ",  b[2], paste0("(",SE[2],")"))
  cat("\n")
  cat("a3:   ",  b[3], paste0("(",SE[3],")"))
  cat("\n")
  cat("b1:   ",  b[4], paste0("(",SE[4],")"))
  cat("\n")
  cat("b2:   ",  b[5], paste0("(",SE[5],")"))
  cat("\n")
  cat("b3:   ", b[6], paste0("(",SE[6],")"))
  cat("\n")
  cat("g11:  ", b[7], paste0("(",SE[7],")"))
  cat("\n")
  cat("g12:  ", b[8], paste0("(",SE[8],")"))
  cat("\n")
  cat("g13:  ", b[9], paste0("(",SE[9],")"))
  cat("\n")
  cat("g22:  ", b[10], paste0("(",SE[10],")"))
  cat("\n")
  cat("g23:  ", b[11], paste0("(",SE[11],")"))
  cat("\n")
  cat("g33:  ", b[12], paste0("(",SE[12],")"))
  cat("\n")
  cat("l1:   ", b[13], paste0("(",SE[13],")"))
  cat("\n")
  cat("l2:   ", b[14], paste0("(",SE[14],")"))
  cat("\n")
  cat("l3:   ", b[15], paste0("(",SE[15],")"))
  cat("\n")
  cat("t11:  ", b[16], paste0("(",SE[16],")"))
  cat("\n")
  cat("t12:  ", b[17], paste0("(",SE[17],")"))
  cat("\n")
  cat("t13:  ", b[18], paste0("(",SE[18],")"))
  cat("\n")
  cat("t21:  ", b[19], paste0("(",SE[19],")"))
  cat("\n")
  cat("t22:  ", b[20], paste0("(",SE[20],")"))
  cat("\n")
  cat("t23:  ", b[21], paste0("(",SE[21],")"))
  cat("\n")
  cat("t31:  ", b[22], paste0("(",SE[22],")"))
  cat("\n")
  cat("t32:  ", b[23], paste0("(",SE[23],")"))
  cat("\n")
  cat("t33:  ", b[24], paste0("(",SE[24],")"))
  cat("\n")
  cat("t41:  ", b[25], paste0("(",SE[25],")"))
  cat("\n")
  cat("t42:  ", b[26], paste0("(",SE[26],")"))
  cat("\n")
  cat("t43:  ", b[27], paste0("(",SE[27],")"))
  cat("\n")
  cat("sig11:", b[28], paste0("(",SE[28],")"))
  cat("\n")
  cat("sig12:", b[29], paste0("(",SE[29],")"))
  cat("\n")
  cat("sig13:", b[30], paste0("(",SE[30],")"))
  cat("\n")
  cat("sig22:", b[31], paste0("(",SE[31],")"))
  cat("\n")
  cat("sig23:", b[32], paste0("(",SE[32],")"))
  cat("\n")
  cat("sig33:", b[33], paste0("(",SE[33],")"))
  cat("\n")
  cat("\n")
  cat("Code started at:", format(genesis, "%a %b %d %X %Y") )
  cat("\n")
  cat("Code finished at:", format(Sys.time(), "%a %b %d %X %Y")
)
}
#end