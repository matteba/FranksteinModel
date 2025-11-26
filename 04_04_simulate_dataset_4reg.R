# This script simulates spatio-temporal abundance data for FUOR distinct ecological regions (A, B, C and D), 
# to test species distribution model (SDM) configurations under complex conditions. 
# The simulation is designed to reflect heterogeneous spatio-temporal dynamics: Region A and C follows 
# a progressive dynamic (i.e., spatio-temporal AR(1) ), 
#Region B follows a persistent dynamic (i.e., fixed spatial field with a linear temporal trend),
#Region D follows a opportunistic synamic (i.e., changing spatial structure every time unit with no interannual correlation) .
#
# A regular grid is generated for each region, and spatial autocorrelation is introduced using an 
# exponential covariance function. 


# Abundance observations are simulated using a Gamma distribution. Observational error is included. 
# he final output is a tidy # data frame combining both regions, including spatial coordinates, time, 
and abundance values.
#
# This synthetic dataset is intended for evaluating the fit and predictive capacity of different SDM configurations 
# (persistent, opportunistic, progressive-Type-IV-AR1 and Type-III) and for building the Frankenstein SDM, which integrates 
# region-specific spatio-temporal structures. This script does not include independent variable selection, which should 
# be performed after identifying the best-fitting configuration for each subregion.
library(MASS)# mvrnorm()
{
  set.seed(123)
  ##Parameter -------------------- ##
  n.rep         <- 30          # number of dataset to simulate
  n.time        <- 6           # time unit
  rho           <- 0.8         # rho parameter of AR(1) – region A
  gamma.shape   <- 80         # (≈ NB ‘size’)
  gamma.scale   <- 10          # mean starting value
  mu0           <- gamma.shape*gamma.scale  # ≈ 200 
  trend.A       <- 0.4        # linear trend in region A AND C
  trend.B       <- 0.3        # linear trend in region B
  trend.C       <- -0.4        # linear trend in region A AND C
  #trend.D       <- -0.5     # linear trend in region D
  sigma2       <-0.6
  sigma2_D      <- 1.2           # variance of spatial field
  phi.B         <- 1.2        # range parameter for region B
  phi_D         <- 0.35       # range parameter for region D
  sd.obs        <- 0.1        # observation error parameter (log-normal)
  sd.obs.B  <- 0.02
  sd.obs.D  <- 0.002
  phi.vec   <- seq(0.8, 0.1, length.out = n.time) #vector of range parmeters for region A
  ## --------------------------------------------------------- ##
  
  
  # Regular grid with cell-centre points
grid.coord <- function(xmin, xmax, ymin, ymax, nx, ny) {
  dx <- (xmax - xmin) / nx
  dy <- (ymax - ymin) / ny
  expand.grid(
    x = seq(xmin + dx/2, xmax - dx/2, length.out = nx),
    y = seq(ymin + dy/2, ymax - dy/2, length.out = ny)
  )
}

# 50 pts per quadrant (5 x 10)
A.xy <- grid.coord(0,   0.5, 0,   0.5, nx = 5,  ny = 10)  # region A
B.xy <- grid.coord(0.5, 1.0, 0,   0.5, nx = 5,  ny = 10)  # region B
C.xy <- grid.coord(0,   0.5, 0.5, 1.0, nx = 5,  ny = 10)  # region C
D.xy <- grid.coord(0.5, 1.0, 0.5, 1.0, nx = 5,  ny = 10)  # region D

n.A <- nrow(A.xy); n.B <- nrow(B.xy)
n.C <- nrow(C.xy); n.D <- nrow(D.xy)
  ## spatial covariance matrix (exponential kernel)
  expCov <- function(coords, phi, sigma2) {
    d <- as.matrix(dist(coords))
    sigma2 * exp(-d / phi)
  }  
  ##Simulation scenario
  one.sim <- function(){
    field.A <- matrix(NA, n.time, n.A)
    for(t in 1:n.time){
      phi.A     <- phi.vec[t]
      Sigma.A.t <- expCov(A.xy, phi = phi.A,sigma2)
      chol.A.t  <- chol(Sigma.A.t)
      eps.A.t   <- rnorm(n.A) %*% chol.A.t
      if(t == 1){
        field.A[t, ] <- eps.A.t
      } else {
        field.A[t, ] <- rho * field.A[t-1, ] + sqrt(1 - rho^2) * eps.A.t + trend.A
      }
    }
    ##Region B: linear trend + fixed spatial field
    Sigma.B <- expCov(B.xy, phi.B,sigma2);  chol.B <- chol(Sigma.B)
    S.B     <- as.numeric(rnorm(n.B) %*% chol.B)      # parte spaziale
    time_seq <- 0:(n.time - 1)
    field.B <- sweep(matrix(rep(S.B, each=n.time), nrow=n.time), 1, trend.B*time_seq, `+`)
    
    #region C
    field.C <- matrix(NA, n.time, n.A)
    for(t in 1:n.time){
      phi.C     <- phi.vec[t]
      Sigma.C.t <- expCov(C.xy, phi = phi.C,sigma2)
      chol.C.t  <- chol(Sigma.C.t)
      eps.C.t   <- rnorm(n.C) %*% chol.C.t
      if(t == 1){
        field.C[t, ] <- eps.C.t
      } else {
        field.C[t, ] <- rho * field.C[t-1, ] + sqrt(1 - rho^2) * eps.C.t + trend.C
      }
    }
    #region D
    Sigma.D  <- expCov(D.xy, phi = phi_D, sigma2 = sigma2_D)
    chol.D   <- chol(Sigma.D)
    
    # Draw two very different base fields once
    wA <- as.numeric(rnorm(n.D) %*% chol.D)
    wB <- as.numeric(rnorm(n.D) %*% chol.D)
    
    # Scale them to be clearly distinct (optional but helpful)
    wA <- wA / sd(wA)
    wB <- wB / sd(wB)
    
    # Tiny within-time jitter so they’re not exactly identical copies
    eps_sd <- 0.05   # keep small; reducing this increases separation in WAIC
    
    field.D <- matrix(NA, n.time, n.D)
    for (t in 1:n.time) {
      base <- if ((t %% 2) == 1) wA else wB   # odd: wA, even: wB
      field.D[t, ] <- base + rnorm(n.D, 0, eps_sd)
    }
    
    
    
    
    ##Gamma distribution Observation from NB + observation error
    obs.A <- obs.B <- obs.C<- obs.D<- vector("list", n.time)
    for(t in 1:n.time){
      mu.A <- mu0 * exp(field.A[t, ] + rnorm(n.A, 0, sd.obs))
      mu.B <- mu0 * exp(field.B[t, ] + rnorm(n.B, 0, sd.obs.B))
      mu.C <- mu0 * exp(field.C[t, ] + rnorm(n.C, 0, sd.obs))
      mu.D <- mu0 * exp(field.D[t, ] + rnorm(n.D, 0, sd.obs.D))
      ## Conteggi: Negative Binomial + 1 (per evitare zeri)
      ## In R: 'size' = gamma.shape, 'mu' = media desiderata
      obs.A[[t]] <- rgamma(n.A, shape=gamma.shape, scale=mu.A/gamma.shape)
      obs.B[[t]] <- rgamma(n.B, shape=gamma.shape, scale=mu.B/gamma.shape)
      obs.C[[t]] <- rgamma(n.C, shape=gamma.shape, scale=mu.C/gamma.shape)
      obs.D[[t]] <- rgamma(n.D, shape=gamma.shape, scale=mu.D/gamma.shape)
    }
    
    ## -------- Assembly of data frame “tidy” -------- ##
    datA <- data.frame(
      region    = "A",
      time      = rep(1:n.time, each = n.A),
      A.xy[rep(1:n.A, n.time), ],
      abundance = unlist(obs.A)
    )
    datB <- data.frame(
      region    = "B",
      time      = rep(1:n.time, each = n.B),
      B.xy[rep(1:n.B, n.time), ],
      abundance = unlist(obs.B)
    )
    datC <- data.frame(
      region    = "C",
      time      = rep(1:n.time, each = n.C),
      C.xy[rep(1:n.C, n.time), ],
      abundance = unlist(obs.C)
    )
    datD <- data.frame(
      region    = "D",
      time      = rep(1:n.time, each = n.D),
      D.xy[rep(1:n.D, n.time), ],
      abundance = unlist(obs.D)
    )
    
   rbind(datA, datB,datC,datD)
  }
  sim.list <- lapply(1:n.rep, function(i) one.sim())
}

