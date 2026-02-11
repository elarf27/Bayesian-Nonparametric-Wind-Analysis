source('DtaTrans.R')

library(nimble)
library(coda)
library(MASS)
library(clue)
library(mvtnorm)

N = nrow(df)
# Set the number of mixture components
K = 6 

# Define the nimble model
WindCode <- nimbleCode({
  # Set the DP concentration parameter
  eta ~ dgamma(1, 1)
  # Specify the stick-breaking construction
  for (k in 1:K){
    v[k] ~ dbeta(1, eta)
  }
  # Compute mixture weights from stick-breaking
  RemPrMass[1] <- 1
  for (k in 2:K){
    RemPrMass[k] <- RemPrMass[k - 1] * (1 - v[k - 1])
  }
  for (k in 1:K){
    w[k] <- v[k] * RemPrMass[k]
  }
  # Normalise weights to sum to 1
  Wsum <- sum(w[1:K])
  for (k in 1:K){
    Wnorm[k] <- w[k] / Wsum
  }
  # Specify Bivariate normal parameters for each regime
  for (k in 1:K){
    # Mean vector for each regime
    mu[k, 1:2] ~ dmnorm(mu0[1:2], prec = R0[1:2, 1:2])
    # Precision matrix with conjugate Wishart prior
    prec[k, 1:2, 1:2] ~ dwish(S[1:2, 1:2], df = 3)
  }
  # Specify the likelihood
  for (t in 1:N){
    z[t] ~ dcat(Wnorm[1:K])
    WindVec[t, 1:2] ~ dmnorm(mu[z[t], 1:2], prec = prec[z[t], 1:2, 1:2])
  }
  
})

# Use K-means clustering for initialisation
set.seed(123)
km <- kmeans(cbind(df$u, df$v), centers = K, nstart = 25)
# Specify prior hyperparameters for regime means
mu0 <- c(mean(df$u), mean(df$v))
R0 <- diag(2) /100 # Vague prior means
# Specify Wishart prior for precision
ObsCov <- cov(cbind(df$u, df$v))
S <- solve(ObsCov) / 2 # Scale matrix for Wishart

# Specify model constant
WindConsts <- list(
  N = N,
  K = K,
  mu0 = mu0,
  R0 = R0,
  S = S
)
# Specify model data
WindData <- list(WindVec = cbind(df$u, df$v))
# Specify initial values
InitsFx <- function(){
  # Initialise regime means near K-means centres
  MuInit <- km$centers + matrix(rnorm(K * 2, 0, 2), K, 2)
  # Initialise precision matrices near observed covariance structure
  PrecInit <- array(0, dim = c(K, 2, 2))
  for (k in 1:K){
    PrecInit[k,,] <- solve(ObsCov) + diag(2) * runif(1, 0.1, 0.5)
  }
  # Return list of initial values
  list(
    eta = 1,
    v = rep(0.5, K),
    mu = MuInit,
    prec = PrecInit,
    z = km$cluster
  )
}
# Generate 3 independent initial states for parallel chains
nchains <- 3
WindInits <- lapply(1:nchains, function(x) InitsFx())

# Build the model
WindModel <- nimbleModel(
  WindCode,
  constants = WindConsts,
  data = WindData,
  inits = WindInits[[1]]
) 
# Compile the model
WindComp <- compileNimble(WindModel)
# Configure MCMC with conjugate samplers enabled
ConfStart <- Sys.time()
WindConf <- configureMCMC(WindModel, useConjugacy = TRUE)
ConfEnd <- Sys.time()
ConfRunTime <- difftime(ConfEnd, ConfStart, units = 'mins')
ConfRunTime
# List default samplers assigned by NIMBLE
#WindConf$printSamplers()
# Monitor key parameters
WindConf$addMonitors(c('mu', 'prec', 'eta', 'w', 'Wnorm'))
# Build MCMC algorithm
WindMCMC <- buildMCMC(WindConf)
# Compile MCMC
CWindMCMC <- compileNimble(WindMCMC, project = WindModel)
# Run MCMC with 3 chains
MCMCStart <- Sys.time()
WindSamples <- runMCMC(
  CWindMCMC,
  niter = 50000,
  nburnin = 10000,
  thin = 10,
  nchains = nchains,
  inits = WindInits,
  setSeed = c(123, 456, 516),
  samplesAsCodaMCMC = TRUE
)
MCMCEnd <- Sys.time()
MCMCRuntime <- difftime(MCMCEnd, MCMCStart, units = 'mins')
MCMCRuntime

# MCMC diagnostics
# Gelman-Rubin convergence diagnostics
GrDiag <- gelman.diag(WindSamples, multivariate = FALSE)
MaxRhat <- max(GrDiag$psrf[,1])
nBadRhat <- sum(GrDiag$psrf[,1] > 1.1)
cat('Gelman-Rubin diagnostics:')
cat('Maximum R-hat:', round(MaxRhat, 3))
cat('Parameters with R-hat > 1.1:', nBadRhat, 'out of', nrow(GrDiag$psrf))

# Effective sample size
ess <- effectiveSize(WindSamples)
cat('Effective sample size:')
cat('Minimum:', round(min(ess)), '')
cat('Median', round(median(ess)), '')
cat('Mean', round(mean(ess)), '')
PoorEss <- sum(ess < 100)
cat('Parameters with ESS < 100', PoorEss, 'out of', length(ess))

# Posterior analysis
CombinedChains <- do.call(rbind, lapply(WindSamples, as.matrix))

# Extract parameters
Wcols <- grep('^w\\[', colnames(CombinedChains))
Mucols <- grep('^mu\\[', colnames(CombinedChains))
Preccols <- grep('^prec\\[', colnames(CombinedChains))

PostW <- colMeans(CombinedChains[, Wcols])

# Extract mean vectors
MuMatrix <- matrix(0, K, 2)
for (k in 1:K){
  MuMatrix[k, 1] <- mean(CombinedChains[, paste0('mu[', k, ', 1]')])
  MuMatrix[k, 2] <- mean(CombinedChains[, paste0('mu[', k, ', 2]')])                                     
}

# Extract covariance matrices
CovMatrices <- array(0, dim = c(K, 2, 2))
for (k in 1:K){
  # Get precision matrix samples
  PrecSamples <- CombinedChains[, paste0('prec[', k, ', ', rep(1:2, each = 2), ', ', rep(1:2, times = 2), ']')]
  # Convert to covariance (average over posterior samples)
  nsamples <- nrow(PrecSamples)
  CovSum <- matrix(0, 2, 2)
  for (i in 1:nsamples){
    PrecMat <- matrix(PrecSamples[i, ], 2, 2)
    if (det(PrecMat) > 1e-10){ # Check for numerical issues
      CovSum <- CovSum + solve(PrecMat)
    }
  }
  CovMatrices[k,,] <- CovSum / nsamples
}

# Create a summary
ClusterSummary <- data.frame(
  cluster = 1:K,
  weight = PostW,
  MuU = MuMatrix[, 1],
  MuV = MuMatrix[, 2],
  SigmaU = sqrt(sapply(1:K, function(k) CovMatrices[k, 1, 1])),
  SigmaV = sqrt(sapply(1:K, function(k) CovMatrices[k, 2, 2])),
  rho = sapply(1:K, function(k) CovMatrices[k, 1, 2] / (sqrt(CovMatrices[k, 1, 1]) * sqrt(CovMatrices[k, 2, 2]))),
  MeanSpd = sqrt(MuMatrix[, 1]^2 + MuMatrix[, 2]^2),
  MeanDirDeg = (atan2(MuMatrix[, 2], MuMatrix[, 1]) * 180 / pi + 360) %% 360
) %>% 
  arrange(desc(weight))
ClusterSummary

# Posterior predictive sampling
set.seed(123)
npred <- 2000
ClusterAssignments <- sample(1:K, npred, replace = TRUE, prob = ClusterSummary$weight)
PredU <- numeric(npred)
PredV <- numeric(npred)
for (i in 1:npred){
  k <- ClusterAssignments[i]
  WindVec <- mvrnorm(1, mu = MuMatrix[k,], Sigma = CovMatrices[k,,])
  PredU[i] <- WindVec[1]
  PredV[i] <- WindVec[2]
}

PredWsp <- sqrt(PredU^2 + PredV^2)
PredWdrRad <- atan2(PredV, PredU)
PredWdrRad[PredWdrRad < 0] <- PredWdrRad[PredWdrRad < 0] + 2 * pi

PredData <- data.frame(
  cluster = ClusterAssignments,
  u = PredU,
  v = PredV,
  WSp = PredWsp,
  WdrRad = PredWdrRad,
  WdrDeg = (PredWdrRad * 180 / pi) %% 360
)

# Posterior predictive checks
# Wind speed
cat('Observed mean =', round(mean(df$Wsp), 3), 'sd =', round(sd(df$Wsp), 3))
cat('Predicted mean=', round(mean(PredData$WSp), 3), 'sd =', round(sd(PredData$WSp), 3))
PredErrSpd <- abs(mean(PredData$WSp) - mean(df$Wsp)) / mean(df$Wsp) * 100
cat('Relative error =', round(PredErrSpd, 1))

# Wind vectors (joint structure)
cat('Observed u mean =', round(mean(df$u), 2), 'sd =', round(sd(df$u), 2))
cat('Predicted u mean=', round(mean(PredData$u), 2), 'sd =', round(sd(PredData$u), 2))
cat('Observed v mean =', round(mean(df$v), 2), 'sd =', round(sd(df$v), 2))
cat('Predicted v mean =', round(mean(PredData$v), 2), 'sd =', round(sd(PredData$v), 2))
cat('Observed cor(u, v) =', round(cor(df$u, df$v), 2))
cat('Predicted cor(u, v) =', round(cor(PredData$u, PredData$v), 2))

# Pivotal relabelling
PivRel <- function(SamplesMCMC, K){
  MxSamplesMCMC <- as.matrix(SamplesMCMC)
  niter <- nrow(MxSamplesMCMC)
  
  PivotMu <- matrix(NA, K, 2)
  PivotW <- numeric(K)
  
  for (k in 1:K){
    # Set parameter names
    Mu1Name <- paste0('mu[', k, ', 1]')
    Mu2Name <- paste0('mu[', k, ', 2]')
    Wname <- paste0('w[', k, ']')
    
    # Compute medians
    PivotMu[k, 1] <- median(MxSamplesMCMC[, Mu1Name])
    PivotMu[k, 2] <- median(MxSamplesMCMC[, Mu2Name])
    PivotW[k] <- median(MxSamplesMCMC[, Wname])
  }
  
  # Initialise storage
  relabeled <- MxSamplesMCMC
  perms <- matrix(NA, niter, K)
  distances <- numeric(niter)
  
  # Find optimal permutation for each iteration
  # Extract current iteration parameters
  for (t in 1:niter){
    Mut <- matrix(NA, K, 2)
    Wt <- numeric(K)
    
    for (k in 1:K){
      Mut[k, 1] <- MxSamplesMCMC[t, paste0('mu[', k, ', 1]')]
      Mut[k, 2] <- MxSamplesMCMC[t, paste0('mu[', k, ', 2]')]
      Wt[k] <- MxSamplesMCMC[t, paste0('w[', k, ']')]
      }
    
    # Compute distance matrix
    D <- matrix(NA, K, K)
    for (i in 1:K){
      for (j in 1:K){
        EuclideanDist <- sqrt(sum((Mut[i, ] - PivotMu[j, ])^2))
        weight <- sqrt(Wt[i] * PivotW[j])
        D[i, j] <- weight * EuclideanDist
      }
    }
    # Solve linear assignment problem
    perm <- solve_LSAP(D)
    perms[t, ] <- as.integer(perm)
    distances[t] <- sum(D[cbind(1:K, perm)])
      
    # Apply permutation to relabel parameters
    for (k in 1:K){
      knew <- perm[k]
      
      # Relabel mu
      relabeled[t, paste0('mu[', knew, ', 1]')] <- 
        MxSamplesMCMC[t, paste0('mu[', k, ', 1]')]
      relabeled[t, paste0('mu[', knew, ', 2]')] <- 
        MxSamplesMCMC[t, paste0('mu[', k, ', 2]')]
        
      # Relabel precision
      for (i in 1:2){
        for (j in 1:2){
          PrecNameNew <- paste0('prec[', knew, ', ', i, ', ', j, ']')
          PrecNameOld <- paste0('prec[', k, ', ', i, ', ', j, ']')
          relabeled[t, PrecNameNew] <- MxSamplesMCMC[t, PrecNameOld]
        }
      }
      
      # Relabel weights
      relabeled[t, paste0('w[', knew, ']')] <- 
        MxSamplesMCMC[t, paste0('w[', k, ']')]
        
      # Relabel normalised weights
      relabeled[t, paste0('Wnorm[', knew, ']')] <- 
        MxSamplesMCMC[t, paste0('Wnorm[', k, ']')]
    }
  }
  
  # Compute diagnostics
  # Count iterations requiring relabelling
  nswitches <- sum(apply(perms, 1, function(p) any(p != 1:K)))
  
  # Compute permutation distribution
  PermStrings <- apply(perms, 1, paste, collapse = '-')
  PermTable <- sort(table(PermStrings), decreasing = TRUE)
    
  # Compute permutation entropy
  PermProbs <- PermTable / niter
  PermEntropy <- -sum(PermProbs * log(PermProbs))
  MaxEntropy <- log(factorial(K))
  NormEntropy <- PermEntropy / MaxEntropy
    
  # Return results
  result <- list(
    samples = mcmc(relabeled),
    diagnostics = list(
      permutations = perms,
      permutation_table = PermTable,
      entropy = PermEntropy,
      max_entropy = MaxEntropy,
      normalised_entropy = NormEntropy,
      n_switches = nswitches,
      pct_switches = 100 * nswitches / niter,
      distances = distances,
      pivot_mu = PivotMu,
      pivot_w = PivotW
      )
    )
  return(result)
}

PivResult <- PivRel(
  SamplesMCMC = CombinedChains,
  K = 6
)
PivResult$diagnostics$entropy
PivResult$diagnostics$normalised_entropy
