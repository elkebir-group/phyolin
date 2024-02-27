

bayesFactor <- function(y, N, muL=0.05, sL=10, muB=0.15, sB=10){
  ####
  # 
  #  y: the number of entries in B' that were flipped from 0 in B to 1 in B'
  #  N: the total number of ones in B' where 
  #  muL (muB): estimated false negative rate under a hypothesis of linear (branching) evolution 
  #  sL (sB): precision parameter controlling dispersion of the single-cel technology
  # 
  #
  #  Note: B' is the matrix containing the resultant bits after flipping 
  #
  ######
  
  alphaL <- muL* sL
  betaL <- sL - alphaL
  alphaB <- muB* sB
  betaB <- sB - alphaB
  
  like_linear <- calc_beta_likelihood(y, N, alphaL, betaL)
  like_branched <- calc_beta_likelihood(y, N, alphaB, betaB)
  
  bf <- exp(like_linear - like_branched)
  return(bf)
  
}

calc_beta_likelihood <- function(y, N, alpha, beta){
  
  prior <- lbeta(alpha, beta)
  posterior <- lbeta(alpha + y, beta + N - y)
  return(posterior - prior)
}



