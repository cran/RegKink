# File name: RegKink.R
# It uses to estimate the regression kink with a state-dependent threshold
# proposed by Yang and Su(2018, JIMF)
# Lixiong Yang, Jen-Je Su. Debt and growth: Is there a constant tipping point?[J].
# Journal of International Money and Finance, 2018, 87:133-143.



# Some useful functions: compute OLS estimates
reg <- function(X,y) {
  X <- qr(X)
  bols <- as.matrix(qr.coef(X,y))
  return(bols)
}

pos.part <- function(x) {
  ps <- x*(x>0)
 return(ps)
}  #positive part

neg.part <- function(x){
  ne <- x*(x<0)
 return(ne)} #negative part



rkt <- function(y,x,z,q, r01, r02, r11, r12,
                stp1, stp2) {
  # model y = b1 (x-r_t)_- + b2 (x-r_t)_+ + b2*z + e
  # r_t=r0+r1 q_1
  # Input:
  #   y: outcome variable
  #   x: threshold dependent variable
  #   z: control variable
  #   q: state variable affecting threshold
  #   stp1:
  #   r01:
  #   r01:
  # gammas0 = seq(r01,r02,by=stp1)	# Grid on Threshold parameter for estimation
  #   stp2:
  #   r11:
  #   r12:
  # gammas1 = seq(r11,r12,by=stp2)	# Grid on Threshold parameter for estimation
  #   L.bt / U.bt: Lower and upper bounds for the threshold parameters gamma.
  #   L.dt / U.dt: Lower and upper bounds for dt
  #   L.gm / U.gm: Lower and upper bounds for gm
  #   tau1: Lower bound for the proportion of regime 1, i.e. (f'gm > 0)
  #   tau2: Upper bound for the proportion of regime 1, i.e. (f'gm > 0)
  #   eta: effective zero
  #   params: parameters for gurobi engine


  # Linear Model
  q1 = q
  x0 = cbind(x,z)
  k0 = ncol(x0)
  kz = ncol(z)
  x00 = solve(crossprod(x0))  #inverse of x0'x0
  bols = x00%*%crossprod(x0,y)
  e0 = y - x0%*%bols
  sse0 = sum(e0^2)
  n = length(y)
  sigols = sse0/n
  v0 = x00%*%crossprod(x0*matrix(e0,n,k0))%*%x00*(n/(n-k0))
  seols = as.matrix(sqrt(diag(v0)))


  # kink Model with a state-dependent threshold
  gammas0 = seq(r01,r02,by=stp1)	# Grid on Threshold parameter for r0
  gammas1 = seq(r11,r12,by=stp2)
  grid0 = length(gammas0)
  grid1 = length(gammas1)

  sse = matrix(0,grid0,grid1)
  k = kz + 3
  for (j in 1:grid0) {
    for (i in 1:grid1){
      rt=gammas0[j] + gammas1[i]*q1
      x1 = cbind(neg.part(x-rt),pos.part(x-rt),z)
      #############
      np <- neg.part(x-rt)
      pp <- pos.part(x-rt)
      nm1 <- length(np[which(np==0)])
      nm2 <- length(pp[which(pp==0)])
      ##############

      e1 = y - x1%*%reg(x1,y)
      sse[j,i] = sum(e1^2)
      # unsatistified if too few observations are in one regime
      ifelse(nm1<0.15*n | nm2<0.15*n,sse0,sum(e1^2))
    }
  }
  gi = which(sse==min(sse),arr.ind=T)
  jmin = gi[1]
  imin = gi[2]
  gammahat0 = gammas0[jmin]
  gammahat1 = gammas1[imin]
  ssemin = sse[jmin,imin]
  rthat = gammahat0 + gammahat1*q1
  x1 = cbind(neg.part(x-rthat),pos.part(x-rthat),z)
  bt = reg(x1,y)
  et = y - x1%*% bt

  betahat = rbind(bt,gammahat0, gammahat1)
  sig = sum(et^2)/n




  return(list(bols=bols, bt=bt, gammahat0=gammahat0, gammahat1=gammahat1, sig=sig))
}







