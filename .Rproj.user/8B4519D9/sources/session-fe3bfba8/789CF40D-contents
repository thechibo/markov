avail_simul <- function(M,init,mat,nmax){
  #
  #it simulates a trajectory of length M and returns
  #the availability function for n=1:nmax
  #M:length of the trajectory
  #init: initial law
  #mat: transition probability matrix
  #nmax: max time for the evaluation of the availability
  #
  avail_estim      <- vector("numeric",length=nmax)
  simul_traj       <- simulmarkov(M,mat)
  trans_prob_estim <- trans_kern_estim_MC(simul_traj)
  #
  for (i in 1:nmax){
    avail_estim[i] <- availability(i,trans_prob_estim,init)
  }
  #
  return(avail_estim)
  #
}

calc_diff_time <- function(n, mylist, mat) {
  #
  #The function returns the difference in the calculation time
  #of the covariance of F^n, for n fixed, between Theorem 4 (Sadek and Limnios)
  #and the new article (George, Samis and Irene)
  #
  #Calculate the covariance of the vector F^n, for every fixed "n", by Sadek and Limnios (Theorem 4)
  #
  trans_prob_theor <- mylist[[2]]
  #
  #
  #pi_bru          <- (trans_prob_theor %^% 100)[1,]          #stationary distribution
  pi_bru           <- mylist[[100]][1,]
  #
  start.time1      <- Sys.time()
  #
  Sigmaf1          <- Sigma(pi_bru,trans_prob_theor,n)        #covariance matrix
  #
  end.time1        <- Sys.time()
  #
  time.taken1      <- end.time1-start.time1
  #
  ########################################################"
  #
  #Calculate the covariance of the vector F^n, for every fixed "n", by Samis, George and Irene, section 3.2.2
  #
  start.time2      <- Sys.time()
  #
  Sigma_pn         <- power_gsi(n,mylist,mat)
  #
  end.time2        <- Sys.time()
  #
  time.taken2      <- end.time2-start.time2
  #
  #
  res <- time.taken2-time.taken1
  #
  return(res)
  #
}

fun_pi <- function(i,u,nstate,init,q,tmax,param){
  #
  res <- 0
  #
  for (r in 1:nstate){
    for (gamma in 0:tmax){
      res <- res + init[r]*elem_riu(r,i,u,gamma,nstate,q,tmax,param)
    }
  }
  #
  return(res)
}
