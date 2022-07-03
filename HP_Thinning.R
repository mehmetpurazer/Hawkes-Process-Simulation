#' @title Hawkes Process Simulation with Thinning Algorithm
#' @description This function simulates the event times of Hawkes Process(HP).
#' @description It uses the thinning algorigthm (Lewis and Shedler 1979)
#'
#' @param T The simulation duration of each run
#' @param lam0  The base intensity function of Hawkes Process
#' @param phi The excitation function of Hawkes Process
#'
#' @return NULL
#' @export  HP_thinning
#' @examples HP_thinning(T=20,lam0 = function(x) {0*x+1},phi = function(x) {3.9*exp(-4*x)})
#' @examples HP_thinning(T=20,lam0 = function(x) {0*x+1},phi = function(x) {2000*(x+87)^(-2.33)})
#'
HP_thinning <- function(T,lam0,phi)
{

  # STEP 1.	Compute maximum value of lam0 (t) between time 0 and ST, maxlam0
  maxlam0     <- optimize(lam0, interval=c(0,T), maximum=TRUE)[[2]]
  # STEP 2.	Compute maximum value of phi(t) between time 0 and infinity, maxphi
  maxphi      <- optimize(phi, interval=c(0,100), maximum=TRUE)[[2]]
  # STEP 3.	lstar = maxlam0		# ... upper bound for the thinning algorithm
  lstar       <- maxlam0
  # STEP 4. 	HawkesEv =  empty set 					# ... Hawkes Process(HP) event set
  HawkesEv    <- vector()
  # STEP 5.	 n = 0 						# . number of events in the set HawkesEv
  n <- 0
  # STEP 6.	curT = 0					# . current time
  curT <- 0
  # STEP 7.	Repeat
  repeat
  {
    # STEP 7.a.	Generate standard Exponential, E
    # STEP 7.b 	t = E\lstar				# . inter event time
    t <- rexp(1)/lstar
    # STEP 7.c.	curT = curT + t				# update current time
    curT <- curT + t
    # STEP 7.d.	If (curT > T), {Go to Step 8}
    if (curT > T) {break}
    # STEP 7.e.	sumphi = 0	# . sum of excitations triggered by the events in the process history
    sumphi <- 0
    # STEP 7.f. If (n > 0)
    if (n > 0)
    {
      # STEP 7.f.i. For i = 1 to n {sumphi = sumphi + phi (curT - HawkesEv[n])}
      for (n in 1:length(HawkesEv))
      {
        sumphi <- sumphi + phi (curT - HawkesEv[n])
      } 
    }
    # STEP 7.g.	curInt = lam0(curT) + sumphi		# .  current intensity
    curInt <- lam0(curT) + sumphi
    # STEP 7.h.	Generate standard Uniform, U
    # STEP 7.i 	If U L.E. curInt\lstar
    if (runif(1, min = 0, max = 1) <= curInt / lstar)
      # STEP 7.i.i.i.	then 
    {
      # STEP 7.i.i.1	n = n + 1
      n <- n + 1
      # STEP 7.i.i.2	HawkesEv[n] = curT	# Add the new event into the HP event set
      HawkesEv[n] <- curT
      # STEP 7.i.i.3.	Update maximum value of lam0(t) between time curT and T, maxlam0
      maxlam0     <- optimize(lam0, interval=c(curT,T), maximum=TRUE)[[2]]
      # STEP 7.i.i.4.	lambda star = maxlam0 + maxphi + sumphi
      lstar <- maxlam0 + maxphi + sumphi
    }
    # STEP 7.i.i.ii.	else 
    else
    {
      # STEP 7.i.i.ii.1 Update maximum value of lam0(t) between time curT and T, maxlam0
      maxlam0     <- optimize(lam0, interval=c(curT,T), maximum=TRUE)[[2]]
      # STEP 7.i.i.ii.2.	lambda star = maxlam0 + sumphi
      lstar <- maxlam0 + sumphi
    }
  }
  # STEP 8.	Return Hawkes Process event set, HawkesEv 
  return(HawkesEv)
}


