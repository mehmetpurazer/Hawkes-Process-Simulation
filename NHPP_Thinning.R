#' @title Non-Homogeneous Possion Process Simulation with Thinning Algorithm
#' @description This function simulates the event times of Non-Homogeneous Possion Process (NHHP).
#' @description It uses Thinning Algorithm (Lewis and Shedler, 1979)
#'
#' @param T The simulation duration of each run
#' @param lam The intensity function of NHPP
#'
#' @return NULL
#' @export  NHPP.thinning
#' @examples NHPP.thinning(T=20,lam = function(x) {0*x+1})
#' @examples NHPP.thinning(T=20,lam = function(x) {2000*(x+87)^(-2.33)})
#' @examples NHPP.thinning(T=20,lam = function(x) {3.9*exp(-4*x)})
#'
NHPP.thinning <- function(T,lam)
{
  # STEP 1.	Compute maximum value of lam (t) between time 0 and T, maxlam
  maxlam     <- optimize(lam, interval=c(0,T), maximum=TRUE)[[2]]
  # STEP 2.	lstar = maxlam		# ... upper bound for the thinning algorithm
  lstar       <- maxlam
  # STEP 3. 	NHPPEv =  empty set 	# ... Non-Homogeneous Poisson Process(NHHP) event set
  NHPPEv    <- vector()
  # STEP 4.	 n = 0 						      # . number of events in the set NHPPEv
  n <- 0
  # STEP 5.	curT = 0					      # . current time
  curT <- 0
  # STEP 6.	Repeat
  repeat
  {
    # STEP 6.a.	Generate standard Exponential, E
    # STEP 6.b 	t = E\lstar				# . inter event time
    t <- rexp(1)/lstar
    # STEP 6.c.	curT = curT + t				# update current time
    curT <- curT + t
    # STEP 6.c.	If (curT > T), {Go to Step 7}
    if (curT > T) {break}
    # STEP 6.d.	curInt = lam(curT)		# .  current intensity
    curInt <- lam(curT)
    # STEP 6.e.	Generate standard Uniform, U
    # STEP 6.f 	If U L.E. curInt\lstar
    if (runif(1, min = 0, max = 1) <= curInt / lstar)
      # STEP 6.f.i.	then
    {
      # STEP 6.f.ii.	n = n + 1
      n <- n + 1
      # STEP 6.f.iii.	NHPP[n] = curT	# Add the new event into the NHHP event set
      NHPPEv[n] <- curT
    }
    # STEP 6.g.	Update maximum value of lam(t) between time curT and T, maxlam
    #maxlam     <- optimize(lam, interval=c(curT,T), maximum=TRUE)[[2]]
    # STEP 6.h	lambda star = maxlam
    #lstar <- maxlam
  }
  # STEP 7.	Return Hawkes Process event set, HawkesEv
  return(NHPPEv)
}

