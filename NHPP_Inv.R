#' @title Non-Homogeneous Possion Process Simulation with Inversion (Transformation)
#' @description This function simulates the event times of Non-Homogeneous Possion Process (NHHP).
#' @description It uses Inversion (Transformation) Algorithm
#'
#' @param T The simulation duration of each run
#' @param lam The intensity function of NHPP
#'
#' @return NULL
#' @export  NHPP.inversion
#' @examples NHPP.inversion(T=20,lam = function(x) {0*x+1})
#' @examples NHPP.inversion(T=20,lam = function(x) {3.9*exp(-4*x)})
#' @examples NHPP.inversion(T=20,lam = function(x) {2000*(x+87)^(-2.33)})
#'
NHPP.inversion <- function(T,lam)
{
  # STEP 1 Compute inverse of cumulative intensity function, LamInv(t)
  LamInv                <-Runuran::pinv.new(pdf=lam, lb=0, ub=T,uresolution=1e-010)
  # Step 2.	Compute total area under the intensity function, lam(t) during the simulation period, T, Lam(T)
  UNU.RAN.object.base   <-Runuran::unuran.details(LamInv,show=FALSE, return.list=TRUE)
  Lam                   <-UNU.RAN.object.base$area.pdf
  # STEP 3 Create the NHPP events
  # STEP 3.a. 	nhppEv =  empty set  		# nhppEv ... NHPP events
  nhppEv <- vector()
  # STEP 3.b. n = 0 			              # n . NHPP events
  n <- 0
  # STEP 3.c.	sExp = 0 		              # sExp . cumulative sum of n standard exponentials
  sExp <- 0
  # STEP 3.d.	Repeat
  repeat
  {
    # STEP 3.d.i.	Generate standard Exponential, E
    # STEP 3.d.ii.	sExp = sExp + E
    sExp <- sExp + rexp(1)
    # STEP 3.d.iii.	If (sExp + E > Lam(T)) {Go to Step 4}
    if (sExp  > Lam) {break}
    # STEP 3.d.iv.	n = n +1
    n <- n + 1
    # STEP 3.d.v.	nhppEv[n] = LamInv(sExp)    # nhppEv[n] . occurrence time of the NHPP event n
    nhppEv[n] <- Runuran::uq(LamInv,sExp/Lam)
  }
  # STEP 4.	Return NHPP event set, nhppEv
  return (nhppEv)
}



