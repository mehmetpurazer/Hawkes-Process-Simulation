#' @title Non-Homogeneous Possion Process Simulation with Inversion (Transformation) combined with Order Statistics
#' @description This function simulates the event times of Non-Homogeneous Possion Process (NHHP).
#' @description It uses Inversion (Transformation) Algorithm
#' @description It is combined with Order Statistics to prevent the events to be rejected
#'
#' @param T The simulation duration of each run
#' @param lam  The intensity function of NHPP
#'
#' @return NULL
#' @export NHPP.inversion.order.statistics
#' @examples NHPP.inversion.order.statistics(T=20,lam = function(x) {0*x+1})
#' @examples NHPP.inversion.order.statistics(T=20,lam = function(x) {3.9*exp(-4*x)})
#' @examples NHPP.inversion.order.statistics(T=20,lam = function(x) {2000*(x+87)^(-2.33)})
#'
NHPP.inversion.order.statistics <- function(T,lam)
{
  # STEP 1 Compute inverse of cumulative intensity function, LamInv(t)
  LamInv                <-Runuran::pinv.new(pdf=lam, lb=0, ub=T,uresolution=1e-010)
  # Step 2.	Compute total area under the intensity function, lam(t) during the simulation period, T, Lam(T)
  UNU.RAN.object.base   <-Runuran::unuran.details(LamInv,show=FALSE, return.list=TRUE)
  Lam                   <-UNU.RAN.object.base$area.pdf
  
  # STEP 3 	Get the number of events, n = Poi(1*Lam)
  n <- rpois(1,Lam)
  
  # STEP 4	Create n items uniformly distributed on (0,Lam), u
  # The items are cumulative sum of standard exponentials
  u <-  runif(n,min = 0, max = Lam)
  
  # STEP 5	Create the NHPP events
  nhppEv <- Runuran::uq(LamInv,u/Lam)
  
  
  # STEP 6 Return NHPP event set sorted, sort (nhppEv)
  return (sort (nhppEv,decreasing = FALSE))
  
  #return(Lam)
}


