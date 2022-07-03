#' @title Hawkes Process Simulation with Composition Algorithm with Runuran
#' @description This function simulates the event times of Hawkes Process(HP).
#' @description It is based on cluster-based representation of Hawkes Process
#' @description The immigrant process and the offspring processes of HP are
#' @description independent Non-Homogeneous Poisson Processes (NHPP).
#' @description The Hawkes Process is the combination of these NHPPs with
#' @description different parent event times.
#' @description Therefore it is a combination algorithm.
#' @description The NHPPs are are simulated using transformation method with order statistics
#' @description The R package RUNURAN is used for the numerical inversion of
#' @description the cumulative intensity function
#'
#' @param n The number of independent Hawkes Process simulations
#' @param T The simulation duration of each run
#' @param lam0  The base intensity function of Hawkes Process
#' @param phi The excitation function of Hawkes Process
#'
#' @return NULL
#' @export  HP_Composition_R
#' @examples HP_Composition_R(n=1,T=20,lam0 = function(x) {0*x+1},phi = function(x) {3.9*exp(-4*x)})
#' @examples HP_Composition_R(n=10,T=20,lam0 = function(x) {0*x+1},phi = function(x) {2000*(x+87)^(-2.33)})
#'
HP_Composition_R <- function(n,T,lam0,phi)
{
  #To ignore the warnings during usage
  options(warn=-1)
  options("getSymbols.warning4.0"=FALSE)
  # STEP 1.	Compute inverse of cumulative base intensity function, Lam0Inv(t)
  # Lam0Inv                <-Runuran::pinv.new(pdf=lam0, lb=0, ub=T,uresolution=1e-010)
  # Step 2.	Compute inverse of cumulative excitation function, PhiInv(t)
  # PhiInv                 <-Runuran::pinv.new(pdf=phi, lb=0, ub=T,uresolution=1e-010)
  # STEP 3.	Compute the cumulative base intensity function, Lam0(t)

  # STEP 4.	Compute the cumulative excitation function, Phi(t)

  HawkesEv <- list()

  # simulates Hawkes Process for n times
  for (re in (1:n))
  {
    # STEP 5.	Create immigrant generation (generation 1)
    # STEP 5.a.	Create the immigrant events,
    imEv <- vector()
    imEv <- NHPP_Simulation_R (0,lam0,T)
    # imEv . occurrence time of the immigrants
    # STEP 5.b.	i = 1		# immigrant events are generation 1
    i <- 1
    # STEP 5.c.	Generation [1] = imEv
    Generation <- list()
    Generation[[i]] <- imEv
    # STEP 6.  Create offspring generations
    repeat
    {
      # STEP 6.a.	if (length(generation[[i]] == 0) break
      if (length(Generation[[i]]) == 0) {break}
      # There is no event in the generation to produce off springs
      # STEP 6.b.	NHPP start times that are the parent event times
      NHPP.start <- vector()
      NHPP.start <- Generation[[i]]
      # STEP 6.c.	Remaining time for offspring production, r.time = T- Parent event time
      r.time <- vector()
      r.time <- T - NHPP.start
      # STEP 6.d.	Offspring event times
      Off.Ev <- list()
      # STEP 6.d.i.	for j = 1 to length(Parent event time)
      # Offspring event times = Offspring event times +
      #          NHPP (Parent event time[j],phi(t),Phi(t),PhiInv(t),simulation duration = r.time[j])

      for (j in 1:length(r.time))
      {
        Off.Ev[[j]] <- NHPP_Simulation_R(NHPP.start[j],phi,r.time[j])
      }

      Off.Ev.1 <- vector()
      Off.Ev.1 <- unlist(Off.Ev)
      # STEP 6.e.if no off springs are produced then break
      if(length(Off.Ev.1) == 0) break
      # STEP 6.f.	Next generation event times
      i <- i+1
      Generation[[i]] <- Off.Ev.1
    }

    # STEP 7.	Return Hawkes Process event set
    # STEP 7.a.	HawkesEv = sort (unlist(Generation))	# HawkesEv ... Hawkes Process events

    # return (HawkesEv <- sort(unlist(Generation),decreasing = FALSE))
    HawkesEv[[re]] <- sort(unlist(Generation),decreasing = FALSE)
  }
  return (HawkesEv)
}

NHPP_Simulation_R <- function (start,intensity, sim.dur)
  #  start            : NHPP start time, parent.event.time
  #  intensity        : intensity function, intensity(t)
  #  cum.intensity    : cumulative intensity function, Intensity(t)
  #  parent.event.time: inverse cumulative intensity function, InvIntensity(t)
  #  sim.dur          : simulation duration, ST

{
  # STEP 1.	Calculate the area under the intensity function during the simulation period,
  area <- vector()

  InvCumIntensityInv          <-Runuran::pinv.new(pdf=intensity, lb=0, ub=sim.dur,uresolution=1e-010)
  UNU.RAN.object.intensity    <-Runuran::unuran.details(InvCumIntensityInv,show=FALSE, return.list=TRUE)
  area                        <-UNU.RAN.object.intensity$area.pdf


  # STEP 2.	Generate the number of events, n from the Poisson distribution with  mean "area"
  n <- vector()
  n <- rpois(1,area)
  # STEP 3.	Create n items uniformly distributed on (0, area), u
  #The items are cumulative sum of standard exponentials
  u <- vector()
  u <- runif(n,min = 0, max = area)
  # STEP 4.	Create the NHPP events, NHPP.Ev = inv.cum.intensity(u)  + nhpp.start
  # NHPP.Ev . occurrence time of the events
  NHPP.Ev         <- vector()
  interevent.time <- vector()
  interevent.time <- Runuran::uq(InvCumIntensityInv,u/area)
  NHPP.Ev         <- interevent.time + start
  return (NHPP.Ev)
}