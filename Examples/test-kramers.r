## test-kramers.r

## This file is part of RSTOC, a system for stochastic
## optimal control and experimental design written in R.

## Copyright (C) 2013 by the authors:

## Giles Hooker <gjh27@cornell.edu>
## Kevin K Lin <klin@math.arizona.edu>
## Bruce Rogers <bruce@math.duke.edu>

## This program is free software; you can redistribute it
## and/or modify it under the terms of the GNU General
## Public License as published by the Free Software
## Foundation; either version 2 of the License, or (at your
## option) any later version.

## This program is distributed in the hope that it will be
## useful, but WITHOUT ANY WARRANTY; without even the
## implied warranty of MERCHANTABILITY or FITNESS FOR A
## PARTICULAR PURPOSE.  See the GNU General Public License
## for more details.

## You should have received a copy of the GNU General Public
## License along with this program; if not, write to the
## Free Software Foundation, Inc., 51 Franklin Street, Fifth
## Floor, Boston, MA 02110-1301 USA.

##--------------------------------------------------------##

## This code showcases the Kramers example.

source('kramers.r')
source('../load.r', chdir=TRUE)

## Dynamical parameters
noise.dyn <- 0.1                      ## noise amplitude
A.true <- default.kramers.params$A    ## true value of A
start = 0                             ## initial condition

## Control and filtering-related parameters
noise.obs <- 0.05                     ## observation noise
obs.dt <- 0.25                        ## observation time interval
nparticles <- 10000                   ## #(particles) for particle filter
ntries <- 4                           ## #(repeated trials)
param.grid <- seq(2.0, 5.0, 0.3)      ## param grid for MLE
kramers.prior <- seq(2.0, 5.0, 0.3)   ## prior interval

## Control-related parameters
cntl <- seq(-10,10,2.0)               ## set of control values
cntl.default <- 0.0               ## set of control values

## Subroutine for computing control policy
kramers.control <- function(objective,
                            T,                     
                            xr=c(-5,5),
                            nx=100,
                            diff.ratio=1.0,        ## keep dt <= diff.ratio * (dx/noise.dyn)^2
                            max.dt=0.01,           ## limit on max dt
                            do.null=FALSE,         ## compute no-control case
                            prior=c()              ## specify prior param range
                            )
{
  ## Derived parameters
  dx <- (xr[[2]]-xr[[1]])/nx
  kramers.dt <- min(max.dt, diff.ratio * dx / noise.dyn^2)
  nsteps <- ceiling(T / kramers.dt)
  nskip <- max(0, ceiling(obs.dt / kramers.dt)-1)

  print(paste('dx=',dx,'nsteps=',nsteps, 'kramers.dt=',kramers.dt))

  ## Compute control policy
  bellman1d(
            function(u,x,A) {
              kramers.drift(u,x, A=A)
            },
            noise.dyn,
            cntl,
            objective,
            xr[[1]], xr[[2]], nx,
            nsteps, kramers.dt,
            infinite.domain=TRUE,
            params=prior,
            skip=nskip
            )
}


## generate sample observations and filter away!
## outputs are log-likelihood estimates for each value in param.grid
## these can be post-processed to obtain the data shown in the paper
T.list <- c(4.0, 30.0)

do.it <- function () {
  for ( T in T.list ) {

    ## compute and plot control policy
    print(system.time( bc <- kramers.control(kramers.fisher, T, prior=kramers.prior)))
    plot(sapply(bc$mc$grid$cells, function(z) { z$xm }),
         bc$cntl[[1]],
         main=paste('T=', T),
         xlab='x',
         ylab='u')

    ## set params
    kramers.dt <- bc$dt
    nsteps <- bc$nsteps
    obs.nsteps <- bc$nsteps.actual
    kramers.obs.dt <- bc$dt.actual
    nskip <- round(kramers.obs.dt / kramers.dt)-1
    obs.times <- (seq(1,nsteps,nskip+1)-1) * kramers.dt
    nparams <- length(param.grid)
    
    for ( i.try in 1:ntries ) {

      ########################################################
      ## Uncontrolled
      
      ## uncontrolled sample path and corresponding observations
      npath <- make.sample.path(kramers.drift, nsteps, kramers.dt, noise.dyn,
                                control.value=cntl.default, start=start, A=A.true)
      nsubpath <- npath[seq(1, nsteps, nskip+1)]
      nfobs <- sapply(nsubpath,function(x) {x[1] + noise.obs[[1]] * rnorm(1)})

      ## controlled sample path and corresponding observations
      cpath <- make.controlled.path(kramers.drift, nsteps, kramers.dt, bc, noise.dyn,
                                    start=start, A=A.true)
      csubpath <- npath[seq(1, nsteps, nskip+1)]
      cfobs <- sapply(csubpath,function(x) {x[1] + noise.obs[[1]] * rnorm(1)})

      ## Full observations w/o control
      print(system.time( lldat <-
                        sde.loglikelihoods(1,
                                           function(t,state,A) {
                                             kramers.drift(cntl.default, state, A=A)
                                           },
                                           kramers.dt,
                                           noise.dyn,
                                           param.grid,
                                           npath)))

      pp(rowMeans(lldat$loglike))
      
      ## Full observations + control
      print(system.time( lldat <-
                        sde.loglikelihoods(1,
                                           function(t,state,A) {
                                             x <- state[1]
                                             u <- bellman1d.control.value(bc, t, x)
                                             kramers.drift(u, x, A=A)
                                           },
                                           kramers.dt,
                                           noise.dyn,
                                           param.grid,
                                           lapply(cpath,function(x) {c(x[1])}))))
        
      pp(rowMeans(lldat$loglike))

      ## Partial observation + control
      print(system.time( cx <- pf.expt(1, 1, nparticles,
                                       bc, kramers.drift,
                                       noise.dyn, noise.obs,
                                       param.grid)))

      pp(rowMeans(cx$loglikes))
      
      ## Static controls
      for ( u in cntl ) {

        ## need new sample paths
        path <- make.sample.path(kramers.drift, nsteps, kramers.dt, noise.dyn,
                                 control.value=u, A=A.true,start=start)
        subpath <- path[seq(1, nsteps, nskip+1)]
        fobs <- sapply(subpath,function(x) {x[1] + noise.obs[[1]] * rnorm(1)})
        
        print(system.time( lldat <-
                          pf.loglikelihoods(1, 1, nparticles,
                                            function(t,state,A) {
                                              kramers.drift(u, state, A=A)
                                            },
                                            kramers.dt,
                                            noise.dyn, noise.obs,
                                            param.grid,
                                            fobs, obs.times=obs.times,
                                            init.ensemble=matrix(start, 1, nparticles))))

        pp(rowMeans(lldat$loglike))
      }
    }
  }
}

print(system.time( do.it()))
