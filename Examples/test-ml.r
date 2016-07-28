## test-ml.r

## This file is part of RSTOC, a system for stochastic
## optimal control and experimental design written in R.

## Copyright (C) 2012 by the authors:

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

## Compute optimal control for the ML example, then use it
## to generate log-likelihoods for sample paths.

source('ml.r')
source('../load.r', chdir=TRUE)

## System parameters
noise.dyn <- c(1.0,0.01)   ## Ermentrout-Terman use sigma_v=1, not sure why
                           ## Also, 1e-3 = 1/sqrt(1e6), roughly #(ion channels)?
gca.true <- 4.41498308
ml.gca.prior <- param.grid <- seq(3,5,0.25)

## Filtering-related parameters
noise.obs <- c(0.2,0.01)
nparticles <- 16000
ntries <- 32

## Control-related parameters
cntl <- c(0,3.5,5)
cntl.default <- 3.5


## Compute control policy
ml.control <- function(objective,
                       T,                     ## all time mesaured in ms
                       xr=c(-80,80),          ## voltage range
                       nx=70,
                       yr=c(0,1),             ## gating variable should stay in [0,1]
                       ny=70,
                       diff.ratio=1.0,        ## keep dt <= diff.ratio * (dx/noise.dyn)^2
                       max.dt=0.5,            ## limit on max dt
                       do.null=FALSE,         ## compute no-control case
                       prior=c()              ## specify prior param range
                       )
{
  ## Derived parameters
  dx <- (xr[[2]]-xr[[1]])/nx
  dy <- (yr[[2]]-yr[[1]])/ny
  ml.dt <- min(max.dt, diff.ratio * min(c(dx,dy) / noise.dyn)^2)
  nsteps <- ceiling(T / ml.dt)

  print(paste('dx=',dx,'dy=',dy,'nsteps=',nsteps, 'ml.dt=',ml.dt))

  ## Compute control policy
  bellman2d(
            function(u,v,w,gca) {
              ml.drift.v(u, v, w, g_ca=gca)
            },
            function(u,v,w,gca) {
              ml.drift.w(u, v, w, g_ca=gca)
            },
            noise.dyn[[1]],
            function(v,w) {
              channel.noise.coeff(v,w,channel.noise.amp=noise.dyn[[2]])
            },
            cntl,
            objective,
            xr[[1]], xr[[2]], nx,
            yr[[1]], yr[[2]], ny,
            nsteps, ml.dt,
            infinite.domain=TRUE,
            params=prior
            )
}


## multiplicative noise
ml.noise.dyn <- function(u, state, ...) {
  c(noise.dyn[[1]],
    channel.noise.coeff(state[[1]], state[[2]],
                        channel.noise.amp=noise.dyn[[2]]))
}

## Generate sample observations and filter away!
T.list <- c(10.0, 40.0)

do.it <- function () {
  for ( T in T.list ) {

    ## compute and plot control policy
    print(system.time( bc <- ml.control(ml.fisher.gca, T, prior=ml.gca.prior)))
    image.grid.fun(bc$mc$grid, bc$cntl[[1]])

    ml.dt <- bc$dt
    nsteps <- bc$nsteps.actual
    
    for ( i.try in 1:ntries ) {

      ## uncontrolled sample path and corresponding observations
      npath <- make.sample.path.mult(2, ml.drift.vector, nsteps, ml.dt, ml.noise.dyn,
                                     control.value=cntl.default, g_ca=gca.true)
      nfobs <- lapply(npath,function(x) { x+rnorm(2,sd=noise.obs) })
      npobs <- sapply(npath,function(x) {x[1]}) + noise.obs[[1]] * rnorm(length(npath))

      ## controlled sample path and corresponding observations
      cpath <- make.controlled.path.mult(2, ml.drift.vector, nsteps, ml.dt, bc,
                                         ml.noise.dyn, g_ca=gca.true)
      cfobs <- lapply(cpath,function(x) { c(x[1],x[2]) + rnorm(2,sd=noise.obs) })
      cpobs <- sapply(cpath,function(x) {x[1]}) + noise.obs[[1]] * rnorm(length(cpath))

      ## Full observations w/o control
      print(system.time( lldat <-
                        sde.loglikelihoods.mult(2,
                                                function(t,state,gca) {
                                                  ml.drift.vector(cntl.default, state, g_ca=gca)
                                                },
                                                ml.dt,
                                                function(state) {
                                                  ml.noise.dyn(cntl.default, state)
                                                },
                                                param.grid,
                                                npath)))

      pp(rowMeans(lldat$loglike))
      
      ## Full observations + control
      print(system.time( lldat <-
                        sde.loglikelihoods.mult(2,
                                                function(t,state,gca) {
                                                  x <- state[1]
                                                  y <- state[2]
                                                  u <- bellman2d.control.value(bc, t, x, y)
                                                  ml.drift.vector(u, c(x,y), g_ca=gca)
                                                },
                                                ml.dt,
                                                function(state) {
                                                  ml.noise.dyn(NULL, state)
                                                },
                                                param.grid,
                                                lapply(cpath,function(x) {c(x[1],x[2])}))))
      
      pp(rowMeans(lldat$loglike))

      ## Partial observations w/o control
      print(system.time( lldat <-
                        pf.loglikelihoods.mult(2, 1, nparticles,
                                               function(t,state,gca) {
                                                 ml.drift.vector(cntl.default, state, g_ca=gca)
                                               },
                                               ml.dt,
                                               function(state) {
                                                 ml.noise.dyn(NULL, state)
                                               },
                                               noise.obs[[1]],
                                               param.grid,
                                               npobs)))

      pp(rowMeans(lldat$loglike))
      
      ## Partial observations + control
      print(system.time( cx <- pf.expt.mult(2, 1, nparticles,
                                            bc,
                                            function(u, state, theta=gca.true) {
                                              ml.drift.vector(u, state, g_ca=theta)
                                            },
                                            function(state) {
                                              ml.noise.dyn(NULL, state)
                                            },
                                            noise.obs[[1]],
                                            param.grid)))

      pp(cx$loglike)

      ## Partial observations + static controls
      for ( u in cntl ) {

        ## need new sample paths
        path <- make.sample.path.mult(2, ml.drift.vector, nsteps, ml.dt, ml.noise.dyn,
                                      control.value=u, g_ca=gca.true)
        pobs <- sapply(path,function(x) {x[1]}) + noise.obs[[1]] * rnorm(length(path))
        print(system.time( lldat <-
                          pf.loglikelihoods.mult(2, 1, nparticles,
                                                 function(t,state,gca) {
                                                   ml.drift.vector(u, state, g_ca=gca)
                                                 },
                                                 ml.dt,
                                                 function(state) {
                                                   ml.noise.dyn(NULL, state)
                                                 },
                                                 noise.obs[[1]],
                                                 param.grid,
                                                 pobs)))

        pp(rowMeans(lldat$loglike))
      }
    }
  }
}

print(system.time( do.it()))
