## test-chemo.r

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

source('chemo.r')
source('../load.r', chdir=TRUE)
options(digits=16)

## System parameters
noise.dyn <- c(0.1,0.1)   

K.true <- 4.4
chemo.K.prior <- seq(3.5,5.5,0.25)
param.grid <- seq(3.5,5.5,0.1)

start = c(-4,2)

## Filtering-related parameters
noise.obs <- c(0.025,0.025)
nparticles <- 10000
ntries <- 32
obs.dt <- 0.5          ## observe twice / day

## Control-related parameters
## T.list <- c(7.0, 15.0, 30.0)
T.list <- c(4.0)
cntl <- c(0.1, 0.3, 0.5, 0.68)
cntl.default <- 0.3 


## Compute control policy
chemo.control <- function(objective,
                       T,                     
                       xr=c(-5,0),           ## bigger range
                       nx=50,
                       yr=c(-2,6),           ## bigger range
                       ny=50,
                       diff.ratio=1.0,        ## keep dt <= diff.ratio * (dx/noise.dyn)^2
                       max.dt=0.01,           ## limit on max dt
                       do.null=FALSE,         ## compute no-control case
                       prior=c()              ## specify prior param range
                       )
{
  ## Derived parameters
  dx <- (xr[[2]]-xr[[1]])/nx
  dy <- (yr[[2]]-yr[[1]])/ny
  chemo.dt <- min(max.dt, diff.ratio * min(c(dx,dy) / noise.dyn)^2)
  nsteps <- ceiling(T / chemo.dt)
  nskip <- max(0, ceiling(obs.dt / chemo.dt)-1)

  print(paste('dx=',dx,'dy=',dy,'nsteps=',nsteps, 'chemo.dt=',chemo.dt))

  ## Compute control policy
  bellman2d(
            function(delt,c,n,K) {
              chemo.drift.c(delt,c,n, K=K)
            },
            function(delt,c,n,K) {
              chemo.drift.n(delt,c,n, K=K)
            },
            noise.dyn[[1]], noise.dyn[[2]],
            cntl,
            objective,
            xr[[1]], xr[[2]], nx,
            yr[[1]], yr[[2]], ny,
            nsteps, chemo.dt,
            infinite.domain=TRUE,
            params=prior,
            skip=nskip
            )
}

## Generate sample observations and filter away!
do.it <- function () {
  for ( T in T.list ) {

    ## compute and plot control policy
    print(system.time( bc <- chemo.control(chemo.fisher, T, prior=chemo.K.prior)))
    image.grid.fun(bc$mc$grid, bc$cntl[[1]])

    chemo.dt <- bc$dt
    nsteps <- bc$nsteps
    obs.nsteps <- bc$nsteps.actual
    chemo.obs.dt <- bc$dt.actual
    nskip <- round(chemo.obs.dt / chemo.dt)-1
    obs.times <- (seq(1,nsteps,nskip+1)-1) * chemo.dt
    nparams <- length(param.grid)
    
    for ( i.try in 1:ntries ) {

      ## uncontrolled sample path and corresponding observations
      npath <- make.sample.path(chemo.drift.vector, nsteps, chemo.dt, noise.dyn,
                                control.value=cntl.default, start=start, K=K.true)
      nsubpath <- npath[seq(1, nsteps, nskip+1)]
      nfobs <- lapply(nsubpath,function(x) { x+rnorm(2,sd=noise.obs) })
      npobs <- sapply(nsubpath,function(x) {x[1] + noise.obs[[1]] * rnorm(1)})

      ## controlled sample path and corresponding observations
      cpath <- make.controlled.path(chemo.drift.vector, nsteps, chemo.dt, bc, noise.dyn,
                                    start=start, K=K.true)
      csubpath <- cpath[seq(1, nsteps, nskip+1)]
      cfobs <- lapply(csubpath,function(x) { c(x[1],x[2]) + rnorm(2,sd=noise.obs) })
      cpobs <- sapply(csubpath,function(x) {x[1] + noise.obs[[1]] * rnorm(1)})

      ## Partial observations w/o control
      print(system.time( lldat <-
                        pf.loglikelihoods(2, 1, nparticles,
                                          function(t,state,K) {
                                            chemo.drift.vector(cntl.default, state, K=K)
                                          },
                                          chemo.dt,
                                          noise.dyn, noise.obs[[1]],
                                          param.grid,
                                          npobs, obs.times=obs.times,
                                          init.ensemble=matrix(start, 2, nparticles))))

      pp(rowMeans(lldat$loglike))
      
      ## Partial observations + control
      print(system.time( cx <- pf.expt(2, 1, nparticles,
                                       bc,
                                       function(u, state, K=default.chemo.params$K) {
                                         chemo.drift.vector(u, state, K=K)
                                       },
                                       noise.dyn, noise.obs[[1]],
                                       param.grid,
                                       start=start)))

      pp(rowMeans(cx$loglikes))

      ## Partial observations + constant controls
      for ( u in cntl ) {

        ## need new sample paths
        path <- make.sample.path(chemo.drift.vector, nsteps, chemo.dt, noise.dyn,
                                 control.value=u, K=K.true,start=start)
        subpath <- path[seq(1, nsteps, nskip+1)]
        pobs <- sapply(subpath,function(x) {x[1] + noise.obs[[1]] * rnorm(1)})
        print(system.time( lldat <-
                          pf.loglikelihoods(2, 1, nparticles,
                                            function(t,state,K) {
                                              chemo.drift.vector(u, state, K=K)
                                            },
                                            chemo.dt,
                                            noise.dyn, noise.obs[[1]],
                                            param.grid,
                                            pobs, obs.times=obs.times,
                                            init.ensemble=matrix(start, 2, nparticles))))

        pp(rowMeans(lldat$loglike))
      }
    }
  }
}

print(system.time( do.it()))
