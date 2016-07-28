## pf.r

## This file is part of RSTOC, a system for stochastic
## optimal control and experimental design written in R.

## Copyright (C) 2012,2013,2014,2015,2016 by the authors:

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

################################
## This code implements a particle filter for SDEs with
## *linear* observations.  The SDE is assumed to be
## time-discretized using Euler, so transition densities are
## gaussian.  Linearity of observations can be exploited to
## provide exact samples, i.e., we can avoid the resampling
## needed in, e.g., bootstrap particle filters.  This is
## both faster and performs better.

## Unlike its predecessor, this code operates in any number
## of dimensions.  It can be used for parameter- as well as
## state-estimation problems.  It can also be used as the
## basis for more complicated algorithms, e.g., MIF.

## The algorithm used here is fully general, but for speed,
## we assume the diffusion matrix is diagonal.  We also
## assume that the observation matrix is a canonical
## projection onto the first DIM.OBS components.

## Note to self: See notebook #15p136-8 for general
## formulae.


################################
## make a new filter object

make.pf <- function(dim.state,  # state space dimension
                    dim.obs,    # observation space dimension
                    n.particles,  # num particles
                    drift0,     # drift=drift(state) or drift(t,state)
                    
                    dyn.noise.amps,   # SDE diffusion coeff matrix
                                      # assumed diagonal ie a vector
                    
                    obs.noise.amps,         # observational noise
                    
                    init.ensemble=FALSE # optional: specify an initial ensemble
                                        # default: all particles start at 0
                    ) {

  ## Store ensemble in a matrix; each state vector is a
  ## column
  if ( is.matrix(init.ensemble)) {
    ensemble <- init.ensemble
  }
  else {
    ensemble <- matrix(0, dim.state, n.particles)
  }

  ## Check if drift function depends on time.
  if ( length(formals(drift0)) == 1 ) {
    drift <- function(t, ensemble) {
      drift0(ensemble)
    }
  }
  else {
    drift <- drift0
  }

  ## Pack up everything in a filter function whose job is
  ## given a new observation, update the ensemble and output
  ## the expected state of the system.

  ## C <- dyn.noise.amps^2 * dt
  ## C.noise <- obs.noise.amps^2 * dt
  C <- dyn.noise.amps^2
  C.noise <- obs.noise.amps^2

  ## coefficients for conditional means
  A <- C[1:dim.obs] / (C[1:dim.obs] + C.noise)
  B <- C.noise / (C[1:dim.obs] + C.noise)

  ## compute conditional variances
  D <- vector(mode='numeric',dim.state)
  D[1:dim.obs] <- 1/(C[1:dim.obs] + C.noise)
  cond.vars <- C - C*D*C
  ## In the general case, this should be a matrix square root
  cond.scales <- sqrt(cond.vars)
  diffusion.scales <- sqrt(C)

  ## Track current time
  t <- 0

  ## Evolve the ensemble forward by one time step, without
  ## taking into account any observations.
  step <- function(new.t, randomness=NULL) {

    if ( new.t <= t ) {
      stop(paste('new.t =', new.t, '<= t =', t))
    }
    else {
      dt <- new.t-t
    }

    if ( is.null(randomness) ) {
      randomness <- rnorm(dim.state*n.particles)
    }

    ensemble <<- ( ensemble +
                  drift(t, ensemble)*dt +
                  dyn.noise.amps*matrix(randomness,dim.state,n.particles)*sqrt(dt))
    t <<- new.t
  }

  ## Given a new observation, update the filter.  Note that
  ## this modifies the ensemble.
  update <- function(new.t, new.obs, randomness=NULL) {

    if ( new.t <= t ) {
      stop(paste('new.t =', new.t, '<= t =', t))
    }
    else {
      dt <- new.t-t
    }

    cond.mean <- ensemble + drift(t, ensemble)*dt
    cond.mean[1:dim.obs,] <- A*new.obs + B*cond.mean[1:dim.obs,]

    if ( is.null(randomness) ) {
      randomness <- rnorm(dim.state*n.particles)
    }

    ensemble <<- cond.mean + cond.scales*matrix(randomness,dim.state,n.particles)*sqrt(dt)
    t <<- new.t
  }

  ## Estimate the contribution of a given observation to the
  ## total log-likelihood, following the supplemental
  ## material for [Ionides et al PNAS 2006].  Note that we
  ## have to resample to get log-likelihood estimates
  ## because we sample exactly; in a bootstrap filter this
  ## can be done more or less for free.
  logprod <-
    if ( dim.obs > 1 ) {
      function(A) {
        log(mean(exp(colSums(log(A)))))
      }
    }
    else {
      function(A) {
        log(mean(A))
      }
    }
  
  loglike <- function(new.t, new.obs, randomness=NULL) {

    if ( new.t <= t ) {
      stop(paste('new.t =', new.t, '<= t =', t))
    }
    else {
      dt <- new.t-t
    }

    if ( is.null(randomness) ) {
      randomness <- rnorm(dim.state*n.particles)
    }

    ens <- ( ensemble + drift(t, ensemble)*dt
            + diffusion.scales * sqrt(dt) * matrix(randomness,
                                                   dim.state, n.particles))

    logprod(dnorm(new.obs, mean=ens[1:dim.obs,], sd=obs.noise.amps))
  }

  list(ensemble=function() {ensemble},
       time=function() {t},
       step=step,
       update=update,
       drift=drift,
       loglike=loglike,
       cond.scales=cond.scales,
       diffusion.scales=diffusion.scales,
       cond.mean.coeffs=c(A,B))
}



## Multiplicative noise version is slower, but interface is
## essentially the same.

make.pf.mult <- function(dim.state,  # state space dimension
                         dim.obs,    # observation space dimension
                         n.particles,  # num particles
                         drift0,     # drift=drift(state) or drift(t,state)
                         
                         noise.dyn,   # SDE diffusion coeff matrix
                                      # assumed diagonal ie a vector
                         
                         obs.noise.amps,         # observational noise
                         
                         init.ensemble=FALSE # optional: specify an initial ensemble
                                        # default: all particles start at 0
                         ) {

  ## Store ensemble in a matrix; each state vector is a
  ## column
  if ( is.matrix(init.ensemble)) {
    ensemble <- init.ensemble
  }
  else {
    ensemble <- matrix(0, dim.state, n.particles)
  }

  ## Check if drift function depends on time.
  if ( length(formals(drift0)) == 1 ) {
    drift <- function(t, ensemble) {
      drift0(ensemble)
    }
  }
  else {
    drift <- drift0
  }

  ## Track current time
  t <- 0

  ## Evolve the ensemble forward by one time step, without
  ## taking into account any observations.
  step <- function(new.t, randomness=NULL) {

    if ( new.t <= t ) {
      stop(paste('new.t =', new.t, '<= t =', t))
    }
    else {
      dt <- new.t-t
    }

    if ( is.null(randomness) ) {
      randomness <- rnorm(dim.state*n.particles)
    }

    ensemble <<- ( ensemble +
                  drift(t, ensemble)*dt +
                  noise.dyn(ensemble)*matrix(randomness,dim.state,n.particles)*sqrt(dt))
    t <<- new.t
  }

  ## Given a new observation, update the filter.  Note that
  ## this modifies the ensemble.
  update <- function(new.t, new.obs, randomness=NULL) {

    if ( new.t <= t ) {
      stop(paste('new.t =', new.t, '<= t =', t))
    }
    else {
      dt <- new.t-t
    }

    ## Pack up everything in a filter function whose job is
    ## given a new observation, update the ensemble and output
    ## the expected state of the system.

    diffusion.scales <- noise.dyn(ensemble)
    C <- diffusion.scales^2
    C.noise <- obs.noise.amps^2

    ## coefficients for conditional means
    A <- C[1:dim.obs] / (C[1:dim.obs] + C.noise)
    B <- C.noise / (C[1:dim.obs] + C.noise)

    ## compute conditional variances
    D <- vector(mode='numeric',dim.state)
    D[1:dim.obs] <- 1/(C[1:dim.obs] + C.noise)
    cond.vars <- C - C*D*C
    
    ## In the general case, this should be a matrix square root
    cond.scales <- sqrt(cond.vars)

    cond.mean <- ensemble + drift(t, ensemble)*dt
    cond.mean[1:dim.obs,] <- A*new.obs + B*cond.mean[1:dim.obs,]

    if ( is.null(randomness) ) {
      randomness <- rnorm(dim.state*n.particles)
    }

    ensemble <<- cond.mean + cond.scales*matrix(randomness,dim.state,n.particles)*sqrt(dt)
    t <<- new.t
  }

  ## Estimate the contribution of a given observation to the
  ## total log-likelihood, following the supplemental
  ## material for [Ionides et al PNAS 2006].  Note that we
  ## have to resample to get log-likelihood estimates
  ## because we sample exactly; in a bootstrap filter this
  ## can be done more or less for free.
  logprod <-
    if ( dim.obs > 1 ) {
      function(A) {
        log(mean(exp(colSums(log(A)))))
      }
    }
    else {
      function(A) {
        log(mean(A))
      }
    }
  
  loglike <- function(new.t, new.obs, randomness=NULL) {

    if ( new.t <= t ) {
      stop(paste('new.t =', new.t, '<= t =', t))
    }
    else {
      dt <- new.t-t
    }

    if ( is.null(randomness) ) {
      randomness <- rnorm(dim.state*n.particles)
    }

    ## Pack up everything in a filter function whose job is
    ## given a new observation, update the ensemble and output
    ## the expected state of the system.

    diffusion.scales <- noise.dyn(ensemble)

    ens <- ( ensemble + drift(t, ensemble)*dt
            + diffusion.scales * sqrt(dt) * matrix(randomness,
                                                   dim.state, n.particles))

    logprod(dnorm(new.obs, mean=ens[1:dim.obs,], sd=obs.noise.amps))
  }

  list(ensemble=function() {ensemble},
       time=function() {t},
       step=step,
       update=update,
       drift=drift,
       loglike=loglike)
}



################################
## Standard filter methods


## Evolve particle filter by one step
pf.step <- function(pf, t, ...) { pf$step(t, ...) }

## Update particle filter with given observation
pf.update <- function(pf, t, obs, ...) { pf$update(t, obs,...) }

## Estimate log-likelihood.  Note that if called
## sequentially, one should always call pf.loglike() before
## calling pf.update().
pf.loglike <- function(pf, t, obs, ...) { pf$loglike(t, obs,...) }

## Estimate mean state and corresponding variance.
pf.means <- function(pf) { pf.reduce(pf,mean) }
pf.vars <- function(pf) { pf.reduce(pf,var) }

## Useful combinator
pf.reduce <- function(pf,f) {
  ens <- pf$ensemble()
  sapply(1:nrow(ens), function(i) { f(ens[i,]) })
}

## Run particle filter PF with timestep DT on the given
## observation sequence; output is a sequence of states.

pf.run <- function(pf, dt, observations, obs.times=NULL, verbose=FALSE) {

  if ( is.null(obs.times)) {
    obs.times <- ((1:length(observations))-1)*dt
  }

  nsteps <- length(observations)-1
  output <- vector(mode='list', nsteps+1)
  output[[1]] <- c(pf.means(pf), pf.vars(pf), 0)

  for ( k in 1:nsteps ) {

    obs <- observations[[k+1]]
    curr.time <- obs.times[[k]]
    obs.time <- obs.times[[k+1]]

    if ( verbose ) {
      cat('# step=', k, '/', nsteps, 'obs=', obs, '\n')
    }

    ## evolve until just before next observation time,
    ## assuming uniform spacing of observation times
    l <- 1
    m <- round((obs.time - curr.time) / dt)
    while ( l < m ) {
      curr.time <- curr.time+dt
      pf.step(pf, curr.time)
      l <- l+1
    }
    
    logl <- pf.loglike(pf, obs.time, obs)  ## this must come before updating
    pf.update(pf, obs.time, obs)
    output[[k+1]] <- c(pf.means(pf), pf.vars(pf), logl)
  }

  output
}

## save a lot more info

pf.run.big <- function(pf, dt, observations, obs.times=NULL, verbose=FALSE) {

  if ( is.null(obs.times)) {
    obs.times <- ((1:length(observations))-1)*dt
  }

  nsteps <- length(observations)-1
  output <- vector(mode='list', round(max(obs.times)/dt))
  output[[1]] <- pf$ensemble()
  count <- 1

  for ( k in 1:nsteps ) {

    obs <- observations[[k+1]]
    curr.time <- obs.times[[k]]
    obs.time <- obs.times[[k+1]]

    if ( verbose ) {
      cat('# step=', k, '/', nsteps, 'obs=', obs, '\n')
    }

    ## evolve until just before next observation time,
    ## assuming uniform spacing of observation times
    l <- 1
    m <- round((obs.time - curr.time) / dt)
    while ( l < m ) {

      curr.time <- curr.time+dt
      pf.step(pf, curr.time)

      if ( verbose ) {
        cat('pf.run.big(): count=', count,
            'pf.time=', pf$time(),
            'curr.time=', curr.time, '\n')
      }

      count <- count+1
      output[[count]] <- pf$ensemble()
      l <- l+1
    }
    
    pf.update(pf, obs.time, obs)

    if ( verbose ) {
      cat('pf.run.big(): count=', count,
          'pf.time=', pf$time(), '\n')
    }

    count <- count+1
    output[[count]] <- pf$ensemble()
  }

  if ( verbose ) {
    cat('pf.run.big(): count=', count, '\n')
  }

  output
}


## Given two particle filters, compare the maximum distance
## between their respective ensembles.

pf.dist <- function(pf1,pf2) {
  ens1 <- pf1$ensemble()
  ens2 <- pf2$ensemble()
  if ( length(ens1) != length(ens2)) {
    warning('ensembles do not have same size')
    Inf
  }
  else {
    max(abs(ens1-ens2))
 }
}
