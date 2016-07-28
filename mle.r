## mle.r

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

## Here we estimate log-likelihoods, and use them to do MLE.

## Notes:

## Derivative-based maximization can be made to work using a
## good ensemble coupling.



################################
## Rewrite of pf.loglikelihoods().  Starting at t=0, does
## the following:

## 1: Generate a short sample path segment, up to the next
## observation time,

## 2: Run ensemble of particle filters on next observation,
## then average the state estimates over all priors

## 3: Set new control value using state estimate

## 4: Update log-likelihood estimates for all parameters

## 5: Goto 1

## Outputs:

## resulting sample path + control values

## log-likelihood values

pf.expt <- function(dim.state, dim.obs, n.particles,
                    bc, drift,
                    noise.dyn, noise.obs,
                    params,
                    start=NULL,
                    avg.ctl=FALSE,
                    verbose=TRUE,
                    ...) {

  ## default init state is 0
  if (is.null(start)) {
    start <- 0*(1:dim.state)
  }

  ## echo useful msgs to stdout
  echo <- function(...) {
    if ( verbose ) {
      cat('## pf.expt():', ..., '\n')
    }
  }

  ## figure out runtime params
  echo('n.params', n.params <- length(params))
  echo('n.steps', n.steps <- bc$nsteps)
  echo('dt', dt <- bc$dt)
  echo('n.obs', n.obs <- bc$nsteps.actual)
  echo('obs.dt', obs.dt <- bc$dt.actual)
  echo('nsubsteps', nsubsteps <- ceiling(obs.dt / dt))
  echo('new dt', dt <- obs.dt / nsubsteps)

  ## initialize control + particle filters
  if ( dim.state == 2 ) {
    ctl <- bellman2d.control.value(bc, 0, start[[1]], start[[2]])
  }
  else if ( dim.state == 1 ) {
    ctl <- bellman1d.control.value(bc, 0, start)
  }
  else {
    stop('## pf.expt() can only do 1 & 2 dims')
  }
  
  pf.list <- lapply(1:n.params,
                    function(pptr) {
                      
                      ## force evaluation
                      param <- params[[pptr]]
                      echo('param(', pptr, ')=', param)

                      make.pf(dim.state, dim.obs, n.particles,
                              function(t, state) {
                                drift(ctl, state, param)
                              },
                              noise.dyn, noise.obs,
                              init.ensemble=matrix(start, dim.state, nparticles),
                              ...)
                    })

  ## storage for output data
  loglikes <- matrix(0, n.params, n.obs)
  path <- list(c(start,ctl))
  obs.list <- list(start[1:dim.obs])
  est.list <- list(start)
  var.list <- list(0*start)
  obs.times <- list(0)

  ## run
  for ( k in 1:n.obs ) {

    ## generate sample path up to the next observation time
    curr.time <- (k-1)*obs.dt
    new.obs.time <- obs.times[[k+1]] <- k*obs.dt
    curr.state <- path[[length(path)]]

    path.seg <- make.sample.path(drift, nsubsteps, dt, noise.dyn,
                                 control.value=ctl,
                                 start=curr.state,
                                 include.control=TRUE)

    path <- c(path, path.seg[2:(nsubsteps+1)])

    ## make a noisy observation
    new.obs <- path.seg[[nsubsteps+1]][1:dim.obs] + rnorm(dim.obs,sd=noise.obs)
    echo('obs', k, 'at', new.obs.time, 'new.obs', new.obs, 'ctl', ctl)
    echo('new path length', length(path))

    ## evolve particle filters up to just before the next
    ## observation time
    if ( nsubsteps>1 ) {
      for ( l in 1:(nsubsteps-1)) {
        w <- rnorm(dim.state * n.particles)
        for ( pf in pf.list ) {
          pf.step(pf, curr.time + l*dt, randomness=w)
        }
      }
    }

    ## grab new observation and compute log-likelihoods
    obs.list[[k+1]] <- new.obs
    w <- rnorm(dim.state * n.particles)
    
    for ( pfptr in 1:length(pf.list) ) {
      pf <- pf.list[[pfptr]]

      ## for reasons having to do with how particle filters
      ## are implemented, this should come before pf.update()
      loglikes[pfptr, k] <- pf.loglike(pf, new.obs.time, new.obs,
                                       randomness=w)
    }

    ## update particle filters and control using new obs
    for ( pf in pf.list ) {
      pf.update(pf, new.obs.time, new.obs, randomness=w)
    }

    new.state.ests <- lapply(pf.list, function(pf) { pf.means(pf) })
    new.state.est <- sapply(1:dim.state,
                            function(i) {
                              mean(sapply(new.state.ests, function(x) { x[[i]] }))
                            })
    est.list[[k+1]] <- new.state.est
    var.list[[k+1]] <- mean(sapply(pf.list, function(pf) { pf.vars(pf) }))

    if ( avg.ctl ) {
      ## compute the average control value over all
      ## estimated states
      if ( dim.state == 2 ) {
        ctl <- mean(sapply(pf.list,
                           function(pf) {
                             ens <- pf$ensemble()
                             mean(bellman2d.control.value(bc, new.obs.time, ens[1,], ens[2,]))
                           }))
      }
      else {
        ctl <- mean(sapply(pf.list,
                           function(pf) {
                             ens <- pf$ensemble()
                             mean(bellman1d.control.value(bc, new.obs.time, ens[1,]))
                           }))
      }
    }
    else {
      ## compute the control value of the average estimated
      ## state
      if ( dim.state == 2 ) {
        ctl <- bellman2d.control.value(bc, new.obs.time,
                                       new.state.est[[1]], new.state.est[[2]])
      }
      else {
        ctl <- bellman1d.control.value(bc, new.obs.time, new.state.est[[1]])
      }
    }
  }

  list(loglikes=loglikes,
       loglike=rowMeans(loglikes),
       path=path,
       times=(0:n.steps)*dt,
       obs=obs.list,
       ests=est.list,
       vars=var.list,
       obs.times=obs.times,
       pf.list=pf.list)
}


## Multiplicative noise version
pf.expt.mult <- function(dim.state, dim.obs, n.particles,
                         bc, drift,
                         noise.dyn, noise.obs,
                         params,
                         start=NULL,
                         avg.ctl=FALSE,
                         verbose=TRUE,
                         ...) {

  ## default init state is 0
  if (is.null(start)) {
    start <- 0*(1:dim.state)
  }

  ## echo useful msgs to stdout
  echo <- function(...) {
    if ( verbose ) {
      cat('## pf.expt.mult():', ..., '\n')
    }
  }

  ## figure out runtime params
  echo('n.params', n.params <- length(params))
  echo('n.steps', n.steps <- bc$nsteps)
  echo('dt', dt <- bc$dt)
  echo('n.obs', n.obs <- bc$nsteps.actual)
  echo('obs.dt', obs.dt <- bc$dt.actual)
  echo('nsubsteps', nsubsteps <- ceiling(obs.dt / dt))
  echo('new dt', dt <- obs.dt / nsubsteps)

  ## initialize control + particle filters
  if ( dim.state == 2 ) {
    ctl <- bellman2d.control.value(bc, 0, start[[1]], start[[2]])
  }
  else if ( dim.state == 1 ) {
    ctl <- bellman1d.control.value(bc, 0, start)
  }
  else {
    stop('## pf.expt.mult() can only do 1 & 2 dims')
  }
  
  pf.list <- lapply(1:n.params,
                    function(pptr) {
                      
                      ## force evaluation
                      param <- params[[pptr]]
                      echo('param(', pptr, ')=', param)

                      make.pf.mult(dim.state, dim.obs, n.particles,
                                   function(t, state) {
                                     drift(ctl, state, param)
                                   },
                                   noise.dyn, noise.obs,
                                   init.ensemble=matrix(start, dim.state, nparticles),
                                   ...)
                    })

  ## storage for output data
  loglikes <- matrix(0, n.params, n.obs)
  path <- list(c(start,ctl))
  obs.list <- list(start[1:dim.obs])
  est.list <- list(start)
  var.list <- list(0*start)
  obs.times <- list(0)

  ## run
  for ( k in 1:n.obs ) {

    ## generate sample path up to the next observation time
    curr.time <- (k-1)*obs.dt
    new.obs.time <- obs.times[[k+1]] <- k*obs.dt
    curr.state <- path[[length(path)]]

    path.seg <- make.sample.path.mult(dim.state, drift, nsubsteps, dt,
                                      function(control, state) {
                                        noise.dyn(state)
                                      },
                                      control.value=ctl,
                                      start=curr.state,
                                      include.control=TRUE)

    path <- c(path, path.seg[2:(nsubsteps+1)])

    ## make a noisy observation
    new.obs <- path.seg[[nsubsteps+1]][1:dim.obs] + rnorm(dim.obs,sd=noise.obs)
    echo('obs', k, 'at', new.obs.time, 'new.obs', new.obs, 'ctl', ctl)
    echo('new path length', length(path))

    ## evolve particle filters up to just before the next
    ## observation time
    if ( nsubsteps>1 ) {
      for ( l in 1:(nsubsteps-1)) {
        w <- rnorm(dim.state * n.particles)
        for ( pf in pf.list ) {
          pf.step(pf, curr.time + l*dt, randomness=w)
        }
      }
    }

    ## grab new observation and compute log-likelihoods
    obs.list[[k+1]] <- new.obs
    w <- rnorm(dim.state * n.particles)
    
    for ( pfptr in 1:length(pf.list) ) {
      pf <- pf.list[[pfptr]]

      ## for reasons having to do with how particle filters
      ## are implemented, this should come before pf.update()
      loglikes[pfptr, k] <- pf.loglike(pf, new.obs.time, new.obs,
                                       randomness=w)
    }

    ## update particle filters and control using new obs
    for ( pf in pf.list ) {
      pf.update(pf, new.obs.time, new.obs, randomness=w)
    }

    new.state.ests <- lapply(pf.list, function(pf) { pf.means(pf) })
    new.state.est <- sapply(1:dim.state,
                            function(i) {
                              mean(sapply(new.state.ests, function(x) { x[[i]] }))
                            })
    est.list[[k+1]] <- new.state.est
    var.list[[k+1]] <- mean(sapply(pf.list, function(pf) { pf.vars(pf) }))

    if ( avg.ctl ) {
      ## compute the average control value over all
      ## estimated states
      if ( dim.state == 2 ) {
        ctl <- mean(sapply(pf.list,
                           function(pf) {
                             ens <- pf$ensemble()
                             mean(bellman2d.control.value(bc, new.obs.time, ens[1,], ens[2,]))
                           }))
      }
      else {
        ctl <- mean(sapply(pf.list,
                           function(pf) {
                             ens <- pf$ensemble()
                             mean(bellman1d.control.value(bc, new.obs.time, ens[1,]))
                           }))
      }
    }
    else {
      ## compute the control value of the average estimated
      ## state
      if ( dim.state == 2 ) {
        ctl <- bellman2d.control.value(bc, new.obs.time,
                                       new.state.est[[1]], new.state.est[[2]])
      }
      else {
        ctl <- bellman1d.control.value(bc, new.obs.time, new.state.est[[1]])
      }
    }
  }

  list(loglikes=loglikes,
       loglike=rowMeans(loglikes),
       path=path,
       times=(0:n.steps)*dt,
       obs=obs.list,
       ests=est.list,
       vars=var.list,
       obs.times=obs.times,
       pf.list=pf.list)
}



################################
## Use particle filter to compute log-likelihoods of
## observations.

## Input: sequence of observations

## Output: parameter estimate

## Derivative-based maximization can be made to work using a
## good ensemble coupling.

## Will need to do some experiments to figure this out...

pf.loglikelihoods <- function(dim.state, dim.obs, n.particles,
                              param.drift, dt,
                              dyn.noise.amps, obs.noise.amps,
                              params,
                              observations,
                              obs.times=NULL,
                              verbose=FALSE,
                              ...) {

  if ( is.null(obs.times)) {
    ## observations DO NOT include observation times (old style)
    pf.loglikelihoods.old(dim.state, dim.obs, n.particles,
                          param.drift, dt,
                          dyn.noise.amps, obs.noise.amps,
                          params,
                          observations,
                          verbose=verbose,
                          ...)
  }
  else {
    ## observation times are provided
    pf.loglikelihoods.new(dim.state, dim.obs, n.particles,
                          param.drift, dt,
                          dyn.noise.amps, obs.noise.amps,
                          params,
                          observations, obs.times,
                          verbose=verbose,
                          ...)
  }    
}

pf.loglikelihoods.new <- function(dim.state, dim.obs, n.particles,
                                  param.drift, dt,
                                  dyn.noise.amps, obs.noise.amps,
                                  params,
                                  observations, obs.times,
                                  verbose=FALSE,
                                  ...) {

  nsteps <- length(observations)-1
  n.params <- length(params)

  if ( verbose ) {
    cat('nsteps=',nsteps,'|obs|=', length(observations), length(obs.times), '\n')
  }
  
  pf.list <- lapply(params,
                    function(param) {
                      param <- param   ## force evaluation
                      make.pf(dim.state, dim.obs, n.particles,
                              function(t, state) {
                                param.drift(t, state, param)
                              },
                              dyn.noise.amps, obs.noise.amps,
                              ...)
                    })

  output <- matrix(0, n.params, nsteps+1)

  for ( k in 1:nsteps ) {

    obs <- observations[[k+1]]
    curr.time <- obs.times[[k]]
    obs.time <- obs.times[[k+1]]

    ## evolve until just before next observation time,
    ## assuming uniform spacing of observation times
    l <- 1
    m <- round((obs.time - curr.time) / dt)
    while ( l < m ) {

      w <- rnorm(dim.state*n.particles) ## standard same-noise coupling
      curr.time <- curr.time+dt

      for ( pf in pf.list ) {
        pf.step(pf, curr.time, randomness=w)
      }

      l <- l+1
    }

    ## now handle the observation
    w <- rnorm(dim.state*n.particles) ## standard same-noise coupling

    if ( verbose ) {
      cat('k=', k, 'obs.time=',obs.time,'\n')
    }

    for ( i in 1:n.params ) {
      pf <- pf.list[[i]]
      output[i,k+1] <- pf.loglike(pf, obs.time, obs, w)  ## this must come before updating
      pf.update(pf, obs.time, obs, w)
    }

    if ( verbose ) {
      cat('#', 'step=', k, '/', nsteps, 'obs=', obs, 'loglike-range=',
          range(output[,k]), '\n')
    }
  }

  list(loglike=output, pf.list=pf.list)
}

pf.loglikelihoods.old <- function(dim.state, dim.obs, n.particles,
                                  param.drift, dt,
                                  dyn.noise.amps, obs.noise.amps,
                                  params,
                                  observations,
                                  verbose=FALSE,
                                  ...) {

  nsteps <- length(observations)-1
  n.params <- length(params)
  
  pf.list <- lapply(params,
                    function(param) {
                      param <- param   ## force evaluation
                      make.pf(dim.state, dim.obs, n.particles,
                              function(t, state) {
                                param.drift(t, state, param)
                              },
                              dyn.noise.amps, obs.noise.amps,
                              ...)
                    })

  output <- matrix(0, n.params, nsteps+1)

  for ( k in 1:nsteps ) {

    obs <- observations[[k+1]]
    w <- rnorm(dim.state*n.particles) ## standard same-noise coupling

    for ( i in 1:n.params ) {
      pf <- pf.list[[i]]
      output[i,k+1] <- pf.loglike(pf, k*dt, obs, w)  ## this must come before updating
      pf.update(pf, k*dt, obs, w)
    }

    if ( verbose ) {
      cat('#', 'step=', k, '/', nsteps, 'obs=', obs, 'loglike-range=',
          range(output[,k]), '\n')
    }
  }

  list(loglike=output, pf.list=pf.list)
}


## Multiplicative noise

pf.loglikelihoods.mult <- function(dim.state, dim.obs, n.particles,
                                   param.drift, dt,
                                   noise.dyn, noise.obs,
                                   params,
                                   observations,
                                   verbose=FALSE,
                                   ...) {

  nsteps <- length(observations)-1
  n.params <- length(params)
  
  pf.list <- lapply(params,
                    function(param) {
                      param <- param   ## force evaluation
                      make.pf.mult(dim.state, dim.obs, n.particles,
                                   function(t, state) {
                                     param.drift(t, state, param)
                                   },
                                   noise.dyn, noise.obs,
                                   ...)
                    })

  output <- matrix(0, n.params, nsteps+1)

  for ( k in 1:nsteps ) {

    obs <- observations[[k+1]]
    w <- rnorm(dim.state*n.particles) ## standard same-noise coupling

    for ( i in 1:n.params ) {
      pf <- pf.list[[i]]
      output[i,k+1] <- pf.loglike(pf, k*dt, obs, w)  ## this must come before updating
      pf.update(pf, k*dt, obs, w)
    }

    if ( verbose ) {
      cat('#', 'step=', k, '/', nsteps, 'obs=', obs, 'loglike-range=',
          range(output[,k]), '\n')
    }
  }

  list(loglike=output, pf.list=pf.list)
}


################################
## Compute log-likelihoods of full-state perfect
## observations.

## Input: sequence of observations

## Output: parameter estimate

sde.loglikelihoods <- function(dim.state,
                               param.drift, dt,
                               dyn.noise.amps,
                               params,
                               observations,
                               verbose=FALSE) {

  nsteps <- length(observations)-1
  n.params <- length(params)
  
  output <- matrix(0, n.params, nsteps+1)

  for ( k in 1:nsteps ) {
    
    t <- (k-1)*dt
    
    for ( i in 1:n.params ) {

      f <- param.drift(t, observations[[k]], params[[i]])
      r <- observations[[k+1]] - observations[[k]] - f*dt
      
      output[i,k+1] <- ( sum((r / dyn.noise.amps)^2)
                        +log(2*pi*prod(dyn.noise.amps)^2 ))
    }

    if ( verbose ) {
      cat('#', 'step=', k, '/', nsteps, 'obs=', obs, 'loglike-range=',
          range(output[,k]), '\n')
    }
  }

  output <- -0.5*output
  list(loglike=output)
}

sde.loglikelihoods.mult <- function(dim.state,
                                    param.drift, dt,
                                    noise.dyn,
                                    params,
                                    observations,
                                    verbose=FALSE) {

  nsteps <- length(observations)-1
  n.params <- length(params)
  
  output <- matrix(0, n.params, nsteps+1)

  for ( k in 1:nsteps ) {
    
    t <- (k-1)*dt
    
    for ( i in 1:n.params ) {

      f <- param.drift(t, observations[[k]], params[[i]])
      r <- observations[[k+1]] - observations[[k]] - f*dt

      sigma <- noise.dyn(observations[[k]])
      
      output[i,k+1] <- sum((r / sigma)^2) + log(2*pi*prod(sigma)^2)
    }

    if ( verbose ) {
      cat('#', 'step=', k, '/', nsteps, 'obs=', obs, 'loglike-range=',
          range(output[,k]), '\n')
    }
  }

  output <- -0.5*output
  list(loglike=output)
}


################################
## Naive MLE: estimate log-likelihoods on a grid of
## parameter values, then interpolate to find max.  See
## comments below for output format.

pf.mle <- function(dim.state, dim.obs, n.particles,
                   param.drift, dt,
                   dyn.noise.amps, obs.noise.amps,
                   params,
                   list.of.observations,
                   verbose=FALSE) {

  if ( max(abs(obs.noise.amps)) == 0 ) {
    mle.perfect(dim.state,
                param.drift, dt,
                dyn.noise.amps,
                params,
                list.of.observations,
                verbose)
  }
  else {
    pf.mle.noisy(dim.state, dim.obs, n.particles,
                 param.drift, dt,
                 dyn.noise.amps, obs.noise.amps,
                 params,
                 list.of.observations,
                 verbose)
  }
}

mle.perfect <- function(dim.state,
                        param.drift, dt,
                        dyn.noise.amps,
                        params,
                        list.of.observations,
                        verbose=FALSE) {

  n.params <- length(params)
  nsteps <- length(list.of.observations)
  llavgs <- vector(mode='numeric', n.params)

  lapply(list.of.observations,
         function(observations) {

           print(system.time(
                             ll <-
                             sde.loglikelihoods(dim.state,
                                                param.drift, dt,
                                                dyn.noise.amps,
                                                params,
                                                observations,
                                                verbose=verbose)))

           ## Using rowMeans() -- instead of rowSums() --
           ## effectively time-averages log-likelihoods
           ## along each sample path; this keeps the
           ## magnitude of log-likelihoods under control.
           llavgs <<- llavgs + rowMeans(ll$loglike) / nsteps
         })

  ## Find best fit of form A + B*theta + C*theta^2
  qfit <- lm( llavgs ~ params + I(params^2))
  A <- qfit$coefficients[[1]]
  B <- qfit$coefficients[[2]]
  C <- qfit$coefficients[[3]]

  list(params=params,    ## Grid of parameters used
       loglikes=llavgs,  ## Log-likelihood of each parameter 
       A=A,              ## Coefficients of interpolating polynomial
       B=B,              ## which has form AX^2 + BX + C
       C=C,
       mle=-B/2.0/C)     ## MLE
}

pf.mle.noisy <- function(dim.state, dim.obs, n.particles,
                         param.drift, dt,
                         dyn.noise.amps, obs.noise.amps,
                         params,
                         list.of.observations,
                         verbose=FALSE) {

  n.params <- length(params)
  nsteps <- length(list.of.observations)
  llavgs <- vector(mode='numeric', n.params)

  lapply(list.of.observations,
         function(observations) {

           print(system.time(
                             ll <-
                             pf.loglikelihoods(dim.state, dim.obs, n.particles,
                                               param.drift, dt,
                                               dyn.noise.amps, obs.noise.amps,
                                               params,
                                               observations,
                                               verbose=verbose)))

           ## Using rowMeans() -- instead of rowSums() --
           ## effectively time-averages log-likelihoods
           ## along each sample path; this keeps the
           ## magnitude of log-likelihoods under control.
           llavgs <<- llavgs + rowMeans(ll$loglike) / nsteps
         })

  ## Find best fit of form A + B*theta + C*theta^2
  qfit <- lm( llavgs ~ params + I(params^2))
  A <- qfit$coefficients[[1]]
  B <- qfit$coefficients[[2]]
  C <- qfit$coefficients[[3]]

  list(params=params,    ## Grid of parameters used
       loglikes=llavgs,  ## Log-likelihood of each parameter 
       A=A,              ## Coefficients of interpolating polynomial
       B=B,              ## which has form AX^2 + BX + C
       C=C,
       mle=-B/2.0/C)     ## MLE
}
