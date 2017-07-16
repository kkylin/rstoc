## paths.r

## This file is part of RSTOC, a system for stochastic
## optimal control and experimental design written in R.

## Copyright (C) 2012,2013,2014,2015,2016,2017 by the authors:

## Giles Hooker <gjh27@cornell.edu>
## Kevin K Lin <klin@math.arizona.edu>
## Bruce Rogers <bruce.warren.rogers@gmail.com>

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

## Generate a bunch of paths, using the given path generator.

## make.paths <- function(npaths, make.path, ...) {
##   paths <- vector(mode="list", npaths)
##   for ( n in 1:npaths ) {
##     paths[[n]] <- make.path(...)
##   }
##   paths
## }

make.sample.path <- function(drift.vector, nsteps, dt, noise.dyn,
                             control.value=0,
                             start=NULL,
                             include.control=FALSE,
                             ...) {

  if (is.null(start)) {
    start <- 0*noise.dyn
  }
  
  path <- vector(mode="list", nsteps+1)
  sqrt.dt <- sqrt(dt)
  ndim <- length(noise.dyn)

  if ( include.control ) {

    path[[1]] <- c(start, control.value)

    for ( step in 1:nsteps ) {

      curr.state <- path[[step]][1:ndim]
      
      path[[step+1]] <- c(curr.state +
                          dt * drift.vector(control.value, curr.state, ...) +
                          sqrt.dt * rnorm(ndim, mean=0, sd=noise.dyn),
                          control.value)
    }
  }
  else {

    path[[1]] <- start

    for ( step in 1:nsteps ) {
      path[[step+1]] <- (path[[step]] +
                         dt * drift.vector(control.value, path[[step]], ...) +
                         sqrt.dt * rnorm(ndim, mean=0, sd=noise.dyn))
    }
  }
  
  path
}


## multiplicativ noise

make.sample.path.mult <- function(ndim, drift.vector, nsteps, dt, noise.dyn,
                                  control.value=0,
                                  start=NULL,
                                  include.control=FALSE,
                                  ...) {

  if (is.null(start)) {
    start <- vector(mode='numeric', ndim)
  }
  
  path <- vector(mode="list", nsteps+1)
  sqrt.dt <- sqrt(dt)

  if ( include.control ) {

    path[[1]] <- c(start, control.value)

    for ( step in 1:nsteps ) {

      curr.state <- path[[step]][1:ndim]
      
      path[[step+1]] <- c(curr.state +
                          dt * drift.vector(control.value, curr.state, ...) +
                          sqrt.dt * rnorm(ndim,
                                          mean=0,
                                          sd=noise.dyn(control.value, curr.state, ...)),
                          control.value)
    }
  }
  else {

    path[[1]] <- start

    for ( step in 1:nsteps ) {
      path[[step+1]] <- (path[[step]] +
                         dt * drift.vector(control.value, path[[step]], ...) +
                         sqrt.dt * rnorm(ndim,
                                         mean=0,
                                         sd=noise.dyn(control.value, path[[step]], ...)))
    }
  }
  
  path
}


## Generate sample paths with given control.  The `control'
## object should be the output of bellman().

make.controlled.path <- function(drift.vector, nsteps, dt, control,
                                 noise.dyn, T0=0, start=NULL, skip=0,
                                 ...) {

  if (is.null(start)) {
    start <- 0*noise.dyn
  }
  
  path <- vector(mode="list", nsteps+1)
  path[[1]] <- c(start,0)
  sqrt.dt <- sqrt(dt)

  if ( length(noise.dyn) == 2 ) {
    for ( step in 1:nsteps ) {

      z <- path[[step]]

      if ( (step-1)%%(skip+1) == 0 ) {
        u <- bellman2d.control.value(control, T0+(step-1)*dt, z[[1]], z[[2]])
      }

      path[[step]][[3]] <- u
      path[[step+1]] <- c(path[[step]][1:2] +
                          dt * drift.vector(u, path[[step]][1:2], ...) +
                          sqrt.dt * noise.dyn * rnorm(2),
                          0)
    }
  }
  else if ( length(noise.dyn) == 1 ) {
    for ( step in 1:nsteps ) {

      z <- path[[step]]
      
      if ( (step-1)%%(skip+1) == 0 ) {
        u <- bellman1d.control.value(control, T0+(step-1)*dt, z[[1]])
      }
      
      path[[step]][[2]] <- u
      path[[step+1]] <- c(path[[step]][1] +
                          dt * drift.vector(u, path[[step]][1], ...) +
                          sqrt.dt * noise.dyn * rnorm(1),
                          0)
    }
  }
  else {
    stop('# make.controlled.path() only supports dimensions 1 and 2')
  }


  path
}

make.controlled.path.mult <- function(ndim, drift.vector, nsteps, dt, control,
                                      noise.dyn, T0=0, start=NULL, skip=0,
                                      ...) {

  if (is.null(start)) {
    start <- vector(mode='numeric', ndim)
  }
  
  path <- vector(mode="list", nsteps+1)
  path[[1]] <- c(start,0)
  sqrt.dt <- sqrt(dt)

  if ( ndim == 2 ) {
    for ( step in 1:nsteps ) {

      z <- path[[step]]

      if ( (step-1)%%(skip+1) == 0 ) {
        u <- bellman2d.control.value(control, T0+(step-1)*dt, z[[1]], z[[2]])
      }

      path[[step]][[3]] <- u
      path[[step+1]] <- c(path[[step]][1:2] +
                          dt * drift.vector(u, path[[step]][1:2], ...) +
                          sqrt.dt * noise.dyn(u, path[[step]][1:2], ...) * rnorm(2),
                          0)
    }
  }
  else if ( ndim == 1 ) {
    for ( step in 1:nsteps ) {

      z <- path[[step]]
      
      if ( (step-1)%%(skip+1) == 0 ) {
        u <- bellman1d.control.value(control, T0+(step-1)*dt, z[[1]])
      }
      
      path[[step]][[2]] <- u
      path[[step+1]] <- c(path[[step]][1] +
                          dt * drift.vector(u, path[[step]][1], ...) +
                          sqrt.dt * noise.dyn(u, path[[step]][1], ...) * rnorm(1),
                          0)
    }
  }
  else {
    stop('# make.controlled.path() only supports dimensions 1 and 2')
  }


  path
}


## Plot path

plot.path <- function(path) {

  x <- sapply(path, function(x) { x[[1]] })
  y <- sapply(path, function(x) { x[[2]] })
  u <- sapply(path, function(x) { x[[3]] })

  ## par(new=FALSE, yaxp=c(-4,4,2))
  ## plot(x, type='l')

  ## par(new=TRUE, yaxp=c(-4,4,2))
  ## plot(y, type='l', col='blue')

  ## par(new=TRUE, yaxp=c(-4,4,2))
  ## plot(u, type='l', col='red')
  plot((x-y)*u, type='l', col='blue')
}  

pathref <- function(path,i,T=FALSE) {
  if ( i > 0 ) {
    sapply(path, function(x) {x[i]})
  }
  else {
    if ( !T ) {
      T <- length(path)
    }

    dt <- T / (length(path)-1)
    
    ( 1:length(path) - 1 ) * dt
  }
}

plot.traj <- function(path, i, j=FALSE, ...) {
  if ( j ) {
    plot.traj2(path, i, j, ...)
  }
  else {
    plot.traj1(path, i, ...)
  }
}

plot.traj1 <- function(path,i,...) {

  if ( length(path[[1]]) > 2 ) {
    cstat <- 'controlled'
  }
  else {
    cstat <- 'uncontrolled'
  }
  
  plot(pathref(path,0,...), pathref(path,i,...),
       xlab='t', ylab=paste('x[', i, ']'))
}

plot.traj2 <- function(path, i, j, ...) {

  if ( length(path[[1]]) > 2 ) {
    cstat <- 'controlled'
  }
  else {
    cstat <- 'uncontrolled'
  }
  
  plot(pathref(path,i), pathref(path,j),
       xlab=paste('x[', i, ']'),
       ylab=paste('x[', j, ']'),
       ...)
}
