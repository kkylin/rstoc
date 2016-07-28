## chemo.r

## This file is part of RSTOC, a system for stochastic
## optimal control and experimental design written in R.

## Copyright (C) 2012,2013,2014,2015,2016 by the authors:

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

## chemostat model

## default parameters
default.chemo.params <- data.frame(NI = 160, rho = 270, X = 0.0027, K = 4.4)


## With default parameters, there is a stable fixed point:
chemo.fp <- function(delta=0.3,
                     K   = default.chemo.params$K,
                     NI  = default.chemo.params$NI,
                     rho = default.chemo.params$rho,
                     X   = default.chemo.params$X) {

  exp.n.fix <- K*delta/(X*rho-delta)
  exp.c.fix <- delta*(NI-exp.n.fix)*(K+exp.n.fix)/exp.n.fix/rho
  log(c(exp.c.fix, exp.n.fix))
}

cat('## fixed point at', chemo.fp(0.3), '\n')


chemo.nullclines <- function(delta,
                             x0=-5, x1=0, nx=100,
                             y0=-2, y1=6, ny=100,
                             K   = default.chemo.params$K,
                             NI  = default.chemo.params$NI,
                             rho = default.chemo.params$rho,
                             X   = default.chemo.params$X) {

  # c'=0
  dx=(x1-x0)/(nx-1)
  cl=seq(x0, x1, dx)
  c.nullcline = lapply(cl, function(x) { c(x, log(delta*K/(X*rho-delta))) })
  
  # n'=0
  dy=(min(y1,log(NI))-y0)/(ny-1)
  nl <- seq(y0, min(y1,log(NI)), dy)
  n.nullcline <- lapply(nl,
                        function(n) {
                          c(log(delta*(NI-exp(n))*(K+exp(n))*exp(-n)/rho), n)
                        })

  list(c.null=c.nullcline,
       n.null=n.nullcline)
}



## ML vector field
chemo.drift.c <- function(delt, c, n,
                          K   = default.chemo.params$K,
                          NI  = default.chemo.params$NI,
                          rho = default.chemo.params$rho,
                          X   = default.chemo.params$X) {
  
  X*rho*exp(n)/(K+exp(n)) - delt
}

chemo.drift.n <- function(delt, c, n,
                          K   = default.chemo.params$K,
                          NI  = default.chemo.params$NI,
                          rho = default.chemo.params$rho,
                          X   = default.chemo.params$X) {
  
  delt*NI*exp(-n) - rho*exp(c)/(K+exp(n)) - delt
}


################################
## Drift

## x and y are assumed to be row vectors

chemo.drift.vector <- function(u, state, ...) {

  num.states <- length(state)
  c <- state[seq(1,num.states,2)]
  n <- state[seq(2,num.states,2)]
  cdot <- chemo.drift.c(u, c, n, ...)
  ndot <- chemo.drift.n(u, c, n, ...)
  
  if ( is.matrix(state)) {
    matrix(c(cdot,ndot), 2, ncol(state), byrow=T)
  }
  else {
    c(cdot,ndot)
  }
}


################################
## Fisher informations

## Given initial condition (X Y), compute the Fisher
## information for the parameter K, using the value U as
## control over the interval [0,DT].

## Note minus sign: bellman() minimizes cost by default,
## whereas we want to maximize Fisher information.



## note negative signs
chemo.fisher.c <- function(delt, c, n,
                           K   = default.chemo.params$K,
                           NI  = default.chemo.params$NI,
                           rho = default.chemo.params$rho,
                           X   = default.chemo.params$X) {
  
  -( X*rho*exp(n)/(K+exp(n))^2 )^2
}

chemo.fisher.n <- function(delt, c, n,
                           K   = default.chemo.params$K,
                           NI  = default.chemo.params$NI,
                           rho = default.chemo.params$rho,
                           X   = default.chemo.params$X) {

  -( rho*exp(c)/(K+exp(n))^2 )^2
}

chemo.fisher.K <- function(delt, c, n, ...) {
  chemo.fisher.n(delt,c,n,...) + chemo.fisher.c(delt,c,n,...)
}

chemo.fisher <- chemo.fisher.K


## FI for the parameter rho
chemo.fisher.rho <- function(delta, c, n,
                             K   = default.chemo.params$K,
                             NI  = default.chemo.params$NI,
                             rho = default.chemo.params$rho,
                             X   = default.chemo.params$X) {
  -((X*exp(n)/(K+exp(n)))^2 + (exp(c)/(K+exp(n)))^2)
}


## FI for the parameter X
chemo.fisher.X <- function(delta, c, n,
                           K   = default.chemo.params$K,
                           NI  = default.chemo.params$NI,
                           rho = default.chemo.params$rho,
                           X   = default.chemo.params$X) {
  -(rho*exp(n)/(K+exp(n)))^2
}

## FI for the parameter NI
chemo.fisher.NI <- function(delta, c, n,
                            K   = default.chemo.params$K,
                            NI  = default.chemo.params$NI,
                            rho = default.chemo.params$rho,
                            X   = default.chemo.params$X) {
  -(delta*exp(-n))^2
}


################################
## Alternative objection functions

## Sensitivity of fixed point to perturbations

chemo.sens <- function(delta, c, n,
                       K   = default.chemo.params$K,
                       NI  = default.chemo.params$NI,
                       rho = default.chemo.params$rho,
                       X   = default.chemo.params$X) {

  a <- matrix(c(0,                            ## dc'/dc
                -rho*exp(c)/(K+exp(n)),       ## dn'/dc
                X*rho*exp(n)/(K+exp(n))       ## dc'/dn
                -X*rho*exp(2*n)/(K+exp(n))^2,
                -delta*NI*exp(-n)             ## dn'/dn
                +rho*exp(c+n)/(K+exp(n))^2),
              2,2)
  
  b <- c(-X*rho*exp(n)/(K+exp(n))^2, rho*exp(c)/(K+exp(n))^2)

  -solve(a,b)
  ## -((solve(a,b))[[1]])^2
}

chemo.fp.sens <- function(delta,
                          K   = default.chemo.params$K,
                          NI  = default.chemo.params$NI,
                          rho = default.chemo.params$rho,
                          X   = default.chemo.params$X) {
  
  fp <- chemo.fp(delta, K, NI, rho, X)
  chemo.sens(delta, fp[[1]], fp[[2]], K, NI, rho, X)
}
