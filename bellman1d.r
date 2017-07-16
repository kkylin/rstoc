## bellman1d.r

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

## OVERVIEW

## Given a controlled SDE of the form

##   dX_t = F(u_t,X_t) dt + b_x dW^1_t

## and a cost functional of the form

##   COST[u] := \int_0^T C(u_t, X_t) dt

## we compute the control that minimizes COST[u].  This is
## done by calling the 2D code.


##--------------------------------------------------------##


bellman1d <- function(F, b_x, u.values, cost,
                      x0, x1, nx,
                      nsteps, dt,
                      skip=0,
                      rescale=TRUE,
                      cutoff=0,
                      infinite.domain=TRUE,
                      debug=FALSE,
                      params=c(),
                      param.weights=c()) {

  bellman2d(function(u, x, y, ...) { F(u,x,...) },
            function(u, x, y, ...) { -y },    ## stub
            b_x,
            1.0,                              ## stub
            u.values,

            function(u, x, y, ...) { cost(u, x, ...) },
            
            x0, x1, nx,
            -1, 1, 1,                         ## stub
            
            nsteps, dt,
            skip=skip,
            rescale=rescale,
            cutoff=0,
            infinite.domain=infinite.domain,
            discretize=discretize2d.cheap,
            debug=debug,
            params=params,
            param.weights=param.weights)
}

bellman1d.control.value <- function(control, t, x) {
  bellman2d.control.value(control, t, x, 0)
}


################################
## Utilities

image.1d.cntl <- function(bc, step, ...) {
  plot(sapply(bc$mc$grid$cells, function(c) { c$xm }),
       bc$cntl[[step]],
       ...)
}

image.1d.cost <- function(bc, step, ...) {
  plot(sapply(bc$mc$grid$cells, function(c) { c$xm }),
       -bc$cost[[step]],
       ...)
}

image.1d.fun <- function(bc, f, ...) {
  x <- sapply(bc$mc$grid$cells, function(c) { c$xm })
  plot(x, sapply(x,f), ...)
}
