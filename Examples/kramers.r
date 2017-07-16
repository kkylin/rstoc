## kramers.r

## This file is part of RSTOC, a system for stochastic
## optimal control and experimental design written in R.

## Copyright (C) 2013,2017 by the authors:

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

## double-well model with adjustable barrier

## default parameters
default.kramers.params <- data.frame(A = 3.84, w = 0.3)

## drift field
kramers.drift <- function(u, x,
                          A   = default.kramers.params$A,
                          w   = default.kramers.params$w) {
  
  4*x - 4*x^3 + A*(x/w)*exp(-0.5*(x/w)^2) + u
}


################################
## Fisher informations

## Given initial condition (X Y), compute the Fisher
## information for time DT, using the value U as control
## over the interval [0,DT].

## Note minus sign: bellman() minimizes cost by default,
## whereas we want to maximize Fisher information.

kramers.fisher <-function(u, x,
                          A   = default.kramers.params$A,
                          w   = default.kramers.params$w) {
  
  -((x/w)*exp(-0.5*(x/w)^2))^2
}
