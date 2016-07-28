## ml.r

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

## ml.r: Morris-Lecar neuron model in the Hopf regime.  See
## [Ermentrout-Terman] for more details.


## default parameters
default.ml.params <- data.frame(I0=0.0,
                                Cm=20.0,
                                g_ca=4.4, g_k=8.0, g_leak=2.0,
                                v_k=-84.0, v_leak=-60.0, v_ca=120.0,
                                phi=0.04,
                                v_1=-1.2, v_2=18.0, v_3=2.0, v_4=30.0,

                                ## Note that by construction, the amplitude is set to 1.
                                ## User scripts should scale this appropriately.
                                channel.noise.amp=1)

## ML vector field

m_inf <-
  with(default.ml.params,
       function(v) {
         0.5*(1.0+(tanh ((v - v_1) / v_2)))
       })
       
w_inf <-
  with(default.ml.params,
       function(v) {
         0.5*(1.0+(tanh ((v - v_3) / v_4)))
       })

tau_w <-
  with(default.ml.params,
       function(v) {
         1.0/(cosh((v - v_3) / 2.0 / v_4))
       })

## Derivatives of the above -- not needed here.

## dm_inf <- function(v) {
##   0.5*((dtanh ((v - v_1) / v_2))/v_2)
## }

## dw_inf <- function(v) {
##   0.5*((dtanh ((v - v_3) / v_4))/v_4)
## }

## dtau_w <- function(v) {
##   (-1.0/(square (cosh ((v - v_3) / 2.0 / v_4))))*
##     ((dcosh ((v - v_3) / 2.0 / v_4))/2.0/v_4)
## }

ml.drift.v <-
  function(u,v,w,
           g_ca=default.ml.params$g_ca,
           v_ca=default.ml.params$v_ca,
           g_k=default.ml.params$g_k,
           v_k=default.ml.params$v_k,
           g_leak=default.ml.params$g_leak,
           v_leak=default.ml.params$v_leak,
           I0=default.ml.params$I0,
           phi=default.ml.params$phi,
           Cm=default.ml.params$Cm) {
  (((-g_ca * (m_inf (v)) * (v - v_ca)) +
    (-g_k * w * (v - v_k))
    + (-g_leak * (v - v_leak)) + I0) / Cm + u) }
  
ml.drift.w <-
  function(u,v,w,
           g_ca=default.ml.params$g_ca,
           v_ca=default.ml.params$v_ca,
           g_k=default.ml.params$g_k,
           v_k=default.ml.params$v_k,
           g_leak=default.ml.params$g_leak,
           v_leak=default.ml.params$v_leak,
           I0=default.ml.params$I0,
           phi=default.ml.params$phi,
           Cm=default.ml.params$Cm) {
    (phi * (((w_inf (v)) - w) / (tau_w (v)))) }
  
## Diffusion coeff for channel noise
channel.noise.coeff <-
  function(v,w,
           g_ca=default.ml.params$g_ca,
           v_ca=default.ml.params$v_ca,
           g_k=default.ml.params$g_k,
           v_k=default.ml.params$v_k,
           g_leak=default.ml.params$g_leak,
           v_leak=default.ml.params$v_leak,
           I0=default.ml.params$I0,
           phi=default.ml.params$phi,
           Cm=default.ml.params$Cm,
           channel.noise.amp=default.ml.params$channel.noise.amp) {

    channel.noise.amp*sqrt(max(0.0 ,
                               phi / tau_w (v) * (w_inf (v) *
                                                  (1.0 - 2.0 * w) + w)))
  }


################################
## Drift

## x and y are assumed to be row vectors

ml.drift.vector <- function(u, state, ...) {
  n <- length(state)
  v <- state[seq(1,n,2)]
  w <- state[seq(2,n,2)]
  vdot <- ml.drift.v(u, v, w, ...)
  wdot <- ml.drift.w(u, v, w, ...)
  if ( is.matrix(state)) {
    matrix(c(vdot,wdot), 2, ncol(state), byrow=T)
  }
  else {
    c(vdot,wdot)
  }
}



################################
## Fisher informations

## Given initial condition (X Y), compute the Fisher
## information for time DT, using the value U as control
## over the interval [0,DT].

## Note minus sign: bellman() minimizes cost by default,
## whereas we want to maximize Fisher information.


## FI for membrane capacitance
ml.fisher.Cm <- function(u,v,w,
                         Cm=default.ml.params$Cm,
                         g_ca=default.ml.params$g_ca,
                         v_ca=default.ml.params$v_ca,
                         g_k=default.ml.params$g_k,
                         v_k=default.ml.params$v_k,
                         g_leak=default.ml.params$g_leak,
                         v_leak=default.ml.params$v_leak,
                         I0=default.ml.params$I0,
                         phi=default.ml.params$phi)
{
  -ml.drift.v(0,v,w)^2
}

## FI for phi
ml.fisher.phi <- function(u,v,w,
                          phi=default.ml.params$phi,
                          g_ca=default.ml.params$g_ca,
                          v_ca=default.ml.params$v_ca,
                          g_k=default.ml.params$g_k,
                          v_k=default.ml.params$v_k,
                          g_leak=default.ml.params$g_leak,
                          v_leak=default.ml.params$v_leak,
                          I0=default.ml.params$I0,
                          Cm=default.ml.params$Cm)
{
  -ml.drift.w(0,v,w)^2
}

## FI for Ca conductance constant
ml.fisher.gca <- function(u,v,w,
                          g_ca=default.ml.params$g_ca,
                          v_ca=default.ml.params$v_ca,
                          g_k=default.ml.params$g_k,
                          v_k=default.ml.params$v_k,
                          g_leak=default.ml.params$g_leak,
                          v_leak=default.ml.params$v_leak,
                          I0=default.ml.params$I0,
                          phi=default.ml.params$phi,
                          Cm=default.ml.params$Cm)
{
  -(m_inf(v) * (v - v_ca))^2
}
