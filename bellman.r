## bellman.r

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

## Given a controlled finite-state Markov chain, this code
## computes the sequence of controls that optimize a given
## cost function.  The controls are assumed to come from a
## finite set, which we think of as a set of integers
## U={1,2,...,U.max}.  The optimal control sequence is
## computed using Bellman's dynamic programming algorithm.

## Currently, there is no end-point cost.  This code is
## specific to Markov chains and has no knowledge of SDEs
## etc.


bellman.mc <- function(trans.matrices, costs, nsteps,
                       average=FALSE ) {

  ## Inputs:
  
  ## TRANS.MATRICES - A finite list of transition
  ##   probability matrices, one for each control value.
  ##   Transition matrices should be represented as sparse
  ##   matrices; see `sparse.r'.  The length of this list
  ##   should be |U|, and each transition matrix should have
  ##   dimensions |S|x|S| (where S = state space).

  ## COSTS - Cost vectors.  This can depend on the control
  ## value or not.  In the first case, a single cost vector
  ## is expected; its length must be |S|.  In the second
  ## case, there should be a cost vector for each control
  ## value, each of which should have length |S|.

  ## NSTEPS - The number of time steps to take.

  ## AVERAGE - Assume the transition matrices depend on a
  ##  parameter, and compute the expected cost by averaging
  ##  over this parameter.  This option, when TRUE, modifies
  ##  the meaning of TRANS.MATRICES and COSTS as follows:
  ##  TRANS.MATRICES is a 2D array of sparse matrices,
  ##  indexed so that TRANS.MATRICES[[k]][[u]] is the
  ##  transition matrix for the kth parameter and the uth
  ##  control value.  Similarly, COSTS[[k]][[s]] is the cost
  ##  associated with the kth parameter and the sth state.

  ##--------------------------------
  ## Parse input flags

  ## To average over parameters or not:
  if ( average ) {
    nparams <- length(trans.matrices)
    ncontrols <- length(trans.matrices[[1]])
    nstates <- trans.matrices[[1]][[1]]$m
  }
  else {
    nparams <- 1
    ncontrols <- length(trans.matrices)
    nstates <- trans.matrices[[1]]$m
  }

  print(paste('nparams=',nparams,'ncontrols=',ncontrols,'nstates=',nstates))

  ## Replicate cost vector if needed, one for each control
  ## value

  if ( is.numeric(costs)) {
    print('bellman(): copying cost vectors')
    costs <- lapply(trans.matrices, function(mat) { costs })
  }
  

  ##--------------------------------
  ## Main loop

  ## Step back in time, building the control sequence one
  ## element at a time.

  ## J = optimal cost
  J.list <- vector(mode="list", nsteps+1)
  J.list[[nsteps+1]] <- vector(mode='numeric', nstates)

  ## F = control
  F.list <- vector(mode="list", nsteps+1)
  F.list[[nsteps+1]] <- vector(mode='numeric', nstates)+1

  for ( step in nsteps:1 ) {

    print(paste('; bellman(): optimizing time step',step-1, 'to', step))
    
    ## grab previous-computed optimal cost
    J <- J.list[[step+1]]

    ## allocate matrices for new time step
    new.J <- vector(mode='numeric', nstates)
    new.F <- vector(mode='numeric', nstates)
    
    ## fill every table entry
    for ( s in 1:nstates ) {

      ## Find best control value, with S as starting state
      ## at timestep STEP.

      ## Implementation Note: If, rather than a finite set U
      ## of control values (i.e., bang-bang control), we
      ## have a continuous interval for U and a nice (e.g.,
      ## quadratic) penalization function, one can skip this
      ## optimization loop and solve for the minimizer
      ## explicitly or numerically.

      min.cost <- Inf
      min.u <- FALSE
      
      for ( u in 1:ncontrols ) {

        new.cost <- 0

        ## Compute cost for this choice of u using
        ## pre-computed transition probabilities and cost
        ## values.

        ## (This case analysis could have been avoided if we
        ## used the format
        ## TRANS.MATRICES[[control-value]][[[param-index]]
        ## instead of the other way around, but the speed
        ## difference is <= 1%.

        if ( average ) {
          for ( k in 1:nparams ) {
            P <- sparse.matrix.row(trans.matrices[[k]][[u]], s)
            new.cost <- new.cost + (sparse.vector.dot(P,J) + costs[[u]][[k]][[s]]) / nparams
          }
        }
        else {
          P <- sparse.matrix.row(trans.matrices[[u]], s)
          new.cost <- sparse.vector.dot(P,J) + costs[[u]][[s]]
        }

        if ( new.cost  < min.cost ) {
          min.cost <- new.cost
          min.u <- u
        }
      }

      ## Save computed best cost & control.
      new.J[[s]] <- min.cost
      new.F[[s]] <- min.u
    }

    ## Save the results for the current time step.
    J.list[[step]] <- new.J
    F.list[[step]] <- new.F
  }

  ## Output everything in a big data structure.
  list(nsteps=nsteps,
       cost=J.list,
       cntl=F.list)
}
