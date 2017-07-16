## bellman2d.r

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

##   dX_t = F(u_t,X_t,Y_t) dt + b_x dW^1_t
##   dY_t = G(u_t,X_t,Y_t) dt + b_y dW^2_t

## and a cost functional of the form

##   COST[u] := \int_0^T C(u_t, X_t, Y_t) dt

## we compute the control that minimizes COST[u].  This is
## done in two steps:

## 1) A finite-state Markov chain approximation is computed
##    using a simple discretization method (see below).

## 2) The standard dynamic programming algorithm is applied
##    to the finite-state Markov chain.

## SPACE DISCRETIZATION

## This code works with rectangular domains R in the plane.
## The rectangle R is subdivided into a regular grid of
## cells of area dx-by-dy.  By default, the code can include
## semi-infinite cells, so that the entire plane (not just
## the chosen rectangle R) is included in the
## discretization.  This can be turned off by the user, in
## which case the approximating Markov chain is computed by
## conditioning on the event that the continuous-time NEVER
## exits R.


##--------------------------------------------------------##





###################################################################
## This is the main function in this file.

bellman2d <- function(F, G, b_x, b_y, u.values, cost,
                      x0, x1, nx,
                      y0, y1, ny,
                      nsteps, dt,
                      skip=0,
                      rescale=TRUE,
                      cutoff=0,
                      infinite.domain=TRUE,
                      discretize=discretize2d.cheap,
                      debug=FALSE,
                      params=c(),
                      param.weights=c()) {
  
  ##--------------------------------
  ## Inputs:

  ## F, G, b_x, b_y - As above.

  ## U.VALUES - List of control values.

  ## COST - A function of u, x, and y that returns the cost
  ##   associated with the point (x,y) and control u.

  ## X0, X1, NX, Y0, Y1, NY - The rectangular domain is
  ##   [X0,X1]X[Y0,Y1].  NX and NY are the number of cells
  ##   in the X and Y directions, respectively.

  ## NSTEPS - Number of time steps.

  ## DT - Size of each time step.

  ## SKIP - skip every SKIP timesteps.  Default=0, i.e., no
  ## skipping.  If SKIP>0, then the actual number of
  ## timesteps taken is CEILING(NSTEPS/(SKIP+1)), with
  ## timestep (SKIP+1)*DT.

  ## RESCALE - Whether to rescale the transition matrix so
  ##   each row sums to 1.  This is provided for debugging
  ##   purposes.  Defaults to TRUE.

  ## CUTOFF - The cutoff probability value for tails of
  ##   transition matrices.  This is a cheap way to sparsify
  ##   transition matrices; as long as the cutoff is
  ##   O(DT^2), the scheme should still converge.  Defaults
  ##   to 0 (no cutoff).

  ## INFINITE.DOMAIN - Whether to discretize all of R^2 or
  ##   just the rectangle R.  Defaults to TRUE (all of R^2).

  ## PARAMS - A list of parameters to pass to F and G.  By
  ## default, this is a null list, which means we assume
  ## F=F(u,x,y) and G=G(u,x,y) ((x,y) = system state, u =
  ## control value).  If a list (THETA0,THETA1,...) is
  ## given, then we assume F=F(u,x,y,theta) and
  ## G=G(u,x,y,theta).


  ##--------------------------------
  ## Output:

  ## A lot of garbge.  Main thing is COST, a list of
  ## matrices containing the expected optimized cost for
  ## each time step and each state, and CNTL, a list of
  ## matrices containing the optimal control value.

  nparams <- length(params)
  ncontrols <- length(u.values)
  nsteps.actual <- ceiling(nsteps/(skip+1))
  dt.actual <- (skip+1)*dt

  if ( nparams > 0 ) {

    onestep.matrices <- vector(mode='list',length(params))
    trans.matrices <- vector(mode='list',length(params))

    if ( length(param.weights) == 0 ) {
      ## no weights supplied: assume uniform distribution
      param.weights <- sapply(params, function(z) {1})
    }
    else if ( length(param.weights) == nparams ) {
      
      ## Use given weights

      ## This scaling will trick bellman.mc() into computing
      ## the weighted average instead of the unweighted
      ## average.
      total.weight <- sum(param.weights)
      param.weights <- param.weights / total.weight * nparams
    }
    else {
      ## programming error?
      stop('number of weights != number of params')
    }

    ## Generate transition matrices
    for ( k in 1:nparams ) {

      theta <- params[[k]]
      print(paste('# bellman2d(): discretizing theta=',theta))

      print(system.time( mc <- discretize(function(u,x,y) {F(u,x,y,theta)},
                                          function(u,x,y) {G(u,x,y,theta)},
                                          b_x, b_y, u.values,
                                          x0, x1, nx, y0, y1, ny,
                                          dt,
                                          rescale=rescale,
                                          cutoff=cutoff,
                                          ## don't do this; see below
                                          ## skip=skip,
                                          infinite.domain=infinite.domain)))

      ## these are needed to compute cost vectors correctly
      ## when skip>0
      onestep.matrices[[k]] <- lapply(mc$mat,
                                      function(mat) {
                                        sparse.matrix.scale(param.weights[[k]], mat)
                                      })

      trans.matrices[[k]] <- lapply(mc$mat,
                                    function(mat) {
                                      sparse.matrix.scale(param.weights[[k]],
                                                          sparse.matrix.pow(mat, skip+1))
                                    })
      
      ## The grid should be the same every time; doesn't
      ## hurt to save it here.
      g <- mc$grid
      dx <- g$dx
      dy <- g$dy
    }
  }
  else {
    ## old-style drift fields: no need to pass parameters
    print(system.time( mc <- discretize(F, G, b_x, b_y, u.values,
                                        x0, x1, nx, y0, y1, ny,
                                        dt,
                                        rescale=rescale,
                                        cutoff=cutoff,
                                        ## don't do this; see below
                                        ## skip=skip,
                                        infinite.domain=infinite.domain)))
    onestep.matrices <- mc$mat
    
    trans.matrices <- lapply(mc$mat,
                             function(mat) {
                               sparse.matrix.pow(mat, skip+1)
                             })
    
    g <- mc$grid
    dx <- g$dx
    dy <- g$dy
  }

  ## Actual grid may be larger, if we asked for infinite domain
  nx <- g$nx
  ny <- g$ny

  ## compute costs, making sure ordering is same as transition matrices
  if ( nparams > 0 ) {

    costs <- vector(mode='list', ncontrols)

    for ( u in 1:ncontrols ) {
      for ( k in 1:nparams ) {

        print(paste('computing cost vectors for control', u, 'param', k))

        theta <- params[[k]]
        costs[[u]][[k]] <- sapply(g$cells,
                                  function(cell) {
                                    cost(u.values[[u]], cell$xm, cell$ym, theta)
                                  })
        
        costs[[u]][[k]] <- geometric.sum(onestep.matrices[[k]][[u]], costs[[u]][[k]], skip)
      }
    }

    ## optimize
    print(system.time(
                      opt.policy <-
                      bellman.mc(trans.matrices, costs, nsteps.actual,
                                 average=TRUE)))
  }
  else {

    costs <- vector(mode='list', ncontrols)

    for ( u in 1:ncontrols ) {

      print(paste('computing cost vectors for control', u))

      costs[[u]] <- sapply(g$cells, function(cell) { cost(u.values[[u]], cell$xm, cell$ym) })
      costs[[u]] <- geometric.sum(onestep.matrices[[u]], costs[[u]], skip)
    }
      
    ## optimize
    print(system.time(
                      opt.policy <-
                      bellman.mc(trans.matrices, costs, nsteps.actual,
                                 average=FALSE)))
  }
  
  print(paste('# bellman2d(): repackaging control policy'))

  cntl <- vector(mode='list',nsteps.actual+1)
  cost <- vector(mode='list',nsteps.actual+1)

  for ( k in 1:(nsteps.actual+1) ) {

    ## print(paste('# bellman2d(): repackaging step', k))

    opt.cntl <- opt.policy$cntl[[k]]
    opt.cost <- opt.policy$cost[[k]]
    ## a <- matrix(0,nx,ny)
    ## b <- matrix(0,nx,ny)
    
    ## for ( ind in 1:length(g$cells) ) {

    ##   cell <- g$cells[[ind]]
    ##   i <- cell$i
    ##   j <- cell$j

    ##   a[i,j] <- u.values[[opt.cntl[[ind]]]]
    ##   b[i,j] <- opt.cost[[ind]]
    ## }

    cntl[[k]] <- matrix(u.values[opt.cntl], nx, ny)
    cost[[k]] <- matrix(opt.cost, nx, ny)
  }


  ## Output computed policy and cost, plus discretization
  ## grid etc.  In debug mode, also output transition
  ## matrices.

  list(nsteps=nsteps,
       dt=dt,
       nsteps.actual=nsteps.actual,
       dt.actual=dt.actual,
       mc={ if ( debug ) {
              mc
            }
            else {
              list(mat=FALSE,grid=mc$grid)
            }},
       x0=x0, x1=x1, nx=nx, dx=g$dx,
       y0=y0, y1=y1, ny=ny, dy=g$dy,
       cost=cost,
       cntl=cntl,
       policy={ if ( debug ) {
                  opt.policy
                }
                else {
                  FALSE
                }})
}


## User-friendly way to use the output of bellman2d().  This
## version is now deprecated in favor of the one below.
bellman2d.control <- function(control, k, x, y) {
  g <- control$mc$grid
  i <- max(1,min(g$nx,length(which(g$x <= x))))
  j <- max(1,min(g$ny,length(which(g$y <= y))))
  control$cntl[[k]][i,j]
}


## A better way to compute control values.  Note that T
## should be a scalar, but X and Y can be vectors.  This
## implementation is a bit more general than it needs to be:
## it leaves room for unequal grids.


## This version is more generally applicable.

bellman2d.control.value.general <- function(control, t, x, y) {
  g <- control$mc$grid
  i <- pmax(1,pmin(g$nx,colSums(outer(g$x, x, '<='))))
  j <- pmax(1,pmin(g$ny,colSums(outer(g$y, y, '<='))))
  k <- 1 + min(control$nsteps.actual,max(0,round(t/control$dt.actual)))
  control$cntl[[k]][(j-1)*g$nx+i]
}

## This version is a lot faster.

bellman2d.control.value <- function(control, t, x, y) {
  g <- control$mc$grid
  i <- pmax(1,pmin(g$nx,1+floor((x-g$x0)/g$dx) + ( g$x0 != g$x[[1]] )))
  j <- pmax(1,pmin(g$ny,1+floor((y-g$y0)/g$dy) + ( g$y0 != g$y[[1]] )))
  k <- 1 + min(control$nsteps.actual,max(0,round(t/control$dt.actual)))
  control$cntl[[k]][(j-1)*g$nx+i]
}



## Old version: using this usually requires wrapping an
## mapply() around it; it's about 1/3 the speed.

## bellman2d.control.value <- function(control, t, x, y) {
##   g <- control$mc$grid
##   i <- max(1,min(g$nx,sum(g$x <= x)))
##   j <- max(1,min(g$ny,sum(g$y <= y)))
##   k <- 1 + min(control$nsteps.actual,max(0,floor(t/control$dt.actual)))
##   control$cntl[[k]][i,j]
## }




## Used in skipping steps
geometric.sum <- function(A, v0, n) {

  v <- v0

  while ( n > 0 ) {
    v <- sparse.matrix.vector.mult(A,v) + v0
    n <- n-1
  }

  v
}


############################################################################

## Here we discretize the SDE, i.e., we compute time-t
## transition matrices for a Markov chain approximation to
## the SDE

## dX = F(u,X,Y) + b_x dW_x
## dY = G(u,X,Y) + b_y dW_y

## using an uniform MxN grid over the rectangle
## [x0,x1]x[y0,y1].  This is done for all control values in
## U.VALUES.

## A few methods were tried:

## 1) A transition-density-approximation method

## 2) A split-operator finite-difference scheme

## 3) A hybrid scheme, where advection is computed by finite
## differences and diffusion by density-approximation

## Advantage of Method 1 is that it is reasonably fast
## because of using R's optimized internal routines, and it
## tolerates (in fact, requires) larger timesteps for
## convergence, i.e., dt >> dx^2.  The finite difference
## scheme is (I think) convergent according to standard
## theory, but needs dt = O(dx^2).

## In the end, a version of a finite difference scheme seems
## to strike the best balance betwen speed, memory usage,
## and convergence behavior: the scheme discretizes the
## diffusion operator by adaptively scaling up the spatial
## grid size, so there are no essential time step
## constraints (but one needs to maintain dt = O(dx^2) to
## get convergence).  So it is now the default; see
## discretize2d.cheap() below.  The other discretization
## routines can be `plugged in' easily; they are defined in
## `bellman2d-extra.r'.

## Usage: inputs are basically the same as bellman2d().  In
## all cases, the function returns a list of transition
## matrices and a grid structure (see make.grid2d()).  The
## grid points correspond to the centers of cells.



## Method 2A: Discretize the Fokker-Planck equation by a
## split-operator finite difference scheme, but with a cheat
## to avoid iterating.  Less accurate, but cheaper.

discretize2d.cheap <- function(F, G, bx_val_or_fun, by_val_or_fun, u.values,
                               x0, x1, nx,
                               y0, y1, ny,
                               dt,
                               rescale=TRUE,
                               cutoff=0,
                               skip=0,
                               infinite.domain=FALSE) {

  ## Set grid size and timestep.
  g <- make.grid2d(x0, x1, nx, y0, y1, ny, infinite.domain)
  dx <- (x1-x0)/nx
  x.ratio <- 1.5/dx
  dy <- (y1-y0)/ny
  y.ratio <- 1.5/dy
  nx <- g$nx   ## actual domain may be bigger if infinite.domain==TRUE
  ny <- g$ny

  ## Support variable diffusion coefficients
  sqrt_dt <- sqrt(dt)

  b_x <- if ( is.function(bx_val_or_fun)) {
    function(x,y) {
      bx_val_or_fun(x,y) * sqrt_dt
    }
  }
  else {
    function(x,y) {
      bx_val_or_fun * sqrt_dt
    }
  }
    
  b_y <- if ( is.function(by_val_or_fun)) {
    function(x,y) {
      by_val_or_fun(x,y) * sqrt_dt
    }
  }
  else {
    function(x,y) {
      by_val_or_fun * sqrt_dt
    }
  }

  ## R matrices use column-major format
  cell.index <- function(i,j) {
    i+(j-1)*nx
  }

  ## Map grid to integers
  safe.set <- function(A, src.ind, i, j, val) {
    if ( val < 0 || val > 1 ) {
      print(paste('safe.set:', i, j, src.ind, val))
      stop('')
    }

    sparse.matrix.inc(A,
                      src.ind,
                      cell.index(max(1,min(i,nx)), max(1,min(j,ny))),
                      val)
  }

  ## Compose advection and diffusion matrices on the fly
  set.matrix.elt <- function(A, src.ind, i, j, xmult, ymult, p, q.horiz, q.vert) {
    
    safe.set(A, src.ind, i-xmult, j,       p * q.horiz)
    safe.set(A, src.ind, i+xmult, j,       p * q.horiz)
    safe.set(A, src.ind, i,       j-ymult, p * q.vert)
    safe.set(A, src.ind, i,       j+ymult, p * q.vert)
    safe.set(A, src.ind, i,       j,       p * ( 1-2*( q.horiz + q.vert )))
  }

  
  if ( infinite.domain ) {
    index.shift <- 2
  }
  else {
    index.shift <- 1
  }

  ## Track how many cells are shifted
  cell.shift.count <- 0
  
  ## Generate separate transition matrices for each control
  ## value.
  generate.trans.matrix <- function(u) {

    print(paste('# u=',u))
    new.mat <- sparse.matrix(nx*ny,nx*ny)

    ## Compute transition operator
    print('## computing matrix elements')

    xmult.min <- ymult.min <- Inf
    xmult.max <- ymult.max <- 0

    for ( cell in g$cells ) {

      ## index of source grid cell
      i <- cell$i
      j <- cell$j
      src.ind <- cell.index(i,j)
      x <- cell$xm
      y <- cell$ym

      ## Advection matrix elements: go forward in time along
      ## characteristic, then interpolate.
      z <- euler_step_2d(F,G,u,x,y,dt)
      zx <- z[[1]]
      zy <- z[[2]]

      zi <- floor((zx-x0)/dx-0.5) + index.shift
      zj <- floor((zy-y0)/dy-0.5) + index.shift

      if ( cell.index(zi,zj) != src.ind ) {
        cell.shift.count <- cell.shift.count+1
      }

      px <- (zx - x0)/dx-0.5 - (zi - index.shift)
      py <- (zy - y0)/dy-0.5 - (zj - index.shift)

      ## Diffusion matrix elements
      bx.dw <- b_x(x,y)
      by.dw <- b_y(x,y)
      xk <- ceiling( bx.dw * x.ratio )
      yk <- ceiling( by.dw * y.ratio )

      xmult.min <- min(xmult.min, xk)
      xmult.max <- max(xmult.max, xk)
      ymult.min <- min(ymult.min, yk)
      ymult.max <- max(ymult.max, yk)

      q.horiz <- 0.5*(bx.dw/xk/dx)^2
      q.vert  <- 0.5*(by.dw/yk/dy)^2

      ## Set matrix elements
      set.matrix.elt(new.mat, src.ind, zi,   zj,   xk, yk, (1-px)*(1-py), q.horiz, q.vert)
      set.matrix.elt(new.mat, src.ind, zi+1, zj,   xk, yk, px*(1-py),     q.horiz, q.vert)
      set.matrix.elt(new.mat, src.ind, zi+1, zj+1, xk, yk, px*py,         q.horiz, q.vert)
      set.matrix.elt(new.mat, src.ind, zi,   zj+1, xk, yk, (1-px)*py,     q.horiz, q.vert)
    }

    print(paste('## xmult=', xmult.min, xmult.max, 'ymult=', ymult.min, ymult.max))

    ## Spit out some interesting stats.
    print(paste('## sparsity:', sparsity(new.mat)))
    print(paste('## trace-frac:', sparse.matrix.trace(new.mat)/nx/ny))
    print(paste('## diag-max:', sparse.matrix.diagmax(new.mat)))

    if ( skip > 0 ) {
      print(paste('## skipping', skip, 'steps between observations'))
      new.mat <- sparse.matrix.pow(new.mat,skip+1)
      print(paste('## sparsity:', sparsity(new.mat)))
      print(paste('## trace-frac:', sparse.matrix.trace(new.mat)/nx/ny))
      print(paste('## diag-max:', sparse.matrix.diagmax(new.mat)))
    }
    
    new.mat
  }

  print(paste('## cell.shift.count:', cell.shift.count, '/', length(g$cells)))

  P.list <- lapply(u.values, generate.trans.matrix)
  list(mat=P.list, grid=g)
}

rk2_step_2d <- function(F, G, u, x, y, dt) {
  kx <- x + 0.5*dt*F(u,x,y)
  ky <- y + 0.5*dt*G(u,x,y)
  c(x+dt*F(u,kx,ky), y+dt*G(u,kx,ky))
}

euler_step_2d <- function(F, G, u, x, y, dt) {
  c(x+dt*F(u,x,y), y+dt*G(u,x,y))
}




########################################################
## Generate a 2D grid

## This function constructs a 2D grid.  Each cell is
## explicitly represented.  A `safe midpoint' is also
## computed.  This is just the usual midpoint if the cell
## has finite area, but a point on a finite edge or corner
## if the cell is semi-infinite.  The `safe midpoint' is a
## convenient place to evaluate quantities like transition
## probabilities and costs.

make.grid2d <- function(x0, x1, nx, y0, y1, ny, infinite.domain=FALSE)
{
  dx <- (x1-x0) / nx
  dy <- (y1-y0) / ny
  X <- x0 + dx*(0:nx)
  Y <- y0 + dy*(0:ny)

  if ( infinite.domain ) {
    X <- c(-Inf,X,Inf)
    Y <- c(-Inf,Y,Inf)
    nx <- nx+2
    ny <- ny+2
  }

  cells <- outer(1:(length(X)-1), 1:(length(Y)-1),
                 function(I,J) {
                   mapply(
                          function(i,j) {

                            x0 <- X[i]
                            x1 <- X[i+1]
                            y0 <- Y[j]
                            y1 <- Y[j+1]

                            list(i=i, j=j,
                                 xm=safe.midpoint(x0,x1),
                                 ym=safe.midpoint(y0,y1),
                                 x0=X[i], x1=X[i+1],
                                 y0=Y[j], y1=Y[j+1])
                          },
                          I, J,
                          SIMPLIFY=FALSE)
                 })
  
  list(x0=x0, x1=x1, nx=nx, dx=dx, x=X,
       y0=y0, y1=y1, ny=ny, dy=dy, y=Y,
       cells=cells)
}

safe.midpoint <- function(x0,x1) {
  if ( abs(x0) == Inf ) {
    x1
  }
  else if ( abs(x1) == Inf ) {
    x0
  }
  else {
    (x1+x0)/2
  }
}


################################
## Utilities

image.2d.cntl <- function(bc, step, ...) {
  image.grid.fun(bc$mc$grid, bc$cntl[[step]], ...)
}

image.2d.cost <- function(bc, step, ...) {
  image.grid.fun(bc$mc$grid, -bc$cost[[step]], ...)
}

image.2d.fun <- function(bc, f, ...) {
  image.grid.fun(bc$mc$grid,
                 sapply(bc$mc$grid$cells,
                        function(c) {
                          f(c$xm, c$ym)
                        }),
                 ...)
}
