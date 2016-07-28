## util.r

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

## Utility functions.

## run.vector.program(prog, input.vector): Run the program
## PROG with the numeric vector INPUT.VECTOR as input, and
## returns a numeric vector as output.  More precisely, the
## entries in INPUT.VECTOR will be fed into PROG, one number
## per line, via stdin.  The results from PROG (from its
## stdout) are parsed and returned as a big vector.  Note
## that the output from PROG must contain ONE NUMBER PER
## LINE, or as.numeric() will fail.

## This was intended to be a easy-to-use interface to code
## written in lower-level languages like C.  Not used.

## run.vector.program <- function(prog, input.vector) {
##   as.vector(sapply(system(prog,
##                           input=sapply(input.vector,toString),
##                           intern=TRUE),
##                    as.numeric))
## }


## pp(X,file): Pretty-print the object X (should be a vector
## or a matrix) as an s-expression, for easier parsing in
## other environments, e.g., Lisp or Matlab.

pp <- function(X,file='',varname='',mode='matlab') {
  if ( mode=='lisp' ) {
    if ( varname=='') {
      pp.lisp(X, '', file=file)
    }
    else {
      cat('(', varname, file=file, append=TRUE)
      pp.lisp(X, '', file=file)
      cat(')', file=file, append=TRUE)
    }
  }
  else {
    if ( varname=='') {
      pp.matlab(X, '', file=file)
      cat('\n',file=file,append=TRUE)
    }
    else {
      cat(varname,'=',file=file,append=TRUE)
      pp.matlab(X, '', file=file)
      cat(';\n',file=file,append=TRUE)
    }
  }
}

pp.lisp <- function(X,prefix,file='') {
  if ( is.list(X) ) {
    cat(prefix,'(\n', file=file, append=TRUE)
    for ( x in X ) {
      pp.lisp(x,paste(prefix,' '), file=file)
    }
    cat(prefix,')\n', file=file, append=TRUE)
  }
  else if ( is.vector(X)) {
    cat(prefix, '(', X, ')\n', file=file, append=TRUE)
  }
  else if ( is.matrix(X)) {
    m <- nrow(X)
    pp.lisp(lapply(1:m, function(i) { X[i,] }), prefix, file=file)
  }
  else {
    cat(prefix, X, file=file, append=TRUE)
  }
}

pp.matlab <- function(X,prefix,file='') {
  if ( is.list(X) ) {
    cat(prefix,'[', file=file, append=TRUE)
    for ( x in X ) {
      pp.matlab(x,paste(prefix,' '),file=file)
      cat('\n', file=file, append=TRUE)
    }
    cat(prefix,']', file=file, append=TRUE)
  }
  else if ( is.vector(X)) {
    cat(prefix, '[', X, ']', file=file, append=TRUE)
  }
  else if ( is.matrix(X)) {
    m <- nrow(X)
    pp.matlab(lapply(1:m, function(i) { X[i,] }), prefix, file=file)
  }
  else {
    cat(prefix, X, file=file, append=TRUE)
  }
}





## Given a grid G and a list A, image it as a matrix.

image.grid.fun <- function(g, A, colorbar=TRUE, ...) {

  nx <- g$nx
  ny <- g$ny

  if ( is.numeric(A)) {
    A <- matrix(A, nx, ny)
  }

  ## for debugging
  ## print(paste(length(g$x), length(g$y), dim(A)))

  xl <- filter.out.inf(g$x,g$dx)
  yl <- filter.out.inf(g$y,g$dy)
  xl <- 0.5 * (xl[1:nx] + xl[2:(nx+1)])
  yl <- 0.5 * (yl[1:ny] + yl[2:(ny+1)])

  if ( colorbar ) {
    filled.contour(xl, yl, A, ...)
  }
  else {
    image(xl, yl, A, col=cm.colors(20), ...)
  }
}

filter.out.inf <- function(X,dx) {
  I <- (abs(X) != Inf)
  lo <- min(X[I])
  hi <- max(X[I])
  X[X==Inf] <- hi+dx
  X[X==-Inf] <- lo-dx
  X
}



## Overlay trajectory on top of grid image
plot.over.image <- function(g, A, traj, i, j) {
  image.grid.fun(g, A, colorbar=FALSE)
  par(new=T)
  ## plot(pathref(traj,i), pathref(traj,j),
  ##      xlim=c(g$x0,g$x1), ylim=c(g$y0,g$y1))
  lines(pathref(traj,i), pathref(traj,j))
}


## Set current device, creating a new figure if necessary.
dev.try.set <- function(i) {
  for ( j in dev.list()) {
    if ( j == i ) {
      dev.set(i)
      return()
    }
  }

  dev.new()
}

## More convenient variant
dev.try.next.pointer <- 2

dev.try.next.reset <- function() {
  dev.try.next.pointer <<- 2
}

dev.try.next <- function() {
  dev.try.set(dev.try.next.pointer)
  dev.try.next.pointer <<- dev.try.next.pointer+1
}

## Close all graphics windows
close.all <- function() {
  for ( d in dev.list()) dev.off(d)
  dev.try.next.reset()
}


## Plot a picture in a new graphics device, and set labels
## etc.  Packaging things up this way makes it easier to use
## the same code for on-screen plotting and for printing.

do.plot <- function(thunk,
                    title='',
                    xlab='Membrane voltage V',
                    ylab='Gating variable W',
                    width=5,
                    height=4,
                    print=FALSE) {
  if ( print ) {
    ## print to file
    postscript(paste(title, '.eps', sep=''),
               width=width, height=height,
               horizontal=FALSE,
               onefile=FALSE,
               paper='special')
    thunk()
    title(main=paste(title, sep=''),
          xlab='Membrane voltage V',
          ylab='Gating variable W')
    dev.off()
  }
  else {
    ## plot to screen
    dev.try.next()
    thunk()
    title(main=paste(title, sep=''),
          xlab=xlab,
          ylab=ylab)
  }
}


## Modify WITH so one can specify defaults.  This does *not*
## work as intended because of R's scoping rules.

## with.defaults <- function(defaults, new.params, ...) {
##   with(defaults, with(new.params, ...))
## }
