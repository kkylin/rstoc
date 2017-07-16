## sparse.r

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

################################
## To do things efficiently (*), we need sparse matrices.  R
## has sparse matrices, for example in the Matrix library.
## But I just couldn't figure out how to make the sparse
## matrix stuff in Matrix work well.  So here are some crude
## routines that are just good enough.

## (*) Primary issue is speed, though using less memory is
##     also nice (and contributes to speed -- less GC).

## TO DO: better memory usage via dataframes (instead of
## lists containing methods)

## Note that it is possible to make this a *lot* faster, by
## using a better representation and better algorithms for
## matrix-multiply.  Idea for new design:

## Rep: Keep a list of nonzero entries, in the form (I,J,A).
## The list should be sorted in lexicographic order, either
## I first or J first.

## Performance: suppose the matrix is NxN and there are, on
## average, k nonzero entries per row.  Then references are
## O(log(k)) (via a bisection search), vs O(k) in the naive
## implementation.  Same for assignment.  Matrix multiply
## can be done in O(k*N^2) instead of the naive O(N^3),
## which -- while not fantastic -- is certainly a
## significant speed-up.

## One problem with this design is that getting a row vector
## is more expensive.

as.sparse.vector <- function(A) {

  ## < 1/3 full: sparse matrix.  Note we use sparse rep for
  ## NULL matrices.
  if ( 3*length(which(A!=0)) <= length(A)) {

    I <- which(A!=0)
    V <- A[I]

    ref <- function(k) {
      if ( any(I==k)) {
        V[which(I==k)]
      }
      else {
        0
      }
    }

    set <- function(k,val) {
      if ( any(I==k)) {
        V[which(I==k)] <<- val
      }
      else {
        n <- length(I)
        I[[n+1]] <<- k
        V[[n+1]] <<- val
      }
    }

    dot <- function(B) {
      if ( is.sparse.vector(B) ) {
        sum(V * sparse.to.dense.vector(B)[I])
      }
      else {
        sum(V * B[I])
      }
    }

    rep <- function() {
      list(indices=I, values=V)
    }

    ## get.indices <- function() {I}
    ## get.values <- function() {V}
    ## list(ref=ref, set=set, dot=dot, indices=get.indices, values=get.values)
    list(ref=ref, set=set, dot=dot, rep=rep)
  }
  ## >= 1/3 full: dense matrix
  else {

    ref <- function(k) { A[[k]] }
    set <- function(k,val) { A[[k]] <<- val }
    dot <- function(B) {
      if ( is.sparse.vector(B) ) {
        sum(A * sparse.to.dense.vector(B))
      }
      else {
        sum(A*B)
      }
    }

    rep <- function() { A }

    ## get.values <- function() {A}
    ## list(ref=ref, set=set, dot=dot, values=get.values)
    list(ref=ref, set=set, dot=dot, rep=rep)
  }
}

sparse.vector.dot <- function(A,B) { A$dot(B) }
sparse.vector.ref <- function(A,k) { A$ref(k) }
sparse.vector.set <- function(A,k,val) { A$set(k,val) }

sparse.vector.nnzero <- function(A) {
  r <- A$rep()
  if ( is.numeric(r)) {
    length(which(r!=0))
  }
  else {
    length(r$indices)
  }
}

sparse.vector.nnzero.indices <- function(A) {
  r <- A$rep()
  if ( is.numeric(r)) {
    which(r!=0)
  }
  else {
    r$indices
  }
}

sparse.vector.nnzero.values <- function(A) {
  r <- A$rep()
  if ( is.numeric(r)) {
    r[which(r!=0)]
  }
  else {
    r$values
  }
}

sparse.vector.sum <- function(A) {
  r <- A$rep()
  if ( is.numeric(r)) {
    sum(r)
  }
  else {
    sum(r$values)
  }
}

sparse.to.dense.vector <- function(A,n=0) {
  r <- A$rep()
  if ( is.numeric(r)) {
    r
  }
  else {

    ## This test is not necessary, but avoids an annoying
    ## error message for empty vectors.
    if ( length(r$indices) > 0 ) {
      n <- max(n,max(r$indices))
    }
    else {
      n <- max(n,0)
    }
    
    v <- vector(mode='numeric', n)

    if ( length(r$indices) > 0 ) {
      for ( k in 1:length(r$indices)) {
        v[[r$indices[[k]]]] <- r$values[[k]]
      }
    }
    
    v
  }
}



## A kludge to tell dense and sparse vectors apart
is.sparse.vector <- function(v) {
  ## length(attributes(v)$names) == 4
  !is.atomic(v)
}
  
  


## Sparse matrices are stored row-wise.

sparse.matrix <- function(m,n) {

  rows <- vector(mode="list", m)

  for ( i in 1:m ) {
    rows[[i]] <- as.sparse.vector(c())
  }

  ref <- function(i,j) {
    if ( 1 <= i && i <= m && 1 <= j && j <= n ) {
      rows[[i]]$ref(j)
    }
    else {
      0
    }
  }
  
  set <- function(i,j,val) {
    if ( 1 <= i && i <= m && 1 <= j && j <= n ) {
      rows[[i]]$set(j,val)
    }
  }
  
  get.row <- function(i) { rows[[i]] }
  set.row <- function(i,new.row) { rows[[i]] <<- new.row }
  rep <- function() { rows }

  list(m=m, n=n, ref=ref, set=set, row=get.row, set.row=set.row, rep=rep)
}

sparse.matrix.dim <- function(A) { c(A$m, A$n) }
sparse.matrix.nrows <- function(A) { A$m }
sparse.matrix.ncols <- function(A) { A$n }
sparse.matrix.ref <- function(A,i,j) { A$ref(i,j) }
sparse.matrix.set <- function(A,i,j,val) { A$set(i,j,val) }
sparse.matrix.inc <- function(A,i,j,val) {
  ## I have no f***ing idea why I have to write it this way
  ## to make it work.  The simpler A$set(i,j,A$ref(i,j)+val)
  ## messes up the data structure!  Are R procedures not
  ## composable?!
  old <- A$ref(i,j)
  A$set(i,j,old+val)
}
sparse.matrix.row <- function(A,i) { A$row(i) }
sparse.matrix.set.row <- function(A,i,new.row) { A$set.row(i,new.row) }


## convert dense matrix to sparse matrix

as.sparse.matrix <- function(Am) {
  m <- nrow(Am)
  n <- ncol(Am)
  A <- sparse.matrix(m,n)

  for ( i in 1:m ) {
    sparse.matrix.set.row(A,i,as.sparse.vector(Am[i,1:n]))
  }

  A
}


## Functions useful mainly for debugging and testing

sparsity <- function(A) {
  sum(sapply(A$rep(), sparse.vector.nnzero)) / A$m / A$n
}

sparse.matrix.row.sums <- function(A) {
  sapply(A$rep(), sparse.vector.sum)
}

sparse.matrix.trace <- function(A) {
  n <- min(sparse.matrix.dim(A))
  sum <- 0

  for ( i in 1:n ) {
    sum <- sum + sparse.matrix.ref(A,i,i)
  }

  sum
}

sparse.matrix.diagmax <- function(A) {
  n <- min(sparse.matrix.dim(A))
  m <- 0

  for ( i in 1:n ) {
    m <- max(m,sparse.matrix.ref(A,i,i))
  }

  m
}

sparse.matrix.range <- function(A) {
  range(sapply(A$rep(),
               function(row) {
                 r <- row$rep()
                 if ( is.numeric(r)) {
                   range(r)
                 }
                 else {
                   range(r$values)
                 }
               }
               ))
}

sparse.matrix.transpose <- function(A) {
  as.sparse.matrix(t(sparse.to.dense.matrix(A)))
}


## sparse.matrix.transpose <- function(A) {
##   d <- sparse.matrix.dim(A)
##   m <- d[[1]]
##   n <- d[[2]]
##   B <- sparse.matrix(n,m)
  
##   for ( i in 1:m ) {
##     for ( j in 1:n ) {
##       a <- sparse.matrix.ref(A,i,j)
##       if ( a != 0 ) {
##         sparse.matrix.set(B,j,i,a)
##       }
##     }
##   }

##   B
## }

sparse.to.dense.matrix <- function(A) {
  m <- sparse.matrix.dim(A)[[1]]
  n <- sparse.matrix.dim(A)[[2]]
  C <- matrix(0,m,n)

  for ( i in 1:m ) {

    v <- sparse.matrix.row(A,i)$rep()
    
    if ( is.sparse.vector(v) ) {

      ## if ( is.atomic(v)) {
      ##   print(attributes(v))
      ## }

      ind <- v$indices
      val <- v$values

      C[i,ind] <- val
    }
    else {
      C[i,1:n] <- v
    }
  }

  C
}


## Scale a matrix

sparse.matrix.scale <- function(c,A) {
  as.sparse.matrix(c * sparse.to.dense.matrix(A))
}


## This version turns out to be the fastest.  Alas, uses quite a bit of memory...

sparse.matrix.mult <- function(A,B) {

  m <- sparse.matrix.nrows(A)
  n <- sparse.matrix.ncols(B)
  Bm <- sparse.to.dense.matrix(B)
  C <- sparse.matrix(m, n)

  for ( i in 1:m ) {

    row <- sparse.matrix.row(A,i)
    I <- sparse.vector.nnzero.indices(row)
    V <- sparse.vector.nnzero.values(row)

    if ( length(I) > 1 ) {
      sparse.matrix.set.row(C, i, as.sparse.vector(V %*% Bm[I,]))
    }
    else if ( length(I) == 1 ) {
      sparse.matrix.set.row(C, i, as.sparse.vector(V * Bm[I,]))
    }
  }

  C
}


## This version is slow as heck, and hogs even more memory.

## sparse.matrix.mult.old <- function(A,B) {
##   Am <- sparse.to.dense.matrix(A)
##   Bm <- sparse.to.dense.matrix(B)
##   as.sparse.matrix(Am %*% Bm)
## }


## This version is much more memory-efficient, but 2-3 times
## slower.

sparse.matrix.mult.small <- function(A,B) {

  m <- sparse.matrix.nrows(A)
  n <- sparse.matrix.ncols(B)
  C <- sparse.matrix(m,n)
  
  for ( i in 1:m ) {

    row <- sparse.matrix.row(A,i)

    mapply(
           function(k,v) {

             target <- sparse.matrix.row(B,k)

             mapply(
                    function(j,w) {
                      sparse.matrix.inc(C, i, j, v*w)
                      NULL
                    },
                    sparse.vector.nnzero.indices(target),
                    sparse.vector.nnzero.values(target),
                    SIMPLIFY=FALSE,
                    USE.NAMES=FALSE)
             NULL
           },
           sparse.vector.nnzero.indices(row),  
           sparse.vector.nnzero.values(row),
           SIMPLIFY=FALSE,
           USE.NAMES=FALSE)
  }

  C
}


## This version is even slower than the above, for some
## reason I cannot figure out.

## sparse.matrix.mult.1 <- function(A,B) {

##   m <- sparse.matrix.nrows(A)
##   n <- sparse.matrix.ncols(B)
##   Cm <- matrix(0,m,n)

##   for ( i in 1:m ) {

##     print(paste('i=',i))

##     row <- sparse.matrix.row(A,i)

##     mapply(
##            function(k,v) {
##              target <- sparse.matrix.row(B,k)
##              ind <- sparse.vector.nnzero.indices(target)
##              val <- sparse.vector.nnzero.values(target)
##              Cm[i,ind] <- Cm[i,ind] + v*val
##            },
##            sparse.vector.nnzero.indices(row),  
##            sparse.vector.nnzero.values(row),
##            SIMPLIFY=FALSE,
##            USE.NAMES=FALSE)
##   }

##   as.sparse.matrix(Cm)
## }



## This version is slower than the other memory-efficient
## ones above.

## sparse.matrix.mult <- function(A,B) {
##   m <- sparse.matrix.dim(A)[[1]]
##   n <- sparse.matrix.dim(B)[[2]]
##   C <- sparse.matrix(m,n)

##   for ( i in 1:m ) {
##     v <- sparse.matrix.row(A,i)
##     ind <- v$rep()$indices
##     val <- v$rep()$values
    
##     print(i)

##     for ( j in 1:n ) {
##       s <- 0
##       for ( k in 1:length(ind)) {
##         s <- s + val[[k]] * sparse.matrix.ref(B,ind[[k]],j)
##       }
      
##       if ( s != 0 ) {
##         sparse.matrix.set(C,i,j,s)
##       }
##     }
##   }

##   C
## }


## This version is also very slow.

## sparse.matrix.mult <- function(A,B) {
  
##   Bt <- sparse.matrix.transpose(B)
##   m <- sparse.matrix.nrows(A)
##   n <- sparse.matrix.ncols(B)

##   C <- outer(1:m, 1:n,
##              function(I,J) {
##                mapply(
##                       function(i,j) {
##                         sparse.vector.dot(sparse.matrix.row(A,i),
##                                           sparse.matrix.row(Bt,j))
##                       },
##                       I, J,
##                       SIMPLIFY=FALSE)
##              })

##   as.sparse.matrix(C)
## }




sparse.matrix.vector.mult <- function(A,v) {

  m <- sparse.matrix.nrows(A)
  n <- sparse.matrix.ncols(A)
  w <- vector(mode='numeric', m)

  for ( i in 1:m ) {

    row <- sparse.matrix.row(A,i)
    I <- sparse.vector.nnzero.indices(row)
    V <- sparse.vector.nnzero.values(row)

    if ( length(I) >= 1 ) {
      w[[i]] <- sum(V * v[I])
    }
  }

  w
}


## Matrix powers.  This appears to be the fastest
## implementation possible with the representation...

sparse.matrix.pow <- function(A,k) {
  if ( k > 1 ) {
    if ( k %% 2 == 0 ) {
      B <- sparse.matrix.pow(A,k/2)
      sparse.matrix.mult(B,B)
    }
    else {
      B <- sparse.matrix.pow(A,k-1)
      sparse.matrix.mult(A,B)
    }
  }
  else {
    A
  }
}

sparse.matrix.pow.cost <- function(n) {
  if ( n > 1 ) {
    if ( n %% 2 == 0 ) {
      sparse.matrix.pow.cost(n/2) + 1
    }
    else {
      sparse.matrix.pow.cost(n-1) + 1
    }
  }
  else {
    0
  }
}

## sparse.matrix.pow <- function(A,k) {
##   if ( k > 1 ) {
##     as.sparse.matrix(matrix.pow(sparse.to.dense.matrix(A),k))
##   }
##   else if ( k == 1 ) {
##     A
##   }
##   else if ( A == 0 ) {
##     diag(nrow(A))
##   }
##   else {
##     stop(paste('invalid sparse-matrix power',k))
##   }
## }

## R doesn't have a matrix power function?

matrix.pow <- function(A,k) {
  if ( k > 1 ) {
    if ( k %% 2 == 0 ) {
      B <- matrix.pow(A,k/2)
      B %*% B
    }
    else {
      matrix.pow(A,k-1) %*% A
    }
  }
  else {
    A
  }
}
