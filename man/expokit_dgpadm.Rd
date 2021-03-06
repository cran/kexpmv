% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kexpmv.R
\name{expokit_dgpadm}
\alias{expokit_dgpadm}
\title{EXPOKIT dgpadm matrix exponentiation on a square matrix}
\usage{
expokit_dgpadm(mat = NULL, t = 15, transpose_needed = TRUE)
}
\arguments{
\item{mat}{an input square matrix.}

\item{t}{time value to exponentiate by.}

\item{transpose_needed}{If TRUE (default), matrix will be transposed.}
}
\description{
This function exponentiates a matrix via the EXPOKIT padm function.
}
\details{
From EXPOKIT:\cr

\code{*     Computes exp(t*H), the matrix exponential of a general matrix in }\cr
\code{*     full, using the irreducible rational Pade approximation to the   }\cr
\code{*     exponential function exp(x) = r(x) = (+/-)( I + 2*(q(x)/p(x)) ), }\cr
\code{*     combined with scaling-and-squaring.                              }\cr

This function is recommended when dealing with small dense matrices. However it can also be
used for large sparse matrices when the infinity norm is approximately >100.
}
\examples{
# Define input matrix to be exponentiated
mat=matrix(c(-0.071207, 0.065573, 0.005634, 0, -0.041206, 0.041206, 0, 0, 0), 
nrow=3, byrow=TRUE) 

# Define value of t 
t=15

# Exponentiate with EXPOKIT's dgpadm
	Pmat = expokit_dgpadm(mat=mat, t=t, transpose_needed=TRUE)
	print(Pmat)


}
\author{
Meabh G. McCurdy \email{mmccurdy01@qub.ac.uk}\cr
Nicholas J. Matzke \email{nickmatzke.ncse@gmail.com}
}
