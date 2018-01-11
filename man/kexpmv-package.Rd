\name{kexpmv-package}
\docType{package}
\title{Matrix exponentiation using Krylov subspace routines}
\description{Matrix Exponential using Krylov subspace routines}
\details{
 This package utilises some of the matrix exponential routines included in
 EXPOKIT (\url{http://www.maths.uq.edu.au/expokit/}), which is software 
 designed to calculate matrix exponentials for both small dense and large sparse
 matrices. The use of Krylov subspace methods implemented within these routines 
 should result in more efficient computations. \cr

 The EXPOKIT software was developed by Roger B. Sidje. Nicholas J. Matzke went on to
 adapt this software for the R package Rexpokit, which enables users to make 
 use of EXPOKIT, a FORTRAN library, within an R environment. Kexpmv uses the 
 foundations of Rexpokit with the aim of making computing matrix exponentials more
 efficient, especially for large sparse matrices which can often be computationally 
 complex.\cr

 Permission to distribute the EXPOKIT source under GPL was obtained from 
 Roger B. Sidje.\cr
 
 EXPOKIT includes functions for exponentiating both small dense matrices and
 large sparse matrices. Fast and efficient matrix exponentiation is needed
 for different application areas, i.e. Markov models. This package allows the
 user to calculate both the matrix exponential in isolation as well as the 
 matrix exponential with the product of a vector, which is essential for multi-state
 Markov models. Both \code{expokit_dgexpv} and \code{expokit_dmexpv} functions can compute both 
 of these calculations with the use of the \code{vector} arguement.\cr

 When the matrix has large dimensions and is sparse in nature, this means the 
 matrix has a high volume of elements equal to zero. Similiar to EXPOKIT this package
 transforms the matrix into Compressed Row Storage (CRS) format before the
 matrix exponentiation is performed. See functions \code{mat2crs} and \code{crs2mat}
 for more details.\cr
 
 \bold{Acknowledgements/sources}\cr \cr \bold{1.} Nicholas
  Matzke \email{nickmatzke.ncse@gmail.com} helped
  greatly with the initial setup of the package.  See his
  \code{Rexpokit} for another R
  implementation of EXPOKIT routines.\cr \cr \bold{2.}
  EXPOKIT, original FORTRAN package, by Roger B. Sidje
  \email{rbs@maths.uq.edu.au}, Department of Mathematics,
  University of Queensland, Brisbane, QLD-4072, Australia,
  (c) 1996-2013 All Rights Reserved. Sidje has given
  permission to include EXPOKIT code in this R package
  under the usual GPL license for CRAN R packages. For the
  full EXPOKIT copyright and license, see
  \code{expokit_copyright.txt} under
  \code{inst/notes}. \cr \cr \bold{3.} The development of this
  package was helped by advice and discussions with my PhD
  supervisors Adele Marshall and Karen Cairns.}
 
\examples{
library(kexpmv)

# Construct a square (nxn) matrix
mat=matrix(c(-0.071207, 0.065573, 0.005634, 0, -0.041206, 0.041206, 0, 0, 0), 
nrow=3, byrow=TRUE) 
 
# Define value for t 
t = 15
 
# Exponentiate this matrix using DGPADM for t = 15 
OutputMat = expokit_dgpadm( mat = mat, t = t, transpose_needed = TRUE)
print(OutputMat)

# Construct column (nx1) vector
v = matrix(0,3,1)
v[1] = 1
 	
# Exponentiate the matrix using DMEXPV for t = 15 with the product of vector v.
OutputMat = expokit_dmexpv( mat = mat, t = t, vector = v, transpose_needed = TRUE,
transform_to_crs = TRUE)

# Print corresponding (nx1) results vector
print(OutputMat$output_probs)

# Print message to determine whether the mxstep value needs increased. If NULL
# then mxstep value is valid.
print(OutputMat$message)

# Exponentiate the matrix using DGEXPV for t = 15.
OutputMat = expokit_dgexpv( mat = mat, t = t, vector = NULL, transpose_needed = TRUE,
transform_to_crs = TRUE)
print(OutputMat$output_mat)
print(OutputMat$message)

# Functions DMEXPV and DGEXPV are very similiar, with the only difference being DMEXPV
# carries out an additional check for Markov chains.}












 

 