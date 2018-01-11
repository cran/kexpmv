# kexpmv

This package utilises some of the matrix exponential routines included in
EXPOKIT (\url{http://www.maths.uq.edu.au/expokit/}), which is software 
designed to calculate matrix exponentials for small dense or large sparse
matrices. This package includes functions to calculate both the matrix
exponential in isolation as well as the product of the matrix exponential
with a vector. 

When the matrix is sparse in nature as well as having large dimensions, this 
means many of the elements in the matrix are zero. Similar to EXPOKIT 
this package transforms the matrices into Compressed Row Storage (CRS) 
format before the matrix exponentiation is performed.



