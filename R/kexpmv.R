
#######################################################
# EXPOKIT-RELATED FUNCTIONS
#######################################################

#' EXPOKIT dgpadm matrix exponentiation on a square matrix
#'
#' This function exponentiates a matrix via the EXPOKIT padm function.
#'
#' From EXPOKIT:\cr
#'
#' \code{*     Computes exp(t*H), the matrix exponential of a general matrix in }\cr
#' \code{*     full, using the irreducible rational Pade approximation to the   }\cr
#' \code{*     exponential function exp(x) = r(x) = (+/-)( I + 2*(q(x)/p(x)) ), }\cr
#' \code{*     combined with scaling-and-squaring.                              }\cr
#'
#' This function is recommended when dealing with small dense matrices. However it can also be
#' used for large sparse matrices when the infinity norm is approximately >100.
#'
#' @param mat an input square matrix.
#' @param t time value to exponentiate by.
#' @param transpose_needed If TRUE (default), matrix will be transposed.
#' @export
#' @author Meabh G. McCurdy \email{mmccurdy01@qub.ac.uk}\cr
#' Nicholas J. Matzke \email{nickmatzke.ncse@gmail.com}
#' @examples
#' # Define input matrix to be exponentiated
#' mat=matrix(c(-0.071207, 0.065573, 0.005634, 0, -0.041206, 0.041206, 0, 0, 0), 
#' nrow=3, byrow=TRUE) 
#' 
#' # Define value of t 
#' t=15
#' 
#' # Exponentiate with EXPOKIT's dgpadm
#' 	Pmat = expokit_dgpadm(mat=mat, t=t, transpose_needed=TRUE)
#' 	print(Pmat)
#' 
#' 
expokit_dgpadm <- function(mat=NULL, t=15, transpose_needed=TRUE)
	{
	defaults = '
	mat=NULL
	t = 15
	transpose_needed=TRUE
	'
	# Check if mat is blank
	if (is.null(mat))
		{
		# Default mat
		cat("\nWARNING: expokit_dgpadm() was provided a mat with value NULL.  Example mat provided instead\n")
		mat=matrix(c(-0.071207, 0, 0, 0.065573, -0.041206, 0, 0.005634, 0.014206, 0), nrow=3, byrow=TRUE)
		}
	# Check if t is blank
	if (is.null(t))
		{
		# Default mat
		stop("\nSTOP ERROR: expokit_dgpadm() was provided a t (time or times list) with value NULL.  \n")
		}
		
	# FOR DGPADM
	ideg = as.integer(6)

	# Order (numrows/numcols) of the matrix
	m = as.integer(nrow(mat))

	# output matrix
	res = double(length=m*m)

	# Transpose matrix, if required
	matvec = mat
	if (transpose_needed == TRUE)
		{
		tmatvec = t(matvec)
		H = as.numeric(tmatvec)
		} else {
		H = as.numeric(matvec)
		}
	
	# (ldh,m):(input) argument matrix
	ldh = m

	# lwsp = length of wsp, the workspace
	# wsp(lwsp):(workspace/output) lwsp .ge. 4*m*m+ideg+1
	lwsp = as.integer(4*m*m+ideg+1)
	wsp = double(length=lwsp)
	
	# ipiv(m)   : (workspace)
	ipiv = integer(length=m)
	
	# iexph:(output) number such that wsp(iexph) points to exp(tH)
	# i.e., exp(tH) is located at wsp(iexph ... iexph+m*m-1)
	iexph = as.integer(0)
	
	# ns:(output) number of scaling-squaring used
	ns = as.integer(0)
	
	# iflag:(output) exit flag
	# 0 - no problem
	# <0 - problem
	iflag = as.integer(0)
	
	# Run the function:
	res <- .C("wrapdgpadm_", as.integer(ideg), as.integer(m), as.double(t), as.double(H), as.integer(ldh), as.double(wsp), as.integer(lwsp), as.integer(ipiv), as.integer(iexph), as.integer(ns), as.integer(iflag))
	
	output = res[[6]]
	output_mat_is = seq(res[[9]], res[[9]]+m*m-1, by=1)
	output_mat = output[output_mat_is]
	output_mat = matrix(output_mat, nrow=m, byrow=FALSE)
	
	return(output_mat)
	}	
	
#' EXPOKIT dmexpv matrix exponentiation on a square matrix
#'
#' This function converts a matrix to CRS format and exponentiates
#' it via the EXPOKIT dmexpv function with the wrapper functions \code{wrapalldmexpv_} 
#' and \code{dmexpv_} around dmexpv. This function can be used to calculate both the 
#' matrix exponential in isolation or the product of the matrix exponential with a 
#' vector. This can be achieved by modifying the \code{vector} variable as shown below.\cr
#'
#' From EXPOKIT:\cr
#' \code{*     The method used is based on Krylov subspace projection}\cr
#' \code{*     techniques and the matrix under consideration interacts only}\cr
#' \code{*     via the external routine 'matvec' performing the matrix-vector}\cr
#' \code{*     product (matrix-free method).}\cr
#' \code{*}\cr
#' \code{*     This is a customised version for Markov Chains. This means that a}\cr
#' \code{*     check is done within this code to ensure that the resulting vector}\cr
#' \code{*     w is a probability vector, i.e., w must have all its components }\cr
#' \code{*     in [0,1], with sum equal to 1. This check is done at some expense}\cr
#' \code{*     and the user may try DGEXPV which is cheaper since it ignores }\cr
#' \code{*     probability constraints.}\cr
#'
#' This check assumes that the transition matrix Q, satisfies Qe = 0 where e is a column 
#' vector of 1's. If this condition does not hold then use the DGEPXV function instead. It 
#' should be noted that both the DMEXPV and DGEXPV functions within EXPOKIT require the 
#' matrix-vector product y = A*x = Q'*x i.e, where A is the transpose of Q. Failure to
#' remember this leads to wrong results.\cr
#'
#' CRS (Compressed Row Storage) format is a compressed format that is
#' required for EXPOKIT's sparse-matrix functions such as DGEXPV and
#' DMEXPV. However this format is not necessary in EXPOKIT's padm-related functions.\cr
#'
#' This function is recommended for large sparse matrices, however the infinity norm of the matrix 
#' proves to be crucial when selecting the most efficient routine. If the infinity norm of the large
#' sparse matrix is approximately >100 may be of benefit to use the \code{expokit_dgpadm} function 
#' for faster computations.
#'
#' @param mat an input square matrix.
#' @param t a time value to exponentiate by.
#' @param vector If NULL (default), the full matrix exponential is returned. However, in 
#' order to fully exploit the efficient speed of EXPOKIT on sparse matrices, this vector argument should
#' be equal to a vector, v. This vector is an n dimensional vector, which in the Markovian case, can
#' represent the starting probabilities.
#' @param transpose_needed If TRUE (default), matrix will be transposed before the exponential operator
#' is performed.
#' @param transform_to_crs If the input matrix is in square format then the matrix will
#' need transformed into CRS format. This is required for EXPOKIT's sparse-matrix functions DMEXPV and 
#' DGEXPV. Default TRUE; if FALSE, then the \code{mat} argument must be a CRS-formatted matrix.
#' @param crs_n If a CRS matrix is input, \code{crs_n} specifies the order (# of rows  or # of columns) of
#' the matrix. Default is NULL.
#' @param anorm The \code{expokit_dmexpv} requires an initial guess at the norm of the matrix. Using the
#' R function \code{\link{norm}} might get slow with large matrices. Default is NULL.
#' @param mxstep The EXPOKIT code performs integration steps and \code{mxstep} is the maximum number of 
#' integration steps performed. Default is 10000 steps; May need to increase this value if function 
#' outputs a warning message.
#' @param tol the tolerance defined for the calculations.
#' @export
#' @author Meabh G. McCurdy \email{mmccurdy01@qub.ac.uk}
#' @examples
#' # Make a square matrix A
#' # Use expokit_dmexpv to calculate both exp(At) and exp(At)v, where t is a 
#' # time value and v is an n dimensional column vector.
#' mat=matrix(c(-0.071207, 0.065573, 0.005634, 0, -0.041206, 0.041206, 0, 0, 0), 
#' nrow=3, byrow=TRUE)
#'
#' # Set the time variable t
#' t=15
#'
#' # Exponentiate with EXPOKIT's dmexpv to obtain the full matrix exponential
#'	OutputMat = expokit_dmexpv(mat=mat, t=t, transpose_needed=TRUE, vector=NULL)
#' 
#' print(OutputMat$output_mat)
#' print(OutputMat$message)
#'
#' # Can also calculate the matrix exponential with the product of a vector.
#' # Create the n dimensional vector
#' v = matrix(0,3,1)
#' v[1] = 1
#'
#' # Exponentiate with EXPOKIT's dmexpv
#'	OutputMat = expokit_dmexpv(mat=mat, t=t, transpose_needed=TRUE, vector=v)
#' 
#' print(OutputMat$output_probs)
#' print(OutputMat$message)
#'
#' # If message is 'NULL' then no error has occured and the number of 
#' # mxsteps defined in the function is acceptable.
#'
expokit_dmexpv <- function(mat=NULL, t=15, vector=NULL, transpose_needed=TRUE, transform_to_crs=TRUE, crs_n=NULL, anorm=NULL, mxstep=10000, tol=as.numeric(0.0000000001))
	{
	defaults = '
	mat=NULL
	t = 15
	vector=NULL
	transpose_needed=TRUE
	transform_to_crs=TRUE
	crs_n=NULL
	anorm=NULL
	mxstep=10000
	tol=0.0000000001
	'
	
	matvec =mat
	mxstep=mxstep
	tol=tol
	
	# Check if mat is blank
	if (is.null(matvec))
		{
		# Default mat
		cat("\nWARNING: expokit_dmexpv() was provided a mat with value NULL.  Example mat provided instead\n")
		matvec = matrix(c(-0.071207, 0, 0, 0.065573, -0.041206, 0, 0.005634, 0.014206, 0), nrow=3, byrow=TRUE)
		}
	# Check if t is blank
	if (is.null(t))
		{
		# Default mat
		stop("\nSTOP ERROR: expokit_dmexpv() was provided a t (time or times list) with value NULL.  \n")
		}	
	
	# Count the number of NON-zeros (nz)
	# and input the matrix size
	if (transform_to_crs == TRUE)
		{
		
		# number of rows in matvec	
		n=nrow(matvec)
		# number of non-zeros
		nz  = sum(matvec != 0)
		
		# Set up vectors needed for CRS format
		ia  = integer(length=n+1)
		ja  = integer(length=nz)
		a   = double(length=nz)	
		
		} else {
		n = crs_n

		# number of non-zeros
		nz  = length(matvec[[3]])
		}
		
	# ideg = degree of polynomial, 6 is usually considered sufficient
	ideg = as.integer(6)
	
	# dimension of Krylov subspace, m=n-1 or m=30 if n>31
	if (n > 31){
		m=30
	} else {
		m=n-1
	}
	
	# v should have as many elements as n; first element = 1 
	if (is.null(vector))
		{
		v=double(n)
		v[1] = 1
		} else {
		v = double(n)
		v = vector
		}
		
	w = double(length=n)
	lwsp = as.integer(n*(m+2)+5*(m+2)^2+ideg+1)
	wsp = double(length=lwsp)
	liwsp = max(m+2, 7)
	iwsp = integer(length=liwsp)
	
	#matvec = matrix(data=Q, nrow=n, byrow=TRUE)
	matvec = matvec
	
	# Tranpose the matrix where required
	if (transform_to_crs == TRUE)
		{
		if (transpose_needed == TRUE)
		{		
		tmatvec = t(matvec)
		} else {
		tmatvec = matvec
		}} else {
		if (transpose_needed == TRUE)
			{
			temp = crs2mat(matvec, n)
			temp2 = t(temp)
			tmatvec = mat2crs(temp2)
			} else {
			tmatvec = matvec
			}
		}
	
	# Calculate the matrix norm
	if ((exists("anorm") == FALSE) || is.null(anorm))
		{
		if (transform_to_crs==FALSE)
		{
				tmpmat1 = crs2mat(tmatvec, n)
				anorm = as.numeric(norm(tmpmat1, type="I"))
				} else {
			anorm = as.numeric(norm(tmatvec, type="I"))
			}
		}
		
	# The itrace flag, if set to 1, results in dmexpv printing some details of
	# the function's run to screen.
	itrace = 0
	iflag = 0
	flag = 0	
	
	# Make the input CRS matrix
	# [1] = row pointer
	# [2] = col index for non zero elements
	# [3] = non zero elements
	
	if (transform_to_crs == TRUE)
		{		
		tmpmat_in_kexpmv_crs_fmt = mat2crs(tmatvec)
		} else {
		tmpmat_in_kexpmv_crs_fmt = tmatvec
		}
		
	# Either way, store the rows/columns in the input variables for FORTRAN
	ia=tmpmat_in_kexpmv_crs_fmt[[1]]
	ja=tmpmat_in_kexpmv_crs_fmt[[2]]
	a=tmpmat_in_kexpmv_crs_fmt[[3]]
	
	# Run the wrapper function	
	if (!is.null(vector)) 
	{
		####################################################################
		# Instead of returning the full matrix exponential, just return the 
		# output probabilities
		###################################################################
		
		# Be sure to input the input probabilities
		v = vector		
		res2 <- .C("mydmexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz), as.integer(mxstep), as.integer(flag))
		output_probs <- matrix(res2[[5]], nrow=n, byrow=FALSE)
		steps <-res2[[19]]
		if (steps == 1) {
		message="Value of t not reached. Need to increase the value of mxstep"
		return(list(message=message, output_probs = output_probs))
		} else{
		return(list(output_probs = output_probs))
		}} else 
	{
		################################################################
		# Return the full matrix exponential
		################################################################
		
		res = double(length=n*n)
		flag2 = double(length=n+1)
		flag3 = double(length=n+1)	
		res2 <- .C("wrapalldmexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz), as.double(res), as.integer(mxstep), as.integer(flag), as.double(flag2), as.double(flag3))
	
		output_mat <- matrix(res2[[18]], nrow=n, byrow=FALSE)
		flag <-res2[[20]]
		if (flag >= 1) {
		message="Value of t not reached. Need to increase the value of mxstep"
		return(list(message=message, output_mat = output_mat))
		} else {
		return(list(output_mat = output_mat))
		}		
	}}

	
#' Convert matrix to CRS format using SparseM function
#'
#' The \code{EXPOKIT}'s \code{expokit_dmexpv} and \code{expokit_dgexpv} functions both deal with sparse matrices.
#' These matrices have a lot of zeros, and can therefore be stored more efficiently by converting the
#' matrix into CRS (Compressed Row Storage) format.\cr
#'
#' In \code{EXPOKIT} and its wrapper functions, a CRS-formatted matrix is input as
#' 3 vectors:\cr
#'
#' ia = row pointer. This vector stores the location in the `a' vector that is the first non-zero element 
#' in a row. \cr
#' ja = column indices for non-zero elements.\cr
#' a = non-zero elements of the matrix.\cr
#' 
#' @param mat A square matrix.
#' @export
#' @author Meabh G. McCurdy \email{mmccurdy01@qub.ac.uk}
#' @examples 
#' # Make a square matrix
#' mat=matrix(c(-0.071207, 0, 0, 0.065573, -0.041206, 0, 0.005634, 0.014206, 0), nrow=3, byrow=TRUE)
#' 
#' # Covert to CRS format
#' mat2crs(mat)
#' print(mat)
#' 
mat2crs <- function(mat)
	{
	defaults = '
	mat = matrix(c(-0.071207, 0, 0, 0.065573, -0.041206, 0, 0.005634, 0.014206, 0), nrow=3, byrow=TRUE)
	'
	
	numrows = nrow(mat)
	numcols = ncol(mat)

	if (numrows != numcols)
		{
		stop("ERROR! mat2crs(mat) says that in mat, nrows != ncols")
		return(NA)
		}

	numcells = numrows ^2
	
	# xy.grid
	
	x = 1:numrows
	y = 1:numrows
	
	# This produces a matrix in crs format
	tmpmat_in_SparseMcrs_fmt = as.matrix.csr(mat)
	tmpcrs = tmpmat_in_SparseMcrs_fmt
  
    # Need the 3 columns: 
  	# [ia] = row pointer
	# [ja] = col index for non zero elements
	# [a] = non zero elements
	ia = tmpcrs@ia
	ja = tmpcrs@ja
	a = tmpcrs@ra
  
	return(list(ia, ja, a))
}

#' Convert a CRS-formatted matrix to standard square format
#'
#' The \code{EXPOKIT}'s DMEPXV and DGEXPV functions both deal with sparse matrices.
#' These matrices have a lot of zeros, and can therefore be stored more efficiently by converting the
#' matrix into CRS (Compressed Row Storage) format.\cr
#'
#' In \code{EXPOKIT} and its wrapper functions, a CRS-formatted matrix is input as
#' 3 vectors: \cr
#'
#' ia = row pointer. This vector stores the location in the `a' vector that is the first non-zero element 
#' in a row. \cr
#' ja = column indices for non-zero elements.\cr
#' a = non-zero elements of the matrix.\cr
#'
#' This function takes a 3-column list CRS matrix and converts it back to standard square format.\cr
#'
#' @param mat a 3 column list including the vectors; ia, ja and a.
#' @param n the dimension of the square (nxn) matrix.
#' @export
#' @author Meabh G. McCurdy \email{mmccurdy01@qub.ac.uk}
#' @examples 
#' # Create a CRS format matrix
#' ia = c(1, 2, 4, 6)
#' ja = c(1, 1, 2, 1, 2)
#' a  = c(-0.071207,  0.065573, -0.041206,  0.005634,  0.041206)
#' crsmat=list(ia, ja, a)
#'
#' # Convert CRS format matrix to square format
#' n = 3
#' mat = crs2mat(crsmat, n)
#' print(mat)
#'
crs2mat <- function(mat, n){
  
	ia = mat[[1]]
	ja = mat[[2]]
	a = mat[[3]]
	
	# Convert crs matrix to square format
	
	nrows=length(ia)-1
	A<-matrix(0,nrows,nrows)                        
	for (i in 1:nrows){                             
		for (j in ia[i]:(ia[i+1]-1)){                   
			A[i,ja[j]]=a[j]                             
		}
	}
	
	# Correct the matrix to account for rows with all elements equal to 0
	
		for (i in 1:n){
		if (ia[i]>(ia[i+1]-1)){
			for (j in 1:n){
			A[i, j]=0
			}
		}
		}
return(A)
}
	
##############################################################################
# 
# NOTE: DGEXPV section.  Same code as dmexpv, but EXPOKIT's DGEXPV should be
# faster than DMEXPV, however DMEXPV runs an accuracy check appropriate for
# Markov chains, which is not done in DGEXPV.
#
##############################################################################

#' EXPOKIT DGEXPV matrix exponentiation on a square matrix
#'
#' This function converts a matrix to CRS format and exponentiates
#' it via the EXPOKIT dgexpv function with the wrapper functions \code{wrapalldgexpv_}  
#' and \code{dgexpv_} around DGEXPV. This function can be used to calculate both the matrix 
#' exponential in isolation or the product of the matrix exponential with a vector. 
#' This can be achieved by modifying the \code{vector} variable as shown below.\cr
#'
#'
#' NOTE: When looking at DGEXPV vs. DMEXPV. According to the EXPOKIT documentation, the
#' DGEXPV routine should be faster than DMEXPV. This is a result of the additional check the DMEXPV
#' routine runs to check whether the output values satisfy the conditions needed for a Markovian 
#' model, which is not done in DGEXPV. \cr
#'
#' From EXPOKIT:\cr
#'
#' \code{*     The method used is based on Krylov subspace projection}\cr
#' \code{*     techniques and the matrix under consideration interacts only}\cr
#' \code{*     via the external routine 'matvec' performing the matrix-vector} \cr
#' \code{*     product (matrix-free method).}\cr
#'
#' It should be noted that both the DMEXPV and DGEXPV functions within EXPOKIT require 
#' the matrix-vector product y = A*x = Q'*x i.e, where A is the transpose of Q. Failure 
#' to remember this leads to wrong results.\cr
#'
#' CRS (Compressed Row Storage) format is a compressed format that is
#' required for EXPOKIT's sparse-matrix functions such as DGEXPV and
#' DMEXPV. However this format is not necessary in EXPOKIT's padm-related functions.\cr 
#' 
#' This function is recommended for large sparse matrices, however the infinity norm of the matrix 
#' proves to be crucial when selecting the most efficient routine. If the infinity norm of the large
#' sparse matrix is approximately >100 may be of benefit to use the \code{expokit_dgpadm} function 
#' for faster computations.
#'
#' @param mat an input square matrix.
#' @param t a time value to exponentiate by.
#' @param vector If NULL (default), the full matrix exponential is returned. However, in 
#' order to fully exploit the efficient speed of EXPOKIT on sparse matrices, this vector argument should
#' be equal to a vector, v. This vector is an n dimensional vector, which in the Markovian case, can
#' represent the starting probabilities.
#' @param transpose_needed If TRUE (default), matrix will be transposed before the exponential operator
#' is performed.
#' @param transform_to_crs If the input matrix is in square format then the matrix will
#' need transformed into CRS format. This is required for EXPOKIT's sparse-matrix functions DMEXPV and 
#' DGEXPV. Default TRUE; if FALSE, then the \code{mat} argument must be a CRS-formatted matrix.
#' @param crs_n If a CRS matrix is input, \code{crs_n} specifies the order (# of rows  or # of columns) of
#' the matrix. Default is NULL.
#' @param anorm The \code{expokit_dgexpv} routine requires an initial guess at the norm of the matrix. Using the
#' R function \code{\link{norm}} might get slow with large matrices. Default is NULL.
#' @param mxstep The EXPOKIT code performs integration steps and \code{mxstep} is the maximum number of 
#' integration steps performed. Default is 10000 steps; May need to increase this value if function 
#' outputs a warning message.
#' @param tol the tolerance defined for the calculations.
#' @export
#' @author Meabh G. McCurdy \email{mmccurdy01@qub.ac.uk}
#' @examples
#' # Make a square (n x n) matrix A
#' # Use expokit_dgexpv to calculate both exp(At) and exp(At)v, where t is a 
#' # time value and v is an n dimensional column vector.
#' mat=matrix(c(-0.071207, 0.065573, 0.005634, 0, -0.041206, 0.041206, 0, 0, 0), 
#' nrow=3, byrow=TRUE)
#'
#' # Set the time variable t
#'  t=15
#'
#' # Exponentiate with EXPOKIT's dgexpv to obtain the full matrix exponential
#'	OutputMat = expokit_dgexpv(mat=mat, t=t, transpose_needed=TRUE, vector=NULL)
#' 
#' print(OutputMat$output_mat)
#' print(OutputMat$message)
#'
#' # Can also calculate the matrix exponential with the product of a vector.
#' # Create the n dimensional vector
#' v = matrix(0,3,1)
#' v[1] = 1
#'
#' # Exponentiate with EXPOKIT's dgexpv
#'	OutputMat = expokit_dgexpv(mat=mat, t=t, transpose_needed=TRUE, vector=v)
#'	
#' print(OutputMat$output_probs)
#' print(OutputMat$message)
#'
#' # If message is 'NULL' then no error has occured and the number of 
#' # mxsteps defined in the function is acceptable.
#'
 
expokit_dgexpv <- function(mat=NULL, t=15, vector=NULL, transpose_needed=TRUE, transform_to_crs=TRUE, crs_n=NULL, anorm=NULL, mxstep=10000, tol=as.numeric(0.0000000001))
	{
	defaults = '
	mat=NULL
	t = 15
	vector=NULL
	transpose_needed=TRUE
	transform_to_crs=TRUE
	crs_n=NULL
	anorm=NULL
	mxstep=10000
	tol=0.0000000001
	'
	
	matvec = mat
	mxstep = mxstep
	tol=tol
	
	# Check if mat is blank
	if (is.null(matvec))
		{
		# Default mat
		cat("\nWARNING: expokit_dgexpv() was provided a mat with value NULL.  Example mat provided instead\n")
		matvec = matrix(c(-0.071207, 0, 0, 0.065573, -0.041206, 0, 0.005634, 0.014206, 0), nrow=3, byrow=TRUE)
		}
	# Check if t is blank
	if (is.null(t))
		{
		# Default mat
		stop("\nSTOP ERROR: expokit_dgexpv() was provided a t (time or times list) with value NULL.  \n")
		}

	
	# Count the number of NON-zeros (nz)
	# and input the matrix size
	if (transform_to_crs == TRUE)
		{
		
		# CRS format
		# number of rows in matvec
		n=nrow(matvec)
		
		# number of non-zeros
		nz  = sum(matvec != 0)
		
		# Set up vectors needed for CRS format
		ia  = integer(length=n+1)
		ja  = integer(length=nz)
		a   = double(length=nz)	
		
		} else {
		n = crs_n
		
		# number of non-zeros
		nz  = length(matvec[[3]])
		}

	# ideg = degree of polynomial, 6 is usually considered sufficient
	ideg = as.integer(6)
	
	# dimension of Krylov subspace, m=n-1 or m=30 if n>31
	if (n > 31){
		m=30
	} else {
		m=n-1
	}
	
	# v should have as many elements as n; first element = 1 (?)
	if (is.null(vector))
		{
		v=double(n)
		v[1] = 1
		} else {
		v = double(n)
		v = vector
		}
	
	w = double(length=n)
	lwsp = as.integer(n*(m+2)+5*(m+2)^2+ideg+1)
	wsp = double(length=lwsp)
	liwsp = max(m+2, 7)
	iwsp = integer(length=liwsp)
	
	#matvec = matrix(data=Q, nrow=n, byrow=TRUE)
	matvec = matvec
	
	# Tranpose the matrix where required
	if (transform_to_crs == TRUE)
		{
		if (transpose_needed == TRUE)
		{		
		tmatvec = t(matvec)
		} else {
		tmatvec = matvec
		}} else {
		if (transpose_needed == TRUE)
			{
			temp = crs2mat(matvec, n)
			temp2 = t(temp)
			tmatvec = mat2crs(temp2)
			} else {
			tmatvec = matvec
			}
		}
			
	# Calculate the matrix norm
	if ((exists("anorm") == FALSE) || is.null(anorm))
		{
		if (transform_to_crs==FALSE)
		{
				tmpmat1 = crs2mat(tmatvec, n)
				anorm = as.numeric(norm(tmpmat1, type="I"))
				} else {
			anorm = as.numeric(norm(tmatvec, type="I"))
			}
		}
		
	# The itrace flag, if set to 1, results in dgexpv printing some details of
	# the function's run to screen.
	itrace = 0
	iflag = 0	
	flag = 0
	
	# Make the input CRS matrix
	# [1] = row pointer
	# [2] = column index of non zero elements
	# [3] = non zero elements
	
	if (transform_to_crs == TRUE)
		{		
		tmpmat_in_kexpmv_crs_fmt = mat2crs(tmatvec)
		} else {
		tmpmat_in_kexpmv_crs_fmt = tmatvec
		}
		
	# Either way, store the rows/columns in the input variables for FORTRAN
	ia = tmpmat_in_kexpmv_crs_fmt[[1]]
	ja = tmpmat_in_kexpmv_crs_fmt[[2]]
	a = tmpmat_in_kexpmv_crs_fmt[[3]]
	
		# Run the wrapper function	
		if (!is.null(vector)) 
		{
		####################################################################
		# Instead of returning the full matrix exponential, just return the 
		# output probabilities
		###################################################################
		
		# Be sure to input the input probabilities
		v = vector		
		res2 <- .C("mydgexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz), as.integer(mxstep), as.integer(flag))
		output_probs <- matrix(res2[[5]], nrow=n, byrow=FALSE)
		steps <-res2[[19]]
		if (steps == 1) {
		message="Value of t not reached. Need to increase the value of mxstep"
		return(list(message=message, output_probs = output_probs))
		} else{
		return(list(output_probs = output_probs))
		}} else 
		{
		################################################################
		# Return the full matrix exponential
		################################################################

		res = double(length=n*n)
		flag2 = double(length=n+1)
		flag3 = double(length=n+1)	
		res2 <- .C("wrapalldgexpv_", as.integer(n), as.integer(m), as.double(t), as.double(v), as.double(w), as.double(tol), as.double(anorm), as.double(wsp), as.integer(lwsp), as.integer(iwsp), as.integer(liwsp), as.integer(itrace), as.integer(iflag), as.integer(ia), as.integer(ja), as.double(a), as.integer(nz), as.double(res), as.integer(mxstep), as.integer(flag), as.double(flag2), as.double(flag3))

		output_mat <- matrix(res2[[18]], nrow=n, byrow=FALSE)
		flag <-res2[[20]]
		if (flag >= 1) {
		message="Value of t not reached. Need to increase the value of mxstep"
		return(list(message=message, output_mat = output_mat))
		} else {
		return(list(output_mat = output_mat))
		}}}	
	
