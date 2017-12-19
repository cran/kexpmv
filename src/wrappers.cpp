extern"C" {
void wrapalldmexpv_(int* n,int* m,double* t,double* v,double* w,double* tol,
	double* anorm,double* wsp,int* lwsp,int* iwsp,int* liwsp,int* itrace,int* iflag,
		int* ia,int* ja,double*a ,int* nz,double* res,int* mxstep,int* flag,double* flag2,double* flag3);

void wrapalldgexpv_(int* n,int* m,double* t,double* v,double* w,double* tol,
	double* anorm,double* wsp,int* lwsp,int* iwsp,int* liwsp,int* itrace,int* iflag,
		int* ia,int* ja,double* a,int* nz,double* res,int* mxstep,int* flag,double* flag2,double* flag3);

void myDMEXPV_(int* n,int* m,double* t,double* v,double* w,double* tol,
	double* anorm,double* wsp,int* lwsp,int* iwsp,int* liwsp,int* itrace,int* iflag,
		int* ia,int* ja,double* a,int* nz,int* mxstep,int* flag );

void myDGEXPV_(int* n,int* m,double* t,double* v,double* w,double* tol,
	double* anorm,double* wsp,int* lwsp,int* iwsp,int* liwsp,int* itrace,int* iflag,
		int* ia,int* ja,double* a,int* nz,int* mxstep,int* flag );
		
void wrapdgpadm_(int* ideg,int* m,double* t,double* H,int* ldh,
	double* wsp,int* lwsp,int* ipiv,int* iexph,int *ns,int *iflag ); 

} 
