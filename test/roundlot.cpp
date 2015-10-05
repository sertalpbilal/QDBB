
#include "portfolio.h"


extern double* pBar;
extern double** DhalfVT;
extern double** VDhalfinv;
extern double* mu;
extern double** Q;
extern int localN;



extern int N_;
extern double Rt_;
extern int k_;
extern int cardinaltype_;
extern double C_;
extern string datafolder_;

extern int PROBLEMCODE; // = 1;


int createRoundlot(MSKtask_t* task)
{
  PROBLEMCODE = 1;
  // if(argc >= 7) {
	// std::cout << "Error (1): Argument number is not correct.";
	// cout << argv[0] << endl;
	// }

  // Create variables
  int N = N_;
  localN = N;
  double R = Rt_;
  double mu_0 = 0.001;
  double C = C_;
  
  int* M = new int[N];
  double* P = new double[N];
  
  // Read data files
  char qname[80];
  char pname[80];
  char muname[80];
  char mname[80];
  strcat(strcpy(qname, "../data/Q/"), datafolder_.c_str());
  strcat(strcpy(pname, "../data/P/"), datafolder_.c_str());
  strcat(strcpy(muname, "../data/mu/"), datafolder_.c_str());
  strcat(strcpy(mname, "../data/M/"), datafolder_.c_str());
  
  //
  printf("%s\n",qname);

  // Read the file
  readDouble2DArray(qname, N, Q);
  readDoubleArray(pname, N, P);
  readDoubleArray(muname, N, mu);
  readIntegerArray(mname, N, M);
  
  // Compute parameters
  double* phat = new double[N];
  for(int i=0; i<N; ++i) {
    phat[i] = P[i]*M[i]/C;
  }
  double* mpbar = new double[N];
  for(int i=0; i<N; ++i) {
    mpbar[i] = phat[i]*(mu_0-mu[i]);
  }
  
  double** Q_hat = new double*[N+1];
  for(int i = 0; i < N+1; ++i) {
    Q_hat[i] = new double[N+1];
  }
  for(int i=0; i<N; ++i) {
    for(int j=0; j<N; ++j) {
      Q_hat[i+1][j+1] = phat[i]*phat[j]*Q[i][j];
    }
    Q_hat[i][0] = 0;
    Q_hat[0][i] = 0;
  }
  
  double pbarzero = mu_0 - R;
  
  double* ct = new double[N+1];
  for(int i = 1; i < N+1; ++i) {
    ct[i] = 0;
  }
  
  ct[0] = 1;
  
  pBar = new double[N+1];
  DhalfVT = new double*[N+1];
  VDhalfinv = new double*[N+1];
    for (int i=0; i<N+1; ++i) {
        DhalfVT[i] = new double[N+1];
        VDhalfinv[i] = new double[N+1];
    }
  
  initEigTransform(N, Q_hat, pBar, DhalfVT, VDhalfinv);
  
  /*cout << "npbar" << endl;
  for (int i=0; i<N+1; ++i) {
		cout << npBar[i] << endl;
    }
  
  cout << "VDHALF" << endl;
  for(int i=0; i<N+1; i++) {
	for(int j=0; j<N+1; j++) {
		cout << VDhalfinv[i][j] << endl;
	}
  }
  cout << "VDHALF" << endl; */
  
  // Set MOSEK Environment
  // MSKenv_t env = NULL;
  // MSKtask_t task = NULL;
  MSKrescodee r;

  // /* Create the mosek environment. */
  // r = MSK_makeenv(&env,NULL);
  // r = MSK_maketask(env,3,N+1,&task);
  
  
  
  // TODO Open this
  //r = MSK_linkfunctotaskstream(*task,MSK_STREAM_LOG,NULL,printtxt);
  //r = MSK_linkfunctotaskstream(*task,MSK_STREAM_ERR,NULL,printtxt);
  
  //double maxtime = 3600;
  
  /* MSK_putcallbackfunc(task,
                          usercallback,
                          (void *) &maxtime); */
  
  r = MSK_appendcons(*task,3);
  r = MSK_appendvars(*task,N+1);
  
  // Give Variable Bounds
  r = MSK_putvarbound(*task,0,MSK_BK_FR,-MSK_INFINITY,+MSK_INFINITY);
  for(int j=0; j<N; j++) {
    r = MSK_putvarbound(*task,
    j+1,           /* Index of variable.*/
    MSK_BK_RA,      /* Bound key.*/
			0,// -0.2/phat[j],      /* Numerical value of lower bound.*/
    1/phat[j]); 
  }
   
  
  // Set Variable Types and Names
  for(int j=0; j<N; j++) {
    r = MSK_putvartype(*task,j+1,MSK_VAR_TYPE_INT);
	stringstream varname;
	varname << "z_" << (j+1);
	r = MSK_putvarname(*task,j+1,varname.str().c_str());
  }
  r = MSK_putvarname(*task,0,"t");
  
  // Give Objective
  r = MSK_putcj(*task,0,1);
  
  // Give Constraints
  // Row 1 - pbar pbarzero
  for(int j=0; j<N; j++) {
    r = MSK_putaij(*task, 1, j+1, mpbar[j]);
  }
  // Row 2 - phat leq 1
  for(int j=0; j<N; j++) {
    r = MSK_putaij(*task, 2, j+1, phat[j]);
  }
  
  // Constraint bounds
  // Row 0
  r = MSK_putconbound(*task, 
  0,           /* Index of constraint.*/ 
  MSK_BK_UP,      /* Bound key.*/ 
  -MSK_INFINITY,      /* Numerical value of lower bound.*/ 
  0);     /* Numerical value of upper bound.*/ 
  // Row 1
  r = MSK_putconbound(*task, 1, MSK_BK_UP, -MSK_INFINITY, pbarzero); 
  // Row 2
  r = MSK_putconbound(*task, 2, MSK_BK_UP, -MSK_INFINITY, 1); 
  
  // Give Quadratic Constraints
  
  // Linear Part of QC
  r = MSK_putaij(*task, 0, 0, -1);
  
  // Q Part of QC
  int* rowindex = new int[(N+1)*N/2];
  int* colindex = new int[(N+1)*N/2];
  double* valindex = new double[(N+1)*N/2];
  int busindex = 0;
  for(int i=1; i<N+1; i++) {
    for(int j=1; j<i+1; j++) {
      rowindex[busindex]=i;
      colindex[busindex]=j;
      valindex[busindex]=2*Q_hat[i][j];
      if (i!=j) {
        //valindex[busindex]-=Q_hat[i][j];
      }
      busindex++;
    }
  }
  
  for(int i=0; i<N+1; i++) {
    for(int j=0; j<N+1; j++) {
      //std::cout << Q_hat[i][j] << "\t";
    }
    //std::cout << "\n";
  }
  
  MSK_putqconk(*task, 0, (N+1)*N/2, rowindex, colindex, valindex);
  
  string solver("MOSEK");
  
  delete[] rowindex;
  delete[] colindex;
  delete[] valindex;
  
  delete[] M;
  delete[] P;
  delete[] phat;
  delete[] mpbar;
  
  for(int i=0; i<N+1; i++) {
    delete[] Q_hat[i];
  }
  delete[] Q_hat;
  delete[] ct;
  
  return 1;
  
}

int deleteRoundlot() {
  
  for(int i=0; i<localN+1; ++i) {
    
    delete[] DhalfVT[i];
    delete[] VDhalfinv[i];
    
  }
  
  delete[] DhalfVT;
  delete[] VDhalfinv;
  delete[] pBar;
  
  
  return 1;
}




