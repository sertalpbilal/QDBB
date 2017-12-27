
#include "portfolio.h"


extern double* pBar;
extern double** DhalfVT;
extern double** VDhalfinv;
extern double* mu;
extern double** Q;
extern int localN;



extern int N_;
extern double Rt_;
extern double risk_;
extern int k_;
extern int cardinaltype_;
extern int objtype_;
extern double C_;
extern string datafolder_;
extern int FILEOUTPUT;
extern int PROBLEMCODE; // = 1;

double* phat;
double** Q_hat;

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
  double mu_0 = 0;
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
  phat = new double[N]; // a
  for(int i=0; i<N; ++i) {
    phat[i] = P[i]*M[i]/C;
  }
  double* mpbar = new double[N]; // mu_hat
  for(int i=0; i<N; ++i) {
    mpbar[i] = phat[i]*(mu_0-mu[i]);
  }
  
  Q_hat = new double*[N+1];
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
    
  MSKrescodee r;
  
  r = MSK_appendcons(*task,3);
  r = MSK_appendvars(*task,N+1);
  
  // Give Variable Bounds
  r = MSK_putvarbound(*task,0,MSK_BK_LO,0,+MSK_INFINITY);
  for(int j=0; j<N; j++) {
    r = MSK_putvarbound(*task,
			j+1,           /* Index of variable.*/
			MSK_BK_RA,     /* Bound key.*/
			0,             // -0.2/phat[j],      /* Numerical value of lower bound.*/
			1/phat[j]);    // Upper bound (1/a_i)
  }
   
  
  // Set Variable Types and Names
  for(int j=0; j<N; j++) {
    r = MSK_putvartype(*task,j+1,MSK_VAR_TYPE_INT);
    stringstream varname;
    varname << "z_" << (j+1);
    r = MSK_putvarname(*task,j+1,varname.str().c_str());
  }
  r = MSK_putvarname(*task,0,"t");
  r = r;

  if(objtype_ == 0) {
  
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
    // Row 1
    r = MSK_putconbound(*task, 1, MSK_BK_UP, -MSK_INFINITY, pbarzero); 
    // Row 2
    r = MSK_putconbound(*task, 2, MSK_BK_UP, -MSK_INFINITY, 1); 
    
  }
  else if(objtype_ == 1) {
    // Give objective (max: m_hat z)
    for (int j=0; j<N; j++) {
      r = MSK_putcj(*task, j+1, -phat[j]*mu[j]); // min -mu_hat z == max mu_hat z
    }
    
    // Constraint 1: a' z <= 1
    for (int j=0; j<N; j++) {
      r = MSK_putaij(*task, 1, j+1, phat[j]);
    }
    r = MSK_putconbound(*task, 1, MSK_BK_UP, -MSK_INFINITY, 1);

    // Constraint 2: t = sigma
    r = MSK_putaij(*task, 2, 0, 1);
    r = MSK_putconbound(*task, 2, MSK_BK_RA, risk_, risk_);

  }

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

  // Row 0
  r = MSK_putconbound(*task, 
		      0,           /* Index of constraint.*/ 
		      MSK_BK_UP,      /* Bound key.*/ 
		      -MSK_INFINITY,      /* Numerical value of lower bound.*/ 
		      0);     /* Numerical value of upper bound.*/ 
  

  if(FILEOUTPUT) {
    char tbuffer[120];
    char tbufferlp[120];
    //MSK_writedata(*task, "result/RoundlotOriginal.mps");
    sprintf(tbuffer, "result/roundlot_n%2d_d%s.mps", N, datafolder_.c_str()); 
    sprintf(tbuffer, "result/roundlot_n%2d_d%s.lp", N, datafolder_.c_str()); 
    //MSK_writedata(*task, "result/Roundlot_n" + N + "_d" + datafolder_.c_str() + ".mps");
    //MSK_writedata(*task, "result/RoundlotOriginal.lp");
    MSK_writedata(*task, tbuffer);
    MSK_writedata(*task, tbufferlp);
  }
  string solver("MOSEK");
  
  delete[] rowindex;
  delete[] colindex;
  delete[] valindex;
  
  delete[] M;
  delete[] P;
  delete[] mpbar;
  
  delete[] ct;
  
  return 1;
  
}

int deleteRoundlot() {
  
  for(int i=0; i<localN+1; ++i) {
    
    delete[] DhalfVT[i];
    delete[] VDhalfinv[i];
    
  }

  for(int i=0; i<N_+1; i++) {
    delete[] Q_hat[i];
  }
  delete[] Q_hat;
  delete[] phat;

  
  delete[] DhalfVT;
  delete[] VDhalfinv;
  delete[] pBar;
  
  
  return 1;
}

int postSolveRoundlot(double* soln) {
  try {
    int N = N_;
    double totalrisk = 0;
    for(int i=0; i<N; i++) {
      for(int j=0; j<N; j++) {
	totalrisk += soln[i]*Q_hat[i+1][j+1]*soln[j];
      }
    }
    printf("Total Risk  : %.3f = %e", totalrisk, totalrisk);
    if(objtype_ == 1) {
      printf(" < %e Desired risk", risk_);
    }
    printf("\n");
    double totalreturn = 0;
    for(int i=0; i<N; i++) {
      //printf("Mu[%d]: %f, phat[%d]: %f, soln[%d]: %f\n",i,mu[i],i,phat[i],i,soln[i]);
      totalreturn += soln[i]*mu[i]*phat[i];
    }
    printf("Total Return: %.3f = %e", totalreturn, totalreturn);
    if(objtype_ == 0) {
      printf(" > %e Desired return",  Rt_);
    }
    printf("\n");
  }
  catch (...) {
    printf("Exception at post solve tasks...");
  }
  return 0;
}


