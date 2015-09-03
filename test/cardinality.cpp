
#include "portfolio.h"

extern double* mu;

int localN;
extern int k;

int PROBLEMCODE = 2;

/*  @short Creates the cardinality constrained portfolio optimization problem
 *  @param[in] task  Problem instance (MOSEK)
 *  @param[in] argv  User input at runtime
 */
int createProblem(MSKtask_t* task, char* argv[])
{

  // Create variables
  int N = atoi(argv[1]);
  localN = N; // number of assets, total vars: 2*N+1 [t x1...xn z1...zn]
  double Rt = atof(argv[9]); // Desired return rate, e.g. 0.06
  double C = atof(argv[8]); // Capital, e.g. 100000;
  k = atoi(argv[10]);

  int cardinaltype = atoi(argv[11]);

  double **Q = new double*[N];
  for(int i = 0; i < N; ++i) {
    Q[i] = new double[N];
  }
  
  mu = new double[N];
  
  // Read the file
  readDouble2DArray("../data/Q", N, Q);
  readDoubleArray("../data/mu", N, mu);
  
  MSKrescodee r;
  
  r = MSK_appendcons(*task,1+N+1+2); //return, x<=z, sum to 1, x'Qx<=t, and z'z <= k
  r = MSK_appendvars(*task,2*N+1);
  
  // Give Variable Bounds
  r = MSK_putvarbound(*task,0,MSK_BK_LO,0,+MSK_INFINITY);
  for(int j=0; j<N; j++) {
    r = MSK_putvarbound(*task, j+1 /* Index of variable.*/, MSK_BK_RA, 0, 1);
    r = MSK_putvarbound(*task, j+N+1, MSK_BK_RA, 0, 1);
  }
   
  // Set Variable Types and Names
  for(int j=0; j<N; j++) {
    // z
    r = MSK_putvartype(*task,N+j+1,MSK_VAR_TYPE_INT);
    stringstream varname;
    varname << "z_" << (j+1);
    r = MSK_putvarname(*task,j+N+1,varname.str().c_str());
    
    // x
    stringstream xname;
    xname << "x_" << (j+1);
    r = MSK_putvarname(*task,j+1,xname.str().c_str());
  }
  r = MSK_putvarname(*task,0,"t");
  
  // Give Objective
  r = MSK_putcj(*task,0,1);
  
  // Give Constraints
  // Row 0: Return Rate
  for(int j=0; j<N; j++) {
    r = MSK_putaij(*task, 0, j+1, mu[j]);
  }
  // RHS
  r = MSK_putconbound(*task, 0, MSK_BK_LO, Rt, MSK_INFINITY); 
  
  // Row 1 to N: x/z relation
  for(int j=0; j<N; j++) {
    r = MSK_putaij(*task, j+1, j+1, 1);
    r = MSK_putaij(*task, j+1, N+j+1, -1);
    r = MSK_putconbound(*task, j+1, MSK_BK_UP, -MSK_INFINITY, 0); 
  }
  
  // Row N+1: sum to 1
  for(int j=0; j<N; j++) {
    r = MSK_putaij(*task, N+1, j+1, 1);
  }
  r = MSK_putconbound(*task, N+1, MSK_BK_FX, 1, 1); 
  
  // Row N+2: x'Qx <= t
  r = MSK_putaij(*task, N+2, 0, -1); // t
  int* rowindex = new int[(N+1)*N/2];
  int* colindex = new int[(N+1)*N/2];
  double* valindex = new double[(N+1)*N/2];
  int busindex = 0;
  for(int i=1; i<N+1; i++) {
    for(int j=1; j<i+1; j++) {
      rowindex[busindex]=i;
      colindex[busindex]=j;
      valindex[busindex]=2*Q[i-1][j-1];
      busindex++;
    }
  }
  MSK_putqconk(*task, N+2, (N+1)*N/2, rowindex, colindex, valindex);
  r = MSK_putconbound(*task, N+2, MSK_BK_UP, -MSK_INFINITY, 0); 
  
  int* zrowindex = new int[N];
  int* zcolindex = new int[N];
  double* zvalindex = new double[N];


  if(cardinaltype==1) {  // Quadratic cardinality constraint
  
    // Row N+3: z'z <= k
    
    int zbusindex = 0;
    for(int i=1; i<N+1; i++) {
      zrowindex[zbusindex]=N+i;
      zcolindex[zbusindex]=N+i;
      zvalindex[zbusindex]=2;
      zbusindex++;
    }
    MSK_putqconk(*task, N+3, N, zrowindex, zcolindex, zvalindex);
    r = MSK_putconbound(*task, N+3, MSK_BK_UP, -MSK_INFINITY, k);

  }
  else if(cardinaltype==0) {  // Linear cardinality constraint
    // Linear 
    for(int j=0; j<N; j++) {
      r = MSK_putaij(*task, N+3, N+j+1, 1);
    }
    // RHS
    r = MSK_putconbound(*task, N+3, MSK_BK_UP, -MSK_INFINITY, k); 
  }
  else if(cardinaltype==2) { // Conic problem

    MSK_appendvars(*task, 2); // v_1 and v_2
    MSK_appendcons(*task, 2); // v_1 = (k-1)/2  and  v2 = (k+1)/2
    
    MSK_putvarname(*task, 2*N+1, "v1");
    MSK_putvarname(*task, 2*N+2, "v2");

    MSK_putvarbound(*task, 2*N+1, MSK_BK_LO, 0, MSK_INFINITY);
    MSK_putvarbound(*task, 2*N+2, MSK_BK_LO, 0, MSK_INFINITY);

    int* csub = new int[N+2];
    csub[0] = 2*N+2;
    csub[1] = 2*N+1;
    for(int i=0; i<N; i++) { csub[i+2] = N+i+1; }
    MSK_appendcone(*task, MSK_CT_QUAD, 0.0, N+2, csub);

    delete[] csub;

    MSK_putaij(*task, N+4, 2*N+1, 1);
    MSK_putaij(*task, N+5, 2*N+2, 1);
    double copyk = k+0;
    MSK_putconbound(*task, N+4, MSK_BK_RA, (copyk-1)/2, (copyk-1)/2);
    MSK_putconbound(*task, N+5, MSK_BK_RA, (copyk+1)/2, (copyk+1)/2);
    

  }

  MSK_toconic (*task);
  MSK_writedata(*task, "CardinalityOriginal.mps");
  MSK_writedata(*task, "CardinalityOriginal.lp");
  string solver("MOSEK");
  
  delete[] rowindex;
  delete[] colindex;
  delete[] valindex;
  delete[] zrowindex;
  delete[] zcolindex;
  delete[] zvalindex;
  for(int i=0; i<N; i++) {
    delete[] Q[i];
  }
  delete[] Q;
  
  return 1;
}

int deleteProblem() {
  
  delete[] mu;
  
  return 1;
}




