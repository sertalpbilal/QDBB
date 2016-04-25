// Single.cpp is a part of QDBB Portfolio test problems
// It have the DCC for single cardinality constraint, i.e. x_i^2 <= z_i

#include "portfolio.h"

extern double* mu;
extern double** Q;

extern int localN;
extern int k;

extern int N_;
extern double Rt_;
extern int k_;
extern int cardinaltype_;
extern int PROBLEMCODE;// = 3;
extern string datafolder_;
extern double integerTolerance_;
extern int FILEOUTPUT;

extern int cutRule_;

/*  @short Creates the single cardinality constrained portfolio optimization problem
 *  @param[in] task  Problem instance (MOSEK)
 *  @param[in] argv  User input at runtime
 */
int createSingleCardinality(MSKtask_t* task)
{
  PROBLEMCODE = 3; // Yay!
  //printf("Working\n");
  localN = N_; // number of assets, total vars: 2*N+1 [t x1...xn z1...zn]
  int N = N_;
  double Rt = Rt_;
  k = k_;

  // Read the file
  char qname[80];
  char muname[80];
  strcat(strcpy(qname, "../data/Q/"), datafolder_.c_str());
  strcat(strcpy(muname, "../data/mu/"), datafolder_.c_str());
 
  readDouble2DArray(qname, N, Q);
  readDoubleArray(muname, N, mu);
  
  MSKrescodee r;
  
  r = MSK_appendcons(*task,1+N+1+1+1); 
  // constraints: 0. return, 1-N x^2<=z, N+1 sum to 1, N+2 x'Qx<=t, and N+3 sum(z) <= k
  r = MSK_appendvars(*task,2*N+1);
  // variables: 0. t, 1-N x, N+1 - 2N z
  if(r) {
    r = r;
  }
  //printf("N: %d\n",N);
  // Give Variable Bounds
  r = MSK_putvarbound(*task,0,MSK_BK_FR,-MSK_INFINITY,+MSK_INFINITY);
  for(int j=0; j<N; j++) {
    r = MSK_putvarbound(*task, j+1, MSK_BK_RA, 0, 1);
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
  
  if(cardinaltype_ == 1) { //Quadratic bound
    // Row 1 to N: x/z relation ---> x_j^2 - z_j <= 0
    int* temprow = new int[1];
    int* tempcol = new int[1];
    double* tempval = new double[1];
    for(int j=0; j<N; j++) {
      // r = MSK_putaij(*task, j+1, j+1, 1);
      MSK_putaij(*task, j+1, N+j+1, -1); // z
      temprow[0] = j+1;
      tempcol[0] = j+1;
      tempval[0] = 2;
      MSK_putqconk(*task, j+1, 1, temprow, tempcol, tempval); // x^2
      MSK_putconbound(*task, j+1, MSK_BK_UP, -MSK_INFINITY, 0); // <= 0 
    }
    delete[] temprow;
    delete[] tempcol;
    delete[] tempval;
  }
  else { // Linear bound  -- cardinaltype_ == 0
    for(int j=0; j<N; j++) {
      r = MSK_putaij(*task, j+1, j+1, 1);
      r = MSK_putaij(*task, j+1, N+j+1, -1);
      r = MSK_putconbound(*task, j+1, MSK_BK_UP, -MSK_INFINITY, 0); 
    }
  }
  
  // Row N+1: sum to 1
  for(int j=0; j<N; j++) {
    r = MSK_putaij(*task, N+1, j+1, 1);
  }
  r = MSK_putconbound(*task, N+1, MSK_BK_RA, 1, 1); 
  
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
  
  // Linear 
  // Cardinality sum(z)<=k
  for(int j=0; j<N; j++) {
    r = MSK_putaij(*task, N+3, N+j+1, 1);
  }
  // RHS
  r = MSK_putconbound(*task, N+3, MSK_BK_UP, -MSK_INFINITY, k); 
  
  //MSK_toconic (*task);
  string solver("MOSEK");
  
  delete[] rowindex;
  delete[] colindex;
  delete[] valindex;
  
  
  
  if(FILEOUTPUT) {
    MSK_writedata(*task, "result/singleOriginal.mps");
  }

  MSKrescodee a =MSK_toconic(*task);
    
  if(FILEOUTPUT) {
    MSK_writedata(*task, "result/singleOriginal2.mps");
  }
  
  return 1;
}

int deleteSingleCardinality() {
  return 1;
}




