
#include "portfolio.h"

double* pBar;
double** DhalfVT; 
double** VDhalfinv;
double* mu;
double Rt = 0; // Return rate
int k;
int localN;
extern int numVars_;
string datafolder_;

// Problem params
int N_;
double Rt_;
double C_;
int k_;
int cardinaltype_;

// This value is overwritten inside createProblem
int PROBLEMCODE = 1; // 1 for roundlot
                     // 2 for cardinality

extern int FILEOUTPUT;

extern "C"
{
  #include <cblas.h>
  void dsyev_(char *jobz, char *uplo, int *n, double *a,
  int *lda, double *w, double *work, int *lwork,
  int *info);
  int dpotri_(char *uplo, int *n, double *a, int *lda , int *info);
}

int createProblem(MSKtask_t* task, int argc, char* argv[]) {
    
  int i = 0;
  for(i=0; i<argc; i++) {
    char* tmp = argv[i];
    if(argv[i][0]!='-') {
      continue;
    }
    if(strcmp(tmp, "-t")==0) { // Problem type
      if(strcmp(argv[i+1],"roundlot")==0)
	PROBLEMCODE = 1;
      else if(strcmp(argv[i+1],"cardinality")==0)
	PROBLEMCODE = 2;
      else { printf("Unknown Problem Type\n"); exit(0); }
    }

    if(strcmp(tmp, "-d")==0) { // Data file
      datafolder_ = argv[i+1];
    }

    if(strcmp(tmp, "-a")==0) { // Number of assets
      N_ = atoi(argv[i+1]);
    }

    if(strcmp(tmp, "-r")==0) // return 
      Rt_ = atof(argv[i+1]);

    if(strcmp(tmp, "-C")==0) // capital
      C_ = atof(argv[i+1]);

    if(strcmp(tmp, "-k")==0) // cardinality
      k_ = atoi(argv[i+1]);

    if(strcmp(tmp, "-ct")==0) { // cardinality type - quadratic or linear
      if(strcmp(argv[i+1],"quadratic")==0) {
	cardinaltype_ = 1;
      }
      else if(strcmp(argv[i+1],"linear")==0) {
	cardinaltype_ = 0;
      }
    }

  }

  //printf("PROBLEMCODE: %d\n", PROBLEMCODE);
  if(PROBLEMCODE==1) {
    createRoundlot(task);
  }
  else if(PROBLEMCODE==2) {
    createCardinality(task);
  }
  else if(PROBLEMCODE==3) {
    // createSingle(task);
  }

  return 1;
  
}

int deleteProblem() {
  if(PROBLEMCODE==1) {
    deleteRoundlot();
  }
  else if(PROBLEMCODE==2) {
    deleteCardinality();
  }

  return 1;

}




/* Note: This function is declared using MSKAPI,
        so the correct calling convention is
        employed. */
static int MSKAPI usercallback(MSKtask_t            task,
MSKuserhandle_t      handle,
MSKcallbackcodee     caller,
MSKCONST MSKrealt  * douinf,
MSKCONST MSKint32t * intinf,
MSKCONST MSKint64t * lintinf)
{
  
  MSKtask_t *originaltask=(MSKtask_t *) handle;
  
  switch ( caller )
  {
    //case MSK_CALLBACK_BEGIN_INTPNT:
    //  printf("Starting interior-point optimizer\n");
    //  break;
  case MSK_CALLBACK_END_INTPNT:
    printf("Interior-point optimizer finished.\n");
    MSKrescodee r;
    r = MSK_putvarbound(*originaltask,
    3,           /* Index of variable.*/
    MSK_BK_RA,      /* Bound key.*/
    1,      /* Numerical value of lower bound.*/
    2); 
    r = r;
    
    break;
  default:
    break;
  }

  
  return ( 0 );
} /* usercallback */

// create problem is moved to roundlot and cardinality


int restofoldmain() {
  
  
  // /** Cut Management starts here **/
  // double* contsol = new double[N+1];
  // double* contsol2 = new double[N+1];
  
  // double* intsol = new double[N+1];
  // double* intsol2 = new double[N+1];
  
  
  // //Turn -off printing					  
  // MSK_putintparam(task,  MSK_IPAR_LOG,  0);
  
  // solve(task,r,1,intsol, solver); // Integer version
  
  
  // solve(task,r,0,contsol,solver); // Cont. relaxation
  
  
  // int asset = nextCut(N, contsol, &usedCuts);
  
  // addDCyC(task, r, asset+1, contsol[asset+1], solver, npBar, DhalfVT, VDhalfinv);
  
  // // cout << "CP8" << endl; 
  
  // solve(task,r,0,contsol2,solver); // Cont. relaxation
  
  // // cout << "CP9" << endl; 
  
  // solve(task,r,1,intsol2, solver); // Integer version
  
  // cout << "\tMIS\t\tCont\t\tContCut\t\tMISCut" << endl;
  // for(int i=0; i<N+1; i++) {
  // cout << i << "\t" << intsol[i] << "\t\t" << contsol[i] << "\t\t" << contsol2[i] << "\t" << intsol2[i] << endl;
  // }
  
  // return 1;
  
  
  // // Solve the Problem
  // MSKrescodee trmcode;
  
  // MSK_writedata(task, "mymodel.lp");
  
  // /* Run optimizer */
  // r = MSK_optimizetrm(task,&trmcode);

  // /* Print a summary containing information
  // about the solution for debugging purposes*/
  // //MSK_solutionsummary (task,MSK_STREAM_MSG);
  
  // /* Print solutions */
  // MSKsolstae solsta;
  // MSK_getsolsta (task,MSK_SOL_ITR,&solsta); 

  // double* xx = new double[N+1];
  // MSK_getxx(task, MSK_SOL_ITG, xx); 
  // printf("Optimal primal solution\n"); 
  // for(int j=0; j<N+1; ++j) 
  // printf("x[%d]: %e\n",j,xx[j]); 
  
  
  
  
  // MSK_deletetask(&task);
  // MSK_deleteenv(&env);
  
  
  // return (r);
  
  return 1;
} /* main */



template<class ENV, class PROB>
int solve(ENV env, PROB & prob, int type, double* solution, string &solver)
{
  if(solver.compare("MOSEK")==0) {
    solveWithMOSEK((MSKtask_t*) env, prob, solution, type);
  }
  return 1;
}


template<class ENV, class PROB>
int addDCyC(ENV env, PROB & prob, int asset, double value, string &solver, double* pBar, double** DhalfVT, double** VDhalfinv)
{
  if(solver.compare("MOSEK")==0) {
    MOSEK_addDCyC((MSKtask_t*) env, prob, asset, value, pBar, DhalfVT, VDhalfinv);
  }
  return 1;
}



void getEigenvalue(double** matrix, double* eig, double** eigenvectors, int size) {
  char jobz = 'V';
  char uplo = 'L';
  int info = 0;
  int lwork = -1;
  int n=size;
  double* L = new double[n*n];
  
  int i=0, j=0;
  
  for(i=0; i<n; ++i) {
    for(j=0; j<n; ++j) {
      L[i*n+j]=matrix[i][j];
    }
  }
  
  double worksize[1];

  dsyev_(&jobz, &uplo, &n, L, &n, eig, worksize, &lwork, &info);
  lwork = (int) worksize[0];
  double work[lwork];
  dsyev_(&jobz, &uplo, &n, L, &n, eig, work, &lwork, &info);

  for(i=0; i<n; ++i) {
    //cout << "Eig " << i << ": " << eig[i] << endl;
  }
  for(i=0; i<n*n; ++i) {
    //cout << "L( " << i << "): " << L[i] << endl;
  }
  
  /* EIGENVECTOR TRANSPOSE */
  
  for(i=0; i<n; ++i) {
    for(j=0; j<n; ++j) {
      eigenvectors[j][i] = L[i*n+j];
    }
  }
  
  delete[] L;

}





void initEigTransform(int N, double** Q_hat, double* pBar, double** DhalfVT, double** VDhalfinv) {

  int i=0; int j=0;

  



  /** Eigenvalue decomposition **/

  double** new_Q_hat = new double*[N+1];
  for (i=0; i<N+1; ++i) {
    new_Q_hat[i] = new double[N+1];
  }
  for(i = 0; i < N+1; ++i){
    for(j = 0; j < N+1; ++j) {
      new_Q_hat[i][j]=Q_hat[i][j];
    }
  }
  
  double* eig = new double[N+1];
  eig[0] = 1;
  double** eigvectors = new double*[N+1];
  for (i=0; i<N+1; ++i) {
    eigvectors[i] = new double[N+1];
  }
  getEigenvalue(new_Q_hat, eig, eigvectors, N+1);

  // DEBUGGING
  /*cout << "DEBUGGING EIGV" << endl;
  for(i = 0; i < N+1; ++i){
    for(j = 0; j < N+1; ++j) {
      cout << eigvectors[i][j] << "\t";
    }
    cout << endl;
  }*/
  
  // DEBUGGING

  double** D = new double*[N+1];
  for(i=0; i<N+1; ++i) {
    D[i] = new double[N+1];
    for(j=0; j<N+1; ++j) {
      if(i==j) {
        D[i][j] = eig[i];
      } else {
        D[i][j] = 0;
      }
    }
  }
  D[0][0] = 1;

  // Set D half times V transpose
  for(i = 0; i < N+1; ++i){
    for(j = 0; j < N+1; ++j) {
      DhalfVT[i][j]=sqrt(D[i][i])*eigvectors[j][i];
    }
  }

  double** invsqD = new double*[N+1];
  for(i=0; i<N+1; ++i) {
    invsqD[i] = new double[N+1];
  }
  for(i = 0; i < N+1; ++i){
    for(j = 0; j < N+1; ++j) {
      if(i==j) {
        //cout << "D[" << i << "]= " << D[i][i] << endl;
        invsqD[i][i] = 1/sqrt(D[i][i]);
      } else {
        invsqD[i][j] = 0;
      }
    }
  }
  
  // DEBUGGING
  /*  for(i = 0; i < N+1; ++i){
    for(j = 0; j < N+1; ++j) {
      cout << invsqD[i][j] << "\t";
    }
    cout << endl;
  }
    */
  // DEBUGGING

  for(i = 0; i < N+1; ++i){
    for(j = 0; j < N+1; ++j) {
      VDhalfinv[i][j] = eigvectors[i][j]*(1/sqrt(D[j][j])); //  eigvectors[asset+1][i]*invsqD[i][j];
    }
  }

  //pBar = new double[N+1];
  for(i=0; i<N+1; ++i) {
    //pBar[i] =  (invsqD[i][i]*eigvectors[i][0]);
    pBar[i] = invsqD[i][i]*eigvectors[0][i]*(-0.5);
  }
  //pBar[0] = -0.5 * pBar[0];
  
  /* // DEBUGGING
  cout << "Pbar DEBUGGING" << endl;
  for(i = 0; i < N+1; ++i){
    
      cout << pBar[i] << endl;
  }
  */
  // DEBUGGING


  for(int i = 0; i < N+1; ++i) {
    delete[] eigvectors[i];
  }
  for(int i=0; i<N+1; ++i) {
    delete[] new_Q_hat[i];
    delete[] invsqD[i];
    delete[] D[i];
  }
  delete[] invsqD;
  delete[] D;
  delete[] new_Q_hat;
  delete[] eigvectors;
  delete[] eig;

}





void MOSEK_addDCyC(MSKtask_t env, MSKrescodee & problem, int asset, double value, double* pBar,  double** DhalfVT,  double** VDhalfinv) {

  

  
  //MSK_writedata(env, "beforeCut.lp");
  //MSK_writedata(env, "beforeCut.mps");
  
  int N = 0;
  MSK_getnumvar(env,&N);
  
  N = N-1;
  
  /*cout << "pbarinfo" << endl;
  for(int i=0; i<N+1; i++) {
    cout << "pbar[" << i << "]: " << pBar[i] << endl;
  }*/
  
  int i=0; int j=0;

  double alpha0 = floor(value);
  double beta0 = ceil(value);
  
  double* matC = new double[N+1];
  
  for(j=0; j<N+1; ++j) {
    //cout << "cut info" << endl;
    //cout << asset << endl;
    //cout << N << endl;
    matC[j] = VDhalfinv[asset][j];
  }


  double normC = 0;
  for(i=0; i<N+1; ++i) {
    normC = normC + matC[i]*matC[i];
  }
  normC = sqrt(normC);

  double alpha = alpha0/normC;
  double beta = beta0/normC;
  

  double tau = -1;
  
  
  double** newP = new double*[N+1];
  for(i = 0; i < N+1; ++i){
    newP[i] = new double[N+1];
    for(j = 0; j < N+1; ++j) {
      newP[i][j] = tau * (matC[i]/normC) * (matC[j]/normC);
      if (i==j && i!=0) {
        newP[i][j] = newP[i][j] + 1; // Add matrix J effect
      }
      //cout << "P(" << i << "," << j << "): " << newP[i][j] << endl;
    }
  }
  
  //cout << "newp" << endl;
  double* newp = new double[N+1];
  for(i = 0; i < N+1; ++i){
    //cout << "pbar(" << i << "): " << pBar[i] << endl;
    newp[i] = pBar[i] - tau * ((alpha+beta)/2) * (matC[i]/normC);
    //cout << newp[i] << endl;
  }
  
  double newro = tau*alpha*beta;
  
  ///cout << ">> Compact Disjunctive Cut Insertion" << endl;

  double** P1 = new double*[N+1];
  double** VDPDV = new double*[N+1];
  for(i=0; i<N+1; ++i) { // First part of VDPDV
    P1[i] = new double[N+1];
    VDPDV[i] = new double[N+1];
    for(j=0; j<N+1; ++j) {
      P1[i][j] = 0;
      VDPDV[i][j] = 0;
    }
  }
  

  double** DVTT = new double*[N+1];
  for(i=0; i<N+1; ++i) {
    DVTT[i] = new double[N+1];
    for(j=0; j<N+1; ++j) {
      //cout << i << "," << j << endl;
      DVTT[i][j] = DhalfVT[j][i];
    }
  }
  

  matMult(DVTT, newP, P1, N+1);
  matMult(P1, DhalfVT, VDPDV, N+1);
  
  
  double* pDV = new double[N+1];

  for(i=0; i<N+1; ++i) {
    pDV[i] = 0;
    for(j=0; j<N+1; ++j) {
      pDV[i] = pDV[i] + newp[j]*DhalfVT[j][i];
    }
  }
  
  int cutidx = 0;
  
  problem = MSK_getnumcon(env,&cutidx); 
  problem = MSK_appendcons(env,1); 
  
  problem = MSK_putconbound(env, 
  cutidx, 
  MSK_BK_UP, 
  -MSK_INFINITY, 
  -newro); 
  
  // Give Quadratic Constraints
  
  // cout << "CP1" << endl;
  
  
  // Linear Part of QC
  for(int j=0; j<N+1; j++) {
    problem = MSK_putaij(env, cutidx, j, 2*pDV[j]);
  }
  
  // cout << "CP2" << endl;  
  
  
  
  // Q Part of QC
  int* rowindex = new int[(N+1)*(N+2)/2];
  int* colindex = new int[(N+1)*(N+2)/2];
  double* valindex = new double[(N+1)*(N+2)/2];
  int busindex = 0;
  for(int i=0; i<N+1; i++) {
    for(int j=0; j<i+1; j++) {
      rowindex[busindex]=i;
      colindex[busindex]=j;
      valindex[busindex]=2*VDPDV[i][j];
      
      busindex++;
    }
    
  }
  // cout << "CP3" << endl;  
  
  
  MSK_putqconk(env, cutidx, (N+1)*(N+2)/2, rowindex, colindex, valindex);
  
  // cout << "CP4" << endl;  
  
  //MSK_writedata(env, "afterCut.lp");
  //MSK_writedata(env, "afterCut.mps");
  
  // CPLEX STUFF
  /* IloExpr quadCut(env);
    for (i = 1; i < N+1; ++i) { // Quadratic Part - Only z variables
      for (j = 1; j < N+1; ++j) {
        quadCut += VDPDV[i][j]*(*var_z)[i-1]*(*var_z)[j-1];
      }
      quadCut += 2*pDV[i]*(*var_z)[i-1];
    }

    for (j= 1; j<N+1; ++j) { // Quadratic Part - t z
      quadCut += VDPDV[0][j]*(*var_t)*(*var_z)[j-1];
      quadCut += VDPDV[j][0]*(*var_t)*(*var_z)[j-1];
    }

    // Quadratic part - t t
    quadCut += VDPDV[0][0]*(*var_t)*(*var_t);

    quadCut += 2*pDV[0]*(*var_t);

    quadCut += 1*newro; */
  
  delete[] matC;
  
  for(int i=0; i<N+1; ++i) {
    delete[] newP[i];
    delete[] P1[i];
    delete[] VDPDV[i];
    delete[] DVTT[i];
  }
  delete[] DVTT;
  delete[] newP;
  delete[] P1;
  delete[] VDPDV;
  delete[] newp;
  delete[] rowindex;
  delete[] colindex;
  delete[] valindex;
}





int addNewCut(MSKtask_t env, int asset, double value, int option) { //, double* pBar,  double** DhalfVT,  double** VDhalfinv) {

  int response = 1;
  
  char tbuffer[80];

  sprintf(tbuffer, "result/beforeCut%d.mps", asset);
  
  //MSK_writedata(env, "beforeCut.lp");
  
  if(FILEOUTPUT)
    MSK_writedata(env, tbuffer); // "beforeCut.mps");
  
  int N = 0;
  MSK_getnumvar(env,&N);
  
  if(PROBLEMCODE==1) { // ROUND LOT
    
    
    N = N-1;
    
    /*cout << "pbarinfo" << endl;
    for(int i=0; i<N+1; i++) {
      cout << "pbar[" << i << "]: " << pBar[i] << endl;
    }*/
    
    int i=0; int j=0;

    double alpha0 = floor(value);
    double beta0 = ceil(value);
    
    double* matC = new double[N+1];
    
    for(j=0; j<N+1; ++j) {
      //cout << "cut info" << endl;
      //cout << asset << endl;
      //cout << N << endl;
      matC[j] = VDhalfinv[asset][j];
    }


    double normC = 0;
    for(i=0; i<N+1; ++i) {
      normC = normC + matC[i]*matC[i];
    }
    normC = sqrt(normC);

    double alpha = alpha0/normC;
    double beta = beta0/normC;
    

    double tau = -1;
    
    
    double** newP = new double*[N+1];
    for(i = 0; i < N+1; ++i){
      newP[i] = new double[N+1];
      for(j = 0; j < N+1; ++j) {
        newP[i][j] = tau * (matC[i]/normC) * (matC[j]/normC);
        if (i==j && i!=0) {
          newP[i][j] = newP[i][j] + 1; // Add matrix J effect
        }
        //cout << "P(" << i << "," << j << "): " << newP[i][j] << endl;
      }
    }
    
    //cout << "newp" << endl;
    double* newp = new double[N+1];
    for(i = 0; i < N+1; ++i){
      //cout << "pbar(" << i << "): " << pBar[i] << endl;
      newp[i] = pBar[i] - tau * ((alpha+beta)/2) * (matC[i]/normC);
      //cout << newp[i] << endl;
    }
    
    double newro = tau*alpha*beta;
    
    //cout << ">> Compact Disjunctive Cut Insertion" << endl;


    double** P1 = new double*[N+1];
    double** VDPDV = new double*[N+1];
    for(i=0; i<N+1; ++i) { // First part of VDPDV
      P1[i] = new double[N+1];
      VDPDV[i] = new double[N+1];
      for(j=0; j<N+1; ++j) {
        P1[i][j] = 0;
        VDPDV[i][j] = 0;
      }
    }
    

    double** DVTT = new double*[N+1];
    for(i=0; i<N+1; ++i) {
      DVTT[i] = new double[N+1];
      for(j=0; j<N+1; ++j) {
        //cout << i << "," << j << endl;
        DVTT[i][j] = DhalfVT[j][i];
      }
    }
    

    matMult(DVTT, newP, P1, N+1);
    matMult(P1, DhalfVT, VDPDV, N+1);
    
    
    double* pDV = new double[N+1];

    for(i=0; i<N+1; ++i) {
      pDV[i] = 0;
      for(j=0; j<N+1; ++j) {
        pDV[i] = pDV[i] + newp[j]*DhalfVT[j][i];
      }
    }
    
    int cutidx = 0;
    
    

    // Calculate if cut is deep
    if(option==1) {

      double* mysolnxx = (double*) malloc((N+1)*sizeof(double));
      MSK_getsolutionslice(env, MSK_SOL_ITR, MSK_SOL_ITEM_XX, 0, N+1, mysolnxx);


      double tvalue = 0;
      for(int i=0; i<N+1; i++) {
        for(int j=0; j<N+1; j++) {
          tvalue = tvalue + VDPDV[i][j]*mysolnxx[i]*mysolnxx[j];
        }
        tvalue = tvalue + 2*pDV[i]*mysolnxx[i];
      }
      //std::cout << "Value is : " << tvalue << ", rho is : " << -newro << std::endl;
      if(abs(tvalue+newro)<1e-8) {
        response = 0;
        free(mysolnxx);
        
        delete[] matC;
        
        for(int i=0; i<N+1; ++i) {
          delete[] newP[i];
          delete[] P1[i];
          delete[] VDPDV[i];
          delete[] DVTT[i];
        }
        delete[] DVTT;
        delete[] newP;
        delete[] P1;
        delete[] VDPDV;
        delete[] newp;
        delete[] pDV;
        
        return 0;
      }
      
      free(mysolnxx);
    }

    
    MSK_getnumcon(env,&cutidx); 
    MSK_appendcons(env,1); 
    
    MSK_putconbound(env, 
    cutidx, 
    MSK_BK_UP, 
    -MSK_INFINITY, 
    -newro); 
    
    // Give Quadratic Constraints
    
    // cout << "CP1" << endl;
    
    
    // Linear Part of QC
    for(int j=0; j<N+1; j++) {
      MSK_putaij(env, cutidx, j, 2*pDV[j]);
    }
    
    // cout << "CP2" << endl;  
    
    
    
    // Q Part of QC
    int* rowindex = new int[(N+1)*(N+2)/2];
    int* colindex = new int[(N+1)*(N+2)/2];
    double* valindex = new double[(N+1)*(N+2)/2];
    int busindex = 0;
    for(int i=0; i<N+1; i++) {
      for(int j=0; j<i+1; j++) {
        rowindex[busindex]=i;
        colindex[busindex]=j;
        valindex[busindex]=2*VDPDV[i][j];
        
        busindex++;
      }
      
    }
    // cout << "CP3" << endl;  
    
    
    MSK_putqconk(env, cutidx, (N+1)*(N+2)/2, rowindex, colindex, valindex);
    
    
    
    
    
    delete[] rowindex;
    delete[] colindex;
    delete[] valindex;
    
    delete[] matC;
    
    for(int i=0; i<N+1; ++i) {
      delete[] newP[i];
      delete[] P1[i];
      delete[] VDPDV[i];
      delete[] DVTT[i];
    }
    delete[] DVTT;
    delete[] newP;
    delete[] P1;
    delete[] VDPDV;
    delete[] newp;
    delete[] pDV;
    
    return 1;
  
  } 
  else if(PROBLEMCODE==2) { // CARDINALITY // TODO each of these functions should go to their own file!
    
    N = numVars_-1; //(N-1)/2; // total number of variables is 2N+1
    printf("N: %d\n", N);

    // For this cut we need to introduce 
    // 1) a new variable
    // 2) A linear map between z_j and variable 
    // 3) And finally cut itself
    

    // DEBUGGING
    // VARIABLE DUMP
    /*
    printf("Variable Dump\n");
    int numvars = 20;
    MSK_getnumvar(env, &numvars);
    printf("NUMVAR: %d\n",numvars);
    double* mysolnxy = (double*) malloc((numvars)*sizeof(double));                
    MSK_getsolutionslice(env, MSK_SOL_ITR, MSK_SOL_ITEM_XX, 0, numvars, mysolnxy);
    char mname[20];
    int maxlen = 0;
    MSK_getmaxnamelen(env, &maxlen);
    for(int i=0; i<numvars; i++) {
      // Print var name, value, lower bound, upper bound
      MSK_getvarname(env, i, 3, mname);
      double xlow = 0;
      double xupp = 0;
      MSKboundkeye *  bk;
      MSK_getvarbound (env, i, bk, &xlow, &xupp);
		       
      printf("%2d %5s, %6.3f, %6.3e, %.3e\n", i, mname, mysolnxy[i],xlow, xupp);
    }
    */
    
    long double kd = k + 0.0;
    long double tau = ( sqrt(1-(1/kd))  - 1 ) * 2 * kd;
    //tau = 0;
    printf("\ntau: %Lf, k: %d\n",tau, k);
    
    int* zrowindex = new int[N];
    int* zcolindex = new int[N];
    double*  zvalindex = new double[N];
    for(int i=0; i<N; i++){zrowindex[i]=0; zcolindex[i]=0; zvalindex[i]=0.0;}
    int zbusindex = 0;
    for(int i=1; i<=N; i++) {
      zrowindex[zbusindex]=N+i;
      zcolindex[zbusindex]=N+i;
      if(i==asset) {
        //printf("\n asset: %d\n",asset);
        zvalindex[zbusindex]=2+2*tau+0.0;
      } else {
        zvalindex[zbusindex]=2;
      }
      zbusindex++;
    }
    
    int cutidx;
    MSK_getnumcon(env,&cutidx);
    MSK_appendcons(env,1);
    
    //printf("CutID: %d, asset: %d\n", cutidx, asset);

    /*
    OLD CONE CODES WERE HERE
    */
    
    
    MSK_putqconk(env, cutidx, N, zrowindex, zcolindex, zvalindex);
    MSK_putaij(env, cutidx, N+asset, - 0.5*tau);
    MSK_putconbound(env, cutidx, MSK_BK_UP, -MSK_INFINITY, k); 
    
    

    

    MSK_toconic(env);

  char txbuffer[80];

  if(FILEOUTPUT)
    sprintf(txbuffer, "result/afterCut%d.mps", asset);


  if(FILEOUTPUT)
    MSK_writedata(env, txbuffer);

  delete[] zrowindex;
    delete[] zcolindex;
    delete[] zvalindex;
    

    //delete[] csub;
    
    
    return 1;
    
    
  }  
  // cout << "CP4" << endl;  
  
  //MSK_writedata(env, "afterCut.lp");
  //MSK_writedata(env, "afterCut.mps");
  
  // CPLEX STUFF
  /* IloExpr quadCut(env);
    for (i = 1; i < N+1; ++i) { // Quadratic Part - Only z variables
      for (j = 1; j < N+1; ++j) {
        quadCut += VDPDV[i][j]*(*var_z)[i-1]*(*var_z)[j-1];
      }
      quadCut += 2*pDV[i]*(*var_z)[i-1];
    }

    for (j= 1; j<N+1; ++j) { // Quadratic Part - t z
      quadCut += VDPDV[0][j]*(*var_t)*(*var_z)[j-1];
      quadCut += VDPDV[j][0]*(*var_t)*(*var_z)[j-1];
    }

    // Quadratic part - t t
    quadCut += VDPDV[0][0]*(*var_t)*(*var_t);

    quadCut += 2*pDV[0]*(*var_t);

    quadCut += 1*newro; */
  
  
  
  

  return response;

}




void matMult( double** mA,  double** mB, double** mC, int sizeN) {
  int i=0, j=0, k=0;
  int n = sizeN;

  for(i=0; i<n; ++i) {
    for(k=0; k<n; ++k) {
      for(j=0; j<n; ++j) {
        mC[i][j] = mC[i][j] + mA[i][k]*mB[k][j];
      }
    }
  }

}


/**
* Generic continuous problem solver for AA
*/
int solveWithMOSEK(MSKtask_t realtask, MSKrescodee & problem, double* solution, int type) {
  
  // MOSEK Interface
  //MSKtask_t task = (MSKtask_t) env;
  //MSKrescodee r = (MSKrescodee) problem;
  
  MSKtask_t task;
  
  //MSKtask_t task;
  //MSK_clonetask(myenv, &task);
  MSK_clonetask(realtask, &task);
  MSKrescodee r(problem);
  // TODO Open this
  //r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printtxt);
  
  //double maxtime = 3600;
  
  MSK_putcallbackfunc(task,
  usercallback,
  (void *) &task);
  
  int N = 0;
  MSK_getnumvar(task,&N);
  
  for(int j=0; j<N-1; j++) {
    if(type==0) {
      r = MSK_putvartype(task,j+1,MSK_VAR_TYPE_CONT);
    } else if(type==1) {
      r = MSK_putvartype(task,j+1,MSK_VAR_TYPE_INT);
    }
  }
  
  cout << r << endl;
  
  MSKrescodee trmcode;
  
  stringstream filename;
  if(type==0)
  filename << "aa_cont_" << (N-1) << ".lp";
  else
  filename << "aa_int_" << (N-1) << ".lp"; 
  

  if(FILEOUTPUT)
    MSK_writedata(task, filename.str().c_str());
  
  // cout << "CP30" << endl;
  try{
    r = MSK_optimizetrm(task,&trmcode);
  } catch(int e) {
    
  }
  
  // cout << "CP31" << endl;
  
  if(type==0) // Continuous QP
  MSK_getxx(task, MSK_SOL_ITR, solution); // Interior Point
  else if(type==1) // Mixed-Integer QP
  MSK_getxx(task, MSK_SOL_ITG, solution); // Integer Solution
  
  // cout << "CP32" << endl;
  
  cout << "Solution: " << endl;
  for(int j=0; j<N+1; j++) {
    cout << solution[j] << endl;
  }
  
  
  return 1;
  
}



int nextCut(int N, int heuType, double* soln, vector< vector<int> > *usedCuts) {
  /** In this part we will check the solution
  *  and by checking usedCuts, we will decide which asset to
  *  write a cut for
  */
  
  int status = -1;

  double* diff = new double[N];
  int i=0;

  if(heuType==0) {

    // Step 1: Sort assets based on their distance to half (.5)
    
    for(i=0; i<N; ++i) {
      diff[i] = fabs(soln[i]-round(soln[i]));
    }

  } else if(heuType==1) {
    // Step 1: Sort assets based on their return values
    for(i=0; i<N; ++i) {
      diff[i] = mu[i];
    }
  }

  // Step 2: Find feasible cuts
  for(i=0; i<N; ++i) { // At most N iteration is needed
    double largest = diff[0];
    int largestIndex = 0;
    int k=0;
    for(k=1; k<N; ++k) { // Step 3: Find the maximum among the remaining candidates
      if(largest<diff[k]) {
        largest = diff[k];
        largestIndex = k;
      }
    }
    if(largest==-1) {
      delete[] diff;
      return -1;
      //return -1; // No cut is feasible
    }

    // Step 4: Check feasibility of cut
    bool feas = true;
    unsigned int j=0;
    for(j=0; j<(*usedCuts)[largestIndex].size(); ++j) {
      if((*usedCuts).at(largestIndex).at(j)==floor(soln[largestIndex+1])) {
        // Not feasible, skip it
        diff[largestIndex] = -1;
        feas=false;
        break;
      }
    }

    if(feas) {
      // Step 5: Found, update cut list
      (*usedCuts).at(largestIndex).push_back(floor(soln[largestIndex+1]));

      //cout << "\n===New cut for asset " << (largestIndex+1) << " for value: " << soln[largestIndex+1] << endl;
      delete[] diff;
      return largestIndex;
    } else {
      // Go on for the next candidate
    }
  }
  
  
  delete[] diff;
  
  return status;
}





void readDoubleArray(string _filename, int N, double* target) {
  ifstream inputFile(_filename.c_str()); //"data/NASDAQ-P-16"); //"data/NASDAQ-mu-16");
  char const row_delim = '\n';
  int i=0;
  string row;
  for (i=0; i<N; ++i) {
    getline(inputFile, row, row_delim);
    target[i] = atof(row.c_str());
  }
}

void readIntegerArray(string _filename, int N, int* target) {
  ifstream inputFile(_filename.c_str()); //"data/M1-P-16");
  char const row_delim = '\n';
  int i=0;
  string row;
  for (i=0; i<N; ++i) {
    getline(inputFile, row, row_delim);
    target[i] = atof(row.c_str());
  }
}

void readDouble2DArray(string _filename, int N, double** target) {
  ifstream inputFile(_filename.c_str()); //"data/NASDAQ-Q-16");
  char const row_delim = '\n';
  char const field_delim = '\t';
  int i=0, j=0;
  string row;
  string col;
  for (i=0; i<N; ++i) {
    getline(inputFile, row, row_delim);
    istringstream ss(row);
    for (j=0; j<N; ++j) {
      getline(ss, col, field_delim);
      target[i][j] = atof(col.c_str());
    }
  }
}


/*
OLD CONE CODES

    // 1) new variable
    int numvar;
    MSK_getnumvar (env, &numvar); 
    MSK_appendvars (env, 1);
    MSK_putvarbound(env, numvar, MSK_BK_LO, 0, MSK_INFINITY);
    MSK_putvarname(env, numvar, "u");
    
    // 2) Linear map
    int numcon;
    MSK_getnumcon(env, &numcon);
    MSK_appendcons(env, 1);
    MSK_putaij(env, numcon, numvar, 1);
    double sqtau = -sqrt(fabs(1+tau));
    MSK_putaij(env, numcon, N+asset, sqtau);
    double linmaprhs = - (sqrt(fabs(1+tau))/(1+tau)) * tau / 2;
    MSK_putconbound(env, numcon, MSK_BK_RA, linmaprhs, linmaprhs);
    
    
    // 3) Conic Constraint DCC 
    int* csub = new int[N];
    int loopindex = 1;
    csub[0] = numvar;
    for(int i=0; i<N; i++) {
      if(i==asset-1) {
	// skip it
      } else {
	csub[loopindex] = N+i+1;
	loopindex++;
      }
    }
    int r = MSK_appendcone(env, MSK_CT_QUAD, 0.0, N, csub);
    
*/



