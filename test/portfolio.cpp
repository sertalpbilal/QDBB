
#include "portfolio.h"

#define printText(flag, ...) if (flag <= OUTLEV) { printf("%d:%*s",flag,2*flag,""); printf(__VA_ARGS__); printf("\n"); }

extern int OUTLEV;

double* pBar;
double** DhalfVT; 
double** VDhalfinv;

double* mu;
double** Q;

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
extern double integerTolerance_;

// This value is overwritten inside createProblem
int PROBLEMCODE = 1; // 1 for roundlot
                     // 2 for cardinality
                     // 3 for single
                     // 4 for diversification

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
    if(strcmp(tmp, "-t")==0 || strcmp(tmp, "-type")==0) { // Problem type
      if(strcmp(argv[i+1],"roundlot")==0)
	PROBLEMCODE = 1;
      else if(strcmp(argv[i+1],"cardinality")==0)
	PROBLEMCODE = 2;
      else if(strcmp(argv[i+1],"single")==0)
	PROBLEMCODE = 3;
      else if(strcmp(argv[i+1],"diverse")==0)
	PROBLEMCODE = 4;
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

  // create mu and Q here
  Q = new double*[N_];
  for(int i = 0; i < N_; ++i) {
    Q[i] = new double[N_];
  }
  mu = new double[N_];


  //printf("PROBLEMCODE: %d\n", PROBLEMCODE);
  if(PROBLEMCODE==1) {
    createRoundlot(task);
  }
  else if(PROBLEMCODE==2) {
    createCardinality(task);
  }
  else if(PROBLEMCODE==3) {
    createSingleCardinality(task);
  }
  else if(PROBLEMCODE==4) {
    createDiverse(task);
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
  else if(PROBLEMCODE==3) {
    deleteSingleCardinality();
  }
  else if(PROBLEMCODE==4) {
    deleteDiverse();
  }

  for(int i=0; i<N_; i++) {
    delete[] Q[i];
  }
  delete[] Q;
  delete[] mu;
  
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

int addNewCut(MSKtask_t env, int asset, double value, int option) { //, double* pBar,  double** DhalfVT,  double** VDhalfinv) {
  
  //printf("Add new cut PROBLEMCODE: %d\n",PROBLEMCODE);
  
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
      if(tvalue+newro < 1e-2) { // TODO Make this a parameter
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
        printText(7,"Cut is NOT deep: %f",tvalue+newro);
        return 0;
      }
      else {
	printText(7,"Cut is deep: %f", tvalue+newro);
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
        //if(i==j) { // diagonal
	//  valindex[busindex] += 2*0.000001; // TODO DEBUGGING Q NOT PSD ERROR
	//}
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
    //printf("N: %d\n", N);

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
    //printf("\ntau: %Lf, k: %d\n",tau, k);
    
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
  else if(PROBLEMCODE == 3) { // SINGLE CARDINAL
    N = numVars_-1;
    //printf("N: %d\n",N);
    int* zrowindex = new int[2];
    int* zcolindex = new int[2];
    double*  zvalindex = new double[2];
    zrowindex[0] = asset;
    zrowindex[1] = N+asset;
    zcolindex[0] = asset;
    zcolindex[1] = N+asset;
    zvalindex[0] = 2;
    zvalindex[1] = -2;
    
    //printf("Working?\n");
    int cutidx;
    MSK_getnumcon(env,&cutidx);
    MSK_appendcons(env,1);
    MSK_putqconk(env, cutidx, 2, zrowindex, zcolindex, zvalindex);
    MSK_putconbound(env, cutidx, MSK_BK_UP, -MSK_INFINITY, 0); 
    

    MSKrescodee r;
    r = MSK_toconic(env);
    r = r;
    //printf("r: %d\n",r);

    char txbuffer[80];
    
    if(FILEOUTPUT)
      sprintf(txbuffer, "result/afterCut%d.mps", asset);
    
    
    if(FILEOUTPUT)
      MSK_writedata(env, txbuffer);
    
    delete[] zrowindex;
    delete[] zcolindex;
    delete[] zvalindex;
    
    return 1;
    
  }
  else if(PROBLEMCODE==4) { // Diversification
    
    int nofAsset = 3;
    int bAsset = 0;
    
    if(asset==1 || asset==2 || asset==3) { // Germany
      bAsset = 1;
    }
    else if(asset==4 || asset==5 || asset ==6) { // France
      bAsset = 4;
    }
    else if(asset==7 || asset==8 || asset ==9) { // Japan
      bAsset = 7;
    }
    else if(asset==10 || asset==11 || asset ==12) { // UK
      bAsset = 10;
    }
    else if(asset==13 || asset==14) { // US
      bAsset = 13;
      nofAsset = 2;
    }
    else if(asset==15 || asset==16 || asset ==17) { // Canada
      bAsset = 15;
    }
    else if(asset==18 || asset==19 || asset ==20) { // Australia
      bAsset = 18;
    }

    int* zrowindex = new int[nofAsset];
    int* zcolindex = new int[nofAsset];
    double*  zvalindex = new double[nofAsset];
    for(int i=0; i<nofAsset; ++i) {
      zrowindex[i] = N+bAsset+i;
      zcolindex[i] = N+bAsset+i;
      if(bAsset+i == asset) { // if selected asset
	zvalindex[i] = -2;
      }
      else {
	zvalindex[i] = 2;
      }
    }

    int cutidx;
    MSK_getnumcon(env,&cutidx);
    MSK_appendcons(env,1);
    MSK_putqconk(env, cutidx, nofAsset, zrowindex, zcolindex, zvalindex);
    MSK_putaij(env, cutidx, N+asset, 2);
    MSK_putconbound(env, cutidx, MSK_BK_UP, -MSK_INFINITY, 1); 
    
    MSK_toconic(env);

    delete[] zrowindex;
    delete[] zcolindex;
    delete[] zvalindex;
    
    
    return 1;
  }
  
  

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

int nextCut(int N, int heuType, double* soln, vector< vector<int> > *usedCuts) {
  /** In this part we will check the solution
  *  and by checking usedCuts, we will decide which asset to
  *  write a cut for
  */
  
  int status = -1;

  double* diff = new double[N];
  int i=0;

  if(heuType==0) { // most fractional

    // Step 1: Sort assets based on their distance to half (.5)
    
    for(i=0; i<N; ++i) {
      diff[i] = fabs(soln[i]-round(soln[i]));
    }

  } else if(heuType==1) { // highest cost
    // Step 1: Sort assets based on their return values
    for(i=0; i<N; ++i) {
      diff[i] = mu[i];
    }
  } 
  else if(heuType==2) { // random
    for(i=0; i<N; ++i) {
      diff[i] = rand() % 300;
    }
  }
  else if(heuType==3) { //bonami dynamic
    for(int i=0; i<N; ++i) {
      double delta_lb = (soln[i]-floor(soln[i]))*(soln[i]-floor(soln[i]))*Q[i][i];
      double delta_ub = (ceil(soln[i])-soln[i])*(ceil(soln[i])-soln[i])*Q[i][i];
      double lower = 0;
      double upper = 0;
      if(delta_lb < delta_ub) {
	lower = delta_lb;
	upper = delta_ub;
      } else {
	lower = delta_ub;
	upper = delta_ub;
      }
      diff[i] = 1*lower + 2*upper;
    }
  }
  else if(heuType==4) { //bonami static
    for(i=0; i<N; ++i) {
      diff[i] = Q[i][i];
    }
  }
  
  // Prevent integer cuts!
  for(i=0; i<N; ++i) {
    if( fabs(soln[i]-round(soln[i])) < integerTolerance_) {
      diff[i] = 0;
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
