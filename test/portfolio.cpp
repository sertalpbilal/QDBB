/* 
  Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.

  File:      qo1.c

  Purpose: To demonstrate how to solve a quadratic optimization
              problem using the MOSEK API.
*/
using namespace std;
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <exception>
#include <math.h>
#include "mosek.h" /* Include the MOSEK definition file. */
#include <vector>

//#define NUMCON 1   /* Number of constraints.             */
//#define NUMVAR 3   /* Number of variables.               */
//#define NUMANZ 3   /* Number of non-zeros in A.           */
//#define NUMQNZ 4   /* Number of non-zeros in Q.           */

template<class ENV, class PROB>
int solve(ENV env, PROB & prob, int type, double* solution, string &solver);

template<class ENV, class PROB>
int addDCyC(ENV env, PROB & prob, int asset, double value, string &solver, double* pBar, double** DhalfVT, double** VDhalfinv);

void matMult( double**  mA,  double**  mB, double** mC, int sizeN);
void initEigTransform(int N, double** Q_hat, double* pBar, double** DhalfVT, double** VDhalfinv);



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

/* static void MSKAPI printstr(void *handle,
MSKCONST char str[])
{
  printf("%s",str);
} */

static void MSKAPI printtxt(void          *info,
                            MSKCONST char *buffer)
{
  printf("%s",buffer); 
} /* printtxt */


extern "C"
{
	#include <cblas.h>
void dsyev_(char *jobz, char *uplo, int *n, double *a,
	      int *lda, double *w, double *work, int *lwork,
	      int *info);
int dpotri_(char *uplo, int *n, double *a, int *lda , int *info);
}

void readDoubleArray(string _filename, int N, double* target);
void readIntegerArray(string _filename, int N, int* target);
void readDouble2DArray(string _filename, int N, double** target);
int nextCut(int N, double* soln, vector< vector<int> > *usedCuts);

double* pBar;
double** DhalfVT;
double** VDhalfinv;
double* mu;

int createProblem(MSKtask_t* task, char* argv[])
{

  // if(argc >= 7) {
	// std::cout << "Error (1): Argument number is not correct.";
	// cout << argv[0] << endl;
	// }

  // Create variables
  int N = atoi(argv[1]);
  double R = atof(argv[9]); //0.02;
  double mu_0 = 0;
  double C = atof(argv[8]); //100000;

  double **Q = new double*[N];
  for(int i = 0; i < N; ++i) {
    Q[i] = new double[N];
  }
  
  int* M = new int[N];
  double* P = new double[N];
  mu = new double[N];
  
  // Read the file
  readDouble2DArray("../data/Q", N, Q);
  readDoubleArray("../data/P", N, P);
  readDoubleArray("../data/mu", N, mu);
  readIntegerArray("../data/M", N, M);
  
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
  r = MSK_linkfunctotaskstream(*task,MSK_STREAM_LOG,NULL,printtxt);
  r = MSK_linkfunctotaskstream(*task,MSK_STREAM_ERR,NULL,printtxt);
  
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
  
  //cout << task << endl;
  
  // vector< vector<int> > usedCuts; // Test structure for keeping old cuts and values
  // for (int i = 0; i < N; i++) { // Initialize array of used cuts
        // vector<int> row; // Create an empty row for each asset
        // usedCuts.push_back(row); // Add the row to the main vector
    // }
  
  r;
  
  //*originalProblem = &task;
  
  
  return 1;
  
}


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


	for(int i = 0; i < N; ++i) {
		delete[] eigvectors[i];
	}
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
    
	MSK_writedata(env, "afterCut.lp");
    MSK_writedata(env, "afterCut.mps");
    
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
    
	delete[] rowindex;
	delete[] colindex;
	delete[] valindex;
}





int addNewCut(MSKtask_t env, int asset, double value, int option) { //, double* pBar,  double** DhalfVT,  double** VDhalfinv) {

	
  
    
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
	return 0;
      }
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
    
	delete[] rowindex;
	delete[] colindex;
	delete[] valindex;

	return 1;

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
  r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printtxt);
  
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
			return -1; // No cut is feasible
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

	return -1;
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



