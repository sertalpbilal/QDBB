
#ifndef PORTH
#define PORTH

#include "../src/qdbb.h"



void readDoubleArray(string _filename, int N, double* target);
void readIntegerArray(string _filename, int N, int* target);
void readDouble2DArray(string _filename, int N, double** target);
int nextCut(int N, double* soln, vector< vector<int> > *usedCuts);


template<class ENV, class PROB>
int solve(ENV env, PROB & prob, int type, double* solution, string &solver);

template<class ENV, class PROB>
int addDCyC(ENV env, PROB & prob, int asset, double value, string &solver, double* pBar1, double** DhalfVT, double** VDhalfinv);

void matMult( double**  mA,  double**  mB, double** mC, int sizeN);
void initEigTransform(int N, double** Q_hat, double* pBar1, double** DhalfVT, double** VDhalfinv);

// Problem specific
int createRoundlot(MSKtask_t* task);
int createCardinality(MSKtask_t* task);
int deleteRoundlot();
int deleteCardinality();


#endif




