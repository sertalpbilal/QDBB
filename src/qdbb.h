
#ifndef qdbb_h
#define qdbb_h


using namespace std;
#include <iostream>
#include <limits>
#include <vector>
#include <string>
#include <cmath>
#include <stdio.h>
#include <stdarg.h>
#include <ctime>
#include "mosek.h"
// from portfolio file
#include <fstream>
#include <sstream>
#include <exception>
#include <math.h>
#include <string.h>
#include <unistd.h>

/* 

TODO: Solution may not be feasible (OPTIMAL), if so, return

 */



struct Node {
  //attributes(const string& name, bool value) : name(name), value(value) {}
    int         ID;
    double      nodeObj;
    double      lowerBound;
    double*     nodeSoln;
    double*     mosekSoln;
    bool        feasible;
    bool        intfeasible;
    bool        eliminated;
    bool        processed;
    MSKtask_t   problem;
    Node*       parent;
    std::vector< std::vector <int> > usedCuts;
    std::vector< std::vector <int> > usedBranches;
    int depth;
    int totalCuts;
};

int startBB(int argc, char* argv[]);

// This function should be defined by user
int createProblem(MSKtask_t* originalProblem, int argc, char* argv[]); // In problem-specific file

// This function should be defined by user
int deleteProblem();                                         // In problem-specific file

int parseInfo(int argc, char* argv[]);
int createNewNode(Node* parent, Node** newNode, int varID, double bound, int lower);
int printTxt(int level, const char* fmt, ...);
int selectNode(Node** activeNode);
int solveLP(Node* aNode);
int isIntFeasible(Node* aNode);
int branch(Node* activeNode);
int cut(Node* activeNode);
int eliminateNodes();
int nextCut(int N, int heuType, double* soln, std::vector< std::vector<int> > *usedCuts);
int getCurrentGap(double* curGap, double* lowerBound, int* lowerNode);
// TODO: Generates new cut without adding it to the problem
int generateNewCut(MSKtask_t env, double* currSoln, double nodeObj, int asset, double value, int option, long double* cdepth);

int addNewCut(MSKtask_t env, int* cutIndex, double* currSoln, double nodeObj, int asset, double value, int option, long double* cdepth, int addToProblem);
int finishBB();
int deleteNode(Node* aNode);
int printToFile(Node* aNode);

#endif
