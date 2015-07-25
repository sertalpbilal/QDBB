
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

/* 

TODO: Solution may not be feasible (OPTIMAL), if so, return

 */



struct Node {
  //attributes(const string& name, bool value) : name(name), value(value) {}
    int         ID;
    double      nodeObj;
    double      lowerBound;
    double*     nodeSoln;
    bool        feasible;
    bool        intfeasible;
    bool        eliminated;
    MSKtask_t   problem;
    Node*       parent;
  std::vector< std::vector <int> > usedCuts;
    int depth;
    int totalCuts;
};


 


int startBB(char* argv[]);
int createProblem(MSKtask_t* originalProblem, char* argv[]);
int createNewNode(Node* parent, Node** newNode, int varID, double bound, int lower);
int printTxt(int level, const char* fmt, ...);
int selectNode(Node** activeNode);
int solveLP(Node* aNode);
int isIntFeasible(Node* aNode);
int branch(Node* activeNode);
int cut(Node* activeNode);
int eliminateNodes();
int nextCut(int N, int heuType, double* soln, std::vector< std::vector<int> > *usedCuts);
int addNewCut(MSKtask_t env, int asset, double value, int option);

#endif
