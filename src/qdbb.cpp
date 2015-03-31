 
#include <iostream>
#include <limits>
#include <vector>
#include <string>
#include <stdio.h>
#include <stdarg.h>
#include "mosek.h"

struct Node {
    int         ID;
    double      nodeObj;
    double      lowerBound;
    bool        feasible;
    MSKtask_t*  problem;
    Node*       parent;
    std::vector< std::vector <int> > usedCuts;
    int depth;
};

MSKtask_t* oProblem = NULL; // Original problem instance
double globalUpperBound_ = std::numeric_limits<double>::infinity();
int branchingRule_ = 0; // 0: most fractional, 1: index-based, 2: value-based
int cutPriority_ = 0; // 0: default, most fractional, 1: best-value
int cutRule_ = 0; // 0: default, no cut, 1: always cut, 2: deep-heuristic-cut
int cutLimit_ = 1; // max number of cuts to add, default is 1
int nodesProcessed_ = 0;
int totalNodes_ = 0;
int numVars_ = 0;
int N = 0; /* number of assets, always 1 less than numVars */
Node* root;
std::vector< Node* > nodeList_;
int printLevel_ = 5;

int startBB();
int createProblem(MSKtask_t* originalProblem);
int createNewNode(Node* parent, Node** newNode, int varID, double bound, int lower);
int printTxt(int level, const char* fmt, ...);
int selectNode(Node** activeNode);

int main(int argc, char* argv[]) {

    std::cout << "Working!" << std::endl;
    root = NULL;
    // if(argc) ....
    numVars_ = 11;
    N = numVars_ - 1;
    
    MSKenv_t env = NULL;
    MSKtask_t task = NULL;
    MSKrescodee r;

    /* Create the mosek environment. */
    r = MSK_makeenv(&env,NULL);
    r = MSK_maketask(env,3,N+1,&task);
    
    oProblem = &task;
    
    startBB();
    
    return 1;
}
 
int startBB() {
    int status = 0;
    createProblem(oProblem);
    
    root = (Node *) malloc(sizeof(Node));
    createNewNode(0, &root, -1 /*var id*/, 0 /* bound */, 0 /* lower */);
    
    while(true) {
        
        printTxt(2, "Active Nodes: %d", nodeList_.size());
        
        if(nodeList_.size() == 0) {
            break;
        }
        
        Node* activeNode;
        selectNode(&activeNode);
        
        solveLP(activeNode);
        
        if(true) {
            break;
        }
    }
    
    return status;
}

int branch() {
    int status = 0;

    // Call newNode twice
    
    // Add them to active list
    
    
    return status;
}

int cut() {
    int status = 0;
    
    
    
    
    
    return status;
}

int selectNode(Node** activeNode) {
    int status = 0;
    int selected = 0;
    int max_depth = 0;
    *activeNode = nodeList_[0];
    for(unsigned int i=0; i < nodeList_.size(); i++) {
        if(nodeList_[i]->depth > max_depth) {
            selected = i;
            *activeNode = nodeList_[i];
        }
    }
    nodeList_.erase(nodeList_.begin()+selected);
    
    printTxt(4,"Node selected: %d (%p)", (*activeNode)->ID, *activeNode);
    
    return status;
}

int createNewNode(Node* parent, Node** newNode, int varID, double bound, int lower) {
    int status = 0;
    
 
    if(!parent) {
    
        
    
        for (int i = 0; i < N; ++i) { // Initialize array of used cuts
            std::vector<int> row; // Create an empty row for each asset
            //(*newNode)->usedCuts.push_back(row); // Add the row to the main vector
        }
        
        
        (*newNode)->ID = 0;
        (*newNode)->nodeObj = 0;
        (*newNode)->lowerBound = 0;
        (*newNode)->feasible = false;
        (*newNode)->problem = oProblem;
        (*newNode)->parent = 0;
        //(*newNode)->usedCuts = usedCuts;
        (*newNode)->depth = 0;
        
        printTxt(3, "Root is created");
        
    } else {
        
        printTxt(3, "New node is created");
        
        // Add branch constraint here
        // Change problem
        MSKtask_t newProblem = NULL;
        MSK_clonetask((parent->problem), &newProblem);
        if(varID!=-1) {
            MSK_chgbound ( 
                newProblem,         //MSKtask_t    task, 
                MSK_ACC_VAR,        //MSKaccmodee  accmode, 
                varID,              //MSKint32t    i, 
                lower,              //MSKint32t    lower, 
                1,                  //MSKint32t    finite, 
                bound               //MSKrealt     value); 
                );
        }
        
        std::vector< std::vector <int> > usedCuts;
        // copy usedCuts
        for (int i = 0; i < N; i++) { // Initialize array of used cuts
            std::vector<int> row = parent->usedCuts[i]; // Create an empty row for each asset
            usedCuts.push_back(row); // Add the row to the main vector
        }
        
        (*newNode)->ID = totalNodes_;
        (*newNode)->nodeObj = 0;
        (*newNode)->lowerBound = parent->lowerBound;
        (*newNode)->feasible = false;
        (*newNode)->problem = &newProblem;
        (*newNode)->parent = parent;
        (*newNode)->usedCuts = usedCuts;
        (*newNode)->depth = parent->depth+1;
        
    
    }
    
    
    nodeList_.push_back(*newNode);
    totalNodes_++;
    
    
    return status;
}


int solveLP(Node* aNode) {
    int status = 0;
    
    MSKvariabletypee* vartype = malloc(numVars_*sizeof(MSKvariabletypee));
    
    
    
    
    return status;
    
}

int printTxt(int level, const char* fmt, ...) {

    if(level<=printLevel_) {
        std::cout << level << ":";
        for(int i=0; i<level; ++i) {
            std::cout << "   ";
        }
        va_list args;
        va_start(args,fmt);
        vprintf(fmt,args);
        va_end(args);
        std::cout << std::endl;
    }
    return 1;
    
}



