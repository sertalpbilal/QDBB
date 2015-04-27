 
#include <iostream>
#include <limits>
#include <vector>
#include <string>
#include <cmath>
#include <stdio.h>
#include <stdarg.h>
#include "mosek.h"

/* 

TODO: Solution may not be feasible (OPTIMAL), if so, return

 */



struct Node {
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
};

static void MSKAPI printstr(void *handle, MSKCONST char str[])
{
  printf("%s",str);
}

MSKtask_t oProblem = NULL; // Original problem instance
double globalUpperBound_ = std::numeric_limits<double>::infinity();
double* bestSoln_;
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
std::vector< Node* > allNodeList_; // Keeps pointer of all nodes for destructor
int printLevel_ = 2;

int startBB();
int createProblem(MSKtask_t* originalProblem);
int createNewNode(Node* parent, Node** newNode, int varID, double bound, int lower);
int printTxt(int level, const char* fmt, ...);
int selectNode(Node** activeNode);
int solveLP(Node* aNode);
int isIntFeasible(Node* aNode);
int branch(Node* activeNode);
int eliminateNodes();

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
    
    oProblem = task;
    
    startBB();
    
    return 1;
}
 
int startBB() {
    int status = 0;
    createProblem(&oProblem);
    
    root = (Node *) malloc(sizeof(Node));
    createNewNode(0, &root, -1 /*var id*/, 0 /* bound */, 0 /* lower */);
    
    while(true) {
        
        printTxt(3, "Active Nodes: %d", nodeList_.size());
        
        if(nodeList_.size() == 0) {
            break;
        }
        
        Node* activeNode;
        selectNode(&activeNode);
        
	if(!activeNode->eliminated) {
	  solveLP(activeNode);
	} else {
	  continue;
	}

	if( activeNode->feasible) { isIntFeasible(activeNode); }

	if( activeNode->feasible && activeNode->intfeasible) {
	  if(activeNode->nodeObj < globalUpperBound_) {
	    globalUpperBound_ = activeNode->nodeObj;
	    bestSoln_ = activeNode->nodeSoln;
	    eliminateNodes();
	    printTxt(1, "Best objective value: %f", globalUpperBound_);
	  }
	} else if(activeNode->feasible && !activeNode->intfeasible) {
	  // branch
	  if(activeNode->nodeObj <= globalUpperBound_) {
	    printTxt(3, "Upper bound: %f, current objective: %f, branching...",globalUpperBound_,activeNode->nodeObj);
	    branch(activeNode);
	  }
	} else {
	  continue;
	}
        
    }

    // Summary report
    printTxt(1, "========== SUMMARY ==========");
    printTxt(1, "Best objective value: %f", globalUpperBound_);
    printTxt(1, "Values: ");
    for(int i=0; i<N; i++) {
      printTxt(1, "\t%d: %f", (i+1), bestSoln_[i]);
    }
    printTxt(1,"Done...");
    return status;
}

int branch(Node* activeNode) {
    int status = 0;

    // find most fractional
    int asset = 0;
    double mostfrac = 0;
    for(int i=0; i<N; ++i) {
      double cValue = activeNode->nodeSoln[i];
      double smallest = fmin(cValue - floor(cValue), ceil(cValue) - cValue);
      if(smallest >= mostfrac) {
	asset = i;
	mostfrac = smallest;
      }
    }

    // Call newNode twice

    // Left - Less than bound
    Node* leftNode = (Node*) malloc(sizeof(Node));
    createNewNode(activeNode /*parent*/, &leftNode, asset /*var id*/,  floor(activeNode->nodeSoln[asset])/* bound */, 0 /* lower */);
    // Right - Greater than bound
    Node* rightNode = (Node*) malloc(sizeof(Node));
    createNewNode(activeNode /*parent*/, &rightNode, asset /*var id*/, ceil(activeNode->nodeSoln[asset]) /* bound */, 1 /* lower */);
    
    
    // Add them to active list
    // Done in create
    
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

int eliminateNodes() {
  for(unsigned int i=0; i<nodeList_.size(); i++) {
     if(nodeList_[i]->lowerBound > globalUpperBound_) {
       nodeList_[i]->eliminated = true;
     }
  }
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
        (*newNode)->intfeasible = false;
	(*newNode)->eliminated = false;
        (*newNode)->problem = oProblem;
        (*newNode)->parent = 0;
        //(*newNode)->usedCuts = usedCuts;
        (*newNode)->depth = 0;
        
        printTxt(3, "Root is created");
        
    } else {
        
        printTxt(3, "New node is created");
        
        // Add branch constraint here
        // Change problem
        MSKtask_t newProblem;
        MSK_clonetask((parent->problem), &newProblem);
	
	//MSK_analyzeproblem(parent->problem, MSK_STREAM_LOG);
	
	//MSK_analyzeproblem(newProblem, MSK_STREAM_LOG);
	
	
        if(varID!=-1) {
            MSK_chgbound ( 
                newProblem,         //MSKtask_t    task, 
                MSK_ACC_VAR,        //MSKaccmodee  accmode, 
                varID+1,            //MSKint32t    i, 
                lower,              //MSKint32t    lower, 
                1,                  //MSKint32t    finite, 
                bound               //MSKrealt     value); 
                );
        }
        
        std::vector< std::vector <int> > usedCuts;
        // copy usedCuts
        for (int i = 0; i < N; i++) { // Initialize array of used cuts
	  //std::vector<int> row = parent->usedCuts[i]; // Create an empty row for each asset
            //usedCuts.push_back(row); // Add the row to the main vector
        }
        
        (*newNode)->ID = totalNodes_;
        (*newNode)->nodeObj = 0;
        (*newNode)->lowerBound = parent->lowerBound;
        (*newNode)->feasible = false;
        (*newNode)->intfeasible = false;
	(*newNode)->eliminated = false;
        //(*newNode)->problem = &newProblem;
        (*newNode)->problem = newProblem;
	(*newNode)->parent = parent;
        //(*newNode)->usedCuts = usedCuts;
        (*newNode)->depth = parent->depth+1;
        
	
    }
    
    
    nodeList_.push_back(*newNode);
    allNodeList_.push_back(*newNode);
    totalNodes_++;
    
    
    return status;
}


int solveLP(Node* aNode) {
    int status = 0;
    MSKtask_t mtask = NULL;
    
    //MSK_analyzeproblem(mtask, MSK_STREAM_LOG);


    MSK_clonetask(aNode->problem, &mtask);
    
    //MSKvariabletypee* vartype = malloc(numVars_*sizeof(MSKvariabletypee));
    for(int i=0; i<N; ++i) {
        MSK_putvartype (mtask, i+1, MSK_VAR_TYPE_CONT);
    }
    //std::cout << "Running so far!" << std::endl;
    //MSK_linkfunctotaskstream (mtask, MSK_STREAM_LOG, NULL, printstr);
    MSK_optimize(mtask);
    
    MSKsolstae solsta;
    MSK_getsolsta (mtask, MSK_SOL_ITR, &solsta); 
    
    MSKrealt primalobj;

    //if(solsta) {
        // Save it, etc..
      
      //}

    switch( solsta ) 
      { 
      case MSK_SOL_STA_OPTIMAL: 
      case MSK_SOL_STA_NEAR_OPTIMAL: 
        { 
          aNode->feasible = true;
	  MSK_getprimalobj(mtask, MSK_SOL_ITR, &primalobj);
	  aNode->nodeObj = primalobj;
	  if(primalobj > aNode->lowerBound) {
	    aNode -> lowerBound = primalobj;
	  }
	  aNode -> nodeSoln = (double*) malloc(N*sizeof(double));
	  MSK_getsolutionslice(mtask, MSK_SOL_ITR, MSK_SOL_ITEM_XX, 1, N+1, aNode->nodeSoln);
	   
	  //for(int j=0; j<N; ++j) {
	  //  std::cout << "Variable " << j << " value: " << aNode->nodeSoln[j] << std::endl;
	  //}
	  
	  // std::cout << primalobj << std::endl;
          printTxt(3, "Node (%d) Optimal Objective: %f", aNode->ID, aNode->nodeObj);
          break; 
        } 
      case MSK_SOL_STA_DUAL_INFEAS_CER: 
      case MSK_SOL_STA_PRIM_INFEAS_CER: 
      case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER: 
      case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
	aNode->feasible = false;
	printTxt(3,"Node (%d) infeasibility certificate found.", aNode->ID); 
	break; 
      default: 
	aNode->feasible = false;
	break; 
      } 
    
    //MSK_solutionsummary (mtask, MSK_STREAM_LOG);
    
    
    return status;
    
}

int isIntFeasible(Node* aNode) {
  bool intfeasible = true;
  for(int i=0; i<N; ++i) {
    
    double intpart;
    if( std::modf(aNode->nodeSoln[i], &intpart) >= 1e-8) {
      intfeasible = false;
      printTxt(5,"Integer infeasibility at asset %d, value: %f",i, aNode->nodeSoln[i]);
    }
  }
  
  aNode->intfeasible = intfeasible;
  intfeasible ? printTxt(4,"Node %d is integer feasible",aNode->ID) : printTxt(4,"Node %d is NOT integer feasible",aNode->ID);
    return 1;
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



