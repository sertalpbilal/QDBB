 
#include "qdbb.h"
 
static void MSKAPI printstr(void *handle, MSKCONST char str[])
{
  printf("%s",str);
}


MSKtask_t oProblem = NULL; // Original problem instance
double globalUpperBound_ = std::numeric_limits<double>::infinity();
double* bestSoln_;
int branchingRule_ = 0; // 0: most fractional, 1: index-based, 2: value-based
int cutPriority_ = 0; // 0: default, most fractional, 1: best-value
int cutRule_ = 1; // 0: default, no cut, 1: always cut, 2: fading-cuts, 3: root-heuristic cut
                  // 4: min depth for cut, 5: only if deep cut
double deepCutThreshold_ = 1e-1; // deep-cut threshold value
int cutSelection_ = 1;
int cutLimit_ = 1; // max number of cuts to add, default is 1
int cutPerIteration_ = 3; // number of cuts to be added at each relaxation
int iterationLimit_ = 1;
int nodesProcessed_ = 0;
int totalNodes_ = 0;
int totalCutsGenerated_ = 0;
int totalCutsApplied_ = 0;
int totalSocoSolved_ = 0; 
int totalFadingCuts_ = 0;
int totalNodeFadingCuts_ = 0;
double objectiveTolerance_ = 1e-4;
int numVars_ = 0;
int N = 0; /* number of assets, always 1 less than numVars */
Node* root;
std::vector< Node* > nodeList_;
std::vector< Node* > allNodeList_; // Keeps pointer of all nodes for destructor
int printLevel_ = 2;

int firstFeasibleObjective_ = 0;

int main(int argc, char* argv[]) {

    std::cout << "QDBB." << std::endl;
    root = NULL;
    

    numVars_ = atoi(argv[1])+1;
    N = numVars_ - 1;

    cutRule_ = atoi(argv[2]);
    cutPriority_ = atoi(argv[3]);
    cutSelection_ = atoi(argv[4]);
    cutLimit_ = atoi(argv[5]);
    cutPerIteration_ = atoi(argv[6]);
    iterationLimit_ = atoi(argv[7]);
    
    MSKenv_t env = NULL;
    MSKtask_t task = NULL;
    MSKrescodee r;

    /* Create the mosek environment. */
    r = MSK_makeenv(&env,NULL);
    r = MSK_maketask(env,3,N+1,&task);
    
    oProblem = task;
    
    startBB(argv);
    
    return 1;
}
 
int startBB(char* argv[]) {
    int status = 0;
    createProblem(&oProblem, argv);
    
    root = new Node; //(Node *) malloc(sizeof(Node));
    createNewNode(0, &root, -1 /*var id*/, 0 /* bound */, 0 /* lower */);
    
    clock_t begin = clock();
    while(true) {
        
        printTxt(3, "Active Nodes: %d", nodeList_.size());
        
        if(nodeList_.size() == 0) {
            break;
        }
        
        Node* activeNode;
        selectNode(&activeNode);
	if(!activeNode->eliminated) {
	  printTxt(4, "Processing the node now...");
	  nodesProcessed_++;
	  solveLP(activeNode);
	  totalSocoSolved_++;
	  if(!activeNode->feasible) {
	    goto endofcutting;
	  }
	  isIntFeasible(activeNode);
	  if(activeNode->intfeasible) {
	    goto endofcutting;
	  }
	  if(cutRule_>0) { // B&C
	    int totalCut = 0;
	    int iterN = 0;
	    //double prevObjective = activeNode->nodeObj;
	    if(cutRule_==1) {
	      while(iterN < iterationLimit_) {
		int iterCut = 0;
		while(iterCut < cutPerIteration_ && activeNode->totalCuts < cutLimit_) {
		  int isNewCut = cut(activeNode);
		  activeNode->totalCuts += isNewCut;
		  totalCutsApplied_ += isNewCut;
		  totalCutsGenerated_++;
		  iterCut++;
		  totalCut++;
		}
		solveLP(activeNode);
		totalSocoSolved_++;
		iterN++;
	      }
	    } else if(cutRule_ == 2) {
	      double prevObjective = activeNode->nodeObj + 2*objectiveTolerance_;
	      double newObjective = activeNode->nodeObj;
	      totalNodeFadingCuts_++;
	      while(newObjective <= prevObjective - objectiveTolerance_) {
		if(activeNode->totalCuts < cutLimit_) {
		  int isNewCut = cut(activeNode);
		  totalFadingCuts_ += isNewCut;
		  activeNode->totalCuts += isNewCut;
		  totalCutsGenerated_ ++;
		  totalCutsApplied_ += isNewCut;
		} else {
		  printTxt(6, "Cut limit reached for node %d",activeNode->ID);
		  break;
		}
		solveLP(activeNode);
		totalSocoSolved_++;
		if(!activeNode->feasible) { break; }
		prevObjective = newObjective;
                newObjective = activeNode -> nodeObj;
	      }
	    } else {
	      printTxt(3, "Undefined cutting rule. Please check readme.");
	    }
	  }

	  //cut(activeNode);
	  //solveLP(activeNode);
	} else {
	  continue;
	}

        endofcutting:

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
    clock_t end = clock();
    double elapsed = double(end-begin) / CLOCKS_PER_SEC;

    // Summary report
    printTxt(1, "========== SUMMARY ==========");
    printTxt(1, "Best objective value: %f", globalUpperBound_);
    printTxt(1, "Values: ");
    for(int i=0; i<N; i++) {
      printTxt(1, "\t%d: %f", (i+1), bestSoln_[i]);
    }
    printTxt(1, "Number of nodes processed: %d", nodesProcessed_);
    printTxt(1, "Number of nodes generated: %d", totalNodes_);
    printTxt(1, "Total SOCO solved: %d", totalSocoSolved_);
    printTxt(1, "Total time elapsed: %f seconds",  elapsed );
    printTxt(1, "Total cuts generated: %d", totalCutsGenerated_);
    printTxt(1, "Total cuts applied: %d", totalCutsApplied_);
    if(cutRule_==2) {
      printTxt(1, "Average effective cuts: %f", (double) totalFadingCuts_ / totalNodeFadingCuts_);
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
    Node* leftNode =  new Node; // (Node*) malloc(sizeof(Node));
    createNewNode(activeNode /*parent*/, &leftNode, asset /*var id*/,  floor(activeNode->nodeSoln[asset])/* bound */, 0 /* lower */);
    // Right - Greater than bound
    Node* rightNode = new Node; //(Node*) malloc(sizeof(Node));
    createNewNode(activeNode /*parent*/, &rightNode, asset /*var id*/, ceil(activeNode->nodeSoln[asset]) /* bound */, 1 /* lower */);
    
    
    // Add them to active list
    // Done in create
    
    return status;
}

int cut(Node* activeNode) {
    int status = 0;
    
    int variableForCut = nextCut(N, cutPriority_, activeNode->nodeSoln, &(activeNode->usedCuts)); 
    
    if(variableForCut >= 0 ) {
      status = addNewCut(activeNode->problem, variableForCut+1, activeNode->nodeSoln[variableForCut], cutSelection_);
    }
    
    return status;
}

int selectNode(Node** activeNode) {
    int status = 0;
    int selected = 0;
    int max_depth = 0;
    *activeNode = nodeList_[0];
    for(unsigned int i=0; i < nodeList_.size(); i++) {
        if(nodeList_[i]->depth > max_depth && i>selected) {
            selected = i;
            *activeNode = nodeList_[i];
	    max_depth = nodeList_[i]->depth;
        }
    }
    nodeList_.erase(nodeList_.begin()+selected);
    
    printTxt(3,"Node selected: %d (%p)", (*activeNode)->ID, *activeNode);
    
    return status;
}

int eliminateNodes() {
  for(unsigned int i=0; i<nodeList_.size(); i++) {
     if(nodeList_[i]->lowerBound >= globalUpperBound_) {
       nodeList_[i]->eliminated = true;
     }
  }

  return 1;
}

int createNewNode(Node* parent, Node** newNode, int varID, double bound, int lower) {
    int status = 0;
    
 
    if(!parent) {
    
        
    
        for (int i = 0; i < N; ++i) { // Initialize array of used cuts
            std::vector<int> row; // Create an empty row for each asset  
            (*newNode)->usedCuts.push_back(row); // Add the row to the main vector
        }
        
        
        (*newNode)->ID = 0;
        (*newNode)->nodeObj = 0;
        (*newNode)->lowerBound = 0;
        (*newNode)->feasible = false;
        (*newNode)->intfeasible = false;
	(*newNode)->eliminated = false;
        (*newNode)->problem = oProblem;
        (*newNode)->parent = 0;
        //(*newNode)->usedCuts;
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
	  std::vector<int> newRow(parent->usedCuts[i]); // Create an empty row for each asset
	  (*newNode)-> usedCuts.push_back(newRow); // Add the row to the main vector
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
        //(*newNode)->usedCuts.swap(usedCuts);
        (*newNode)->depth = parent->depth+1;
	(*newNode)->totalCuts = parent->totalCuts;
        
	
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

    //printTxt(6,"Task is preparing to be cloned");


    MSK_clonetask(aNode->problem, &mtask);
    
    printTxt(6,"Task is cloned");


    //MSKvariabletypee* vartype = malloc(numVars_*sizeof(MSKvariabletypee));
    for(int i=0; i<N; ++i) {
        MSK_putvartype (mtask, i+1, MSK_VAR_TYPE_CONT);
    }

    //MSK_toconic(mtask);

    //std::cout << "Running so far!" << std::endl;
    //MSK_linkfunctotaskstream (mtask, MSK_STREAM_LOG, NULL, printstr);
    //printTxt(6,"Ready to solve with MOSEK");
    MSK_optimize(mtask);
    printTxt(6,"Problem is solved with  MOSEK");
   
    MSK_solutionsummary(mtask,MSK_STREAM_LOG);
    
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
	printTxt(3,"Node status unknown, node is pruned.");
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
        std::cout  << level << ":";
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



