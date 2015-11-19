
#include "qdbb.h"

//#define OUTLEV 5

int OUTLEV = 1;
int FILEOUTPUT = 0;

#define printText(flag, ...) if (flag <= OUTLEV) { printf("%d:%*s",flag,2*flag,""); printf(__VA_ARGS__); printf("\n"); }


/*static void MSKAPI printstr(void *handle, MSKCONST char str[]) 
{ 
  printf("%s",str); 
  } */



MSKtask_t oProblem = NULL; // Original problem instance
double globalUpperBound_ = std::numeric_limits<double>::infinity();
int finiteUpperBound_ = 0;
double* bestSoln_;
int branchingRule_ = 0; // 0: most fractional, 1: highest cost, 2: random, 3:bonami 
int cutPriority_ = 0; // 0: default, most fractional, 1: best-value
int cutRule_ = 1; // 0: default, no cut, 1: always cut, 2: fading-cuts, 3: root-heuristic cut
// 4: min depth for cut, 5: only if deep cut
int searchRule_ = 0; // 0: default, depth first, left, 1: depth first right, 2: breadth first, 3: best (lower bound)
double deepCutThreshold_ = 0; // deep-cut threshold value
double bestImprovement_ = 0;
double totalImprovement_ = 0;
int controlValve = 50;
int cutSelection_ = 1;
int cutLimit_ = 100; // max number of cuts to add, default is 1
int cutPerIteration_ = 3; // number of cuts to be added at each relaxation
int iterationLimit_ = 1;
int minCutDepth_ = 0;
int maxCutDepth_ = 100;
int nodesProcessed_ = 0;
int totalNodes_ = 0;
int totalCutsGenerated_ = 0;
int totalCutsApplied_ = 0;
int totalSocoSolved_ = 0; 
int totalFadingCuts_ = 0;
int totalNodeFadingCuts_ = 0;
int bestNodeNumber_ = 0;
double objectiveTolerance_ = 1e-3;
double integerTolerance_ = 1e-5;
int numVars_ = 0;
int N = 0; /* number of assets, always 1 less than numVars */
/* TODO Change N to number of variables! */
Node* root;
Node* bestNode;
std::vector< Node* > nodeList_;
std::vector< Node* > allNodeList_; // Keeps pointer of all nodes for destructor

extern int PROBLEMCODE;
extern double* mu;
extern double** Q;

int firstFeasibleObjective_ = 0;

int main(int argc, char* argv[]) {

  std::cout << "QDBB." << std::endl;
  root = NULL;
  srand(42); // Magic number
    
  MSKenv_t env = NULL;
  MSKtask_t task = NULL;

  /* Create the mosek environment. */
  MSK_makeenv(&env,NULL);
  MSK_maketask(env,0,0,&task);
  
  oProblem = task;

  parseInfo(argc,argv);
  
  startBB(argc, argv);
  
  finishBB();
  
  //MSK_deletetask (&task); 
  MSK_deleteenv (&env);
  
  return 1;
}

int parseInfo(int argc, char* argv[]) {

  int i = 0;
  
  for(i=0; i<argc; i++) {
    char* tmp = argv[i];
    if(argv[i][0]!='-') {
      continue;
    }
    //printText(3,"Option %s: %s\n", argv[i],argv[i+1]);
    
    if(strcmp(tmp, "-v")==0 || strcmp(tmp, "-a")==0) { // Number of assets
      N = atoi(argv[i+1]);
      numVars_ = N + 1;
    }

    if(strcmp(tmp, "-o")==0)
      OUTLEV = atoi(argv[i+1]);

    if(strcmp(tmp, "-f")==0)
      FILEOUTPUT = atoi(argv[i+1]);

    if(strcmp(tmp, "-b")==0) { // Branching rule
      if(strcmp(argv[i+1],"mf")==0) // most fractional
	branchingRule_ = 0;
      else if(strcmp(argv[i+1],"hc")==0) // highest cost
	branchingRule_ = 1;
      else if(strcmp(argv[i+1],"random")==0) // random
	branchingRule_ = 2;
      else if(strcmp(argv[i+1],"bonami")==0) // bonami
	branchingRule_ = 3;
      else if(strcmp(argv[i+1],"hvar")==0) // bonami-1 static
	branchingRule_ = 4;      
    }

    if(strcmp(tmp, "-c")==0) { // Cut rule
      if(strcmp(argv[i+1],"mf")==0) // most fractional
	cutPriority_ = 0;
      else if(strcmp(argv[i+1],"hc")==0) // highest cost
	cutPriority_ = 1;
      else if(strcmp(argv[i+1],"random")==0) // random
	cutPriority_ = 2;
      else if(strcmp(argv[i+1],"bonami")==0) // bonami
	cutPriority_ = 3;
      else if(strcmp(argv[i+1],"hvar")==0) // bonami-1 static
	cutPriority_ = 4;
    }

    if(strcmp(tmp, "-s")==0) { // Search rule
      if(strcmp(argv[i+1],"df0")==0) // Depth first - left
	searchRule_ = 0;
      else if(strcmp(argv[i+1],"df1")==0) // Depth first - right
	searchRule_ = 1;
      else if(strcmp(argv[i+1], "breadth")==0) // breadth-first
	searchRule_ = 2;
      else if(strcmp(argv[i+1], "best")==0) // best lower bound
	searchRule_ = 3;
    }

    if(strcmp(tmp, "-l")==0) // cut limit
      cutLimit_ = atoi(argv[i+1]);

    if(strcmp(tmp, "-i")==0) // cut generation iteration
      iterationLimit_ = atoi(argv[i+1]);

    if(strcmp(tmp, "-p")==0) // Cut per iteration
      cutPerIteration_ = atoi(argv[i+1]);

    if(strcmp(tmp, "-x")==0) { // Termination
      if(strcmp(argv[i+1],"0")==0) // No Cut
	cutRule_ = 0;
      if(strcmp(argv[i+1],"1")==1) // Always cuts!
	cutRule_ = 1;
      else if(strcmp(argv[i+1],"2")==0) // Fading cuts
	cutRule_ = 2;
      else if(strcmp(argv[i+1],"3")==0) // Special: all cuts
	cutRule_ = 5;
    }
    
    if(strcmp(tmp, "-mind")==0) // Min depth for cut generation
      minCutDepth_ = atoi(argv[i+1]);
    if(strcmp(tmp, "-maxd")==0) // Min depth for cut generation
      maxCutDepth_ = atoi(argv[i+1]);
    
    if(strcmp(tmp, "-dct")==0) // Deep cut threshold
      deepCutThreshold_ = atof(argv[i+1]);

  }

  printText(1,"Input= a: %d, branch: %d, cut: %d, iter: %d, cutper: %d, term: %d", N, branchingRule_, cutPriority_, iterationLimit_, cutPerIteration_, cutRule_);
  
  return 1;
}

int startBB(int argc, char* argv[]) {
  int status = 0;
  createProblem(&oProblem, argc, argv);
  
  root = new Node; //(Node *) malloc(sizeof(Node));
  createNewNode(0, &root, -1 /*var id*/, 0 /* bound */, 0 /* lower */);
  
  if(OUTLEV==1) {
    printText(1, "Update      Node    Best Obj.    Current Gap   Lowest LB  Node  ");
    printText(1, "---------   -----   ----------   -----------   ---------- ------");
  }      

  clock_t begin = clock();
  while(true) {
    
    printText(3, "Active Nodes: %ld", nodeList_.size());
    
    if(nodeList_.size() == 0) {
      break;
    }
    
    Node* activeNode;
    selectNode(&activeNode);
    if(!activeNode->eliminated) {
      printText(4, "Processing the node now...");
      activeNode->processed = true;
      nodesProcessed_++;
      solveLP(activeNode);
      totalSocoSolved_++;
      if(!activeNode->feasible) {
        printText(4, "Node's (%d) relaxation is infeasible. Pruning...", activeNode->ID);
        goto finaldecision;
      }
      isIntFeasible(activeNode);
      if(activeNode->intfeasible) {
        goto finaldecision;
      }
      if(activeNode->nodeObj > globalUpperBound_ + 1e-12) {
        printText(3, "Node (%d) has a higher objective than upper bound, pruning...", activeNode->ID);
        goto finaldecision;
      }
      if(cutRule_>0 && activeNode->totalCuts < cutLimit_ && activeNode->depth >= minCutDepth_ && activeNode->depth <= maxCutDepth_) { // B&C
        int totalCut = 0;
        int iterN = 0;
        //double prevObjective = activeNode->nodeObj;
        double firstObjective = activeNode->nodeObj; // Objective at relaxation
	int firstNCut = totalCutsApplied_;
        if(cutRule_==1) {
          while(iterN < iterationLimit_) {
            printText(4, "Node (%d) cut iteration %d of %d", activeNode->ID, iterN+1, iterationLimit_);
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
	    // IMPROVEMENT
	    if(totalCutsApplied_ > firstNCut) {
	      long double objImprovement = activeNode->nodeObj -  firstObjective;
	      if(objImprovement  > bestImprovement_) {
		bestImprovement_ = objImprovement;
	      }
	      totalImprovement_ += objImprovement;
	      firstObjective = activeNode->nodeObj;
	      firstNCut = totalCutsApplied_;
	    }
	    // ENDOFIMPROVEMENT
	    totalSocoSolved_++;
            iterN++;
          }
        }
        else if(cutRule_ == 2) { // Fading Cuts, -x 1
          double prevObjective = activeNode->nodeObj + 2*objectiveTolerance_;
          double newObjective = activeNode->nodeObj;
          totalNodeFadingCuts_++;
	  int fadingIter = 0;
          while(newObjective <= prevObjective - objectiveTolerance_ && fadingIter < controlValve ) {
            if(activeNode->totalCuts < cutLimit_) {
              int isNewCut = cut(activeNode);
              totalFadingCuts_ += isNewCut;
              activeNode->totalCuts += isNewCut;
              totalCutsGenerated_ ++;
              totalCutsApplied_ += isNewCut;
            } else {
	      //printf("Tots: %d, cut lim: %d\n",activeNode->totalCuts,cutLimit_);
              printText(6, "Cut limit reached for node %d",activeNode->ID);
              break;
            }
            solveLP(activeNode);
            totalSocoSolved_++;
	    // IMPROVEMENT
	    if(totalCutsApplied_ > firstNCut) {
	      long double objImprovement = activeNode->nodeObj -  firstObjective;
	      if(objImprovement  > bestImprovement_) {
		bestImprovement_ = objImprovement;
	      }
	      totalImprovement_ += objImprovement;
	      firstObjective = activeNode->nodeObj;
	      firstNCut = totalCutsApplied_;
	    }
	    // ENDOFIMPROVEMENT
	    if(!activeNode->feasible) { break; }
            prevObjective = newObjective;
            newObjective = activeNode -> nodeObj;
	    fadingIter++;
          }
        } 
        else if(cutRule_ == 5) { // all cuts at root node!   -x 2
          
          if(activeNode->ID==0) { // If it is root node
          
            for(int i=0; i<N; i++) {
              int isNewCut = cut(activeNode);
              activeNode->totalCuts += isNewCut;
              totalCutsGenerated_ ++;
              totalCutsApplied_ += isNewCut;
            }
            solveLP(activeNode);
	    // IMPROVEMENT
	    if(totalCutsApplied_ > firstNCut) {
	      long double objImprovement = activeNode->nodeObj -  firstObjective;
	      if(objImprovement  > bestImprovement_) {
		bestImprovement_ = objImprovement;
	      }
	      totalImprovement_ += objImprovement;
	      firstObjective = activeNode->nodeObj;
	      firstNCut = totalCutsApplied_;
	    }
	    // ENDOFIMPROVEMENT
	    totalSocoSolved_++;
            
          }
          
	  // no cut for other nodes
          
        }
	else if(cutRule_ == 6) { // Order cuts on violation choose fixed num
	  // TODO
	}
	else {
          printText(3, "Undefined cutting rule. Please check readme.");
        }
        
      }

      //cut(activeNode);
      //solveLP(activeNode);
    } else {
      printText(3, "Node %d is already eliminated, Current LB: %e, Global UB: %e", activeNode->ID, activeNode->lowerBound, globalUpperBound_);
      continue;
    }

    
    if( activeNode->feasible) { isIntFeasible(activeNode); }

finaldecision:

    if( activeNode->feasible && activeNode->intfeasible) { // better objective function
      if(activeNode->nodeObj < globalUpperBound_) {
        globalUpperBound_ = activeNode->nodeObj;
	finiteUpperBound_ = 1;
        bestSoln_ = activeNode->nodeSoln;
	bestNode = activeNode;
	double curGap = 0;
	double lowBnd = 0;
	int lowNode = 0;
	getCurrentGap(&curGap, &lowBnd, &lowNode);
	if(lowBnd > globalUpperBound_) {
	  curGap = 0;
	  lowBnd = globalUpperBound_;
	  lowNode = activeNode->ID;
	}
        printText(1, "New bound   %5d   %10.3e     %5.2f %%     %10.3e (%d)", activeNode->ID, globalUpperBound_, curGap, lowBnd, lowNode);
        bestNodeNumber_ = activeNode->ID;
        eliminateNodes();
	if(FILEOUTPUT)
	  MSK_writedata(activeNode->problem, "result/bestnode.mps");
      }
    } else if(activeNode->feasible && !activeNode->intfeasible) {
      // branch
      if(activeNode->nodeObj <= globalUpperBound_) {
        printText(3, "Upper bound: %e, current objective: %e, branching...",globalUpperBound_,activeNode->nodeObj);
        branch(activeNode);
      }
    } else {
      continue;
    }
    
  }
  clock_t end = clock();
  double elapsed = double(end-begin) / CLOCKS_PER_SEC;

  // Summary report
  printText(1, "========== SUMMARY ==========");
  if(globalUpperBound_== std::numeric_limits<double>::infinity()) {
    printText(1, "Problem is infeasible :(");
    return 0;
  }
  
  
  printText(1, "Best objective value: %e", globalUpperBound_);
  printText(2, "Values: ");
  for(int i=0; i<N; i++) {
    if(bestSoln_[i] < integerTolerance_) {
      printText(3, "  %2d: %.12f\t%.8f", (i+1), bestSoln_[i], bestNode->mosekSoln[i+1]);
    } else {
      printText(2, "  %2d: %.12f\t%.8f", (i+1), bestSoln_[i], bestNode->mosekSoln[i+1]);
    }
  }
  for(int i=0; i<2*N+1; i++) {
    //printf("Var[%d]:\t%.6f\n", i, bestNode->mosekSoln[i+1]);
  }
  printText(1, "Number of nodes processed: %d", nodesProcessed_);
  printText(1, "Number of nodes generated: %d", totalNodes_);
  printText(1, "Total SOCO solved: %d", totalSocoSolved_);
  printText(1, "Total time elapsed: %f seconds",  elapsed );
  printText(1, "Total cuts generated: %d", totalCutsGenerated_);
  printText(1, "Total cuts applied: %d", totalCutsApplied_);
  if(cutRule_==2)
    printText(1, "Average effective cuts: %f", (double) totalFadingCuts_ / totalNodeFadingCuts_);
  printText(1, "Total objective improvement: %e", totalImprovement_);
  printText(1, "Best objective improvement: %e", bestImprovement_);  
  printText(1,"Optimal node: %d", bestNodeNumber_);
  printText(1,"Optimal node depth: %d", bestNode->depth);
  printText(1,"Optimal value: %e", globalUpperBound_);
  printText(1,"Done...");
  return status;
}

int branch(Node* activeNode) {
  int status = 0;
  int asset = 0;

  if(branchingRule_ == 0) { // most fractional branching
    double initval = activeNode->nodeSoln[0];
    double mostfrac = fmin(initval - floor(initval), ceil(initval) - initval);
    for(int i=1; i<N; ++i) {
      double cValue = activeNode->nodeSoln[i];
      double smallest = fmin(cValue - floor(cValue), ceil(cValue) - cValue);
      if(smallest > mostfrac) {
	asset = i;
	mostfrac = smallest;
      }
    }
  }
  else if (branchingRule_ == 1) { // highest cost branching
    double *cost = mu;
    asset = -1;
    double hc = -1000000; //cost[0];
    for(int i=0; i<N; ++i) {
      //MSKvariabletypee vartype;
      //MSK_getvartype(activeNode->problem, i+1, &vartype);
      //if(vartype == MSK_VAR_TYPE_INT) {
      //printf("Variable %d, value: %f, round: %f, diff: %f, hc: %f\n", i, activeNode->nodeSoln[i], round(activeNode->nodeSoln[i]), fabs(activeNode->nodeSoln[i+1] - round(activeNode->nodeSoln[i+1])), hc);
	if( fabs(activeNode->nodeSoln[i] - round(activeNode->nodeSoln[i])) > integerTolerance_) {
	  if( cost[i] > hc ) {
	    asset = i;
	    hc = cost[i];
	  }
	  //}
    
	}
    }
  }
  else if (branchingRule_ == 2) { // random
    int i=0;
    while(true) {
      i = rand() % N;
      if( fabs(activeNode->nodeSoln[i] - round(activeNode->nodeSoln[i])) > integerTolerance_) {
	asset = i;
	break;
      }
    }
    
    
  } 
  else if (branchingRule_ == 3) { // bonami dynamic
    double hdiff = -std::numeric_limits<double>::infinity();
    asset = 0;
    double* soln = activeNode->nodeSoln;
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
      double score = 1*lower + 2*upper;
      //printf("[%d] score: %f\n",i,score);
      if( fabs(activeNode->nodeSoln[i] - round(activeNode->nodeSoln[i])) > integerTolerance_) {
	if(score > hdiff) {
	  hdiff = score;
	  asset = i;
	}
      }
    }
    
    // not implemented yet!
    
  }
  else if (branchingRule_ == 4) { // bonami static - hvar
    asset = 0;
    double hvar = -std::numeric_limits<double>::infinity();
    for(int i=0; i<N; ++i) {
      if( fabs(activeNode->nodeSoln[i] - round(activeNode->nodeSoln[i])) > integerTolerance_) {
	if(Q[i][i] > hvar) {
	  asset = i;
	  hvar = Q[i][i];
	}
      }
    }
  }
  
  printText(2,"Branching at [Node %04d], on var(%02d), <=%2.0f [Node %04d] and >=%2.0f [Node %04d], Obj: %.5f",activeNode->ID, asset, floor(activeNode->nodeSoln[asset]), totalNodes_, ceil(activeNode->nodeSoln[asset]), totalNodes_+1, activeNode->nodeObj);
  
  // Call newNode twice

  // Left - Less than bound
  Node* leftNode =  new Node; // (Node*) malloc(sizeof(Node));
  createNewNode(activeNode /*parent*/, &leftNode, asset /*var id*/,  floor(activeNode->nodeSoln[asset])/* bound */, 0 /* lower */);
  printText(3,"New node (%d) is child of (%d) with var(%d)<=%.2f",totalNodes_-1, activeNode->ID, asset, floor(activeNode->nodeSoln[asset])); 
  // Right - Greater than bound
  Node* rightNode = new Node; //(Node*) malloc(sizeof(Node));
  createNewNode(activeNode /*parent*/, &rightNode, asset /*var id*/, ceil(activeNode->nodeSoln[asset]) /* bound */, 1 /* lower */);
  printText(3,"New node (%d) is child of (%d) with var(%d)>=%.2f",totalNodes_-1, activeNode->ID, asset, ceil(activeNode->nodeSoln[asset]));
  
  // Add them to active list
  // Done in create

  
  printToFile(leftNode);
  printToFile(rightNode);
  
  
  return status;
}

int cut(Node* activeNode) {
  int status = 0;
  
  int variableForCut = nextCut(N, cutPriority_, activeNode->nodeSoln, &(activeNode->usedCuts)); 
  if(variableForCut >= 0 ) {

    double currsoln = activeNode->nodeSoln[variableForCut];
    if(PROBLEMCODE==1) { // TODO Remove this if condition in future - problem specific
      if(fabs(currsoln-round(currsoln))<1e-6) { // TODO make this one a parameter
      printText(3,"Node (%d) variable %d is close to integer %.3f, cut is not generated",activeNode->ID,variableForCut,currsoln);
      return 0;
    } 
    }



    status = addNewCut(activeNode->problem, variableForCut+1, activeNode->nodeSoln[variableForCut], cutSelection_);
    if(status>0) {
      printText(3, "Node (%d) New cut for variable %d is added, value: %f", activeNode->ID, variableForCut, activeNode->nodeSoln[variableForCut]);
    } else {
      printText(3, "Node (%d) Cut generation failed, variable %d, value: %f", activeNode->ID, variableForCut, activeNode->nodeSoln[variableForCut]);
    }
  }
  
  return status;
}

int selectNode(Node** activeNode) {
  int status = 0;
  unsigned int selected = 0;
  int max_depth = -1;

  double lowestlb = nodeList_[0]->lowerBound;
    for(unsigned int i=1; i < nodeList_.size(); i++) {
      if(nodeList_[i]->lowerBound < lowestlb) {
	lowestlb = nodeList_[i]->lowerBound;
      }
    }
    double redline = 10*fabs(lowestlb);
    //printf("Lowest LB: %.5f, Redline: %.5f\n",lowestlb, redline);

  if(searchRule_ == 0) {  // df0 depth-first left
    *activeNode = nodeList_[0];
    max_depth = nodeList_[0]->depth;
    for(unsigned int i=0; i < nodeList_.size(); i++) {
      if(nodeList_[i]->depth > max_depth && nodeList_[i]->lowerBound < redline) {
	selected = i;
	*activeNode = nodeList_[i];
	max_depth = nodeList_[i]->depth;
      }
    }
  }
  else if(searchRule_ == 1) { // df1 depth-first right
    *activeNode = nodeList_[0];
    max_depth = nodeList_[0]->depth;
    for(unsigned int i=1; i < nodeList_.size(); i++) {
      if(nodeList_[i]->depth >= max_depth && nodeList_[i]->lowerBound < redline) {
	selected = i;
	*activeNode = nodeList_[i];
	max_depth = nodeList_[i]->depth;
      }
    }
  }
  else if(searchRule_ == 2) { // breadth
    *activeNode = nodeList_[0];
    selected = 0;
  }
  else if(searchRule_ == 3) { // best lower bound
    *activeNode = nodeList_[0];
    double bestbound = nodeList_[0]->lowerBound;
    for(unsigned int i=1; i < nodeList_.size(); i++) {
      if(nodeList_[i]->lowerBound > bestbound && nodeList_[i]->lowerBound < redline) {
	selected = i;
	*activeNode = nodeList_[i];
	bestbound = nodeList_[i]->lowerBound;
      }
    }
  }
  
  nodeList_.erase(nodeList_.begin()+selected);
  
  printText(3,"Node selected: %d (%p)", (*activeNode)->ID, *activeNode);
  
  return status;
}

int eliminateNodes() {
  int totaln = 0;
  for(unsigned int i=0; i<nodeList_.size(); i++) {
    if(nodeList_[i]->lowerBound >= globalUpperBound_) {
      nodeList_[i]->eliminated = true;
      totaln++;
    }
  }
  printText(3, "%d nodes are eliminated due to new upper bound", totaln);

  return 1;
}

int createNewNode(Node* parent, Node** newNode, int varID, double bound, int lower) {
  int status = 0;
  
  // No parent
  if(!parent) {
    
    
    
    for (int i = 0; i < N; ++i) { // Initialize array of used cuts
      std::vector<int> row; // Create an empty row for each asset  
      (*newNode)->usedCuts.push_back(row); // Add the row to the main vector
      std::vector<int> row2; // Empty row for branch
      (*newNode)->usedBranches.push_back(row2); // TODO (I'm not using this right now)
    }
    
    
    (*newNode)->ID = 0;
    (*newNode)->nodeObj = 0;
    (*newNode)->lowerBound = 0;
    (*newNode)->feasible = false;
    (*newNode)->intfeasible = false;
    (*newNode)->eliminated = false;
    (*newNode)->problem = oProblem;
    (*newNode)->parent = 0;
    (*newNode)->totalCuts = 0;
    //(*newNode)->usedCuts;
    (*newNode)->depth = 0;
    (*newNode) -> nodeSoln = (double*) malloc(N*sizeof(double));
    (*newNode) -> mosekSoln = (double*) malloc((2*N+1)*sizeof(double));
    
    printText(1, "Root is created, Branch and Conic Cut has started.");

    printToFile(*newNode);
    
  } else {
    
    //printText(3, "New node is created");
    
    // Add branch constraint here
    // Change problem
    MSKtask_t newProblem;
    MSK_clonetask((parent->problem), &newProblem);
    
    //MSK_analyzeproblem(parent->problem, MSK_STREAM_LOG);
    
    //MSK_analyzeproblem(newProblem, MSK_STREAM_LOG);
    int corrector = 0;
    if(PROBLEMCODE==2 || PROBLEMCODE==3 || PROBLEMCODE==4) {
      corrector = N;
      //if(bound!=0)
      //bound = bound + (1-2*lower)*(1e-6);
    }
    
    MSKboundkeye exbk;
    MSKrealt exlbound;
    MSKrealt exubound;
    //bound = bound; // + 1e-8
    // printf("Bound: %.8f\n",bound);
    if(varID!=-1) {
      MSK_getvarbound ( 
				   newProblem, 
				   varID+corrector+1, 
				   &exbk, 
				   &exlbound, 
				   &exubound);
      //printf("Variable %d, ex lb: %f, ub: %f\n",varID, exlbound, exubound);
      MSK_chgbound ( 
      newProblem,         //MSKtask_t    task, 
      MSK_ACC_VAR,        //MSKaccmodee  accmode, 
      varID+corrector+1,            //MSKint32t    i, 
      lower,              //MSKint32t    lower, 
      1,                  //MSKint32t    finite, 
      bound               //MSKrealt     value); 
      );
      
      //MSK_putbound(newProblem, MSK_ACC_VAR, varID+corrector+1, bk, lbound, ubound); -- DO NOT USE THIS ONE
      MSK_getvarbound ( 
				   newProblem, 
				   varID+corrector+1, 
				   &exbk, 
				   &exlbound, 
				   &exubound);
      //printf("Variable %d, ex lb: %f, ub: %f\n",varID, exlbound, exubound);
      
    }
    else {
      printf("Branching error: -1\n");
    }
    
    std::vector< std::vector <int> > usedCuts;
    // copy usedCuts
    for (int i = 0; i < N; i++) { // Initialize array of used cuts
      std::vector<int> newRow(parent->usedCuts[i]); // Create an empty row for each asset
      (*newNode)-> usedCuts.push_back(newRow); // Add the row to the main vector
    }
    
    std::vector< std::vector <int> > usedBranches;
    for(int i=0; i<N; i++) {
      std::vector<int> newRow(parent->usedBranches[i]); // empty row
      (*newNode)-> usedBranches.push_back(newRow);
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
    (*newNode) -> nodeSoln = (double*) malloc(N*sizeof(double));
    (*newNode) -> mosekSoln = (double*) malloc((2*N+1)*sizeof(double));
    
  }
  
  
  nodeList_.push_back(*newNode);
  allNodeList_.push_back(*newNode);
  totalNodes_++;
  
  
  return status;
}

int getCurrentGap(double* curGap, double* lowerBound, int* lowerNode) {
  // Loop thru node list
  // if not eliminated, not processed and has not a child
  // if lower bound is less than current lb, update it
  
  int found = 0;
  double lb = std::numeric_limits<double>::infinity();
  for(unsigned int i=0; i < nodeList_.size(); i++) {
    //printText(1,"%d\t%e",nodeList_[i]->ID, nodeList_[i]->lowerBound);
    if(!nodeList_[i]->eliminated && !nodeList_[i]->processed) {
      if(nodeList_[i]->lowerBound < lb) {
	lb = nodeList_[i]->lowerBound;
	*lowerNode = nodeList_[i]->ID;
	found = 1;
      }
    }
  }

  if(found) {
    *lowerBound = lb;
    *curGap = (globalUpperBound_ - lb) / (globalUpperBound_) * 100;
  }
  else {
    *lowerBound = globalUpperBound_;
    *curGap = 0;
    *lowerNode = 0;
  }

  return 1;

}

int solveLP(Node* aNode) {
  int status = 0;
  
  MSKtask_t originaltask = aNode->problem;

  MSKtask_t mtask;
  //MSK_clonetask(originaltask, &mtask);
  mtask = originaltask;
  
  MSK_putintparam(mtask, MSK_IPAR_MIO_MODE, MSK_MIO_MODE_IGNORED); // make it LP
  MSK_putintparam(mtask, MSK_IPAR_PRESOLVE_LINDEP_USE, MSK_OFF);
  MSK_putintparam(mtask, MSK_IPAR_PRESOLVE_USE, MSK_PRESOLVE_MODE_OFF);
  MSK_putintparam(mtask, MSK_IPAR_INTPNT_REGULARIZATION_USE, MSK_OFF);
  MSK_putintparam(mtask, MSK_IPAR_INTPNT_SCALING, MSK_SCALING_NONE);
  MSK_putintparam(mtask, MSK_IPAR_NUM_THREADS, 1);
  MSK_putintparam(mtask, MSK_IPAR_INTPNT_BASIS, 0);
  //MSK_putdouparam(mtask, MSK_DPAR_INTPNT_CO_TOL_INFEAS, 1e-16);
  //MSK_putdouparam(mtask, MSK_DPAR_INTPNT_CO_TOL_PFEAS, 1e-16);
  MSK_putdouparam(mtask, MSK_DPAR_INTPNT_TOL_PFEAS, 1e-16);
  //MSK_putdouparam(mtask, MSK_DPAR_INTPNT_CO_TOL_DFEAS, 1e-16);
  //MSK_putdouparam(mtask, MSK_DPAR_INTPNT_CO_TOL_REL_GAP, 1e-16);
  MSK_putintparam(mtask, MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_CONIC);
  MSK_putintparam(mtask, MSK_IPAR_INTPNT_MAX_ITERATIONS, 30000);
  MSK_putintparam(mtask, MSK_IPAR_BI_IGNORE_MAX_ITER, MSK_ON);
  //MSK_linkfunctotaskstream (mtask, MSK_STREAM_LOG, NULL, printstr);
  //MSK_putdouparam(mtask, MSK_DPAR_CHECK_CONVEXITY_REL_TOL, 100000);
  
  // QP PARAMS!
  MSK_putdouparam(mtask, MSK_DPAR_INTPNT_TOL_DFEAS, 1e-16);
  MSK_putdouparam(mtask, MSK_DPAR_INTPNT_NL_TOL_DFEAS, 1e-16);
  MSK_putdouparam(mtask, MSK_DPAR_INTPNT_NL_TOL_PFEAS, 1e-16);


  int reattempted = 1;
  //solutionpoint:

  MSKrescodee trmcode;
  MSKrescodee r; 
  r = MSK_optimizetrm(mtask,&trmcode);
  printText(6,"Problem is solved with  MOSEK");
  r = r;
  MSK_solutionsummary(mtask,MSK_STREAM_LOG);
  
  MSKsolstae solsta;
  MSKrescodee obsolsta = MSK_getsolsta (mtask, MSK_SOL_ITR, &solsta); 
  MSKrealt primalobj;

  //printf("TRMCODE: %d\n", trmcode);
  //printf("SOLSTA: %d\n", solsta);
  if(trmcode==10006) {
    solsta = MSK_SOL_STA_OPTIMAL;
    //printf("Stall!\n");
  }
  
  if (obsolsta == MSK_RES_OK) {
  switch( solsta ) 
  { 
  case MSK_SOL_STA_OPTIMAL: 
  case MSK_SOL_STA_NEAR_OPTIMAL: 
  case MSK_SOL_STA_PRIM_FEAS:
  case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
  case MSK_SOL_STA_INTEGER_OPTIMAL:
  case MSK_SOL_STA_NEAR_INTEGER_OPTIMAL:
     
      aNode->feasible = true;
      MSK_getprimalobj(mtask, MSK_SOL_ITR, &primalobj);
      aNode->nodeObj = primalobj;
      if(primalobj > aNode->lowerBound) {
        aNode -> lowerBound = primalobj;
      }
      //aNode -> nodeSoln = (double*) malloc(N*sizeof(double));
      if(PROBLEMCODE==1) {
        MSK_getsolutionslice(mtask, MSK_SOL_ITR, MSK_SOL_ITEM_XX, 1, N+1, aNode->nodeSoln);
      } else if(PROBLEMCODE==2) {
        MSK_getsolutionslice(mtask, MSK_SOL_ITR, MSK_SOL_ITEM_XX, N+1, 2*N+1, aNode->nodeSoln);
      } else if(PROBLEMCODE==3 || PROBLEMCODE==4) {
	MSK_getsolutionslice(mtask, MSK_SOL_ITR, MSK_SOL_ITEM_XX, N+1, 2*N+1, aNode->nodeSoln);
      }
      MSK_getsolutionslice(mtask, MSK_SOL_ITR, MSK_SOL_ITEM_XX, 0, 2*N+1, aNode->mosekSoln);
      
      // std::cout << primalobj << std::endl;
      printText(3, "Node (%d) Optimal Objective: %f", aNode->ID, aNode->nodeObj);
      break; 
    
  case MSK_SOL_STA_DUAL_INFEAS_CER: 
    printText(3,"Node (%d) dual infeasibility certificate found.", aNode->ID); 
  case MSK_SOL_STA_PRIM_INFEAS_CER: 
    //case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER: 
    //case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
    aNode->feasible = false;
    printText(3,"Node (%d) primal infeasibility certificate found.", aNode->ID); 
    char fname[80];
    sprintf(fname, "../test/result/infnode%d.mps", aNode->ID);
    
    if(FILEOUTPUT)
      MSK_writedata(mtask, fname);
    break; 
  case MSK_SOL_STA_UNKNOWN: 
  default:  // most probably solution is basic
    aNode->feasible = false;
    if(reattempted) {
      printText(3,"Node status unknown %d, term code: %d, pruning, attempt failed!",solsta,trmcode);
    }
    else {
      printText(4,"Node status unknown %d, term code: %d, reattempting with higher precision", solsta, trmcode);
      MSK_putdouparam(mtask, MSK_DPAR_INTPNT_CO_TOL_INFEAS, 1e-16);
      MSK_putdouparam(mtask, MSK_DPAR_INTPNT_CO_TOL_PFEAS, 1e-16);
      MSK_putdouparam(mtask, MSK_DPAR_INTPNT_CO_TOL_DFEAS, 1e-16);
      MSK_putdouparam(mtask, MSK_DPAR_INTPNT_CO_TOL_REL_GAP, 1e-16);
  
      reattempted = 1;
      //goto solutionpoint;
    }
    
    char funame[80];
    sprintf(funame, "../test/result/unknown%d.mps", aNode->ID);
    if(FILEOUTPUT)
      MSK_writedata(mtask, funame);
    break;
    
    

  }
  } else {
    printText(3,"Solution status cannot be obtained, not OK");
  } 

  //
  
  //MSK_solutionsummary (mtask, MSK_STREAM_LOG);
  
  MSK_putintparam(mtask, MSK_IPAR_MIO_MODE, MSK_MIO_MODE_SATISFIED);
  
  return status;
  
}

int isIntFeasible(Node* aNode) {
  bool intfeasible = true;
  for(int i=0; i<N; ++i) {
    
    
    
    //if( std::modf(aNode->nodeSoln[i], &intpart) >= 1e-5) {
    if(fabs(round(aNode->nodeSoln[i])-aNode->nodeSoln[i]) > integerTolerance_) {
      intfeasible = false;
      printText(5,"Integer infeasibility at variable %d, value: %f",i, aNode->nodeSoln[i]);
    }
    printText(6,"Proximity to integer (%d): %.6f, value: %.6f",i+1,fabs(round(aNode->nodeSoln[i])-aNode->nodeSoln[i]), aNode->nodeSoln[i]);
  }
  
  aNode->intfeasible = intfeasible;
  if(intfeasible) {
    printText(4,"Node %d is integer feasible with objective: %f",aNode->ID, aNode->nodeObj);
  } else {
    printText(4,"Node %d is NOT integer feasible",aNode->ID);
  }
  //intfeasible ?  : 
  return 1;
}

//int printText(int level, const char* fmt, ...) {
//
//  if(level<=printLevel_) {
//    std::cout  << level << ":";
//    for(int i=0; i<level; ++i) {
//      std::cout << "   ";
//    }
//    va_list args;
//    va_start(args,fmt);
//    vprintf(fmt,args);
//    va_end(args);
//    std::cout << std::endl;
//  }
//  return 1;
//  
//}

int finishBB() {
  
  for(unsigned int i=0; i < allNodeList_.size(); i++) {
    deleteNode(allNodeList_[i]);
  }
  
  deleteProblem();
  
  return 1;
}

int deleteNode(Node* aNode) {
  
  free(aNode -> nodeSoln);
  free(aNode -> mosekSoln);
  MSK_deletetask(&(aNode->problem)); 
  
  delete aNode;
  
  
  
  return 1;
}

int printToFile(Node* aNode) {
  
  // Get problem info
  // Print it to file

  int fno = aNode->ID;
  char fname[80];
  sprintf(fname, "../test/result/node%d.mps", fno);
  //MSK_toconic(&(aNode->problem));

  if(FILEOUTPUT)  
    MSK_writedata(aNode->problem, fname);

  return 1;
  
}








