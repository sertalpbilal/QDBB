
#include "qdbb.h"

//#define OUTLEV 5

int OUTLEV = 1;
int FILEOUTPUT = 0;

#define printText(flag, ...) if (flag <= OUTLEV) { printf("%d:%*s",flag,2*flag,""); printf(__VA_ARGS__); printf("\n"); fflush(stdout); }

//static void MSKAPI printstr(void *handle,
//			    MSKCONST char str[])
//{
//  printf("%s",str);
//}


MSKtask_t oProblem = NULL; // Original problem instance
double globalUpperBound_ = std::numeric_limits<double>::infinity();
int finiteUpperBound_ = 0;
double* bestSoln_;
int branchingRule_ = 0; // 0: most fractional, 1: highest cost, 2: random, 3:bonami, 4: hvar 
int cutPriority_ = 1; // 0: default, most fractional (MF), 1: highest cost (HC), 2: random, 3: bonami(PR), 4: highest variance (IR)
int cutRule_ = 2; // 0: no cut, 1: always cut, 2: default, fading-cuts, 3: root-heuristic cut
// 4: min depth for cut, 5: only if deep cut
int searchRule_ = 1; // 0: depth first, left, 1: default, depth first right, 2: breadth first, 3: best (lower bound)
double deepCutThreshold_ = 0; // deep-cut threshold value (percentage)
double bestImprovement_ = 0;
double totalImprovement_ = 0;
int controlValve = 20;
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
int dryRun_ = 0;
double objectiveTolerance_ = 1e-4; //1e-3;
double integerTolerance_ = 1e-5;//1e-5;
double timeLimit_ = 14400;
int terminationReason_ = 0;
int numVars_ = 0;
int N = 0; /* number of assets, always 1 less than numVars */
/* TODO Change N to number of variables! */
Node* root;
Node* bestNode;
std::vector< Node* > nodeList_;
std::vector< Node* > allNodeList_; // Keeps pointer of all nodes for destructor
std::vector<long double*> cutImprovements_;
long double cdepth_; // depth of latest cyt
int cvariable_; // variable that cut is generated on
double cvarvalue_; // value of the variable at the time cut is generated
double branchingVarValue_ = 0;

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
      else if(strcmp(argv[i+1], "dfc")==0) // Depth first - closest
	searchRule_ = 4;
    }

    if(strcmp(tmp, "-l")==0) // cut limit
      cutLimit_ = atoi(argv[i+1]);

    if(strcmp(tmp, "-i")==0) // cut generation iteration
      iterationLimit_ = atoi(argv[i+1]);

    if(strcmp(tmp, "-p")==0) // Cut per iteration
      cutPerIteration_ = atoi(argv[i+1]);

    if(strcmp(tmp, "-x")==0) { // Termination
      if(strcmp(argv[i+1],"0")==0) // No Cut (BB)
	{cutRule_ = 0;}
      else if(strcmp(argv[i+1],"1")==0) // Always cuts! (BCC-F)
	{cutRule_ = 1;}
      else if(strcmp(argv[i+1],"2")==0) // Fading cuts (BCC-I)
	{cutRule_ = 2;}
      else if(strcmp(argv[i+1],"3")==0) // Special: all cuts (BCC-R)
	{cutRule_ = 5;}
      else if(strcmp(argv[i+1],"4")==0) // Max violation order (BCC-D)
	{cutRule_ = 6;}
      else if(strcmp(argv[i+1],"10")==0) // Pure MOSEK (MOSEK)
	{cutRule_ = 10;}
      else if(strcmp(argv[i+1],"11")==0) // Mosek with cuts at root (MOSEK-R)
	{cutRule_ = 11;}
    }
    
    if(strcmp(tmp, "-mind")==0) // Min depth for cut generation
      minCutDepth_ = atoi(argv[i+1]);
    if(strcmp(tmp, "-maxd")==0) // Min depth for cut generation
      maxCutDepth_ = atoi(argv[i+1]);
    
    if(strcmp(tmp, "-dct")==0) // Deep cut threshold
      deepCutThreshold_ = atof(argv[i+1]);

    if(strcmp(tmp, "-timelimit")==0 || strcmp(tmp, "-tl")==0) // Time limit
      timeLimit_ = atof(argv[i+1]);

    if(strcmp(tmp, "-dryrun") == 0 || strcmp(tmp, "-dr")==0)
      dryRun_ = atoi(argv[i+1]);
    
  }

  printText(1,"Input= a: %d, branch: %d, cut: %d, iter: %d, cutper: %d, term: %d", N, branchingRule_, cutPriority_, iterationLimit_, cutPerIteration_, cutRule_);
  
  return 1;
}

int startBB(int argc, char* argv[]) {


  
  int status = 0;
  createProblem(&oProblem, argc, argv);

  if(PROBLEMCODE==2) {
    integerTolerance_ = 1e-4;
  }

  
  root = new Node; //(Node *) malloc(sizeof(Node));
  createNewNode(0, &root, -1 /*var id*/, 0 /* bound */, 0 /* lower */);
  
  if(OUTLEV==1) {
    printText(1, "Update      Node    Best Obj.    Current Gap   Lowest LB   Node   Cuts Appl / Gen");
    printText(1, "---------   -----   ----------   -----------   ---------- ------  ---------------");
  }      

  clock_t begin = clock();
  while(true) {
    
    printText(3, "Active Nodes: %ld", nodeList_.size());
    
    if(nodeList_.size() == 0) {
      break;
    }
    
    Node* activeNode;
    clock_t current = clock();
    double elapsed = double(current-begin) / CLOCKS_PER_SEC;
    if(elapsed > timeLimit_) { terminationReason_ = 2; goto summary_report; }
    selectNode(&activeNode);
    if (dryRun_==1) {
      printText(1, "Dry run, printing root node to file...");
      printToFile(activeNode);
      goto summary_report;
    }
    if(!activeNode->eliminated) {
      printText(4, "Processing the node now...");
      activeNode->processed = true;
      nodesProcessed_++;
      solveLP(activeNode,1);
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
      //std::cout << minCutDepth_ << " : " << activeNode->depth << std::endl;
      if( (cutRule_>0 ) && (activeNode->totalCuts < cutLimit_) && (activeNode->depth >= minCutDepth_) && (activeNode->depth <= maxCutDepth_) ) { // B&C
        int totalCut = 0;
        int iterN = 0;
        //double prevObjective = activeNode->nodeObj;
        double firstObjective = activeNode->nodeObj; // Objective at relaxation
	int firstNCut = totalCutsApplied_;
        if(cutRule_==1) { // Always cut if possible
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
            solveLP(activeNode,1);
	    // IMPROVEMENT
	    long double objImprovement = activeNode->nodeObj -  firstObjective;
	    if(totalCutsApplied_ > firstNCut) {
	      if(objImprovement  > bestImprovement_) {
		bestImprovement_ = objImprovement;
	      }
	      if(totalCutsApplied_ == firstNCut+1) { // Exactly 1 cut
		long double* impr = (long double*) malloc(6*sizeof(long double));
		impr[0] = cdepth_; // depth
		impr[1] = activeNode->depth; // tree depth
		impr[2] = objImprovement; // obj improvement
		impr[3] = activeNode->nodeObj; // current objective
		impr[4] = cvariable_; // asset
		impr[5] = cvarvalue_; //value
		cutImprovements_.push_back(impr);
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
        else if(cutRule_ == 2) { // Fading Cuts, -x 2
          double prevObjective = 0; //activeNode->nodeObj - 2*objectiveTolerance_;
          double newObjective = 1; //activeNode->nodeObj;
          totalNodeFadingCuts_++;
	  int fadingIter = 0;
          while(newObjective >= prevObjective + objectiveTolerance_ && fadingIter < controlValve ) {
	    prevObjective = activeNode->nodeObj;
            if(activeNode->totalCuts < cutLimit_) {
              int isNewCut = cut(activeNode);
	      totalFadingCuts_ += isNewCut;
              activeNode->totalCuts += isNewCut;
              totalCutsGenerated_ ++;
              totalCutsApplied_ += isNewCut;
            } else {
	      //printf("Tots: %d, cut lim: %d\n",activeNode->totalCuts,cutLimit_);
              printText(5, "Cut limit reached for node %d",activeNode->ID);
              break;
            }
	    printText(4, "Mid-loop");
            solveLP(activeNode,1);
            newObjective = activeNode -> nodeObj;
	    totalSocoSolved_++;
	    if(!activeNode->feasible) { break; }
            fadingIter++;
	    printText(6, "Diff in objectives: %.5e", newObjective-prevObjective);
	    printText(6, "Iteration: %d / %d", fadingIter, controlValve);
	    double objChange = newObjective-prevObjective;
	    totalImprovement_ += objChange;
	    if (objChange > bestImprovement_) {
	      bestImprovement_ = objChange;
	    }
	    
          }
        } 
        else if(cutRule_ == 5) { // all cuts at root node!   -x 3
          
          if(activeNode->ID==0) { // If it is root node
          
            for(int i=0; i<N; i++) {
              int isNewCut = cut(activeNode);
              activeNode->totalCuts += isNewCut;
              totalCutsGenerated_ ++;
              totalCutsApplied_ += isNewCut;
            }
            solveLP(activeNode,1);
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
	  // Step 0: For each asset, generate cut, and then add top -p of them
	  // Cutting order is not important for this one (-c)
	  iterN = 0;
	  while(iterN < iterationLimit_) {
            printText(4, "Node (%d) cut iteration %d of %d", activeNode->ID, iterN+1, iterationLimit_);
	    int isNewCut = 0;
	    if(activeNode->totalCuts < cutLimit_) {
	      isNewCut = cut(activeNode); // could be more than 1
              activeNode->totalCuts += isNewCut;
              totalCutsApplied_ += isNewCut;
              totalCutsGenerated_ += N;
              totalCut += N;
	    }
	    else {
	      printText(4, "Cut limit is reached for node, terminating...");
	      break;
	    }
	    if (isNewCut == 0) {
	      printText(4, "No improvement after cut, terminating...");
	      break;
	    }
	    double prevObj = activeNode->nodeObj;
            solveLP(activeNode,1);
	    double newObj = activeNode->nodeObj;
	    double objChange = newObj - prevObj;
	    totalImprovement_ += objChange;
	    if (objChange > bestImprovement_) {
	      bestImprovement_ = objChange;
	    }
	    totalSocoSolved_++;
            iterN++;
          }
	  
	}
	else if(cutRule_ == 10) {
	  // pure MOSEK
	  if(activeNode->ID==0) {
	    solveLP(activeNode,0); // Not a relaxation
	    // Parse info
	    // Node: MSK_IINF_MIO_NUM_RELAX
	    int totnode = 0, pronode = 0;
	    MSK_getintinf (activeNode->problem,  MSK_IINF_MIO_NUM_RELAX , &pronode);
	    MSK_getintinf (activeNode->problem,  MSK_IINF_MIO_NUM_BRANCH, &totnode);
	    //printf("%d\n",totnode);
	    nodesProcessed_ = pronode;
	    totalNodes_ = totnode;
	    totalSocoSolved_ = pronode;
	  }
	}
	else if(cutRule_ == 11) {
	  // MOSEK with DCCs at ROOT NODE
	  if(activeNode->ID==0) { // If it is root node
          
            for(int i=0; i<N; i++) {
              int isNewCut = cut(activeNode);
              activeNode->totalCuts += isNewCut;
              totalCutsGenerated_ ++;
              totalCutsApplied_ += isNewCut;
            }
            solveLP(activeNode,0); // Not a relaxation
	    int totnode = 0, pronode = 0;
	    MSK_getintinf (activeNode->problem,  MSK_IINF_MIO_NUM_RELAX , &pronode);
	    MSK_getintinf (activeNode->problem,  MSK_IINF_MIO_NUM_BRANCH, &totnode);
	    //printf("%d\n",totnode);
	    nodesProcessed_ = pronode;
	    totalNodes_ = totnode;
	    totalSocoSolved_ = pronode;
          }
	}
	else {
          printText(3, "Undefined cutting strategy. Please check readme.");
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
        printText(1, "New bound   %5d   %10.3e     %5.2f %%     %10.3e  %4d   %4d  /  %4d ", activeNode->ID, globalUpperBound_, curGap, lowBnd, lowNode, totalCutsApplied_, totalCutsGenerated_);
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

 summary_report:  
  clock_t end = clock();
  double elapsed = double(end-begin) / CLOCKS_PER_SEC;


  // Summary report
  printText(1, "========== SUMMARY ==========");
  if(globalUpperBound_== std::numeric_limits<double>::infinity() && terminationReason_ == 0) {
    printText(1, "Problem is infeasible.");
    return 0;
  }
  else if(terminationReason_ == 2) {
    printText(1, "Time limit has been hit.");
  }
  
  
  printText(1, "Best objective value: %e", globalUpperBound_);
  printText(2, "Values: ");
  for(int i=0; i<N; i++) {
    if(bestSoln_[i] < integerTolerance_) {
      printText(3, "  %2d: %.6e\t%.6e", (i+1), bestSoln_[i], bestNode->mosekSoln[i+1]);
    } else {
      printText(2, "  %2d: %.6e\t%.6e", (i+1), bestSoln_[i], bestNode->mosekSoln[i+1]);
    }
  }
  // RISK and RETURN values for portfolio problems
  int numcons = 0;
  MSK_getmaxnumcon(bestNode->problem, &numcons);
  printText(2, "Number of constraints: %d", numcons);
  for(int i=0; i<numcons; i++) {
    MSKstakeye c_sk;
    MSKrealt c_sx;
    MSKrealt c_lb;
    MSKrealt c_ub;
    MSKrealt c_sn;
    char c_string[120];
    MSK_getsolutioni(bestNode->problem, MSK_ACC_CON, i, MSK_SOL_ITR, &c_sk, &c_sx, &c_lb, &c_ub, &c_sn);
    MSK_sktostr(bestNode->problem, c_sk, c_string);
    printText(4, "Constraint [%d]: %.6f - %s (sl: %.6f, su: %.6f, sn: %.6f)", i, c_sx, c_string, c_lb, c_ub, c_sn);
  }
  
  for(int i=0; i<2*N+1; i++) {
    printText(7,"MOSEK Var[%d]:\t%.6f", i, bestNode->mosekSoln[i+1]);
  }
  printText(1, "Number of nodes processed: %d", nodesProcessed_);
  printText(1, "Number of nodes generated: %d", totalNodes_);
  printText(1, "Total SOCO solved: %d", totalSocoSolved_);
  printText(1, "Total time elapsed: %f seconds",  elapsed );
  printText(1, "Time per SOCO subproblem: %f", ((0.0+elapsed)/totalSocoSolved_)  );
  printText(1, "Total cuts generated: %d", totalCutsGenerated_);
  printText(1, "Total cuts applied: %d", totalCutsApplied_);
  if(cutRule_==2)
    printText(1, "Average effective cuts: %f", (double) totalFadingCuts_ / totalNodeFadingCuts_);
  printText(1, "Total objective improvement: %e", totalImprovement_);
  printText(1, "Best objective improvement: %e", bestImprovement_);  
  printText(1,"Optimal node: %d", bestNodeNumber_);
  printText(1,"Optimal node depth: %d", bestNode->depth);
  printText(1,"Optimal value: %e", globalUpperBound_);
  double finalgap_;
  double t1;
  int t2;
  getCurrentGap(&finalgap_, &t1, &t2);
  printText(1,"Final gap: %5.2f %%", finalgap_);
  for(unsigned int i=0; i<cutImprovements_.size(); i++) {
    printText(3, "Cut[ %4d ], Depth: %8Le, Tree Depth: %Lf, Improvement: %8Le, CurrentObj: %8Le, Variable: %Lf, Value: %Lf",(i+1), cutImprovements_[i][0], cutImprovements_[i][1], cutImprovements_[i][2], cutImprovements_[i][3], cutImprovements_[i][4], cutImprovements_[i][5]);
  }
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

  branchingVarValue_ = activeNode->nodeSoln[asset];
  printText(2,"Branching at [Node %04d], on var(%02d), <=%2.0f [Node %04d] and >=%2.0f [Node %04d], Obj: %.5f, Cur.Var.Val: %e, ActiveNodes: %ld",activeNode->ID, asset, floor(activeNode->nodeSoln[asset]), totalNodes_, ceil(activeNode->nodeSoln[asset]), totalNodes_+1, activeNode->nodeObj, activeNode->nodeSoln[asset],nodeList_.size());
  
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
  
  // Type 1: Generate in order specified by -c option, for x=1,2,3
  if(cutRule_ < 6 || cutRule_ == 11) {
    int variableForCut = nextCut(N, cutPriority_, activeNode->nodeSoln, &(activeNode->usedCuts));
    if(variableForCut >= 0 ) {
      
      double currsoln = activeNode->nodeSoln[variableForCut];
      if(PROBLEMCODE==1) { // TODO Remove this if condition in future - problem specific
	if(fabs(currsoln-round(currsoln))<1e-6) { // TODO make this one a parameter
	  printText(3,"Node (%d) variable %d is close to integer %.3f, cut is not generated",activeNode->ID,variableForCut,currsoln);
	  return 0;
	} 
      }
      int cutIndex = 0;
      int addToProblem  = 1;
      status = addNewCut(activeNode->problem, &cutIndex, activeNode->nodeSoln, activeNode->nodeObj, variableForCut+1, activeNode->nodeSoln[variableForCut], cutSelection_, &cdepth_, addToProblem);
      if(status>0) {
	printText(3, "Node (%d) New cut for variable %d is added, value: %f", activeNode->ID, variableForCut, activeNode->nodeSoln[variableForCut]);
	cvariable_ = variableForCut;
	cvarvalue_ = activeNode->nodeSoln[variableForCut];
      } else {
	printText(3, "Node (%d) Cut generation failed, variable %d, value: %f", activeNode->ID, variableForCut, activeNode->nodeSoln[variableForCut]);
      }
    }
    else {
      printText(3,"No cuts can be added, all integer cuts are active for current values...");
    }
  }
  else if(cutRule_ == 6) { // Generate all possible cuts (N) and order them by violation, add top -p of them if possible
    int i = 0;
    double* cutDepthIndex = (double*) malloc(sizeof(double)*N);
    for(i=0; i<N; i++) {
      cdepth_ = -1e-8;
      double currsoln = activeNode->nodeSoln[i];
      if(fabs(currsoln-round(currsoln))< integerTolerance_) {
	printText(3,"Node (%d) variable %d is close to integer %.3f, cut is not generated",activeNode->ID,i,currsoln);
      }
      else {
	int response = 0;
	int addToProblem = 0; // DO NOT add new cut to the problem, yet
	int cutIndex = 0;
	
	response = addNewCut(activeNode->problem, &cutIndex, activeNode->nodeSoln, activeNode->nodeObj, i+1, activeNode->nodeSoln[i], cutSelection_, &cdepth_, addToProblem);

	if(response > 0) {
	  printText(4, "Node (%d) New cut for variable %d is generated, value: %f", activeNode->ID, i, activeNode->nodeSoln[i]);
	}
	else {
	  printText(4, "Node (%d) New cut for variable %d is failed, value: %f", activeNode->ID, i, activeNode->nodeSoln[i]);
	}

      }
      cutDepthIndex[i] = cdepth_;
      
    }

    // All cuts are generated, now order and add
    int effCutIter = 0;
    while(effCutIter < cutPerIteration_) {
      int biggestI = 0;
      double biggestV = cutDepthIndex[0];
      for(i=0; i<N; i++) {
	if(cutDepthIndex[i] >= biggestV) {
	  biggestI = i;
	  biggestV = cutDepthIndex[i];
	}
      }
      printText(6, "Max depth: %e, index: %d", biggestV, biggestI);
      cutDepthIndex[biggestI] = -1e-8; // Prevent it being added again
      if(biggestV >= 0) { // valid cut
	int response = 0;
	int addToProblem = 1;
	int cutIndex = 0;
	response = addNewCut(activeNode->problem, &cutIndex, activeNode->nodeSoln, activeNode->nodeObj, biggestI+1, activeNode->nodeSoln[biggestI], cutSelection_, &cdepth_, addToProblem);
	status += response;
	printText(3, "Node (%d) New cut for variable %d is added, value: %f", activeNode->ID, biggestI, activeNode->nodeSoln[biggestI]);
      }
      else {
	break; // break while loop as no cut is effective anymore
      }
    }
    
    free(cutDepthIndex);
    
  }
  
  return status;
}

int selectNode(Node** activeNode) {
  int status = 0;
  unsigned int selected = 0;
  int max_depth = -1;

  double lowestlb = nodeList_[0]->lowerBound;
  for(unsigned int i=0; i < nodeList_.size(); i++) {
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
  else if(searchRule_ == 4) { // dfc depth-first closest
    *activeNode = nodeList_[0];
    max_depth = nodeList_[0]->depth;
    for(unsigned int i=1; i < nodeList_.size(); i++) {
      if ( (nodeList_[i]->depth >  max_depth && nodeList_[i]->lowerBound < redline) ||
	   (nodeList_[i]->depth >= max_depth && nodeList_[i]->lowerBound < redline && branchingVarValue_ - floor(branchingVarValue_) >= 0.5) ) {
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
      if(nodeList_[i]->lowerBound < bestbound && nodeList_[i]->lowerBound < redline) {
 	//printText(0,"    Better bound, node[%d] = %f < %f",i,nodeList_[i]->lowerBound,bestbound);
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

  if (totaln>0) {
    for(unsigned int i=nodeList_.size(); i>=1; i--) {
      if(nodeList_[i-1]->eliminated) {
	printText(4, "Node [%d] (pos:%d, obj:%e) is removed, ActiveNodes: %ld", nodeList_[i-1]->ID, i-1, nodeList_[i-1]->lowerBound, nodeList_.size());
	nodeList_.erase(nodeList_.begin() + (i-1));
      }
    }
  }
  
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
    (*newNode)->lowerBound = -999999999;
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

      if(PROBLEMCODE==3) { // HOTFIX FOR SINGLE BOUND CONSTRAINTS DUE TO ERRORS
	if(lower==0) {
	  MSK_chgbound ( 
			newProblem,         //MSKtask_t    task, 
			MSK_ACC_VAR,        //MSKaccmodee  accmode, 
			varID+1,            //MSKint32t    i, 
			lower,              //MSKint32t    lower, 
			1,                  //MSKint32t    finite, 
			bound               //MSKrealt     value); 
			 );
	}
      }
      
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
    printText(7,"%d\t%e",nodeList_[i]->ID, nodeList_[i]->lowerBound);
    if(!nodeList_[i]->eliminated && !nodeList_[i]->processed) {
      if(nodeList_[i]->lowerBound < lb) {
	lb = nodeList_[i]->lowerBound;
	*lowerNode = nodeList_[i]->ID;
	found = 1;
      }
    }
  }

  printText(6,"Lowest lb found: %d, Node: %d, Value: %e\n", found, *lowerNode, lb);

  if(found) {
    *lowerBound = lb;
    *curGap = (globalUpperBound_ - lb) / fabs(globalUpperBound_) * 100;
  }
  else {
    *lowerBound = globalUpperBound_;
    *curGap = 0;
    *lowerNode = 0;
  }

  return 1;

}

int solveLP(Node* aNode, int relax) {
  int status = 0;
  
  MSKtask_t originaltask = aNode->problem;

  MSKtask_t mtask;
  //MSK_clonetask(originaltask, &mtask);
  mtask = originaltask;

  // 0 MSK_putintparam(mtask, MSK_IPAR_NUM_THREADS, 1);
  
  if(relax==1) {

    
    MSK_putintparam(mtask, MSK_IPAR_MIO_MODE, MSK_MIO_MODE_IGNORED); // make it LP
    MSK_putintparam(mtask, MSK_IPAR_INTPNT_REGULARIZATION_USE, MSK_OFF);
    MSK_putintparam(mtask, MSK_IPAR_INTPNT_SCALING, MSK_SCALING_NONE);
    MSK_putintparam(mtask, MSK_IPAR_INTPNT_BASIS, 0);
    MSK_putintparam(mtask, MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_CONIC);

    //MSK_putintparam(mtask, MSK_IPAR_PRESOLVE_LINDEP_USE, MSK_OFF);
    MSK_putintparam(mtask, MSK_IPAR_PRESOLVE_USE, MSK_PRESOLVE_MODE_OFF);
    // xxMSK_putdouparam(mtask, MSK_DPAR_INTPNT_CO_TOL_INFEAS, 1e-16);
    // xxMSK_putdouparam(mtask, MSK_DPAR_INTPNT_CO_TOL_PFEAS, 1e-16);
    //MSK_putdouparam(mtask, MSK_DPAR_INTPNT_TOL_PFEAS, 1e-16);
    //  xxMSK_putdouparam(mtask, MSK_DPAR_INTPNT_CO_TOL_DFEAS, 1e-16);
    //  xxMSK_putdouparam(mtask, MSK_DPAR_INTPNT_CO_TOL_REL_GAP, 1e-16);
    //  xxMSK_putdouparam(mtask, MSK_DPAR_CHECK_CONVEXITY_REL_TOL, 100000);


    
    MSK_putintparam(mtask, MSK_IPAR_INTPNT_MAX_ITERATIONS, 30000);
    MSK_putintparam(mtask, MSK_IPAR_BI_IGNORE_MAX_ITER, MSK_ON);

    
    // QP PARAMS!
    MSK_putdouparam(mtask, MSK_DPAR_INTPNT_TOL_DFEAS, 1e-16);
    MSK_putdouparam(mtask, MSK_DPAR_INTPNT_NL_TOL_DFEAS, 1e-16);
    MSK_putdouparam(mtask, MSK_DPAR_INTPNT_NL_TOL_PFEAS, 1e-16);
    //if(aNode->ID == 0) {
    //MSK_linkfunctotaskstream (mtask, MSK_STREAM_LOG, NULL, printstr);
    //}
  }
  else {
    char funame[80];
    sprintf(funame, "../test/result/int%d.mps", aNode->ID);
    MSK_writedata(mtask, funame);
    MSK_readdata (mtask, funame);
    MSK_putintparam(mtask, MSK_IPAR_MIO_MODE, MSK_MIO_MODE_SATISFIED);
    MSK_putintparam(mtask, MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_FREE);
    double dtol = 1e-6;
    MSK_putdouparam(mtask, MSK_DPAR_MIO_NEAR_TOL_ABS_GAP, dtol);
    MSK_putdouparam(mtask, MSK_DPAR_MIO_TOL_ABS_GAP, dtol);
    MSK_putdouparam(mtask, MSK_DPAR_MIO_TOL_ABS_RELAX_INT,dtol);
    MSK_putdouparam(mtask, MSK_DPAR_MIO_TOL_FEAS, dtol);
    //MSK_putdouparam(mtask, MSK_DPAR_MIO_TOL_REL_RELAX_INT, dtol);       // commented for 8!
    //MSK_linkfunctotaskstream (mtask, MSK_STREAM_LOG, NULL, printstr);
  }

  // MOSEK PRECISION
  //MSK_putdouparam(mtask, MSK_DPAR_INTPNT_CO_TOL_PFEAS, 1e-8);
  //MSK_putdouparam(mtask, MSK_DPAR_INTPNT_CO_TOL_DFEAS, 1e-16);
  //MSK_putdouparam(mtask, MSK_DPAR_INTPNT_CO_TOL_REL_GAP, 1e-16);
  //MSK_putdouparam(mtask, MSK_DPAR_INTPNT_TOL_INFEAS, 1e-16);
  //MSK_putdouparam(mtask, MSK_DPAR_INTPNT_CO_TOL_MU_RED, 1e-16);
  //MSK_putintparam(mtask, MSK_IPAR_NUM_THREADS, 1);
  

  int reattempted = 0;
  solutionpoint:
 
  MSKrescodee trmcode;
  MSKrescodee r;
  printToFile(aNode);
  clock_t mosek_start = clock();
  r = MSK_optimizetrm(mtask,&trmcode);
  printText(3,"Term code: %d", trmcode);
  clock_t mosek_finish = clock();
  double elapsed = double(mosek_finish-mosek_start) / CLOCKS_PER_SEC;
  printText(6,"Rel. sol time %.3e",elapsed);
  printText(6,"Problem is solved with  MOSEK");
  r = r;
  MSK_solutionsummary(mtask,MSK_STREAM_LOG);
  
  MSKsolstae solsta;
  MSKprostae prosta;
  MSKsoltypee whichsol;
  if(relax) {  whichsol = MSK_SOL_ITR; }
  else      {  whichsol = MSK_SOL_ITG; }
  MSKrescodee obsolsta = MSK_getsolsta (mtask, whichsol, &solsta);
  MSK_getprosta(mtask, whichsol, &prosta);
  char solsta_str[50];
  char prosta_str[80];
  MSK_solstatostr(mtask, solsta, solsta_str);
  MSK_prostatostr(mtask, prosta, prosta_str);
  printText(3,"Problem status: %d %s", prosta, prosta_str);
  printText(3,"Solution status: %d %d %s", solsta, trmcode, solsta_str);
  MSKrealt primalobj;

  // hack for mosek problems
  //  if(solsta==MSK_SOL_STA_UNKNOWN && relax==0) {
    
  //  goto resolve;
  //}

  //printf("TRMCODE: %d\n", trmcode);
  //printf("SOLSTA: %d\n", solsta);
  if(trmcode==10006) {
    solsta = MSK_SOL_STA_OPTIMAL;
    //printf("Stall!\n");
  }

  switch(prosta) {
  case MSK_PRO_STA_PRIM_AND_DUAL_FEAS:
  case MSK_PRO_STA_PRIM_FEAS:
  case MSK_PRO_STA_NEAR_PRIM_FEAS:
    goto cpsolutionfeas;
  case MSK_PRO_STA_PRIM_INFEAS:
  case MSK_PRO_STA_PRIM_AND_DUAL_INFEAS:
  case MSK_PRO_STA_ILL_POSED:
  case MSK_PRO_STA_PRIM_INFEAS_OR_UNBOUNDED:
    goto unknowncase;
  default:
    printText(6,"Problem status: %s", prosta_str);
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
      cpsolutionfeas:
	printText(3,"Node status %d, term code: %d",solsta,trmcode);
	aNode->feasible = true;
	MSK_getprimalobj(mtask, whichsol, &primalobj);
	aNode->nodeObj = primalobj;
	if(primalobj > aNode->lowerBound) {
	  aNode -> lowerBound = primalobj;
	}
	//aNode -> nodeSoln = (double*) malloc(N*sizeof(double));
	if(PROBLEMCODE==1) {
	  MSK_getsolutionslice(mtask, whichsol, MSK_SOL_ITEM_XX, 1, N+1, aNode->nodeSoln);
	} else if(PROBLEMCODE==2) {
	  MSK_getsolutionslice(mtask, whichsol, MSK_SOL_ITEM_XX, N+1, 2*N+1, aNode->nodeSoln);
	} else if(PROBLEMCODE==3 || PROBLEMCODE==4) {
	  MSK_getsolutionslice(mtask, whichsol, MSK_SOL_ITEM_XX, N+1, 2*N+1, aNode->nodeSoln);
	}
	MSK_getsolutionslice(mtask, whichsol, MSK_SOL_ITEM_XX, 0, 2*N+1, aNode->mosekSoln);
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
	if(reattempted) {
	  printText(3,"Node status unknown %d, term code: %d, pruning, attempt failed!",solsta,trmcode);
	  try {
	    printText(3,"Trying to recover solution...");
	    MSK_getprimalobj(mtask, whichsol, &primalobj);
	    printText(3,"Recovered objective: %e", primalobj);
	    double* currsoln = (double*) malloc(N*sizeof(double));
	    MSK_getsolutionslice(mtask, whichsol, MSK_SOL_ITEM_XX, 1, N+1, currsoln);
	    double sumofsoln = 0;
	    for (int i=0; i<N; i++) {
	      //printf("X[%d]: %e\n", i, currsoln[i]);
	      sumofsoln += fabs(currsoln[i]);
	    }
	    free(currsoln);
	    if(sumofsoln < 1e-2 && PROBLEMCODE==1) {
	      printText(3,"Sum of solution is %e, result is unknown...", sumofsoln);
	      goto unknowncase;
	    }
	    if(fabs(primalobj) >= 1e-3) {
	      printText(3, "Non-zero objective, problem could be feasible. Continuing to get the solution");
	      goto cpsolutionfeas;
	    }
	    
	  }
	  catch (...) {
	    printText(3, "Exception... continuing with prune operation...");
	  }
	}
	else {
	  printText(4,"Node status unknown %d, term code: %d, reattempting with lower precision", solsta, trmcode);
	  MSK_putdouparam(mtask, MSK_DPAR_INTPNT_CO_TOL_INFEAS, 1e-3);
	  MSK_putdouparam(mtask, MSK_DPAR_INTPNT_CO_TOL_PFEAS, 1e-3);
	  MSK_putdouparam(mtask, MSK_DPAR_INTPNT_CO_TOL_DFEAS, 1e-3);
	  MSK_putdouparam(mtask, MSK_DPAR_INTPNT_CO_TOL_REL_GAP, 1e-3);
  
	  reattempted = 1;
	  goto solutionpoint;
	}
	
	unknowncase:
	aNode->feasible = false;
	
	if(FILEOUTPUT) {
	  char funame[80];
	  sprintf(funame, "../test/result/unknown%d.mps", aNode->ID);
	  MSK_writedata(mtask, funame);
	}
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
      printText(5,"Integer infeasibility at variable %d, value: %e, inttol: %e",i, aNode->nodeSoln[i], integerTolerance_);
    }
    printText(6,"Proximity to integer (%d - asset %d -): %.6f, value: %.6f",i,i+1,fabs(round(aNode->nodeSoln[i])-aNode->nodeSoln[i]), aNode->nodeSoln[i]);
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

  postSolve(bestSoln_);

  for(unsigned int i=0; i < allNodeList_.size(); i++) {
    deleteNode(allNodeList_[i]);
  }

  
  deleteProblem();

  for(unsigned int i=0; i < cutImprovements_.size(); i++) {
    free(cutImprovements_[i]);
  }
  
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
  char fname2[80];
  sprintf(fname, "../test/result/node%d.mps", fno);
  sprintf(fname2, "../test/result/node%d.lp", fno);

  //MSK_toconic(&(aNode->problem));

  if(FILEOUTPUT) { 
    MSK_writedata(aNode->problem, fname);
    MSK_writedata(aNode->problem, fname2);
  }
  
  return 1;
  
}

/**
 * @short Solves primal and dual rounding problems and returns solution
 * @param[in] aNode  Node whose rounding will be solved
 * @return Success code (0) or Failure (1)
 */
int solveRounding(Node* aNode) {
  
  // Step 0: Initialize memory and variables
  
  // Step 1: Grab cone structure, optimal solution  and find conic classes
  
  // Step 2: Get eigenvectors
  
  // Step 3: Solve primal rounding problem
  
  // Step 4: Solve dual rounding problem
  
  // Step 5: Convert results back into original space
  
  
  // Return success code!
  return 0;
  
}

/**
 * @short Detects structures in a given conic solution
 * @param[in] aNode  Subject to change
 * @return Success code (0) or Failure (1)
 */
int getConicClasses(Node* aNode) {




  return 0;
}








