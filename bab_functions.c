#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>       
#include <math.h>
#include <sys/time.h>

#include "biqbin.h"
  
/* defined in heap.c */  
extern int BabPbSize;           
extern BabSolution *BabSol;

extern FILE *output;    
extern BiqBinParameters params;
extern Problem *SP;
extern Problem *PP;
extern double root_bound;
extern double TIME;
extern int stopped;


/* initialize global lower bound to 0 and global solution vector to zero */
void initializeBabSolution(void) { 

    BabSolution bs;

    for (int i = 0; i < BabPbSize; ++i) {
        bs.X[i] = 0;
    }

    Bab_LBInit(0, &bs);
}


/************** Initialization: root node and priority queue **************/
void Init_PQ(void) {
    
    extern BabNode *BabRoot;

    // Create the root node
    BabRoot = newNode(NULL);

    // increase number of evaluated nodes
    Bab_incEvalNodes();

    // Evaluate root node: compute upper and lower bound 
    root_bound = Evaluate(BabRoot, SP, PP);

    // save upper bound
    BabRoot->upper_bound = root_bound; 

    /* insert node into the priority queue or prune */
    // NOTE: optimal solution has INTEGER value, i.e. add +1 to lower bound
    if (Bab_LBGet() + 1.0 < BabRoot->upper_bound) {      
        Bab_PQInsert(BabRoot); 
    }
    else {
        // otherwise, upper_bound <= BabLB, so we can prune
        free(BabRoot);
    }
}


/* 
 * Bab function which initializes the problem, allocate the structures 
 * and evaluate the root node.
 */
void Bab_Init(int argc, char **argv) {

    // Start the timer
    TIME = time_wall_clock();

    // Process the command line arguments
    processCommandLineArguments(argc, argv);

    // Provide B&B with an initial solution
    initializeBabSolution();

    // Allocate the memory
    allocMemory();

    // Initialize root node
    Init_PQ(); 
}

/* NOTE: int *sol in functions evaluateSolution and updateSolution have length BabPbSize
 * -> to get objecive multiple with Laplacian that is stored in upper left corner of SP->L
 */
double evaluateSolution(int *sol) {

    double val = 0.0;
    
    for (int i = 0; i < BabPbSize; ++i) {
        for (int j = 0; j < BabPbSize; ++j) {
            val += SP->L[j + i * SP->n] * sol[i] * sol[j];
        }
    }

    return val;
}


/*
 * Only this function can update best solution and value.
 * Returns 1 if success.
 */
int updateSolution(int *x) {
    
    int solutionAdded = 0;
    double sol_value;
    BabSolution solx;

    // Copy x into solx --> because Bab_LBUpd needs BabSolution and not int*
    for (int i = 0; i < BabPbSize; ++i) {
      solx.X[i] = x[i];
    }

    sol_value = evaluateSolution(x); // computes objective value of solx

    /* If new solution is better than the global solution, 
     * then update and print the new solution. */
    if (Bab_LBUpd(sol_value, &solx)) {
        solutionAdded = 1;
        printf("Node %d Feasible solution %.0lf\n", Bab_numEvalNodes(), Bab_LBGet() );
        fprintf(output,"Node %d Feasible solution %.0lf\n", Bab_numEvalNodes(), Bab_LBGet() );
    }
    
    return solutionAdded;
}


/* 
 * The function generates two new children of the current node of the 
 * branch-and-bound tree.
 * It also evaluates new generated nodes, i.e. computes upper and lower bound.
 */
void Bab_GenChild(BabNode *node) {

    BabNode *child_node;

    // If the algorithm stops before finding the optimal solution, search in the 
    // nodes queue for the worst upper bound. 
    if (stopped || params.root || (params.time_limit > 0 && (time_wall_clock() - TIME) > params.time_limit) ) {
        
        // signal to printFinalOutput that algorithm stopped early
        stopped = 1;        

        free(node);
        return;
    }

    // If it's a leaf of the B&B tree, add the solution and don't branch
    if (isLeafNode(node)) {
        updateSolution(node->sol.X);
        free(node);
        return;
    }

    // Determine the variable x[ic] to branch on
    int ic = getBranchingVariable(node);

    if (params.detailedOutput) {
        fprintf(output, "Branching on x[%d] = %.2f\n", ic, node->fracsol[ic]);
    }    

    // add two nodes to the search tree
    for (int xic = 0; xic <= 1; ++xic) { 

        // Create a new child node from the parent node
        child_node = newNode(node);

        // split on node ic
        child_node->xfixed[ic] = 1;
        child_node->sol.X[ic] = xic;

        if (params.detailedOutput) {
            fprintf(output, "Fixing x[%d] = %d\n", ic, xic);
        }
            
        //increment the number of explored nodes
        Bab_incEvalNodes();

        // If it's a leaf of the B&B tree, add the solution and don't branch
        if (isLeafNode(child_node)) {
            updateSolution(child_node->sol.X);
            free(child_node);
        }
        else {

            /* compute upper bound (SDP bound) and lower bound (via heuristic) for this node */
            child_node->upper_bound = Evaluate(child_node, SP, PP);

            /* if BabLB + 1.0 < child_node->upper_bound, 
             * then we must branch since there could be a better feasible 
             * solution in this subproblem
             */
            if (Bab_LBGet() + 1.0 < child_node->upper_bound) {
                /* insert node into the priority queue */
                Bab_PQInsert(child_node); 
            }
            else {
                // otherwise, upper_bound <= BabLB, so we can prune
                free(child_node);
            }
        }

    } // end for xic

    // free parent node
    free(node);
}


/* timer */
double time_wall_clock(void)  {
 
    struct timeval timecheck;
    gettimeofday(&timecheck, NULL);
    return timecheck.tv_sec + timecheck.tv_usec * 1e-6;

}


/* print solution 0-1 vector */
void printSolution(FILE *file) {

    fprintf(file, "Solution = ( ");
    for (int i = 0; i < BabPbSize; ++i) {
        if (BabSol->X[i] == 1) {
            fprintf(file, "%d ", i + 1);
        }
    }
    fprintf(file, ")\n");
}


/* print final output */
void printFinalOutput(FILE *file, int num_nodes) {

    // Best solution found
    double best_sol = Bab_LBGet();

    fprintf(file, "\nNodes = %d\n", num_nodes);
    
    // normal termination
    if (!stopped) {
        fprintf(file, "Root node bound = %.3lf\n", root_bound);
        fprintf(file, "Maximum value = %.0lf\n", best_sol);
        printSolution(file);
        
    } else { // B&B stopped early
        if (params.root) {
            fprintf(file, "Root node bound = %.3lf\n", root_bound);
            fprintf(file, "Best value = %.0lf\n", best_sol);
        }
        else { /* time limit reached */
            fprintf(file, "TIME LIMIT REACHED.\n");
            fprintf(file, "Root node bound = %.3lf\n", root_bound);
            fprintf(file, "Best value = %.0lf\n", best_sol);
        }    
        printSolution(file);
    }

    fprintf(file, "Wall clock time = %.2f s\n\n", time_wall_clock() - TIME);
}


/* Bab function called at the end of the execution.
 * This function prints output: number of 
 * nodes evaluated, the best solution found, the nodes best bound, the wall clock time.
 * It is also used to free the memory allocated by the program.
 */
void Bab_End(void) {

    /* Print results to the standard output and to the output file */
    printFinalOutput(stdout,Bab_numEvalNodes());
    printFinalOutput(output,Bab_numEvalNodes());

    freeMemory();

    fclose(output);
}


/*
 * getBranchingVariable function used in the Bab_GenChild routine to determine
 * which variable x[ic] to branch on.
 *
 * node: the current node of the branch-and-bound search tree
 */
int getBranchingVariable(BabNode *node) {

    int ic = -1;  // x[ic] is the variable to branch on
    double maxValue, minValue;

    /* 
     * Choose the branching variable x[ic] based on params.branchingStrategy
     */
    if (params.branchingStrategy == LEAST_FRACTIONAL) {
        // Branch on the variable x[ic] that has the least fractional value
        maxValue = -BIG_NUMBER;
        for (int i = 0; i < BabPbSize; ++i) {
            if (!(node->xfixed[i]) && fabs(0.5 - node->fracsol[i]) > maxValue) {
                ic = i;
                maxValue = fabs(0.5 - node->fracsol[ic]);
            }
        }
    }
    else if (params.branchingStrategy == MOST_FRACTIONAL) {
        // Branch on the variable x[ic] that has the most fractional value
        minValue = BIG_NUMBER;
        for (int i = 0; i < BabPbSize; ++i) {
            if (!(node->xfixed[i]) && fabs(0.5 - node->fracsol[i]) < minValue) {
                ic = i;
                minValue = fabs(0.5 - node->fracsol[ic]);
            }
        }
    }
    else {
        fprintf(stderr, "Error: Wrong value for params.branchingStrategy\n");
        exit(1);
    }

    return ic;
}


/* Count the number of fixed variables */
int countFixedVariables(BabNode *node) {
    
    int numFixedVariables = 0;

    for (int i = 0; i < BabPbSize; ++i) {
        if (node->xfixed[i]) {
            ++numFixedVariables;
        }
    }

    return numFixedVariables;
}

/* Determine if node is a leaf node by counting the number of fixed variables */
int isLeafNode(BabNode *node) {
    return (countFixedVariables(node) == BabPbSize);
}
