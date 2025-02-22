#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "biqbin.h"

extern double *X;
extern double *Z;       // stores Cholesky decomposition: X = ZZ^T

double runHeuristic(Problem *P0, Problem *P, BabNode *node, int *x) {

    // Problem *P0 ... the original problem
    // Problem *P  ... the current subproblem
    // int *x      ... current best feasible solution 

    int n = P->n;
    int inc = 1;
    char UPLO = 'L';
    int info = 0;
    int nn = n * n;
    double heur_val;
    // Z = X
    dcopy_(&nn, X, &inc, Z, &inc);

    // compute Cholesky factorization
    dpotrf_(&UPLO, &n, Z, &n, &info);

    // set lower triangle of Z to zero
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < i; ++j)
            Z[j + i*n] = 0.0;

    // Goemans-Williamson heuristic
    heur_val = GW_heuristic(P0, P, node, x, P0->n);

    return heur_val;    


#if 0
    int n = P->n;
    int N = P0->n - 1; // BabPbSize
    int nn = n * n;
    int inc = 1;
    char UPLO = 'L';
    int info = 0;
    double heur_val;

    double xh[n];       // is used for convex combination with matrix X (n x n)
    int temp_x[N];      // stores xh + some variables are fixed in {0,1} model
    int index;
    
    // generate first random cut vector {-1,1}^n
    for (int i = 0; i < n; ++i)
        xh[i] = 2 * (rand() % 2) - 1; 

    // compute its objective value (store in temp_x and transform to {0,1})
    index = 0;
    for (int i = 0; i < N; ++i) {

        if (node->xfixed[i]) 
            temp_x[i] = node->sol.X[i];

        else {
            temp_x[i] = (xh[index] + 1) / 2.0;
            ++index;
        }
    }

    double fh = evaluateSolution(temp_x);

    int done = 0;
    double constant;    // scalar in convex combiantion
    double alpha;

    // Z = X
    dcopy_(&nn, X, &inc, Z, &inc);

    while (done < 2) {

        ++done;
        
        // compute Cholesky factorization
        dpotrf_(&UPLO, &n, Z, &n, &info);

        if (info != 0) {
            fprintf(stderr, "%s: Problem with Cholesky factorization \
                (line: %d).\n", __func__, __LINE__);
            exit(1);
        }

        // set lower triangle of Z to zero
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < i; ++j)
                Z[j + i*n] = 0.0;

        // Goemans-Williamson heuristic
        heur_val = GW_heuristic(P0, P, node, x, P0->n);

        if (heur_val > fh) {

            done = 0;
            fh = heur_val;

            // copy global cut vector x into xh
            // NOTE: skip fixed vertices
            index = 0;
            for (int i = 0; i < N; ++i) {

                if (!node->xfixed[i]) {
                    xh[index] = 2 * x[i] - 1;
                    ++index;
                }
            }
            xh[n-1] = -1.0;                  // last vertex in original is fixed to 0
                      
        }

        constant = 0.3 + 0.6 * ( (double)rand()/(double)(RAND_MAX) );

        // Z = (1-constant)*X + constant* xh *xh'
        alpha = 1.0 - constant;
        dcopy_(&nn, X, &inc, Z, &inc);
        dscal_(&nn, &alpha, Z, &inc);
        alpha = constant;
        dsyr_(&UPLO, &n, &alpha, xh, &inc, Z, &n);

    }

    return heur_val;
    #endif
}


/* Goemans-Williamson random hyperplane heuristic */
double GW_heuristic(Problem *P0, Problem *P, BabNode *node, int *x, int num) {

    // Problem *P0 ... the original problem
    // Problem *P  ... the current subproblem
    //         num ... number of random hyperplanes

    int index;
    int N = P->n;

    // (local) temporary vector of size X
    int temp_x[N];                    

    // (global) temporary vector of size BabPbSize to store heuristic solutions
    int sol[P0->n - 1];                 

    double sca;                         // dot product of random vector v and col of Z
    double best = -BIG_NUMBER;          // best lower bound found
    double v[N];                        // defines random hyperplane v   


    for (int count = 0; count < num; ++count) {

        // compute random hyperplane v
        for (int i = 0; i < N; ++i) 
            v[i] = ( (double)rand() / (double)(RAND_MAX) ) - 0.5;

        // compute cut temp_x generated by hyperplane v
        index = 0;
        for (int i = 0; i < N; ++i) {

                sca = 0.0;
                for (int j = 0; j < N; ++j)
                    sca += v[j] * Z[j * N + index];

                if (sca < 0) {

                    temp_x[i] = -1;
                    
                }
                else {

                    temp_x[i] = 1;
                    
                }

                ++index;            
        }

        // improve feasible solution through 1-opt
        mc_1opt(temp_x, P);

        // store local cut temp_x into global cut sol
        index = 0;
        for (int i = 0; i < P0->n-1; ++i) {
            if (node->xfixed[i]) 
                sol[i] = node->sol.X[i];
            else {
                sol[i] = (temp_x[index]+1)/2;
                ++index;
            }
        }

        update_best(x, sol, &best, P0);
      
    }

    return best;
}


/*
 * Performs a simple local search starting from the given feasible solution x. 
 * Returns a feasible solution x that is locally optimal.
 * The objective value of x is returned.
 */
// NOTE: this function is working in {-1,1} model!
double mc_1opt(int *x, Problem *P) {

    int N = P->n;

    double *Lx, *d, *delta;
    int *I;
    alloc_vector(Lx, N, double);
    alloc_vector(d, N, double);
    alloc_vector(delta, N, double);
    alloc_vector(I, N, int);


    // Lx = L*x
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            Lx[i] += P->L[j + i*N] * x[j];

    // d = diag(L);
    // cost = x'*Lx
    // delta = d - x.*Lx
    double cost = 0.0;
    
    for (int i = 0; i < N; ++i) {
        d[i] = P->L[i + i * N];
        cost += x[i] * Lx[i];
        delta[i] = d[i] - x[i] * Lx[i];
    }

    // [best, i] = max(delta);
    double best = -BIG_NUMBER;
    int index = 0;

    for (int i = 0; i < N; ++i) {

        if (delta[i] > best) {
            best = delta[i];
            index = i;
        }

    }

    int num_I;      // number of elements in I

    /*** main loop ***/
    while (best > 0.001) {


        // I = find(L(:,index))
        num_I = 0;
        for (int j = 0; j < N; ++j) {
            
            if ( fabs(P->L[index + N * j]) > 0.001 ) { // add to I
                I[num_I] = j;
                ++num_I;
            }
        
        }

        if (x[index] > 0) { // Lx(I) = Lx(I)  - 2 *L(I,index);
            for (int i = 0; i < num_I; ++i) {
                Lx[I[i]] -= 2 * P->L[index + I[i] * N];
            }
        }
        else { // Lx(I) = Lx(I)  + 2 *L(I,index);
            for (int i = 0; i < num_I; ++i) {
                Lx[I[i]] += 2 * P->L[index + I[i] * N];
            }
        }

        // update new cut: x(index) = -x(index) 
        x[index] *= -1;

        // update weight of cut: cost = cost + 4*best
        cost += 4 * best;

        // update new differences: delta = d - x.*Lx
        for (int i = 0; i < N; ++i) {
            delta[i] = d[i] - x[i] * Lx[i];
        }

        // find new champion: [best, i] = max(delta) 
        best = -BIG_NUMBER;
        index = 0;

        for (int i = 0; i < N; ++i) {

            if (delta[i] > best) {
                best = delta[i];
                index = i;
            }

        }

    }

    free(Lx);
    free(d);
    free(delta);
    free(I);

    return cost;
}


/*
 * Given the current best solution, xbest, and a new solution, xnew, determines
 * the objective value of xnew, then replaces xbest with xnew if
 * xnew is better. Also updates the best objective value, best.
 */
int update_best(int *xbest, int *xnew, double *best, Problem *P0) {

    int success = 0;
    int N = P0->n - 1; // N = BabPbSize

    double heur_val = evaluateSolution(xnew);

    if ( *best < heur_val ) {
        memcpy(xbest, xnew, sizeof(int) * N);
        *best = heur_val;
        success = 1;
    }

    return success;
}

