#include <stdio.h>
#include <stdlib.h>

#include "biqbin.h"  

#define HEAP_SIZE 10000000
extern Heap *heap;

                      
int main(int argc, char **argv) {

    // Seed the random number generator
    srand(2024);

    BabNode *node;

    /*** allocate priority queue ***/
    heap = Init_Heap(HEAP_SIZE);

    /*
     * reads graph and params, initializes B&B solution,
     * evaluates root node and places it in priority queue
     * if not able to prune
     */
    Bab_Init(argc, argv);

    while (!isPQEmpty()) {
        node = Bab_PQPop();
        Bab_GenChild(node);
    }

    /* prints solution and frees memory */
    Bab_End();

    free(heap->data);
    free(heap);

    return 0;
    
}
