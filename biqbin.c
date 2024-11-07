#include "biqbin.h"  


#define HEAP_SIZE 10000000

extern Heap *heap;
extern FILE *output;
extern BiqBinParameters params;

                      
int compute(InputData input_data, BiqBinParameters biqbin_parameters) {

    srand(2024);
    params = biqbin_parameters;

    if (output == NULL) {
        open_output_file(input_data.name);
    }

    post_process_BQP_input(input_data);
    print_parameters(params);
 
    BabNode *node;

    /*** allocate priority queue ***/
    heap = Init_Heap(HEAP_SIZE);

    Bab_Init();

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