#include "biqbin.h"  


extern FILE *output;
                      
int main(int argc, char **argv) {
    InputData input_data;
    BiqBinParameters biqbin_parameters;

    // Process the command line arguments
    input_data = processCommandLineArguments(argc, argv);
    
    /*** Read the parameters from a user file ***/
    biqbin_parameters = readParameters(argv[2]);

    compute(input_data, biqbin_parameters);

    free(input_data.A);
    free(input_data.b);
    free(input_data.F);
    free(input_data.c);

    return 0;
    
}
