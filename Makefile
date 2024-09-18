#### BiqBin makefile ####

# Directories
OBJ = Obj

# Compiler: other options (linux users) --> gcc, icc
CC = gcc

# NOTE: -framework Accelerate is for MAC, Linux users set to -lopenblas -lm (or use intel mkl)
LINALG 	 = -lopenblas -lm 
OPTI     = -O3 

# binary
BINS =  biqbin

# BiqBin objects
BBOBJS = $(OBJ)/bundle.o $(OBJ)/allocate_free.o $(OBJ)/bab_functions.o \
	 $(OBJ)/bounding.o $(OBJ)/cutting_planes.o \
         $(OBJ)/evaluate.o $(OBJ)/heap.o $(OBJ)/ipm_mc_pk.o \
         $(OBJ)/heuristic.o $(OBJ)/main.o $(OBJ)/operators.o \
         $(OBJ)/process_input.o $(OBJ)/qap_simulated_annealing.o

# All objects
OBJS = $(BBOBJS)

CFLAGS = $(OPTI) -Wall -W -pedantic 


#### Rules ####

.PHONY : all clean

# Default rule is to create all binaries #
all: $(BINS)

# Rules for binaries #
$(BINS) : $(OBJS)
	$(CC) -o $@ $^ $(INCLUDES) $(LIB) $(OPTI) $(LINALG)  


# BiqBin code rules 
$(OBJ)/%.o : %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<


# Clean rule #
clean :
	rm $(BINS) $(OBJS)
