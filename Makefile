#### BiqBin makefile ####

# container image name
IMAGE ?= serial-biqbin-general-bqp
# container image tag
TAG ?= 1.0.0
DOCKER_BUILD_PARAMS ?=
 
# Directories
OBJ = obj

# Compiler: other options (linux users) --> CC=icc make
CC ?= gcc

# NOTE: -framework Accelerate is for MAC, Linux users set to -lopenblas -lm (or use intel mkl)
LINALG 	 = -lopenblas -lm 
OPTI     = -O3 

# binary
BINS =  biqbin

# test command
TEST = ./test.sh \
	$(BINS) \
	test/test.dat \
	test/test.dat-expected_output \
	test/params

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

.PHONY : all clean test


# Default rule is to create all binaries #
all: $(BINS)

test: all
	$(TEST)

docker: 
	docker build $(DOCKER_BUILD_PARAMS) --progress=plain -t $(IMAGE):$(TAG)  . 

docker-clean: 
	docker rmi -f $(IMAGE):$(TAG) 

docker-test:
	docker run --rm $(IMAGE):$(TAG) $(TEST)

# Rules for binaries #
$(BINS) : $(OBJS)
	$(CC) -o $@ $^ $(INCLUDES) $(LIB) $(OPTI) $(LINALG)  


# BiqBin code rules 
$(OBJ)/%.o : %.c | $(OBJ)/
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<

$(OBJ)/:
	 mkdir -p $@


# Clean rule #
clean :
	rm -rf $(BINS) $(OBJ)
