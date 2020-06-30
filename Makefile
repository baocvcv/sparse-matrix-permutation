# Target
TARGET := OPT_PERM

# Options
CC := $(if $(shell command -v mpiicpc), mpiicpc, mpic++)
CC_FLAG := -std=c++11 -O3 -Wall -Wno-unused-result -fopenmp -pthread
CC_OBJ_FLAG := $(CC_FLAG) -c
DBG_FLAG := -DDEBUG
LIB := -lcholmod \
      -lsuitesparseconfig -lccolamd -lcolamd \
      -lcamd -lamd -lmetis -lrt -lopenblas
#LIB_DIR := -L./libs/OpenBLAS/lib -L./libs/SuiteSparse/lib
#INC_DIR := -I./src -I./libs/OpenBLAS/include -I./libs/SuiteSparse/include
#LIB_DIR_ENV := ./libs/OpenBLAS/lib:./libs/SuiteSparse/lib
#INC_DIR_ENV := ./libs/OpenBLAS/include:./libs/SuiteSparse/include
RUN_CMD := $(if $(shell command -v srun), srun, mpiexec)
LD_LIBRARY_PATH := /usr/local/lib:./libs/SuiteSparse/lib:./libs/OpenBLAS/lib

# Files
DBG_SUFFIX := _dbg
TARGET_DBG := $(addsuffix $(DBG_SUFFIX), $(TARGET))
SRC := $(basename $(wildcard src/*.c*))
OBJ := $(addprefix obj/, $(notdir $(SRC)))

# Rules
default: all

bin/$(TARGET): $(addsuffix .o, $(OBJ))
	$(CC) $(CC_FLAG) $(LIB) $(LIB_DIR) $(INC_DIR) -o $@ $^
	
obj/%.o: src/%.cpp
	CPLUS_INCLUDE_PATH=$(INC_DIR_ENV) $(CC) $(CC_OBJ_FLAG) $(LIB) $(LIB_DIR) $(INC_DIR) -o $@ $<

bin/$(TARGET_DBG): $(addsuffix $(DBG_SUFFIX).o, $(OBJ))
	$(CC) $(CC_FLAG) $(DBG_FLAG) $(LIB) $(LIB_DIR) $(INC_DIR) -o $@ $^

obj/%$(DBG_SUFFIX).o: src/%.cpp
	CPLUS_INCLUDE_PATH=$(INC_DIR_ENV) $(CC) $(CC_OBJ_FLAG) $(DBG_FLAG) $(LIB) $(LIB_DIR) $(INC_DIR) -o $@ $<

# Phony
.PHONY: all
all: mkdirs bin/$(TARGET)

.PHONY: debug
debug: mkdirs bin/$(TARGET_DBG)

.PHONY: clean
clean:
	rm -f obj/*.o bin/$(TARGET) bin/$(TARGET_DBG)

.PHONY: run
run: mkdirs bin/$(TARGET)
	$(RUN_CMD) -n $(N) bin/$(TARGET) -f cases/1138_bus.mtx -n $(n) -s $(s)

.PHONY: run_d
run_d: mkdirs bin/$(TARGET)
	LD_LIBRARY_PATH=$(LD_LIBRARY_PATH) $(RUN_CMD) -n 2 bin/$(TARGET) -f cases/1138_bus.mtx -n 4 -s 5

.PHONY: mkdirs
mkdirs:
	@mkdir -p bin
	@mkdir -p out
	@mkdir -p obj
