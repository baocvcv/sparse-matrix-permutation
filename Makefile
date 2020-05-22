# Target
TARGET := OPT_PERM

# Options
CC := $(if $(shell command -v mpiicpc), mpiicpc, mpic++)
CC_FLAG := -std=c++11 -O3 -Wall -Wno-unused-result -fopenmp -pthread
CC_OBJ_FLAG := $(CCFLAG) -c
DBG_FLAG := -DDEBUG
LIB := -lcholmod \
      -lsuitesparseconfig -lccolamd -lcolamd \
      -lcamd -lamd -lmetis -lrt -lopenblas
LIB_DIR := #-L/home/yangming/Lib -lopenblas
INC_DIR := #-I./include -I/home/yangming/Include
RUN_CMD := $(if $(shell command -v srun), srun, mpiexec)

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
	$(CC) $(CC_OBJ_FLAG) $(LIB) $(LIB_DIR) $(INC_DIR) -o $@ $<

bin/$(TARGET_DBG): $(addsuffix $(DBG_SUFFIX).o, $(OBJ))
	$(CC) $(CC_FLAG) $(DBG_FLAG) $(LIB) $(LIB_DIR) $(INC_DIR) -o $@ $^

obj/%$(DBG_SUFFIX).o: src/%.cpp
	$(CC) $(CC_OBJ_FLAG) $(DBG_FLAG) $(LIB) $(LIB_DIR) $(INC_DIR) -o $@ $<

# Phony
.PHONY: all
all: bin/$(TARGET)

.PHONY: debug
debug: bin/$(TARGET_DBG)

.PHONY: clean
clean:
	rm -f obj/*.o bin/$(TARGET) bin/$(TARGET_DBG)

.PHONY: run
run: bin/$(TARGET)
	$(RUN_CMD) -n $(N) bin/$(TARGET) -f cases/1138_bus.mtx -n $(n) -s $(s)

.PHONY: run_d
run_d: bin/$(TARGET)
	$(RUN_CMD) -n $(N) bin/$(TARGET) -f cases/1138_bus.mtx -n 8 -s 1