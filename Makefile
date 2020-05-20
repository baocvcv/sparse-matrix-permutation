EXE = OPT_PERM

CC = g++ -std=c++11
OPTFLAG= -O3
COMPILE_FLAGS = -g -Wall -Wno-unused-result -fopenmp
lib = -pthread \
      -lcholmod \
      -lsuitesparseconfig -lccolamd -lcolamd \
      -lcamd -lamd -lmetis -lrt -lopenblas -lopenblas
MYINCDIR = ./include

# TODO: change this
LIBS = #-L/home/yangming/Lib -lopenblas
PATH_INCLUDE=#-I/home/yangming/Include
INCDIR =  -I$(MYINCDIR)
CFLAGS = $(OPTFLAG) $(COMPILE_FLAGS) $(INCDIR)
COMPILE = $(CC) $(CFLAGS) -c

SUB= main.o

SUBJECTS=$(SUB:%.o=obj/%.o)

ALL: bin/$(EXE)

bin/$(EXE): $(SUBJECTS)
	$(CC) $(CFLAGS) -o  bin/$(EXE) $(SUBJECTS) $(LIBS) $(lib)
	
obj/%.o: src/%.cpp
	$(COMPILE) -o  $@ $(PATH_INCLUDE) $< $(MYINCLDIR)

clean:
	\rm -f ./obj/*.o  bin/OPT_PERM
