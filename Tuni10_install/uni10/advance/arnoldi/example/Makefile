EXU := cust1_arnoldi
INC := $(TUNI10_ROOT)/include
LIB := $(TUNI10_ROOT)/lib
CC := g++
FLAGS:= -m64 -std=c++11 -DCPU -DLAPACK
MPO_PATH:= ../../../example/nsy_mpo
OBJ := $(MPO_PATH)/mpo.o 

all: $(OBJ) $(EXU)

$(OBJ): %.o: %.cpp
	$(CC) -c -I$(INC) -L$(LIB) $(FLAGS) $< -o $@ -luni10

$(EXU): $(EXU).cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -luni10 -lblas -llapack -latlas

clean:
	rm -f *.o *.e asd $(EXU)
