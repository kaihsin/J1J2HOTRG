EXU := default_lanczos
INC := $(TUNI10_ROOT)/include
LIB := $(TUNI10_ROOT)/lib
CC := icpc
FLAGS:= -m64 -std=c++11 -DCPU -DLAPACK
MPO_PATH:= ../../../example/nsy_mpo
OBJ := $(MPO_PATH)/mpo.o 
#$(LIB)/libuni10.a

all: $(OBJ) $(EXU)

$(OBJ): %.o: %.cpp
	$(CC) -c -I$(INC) -L$(LIB) $(FLAGS) $< -o $@ -luni10

$(EXU): $(EXU).cpp
	$(CC) -I$(INC) -L$(LIB) $(FLAGS) -o $@ $< $(OBJ) -luni10 -mkl=parallel

clean:
	rm -f *.o *.e asd $(EXU)
