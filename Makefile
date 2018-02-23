uni10v2ROOT := /home/kaywu/uni10_20

CC := icpc
CCFLAGS := -std=c++11 -O3 -mkl=parallel -DUNI_CPU -DUNI_LAPACK -DUNI_MKL
LDFLAGS := -L$(uni10v2ROOT)/lib -I$(uni10v2ROOT)/include
OBJ := $(uni10v2ROOT)/lib/libuni10.so.2

all : j1j2hotrg

j1j2hotrg : j1j2hotrg.o parser.o utils.o
	$(CC) $(CCFLAGS) $^ -o $@ $(OBJ) $(LDFLAGS)

j1j2hotrg.o : j1j2hotrg.cpp
	$(CC) -c $(CCFLAGS) $(LDFLAGS) $< -o $@ $(OBJ)

parser.o : Parser.cpp Parser.hpp
	$(CC) -c $(CCFLAGS) $(LDFLAGS) $< -o $@ $(OBJ)

utils.o : Utils.cpp Utils.hpp
	$(CC) -c $(CCFLAGS) $(LDFLAGS) $< -o $@ $(OBJ)

.phony : clean

clean :
	rm *.o j1j2hotrg






