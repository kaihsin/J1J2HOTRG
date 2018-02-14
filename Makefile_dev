uni10v2ROOT := /home/kaywu/uni10_20

CC := icpc
CCFLAGS := -std=c++11 -g -DUNI_CPU -DUNI_LAPACK -DUNI_MKL
LDFLAGS := -I$(uni10v2ROOT)/include -L$(uni10v2ROOT)/lib -luni10
OBJ := $(uni10v2ROOT)/lib/libuni10.so.2

all : j1j2hotrg

j1j2hotrg : j1j2hotrg.o parser.o utils.o
	$(CC) $(CCFLAGS) $(LDFLAGS) $^ -o $@

j1j2hotrg.o : j1j2hotrg.cpp
	$(CC) -c $(CCFLAGS) $(LDFLAGS) $< -o $@

parser.o : Parser.cpp Parser.hpp
	$(CC) -c $(CCFLAGS) $(LDFLAGS) $< -o $@

utils.o : Utils.cpp Utils.hpp
	$(CC) -c $(CCFLAGS) $(LDFLAGS) $< -o $@

.phony : clean

clean :
	rm *.o j1j2hotrg






