CC := g++
UNIFLAGS := -DUNI_CPU -DUNI_LAPACK
CXXFLAGS := -O3 -std=c++11
LDFLAGS := -lgtest -lgtest_main -lpthread -luni10_lapack_cpu
INC := $(UNI10_ROOT_CPU)/include
LIB := $(UNI10_ROOT_CPU)/lib

all:
	$(CC) -I$(INC) -L$(LIB) $(CXXFLAGS) $(UNIFLAGS) runUni10Test.cpp -o runUni10Test $(LDFLAGS)

clean:
	rm -f runUni10Test
