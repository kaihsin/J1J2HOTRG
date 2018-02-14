### BLAS configuration ###
ifneq ($(shell readelf --dynamic $(TUNI10_ROOT)/lib/libuni10.so | grep atlas),)
	CXX := g++
	BLASFLAGS := -lblas -llapack -latlas
else ifneq ($(shell readelf --dynamic $(TUNI10_ROOT)/lib/libuni10.so | grep mkl),)
	CXX := icpc
	BLASFLAGS := -lmkl_rt
endif

### Set directories ###
SRC_DIR := ./src
BUILD_DIR := ./build
INC_DIRS := $(TUNI10_ROOT)/include ./include
LIB_DIRS := $(TUNI10_ROOT)/lib

### Compiler and Linker configuration ###
LIB := $(foreach dir,$(LIB_DIRS),-L$(dir)) $(BLASFLAGS) -luni10
INC := $(foreach dir,$(INC_DIRS),-I$(dir))
CXXFLAGS := $(INC) -m64 -std=c++11 -O0 -pipe
LDFLAGS := $(LIB) $(INC) -O0

### All testing tools ###
TOOLS := testing_qr_d testing_qr_z \
	testing_rq_d testing_rq_z \
	testing_ql_d testing_ql_z \
	testing_lq_d testing_lq_z \
	testing_qdr_d testing_qdr_z \
	testing_ldq_d testing_ldq_z \
	testing_qdr_cpivot_d testing_qdr_cpivot_z \
	testing_svd_d testing_svd_z \
	testing_inverse_d testing_inverse_z

### Targets ###
all: $(TOOLS)

$(TOOLS): %: $(BUILD_DIR)/%.o $(BUILD_DIR)/testing_timer.o
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(LDFLAGS) $^ -o $(BUILD_DIR)/$@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -rf $(BUILD_DIR)
