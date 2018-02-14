# The name of binary object.
EXU := run_5pess

# Search the source files in under paths.
vpath %.cpp ../../../operator ../../../hamiltonian ../../../../common 5pess_tools

# The file names of source codes.
src := operator.cpp hamiltonian.cpp 5pess_tools.cpp 5pess_kagome.cpp common_tools.cpp


CONFIG_FILE := make.inc
ifeq ($(wildcard $(CONFIG_FILE)),)
$(error $(CONFIG_FILE) not found.)
endif
include $(CONFIG_FILE)

UNI10_INC := $(UNI10_ROOT)/include
UNI10_LIB := $(UNI10_ROOT)/lib
CXXFLAGS += $(UNI10CXXFLAGS)

UNI10LDLIBRARY := $(addprefix -luni10_, $(UNI10_VERSION))

OBJ_DIR := objects_$(UNI10_VERSION)
OBJECTS := $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(src))

$(EXU): $(EXU).cpp | $(OBJECTS)
	$(CXX) -I$(UNI10_INC) -L$(UNI10_LIB) $(CXXFLAGS) $< -o $@ $(OBJECTS) $(UNI10LDLIBRARY)

$(OBJECTS): | $(OBJ_DIR)

$(OBJ_DIR):
	@mkdir -p $@

$(OBJ_DIR)/%.o: %.cpp
	$(CXX) -I$(UNI10_INC) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(EXU) $(OBJ_DIR)

print-%  : ; @echo $* = $($*)
