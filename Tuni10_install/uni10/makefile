PROJECT := uni10
CONFIG_FILE := make.inc

##########################################

## Explicitly check for the config file, otherwise make -k will proceed anyway.
ifeq ($(wildcard $(CONFIG_FILE)),)
$(error $(CONFIG_FILE) not found.)
endif
include $(CONFIG_FILE)

## Identify the version name
VERSION := $(LNPACKAGE)_$(CALARCH)
UNI10_VERSION := $(addprefix $(PROJECT)_, $(VERSION))

## The direction to save the objects.

UNI10_INSTALL_PATH := uni10/
ifneq ($(UNI10_INSTALL_PREFIX),)
	UNI10_INSTALL_PATH := $(UNI10_INSTALL_PREFIX)/
endif

SRC_DIRS := $(shell find * -type d -exec bash -c "find {} -maxdepth 1 \
	\( -name '*.cpp' -o -name '*.cu' \) | grep -q ." \; -print)

INC_DIRS := $(shell find * -type d -exec bash -c "find {} -maxdepth 1 \
	\( -name '*.h' -o -name '*.hpp' \) | grep -q ." \; -print)

## Names of the static and dynamic shared libraries.

UNI10_LIB_PATH := $(addprefix $(UNI10_INSTALL_PATH),lib)
UNI10_INC_PATH := $(addprefix $(UNI10_INSTALL_PATH),include)
STATIC_NAME    := $(addprefix $(UNI10_LIB_PATH)/,lib$(UNI10_VERSION).a)
DYNAMIC_NAME   := $(addprefix $(UNI10_LIB_PATH)/,lib$(UNI10_VERSION).so)
UNI10_GTEST_PATH := $(addprefix $(UNI10_INSTALL_PATH),gtest)

##############################
# Find all source files
##############################

# CXX_SRCS are the source files excluding the test ones.
#
CXX_SRCS :=
#Find the auxiliary files
CXX_SRCS += $(shell find src/$(PROJECT) -maxdepth 1 -name "*.cpp")
CXX_SRCS += $(shell find src/$(PROJECT)/uni10_env_info/uni10_$(VERSION) -name "*$(VERSION).cpp")
# Find UNI10_API FILES
CXX_SRCS += $(shell find src/$(PROJECT)/uni10_api -name "*.cpp")
# Find the others sources of the specified versions.
CXX_SRCS += $(shell find src/$(PROJECT)/uni10_$(VERSION) -name "*.cpp")
# EXAMPLE_SRCS are the source files for the example binaries
EXAMPLE_SRCS := $(shell find example -name "*.cpp")
# TESTINGS_SRCS are the source files for the unit tests binaries
TESTINGS_SRCS := $(shell find testing -name "*.cpp")
# API_TESTINGS_SRCS are the source files for the unit tests binaries
GTEST_SRCS := $(shell find gtest -name "*.cpp")

# CU_SRCS
CU_SRCS :=
ifeq ($(CALARCH), gpu)
	#Find the auxiliary files
	CU_SRCS += $(shell find src/$(PROJECT)/uni10_env_info/uni10_$(VERSION) -name "*.cu")
	# Find the others sources of the specified versions.
	CU_SRCS += $(shell find src/$(PROJECT)/uni10_$(VERSION) -name "*.cu")
endif

ifeq ($(GTEST), 1)
	CXX_SRCS += $(GTEST_SRCS)
endif

##############################
# Derive generated files
##############################

# CXX objects
CXX_OBJS := $(addprefix $(UNI10_INSTALL_PATH), ${CXX_SRCS:.cpp=.o})

# CU objects
CU_OBJS := 
ifeq ($(CALARCH), gpu)
	CU_OBJS = $(addprefix $(UNI10_INSTALL_PATH), ${CU_SRCS:.cu=.o})
endif

#EXAMPLE_OBJS := $(addprefix $(UNI10_INSTALL_PATH)/, ${EXAMPLE_SRCS:.cpp=.o})
#EXAMPLE_BINS := ${EXAMPLE_OBJS:.o=.bin}
#
#TESTING_OBJS := $(addprefix $(UNI10_INSTALL_PATH)/, ${EXAMPLE_SRCS:.cpp=.o})
#TESTING_BINS := ${EXAMPLE_OBJS:.o=.bin}

GTEST_BIN := $(addprefix $(UNI10_INSTALL_PATH), ${GTEST_SRCS:.cpp=.bin})

##############################
# Derive compiler warning dump locations
##############################

WARNS_EXT := warnings.txt

CXX_WARNS := $(addprefix $(UNI10_INSTALL_PATH), ${CXX_SRCS:.cpp=.o.$(WARNS_EXT)})

#EXAMPLE_WARNS := $(addprefix $(UNI10_INSTALL_PATH)/, ${EXAMPLE_SRCS:.cpp=.o.$(WARNS_EXT)})
#TESTING_WARNS := $(addprefix $(UNI10_INSTALL_PATH)/, ${EXAMPLE_SRCS:.cpp=.o.$(WARNS_EXT)})

##############################
# Set build directories
##############################

INCLUDE_DIRS += $(UNI10_INC_PATH) ./src ./include

ALL_BUILD_DIRS := $(sort $(UNI10_INSTALL_PREFIX) $(UNI10_LIB_PATH) \
	$(addprefix $(UNI10_INSTALL_PATH), $(SRC_DIRS)))

##############################
# Configure build
##############################
#
# Determine platform
UNAME := $(shell uname -s)
ifeq ($(UNAME), Linux)
	LINUX := 1
else ifeq ($(UNAME), Darwin)
	OSX := 1
endif

# Linux
ifeq ($(LINUX), 1)
	COMMON_FLAGS += -DLINUX
	ifeq ($(CXX), g++)
		CXX = /usr/bin/g++
		GCCVERSION := $(shell $(CXX) -dumpversion | cut -f1,2 -d.)
		# older versions of gcc are too dumb to build boost with -Wuninitalized
		ifeq ($(shell echo | awk '{exit $(GCCVERSION) < 4.6;}'), 1)
			WARNINGS  += -Wno-uninitialized
		endif
	endif
	ifeq ($(CXX), icpc)
		CXX = /opt/intel/bin/icpc
	endif
endif

INCLUDE_DIRS += $(BLAS_INCLUDE)
LIBRARY_DIRS += $(BLAS_LIB)
ifeq ($(CALARCH), gpu)
	INCLUDE_DIRS += $(CUDA_INCLUDE)
	LIBRARY_DIRS += $(CUDA_LIB)
endif
#
LIBRARY_DIRS += $(UNI10_LIB_PATH)
#
# Automatic dependency generation (nvcc is handled separately)
CXXFLAGS += -MMD -MP

# Complete build flags.
COMMON_FLAGS += $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))
CXXFLAGS     += $(COMMON_FLAGS) $(WARNINGS)

LDFLAGS += $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) \
		$(foreach library,$(LDLIBRARIES),-l$(library))

GTESTLDLIBS = gtest gtest_main pthread
GTESTLDFLAGS += $(foreach gtestlib,$(GTESTLDLIBS),-l$(gtestlib))

GPU_TARGET ?= Fermi Kepler

NV_SM    :=
NV_COMP  :=

ifeq ($(CALARCH), gpu)

	# NVCC options for the different cards
	# First, add smXX for architecture names
	NVCCFLAGS += $(UNI10CXXFLAGS)	

	ifneq ($(findstring Fermi, $(GPU_TARGET)),)
		GPU_TARGET += sm20
 	endif
		ifneq ($(findstring Kepler, $(GPU_TARGET)),)
	GPU_TARGET += sm30 sm35
	endif
		ifneq ($(findstring Maxwell, $(GPU_TARGET)),)
	GPU_TARGET += sm50 sm52
	endif
		ifneq ($(findstring Pascal, $(GPU_TARGET)),)
	GPU_TARGET += sm60 sm61
	endif

	ifneq ($(findstring sm10, $(GPU_TARGET)),)
	    # sm10 is no longer supported by CUDA 6.x nvcc
	    #NV_SM    += -gencode arch=compute_10,code=sm_10
	    #NV_COMP  := -gencode arch=compute_10,code=compute_10
	    $(warning CUDA arch 1.x is no longer supported by CUDA >= 6.x and MAGMA >= 2.0)
	endif
	ifneq ($(findstring sm13, $(GPU_TARGET)),)
	    #NV_SM    += -gencode arch=compute_13,code=sm_13
	    #NV_COMP  := -gencode arch=compute_13,code=compute_13
	    $(warning CUDA arch 1.x is no longer supported by CUDA >= 6.x and MAGMA >= 2.0)
	endif
	ifneq ($(findstring sm20, $(GPU_TARGET)),)
	    NV_SM    += -gencode arch=compute_20,code=sm_20
	    NV_COMP  := -gencode arch=compute_20,code=compute_20
	endif
	ifneq ($(findstring sm30, $(GPU_TARGET)),)
	    NV_SM    += -gencode arch=compute_30,code=sm_30
	    NV_COMP  := -gencode arch=compute_30,code=compute_30
	endif
	ifneq ($(findstring sm35, $(GPU_TARGET)),)
	    NV_SM    += -gencode arch=compute_35,code=sm_35
	    NV_COMP  := -gencode arch=compute_35,code=compute_35
	endif
	ifneq ($(findstring sm50, $(GPU_TARGET)),)
	    NV_SM    += -gencode arch=compute_50,code=sm_50
	    NV_COMP  := -gencode arch=compute_50,code=compute_50
	endif
	ifneq ($(findstring sm52, $(GPU_TARGET)),)
	    NV_SM    += -gencode arch=compute_52,code=sm_52
	    NV_COMP  := -gencode arch=compute_52,code=compute_52
	endif
	ifneq ($(findstring sm60, $(GPU_TARGET)),)
	    NV_SM    += -gencode arch=compute_60,code=sm_60
	    NV_COMP  := -gencode arch=compute_60,code=compute_60
	endif
	ifneq ($(findstring sm61, $(GPU_TARGET)),)
	    NV_SM    += -gencode arch=compute_61,code=sm_61
	    NV_COMP  := -gencode arch=compute_61,code=compute_61
	endif
	ifeq ($(NV_COMP),)
	    $(error GPU_TARGET, currently $(GPU_TARGET), must contain one or more of Fermi, Kepler, Maxwell, Pascal, or sm{20,30,35,50,52,60,61}. Please edit your make.inc file)
	endif
	NVCCFLAGS += $(NV_SM) $(NV_COMP)

endif

# 'superclean' target recursively* deletes all files ending with an extension
# in $(SUPERCLEAN_EXTS) below.  This may be useful if you've built older
# versions of Caffe that do not place all generated files in a location known
# to the 'clean' target.
#
# 'supercleanlist' will list the files to be deleted by make superclean.
#
# * Recursive with the exception that symbolic links are never followed, per the
# default behavior of 'find'.
SUPERCLEAN_EXTS := .so .a .o

##############################
# Define build targets
##############################

TARGET := lib include

ifeq ($(GTEST), 1)
	TARGET += gtest
endif

all: $(TARGET)

lib: $(STATIC_NAME) $(DYNAMIC_NAME)

gtest: $(GTEST_BIN)

include: $(UNI10_INSTALL_PATH)
	@ echo CP -r include $<
	$(Q)cp -r include $<
	@ echo CP -r example $<
	$(Q)cp -r example $<
	@ echo CP -r advance $<
	$(Q)cp -r advance $<
	@ echo CP -r testing $<
	$(Q)cp -r testing $<

$(ALL_BUILD_DIRS):
	@ mkdir -p $@

$(DYNAMIC_NAME): $(CXX_OBJS) $(CU_OBJS) | $(UNI10_LIB_PATH)
	@ echo LD -o $@
	$(Q)$(CXX) -shared -o $@ $(CXX_OBJS) $(CU_OBJS) $(LINKFLAGS) $(LDFLAGS) $(DYNAMIC_FLAGS)

$(STATIC_NAME): $(CXX_OBJS) $(CU_OBJS) | $(UNI10_LIB_PATH)
	@ echo AR -o $@
	$(Q)ar rcs $@ $(CXX_OBJS) $(CU_OBJS)

$(UNI10_INSTALL_PATH)%.o: %.cpp | $(ALL_BUILD_DIRS)
	@ echo CXX $<
	$(Q)$(CXX) $< $(CXXFLAGS) -c -o $@ 2> $@.$(WARNS_EXT) \
		|| (cat $@.$(WARNS_EXT); exit 1)
	@ cat $@.$(WARNS_EXT)

$(UNI10_INSTALL_PATH)%.o: %.cu | $(ALL_BUILD_DIRS)
	@ echo nvcc $<
	$(Q)nvcc $< $(NVCCFLAGS) $(COMMON_FLAGS) -c -o $@ 2> $@.$(WARNS_EXT) \
		|| (cat $@.$(WARNS_EXT); exit 1)
	@ cat $@.$(WARNS_EXT)

$(GTEST_BIN): $(GTEST_SRCS)
	@ echo CXX $< -o $@
	$(Q)$(CXX) $(CXXFLAGS) $< -o $@ $(STATIC_NAME) $(LDFLAGS) $(GTESTLDFLAGS)

clean:
	@- $(RM) -rf $(ALL_BUILD_DIRS)
	@- $(RM) -rf docs/latex docs/html docs/man

clean_header:
	@- $(RM) -rf $(UNI10_INC_PATH)

clean_doxygen:
	@- $(RM) -rf docs/latex docs/html docs/man

print-%  : ; @echo $* = $($*)
