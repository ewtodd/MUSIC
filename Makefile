# MUSIC monorepo build.
#
# The flake invokes the Makefile with correct args; just use nix build with the 
# chosen dataset
#
# Per-dataset outputs are self-contained: binaries in analysis/<DATASET>/bin,
# objects + libmusic.a in analysis/<DATASET>/build. tooling/src compiles once
# into a static libmusic.a; each main links it with --gc-sections. The GPU lib
# (tooling/gpu/libgpuaccel.so) is optional: when built, it is dlopen'd at
# runtime via an absolute path injected as -DMUSIC_GPU_LIB, and GpuAccel falls
# back to the CPU sort if the lib (or a GPU) is missing.
#
# Dataset config: tooling/include/Constants.hpp defines the DatasetConfig class.
# Each analysis provides analysis/<DATASET>/config/Constants.cpp which creates
# and initializes the global instance. Editing the .cpp rebuilds all objects.

CXX         := g++
ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS   := $(shell root-config --glibs)

# The nix gcc-wrapper strips -march=native unless this is cleared.
unexport NIX_ENFORCE_NO_NATIVE
export NIX_ENFORCE_NO_NATIVE :=

# ---- dataset selection ----
DATASET ?= $(MUSIC_DATASET)
ifeq ($(strip $(DATASET)),)
$(error MUSIC_DATASET not set. Enter a dataset dev shell (e.g. `nix develop .#87Rb`) or pass DATASET=<iso>)
endif
DATASET_DIR := $(abspath analysis/$(DATASET))
ifeq ($(wildcard $(DATASET_DIR)/config/Constants.cpp),)
$(error No config at analysis/$(DATASET)/config/Constants.cpp for DATASET=$(DATASET))
endif

# ---- tree ----
TOOLING   := tooling
INC_DIR   := $(TOOLING)/include
SRC_DIR   := $(TOOLING)/src
MAIN_DIR  := $(TOOLING)/mains
GPU_DIR   := $(TOOLING)/gpu
CFG_DIR   := $(DATASET_DIR)/config
BIN_DIR   := $(DATASET_DIR)/bin
BUILD_DIR := $(DATASET_DIR)/build

GIT_HASH := $(shell git rev-parse --short HEAD 2>/dev/null || echo unknown)
GPU_LIB  ?= $(abspath $(GPU_DIR)/libgpuaccel.so)
# Runtime paths that can be overridden by nix build (default: same as build paths)
DATASET_DIR_OUT ?= $(DATASET_DIR)
GPU_LIB_OUT     ?= $(GPU_LIB)

CXXFLAGS  := -O3 -g -Wall -Wno-unused-variable -fPIC -std=c++17 \
             -march=native -mtune=native \
             -ffunction-sections -fdata-sections \
             $(ROOT_CFLAGS) -I$(INC_DIR) \
             '-DR__ADD_INCLUDE_PATH(...)=' \
             -DMUSIC_DATASET_NAME='"$(DATASET)"' \
             -DMUSIC_DATASET_DIR='"$(DATASET_DIR_OUT)"' \
             -DMUSIC_GIT_HASH='"$(GIT_HASH)"' \
             -DMUSIC_GPU_LIB='"$(GPU_LIB_OUT)"' \
             -MMD -MP
LDFLAGS   := -Wl,--gc-sections $(ROOT_LIBS) -lSpectrum -lMinuit \
             -l:libanalysis-utils.so -ldl -lpthread

LIB      := $(BUILD_DIR)/libmusic.a
SRC_OBJS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(wildcard $(SRC_DIR)/*.cpp))

# Si-detector calibration / stopping-power library objects. They #include the
# per-dataset SiCalibConstants.hpp, so they only compile for a dataset that
# provides one; otherwise drop them from the library (the matching binaries are
# gated the same way below). This is what lets `make DATASET=87Rb` build.
SI_SRC_OBJS := $(addprefix $(BUILD_DIR)/,SiFits.o SiCalibration.o \
               CalcStoppingPower.o PlotStoppingPower.o)
ifeq ($(wildcard $(CFG_DIR)/SiCalibConstants.hpp),)
SRC_OBJS := $(filter-out $(SI_SRC_OBJS),$(SRC_OBJS))
endif

# Per-dataset config instance: compiled from analysis/<DATASET>/config/Constants.cpp
DATASET_CFG_SRC := $(CFG_DIR)/Constants.cpp
DATASET_CFG_OBJ := $(BUILD_DIR)/dataset_constants.o
SRC_OBJS += $(DATASET_CFG_OBJ)

# Binary "foo-bar" is built from tooling/mains/main_foo_bar.cpp (dashes->underscores).
BINS := pipeline calibrate-beam strip-sum-scatter split-sol preprocess-sol

BIN_PATHS := $(addprefix $(BIN_DIR)/,$(BINS))

.PHONY: all gpu clean
# The GPU lib needs the CUDA toolchain (nvcc), which only the CUDA shell/build
# provides; without it, skip the GPU lib and rely on the CPU fallback. Override
# with `make WITH_GPU=1` (or =0) to force either way.
WITH_GPU ?= $(shell command -v nvcc >/dev/null 2>&1 && echo 1 || echo 0)
ifeq ($(WITH_GPU),1)
all: gpu $(BIN_PATHS)
else
all: $(BIN_PATHS)
endif

gpu:
	@if ! command -v nvcc >/dev/null 2>&1; then \
	  echo "nvcc not found; GPU lib cannot be built (use WITH_GPU=0 or a CUDA shell)"; \
	  exit 1; \
	fi
	$(MAKE) -C $(GPU_DIR)

# Every object explicitly depends on the dataset config source. Editing
# Constants.cpp rebuilds everything -- no make clean needed.
CFG_SRC := $(DATASET_CFG_SRC)
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(CFG_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $(abspath $<)

$(DATASET_CFG_OBJ): $(DATASET_CFG_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $(abspath $<)

$(BUILD_DIR)/%.o: $(MAIN_DIR)/%.cpp $(CFG_SRC) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $(abspath $<)

$(LIB): $(SRC_OBJS)
	ar rcs $@ $^

$(BUILD_DIR) $(BIN_DIR):
	mkdir -p $@

# Each binary: its main object + libmusic.a (--gc-sections drops unused code).
# Rules are generated from BINS: binary "foo-bar" <- build/main_foo_bar.o.
define LINK_RULE
$(BIN_DIR)/$(1): $(BUILD_DIR)/main_$(subst -,_,$(1)).o $(LIB) | $(BIN_DIR)
	$$(CXX) -o $$@ $$< $(LIB) $(LDFLAGS)
.PHONY: $(1)
$(1): $(BIN_DIR)/$(1)
endef
$(foreach b,$(BINS),$(eval $(call LINK_RULE,$(b))))

-include $(SRC_OBJS:.o=.d)
-include $(wildcard $(BUILD_DIR)/main_*.d)

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)
	$(MAKE) -C $(GPU_DIR) clean
