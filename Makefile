# MUSIC monorepo build.
#
# One copy of the tooling (tooling/), built against one dataset's config
# (analysis/<DATASET>/config/Constants.hpp). The dataset is chosen by the
# MUSIC_DATASET env var (set by the per-dataset nix dev shell), or on the
# command line: `make DATASET=87Rb`.
#
#   nix develop .#87Rb && make      # build everything for 87Rb
#   make DATASET=37Cl pipeline      # build one binary for 37Cl
#   make clean                      # remove this dataset's build artifacts
#
# Per-dataset outputs are self-contained: binaries in analysis/<DATASET>/bin,
# objects + libmusic.a in analysis/<DATASET>/build. tooling/src compiles once
# into a static libmusic.a; each main links it with --gc-sections. The GPU lib
# (tooling/gpu/libgpuaccel.so) is built once and dlopen'd at runtime via an
# absolute path injected as -DMUSIC_GPU_LIB.

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
ifeq ($(wildcard $(DATASET_DIR)/config/Constants.hpp),)
$(error No config at analysis/$(DATASET)/config/Constants.hpp for DATASET=$(DATASET))
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
GPU_LIB  := $(abspath $(GPU_DIR)/libgpuaccel.so)

CXXFLAGS  := -O3 -g -Wall -Wno-unused-variable -fPIC -std=c++17 \
             -march=native -mtune=native \
             -ffunction-sections -fdata-sections \
             $(ROOT_CFLAGS) -I$(INC_DIR) -I$(CFG_DIR) \
             '-DR__ADD_INCLUDE_PATH(...)=' \
             -DMUSIC_DATASET_NAME='"$(DATASET)"' \
             -DMUSIC_DATASET_DIR='"$(DATASET_DIR)"' \
             -DMUSIC_GIT_HASH='"$(GIT_HASH)"' \
             -DMUSIC_GPU_LIB='"$(GPU_LIB)"' \
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

# Binary "foo-bar" is built from tooling/mains/main_foo_bar.cpp (dashes->underscores).
BINS := pipeline calibrate-beam traces delta-e-scatter \
        strip-sum-scatter diag-timing diag-events diag-subfile-drift strip-scatter-overlay

# Si-detector calibration / stopping-power binaries. These need a per-dataset
# Si config (analysis/<DATASET>/config/SiCalibConstants.hpp); datasets without
# one (e.g. 87Rb) simply don't build them, instead of failing the build.
SI_BINS := si-fits si-calibration calc-stopping-power plot-stopping-power
ifneq ($(wildcard $(CFG_DIR)/SiCalibConstants.hpp),)
BINS += $(SI_BINS)
endif

BIN_PATHS := $(addprefix $(BIN_DIR)/,$(BINS))

.PHONY: all gpu clean
all: gpu $(BIN_PATHS)

gpu:
	$(MAKE) -C $(GPU_DIR)

# Every object explicitly depends on the dataset config headers. The -MMD auto-
# deps don't reliably re-trigger on the per-dataset config headers (they live on
# the -I path, outside the source tree), so editing Constants.hpp (or the Si
# config) wouldn't rebuild and you'd silently run a stale binary. Making them
# explicit prerequisites means a config change rebuilds everything -- no make
# clean. SiCalibConstants.hpp is wildcard'd so datasets without it still build.
# (Absolute source path so __FILE__ expands to one; ProjectRootOf relies on it.)
CFG_HEADERS := $(CFG_DIR)/Constants.hpp $(wildcard $(CFG_DIR)/SiCalibConstants.hpp)
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(CFG_HEADERS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $(abspath $<)

$(BUILD_DIR)/%.o: $(MAIN_DIR)/%.cpp $(CFG_HEADERS) | $(BUILD_DIR)
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
