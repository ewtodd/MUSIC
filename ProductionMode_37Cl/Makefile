# Build standalone executables for each entry point.
#   make            # build gpu/libgpuaccel.so and all binaries
#   make pipeline   # build just one binary
#   make clean      # remove build artifacts
#
# Binaries land in ./bin/, object/dep files in ./build/. The GPU library
# (gpu/libgpuaccel.so) is built first and dlopen'd at runtime by
# GpuAccel::Init() — its path comes from LD_LIBRARY_PATH (the nix shellHook
# adds ./gpu).
#
# Nix shell users: -I and -L paths for ROOT and analysis-utilities are
# injected by the gcc wrapper via NIX_CFLAGS_COMPILE / NIX_LDFLAGS, so we
# only need to name libraries here.

CXX         := g++
ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS   := $(shell root-config --glibs)

# The nix gcc-wrapper strips -march=native unless this is cleared.
unexport NIX_ENFORCE_NO_NATIVE
export NIX_ENFORCE_NO_NATIVE :=

INC_DIR   := include
SRC_DIR   := src
MAIN_DIR  := mains
BIN_DIR   := bin
BUILD_DIR := build

CXXFLAGS  := -O3 -g -Wall -Wno-unused-variable -fPIC -std=c++17 \
             -march=native -mtune=native \
             -ffunction-sections -fdata-sections \
             $(ROOT_CFLAGS) -I$(INC_DIR) \
             '-DR__ADD_INCLUDE_PATH(...)=' \
             -MMD -MP
LDFLAGS   := -Wl,--gc-sections $(ROOT_LIBS) -lSpectrum \
             -l:libanalysis-utils.so -ldl -lpthread

# Every binary links the full set of library objects; -Wl,--gc-sections drops
# whatever a given entry point does not use.
SRC_OBJS  := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(wildcard $(SRC_DIR)/*.cpp))

BINS := $(BIN_DIR)/pipeline \
        $(BIN_DIR)/calibrate-beam \
        $(BIN_DIR)/traces \
        $(BIN_DIR)/delta-e-scatter \
        $(BIN_DIR)/diag-timing \
        $(BIN_DIR)/diag-events \
        $(BIN_DIR)/diag-subfile-drift

.PHONY: all gpu clean
all: gpu $(BINS)

gpu:
	$(MAKE) -C gpu

# Absolute source path so __FILE__ expands to one too — Paths::ProjectRootOf
# walks two DirNames up (src/ or mains/ -> project root) to find plots/,
# root_files/, and the log file.
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $(abspath $<)

$(BUILD_DIR)/%.o: $(MAIN_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $(abspath $<)

$(BUILD_DIR) $(BIN_DIR):
	mkdir -p $@

$(BIN_DIR)/pipeline:           $(BUILD_DIR)/main_pipeline.o          $(SRC_OBJS) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)
$(BIN_DIR)/calibrate-beam:     $(BUILD_DIR)/main_calibrate_beam.o    $(SRC_OBJS) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)
$(BIN_DIR)/traces:             $(BUILD_DIR)/main_traces.o            $(SRC_OBJS) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)
$(BIN_DIR)/delta-e-scatter:    $(BUILD_DIR)/main_delta_e_scatter.o   $(SRC_OBJS) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)
$(BIN_DIR)/diag-timing:        $(BUILD_DIR)/main_diag_timing.o       $(SRC_OBJS) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)
$(BIN_DIR)/diag-events:        $(BUILD_DIR)/main_diag_events.o       $(SRC_OBJS) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)
$(BIN_DIR)/diag-subfile-drift: $(BUILD_DIR)/main_diag_subfile_drift.o $(SRC_OBJS) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

# Convenience shortcuts so `make pipeline` etc. still work.
.PHONY: pipeline calibrate-beam traces delta-e-scatter \
        diag-timing diag-events diag-subfile-drift
pipeline:           $(BIN_DIR)/pipeline
calibrate-beam:     $(BIN_DIR)/calibrate-beam
traces:             $(BIN_DIR)/traces
delta-e-scatter:    $(BIN_DIR)/delta-e-scatter
diag-timing:        $(BIN_DIR)/diag-timing
diag-events:        $(BIN_DIR)/diag-events
diag-subfile-drift: $(BIN_DIR)/diag-subfile-drift

-include $(SRC_OBJS:.o=.d)
-include $(BUILD_DIR)/main_pipeline.d $(BUILD_DIR)/main_calibrate_beam.d \
         $(BUILD_DIR)/main_traces.d $(BUILD_DIR)/main_delta_e_scatter.d \
         $(BUILD_DIR)/main_diag_timing.d $(BUILD_DIR)/main_diag_events.d \
         $(BUILD_DIR)/main_diag_subfile_drift.d

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)
	$(MAKE) -C gpu clean
