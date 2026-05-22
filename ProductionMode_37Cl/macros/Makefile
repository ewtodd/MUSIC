# Build standalone executables for each entry point.
# Usage:
#   make            # build libgpuaccel.so and all binaries
#   make pipeline   # build just one binary
#   make clean      # remove build artifacts
#
# Binaries land in ./bin/, object/dep files in ./build/. The GPU library
# (../gpu/libgpuaccel.so) is built first and dlopen'd at runtime by
# GpuAccel::Init() — link paths come from LD_LIBRARY_PATH (the nix
# shellHook adds ../gpu).
#
# Nix shell users: -I and -L paths for ROOT and analysis-utilities are
# injected by the gcc wrapper via NIX_CFLAGS_COMPILE / NIX_LDFLAGS, so
# we only need to name libraries here.

CXX       := g++
ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS   := $(shell root-config --glibs)

# The nix gcc-wrapper strips -march=native unless this is cleared.
# Make exports the empty value to every recipe sub-shell.
unexport NIX_ENFORCE_NO_NATIVE
export NIX_ENFORCE_NO_NATIVE :=

BIN_DIR   := bin
BUILD_DIR := build

CXXFLAGS  := -O3 -g -Wall -Wno-unused-variable -fPIC -std=c++17 \
             -march=native -mtune=native \
             -ffunction-sections -fdata-sections \
             $(ROOT_CFLAGS) -I. -I../diagnostics \
             '-DR__ADD_INCLUDE_PATH(...)=' \
             -MMD -MP
LDFLAGS   := -Wl,--gc-sections $(ROOT_LIBS) \
             -l:libanalysis-utils.so -ldl -lpthread

# Shared library objects linked into every binary.
LIB_OBJS  := $(BUILD_DIR)/Timing.o $(BUILD_DIR)/EventBuilder.o \
             $(BUILD_DIR)/GpuAccel.o

# Entry-point translation units (each contributes one binary).
PIPELINE_OBJ := $(BUILD_DIR)/Pipeline.o
CATHODE_OBJ  := $(BUILD_DIR)/CathodeAnodeRatio.o
REGION_OBJ   := $(BUILD_DIR)/RegionTraces.o
CALIB_OBJ    := $(BUILD_DIR)/CalibrateBeam.o
TRACES_OBJ   := $(BUILD_DIR)/Traces.o

DIAG_TIM_OBJ  := $(BUILD_DIR)/DiagnoseTiming.o
DIAG_EVT_OBJ  := $(BUILD_DIR)/DiagnoseEvents.o
DIAG_GAIN_OBJ := $(BUILD_DIR)/DiagnoseGainDrift.o
DIAG_SUB_OBJ  := $(BUILD_DIR)/DiagnoseSubfileDrift.o
DIAG_RATE_OBJ := $(BUILD_DIR)/DiagnoseSubfileRate.o

MAIN_OBJS := $(BUILD_DIR)/main_pipeline.o \
             $(BUILD_DIR)/main_cathode_anode.o \
             $(BUILD_DIR)/main_region_traces.o \
             $(BUILD_DIR)/main_calibrate_beam.o \
             $(BUILD_DIR)/main_traces.o \
             $(BUILD_DIR)/main_diag_timing.o \
             $(BUILD_DIR)/main_diag_events.o \
             $(BUILD_DIR)/main_diag_gain_drift.o \
             $(BUILD_DIR)/main_diag_subfile_drift.o \
             $(BUILD_DIR)/main_diag_subfile_rate.o

BINS := $(BIN_DIR)/pipeline $(BIN_DIR)/cathode-anode $(BIN_DIR)/region-traces \
        $(BIN_DIR)/calibrate-beam $(BIN_DIR)/traces \
        $(BIN_DIR)/diag-timing $(BIN_DIR)/diag-events \
        $(BIN_DIR)/diag-gain-drift $(BIN_DIR)/diag-subfile-drift \
        $(BIN_DIR)/diag-subfile-rate

.PHONY: all gpu clean clean-aclic

all: gpu $(BINS)

gpu:
	$(MAKE) -C ../gpu

# Source -> object rules, one per source directory. Absolute source path
# so __FILE__ expands to one too — Paths::ProjectRootOf(__FILE__) walks two
# DirNames up to find plots/, root_files/, and the log file, and silently
# collapses to CWD with a relative path.
$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $(abspath $<)

$(BUILD_DIR)/%.o: ../diagnostics/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $(abspath $<)

$(BUILD_DIR)/%.o: mains/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $(abspath $<)

$(BUILD_DIR) $(BIN_DIR):
	mkdir -p $@

# Per-binary link rules.
$(BIN_DIR)/pipeline: $(BUILD_DIR)/main_pipeline.o $(PIPELINE_OBJ) $(LIB_OBJS) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/cathode-anode: $(BUILD_DIR)/main_cathode_anode.o $(CATHODE_OBJ) $(LIB_OBJS) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/region-traces: $(BUILD_DIR)/main_region_traces.o $(REGION_OBJ) $(LIB_OBJS) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/calibrate-beam: $(BUILD_DIR)/main_calibrate_beam.o $(CALIB_OBJ) $(LIB_OBJS) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/traces: $(BUILD_DIR)/main_traces.o $(TRACES_OBJ) $(LIB_OBJS) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/diag-timing: $(BUILD_DIR)/main_diag_timing.o $(DIAG_TIM_OBJ) $(LIB_OBJS) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/diag-events: $(BUILD_DIR)/main_diag_events.o $(DIAG_EVT_OBJ) $(LIB_OBJS) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/diag-gain-drift: $(BUILD_DIR)/main_diag_gain_drift.o $(DIAG_GAIN_OBJ) $(LIB_OBJS) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/diag-subfile-drift: $(BUILD_DIR)/main_diag_subfile_drift.o $(DIAG_SUB_OBJ) $(LIB_OBJS) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/diag-subfile-rate: $(BUILD_DIR)/main_diag_subfile_rate.o $(DIAG_RATE_OBJ) $(LIB_OBJS) | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LDFLAGS)

# Convenience shortcuts so `make pipeline` etc. still work.
.PHONY: pipeline cathode-anode region-traces calibrate-beam traces diag-timing diag-events diag-gain-drift diag-subfile-drift diag-subfile-rate
pipeline:           $(BIN_DIR)/pipeline
cathode-anode:      $(BIN_DIR)/cathode-anode
region-traces:      $(BIN_DIR)/region-traces
calibrate-beam:     $(BIN_DIR)/calibrate-beam
traces:             $(BIN_DIR)/traces
diag-timing:        $(BIN_DIR)/diag-timing
diag-events:        $(BIN_DIR)/diag-events
diag-gain-drift:    $(BIN_DIR)/diag-gain-drift
diag-subfile-drift: $(BIN_DIR)/diag-subfile-drift
diag-subfile-rate:  $(BIN_DIR)/diag-subfile-rate

ALL_OBJS := $(LIB_OBJS) $(PIPELINE_OBJ) $(CATHODE_OBJ) $(REGION_OBJ) \
            $(CALIB_OBJ) $(TRACES_OBJ) \
            $(DIAG_TIM_OBJ) $(DIAG_EVT_OBJ) $(DIAG_GAIN_OBJ) $(DIAG_SUB_OBJ) \
            $(DIAG_RATE_OBJ) $(MAIN_OBJS)
-include $(ALL_OBJS:.o=.d)

clean: clean-aclic
	rm -rf $(BUILD_DIR) $(BIN_DIR)
	$(MAKE) -C ../gpu clean

clean-aclic:
	rm -f *_C.so *_C.d *_C_ACLiC_dict_rdict.pcm \
	      *_cpp.so *_cpp.d *_cpp_ACLiC_dict_rdict.pcm \
	      *_cxx.so *_cxx.d *_cxx_ACLiC_dict_rdict.pcm \
	      AutoDict_*
