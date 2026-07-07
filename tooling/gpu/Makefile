NVCC := nvcc

# Compute capability. 
GPU_ARCH ?= $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d '.')
ifeq ($(strip $(GPU_ARCH)),)
GPU_ARCH := 89
endif

NVCC_FLAGS := -shared -Xcompiler -fPIC --expt-extended-lambda -O3 \
              -arch=sm_$(GPU_ARCH)

all: libgpuaccel.so

libgpuaccel.so: gpu_sort.cu gpu_interface.h
	$(NVCC) $(NVCC_FLAGS) -o $@ gpu_sort.cu 

clean:
	rm -f libgpuaccel.so

.PHONY: all clean
