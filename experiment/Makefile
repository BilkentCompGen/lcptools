FASTA := ../data/chm13v2.0.fa

# platform
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
SOEXT = dylib
LDPATH_VAR = DYLD_LIBRARY_PATH
TIME := /usr/bin/time -l
else
SOEXT = so
LDPATH_VAR = LD_LIBRARY_PATH
TIME := /usr/bin/time -v
endif

ifeq ($(firstword $(MAKECMDGOALS)),faMin)
  TARGET := faMin
else ifeq ($(firstword $(MAKECMDGOALS)),faCMin)
  TARGET := faCMin
else ifeq ($(firstword $(MAKECMDGOALS)),faSync)
  TARGET := faSync
endif

BRANCH := $(shell git rev-parse --abbrev-ref HEAD)
HASH := $(shell git rev-parse HEAD)

ifneq ($(TARGET),)
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  $(eval $(RUN_ARGS):;@:)
  ifeq ($(TARGET),faSync)
    RUN_ARGS_JOINED := $(word 1, $(RUN_ARGS))-$(word 2, $(RUN_ARGS))-$(word 3, $(RUN_ARGS))
  else
    RUN_ARGS_JOINED := $(word 1, $(RUN_ARGS))-$(word 2, $(RUN_ARGS))
  endif
endif

# Programs
# fa
FACMIN := facmin
FALCPMEM := falcp-mem
FALCP := falcp
FAMIN := famin
FASYNC := fasync

# Directories
CURRENT_DIR := $(shell pwd)
DATA_DIR := $(CURRENT_DIR)/../data
EXECUTABLE_DIR := $(CURRENT_DIR)/../bin
OUT_DIR := $(CURRENT_DIR)/../out
BRANCH_DIR := $(OUT_DIR)/$(BRANCH)
INCLUDE_DIR ?= $(CURRENT_DIR)/../include
LIB_DIR ?= $(CURRENT_DIR)/../lib

# Extension
CXX := .cpp

# Compiler
GXX := g++
CXXFLAGS = -std=c++11 -O3 -Wall -Wextra -Wpedantic
LCPTOOLS_FLAGS = -L$(LIB_DIR) -llcptools -Wl,-rpath,$(LIB_DIR)

.PHONY: all mkdir_bin mkdir_out check_lcptools check_fasta clean faLcp faLcpMem faMin faCMin faSync

all: mkdir_bin mkdir_out

mkdir_bin:
	@echo "Checking if $(EXECUTABLE_DIR) exists..."
	@mkdir -p $(EXECUTABLE_DIR)

mkdir_out:
	@echo "Checking if $(OUT_DIR)/$(BRANCH) exists..."
	@mkdir -p $(OUT_DIR)/$(BRANCH)

check_lcptools:
	@if [ ! -d "$(INCLUDE_DIR)" ]; then \
		echo "$(INCLUDE_DIR) does not exist. Please compile lcptools to headers to $(INCLUDE_DIR)"; \
		exit 1; \
	fi
	@if [ ! -d "$(LIB_DIR)" ]; then \
		echo "$(LIB_DIR) does not exist. Please compile lcptools and put library to $(LIB_DIR)"; \
		exit 1; \
	fi
	@if [ ! -f "$(LIB_DIR)/liblcptools.a" ] && [ ! -f "$(LIB_DIR)/liblcptools.$(SOEXT)" ]; then \
		echo "No liblcptools.a or liblcptools.$(SOEXT) in $(LIB_DIR)."; \
		if ls $(LIB_DIR)/liblcptools.* >/dev/null 2>&1; then \
			echo "Found instead:"; ls $(LIB_DIR)/liblcptools.*; \
			echo "This is $(UNAME_S); a shared lcptools here must be liblcptools.$(SOEXT)."; \
		fi; \
		echo "Build it with: make -C .. install PREFIX=.."; \
		exit 1; \
	fi

check_fasta:
	@if [ -z "$(FASTA)" ]; then \
        echo "FASTA is not set or is empty. Exiting"; \
		exit 1; \
    fi
	@if [ ! -f $(FASTA) ]; then \
		echo "$(FASTA) does not exist. Exiting."; \
		exit 1; \
	fi

clean:
	@if [ -d "$(OUT_DIR)" ]; then \
		echo "Removing $(OUT_DIR)"; \
		rm -r $(OUT_DIR); \
	fi
	@if [ -d "$(EXECUTABLE_DIR)" ]; then \
		echo "Removing $(EXECUTABLE_DIR)"; \
		rm -r $(EXECUTABLE_DIR); \
	fi
	@echo "Cleaned."

######################################################################
# LCP
######################################################################
faLcp: mkdir_bin mkdir_out check_fasta check_lcptools
	@echo "Compiling $(FALCP)$(CXX)"
	$(GXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -c $(FALCP)$(CXX)
	@echo "Linking lcptoolsS library and binary files"
	$(GXX) $(CXXFLAGS) -o $(FALCP) $(FALCP).o $(LCPTOOLS_FLAGS)
	@rm -f $(FALCP).o
	@echo "Moving $(FALCP) to $(EXECUTABLE_DIR)"
	@mv -f $(FALCP) $(EXECUTABLE_DIR)
	@echo "Preprocessing $(FASTA). Output will be put into $(OUT_DIR)/$(BRANCH)/$(FALCP)-output.txt"
	$(TIME) $(EXECUTABLE_DIR)/$(FALCP) $(FASTA) > $(OUT_DIR)/$(BRANCH)/$(FALCP)-output.txt 2>&1

faLcpMem: mkdir_bin mkdir_out check_fasta check_lcptools
	@echo "Compiling $(FALCPMEM)$(CXX)"
	$(GXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -c $(FALCPMEM)$(CXX)
	@echo "Linking lcptoolsS library and binary files"
	$(GXX) $(CXXFLAGS) -o $(FALCPMEM) $(FALCPMEM).o $(LCPTOOLS_FLAGS)
	@rm -f $(FALCPMEM).o
	@echo "Moving $(FALCPMEM) to $(EXECUTABLE_DIR)"
	@mv -f $(FALCPMEM) $(EXECUTABLE_DIR)
	@echo "Preprocessing $(FASTA). Output will be put into $(OUT_DIR)/$(BRANCH)/$(FALCPMEM)-output.txt"
	$(TIME) $(EXECUTABLE_DIR)/$(FALCPMEM) $(FASTA) > $(OUT_DIR)/$(BRANCH)/$(FALCPMEM)-output.txt 2>&1

######################################################################
# MINIMIZER
######################################################################
faMin: mkdir_bin mkdir_out
	@echo "Compiling $(FAMIN)$(CXX)"
	$(GXX) $(CXXFLAGS) -c $(FAMIN)$(CXX)
	$(GXX) $(CXXFLAGS) -o $(FAMIN) $(FAMIN).o
	@rm -f $(FAMIN).o
	@echo "Moving $(FAMIN) to $(EXECUTABLE_DIR)"
	@mv -f $(FAMIN) $(EXECUTABLE_DIR)
	@echo "Preprocessing $(FASTA). Output will be put into $(OUT_DIR)/$(BRANCH)/$(FAMIN)-output.txt"
	$(TIME) $(EXECUTABLE_DIR)/$(FAMIN) $(FASTA) $(RUN_ARGS) > $(OUT_DIR)/$(BRANCH)/$(FAMIN)-$(RUN_ARGS_JOINED)-output.txt 2>&1

######################################################################
# Canonical kmer - MINIMIZER
######################################################################
faCMin: mkdir_bin mkdir_out
	@echo "Compiling $(FACMIN)$(CXX)"
	$(GXX) $(CXXFLAGS) -c $(FACMIN)$(CXX)
	$(GXX) $(CXXFLAGS) -o $(FACMIN) $(FACMIN).o
	@rm -f $(FACMIN).o
	@echo "Moving $(FACMIN) to $(EXECUTABLE_DIR)"
	@mv -f $(FACMIN) $(EXECUTABLE_DIR)
	@echo "Preprocessing $(FASTA). Output will be put into $(OUT_DIR)/$(BRANCH)/$(FACMIN)-output.txt"
	$(TIME) $(EXECUTABLE_DIR)/$(FACMIN) $(FASTA) $(RUN_ARGS) > $(OUT_DIR)/$(BRANCH)/$(FACMIN)-$(RUN_ARGS_JOINED)-output.txt 2>&1

######################################################################
# SYNCMER
######################################################################
faSync: mkdir_bin mkdir_out check_fasta
	@echo "Compiling $(FASYNC)$(CXX)"
	$(GXX) $(CXXFLAGS) -c $(FASYNC)$(CXX)
	$(GXX) $(CXXFLAGS) -o $(FASYNC) $(FASYNC).o
	@rm -f $(FASYNC).o
	@echo "Moving $(FASYNC) to $(EXECUTABLE_DIR)"
	@mv -f $(FASYNC) $(EXECUTABLE_DIR)
	@echo "Preprocessing $(FASTA). Output will be put into $(OUT_DIR)/$(BRANCH)/$(FASYNC)-output.txt"
	$(TIME) $(EXECUTABLE_DIR)/$(FASYNC) $(FASTA) $(RUN_ARGS) > $(OUT_DIR)/$(BRANCH)/$(FASYNC)-$(RUN_ARGS_JOINED)-output.txt 2>&1
