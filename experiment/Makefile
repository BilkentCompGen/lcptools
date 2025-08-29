FASTA := ../data/chm13v2.0.fa
FASTA2 := ../data/chm13v2.0.chr1.fa
FASTQ := ../data/chm13v2.0.chr1.fq
MAF := ../data/chm13v2.0.chr1.maf

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
# fq - alignment exp
FQCMIN := fqcmin
FQLCP := fqlcp
FQMIN := fqmin
FQSYNC := fqsync
# fq2 - hit experiment
FQLCP2 := fqlcp2
FQMIN2 := fqmin2
FQSYNC2 := fqsync2

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
TIME := /usr/bin/time -v

.PHONY: all

all: mkdir_bin mkdir_out mkdir_branch

mkdir_bin:
	@echo "Checking if $(EXECUTABLE_DIR) exists..."
	@if [ ! -d "$(EXECUTABLE_DIR)" ]; then \
		echo "$(EXECUTABLE_DIR) does not exist. Creating..."; \
		mkdir $(EXECUTABLE_DIR); \
	else \
		echo "$(EXECUTABLE_DIR) already exists."; \
	fi

mkdir_out:
	@echo "Checking if $(OUT_DIR) exists..."
	@if [ ! -d "$(OUT_DIR)" ]; then \
		echo "$(OUT_DIR) does not exist. Creating..."; \
		mkdir $(OUT_DIR); \
	else \
		echo "$(OUT_DIR) already exists."; \
	fi
	@if [ ! -d "$(OUT_DIR)/$(BRANCH)" ]; then \
		echo "$(OUT_DIR)/$(BRANCH) does not exist. Creating..."; \
		mkdir $(OUT_DIR)/$(BRANCH); \
	else \
		echo "$(OUT_DIR)/$(BRANCH) already exists."; \
	fi

check_lcptools:
	@if [ ! -d "$(INCLUDE_DIR)" ]; then \
		echo "$(INCLUDE_DIR) does not exist. Please compile lcptools to headers to $(INCLUDE_DIR)"; \
		exit 1; \
	fi
	@if [ ! -d "$(LIB_DIR)" ]; then \
		echo "$(LIB_DIR) does not exist. Please compile lcptools and put library to $(LIB_DIR)"; \
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
faLcp: mkdir_bin mkdir_out check_lcptools
	@echo "Compiling $(FALCP)$(CXX)"
	$(GXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -c $(FALCP)$(CXX)
	@echo "Linking lcptoolsS library and binary files"
	$(GXX) $(CXXFLAGS) -o $(FALCP) $(FALCP).o -L$(LIB_DIR) -llcptools -Wl,-rpath,$(LIB_DIR)
	rm $(FALCP).o
	@echo "Moving $(FALCP) to $(EXECUTABLE_DIR)"
	mv $(FALCP) $(EXECUTABLE_DIR)
	@echo "Preprocessing $(FASTA). Output will be put into $(OUT_DIR)/$(BRANCH)/$(FALCP)-output.txt"
	$(TIME) $(EXECUTABLE_DIR)/$(FALCP) $(FASTA) > $(OUT_DIR)/$(BRANCH)/$(FALCP)-output.txt 2>&1

faLcpMem: mkdir_bin mkdir_out check_fasta check_lcptools
	@echo "Compiling $(FALCPMEM)$(CXX)"
	$(GXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -c $(FALCPMEM)$(CXX)
	@echo "Linking lcptoolsS library and binary files"
	$(GXX) $(CXXFLAGS) -o $(FALCPMEM) $(FALCPMEM).o -L$(LIB_DIR) -llcptools -Wl,-rpath,$(LIB_DIR)
	rm $(FALCPMEM).o
	@echo "Moving $(FALCPMEM) to $(EXECUTABLE_DIR)"
	mv $(FALCPMEM) $(EXECUTABLE_DIR)
	@echo "Preprocessing $(FASTA). Output will be put into $(OUT_DIR)/$(BRANCH)/$(FALCPMEM)-output.txt"
	$(TIME) $(EXECUTABLE_DIR)/$(FALCPMEM) $(FASTA) > $(OUT_DIR)/$(BRANCH)/$(FALCPMEM)-output.txt 2>&1

fqLcp: mkdir_bin mkdir_out
	@echo "Compiling $(FQLCP)$(CXX)"
	$(GXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -c $(FQLCP)$(CXX)
	$(GXX) $(CXXFLAGS) -o $(FQLCP) $(FQLCP).o -L$(LIB_DIR) -llcptools -Wl,-rpath,$(LIB_DIR)
	rm $(FQLCP).o
	@echo "Moving $(FQLCP) to $(EXECUTABLE_DIR)"
	mv $(FQLCP) $(EXECUTABLE_DIR)
	@echo "Preprocessing $(FASTA). Output will be put into $(OUT_DIR)/$(BRANCH)/$(FQLCP)-output.txt"
	$(TIME) $(EXECUTABLE_DIR)/$(FQLCP) $(MAF) $(FASTQ) > $(OUT_DIR)/$(BRANCH)/$(FQLCP)-output.txt 2>&1

fqLcp2: mkdir_bin mkdir_out
	@echo "Compiling $(FQLCP2)$(CXX)"
	$(GXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -c $(FQLCP2)$(CXX)
	$(GXX) $(CXXFLAGS) -o $(FQLCP2) $(FQLCP2).o -L$(LIB_DIR) -llcptools -Wl,-rpath,$(LIB_DIR)
	rm $(FQLCP2).o
	@echo "Moving $(FQLCP2) to $(EXECUTABLE_DIR)"
	mv $(FQLCP2) $(EXECUTABLE_DIR)
	@echo "Preprocessing $(FASTA). Output will be put into $(OUT_DIR)/$(BRANCH)/$(FQLCP2)-output.txt"
	$(TIME) $(EXECUTABLE_DIR)/$(FQLCP2) $(FASTA2) $(MAF) $(FASTQ) > $(OUT_DIR)/$(BRANCH)/$(FQLCP2)-output.txt 2>&1

######################################################################
# MINIMIZER
######################################################################
faMin: mkdir_bin mkdir_out
	@echo "Compiling $(FAMIN)$(CXX)"
	$(GXX) $(CXXFLAGS) -c $(FAMIN)$(CXX)
	$(GXX) $(CXXFLAGS) -o $(FAMIN) $(FAMIN).o
	rm $(FAMIN).o
	@echo "Moving $(FAMIN) to $(EXECUTABLE_DIR)"
	mv $(FAMIN) $(EXECUTABLE_DIR)
	@echo "Preprocessing $(FASTA). Output will be put into $(OUT_DIR)/$(BRANCH)/$(FAMIN)-output.txt"
	$(TIME) $(EXECUTABLE_DIR)/$(FAMIN) $(FASTA) $(RUN_ARGS) > $(OUT_DIR)/$(BRANCH)/$(FAMIN)-$(RUN_ARGS_JOINED)-output.txt 2>&1

fqMin: mkdir_bin mkdir_out
	@echo "Compiling $(FQMIN)$(CXX)"
	$(GXX) $(CXXFLAGS) -c $(FQMIN)$(CXX)
	$(GXX) $(CXXFLAGS) -o $(FQMIN) $(FQMIN).o
	rm $(FQMIN).o
	@echo "Moving $(FQMIN) to $(EXECUTABLE_DIR)"
	mv $(FQMIN) $(EXECUTABLE_DIR)
	@echo "Preprocessing $(FASTA). Output will be put into $(OUT_DIR)/$(BRANCH)/$(FQMIN)-output.txt"
	$(TIME) $(EXECUTABLE_DIR)/$(FQMIN) $(MAF) $(FASTQ) > $(OUT_DIR)/$(BRANCH)/$(FQMIN)-output.txt 2>&1

fqMin2: mkdir_bin mkdir_out
	@echo "Compiling $(FQMIN2)$(CXX)"
	$(GXX) $(CXXFLAGS) -c $(FQMIN2)$(CXX)
	$(GXX) $(CXXFLAGS) -o $(FQMIN2) $(FQMIN2).o
	rm $(FQMIN2).o
	@echo "Moving $(FQMIN2) to $(EXECUTABLE_DIR)"
	mv $(FQMIN2) $(EXECUTABLE_DIR)
	@echo "Preprocessing $(FASTA). Output will be put into $(OUT_DIR)/$(BRANCH)/$(FQMIN2)-output.txt"
	$(TIME) $(EXECUTABLE_DIR)/$(FQMIN2) $(FASTA2) $(MAF) $(FASTQ) > $(OUT_DIR)/$(BRANCH)/$(FQMIN2)-output.txt 2>&1

######################################################################
# Canonical kmer - MINIMIZER
######################################################################
faCMin: mkdir_bin mkdir_out
	@echo "Compiling $(FACMIN)$(CXX)"
	$(GXX) $(CXXFLAGS) -c $(FACMIN)$(CXX)
	$(GXX) $(CXXFLAGS) -o $(FACMIN) $(FACMIN).o
	rm $(FACMIN).o
	@echo "Moving $(FACMIN) to $(EXECUTABLE_DIR)"
	mv $(FACMIN) $(EXECUTABLE_DIR)
	@echo "Preprocessing $(FASTA). Output will be put into $(OUT_DIR)/$(BRANCH)/$(FACMIN)-output.txt"
	$(TIME) $(EXECUTABLE_DIR)/$(FACMIN) $(FASTA) $(RUN_ARGS) > $(OUT_DIR)/$(BRANCH)/$(FACMIN)-$(RUN_ARGS_JOINED)-output.txt 2>&1

fqCMin: mkdir_bin mkdir_out
	@echo "Compiling $(FQCMIN)$(CXX)"
	$(GXX) $(CXXFLAGS) -c $(FQCMIN)$(CXX)
	$(GXX) $(CXXFLAGS) -o $(FQCMIN) $(FQCMIN).o
	rm $(FQCMIN).o
	@echo "Moving $(FQCMIN) to $(EXECUTABLE_DIR)"
	mv $(FQCMIN) $(EXECUTABLE_DIR)
	@echo "Preprocessing $(FASTA). Output will be put into $(OUT_DIR)/$(BRANCH)/$(FQCMIN)-output.txt"
	$(TIME) $(EXECUTABLE_DIR)/$(FQCMIN) $(MAF) $(FASTQ) > $(OUT_DIR)/$(BRANCH)/$(FQCMIN)-output.txt 2>&1
######################################################################
# SYNCMER
######################################################################
faSync: mkdir_bin mkdir_out check_fasta
	@echo "Compiling $(FASYNC)$(CXX)"
	$(GXX) $(CXXFLAGS) -c $(FASYNC)$(CXX)
	$(GXX) $(CXXFLAGS) -o $(FASYNC) $(FASYNC).o -L$(LIB_DIR) -llcptools -Wl,-rpath,$(LIB_DIR)
	rm $(FASYNC).o
	@echo "Moving $(FASYNC) to $(EXECUTABLE_DIR)"
	mv $(FASYNC) $(EXECUTABLE_DIR)
	@echo "Preprocessing $(FASTA). Output will be put into $(OUT_DIR)/$(BRANCH)/$(FASYNC)-output.txt"
	$(TIME) $(EXECUTABLE_DIR)/$(FASYNC) $(FASTA) $(RUN_ARGS) > $(OUT_DIR)/$(BRANCH)/$(FASYNC)-$(RUN_ARGS_JOINED)-output.txt 2>&1

fqSync: mkdir_bin mkdir_out
	@echo "Compiling $(FQSYNC)$(CXX)"
	$(GXX) $(CXXFLAGS) -c $(FQSYNC)$(CXX)
	$(GXX) $(CXXFLAGS) -o $(FQSYNC) $(FQSYNC).o
	rm $(FQSYNC).o
	@echo "Moving $(FQSYNC) to $(EXECUTABLE_DIR)"
	mv $(FQSYNC) $(EXECUTABLE_DIR)
	@echo "Preprocessing $(FASTA). Output will be put into $(OUT_DIR)/$(BRANCH)/$(FQSYNC)-output.txt"
	$(TIME) $(EXECUTABLE_DIR)/$(FQSYNC) $(MAF) $(FASTQ) > $(OUT_DIR)/$(BRANCH)/$(FQSYNC)-output.txt 2>&1

fqSync2: mkdir_bin mkdir_out
	@echo "Compiling $(FQSYNC2)$(CXX)"
	$(GXX) $(CXXFLAGS) -c $(FQSYNC2)$(CXX)
	$(GXX) $(CXXFLAGS) -o $(FQSYNC2) $(FQSYNC2).o
	rm $(FQSYNC2).o
	@echo "Moving $(FQSYNC2) to $(EXECUTABLE_DIR)"
	mv $(FQSYNC2) $(EXECUTABLE_DIR)
	@echo "Preprocessing $(FASTA). Output will be put into $(OUT_DIR)/$(BRANCH)/$(FQSYNC2)-output.txt"
	$(TIME) $(EXECUTABLE_DIR)/$(FQSYNC2) $(FASTA2) $(MAF) $(FASTQ) > $(OUT_DIR)/$(BRANCH)/$(FQSYNC2)-output.txt 2>&1
