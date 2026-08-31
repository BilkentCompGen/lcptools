# g++
CXX = gcc
CXXFLAGS = -O3 -Wall -Wextra -Wpedantic
CXXEXTRA = -fPIC

# compile-time configuration
LABEL ?= 32
POS ?= 32
DCT ?= 1
CORE ?= var

ifeq ($(filter 0 32 64,$(LABEL)),)
$(error LABEL must be 0, 32 or 64, got '$(LABEL)')
endif

ifeq ($(filter 0 32 64,$(POS)),)
$(error POS must be 0, 32 or 64, got '$(POS)')
endif

ifeq ($(filter var fixed,$(CORE)),)
$(error CORE must be var or fixed, got '$(CORE)')
endif

ifeq ($(CORE),fixed)
FIXED_CORE = 1
FACTOR = 1
else
FIXED_CORE = 0
FACTOR = 1.5
endif

# Optional name suffix, empty by default so a plain build installs exactly what it always has
SUFFIX ?=

# archiver and flags
AR = ar
ARFLAGS = rcs

# platform
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
SOEXT = dylib
LDPATH_VAR = DYLD_LIBRARY_PATH
else
SOEXT = so
LDPATH_VAR = LD_LIBRARY_PATH
endif

# variables
SRC = encoding.c core.c lps.c
CONFIG_HDR = config.h
CONFIG_IN = $(CONFIG_HDR).in
HDR = $(SRC:.c=.h) $(CONFIG_HDR)
OBJ_STATIC = $(SRC:.c=_s.o)
OBJ_DYNAMIC = $(SRC:.c=_d.o)

# test files
TEST_DIR = tests
TESTS = $(patsubst $(TEST_DIR)/%.cpp,%,$(wildcard $(TEST_DIR)/*.cpp))

# library names
LIB_NAME = lcptools$(SUFFIX)
STATIC = lib$(LIB_NAME).a
DYNAMIC = lib$(LIB_NAME).$(SOEXT)

VARIANT = $(CORE)-l$(LABEL)-p$(POS)-d$(DCT)

PREFIX ?= /usr/local
ABS_PREFIX := $(realpath $(PREFIX))
INCLUDE_DIR = $(ABS_PREFIX)/include
LIB_DIR = $(ABS_PREFIX)/lib

# flags for linking the shared library, once LIB_DIR is known
ifeq ($(UNAME_S),Darwin)
SHAREDFLAGS = -dynamiclib -install_name $(LIB_DIR)/$(DYNAMIC)
else
SHAREDFLAGS = -shared -Wl,-soname,$(DYNAMIC)
endif

.PHONY: all clean install uninstall test

install: clean $(STATIC) $(DYNAMIC)
	@mkdir -p $(INCLUDE_DIR) $(LIB_DIR)
	@cp $(HDR) $(INCLUDE_DIR)
	@rm -f *.o $(CONFIG_HDR)
	@echo "Installed $(STATIC) and $(DYNAMIC) to $(LIB_DIR) (variant $(VARIANT))"
	@echo "[[WARNING]]! Please make sure that $(LIB_DIR) included in $(LDPATH_VAR) if you want to include dynamic library";

uninstall:
	@rm -f $(LIB_DIR)/$(STATIC)
	@rm -f $(LIB_DIR)/lib$(LIB_NAME).so $(LIB_DIR)/lib$(LIB_NAME).dylib
	@for hdr in $(HDR); do \
		echo "Removing $(INCLUDE_DIR)/$$hdr;"; \
		rm -f $(INCLUDE_DIR)/$$hdr; \
	done

clean:
	@rm -f $(TEST_DIR)/*.o 
	@if [ -f "$(LIB_DIR)/$(STATIC)" ]; then \
		echo "rm $(LIB_DIR)/$(STATIC)"; \
		rm -f $(LIB_DIR)/$(STATIC) || \
			{ \
				echo "Couldn't remove $(LIB_DIR)/$(STATIC)"; \
				exit 1; \
			}; \
	fi
	@for so in $(LIB_DIR)/lib$(LIB_NAME).so $(LIB_DIR)/lib$(LIB_NAME).dylib; do \
		if [ -f "$$so" ]; then \
			echo "rm $$so"; \
			rm -f "$$so" || { echo "Couldn't remove $$so"; exit 1; }; \
		fi; \
	done
	@rm -f $(OBJ_STATIC)
	@rm -f $(OBJ_DYNAMIC)
	@rm -f $(CONFIG_HDR)
	@rm -f $(INCLUDE_DIR)/*

$(CONFIG_HDR): $(CONFIG_IN) Makefile
	@echo "Generating $@ (LABEL=$(LABEL) POS=$(POS) DCT=$(DCT) CORE=$(CORE))"
	@sed -e 's/@LABEL@/$(LABEL)/' \
	     -e 's/@POS@/$(POS)/' \
	     -e 's/@DCT@/$(DCT)/' \
	     -e 's/@FIXED_CORE@/$(FIXED_CORE)/' \
	     -e 's/@FACTOR@/$(FACTOR)/' $< > $@

# target for static library
$(STATIC): $(OBJ_STATIC)
	$(AR) $(ARFLAGS) $@ $^
	@rm -f $(OBJ_STATIC)
	@mkdir -p $(LIB_DIR)
	@echo "mv $@ $(LIB_DIR)"
	@mv $@ $(LIB_DIR) || \
		{ \
			echo "Couldn't move $@ to $(LIB_DIR)"; \
			exit 1; \
		}

# target for dynamic library
$(DYNAMIC): $(OBJ_DYNAMIC)
	$(CXX) $(SHAREDFLAGS) -o $@ $^ -lm -lpthread
	@rm -f $(OBJ_DYNAMIC)
	@mkdir -p $(LIB_DIR)
	@echo "mv $@ $(LIB_DIR)"
	@mv $@ $(LIB_DIR) || \
		{ \
			echo "Couldn't move $@ to $(LIB_DIR)"; \
			exit 1; \
		}

# rule to compile .c files to .o files for static library
%_s.o: %.c $(CONFIG_HDR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# rule to compile .c files to .o files for dynamic library
%_d.o: %.c $(CONFIG_HDR)
	$(CXX) $(CXXFLAGS) $(CXXEXTRA) -c $< -o $@

# run all tests
test:
	@echo "Running tests (variant $(VARIANT))..."
	@for test in $(TESTS); do \
		echo "Compiling $$test.cpp..."; \
		g++ $(CXXFLAGS) -I$(INCLUDE_DIR) -o tests/$$test tests/$$test.cpp -L$(LIB_DIR) -l$(LIB_NAME) -Wl,-rpath,$(LIB_DIR); \
		if [ $$? -ne 0 ]; then \
			echo "Compilation failed for $$test.c"; \
			exit 1; \
		fi; \
		tests/$$test || exit 1; \
		rm -f tests/$$test; \
	done
	@echo "All tests passed."
