# g++
CXX = gcc
CXXFLAGS = -O3 -Wall -Wextra -Wpedantic
CXXEXTRA = -fPIC

# archiver and flags
AR = ar
ARFLAGS = rcs

# variables
SRC = encoding.c core.c lps.c
HDR = $(SRC:.c=.h)
OBJ_STATIC = $(SRC:.c=_s.o)
OBJ_DYNAMIC = $(SRC:.c=_d.o)

# test files
TEST_DIR = tests
TESTS = $(patsubst $(TEST_DIR)/%.cpp,%,$(wildcard $(TEST_DIR)/*.cpp))

# library names
LIB_NAME = lcptools
STATIC = lib$(LIB_NAME).a
DYNAMIC = lib$(LIB_NAME).so

PREFIX ?= /usr/local
ABS_PREFIX := $(realpath $(PREFIX))
INCLUDE_DIR = $(ABS_PREFIX)/include
LIB_DIR = $(ABS_PREFIX)/lib

.PHONY: all clean install uninstall test

install: clean $(STATIC) $(DYNAMIC) lcptools
	mkdir -p $(INCLUDE_DIR)
	rm -f *.o
	cp $(HDR) $(INCLUDE_DIR)
	@echo "[[WARNING]]! Please make sure that $(LIB_DIR) included in LD_LIBRARY_PATH if you want to include dynamic library";

uninstall:
	rm -f lcptools;
	rm -f $(LIB_DIR)/$(STATIC)
	rm -f $(LIB_DIR)/$(DYNAMIC)
	@for hdr in $(HDR); do \
		echo "Removing $(INCLUDE_DIR)/$$hdr;"; \
		rm -f $(INCLUDE_DIR)/$$hdr; \
	done

clean:
	rm -f $(TEST_DIR)/*.o 
	@echo "rm $(LIB_DIR)/$(STATIC)";
	@if [ -f "$(LIB_DIR)/$(STATIC)" ]; then \
		rm -f $(LIB_DIR)/$(STATIC) || \
			{ \
				echo "Couldn't remove $(LIB_DIR)/$(STATIC)"; \
				exit 1; \
			}; \
	fi
	@echo "rm $(LIB_DIR)/$(DYNAMIC)";
	@if [ -f "$(LIB_DIR)/$(DYNAMIC)" ]; then \
		rm -f $(LIB_DIR)/$(DYNAMIC) || \
			{ \
				echo "Couldn't remove $(LIB_DIR)/$(DYNAMIC)"; \
				exit 1; \
			}; \
	fi
	rm -f $(OBJ_STATIC) 
	rm -f $(OBJ_DYNAMIC) 
	rm -f lcptools

# target for static library
$(STATIC): $(OBJ_STATIC)
	$(AR) $(ARFLAGS) $@ $^
	rm -f $(OBJ_STATIC)
	mkdir -p $(LIB_DIR)
	@mv $@ $(LIB_DIR) || \
		{ \
			echo "Couldn't move $@ to $(LIB_DIR)"; \
			exit 1; \
		}

# target for dynamic library
$(DYNAMIC): $(OBJ_DYNAMIC)
	$(CXX) -shared -o $@ $^ -lm
	rm -f $(OBJ_DYNAMIC)
	mkdir -p $(LIB_DIR)
	@mv $@ $(LIB_DIR) || \
		{ \
			echo "Couldn't move $@ to $(LIB_DIR)"; \
			exit 1; \
		}

# target to compile lcptools executable
lcptools: $(SRC) $(HDR)
	rm -f $@
	$(CXX) $(CXXFLAGS) -I$(INCLUDE_DIR) -c $@.c -o $@.o
	$(CXX) $(CXXFLAGS) -o $@ $@.o -L$(LIB_DIR) -llcptools -Wl,-rpath,$(LIB_DIR)
	chmod +x $@

# rule to compile .c files to .o files for static library
%_s.o: %.c
	$(CXX) $(CXXFLAGS) -c $< -o $@

# rule to compile .c files to .o files for dynamic library
%_d.o: %.c
	$(CXX) $(CXXFLAGS) -c $(CXXEXTRA) $< -o $@

# run all tests
test:
	@echo "Running tests..."
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
