RES_DIR = resources
DOC_DIR = docs
TEST_DIR = tests
GCOV_DIR = coverage
CC = g++
C_FLAGS = -std=c++17 -Wall -Wextra -Werror -O3
LD_FLAGS =
CPP_SOURCES = $(wildcard *.cc)
CPP_OBJECTS = $(patsubst %.cc,$(RES_DIR)/%.o,$(CPP_SOURCES))
TEST_SOURCES = $(wildcard tests/*.cc)
TEST_OBJECTS = $(patsubst %.cc,$(RES_DIR)/%.o,$(TEST_SOURCES))

all: test

ifeq ($(MAKECMDGOALS), gcov_report)
C_FLAGS += --coverage
LD_FLAGS += --coverage
endif

test: $(RES_DIR) $(TEST_OBJECTS)
	$(CC) $(LD_FLAGS) $(TEST_OBJECTS) -o test -lgtest
	./test

$(RES_DIR):
	mkdir -p $(RES_DIR)
	mkdir -p $(RES_DIR)/$(TEST_DIR)

$(RES_DIR)/%.o: %.cc
	$(CC) $(C_FLAGS) -c $< -o $@

clean:
	rm -rf $(RES_DIR) $(DOC_DIR) $(GCOV_DIR) test

gcov_report: test
	mkdir -p $(GCOV_DIR)
	lcov -c --directory . --output-file $(GCOV_DIR)/coverage.info --ignore-errors mismatch \
	--include "$(pwd)/matrix.tpp" --include "$(pwd)/matrix.h"
	genhtml $(GCOV_DIR)/coverage.info --output-directory $(GCOV_DIR)/html

vg: test
	valgrind --track-fds=yes --trace-children=yes --track-origins=yes --leak-check=full --show-leak-kinds=all -s ./test

.PHONY: all test clean gcov_report vg