CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -O2 -Isrc $(shell pkg-config --cflags htslib)
LDLIBS := $(shell pkg-config --libs htslib) -lz

BUILD_DIR := build
TEST_BUILD := tests/build
BIN := dfn_bin

.PHONY: all clean test

all: $(BIN)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(TEST_BUILD):
	mkdir -p $(TEST_BUILD)

$(BUILD_DIR)/main.o: src/main.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BIN): $(BUILD_DIR)/main.o
	$(CXX) $(BUILD_DIR)/main.o -o $@ $(LDLIBS)

test:
	@echo "no tests yet"

clean:
	rm -rf $(BUILD_DIR) $(TEST_BUILD) $(BIN)
