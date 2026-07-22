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

$(BUILD_DIR)/cli_args.o: src/cli_args.cpp src/cli_args.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BIN): $(BUILD_DIR)/main.o $(BUILD_DIR)/cli_args.o
	$(CXX) $(BUILD_DIR)/main.o $(BUILD_DIR)/cli_args.o -o $@ $(LDLIBS)

$(TEST_BUILD)/test_cli_args: tests/test_cli_args.cpp $(BUILD_DIR)/cli_args.o | $(TEST_BUILD)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD_DIR)/bins.o: src/bins.cpp src/bins.hpp src/cli_args.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_BUILD)/test_bins: tests/test_bins.cpp $(BUILD_DIR)/bins.o $(BUILD_DIR)/cli_args.o | $(TEST_BUILD)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

TESTS := $(TEST_BUILD)/test_cli_args $(TEST_BUILD)/test_bins $(TEST_BUILD)/test_barcode_index

.PHONY: test
test: $(TESTS)
	@for t in $(TESTS); do echo "=== $$t ==="; $$t || exit 1; done

$(BUILD_DIR)/barcode_index.o: src/barcode_index.cpp src/barcode_index.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_BUILD)/test_barcode_index: tests/test_barcode_index.cpp $(BUILD_DIR)/barcode_index.o | $(TEST_BUILD)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

clean:
	rm -rf $(BUILD_DIR) $(TEST_BUILD) $(BIN)
