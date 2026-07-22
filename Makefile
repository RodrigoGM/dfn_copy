CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -O2 -Isrc $(shell pkg-config --cflags htslib)
LDLIBS := $(shell pkg-config --libs htslib) -lz

BUILD_DIR := build
TEST_BUILD := tests/build
BIN := dfn_copy

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

$(BIN): $(BUILD_DIR)/main.o $(BUILD_DIR)/cli_args.o $(BUILD_DIR)/bins.o \
        $(BUILD_DIR)/barcode_index.o $(BUILD_DIR)/counts_matrix.o \
        $(BUILD_DIR)/read_filter.o $(BUILD_DIR)/fragment_pairing.o \
        $(BUILD_DIR)/discordant_writer.o
	$(CXX) $^ -o $@ $(LDLIBS)

$(TEST_BUILD)/test_cli_args: tests/test_cli_args.cpp $(BUILD_DIR)/cli_args.o | $(TEST_BUILD)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD_DIR)/bins.o: src/bins.cpp src/bins.hpp src/cli_args.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_BUILD)/test_bins: tests/test_bins.cpp $(BUILD_DIR)/bins.o $(BUILD_DIR)/cli_args.o | $(TEST_BUILD)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

TESTS := $(TEST_BUILD)/test_cli_args $(TEST_BUILD)/test_bins $(TEST_BUILD)/test_barcode_index $(TEST_BUILD)/test_counts_matrix $(TEST_BUILD)/test_read_filter $(TEST_BUILD)/test_fragment_pairing $(TEST_BUILD)/test_discordant_writer

.PHONY: test
test: $(TESTS)
	@for t in $(TESTS); do echo "=== $$t ==="; $$t || exit 1; done

$(BUILD_DIR)/barcode_index.o: src/barcode_index.cpp src/barcode_index.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_BUILD)/test_barcode_index: tests/test_barcode_index.cpp $(BUILD_DIR)/barcode_index.o | $(TEST_BUILD)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD_DIR)/counts_matrix.o: src/counts_matrix.cpp src/counts_matrix.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_BUILD)/test_counts_matrix: tests/test_counts_matrix.cpp $(BUILD_DIR)/counts_matrix.o | $(TEST_BUILD)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD_DIR)/read_filter.o: src/read_filter.cpp src/read_filter.hpp src/cli_args.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_BUILD)/test_read_filter: tests/test_read_filter.cpp $(BUILD_DIR)/read_filter.o $(BUILD_DIR)/cli_args.o | $(TEST_BUILD)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD_DIR)/fragment_pairing.o: src/fragment_pairing.cpp src/fragment_pairing.hpp src/cli_args.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_BUILD)/test_fragment_pairing: tests/test_fragment_pairing.cpp $(BUILD_DIR)/fragment_pairing.o $(BUILD_DIR)/cli_args.o | $(TEST_BUILD)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD_DIR)/discordant_writer.o: src/discordant_writer.cpp src/discordant_writer.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_BUILD)/test_discordant_writer: tests/test_discordant_writer.cpp $(BUILD_DIR)/discordant_writer.o | $(TEST_BUILD)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

clean:
	rm -rf $(BUILD_DIR) $(TEST_BUILD) $(BIN)
