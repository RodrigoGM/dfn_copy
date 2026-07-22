CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -O2 -Isrc $(shell pkg-config --cflags htslib)
LDLIBS := $(shell pkg-config --libs htslib) -lz

BUILD_DIR := build
TEST_BUILD := tests/build
BIN := dfn_copy
CBS_BIN := dfn_cbs

.PHONY: all clean test

all: $(BIN) $(CBS_BIN)

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

$(BUILD_DIR)/cbs_main.o: src/cbs_main.cpp | $(BUILD_DIR)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc -c $< -o $@

$(BUILD_DIR)/cbs_args.o: src/cbs_args.cpp src/cbs_args.hpp | $(BUILD_DIR)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc -c $< -o $@

$(TEST_BUILD)/test_cbs_args: tests/test_cbs_args.cpp $(BUILD_DIR)/cbs_args.o | $(TEST_BUILD)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc $^ -o $@

$(BUILD_DIR)/bin_gc.o: src/bin_gc.cpp src/bin_gc.hpp | $(BUILD_DIR)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc -c $< -o $@

$(TEST_BUILD)/test_bin_gc: tests/test_bin_gc.cpp $(BUILD_DIR)/bin_gc.o | $(TEST_BUILD)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc $^ -o $@

$(BUILD_DIR)/raw_counts_loader.o: src/raw_counts_loader.cpp src/raw_counts_loader.hpp | $(BUILD_DIR)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc -c $< -o $@

$(TEST_BUILD)/test_raw_counts_loader: tests/test_raw_counts_loader.cpp $(BUILD_DIR)/raw_counts_loader.o | $(TEST_BUILD)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc $^ -o $@ -lz

$(BUILD_DIR)/cell_filter.o: src/cell_filter.cpp src/cell_filter.hpp src/raw_counts_loader.hpp | $(BUILD_DIR)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc -c $< -o $@

$(TEST_BUILD)/test_cell_filter: tests/test_cell_filter.cpp $(BUILD_DIR)/cell_filter.o | $(TEST_BUILD)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc $^ -o $@

$(BUILD_DIR)/lowess.o: src/lowess.cpp src/lowess.hpp | $(BUILD_DIR)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc -c $< -o $@

$(TEST_BUILD)/test_lowess: tests/test_lowess.cpp $(BUILD_DIR)/lowess.o | $(TEST_BUILD)
	$(CXX) -std=c++17 -Wall -Wextra -O2 -Isrc $^ -o $@

$(CBS_BIN): $(BUILD_DIR)/cbs_main.o $(BUILD_DIR)/cbs_args.o $(BUILD_DIR)/bin_gc.o $(BUILD_DIR)/raw_counts_loader.o $(BUILD_DIR)/cell_filter.o $(BUILD_DIR)/lowess.o
	$(CXX) $^ -o $@ -lz

$(TEST_BUILD)/test_cli_args: tests/test_cli_args.cpp $(BUILD_DIR)/cli_args.o | $(TEST_BUILD)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

$(BUILD_DIR)/bins.o: src/bins.cpp src/bins.hpp src/cli_args.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_BUILD)/test_bins: tests/test_bins.cpp $(BUILD_DIR)/bins.o $(BUILD_DIR)/cli_args.o | $(TEST_BUILD)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)

TESTS := $(TEST_BUILD)/test_cbs_args $(TEST_BUILD)/test_cli_args $(TEST_BUILD)/test_bins $(TEST_BUILD)/test_barcode_index $(TEST_BUILD)/test_counts_matrix $(TEST_BUILD)/test_read_filter $(TEST_BUILD)/test_fragment_pairing $(TEST_BUILD)/test_discordant_writer $(TEST_BUILD)/test_bin_gc $(TEST_BUILD)/test_raw_counts_loader $(TEST_BUILD)/test_cell_filter $(TEST_BUILD)/test_lowess

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
	rm -rf $(BUILD_DIR) $(TEST_BUILD) $(BIN) $(CBS_BIN)
