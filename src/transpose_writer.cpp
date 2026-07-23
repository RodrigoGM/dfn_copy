#include "transpose_writer.hpp"
#include <zlib.h>
#include <sstream>
#include <stdexcept>

void write_transposed_gz(const std::string& path,
                          const std::vector<float>& scratch_data,
                          size_t num_cells,
                          size_t num_bins,
                          const std::vector<std::string>& bin_names,
                          const std::vector<std::string>& barcode_names) {
    if (bin_names.size() != num_bins) {
        throw std::runtime_error("transpose writer: bin_names size " +
                                  std::to_string(bin_names.size()) +
                                  " does not match num_bins " + std::to_string(num_bins));
    }
    if (barcode_names.size() != num_cells) {
        throw std::runtime_error("transpose writer: barcode_names size " +
                                  std::to_string(barcode_names.size()) +
                                  " does not match num_cells " + std::to_string(num_cells));
    }
    if (scratch_data.size() != num_cells * num_bins) {
        throw std::runtime_error("transpose writer: scratch_data size " +
                                  std::to_string(scratch_data.size()) +
                                  " does not match num_cells * num_bins (" +
                                  std::to_string(num_cells * num_bins) + ")");
    }

    gzFile out = gzopen(path.c_str(), "wb");
    if (!out) {
        throw std::runtime_error("cannot open output file for writing: " + path);
    }

    std::ostringstream header;
    header << "bin";
    for (const auto& name : barcode_names) header << '\t' << name;
    header << '\n';
    std::string header_str = header.str();
    gzwrite(out, header_str.data(), static_cast<unsigned>(header_str.size()));

    for (size_t bin = 0; bin < num_bins; ++bin) {
        std::ostringstream row;
        row << bin_names[bin];
        for (size_t cell = 0; cell < num_cells; ++cell) {
            row << '\t' << scratch_data[cell * num_bins + bin];
        }
        row << '\n';
        std::string row_str = row.str();
        gzwrite(out, row_str.data(), static_cast<unsigned>(row_str.size()));
    }

    gzclose(out);
}
