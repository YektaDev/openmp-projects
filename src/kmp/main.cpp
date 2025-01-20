#include <omp.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <numeric> // std::accumulate
#include <algorithm> // std::min & std::max

using namespace std;

void program() {
    #pragma omp parallel default(none)
    {
        // Your parallel code here...
        // Make sure to use appropriate clauses (shared, private, firstprivate)
        // Example (assuming you have a loop to parallelize):
        #pragma omp for
        for (int i = 0; i < 1000000; ++i) {
            // ... do some work ...
        }
    }
}

int measure(const char *const argv[]) {
    const int num_threads = std::stoi(argv[1]);
    const std::string output_file = argv[2];

    omp_set_num_threads(num_threads);

    constexpr int warmups = 5;
    constexpr int runs = 10;

    std::vector<double> execution_times;

    // Warm-up
    for (int i = 0; i < warmups; ++i) {
        program();
    }

    // Measurement
    for (int i = 0; i < runs; ++i) {
        const double start_time = omp_get_wtime();
        program();
        const double end_time = omp_get_wtime();
        execution_times.push_back(end_time - start_time);
    }

    // Calculate statistics
    const double sum = std::accumulate(execution_times.begin(), execution_times.end(), 0.0);
    const double average = sum / execution_times.size();
    const double min_time = *ranges::min_element(execution_times);
    const double max_time = *ranges::max_element(execution_times);
    const double sq_sum = std::inner_product(execution_times.begin(), execution_times.end(), execution_times.begin(), 0.0);
    const double stdev = std::sqrt(sq_sum / execution_times.size() - average * average);

    // Write output to file
    if (std::ofstream outfile(output_file, std::ios::app); outfile.is_open()) {
        outfile << num_threads << " "
                << std::fixed << std::setprecision(6) << average << " "
                << std::fixed << std::setprecision(6) << min_time << " "
                << std::fixed << std::setprecision(6) << max_time << " "
                << std::fixed << std::setprecision(6) << stdev << "\n";
        outfile.close();
    } else {
        std::cerr << "Error: Could not open output file " << output_file << std::endl;
        return 1;
    }
    return 0;
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <num_threads> <output_file>\n";
        return 1;
    }
    return measure(argv);
}
