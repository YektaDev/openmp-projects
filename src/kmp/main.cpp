#include <omp.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <numeric> // std::accumulate
#include <algorithm> // std::min & std::max

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

int measure(char* argv[]) {
    int num_threads = std::stoi(argv[1]);
    std::string output_file = argv[2];

    omp_set_num_threads(num_threads);

    const int warmups = 5; // Number of warm-up runs
    const int runs = 10; // Number of measurement runs

    std::vector<double> execution_times;

    // Warm-up phase
    for (int i = 0; i < warmups; ++i) {
        program();
    }

    // Measurement phase
    for (int i = 0; i < runs; ++i) {
        double start_time = omp_get_wtime();
        program();
        double end_time = omp_get_wtime();
        execution_times.push_back(end_time - start_time);
    }

    // Calculate statistics
    double sum = std::accumulate(execution_times.begin(), execution_times.end(), 0.0);
    double average = sum / execution_times.size();
    double min_time = *std::min_element(execution_times.begin(), execution_times.end());
    double max_time = *std::max_element(execution_times.begin(), execution_times.end());
    double sq_sum = std::inner_product(execution_times.begin(), execution_times.end(), execution_times.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / execution_times.size() - average * average);

    // Write the output to the specified file
    std::ofstream outfile(output_file, std::ios::app);
    if (outfile.is_open()) {
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

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <num_threads> <output_file>\n";
        return 1;
    }
    return measure(argv);
}