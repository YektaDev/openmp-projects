#include <omp.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <num_threads> <output_file>\n";
        return 1;
    }

    int num_threads = std::stoi(argv[1]);
    std::string output_file = argv[2];

    omp_set_num_threads(num_threads);

    double start_time = omp_get_wtime();

#pragma omp parallel default(none)
    {
        for (int i = 0; i < 1000000; ++i) {
            double temp = i * 0.5;
        }
        // std::cout << "Hello from thread " << omp_get_thread_num() << " of " << omp_get_num_threads() << "\n";
    }

    double end_time = omp_get_wtime();
    double execution_time = end_time - start_time;

    std::ofstream outfile(output_file, std::ios::app);
    outfile << num_threads << " " << std::fixed << std::setprecision(6) << execution_time << "\n";
    outfile.close();

    return 0;
}
