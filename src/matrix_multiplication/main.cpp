#include <omp.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <numeric> // accumulate
#include <algorithm> // min & max

using namespace std;

// Function to multiply two matrices
vector<vector<int> > multiplyMatrices(const vector<vector<int> > &matrixA, const vector<vector<int> > &matrixB) {
    const int rowsA = matrixA.size();
    const int colsA = matrixA[0].size();
    const int rowsB = matrixB.size();
    const int colsB = matrixB[0].size();

    // Check for compatibility
    if (colsA != rowsB) {
        cout << "Error: Matrices are incompatible for multiplication." << endl;
        return {}; // Return an empty matrix to indicate an error
    }

    // Initialize the result matrix with zeros
    vector<vector<int> > resultMatrix(rowsA, vector<int>(colsB, 0));

    // Perform the matrix multiplication
    for (int i = 0; i < rowsA; ++i) {
        for (int j = 0; j < colsB; ++j) {
            for (int k = 0; k < colsA; ++k) {
                resultMatrix[i][j] += matrixA[i][k] * matrixB[k][j];
            }
        }
    }

    return resultMatrix;
}

void program() {
    vector<vector<int> > matrixA = {{1, 2, 3}, {4, 5, 6}};
    vector<vector<int> > matrixB = {{7, 8}, {9, 10}, {11, 12}};
    vector<vector<int> > result = multiplyMatrices(matrixA, matrixB);

    // Print the result
    if (!result.empty()) {
        cout << "Resultant Matrix:" << endl;
        for (const auto &row: result) {
            for (int val: row) {
                cout << val << " ";
            }
            cout << endl;
        }
    }
}

int measure(const int num_threads, const string &output_file) {
    omp_set_num_threads(num_threads);

    constexpr int warmups = 5;
    constexpr int runs = 10;

    vector<double> exec_times;

    // Warm-up
    for (int i = 0; i < warmups; i++) {
        program();
    }

    // Measurement
    for (int i = 0; i < runs; i++) {
        const double start_time = omp_get_wtime();
        program();
        const double end_time = omp_get_wtime();
        exec_times.push_back(end_time - start_time);
    }

    // Calculate statistics
    const double sum = accumulate(exec_times.begin(), exec_times.end(), 0.0);
    const double average = sum / exec_times.size();
    const double min_time = *ranges::min_element(exec_times);
    const double max_time = *ranges::max_element(exec_times);
    const double sq_sum = inner_product(exec_times.begin(), exec_times.end(), exec_times.begin(), 0.0);
    const double stdev = sqrt(sq_sum / exec_times.size() - average * average);

    // Write output to file
    if (ofstream outfile(output_file, ios::app); outfile.is_open()) {
        outfile << num_threads << " "
                << fixed << setprecision(6) << average << " "
                << fixed << setprecision(6) << min_time << " "
                << fixed << setprecision(6) << max_time << " "
                << fixed << setprecision(6) << stdev << "\n";
        outfile.close();
    } else {
        cerr << "Error: Could not open output file " << output_file << endl;
        return 1;
    }
    return 0;
}

int main(const int argc, char *argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <num_threads> <output_file>\n";
        return 1;
    }

    try {
        const int num_threads = std::stoi(argv[1]);
        const string output_file = argv[2];
        return measure(num_threads, output_file);
    } catch (const std::exception &e) {
        cerr << "Error: Invalid input arguments: " << e.what() << endl;
        return 1;
    }
}
