#include <omp.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <numeric> // accumulate
#include <algorithm> // min & max
#include <random> // mt19937

using namespace std;

vector<vector<int> > matrixA;
vector<vector<int> > matrixB;

vector<vector<int> > multiplyMatrices(const vector<vector<int> > &matrixA, const vector<vector<int> > &matrixB) {
    const int rowsA = matrixA.size();
    const int colsA = matrixA[0].size();
    const int rowsB = matrixB.size();
    const int colsB = matrixB[0].size();
    if (colsA != rowsB) {
        throw runtime_error("Error: Matrices are incompatible for multiplication. "
            "Columns of the first matrix must equal rows of the second.");
    }

    vector<vector<int> > resultMatrix(rowsA, vector<int>(colsB, 0));

    //   - Parallelize the outer loop using OpenMP.
    //   - The 'i' loop iterates over the rows of 'resultMatrix'.
    //   - Each thread will handle a different set of rows.
    //   - `schedule(static)` divides the iterations among threads in a static manner.
    //   - `private(j, k)` ensures that each thread has its own private copies of the loop indices
    //      'j' and 'k', preventing race conditions.
#pragma omp parallel for schedule(static)
    for (int i = 0; i < rowsA; i++) {
        for (int j = 0; j < colsB; j++) {
            for (int k = 0; k < colsA; k++) {
                resultMatrix[i][j] += matrixA[i][k] * matrixB[k][j];
            }
        }
    }

    return resultMatrix;
}

vector<vector<int> > generateMatrix(int rows, int cols, int seed) {
    // A random number generator engine with the given seed
    mt19937 generator(seed);
    // A uniform distribution between a suitable range (e.g., 1 to 100)
    uniform_int_distribution<int> distribution(1, 100);
    vector<vector<int> > matrix(rows, vector<int>(cols));
    // Fill the matrix with random numbers
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = distribution(generator);
        }
    }
    return matrix;
}

void program() {
    // Call multiplyMatrices() within a try-catch block to handle potential errors
    try {
        const vector<vector<int> > result = multiplyMatrices(matrixA, matrixB);

        if (const size_t elements = result.size() * result[0].size(); elements < 100) {
            // Print the result if no exception was thrown
            cout << "Resultant Matrix:" << endl;
            for (const auto &row: result) {
                for (const int val: row) {
                    cout << val << " ";
                }
                cout << endl;
            }
        } else {
            cout << "Resultant Matrix Size: " << elements << endl;
        }
    } catch (const runtime_error &error) {
        cerr << error.what() << endl;
    }
}

int measure(const int num_threads, const string &output_file) {
    omp_set_num_threads(num_threads);

    constexpr int warmups = 3;
    constexpr int runs = 8;

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

    matrixA = generateMatrix(2048, 64, 20482048);
    matrixB = generateMatrix(64, 4096, 40964096);

    try {
        const int num_threads = stoi(argv[1]);
        const string output_file = argv[2];
        return measure(num_threads, output_file);
    } catch (const exception &e) {
        cerr << "Error: Invalid input arguments: " << e.what() << endl;
        return 1;
    }
}
