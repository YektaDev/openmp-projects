#include <omp.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <numeric> // accumulate
#include <algorithm> // min & max

// The Initial KMP Implementation:
// https://github.com/Kumar-laxmi/Algorithms/blob/main/C%2B%2B/Pattern-Matching/KMP_Algorithm.cpp

using namespace std;

const string txt =
        "qwertyuiopasdfghjklzxcvbnmqwertyuiopasdfghjklzxcvbnmqwertyuiopasdfghjklzxcvbnmqwertyuiopasdfghjklzxcvbnm________";
const string longTxt = txt + txt + txt + txt + txt + txt + txt;
const string veryLongText = longTxt + longTxt + longTxt + longTxt + longTxt + longTxt + longTxt + longTxt + longTxt;
const string pat =
        "qwertyuiopasdfghjklzxcvbnmqwertyuiopasdfghjklzxcvbnmqwertyuiopasdfghjklzxcvbnmqwertyuiopasdfghjklzxcvbnm";

// Fills lps array for given pattern pat[0..m-1]
void computeLpsArray(const string &pat, vector<int> &lps, const int m) {
    // length of the previous longest prefix suffix
    int len = 0;

    lps[0] = 0; // lps[0] is always 0

    // The loop calculates lps[i] for i = 1 to m-1
    int i = 1;
    while (i < m) {
        if (pat[i] == pat[len]) {
            len++;
            lps[i] = len;
            i++;
        } else {
            // pat[i] != pat[len]
            // The idea is similar to the search step.
            if (len != 0) {
                len = lps[len - 1];
                // Also, note that we do not increment i here
            } else {
                // len == 0
                lps[i] = 0;
                i++;
            }
        }
    }
}

// Finds occurrences of pat in txt
vector<int> kmpSearch(const string &txt, const string &pat) {
    const int m = pat.size();
    const int n = txt.size();

    // Create lps array that will hold the longest prefix suffix values for pattern
    vector<int> lps(m);

    // Preprocess the pattern (calculate lps array)
    computeLpsArray(pat, lps, m);

    vector<int> all_occurrences;

    #pragma omp parallel
    {
        vector<int> local_occurrences;
        const int num_threads = omp_get_num_threads();
        const int thread_id = omp_get_thread_num();

        // Calculate chunk size and overlap
        const int chunk_size = n / num_threads;
        const int overlap = m - 1;

        // Calculate start and end indices for each thread
        int start = thread_id * chunk_size;
        int end = (thread_id == num_threads - 1) ? n : start + chunk_size + overlap;

        // Adjust start for threads other than the first to avoid missing occurrences at chunk boundaries
        if (thread_id != 0) {
            start -= overlap;
        }

        int i = start; // index for txt
        int j = 0;     // index for pat

        while (i < end && n - i >= m - j) {
            if (pat[j] == txt[i]) {
                j++;
                i++;
            }
            if (j == m) {
                const int occurrence_index = i - j;
                // Check if the occurrence is within the thread's main chunk (not just the overlap)
                if (occurrence_index >= thread_id * chunk_size && occurrence_index < (thread_id + 1) * chunk_size) {
                    local_occurrences.push_back(occurrence_index);
                }
                j = lps[j - 1];
            }
            // Mismatch after j matches
            else if (i < end && pat[j] != txt[i]) {
                // Do not match lps[0..lps[j-1]] characters, they will match anyway
                if (j != 0) j = lps[j - 1];
                else i = i + 1;
            }
        }

        // Merge local occurrences into the global vector
        #pragma omp critical
        {
            all_occurrences.insert(all_occurrences.end(), local_occurrences.begin(), local_occurrences.end());
        }
    }

    sort(all_occurrences.begin(), all_occurrences.end());
    return all_occurrences;
}

void program() {
    const vector<int> occurrences = kmpSearch(longTxt, pat);
    cout << "Found pattern at indices: ";
    for (int i = 0; i < occurrences.size(); ++i) {
        cout << occurrences[i] << (i == occurrences.size() - 1 ? "" : ", ");
    }
    cout << endl;
}

int measure(const int num_threads, const string &output_file) {
    omp_set_num_threads(num_threads);

    constexpr int warmups = 3;
    constexpr int runs = 8;

    vector<double> exec_times;

    // Warm-up
    for (int i = 0; i < warmups; i++) {
        cout << "Warm-up Round " << i + 1 << "/" << warmups << " with " << num_threads << " threads\n";
        program();
    }

    // Measurement
    for (int i = 0; i < runs; i++) {
        cout << "Round " << i + 1 << "/" << runs << " with " << num_threads << " threads\n";
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
        const int num_threads = stoi(argv[1]);
        const string output_file = argv[2];
        return measure(num_threads, output_file);
    } catch (const exception &e) {
        cerr << "Error: Invalid input arguments: " << e.what() << endl;
        return 1;
    }
}
