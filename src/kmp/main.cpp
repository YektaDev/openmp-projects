#include <omp.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <numeric> // accumulate
#include <algorithm> // min & max

using namespace std;

const string txt = "qwertyuiopasdfghjklzxcvbnmqwertyuiopasdfghjklzxcvbnmqwertyuiopasdfghjklzxcvbnmqwertyuiopasdfghjklzxcvbnm________";
const string longTxt = txt + txt + txt + txt + txt + txt + txt;
const string pat = "qwertyuiopasdfghjklzxcvbnmqwertyuiopasdfghjklzxcvbnmqwertyuiopasdfghjklzxcvbnmqwertyuiopasdfghjklzxcvbnm";

// Fills lps array for given pattern pat[0..m-1]
void computeLPSArray(const string &pat, vector<int> &lps, const int m) {
    // length of the previous longest prefix suffix
    int len = 0;

    lps[0] = 0; // lps[0] is always 0

    // the loop calculates lps[i] for i = 1 to m-1
    int i = 1;
    while (i < m) {
        if (pat[i] == pat[len]) {
            len++;
            lps[i] = len;
            i++;
        } else // (pat[i] != pat[len])
        {
            //  The idea is similar to search step.
            if (len != 0) {
                len = lps[len - 1];
                // Also, note that we do not increment i here
            } else // if (len == 0)
            {
                lps[i] = 0;
                i++;
            }
        }
    }
}

// Finds occurrences of pat in txt
void KMPSearch(const string &txt, const string &pat) {
    const int m = pat.size();
    const int n = txt.size();

    // create lps array that will hold the longest prefix suffix values for pattern
    vector<int> lps(m);

    // Preprocess the pattern (calculate lps array)
    computeLPSArray(pat, lps, m);

    int i = 0; // index for txt
    int j = 0; // index for pat

    while (n - i >= m - j) {
        if (pat[j] == txt[i]) {
            j++;
            i++;
        }
        if (j == m) {
            cout << "Found pattern at index " << i - j << endl;
            j = lps[j - 1];
        }
        // mismatch after j matches
        else if (i < n && pat[j] != txt[i]) {
            // Do not match lps[0..lps[j-1]] characters, they will match anyway
            if (j != 0) {
                j = lps[j - 1];
            } else {
                i = i + 1;
            }
        }
    }
}

void program() {
    KMPSearch(longTxt, pat);
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
