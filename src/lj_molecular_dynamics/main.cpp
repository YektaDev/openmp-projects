#include <omp.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <numeric> // accumulate
#include <algorithm> // min & max
#include <cmath>

using namespace std;

// Structure to represent a particle
struct Particle {
    double x, y; // Position
    double vx, vy; // Velocity
    double fx, fy; // Force
};

// Function to calculate Lennard-Jones potential
double lennardJonesPotential(const double r, const double epsilon, const double sigma) {
    const double r6 = pow(sigma / r, 6);
    const double r12 = r6 * r6;
    return 4 * epsilon * (r12 - r6);
}

// Function to calculate distance between two particles
double calculateDistance(const Particle &p1, const Particle &p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

// Function to calculate forces between particles
void calculateForces(vector<Particle> &particles, double epsilon, double sigma) {
    const int num_particles = particles.size();

#pragma omp parallel for
    for (int i = 0; i < num_particles; ++i) {
        particles[i].fx = 0.0;
        particles[i].fy = 0.0;
        for (int j = 0; j < num_particles; ++j) {
            if (i != j) {
                const double dist = calculateDistance(particles[i], particles[j]);
                const double r_inv = 1.0 / dist;
                const double r6_inv = pow(r_inv * sigma, 6);
                const double r12_inv = r6_inv * r6_inv;
                const double force_magnitude = 24 * epsilon * r_inv * (2 * r12_inv - r6_inv);

                const double dx = particles[i].x - particles[j].x;
                const double dy = particles[i].y - particles[j].y;
                particles[i].fx += force_magnitude * dx * r_inv;
                particles[i].fy += force_magnitude * dy * r_inv;
            }
        }
    }
}

// Function to update particle positions and velocities
void updateParticles(vector<Particle> &particles, double dt) {
#pragma omp parallel for
    for (auto &particle: particles) {
        // Update velocity
        particle.vx += particle.fx * dt;
        particle.vy += particle.fy * dt;
        // Update position
        particle.x += particle.vx * dt;
        particle.y += particle.vy * dt;
    }
}

void program() {
    // Simulation parameters
    constexpr int num_particles = 100;
    constexpr double epsilon = 1.0;
    constexpr double sigma = 1.0;
    constexpr double dt = 0.001;
    constexpr int num_steps = 2000;

    // Initialize particles
    vector<Particle> particles(num_particles);
    for (int i = 0; i < num_particles; ++i) {
        particles[i].x = (double) rand() / RAND_MAX;
        particles[i].y = (double) rand() / RAND_MAX;
        particles[i].vx = 0.0;
        particles[i].vy = 0.0;
        particles[i].fx = 0.0;
        particles[i].fy = 0.0;
    }

    // Simulation loop
    for (int step = 0; step < num_steps; ++step) {
        calculateForces(particles, epsilon, sigma);
        updateParticles(particles, dt);
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

    try {
        const int num_threads = stoi(argv[1]);
        const string output_file = argv[2];
        return measure(num_threads, output_file);
    } catch (const exception &e) {
        cerr << "Error: Invalid input arguments: " << e.what() << endl;
        return 1;
    }
}
