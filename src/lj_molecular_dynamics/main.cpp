/* ---------------------------------------------------------------------------------------------------------------------
 * MD.c - a simple molecular dynamics program for simulating real gas properties of Lennard-Jones particles.
 *
 *   Copyright (C) 2016  Jonathan J. Foley IV, Chelsea Sweet, Oyewumi Akinfenwa
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 *   Electronic Contact:  foleyj10@wpunj.edu
 *   Mail Contact:   Prof. Jonathan Foley
 *                   Department of Chemistry, William Paterson University
 *                   300 Pompton Road
 *                   Wayne NJ 07470
 * ---------------------------------------------------------------------------------------------------------------------
 * REBORN, Refactored, and Parallelized Using OpenMP - Ali Khaleqi Yekta - 2025 - Disclaimer: I don't like this code.
 * ---------------------------------------------------------------------------------------------------------------------
 */

#include <omp.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <numeric> // accumulate
#include <algorithm> // min & max
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

using namespace std;

constexpr int NUM_ITERATIONS = 1000;

// Number of particles
constexpr int N = 216;

// Lennard-Jones parameters in natural units!
constexpr double sigma = 1.;
constexpr double epsilon = 1.;
constexpr double m = 1.;
constexpr double kB = 1.;

constexpr double NA = 6.022140857e23;
constexpr double kBSI = 1.38064852e-23; // m^2*kg/(s^2*K)

// Size of box, which will be specified in natural units
double L;

// Initial Temperature in Natural Units
constexpr double Tinit = 0.728; // Example: 0.728 in natural units corresponds to 103.45 K for Argon
// Vectors!
//
constexpr int MAXPART = 5001;
// Position
double r[MAXPART][3];
// Velocity
double v[MAXPART][3];
// Acceleration
double a[MAXPART][3];
// Force
double F[MAXPART][3];

// atom type
constexpr char atype[10] = "Ar";

// Function prototypes
// initialize positions on simple cubic lattice, also calls function to initialize velocities
void initialize();

// update positions and velocities using Velocity Verlet algorithm
// print particle coordinates to file for rendering via VMD or other animation software
// return 'instantaneous pressure'
double VelocityVerlet(double dt, int iter);

// Compute Force using F = -dV/dr
// solve F = ma for use in Velocity Verlet
void computeAccelerations();

// Numerical Recipes function for generation gaussian distribution
double gaussdist();

// Initialize velocities according to user-supplied initial Temperature (Tinit)
void initializeVelocities();

// Compute total potential energy from particle coordinates
double Potential();

// Compute mean squared velocity from particle velocities
double MeanSquaredVelocity();

// Compute total kinetic energy from particle mass and velocities
double Kinetic();

void program(int NumTime) {
    // variable delcarations
    int i;
    double dt, Vol, Temp, Press, Pavg, Tavg, rho;
    double VolFac, TempFac, PressFac, timefac;
    double KE, PE, mvs, gc, Z;
    char prefix[1000], tfn[1000], ofn[1000], afn[1000];

    // Define constants based on the chosen atom type
    if (strcmp(atype, "He") == 0) {
        VolFac = 1.8399744000000005e-29;
        PressFac = 8152287.336171632;
        TempFac = 10.864459551225972;
        timefac = 1.7572698825166272e-12;
    } else if (strcmp(atype, "Ne") == 0) {
        VolFac = 2.0570823999999997e-29;
        PressFac = 27223022.27659913;
        TempFac = 40.560648991243625;
        timefac = 2.1192341945685407e-12;
    } else if (strcmp(atype, "Ar") == 0) {
        VolFac = 3.7949992920124995e-29;
        PressFac = 51695201.06691862;
        TempFac = 142.0950000000000;
        timefac = 2.09618e-12;
    } else if (strcmp(atype, "Kr") == 0) {
        VolFac = 4.5882712000000004e-29;
        PressFac = 59935428.40275003;
        TempFac = 199.1817584391428;
        timefac = 8.051563913585078e-13;
    } else if (strcmp(atype, "Xe") == 0) {
        VolFac = 5.4872e-29;
        PressFac = 70527773.72794868;
        TempFac = 280.30305642163006;
        timefac = 9.018957925790732e-13;
    } else {
        // Default to Argon
        VolFac = 3.7949992920124995e-29;
        PressFac = 51695201.06691862;
        TempFac = 142.0950000000000;
        timefac = 2.09618e-12;
    }

    // Set simulation parameters
    strcpy(prefix, "argon_md");
    strcpy(tfn, prefix);
    strcat(tfn, "_traj.xyz");
    strcpy(ofn, prefix);
    strcat(ofn, "_output.txt");
    strcpy(afn, prefix);
    strcat(afn, "_average.txt");

    rho = 40; // Example: 40 moles/m^3
    Vol = N / (rho * NA);
    Vol /= VolFac;
    L = pow(Vol, (1. / 3));

    if (strcmp(atype, "He") == 0) {
        // dt in natural units of time s.t. in SI it is 5 f.s. for all other gasses
        dt = 0.2e-14 / timefac;
    } else {
        dt = 0.5e-14 / timefac;
    }

    // Put all the atoms in simple crystal lattice and give them random velocities
    // that corresponds to the initial temperature we have specified
    initialize();

    // Based on their positions, calculate the ininial intermolecular forces
    // The accellerations of each particle will be defined from the forces and their
    // mass, and this will allow us to update their positions via Newton's law
    computeAccelerations();

    // We want to calculate the average Temperature and Pressure for the simulation
    // The variables need to be set to zero initially
    Pavg = 0;
    Tavg = 0;

    int tenp = floor(NumTime / 10);
    printf("  PERCENTAGE OF CALCULATION COMPLETE:\n  [");
    for (i = 0; i < NumTime + 1; i++) {
        // This just prints updates on progress of the calculation for the users convenience
        if (i == tenp) printf(" 10 |");
        else if (i == 2 * tenp) printf(" 20 |");
        else if (i == 3 * tenp) printf(" 30 |");
        else if (i == 4 * tenp) printf(" 40 |");
        else if (i == 5 * tenp) printf(" 50 |");
        else if (i == 6 * tenp) printf(" 60 |");
        else if (i == 7 * tenp) printf(" 70 |");
        else if (i == 8 * tenp) printf(" 80 |");
        else if (i == 9 * tenp) printf(" 90 |");
        else if (i == 10 * tenp) printf(" 100 ]\n");
        fflush(stdout);

        // This updates the positions and velocities using Newton's Laws
        // Also computes the Pressure as the sum of momentum changes from wall collisions / timestep
        // which is a Kinetic Theory of gasses concept of Pressure
        Press = VelocityVerlet(dt, i + 1);
        Press *= PressFac;

        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // Now we would like to calculate somethings about the system:
        // Instantaneous mean velocity squared, Temperature, Pressure
        // Potential, and Kinetic Energy
        // We would also like to use the IGL to try to see if we can extract the gas constant
        mvs = MeanSquaredVelocity();
        KE = Kinetic();
        PE = Potential();

        // Temperature from Kinetic Theory
        Temp = m * mvs / (3 * kB) * TempFac;

        // Instantaneous gas constant and compressibility - not well defined because
        // pressure may be zero in some instances because there will be zero wall collisions,
        // pressure may be very high in some instances because there will be a number of collisions
        gc = NA * Press * (Vol * VolFac) / (N * Temp);
        Z = Press * (Vol * VolFac) / (N * kBSI * Temp);

        Tavg += Temp;
        Pavg += Press;
    }
}

void initialize() {
    int n, p, i, j, k;
    double pos;

    // Number of atoms in each direction
    n = int(ceil(pow(N, 1.0 / 3)));

    // spacing between atoms along a given direction
    pos = L / n;

    // index for number of particles assigned positions
    p = 0;
    // initialize positions
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                if (p < N) {
                    r[p][0] = (i + 0.5) * pos;
                    r[p][1] = (j + 0.5) * pos;
                    r[p][2] = (k + 0.5) * pos;
                }
                p++;
            }
        }
    }

    // Call function to initialize velocities
    initializeVelocities();

    /***********************************************
  *   Uncomment if you want to see what the initial positions and velocities are
    printf("  Printing initial positions!\n");
    for (i=0; i<N; i++) {
      printf("  %6.3e  %6.3e  %6.3e\n",r[i][0],r[i][1],r[i][2]);
    }

    printf("  Printing initial velocities!\n");
    for (i=0; i<N; i++) {
      printf("  %6.3e  %6.3e  %6.3e\n",v[i][0],v[i][1],v[i][2]);
    }
    */
}

// Function to calculate the averaged velocity squared
double MeanSquaredVelocity() {
    double vx2 = 0;
    double vy2 = 0;
    double vz2 = 0;
    double v2 = 0;
    int i;

#pragma omp parallel for private(i) reduction(+:vx2,vy2,vz2)
    for (i = 0; i < N; i++) {
        vx2 = vx2 + v[i][0] * v[i][0];
        vy2 = vy2 + v[i][1] * v[i][1];
        vz2 = vz2 + v[i][2] * v[i][2];
    }
    v2 = (vx2 + vy2 + vz2) / N;

    //printf("  Average of x-component of velocity squared is %f\n",v2);
    return v2;
}

// Function to calculate the kinetic energy of the system
double Kinetic() {
    //Write Function here!

    double v2_temp, kin = 0.;
    int i, j;

#pragma omp parallel for private(i,j,v2_temp) reduction(+:kin)
    for (i = 0; i < N; i++) {
        v2_temp = 0.;
        for (j = 0; j < 3; j++) {
            v2_temp += v[i][j] * v[i][j];
        }
        kin += m * v2_temp / 2.;
    }

    //printf("  Total Kinetic Energy is %f\n",N*mvs*m/2.);
    return kin;
}

// Function to calculate the potential energy of the system
double Potential() {
    double quot, r2, rnorm, term1, term2, Pot = 0.;
    int i, j, k;

#pragma omp parallel for private(i,j,k,r2,rnorm,quot,term1,term2) reduction(+:Pot)
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            if (j != i) {
                r2 = 0.;
                for (k = 0; k < 3; k++) {
                    r2 += (r[i][k] - r[j][k]) * (r[i][k] - r[j][k]);
                }
                rnorm = sqrt(r2);
                quot = sigma / rnorm;
                term1 = pow(quot, 12.);
                term2 = pow(quot, 6.);

                Pot += 4 * epsilon * (term1 - term2);
            }
        }
    }

    return Pot;
}

// Uses the derivative of the Lennard-Jones potential to calculate
// the forces on each atom. Then uses a = F/m to calculate the
// accelleration of each atom.
void computeAccelerations() {
    int i, j, k;
    double f, rSqd;
    double rij[3]; // position of i relative to j

#pragma omp parallel for private(i,k)
    for (i = 0; i < N; i++) {
        // set all accelerations to zero
        for (k = 0; k < 3; k++) {
            a[i][k] = 0;
        }
    }
#pragma omp parallel for private(i,j,k,rSqd,rij,f)
    for (i = 0; i < N - 1; i++) {
        // loop over all distinct pairs i,j
        for (j = i + 1; j < N; j++) {
            // initialize r^2 to zero
            rSqd = 0;

            for (k = 0; k < 3; k++) {
                // component-by-componenent position of i relative to j
                rij[k] = r[i][k] - r[j][k];
                // sum of squares of the components
                rSqd += rij[k] * rij[k];
            }

            // From derivative of Lennard-Jones with sigma and epsilon set equal to 1 in natural units!
            f = 24 * (2 * pow(rSqd, -7) - pow(rSqd, -4));
            for (k = 0; k < 3; k++) {
                // from F = ma, where m = 1 in natural units!
#pragma omp atomic
                a[i][k] += rij[k] * f;
#pragma omp atomic
                a[j][k] -= rij[k] * f;
            }
        }
    }
}

// returns sum of dv/dt*m/A (aka Pressure) from elastic collisions with walls
double VelocityVerlet(double dt, int iter) {
    int i, j, k;

    double psum = 0.;

    // Compute accelerations from forces at current position
    computeAccelerations();
    // Update positions and velocity with current velocity and acceleration
    //printf("  Updated Positions!\n");
#pragma omp parallel for private(i,j)
    for (i = 0; i < N; i++) {
        for (j = 0; j < 3; j++) {
            r[i][j] += v[i][j] * dt + 0.5 * a[i][j] * dt * dt;

            v[i][j] += 0.5 * a[i][j] * dt;
        }
        //printf("  %i  %6.4e   %6.4e   %6.4e\n",i,r[i][0],r[i][1],r[i][2]);
    }
    // Update accellerations from updated positions
    computeAccelerations();
    // Update velocity with updated acceleration
#pragma omp parallel for private(i,j)
    for (i = 0; i < N; i++) {
        for (j = 0; j < 3; j++) {
            v[i][j] += 0.5 * a[i][j] * dt;
        }
    }

    // Elastic walls
#pragma omp parallel for private(i,j) reduction(+:psum)
    for (i = 0; i < N; i++) {
        for (j = 0; j < 3; j++) {
            if (r[i][j] < 0.) {
                v[i][j] *= -1.; //- elastic walls
                psum += 2 * m * fabs(v[i][j]) / dt; // contribution to pressure from "left" walls
            }
            if (r[i][j] >= L) {
                v[i][j] *= -1.; //- elastic walls
                psum += 2 * m * fabs(v[i][j]) / dt; // contribution to pressure from "right" walls
            }
        }
    }

    return psum / (6 * L * L);
}

void initializeVelocities() {
    int i, j;

#pragma omp parallel for private(i,j)
    for (i = 0; i < N; i++) {
        for (j = 0; j < 3; j++) {
            // Pull a number from a Gaussian Distribution
            v[i][j] = gaussdist();
        }
    }

    // Vcm = sum_i^N  m*v_i/  sum_i^N  M
    // Compute center-of-mas velocity according to the formula above
    double vCM[3] = {0, 0, 0};

#pragma omp parallel for private(i,j) reduction(+:vCM[:3])
    for (i = 0; i < N; i++) {
        for (j = 0; j < 3; j++) {
            vCM[j] += m * v[i][j];
        }
    }

    for (i = 0; i < 3; i++) vCM[i] /= N * m;

    // Subtract out the center-of-mass velocity from the
    // velocity of each particle... effectively set the
    // center of mass velocity to zero so that the system does
    // not drift in space!
#pragma omp parallel for private(i,j)
    for (i = 0; i < N; i++) {
        for (j = 0; j < 3; j++) {
            v[i][j] -= vCM[j];
        }
    }

    // Now we want to scale the average velocity of the system
    // by a factor which is consistent with our initial temperature, Tinit
    double vSqdSum = 0., lambda;
#pragma omp parallel for private(i,j) reduction(+:vSqdSum)
    for (i = 0; i < N; i++) {
        for (j = 0; j < 3; j++) {
            vSqdSum += v[i][j] * v[i][j];
        }
    }

    lambda = sqrt(3 * (N - 1) * Tinit / vSqdSum);

#pragma omp parallel for private(i,j)
    for (i = 0; i < N; i++) {
        for (j = 0; j < 3; j++) {
            v[i][j] *= lambda;
        }
    }
}

// Numerical recipes Gaussian distribution number generator
double gaussdist() {
    static bool available = false;
    static double gset;
    double fac, rsq, v1, v2;
    if (!available) {
        do {
            v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);

        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        available = true;

        return v2 * fac;
    } else {
        available = false;
        return gset;
    }
}

int measure(const int num_threads, const string &output_file) {
    omp_set_num_threads(num_threads);

    constexpr int warmups = 3;
    constexpr int runs = 8;

    vector<double> exec_times;

    // Warm-up
    for (int i = 0; i < warmups; i++) {
        cout << "Warm-up Round " << i + 1 << "/" << warmups << " with " << num_threads << " threads\n";
        program(NUM_ITERATIONS);
    }

    // Measurement
    for (int i = 0; i < runs; i++) {
        cout << "Round " << i + 1 << "/" << runs << " with " << num_threads << " threads\n";
        const double start_time = omp_get_wtime();
        program(NUM_ITERATIONS);
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

int main(int argc, char *argv[]) {
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
