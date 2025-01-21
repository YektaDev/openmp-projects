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

using namespace std;

// Number of particles
constexpr int N = 216;

// Lennard-Jones parameters in natural units!
double sigma = 1.;
double epsilon = 1.;
double m = 1.;
double kB = 1.;

double NA = 6.022140857e23;
double kBSI = 1.38064852e-23; // m^2*kg/(s^2*K)

// Size of box, which will be specified in natural units
double L;

// Initial Temperature in Natural Units
double Tinit; //2;
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
char atype[10];

// Function prototypes
// initialize positions on simple cubic lattice, also calls function to initialize velocities
void initialize();

// update positions and velocities using Velocity Verlet algorithm
// print particle coordinates to file for rendering via VMD or other animation software
// return 'instantaneous pressure'
double VelocityVerlet(double dt);

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

// Forward declaration of the function
int write_output(double Pavg, double Tavg, double Z, double gc, double Vol, int N, double timefac, double dt,
                 const string &ofn);

void program() {
    char prefix[1000], tfn[1000], ofn[1000], afn[1000];

    strcpy(prefix, "md");
    strcpy(tfn, prefix);
    strcat(tfn, "_traj.xyz");
    strcpy(ofn, prefix);
    strcat(ofn, "_output.txt");
    strcpy(afn, prefix);
    strcat(afn, "_average.txt");

    constexpr double VolFac = 3.7949992920124995e-29;
    constexpr double PressFac = 51695201.06691862;
    constexpr double TempFac = 142.0950000000000;
    constexpr double timefac = 2.09618e-12;
    strcpy(atype, "Ar");

    // Initial temprature of the gas in Kelvin
    Tinit = 300;
    // Convert initial temperature from kelvin to natural units
    Tinit /= TempFac;

    // The number of density in moles.
    // Ideal gas: 40 moles/m^3
    // Liquid Argon At 1ATM and 87K: ~35000 moles/m^3
    constexpr double rho = 40;
    const double Vol = (N / (rho * NA)) / VolFac;

    // Limiting N to MAXPART for practical reasons
    if constexpr (N >= MAXPART) {
        printf("\n\n  MAXIMUM NUMBER OF PARTICLES IS %i\n\n  PLEASE ADJUST YOUR INPUT FILE ACCORDINGLY \n\n", MAXPART);
        exit(0);
    }

    // Check to see if the volume makes sense - is it too small?
    // Remember VDW radius of the particles is 1 natural unit of length
    // and volume = L*L*L, so if V = N*L*L*L = N, then all the particles
    // will be initialized with an interparticle separation equal to 2xVDW radius
    if (Vol < N) {
        printf("\n\n\n  YOUR DENSITY IS VERY HIGH!\n\n");
        printf("  THE NUMBER OF PARTICLES IS %i AND THE AVAILABLE VOLUME IS %f NATURAL UNITS\n", N, Vol);
        printf("  SIMULATIONS WITH DENSITY GREATER THAN 1 PARTCICLE/(1 Natural Unit of Volume) MAY DIVERGE\n");
        printf("  PLEASE ADJUST YOUR INPUT FILE ACCORDINGLY AND RETRY\n\n");
        exit(0);
    }
    // Vol = L*L*L;
    // Length of the box in natural units:
    L = pow(Vol, (1. / 3));

    // dt in natural units of time s.t. in SI it is 5 f.s. for all other gasses
    constexpr double dt = 0.5e-14 / timefac;
    // We will run the simulation for NumTime timesteps.
    // The total time will be NumTime*dt in natural units
    // And NumTime*dt multiplied by the appropriate conversion factor for time in seconds
    constexpr int NumTime = 2000;

    // Put all the atoms in simple crystal lattice and give them random velocities
    // that corresponds to the initial temperature we have specified
    initialize();

    // Based on their positions, calculate the ininial intermolecular forces
    // The accellerations of each particle will be defined from the forces and their
    // mass, and this will allow us to update their positions via Newton's law
    computeAccelerations();

    double Pavg = 0;
    double Tavg = 0;

    // Open ofn here
    ofstream ofp(ofn);
    if (!ofp.is_open()) {
        cerr << "Error: Could not open output file " << ofn << endl;
        return;
    }
    ofp <<
            "   time (s)            T(t) (K)             P(t) (Pa)           Kinetic En. (n.u.)     Potential En. (n.u.) Total En. (n.u.)\n";

    for (int i = 0; i < NumTime + 1; i++) {
        double Press = VelocityVerlet(dt);
        Press *= PressFac;

        // Now we would like to calculate somethings about the system:
        // Instantaneous mean velocity squared, Temperature, Pressure
        // Potential, and Kinetic Energy
        // We would also like to use the IGL to try to see if we can extract the gas constant
        const double mvs = MeanSquaredVelocity();
        const double KE = Kinetic();
        const double PE = Potential();

        // Temperature from Kinetic Theory
        const double Temp = m * mvs / (3 * kB) * TempFac;

        Tavg += Temp;
        Pavg += Press;

        // Write instantaneous values to ofp
        ofp << "  " << setw(8) << fixed << setprecision(4) << i * dt * timefac
                << "  " << setw(20) << fixed << setprecision(8) << Temp
                << "  " << setw(20) << fixed << setprecision(8) << Press
                << " " << setw(20) << fixed << setprecision(8) << KE
                << "  " << setw(20) << fixed << setprecision(8) << PE
                << "  " << setw(20) << fixed << setprecision(8) << KE + PE << "\n";
    }
    ofp.close();

    // Because we have calculated the instantaneous temperature and pressure,
    // we can take the average over the whole simulation here
    Pavg /= NumTime;
    Tavg /= NumTime;
    const double Z = Pavg * (Vol * VolFac) / (N * kBSI * Tavg);
    const double gc = NA * Pavg * (Vol * VolFac) / (N * Tavg);
    write_output(Pavg, Tavg, Z, gc, Vol, N, timefac, dt, ofn);
}

void initialize() {
    // Number of atoms in each direction
    const int n = static_cast<int>(ceil(pow(N, 1.0 / 3)));
    // spacing between atoms along a given direction
    const double pos = L / n;
    // index for number of particles assigned positions
    int p = 0;
    // initialize positions
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
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
}

// Function to calculate the averaged velocity squared
double MeanSquaredVelocity() {
    double vx2 = 0;
    double vy2 = 0;
    double vz2 = 0;

#pragma omp parallel for reduction(+:vx2,vy2,vz2)
    for (int i = 0; i < N; i++) {
        vx2 = vx2 + v[i][0] * v[i][0];
        vy2 = vy2 + v[i][1] * v[i][1];
        vz2 = vz2 + v[i][2] * v[i][2];
    }
    // Average of x-component of velocity squared
    return (vx2 + vy2 + vz2) / N;
}

// Function to calculate the kinetic energy of the system
double Kinetic() {
    double kin = 0.;
#pragma omp parallel for reduction(+:kin) private(v2)
    for (int i = 0; i < N; i++) {
        double v2 = 0.;
        for (int j = 0; j < 3; j++) {
            v2 += v[i][j] * v[i][j];
        }
        kin += m * v2 / 2.;
    }

    //printf("  Total Kinetic Energy is %f\n",N*mvs*m/2.);
    return kin;
}

// Function to calculate the potential energy of the system
double Potential() {
    double Pot = 0.;
#pragma omp parallel for reduction(+:Pot) private(r2, rnorm, quot, term1, term2)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (j != i) {
                double r2 = 0.;
                for (int k = 0; k < 3; k++) {
                    r2 += (r[i][k] - r[j][k]) * (r[i][k] - r[j][k]);
                }
                const double rnorm = sqrt(r2);
                const double quot = sigma / rnorm;
                const double term1 = pow(quot, 12.);
                const double term2 = pow(quot, 6.);

                Pot += 4 * epsilon * (term1 - term2);
            }
        }
    }

    return Pot;
}

//  Uses the derivative of the Lennard-Jones potential to calculate
//  the forces on each atom. Then uses a = F/m to calculate the
//  accelleration of each atom.
void computeAccelerations() {
    int i, k;
    double rij[3]; // position of i relative to j

    for (i = 0; i < N; i++) {
        // set all accelerations to zero
        for (k = 0; k < 3; k++) {
            a[i][k] = 0;
        }
    }
#pragma omp parallel for private(j, rSqd, rij, f)
    for (i = 0; i < N - 1; i++) {
        // loop over all distinct pairs i,j
        for (int j = i + 1; j < N; j++) {
            // initialize r^2 to zero
            double rSqd = 0;

            for (k = 0; k < 3; k++) {
                // component-by-componenent position of i relative to j
                rij[k] = r[i][k] - r[j][k];
                // sum of squares of the components
                rSqd += rij[k] * rij[k];
            }

            // From derivative of Lennard-Jones with sigma and epsilon set equal to 1 in natural units!
            const double f = 24 * (2 * pow(rSqd, -7) - pow(rSqd, -4));
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
double VelocityVerlet(const double dt) {
    double psum = 0.;

    // Compute accelerations from forces at current position
    computeAccelerations();
    // Update positions and velocity with current velocity and acceleration
    //printf("  Updated Positions!\n");
#pragma omp parallel for private(j)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 3; j++) {
            r[i][j] += v[i][j] * dt + 0.5 * a[i][j] * dt * dt;
            v[i][j] += 0.5 * a[i][j] * dt;
        }
        //printf("  %i  %6.4e   %6.4e   %6.4e\n",i,r[i][0],r[i][1],r[i][2]);
    }
    // Update accellerations from updated positions
    computeAccelerations();
    // Update velocity with updated acceleration
#pragma omp parallel for private(j)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 3; j++) {
            v[i][j] += 0.5 * a[i][j] * dt;
        }
    }

    // Elastic walls
#pragma omp parallel for reduction(+:psum) private(j)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 3; j++) {
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

    for (i = 0; i < N; i++) {
        for (j = 0; j < 3; j++) {
            // Pull a number from a Gaussian Distribution
            v[i][j] = gaussdist();
        }
    }

    // Vcm = sum_i^N  m*v_i/  sum_i^N  M
    // Compute center-of-mas velocity according to the formula above
    double vCM[3] = {0, 0, 0};

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
    for (i = 0; i < N; i++) {
        for (j = 0; j < 3; j++) {
            v[i][j] -= vCM[j];
        }
    }

    // Now we want to scale the average velocity of the system
    // by a factor which is consistent with our initial temperature, Tinit
    double vSqdSum, lambda;
    vSqdSum = 0.;
    for (i = 0; i < N; i++) {
        for (j = 0; j < 3; j++) {
            vSqdSum += v[i][j] * v[i][j];
        }
    }

    lambda = sqrt(3 * (N - 1) * Tinit / vSqdSum);

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
    if (!available) {
        double rsq, v1, v2;
        do {
            v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);

        double fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        available = true;

        return v2 * fac;
    } else {
        available = false;
        return gset;
    }
}

int write_output(double Pavg, double Tavg, double Z, double gc, double Vol, int N, double timefac, double dt,
                 const string &ofn) {
    constexpr double VolFac = 3.7949992920124995e-29;
    constexpr int NumTime = 2000;
    if (ofstream outfile("simulation_output.txt"); outfile.is_open()) {
        outfile << "\n  AVERAGE TEMPERATURE (K):                 " << fixed << setprecision(5) << Tavg;
        outfile << "\n  AVERAGE PRESSURE  (Pa):                  " << fixed << setprecision(5) << Pavg;
        outfile << "\n  PV/nT (J * mol^-1 K^-1):                 " << fixed << setprecision(5) << gc;
        outfile << "\n  PERCENT ERROR of pV/nT AND GAS CONSTANT: " << fixed << setprecision(5)
                << 100 * fabs(gc - 8.3144598) / 8.3144598;
        outfile << "\n  THE COMPRESSIBILITY (unitless):          " << fixed << setprecision(5) << Z;
        outfile << "\n  TOTAL VOLUME (m^3):                      " << fixed << setprecision(5) << Vol * VolFac;
        outfile << "\n  NUMBER OF PARTICLES (unitless):          " << N;
        outfile << "\n  TOTAL TIME (s):                          " << fixed << setprecision(4) << NumTime * dt * timefac
                << "\n";
        outfile.close();
        cout << "\n  AVERAGE TEMPERATURE (K):                 " << fixed << setprecision(5) << Tavg;
        cout << "\n  AVERAGE PRESSURE  (Pa):                  " << fixed << setprecision(5) << Pavg;
        cout << "\n  PV/nT (J * mol^-1 K^-1):                 " << fixed << setprecision(5) << gc;
        cout << "\n  PERCENT ERROR of pV/nT AND GAS CONSTANT: " << fixed << setprecision(5)
                << 100 * fabs(gc - 8.3144598) / 8.3144598;
        cout << "\n  THE COMPRESSIBILITY (unitless):          " << fixed << setprecision(5) << Z;
        cout << "\n  TOTAL VOLUME (m^3):                      " << fixed << setprecision(5) << Vol * VolFac;
        cout << "\n  NUMBER OF PARTICLES (unitless):          " << N;
        cout << "\n  TOTAL TIME (s):                          " << fixed << setprecision(4) << NumTime * dt * timefac
                << "\n";
        cout << "\n  TO ANALYZE INSTANTANEOUS DATA ABOUT YOUR MOLECULE, OPEN THE FILE \n  '" << ofn
                << "' WITH YOUR FAVORITE TEXT EDITOR OR IMPORT THE DATA INTO EXCEL\n";
    } else {
        cerr << "Error: Could not open output file simulation_output.txt" << endl;
        return 1;
    }
    return 0;
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
        return measure(num_threads, output_file);;
    } catch (const exception &e) {
        cerr << "Error: Invalid input arguments: " << e.what() << endl;
        return 1;
    }
}
