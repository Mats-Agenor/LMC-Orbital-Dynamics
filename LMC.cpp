#include <iostream>
#include <vector>
#include <cmath>
#include <mpi.h>
#include <fstream>
#include <iomanip>
#include <stdexcept>

// Constants
const double G = 4.302e-3;       // Gravitational constant (kpc/Msun)(km/s)^2
const double M_MW = 1.5e12;      // Milky Way mass (Msun)
const double M_LMC = 1e11;       // LMC mass (Msun)
const double a_LMC = 20.0;       // LMC scale radius (kpc)

// NFW & Disk parameters (updated to MW-like values)
const double rho0_NFW = 0.012e9; // NFW scale density (Msun/kpc^3)
const double rs_NFW = 16.0;      // NFW scale radius (kpc)
const double Sigma0_disk = 0.8e9; // Disk central surface density (Msun/kpc^2)
const double Rd_disk = 2.5;      // Disk scale radius (kpc)

// Safe Bessel functions (with asymptotic expansions)
inline double bessel_i0_safe(double x) {
    if (x > 700.0) return exp(x) / sqrt(2 * M_PI * x);
    return std::cyl_bessel_i(0, x);
}

inline double bessel_k0_safe(double x) {
    if (x > 700.0) return exp(-x) * sqrt(M_PI / (2 * x));
    return std::cyl_bessel_k(0, x);
}

inline double bessel_i1_safe(double x) {
    if (x > 700.0) return exp(x) / sqrt(2 * M_PI * x);
    return std::cyl_bessel_i(1, x);
}

inline double bessel_k1_safe(double x) {
    if (x > 700.0) return exp(-x) * sqrt(M_PI / (2 * x));
    return std::cyl_bessel_k(1, x);
}

// Potential functions (with stability checks)
double potential_NFW(double r) {
    if (r <= 1e-6) return -4 * M_PI * G * rho0_NFW * rs_NFW * rs_NFW;
    return -4 * M_PI * G * rho0_NFW * pow(rs_NFW, 3) * log(1.0 + r/rs_NFW) / r;
}

double potential_disk(double R) {
    double x = R / (2.0 * Rd_disk);
    double I0 = bessel_i0_safe(x);
    double K0 = bessel_k0_safe(x);
    return -M_PI * G * Sigma0_disk * R * (I0 * K0);
}

double potential_Hernquist(double r) {
    return -G * M_LMC / (r + a_LMC);
}

// NFW density profile
double density_NFW(double r) {
    return rho0_NFW / ((r/rs_NFW) * pow(1.0 + r/rs_NFW, 2));
}

// Disk density profile
double density_disk(double R) {
    return Sigma0_disk * exp(-R/Rd_disk) / (2.0 * Rd_disk);
}

// Compute acceleration (with safeguards)
void compute_acceleration(const double r[3], const double v[3], double a[3], double lnLambda) {
    double r_mag = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    double R_cyl = sqrt(r[0]*r[0] + r[1]*r[1]);
    double v_mag = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);

    // NFW gradient (regularized)
    double term = r_mag / rs_NFW;
    double dPhi_NFW_dr = 4 * M_PI * G * rho0_NFW * pow(rs_NFW, 3) * 
                        (log(1.0 + term) - term/(1.0 + term)) / (r_mag * r_mag * r_mag);

    // Disk gradient (softened)
    double x = R_cyl / (2.0 * Rd_disk);
    double I0 = bessel_i0_safe(x);
    double K0 = bessel_k0_safe(x);
    double I1 = bessel_i1_safe(x);
    double K1 = bessel_k1_safe(x);
    double dPhi_disk_dR = -M_PI * G * Sigma0_disk * (I0*K0 - I1*K1);

    // Hernquist gradient
    double dPhi_Hernquist_dr = G * M_LMC / pow(r_mag + a_LMC, 2);

    // Total gravitational acceleration
    for (int i = 0; i < 3; i++) {
        a[i] = -dPhi_NFW_dr * r[i];
        if (i < 2) a[i] += -dPhi_disk_dR * r[i]/std::max(R_cyl, 1e-6); // Softening
        a[i] += -dPhi_Hernquist_dr * r[i]/r_mag;
    }

    // Dynamical friction (with velocity floor)
    double rho_disk = density_disk(R_cyl);
    double rho_total = density_NFW(r_mag) + rho_disk;
    double v_mag_safe = std::max(v_mag, 10.0); // Avoid v→0 divergence
    double F_DF_mag = 4 * M_PI * pow(G * M_LMC, 2) * rho_total * lnLambda / pow(v_mag_safe, 3);

    for (int i = 0; i < 3; i++) {
        a[i] += -F_DF_mag * v[i]/(v_mag_safe * M_LMC);
    }
}

// Adaptive Leapfrog integrator
void leapfrog(double r[3], double v[3], double a[3], double& dt, double lnLambda, int step, int rank) {
    const double max_dt = 1e6; // 1 Myr (maximum allowed dt)
    const double min_dt = 1e3; // 1 kyr (minimum allowed dt)
    double r_old[3] = {r[0], r[1], r[2]};
    double v_old[3] = {v[0], v[1], v[2]};

    // First half-step velocity update
    for (int i = 0; i < 3; i++) {
        v[i] += 0.5 * dt * a[i];
    }

    // Full position step
    for (int i = 0; i < 3; i++) {
        r[i] += dt * v[i];
    }

    // Calculate new acceleration with updated position
    compute_acceleration(r, v, a, lnLambda);

    // Second half-step velocity update
    for (int i = 0; i < 3; i++) {
        v[i] += 0.5 * dt * a[i];
    }

    // Check for instability
    double r_mag = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    if (r_mag > 1000.0 || std::isnan(r_mag)) {
        // Roll back and reduce dt
        for (int i = 0; i < 3; i++) {
            r[i] = r_old[i];
            v[i] = v_old[i];
        }
        dt = std::max(dt * 0.5, min_dt);
        if (rank == 0 && step % 100 == 0) {
            std::cerr << "WARNING: Reducing dt to " << dt << " years\n";
        }
    } else if (step % 100 == 0 && dt < max_dt) {
        dt = std::min(dt * 1.1, max_dt); // Slowly increase dt
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        std::cout << "\n=== LMC Orbit Simulation ===\n";
        std::cout << "Running with " << size << " MPI processes\n";
        std::cout << "Simulation time: 5 Gyr\n";
    }

    // Initial conditions (realistic LMC values)
    double r[3] = {-0.5, -41.2, -27.2};  // Initial position (kpc)
    double v[3] = {-86, -268, 252};       // Initial velocity (km/s)
    double lnLambda = 3.0 + 0.1 * rank;   // lnΛ ∈ [3.0, 4.9] for 20 processes

    // Time parameters (5 Gyr simulation)
    double dt = 1e4;                      // Initial time step (10 kyr) in years
    double t_end = 5e9;                   // Total time (5 Gyr) in years
    int output_interval = 1000;           // Output every 1000 steps

    // Output file
    std::string filename = "lmc_orbit_rank_" + std::to_string(rank) + ".txt";
    std::ofstream outfile(filename);
    outfile << std::scientific << std::setprecision(6);
    outfile << "# Time(Myr) X(kpc) Y(kpc) Z(kpc) Vx(km/s) Vy(km/s) Vz(km/s) Distance(kpc) Energy(km^2/s^2)\n";

    // Initial energy calculation
    double r_mag = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    double R_cyl = sqrt(r[0]*r[0] + r[1]*r[1]);
    double energy = 0.5 * (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) + 
                   potential_NFW(r_mag) + 
                   potential_disk(R_cyl) + 
                   potential_Hernquist(r_mag);
    outfile << 0.0 << " " << r[0] << " " << r[1] << " " << r[2] << " "
            << v[0] << " " << v[1] << " " << v[2] << " " << r_mag << " " << energy << "\n";

    if (rank == 0) {
        std::cout << "\nInitial conditions:\n";
        std::cout << "Position: [" << r[0] << ", " << r[1] << ", " << r[2] << "] kpc\n";
        std::cout << "Velocity: [" << v[0] << ", " << v[1] << ", " << v[2] << "] km/s\n";
        std::cout << "Initial energy: " << energy << " km²/s²\n";
        std::cout << "\nStarting simulation...\n";
    }

    // Calculate initial acceleration
    double a[3];
    compute_acceleration(r, v, a, lnLambda);

    // Main integration loop
    double current_time = 0.0;
    int step = 0;
    while (current_time < t_end) {
        leapfrog(r, v, a, dt, lnLambda, step, rank);
        current_time += dt;
        step++;

        // Calculate energy and distance
        r_mag = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        R_cyl = sqrt(r[0]*r[0] + r[1]*r[1]);
        energy = 0.5 * (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) + 
                potential_NFW(r_mag) + 
                potential_disk(R_cyl) + 
                potential_Hernquist(r_mag);

        // Write output at regular intervals
        if (step % output_interval == 0) {
            outfile << current_time / 1e6 << " " << r[0] << " " << r[1] << " " << r[2] << " "
                   << v[0] << " " << v[1] << " " << v[2] << " " << r_mag << " " << energy << "\n";
            outfile.flush();
        }

        // Early termination if escaping
        if (r_mag > 1000.0) {
            if (rank == 0) std::cout << "LMC escaped beyond 1000 kpc at time " << current_time/1e6 << " Myr\n";
            break;
        }
    }

    outfile.close();
    if (rank == 0) {
        std::cout << "\nSimulation complete. Results saved to lmc_orbit_rank_*.txt\n";
        std::cout << "Final time reached: " << current_time/1e6 << " Myr\n";
        std::cout << "Final position: [" << r[0] << ", " << r[1] << ", " << r[2] << "] kpc\n";
        std::cout << "Final velocity: [" << v[0] << ", " << v[1] << ", " << v[2] << "] km/s\n";
        std::cout << "Final energy: " << energy << " km²/s²\n";
    }

    MPI_Finalize();
    return 0;
}
