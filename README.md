# LMC-Orbital-Dynamics
------------
Compile with: mpic++ -O3 -std=c++17 LMC.cpp -o LMC  
Run with: mpirun -np 20 ./LMC  

On a cluster, using 20 MPI processes, the simulation takes approximately 5 minutes to complete.

------------

This C++ code simulates the orbit of the Large Magellanic Cloud (LMC) around the Milky Way galaxy over a timescale of several billion years. It is designed to run in parallel using MPI, which allows the code to distribute different parameter combinations across multiple processes. The main goal is to study the dynamical evolution of the LMC while considering realistic gravitational potentials and dynamical friction.

### Step-by-step Explanation:

1. **Parallelization Setup (MPI):**  
   The code uses MPI to distribute simulations across different processors. Each process explores different initial conditions by varying parameters such as position and velocity within a specified grid.

2. **Time Integration Settings:**  
   The total integration time is set to 5 Gyr (Gigayears), and the timestep is set to 0.01 Gyr (10 Myr). A Leapfrog integrator is used due to its symplectic nature, which conserves energy reasonably well over long timescales. The integrator is implemented in an adaptive form to increase accuracy.

3. **Potential Models Used:**
   - **Disk Potential:**  
     A Miyamoto-Nagai disk potential is used to represent the stellar disk of the Milky Way. This potential provides a flattened disk-like gravitational field, suitable for modeling disk galaxies.
   - **Dark Matter Halo Potential:**  
     A Navarro-Frenk-White (NFW) profile is used to describe the Milky Way's dark matter halo. The NFW profile is widely used in cosmological simulations and matches the density distribution of halos formed in cold dark matter models.
   - **Bulge Potential:**  
     A Hernquist potential is employed to represent the central bulge of the Milky Way. The Hernquist profile closely mimics the observed properties of spheroidal components in galaxies.

4. **Dynamical Friction:**  
   The code includes a dynamical friction term based on Chandrasekhar's formula. This force accounts for the gravitational drag experienced by the LMC as it moves through the Milky Way's dark matter halo, which slows down its motion over time.

5. **Initial Conditions and Parameter Scan:**  
   Each MPI process is responsible for a different combination of initial velocities and positions, allowing the user to scan a grid of initial conditions to find a suitable orbit for the LMC.

6. **Energy and Distance Tracking:**  
   Throughout the simulation, the code saves the energy and galactocentric distance of the LMC at each timestep. These are written to individual output files for further analysis and plotting.

7. **Output:**  
   For each simulation, two files are generated:
   - One containing the galactocentric distance of the LMC as a function of time.
   - Another containing the total energy of the system (kinetic + potential) over time.

### Final Comments:

Unfortunately, all outputs are constant, which physically should not happen. Despite many attempts, I was unable to correct the issue. If anyone accessing this repository knows how to fix the code or wants to give it a try, feel free to do so!
