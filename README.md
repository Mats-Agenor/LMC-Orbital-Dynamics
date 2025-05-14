# LMC-Orbital-Dynamics
------------
Compile with: mpic++ -O3 -std=c++17 LMC.cpp -o LMC
Run with: mpirun -np 20 ./LMC

On a cluster, using 20 MPI process, the simulation takes approximately 5 minutes to complete.

------------

This C++ code simulates the orbit of the Large Magellanic Cloud (LMC) around the Milky Way on a timescale of a few billion years. It is designed to be executed in parallel using MPI, which allows the code to distribute different combinations of parameters among multiple threads. The main goal is to study the dynamical evolution of the LMC, considering realistic gravitational potentials and dynamic friction.

### Step-by-step explanation:

1. **Parallelization Configuration (MPI):**
The code uses MPI to distribute simulations among different processors. Each process should explore different initial conditions, varying parameters such as position and velocity within a specified grid.

2. **Temporal Integration Settings:**
The total integration time is set to 5 Gyr (Gigayears) and the time step is set to 0.01 Gyr (10 Myr). A Leapfrog integrator is used because of its symplectic nature, which conserves energy reasonably well over long time scales. The integrator is implemented adaptively to increase accuracy.

3. **Potential Models Used:**
- **Disk Potential:**
An exponential disk potential is used to represent the stellar disk of the Milky Way. This potential provides a flattened gravitational field, suitable for modeling disk galaxies.
- **Dark Matter Halo Potential:**
A Navarro-Frenk-White (NFW) profile is used to describe the dark matter halo of the Milky Way. The NFW profile is widely used in cosmological simulations and corresponds to the density distribution of halos formed in cold dark matter models.
- **Bulge Potential:**
A Hernquist potential is employed to represent the Large Magellanic Cloud. The Hernquist profile closely mimics the observed properties of spheroidal components in galaxies, in particular, bulges and halos.

4. **Dynamic Friction:**
The code includes a dynamic friction term based on the Chandrasekhar formula. This force is responsible for the gravitational drag experienced by the LMC as it moves through the Milky Way's dark matter halo, which slows its motion over time.

5. **Energy and Distance Tracking:**
Throughout the simulation, the code saves the energy and galactocentric distance of the LMC at each time step. This information is written to individual output files for later analysis and plotting.

### Final Comments:

For each mpi rank, a file is generated containing the galactocentric distance of the LMC as a function of time and the total energy of the system (kinetic + potential) over time. My idea was to obtain the orbital energy loss and compare it with observations and N-body simulations. Unfortunately, all the outputs are constant, which physically should not happen. Despite many attempts, I was unable to fix the problem. If anyone accessing this repository knows how to fix the code or wants to try, feel free to do so!
