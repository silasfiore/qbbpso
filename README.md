# Quaternion Bare Bones Particle Swarm Optimizer
> A PSO to determine the parameters of torque-free attitude profiles that minimize the mean angular residual of vector measurements
## Features
* Text-based input of measurements and output of all indivdual best positions of the particles
* Analytic attitude propagation for symmetric bodies, fixed-step numerical integration for the general case
* Only the initial attitude and the inertial angular momentum normalized by the intermediate moment of inertia are determined, the inertia constants c1,c2, and c3 are assumed to be known and a prior distribution of the magnitude of the normalized angular momentum vector is used.
## Installation & Build
This project uses a standard `Makefile` for compilation. It does not require any external dependencies outside of a standard C compiler and the system math library (`-lm`).
### Prerequisites
Ensure you have a C99-compatible compiler (like `gcc` or `clang`) and `make` installed on your system:
* **Linux (Ubuntu/Debian):** `sudo apt install build-essential`
* **macOS:** `xcode-select --install`
### Building the Project
1. **Clone the repository:**
   ```bash
   git clone https://github.com/silasfiore/qbbpso.git
   cd qbbpso
   ```
2. **Compile the project using the provided Makefile:**
   ```bash
   make
   ```
## Usage
This project compiles into a command-line utility. It requires an input text file and several arguments and generates an output text file.
### Command Syntax

```bash
./qbbpso <inputfile> <magnitude_l> <sigma_l> <c1> <c2> <c3> [-t <numerical integration time step> -n <number of particles> -i <number of iterations> -o <outputfile>]
```
The normalized angular momentum is defined by the initial angular momentum and can be computed using 
the initial angular velocity in the body frame using

$$\mathbf{l} = \frac{1}{I_2}\text{DCM}(\mathbf{q}_0)\begin{bmatrix} I_1 & 0 & 0\\\\ 0 & I_2 & 0 \\\\ 0 & 0 & I_3 \end{bmatrix} \boldsymbol{\omega}_{b,0}$$ 

where 
$$I_1 \leq I_2 \leq I_3$$ 

Its norm has the unit radians / second and for a freely precessing top, it is the angular rate of the precession. The inertia constants are defined as 

$$c_1 = \tfrac{I_2 - I_3}{I_1},
c_2 = \tfrac{I_3 - I_1}{I_2}, \text{ and }
c_3 = \tfrac{I_1 - I_2}{I_3}$$

The number of particles and number of iterations control how many particles are used to
explore the search space and how many times they are resampled. The input file has the format

```text
t rx ry rz bx by bz
```
and the rows should be sorted in increasing order of the relative time t which is in seconds.
The unit vector b is the direction in the body frame that was aligned with the inertial direction r
at time t.
The output file has the format 
```text
mean_residual qw qx qy qz nx ny nz m
```
where the mean residual is the average angular error of the attitude profile for all the measurements.
The attitude quaternion q0 = [qw, qx, qy,qz] define the attitude at t = 0. And the normalized 
angular momentum l = [nx,ny,nz] * m uniquely defines the attitude as a function of time 
under the assumption of torque-free motion and using the known inertia constants.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
## References & Acknowledgments
This project utilizes several specialized mathematical algorithms, third-party routines, and academic research papers:
### Papers
* **Quaternion Integration:** Implemented via the Crouch-Grossman method to preserve geometric constraints. 
  * *Andrle, M. S., & Crassidis, J. L. (2013). Geometric integration of quaternions. Journal of Guidance, Control, and Dynamics, 36(6), 1762-1767.*
* **Directional Sampling (3D vMF):** Numerically stable sampling on $S^2$.
  * *Jakob, W. (2012). Numerically stable sampling of the von Mises-Fisher distribution on $S^2$ (and other tricks). Interactive Geometry Lab, ETH Zürich, Tech. Rep, 6.*
* **Directional Sampling (4D vMF):** Simulation of the von Mises-Fisher distribution on $S^3$.
  * *Wood, A. T. (1994). Simulation of the von Mises Fisher distribution. Communications in statistics-simulation and computation, 23(1), 157-164.*
* **Uniform Sphere Point Picking:** Generated via the polar rejection method.
  * *Marsaglia, G. (1972). Choosing a Point from the Surface of a Sphere. The Annals of Mathematical Statistics, 43(2), 645-646.*
### Codebases & Implementation Sources
* **3D Math Library:** Base vector/matrix routines adapted and modified from the open-source [cglm](https://github.com/recp/cglm) library (MIT License).
* **Thompson Sampling:** Logic adapted and modified from James McCaffrey's article ["Test Run - Thompson Sampling Using C#"](https://learn.microsoft.com/en-us/archive/msdn-magazine/2018/february/test-run-thompson-sampling-using-csharp) (MSDN Magazine, Feb 2018).