# Reference
If you use IBM_TENO for academic applications, please cite our paper. (https://doi.org/10.1016/j.ast.2025.111073)

# IBM_TENO

**IBM_TENO** is a GPU-accelerated solver for simulating supersonic flows using a combination of the Immersed Boundary Method (IBM) and Targeted Essentially Non-Oscillatory (TENO) schemes. This solver is implemented in Fortran and is designed to handle complex geometries efficiently on CUDA-enabled GPUs.

## Features

- GPU parallelized using NVIDIA CUDA Fortran
- Combines Immersed Boundary Method with high-order TENO schemes
- Capable of solving supersonic flow over complex geometries

## Included Test Cases

The repository contains three `.f90` source files, each corresponding to a specific benchmark case:

- **`Cylinder.f90`**: Supersonic flow around a single cylinder at Mach 3
- **`Multi_Cylinder.f90`**: Supersonic flow around multiple cylinders at Mach 2
- **`Schardin.f90`**: Simulation of Schardin's problem (shock interaction with a wedge)

## Compilation

The solver is written in CUDA Fortran and can be compiled using the `nvfortran` compiler from the NVIDIA HPC SDK.

### Compile Command:

```bash
nvfortran -cuda -O3 <source_file>.f90 -o <output_binary>.out
```

For example:

```bash
nvfortran -cuda -O3 Cylinder.f90 -o Cylinder.out
```

# Contributing
We welcome any contributions and feedback that can help improve MHD_GDDC_GPU. If you would like to contribute to the tool, please contact the maintainers or open an Issue in the repository or a thread in Discussions. Pull Requests are encouraged, but please propose or discuss the changes in the associated Issue beforehand.

# Licencing
Please refer to the licence file.
