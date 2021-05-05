# Eratosthenes Sieve Based on MPI
[![](https://img.shields.io/github/workflow/status/Matsubarai/UESTC-MPIEratosthenesSieve/CMake)](https://github.com/Matsubarai/UESTC-MPIEratosthenesSieve/actions/workflows/cmake.yml)
![](https://img.shields.io/badge/platform-linux%20%7C%20macos-lightgrey)

Automatically build and test the MPI implementations of Eratosthenes Sieve. 

For Distributed Parallel Computing Spring 2021 in UESTC.

# Environment
- OS: macOS 10.15+/CentOS 8

- CPU: Intel Core i5-8259U (L3 Cache 6MBx1, L2 Cache 256KBx4)

# Dependency
- OpenMPI 4.0.3+

- CMake 3.11+

- Clang/GCC

# Usage
```shell
mkdir build && cd build
cmake ..
make
ctest
```
