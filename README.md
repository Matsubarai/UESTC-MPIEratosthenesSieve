# Eratosthenes Sieve Based on MPI
Automatically build and test the MPI implementations of Eratosthenes Sieve. 

For Distributed Parallel Computing 2021 Spring in UESTC.

# Environment
OS: macOS 11.2.2/CentOS 8.3.2011

CPU: Intel Core i5-8259U (L3 Cache 6MBx1, L2 Cache 256KBx4)

# Dependency
OpenMPI 4.0.3+

CMake 3.11+

C/C++ Compiler

# Usage
```shell
mkdir build && cd build
cmake ..
make
ctest
```
