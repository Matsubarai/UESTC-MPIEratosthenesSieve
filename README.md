# UESTC-MPIEratosthenesSieve
Automatically build and test the MPI implementations of Eratosthenes sieve. Distributed Parallel Computing 2021 Spring in UESTC.

# Environment
macOS 11.2.2/CentOS 8.3.2011
OpenMPI 4.1.0
Intel Core i5-8259U (L3 Cache 6MBx1, L2 Cache 256KBx4)

# Usage
```shell
mkdir build && cd build
cmake ..
make
ctest
```
