SI-cpp
======
Self-Incompatibility Simulation (c++ version)

A spatially explicit individual-based model of diploid plant population that reproduce according to five different models of self-incompatibility.

Author
------
Tara Furstenau  
Biodesign Institute  
Center for Evolutionary Medicine and Informatics  
Arizona State University  
[Website](http://tfursten@github.io)  

Compiling from Source Code
--------------------------
SI-cpp requires [CMake 2.8](http://www.cmake.org/) to build from source. 

1. Download the source code.  
2. Decompress the tar-bzip archive  
  ```
  tar xvzf SI-cpp-*.tar.bz2
  ```
3. Change to the build directory.  
  ```
  cd SI-cpp-*/build
  ```
4. Run the CMake build system.  
  ```
  cmake ..
  ```
5.Compile  
  ```
  make
  ```


Dependencies
-------------
The Boost c++ Library is required for compilation and usage.

Run
----
Usage:
```
./si config.txt
./si --help
```
