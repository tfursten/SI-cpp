SI-cpp
======
Self-Incompatibility Simulation (c++ version)

A spatially explicit individual-based model of a diploid plant population that reproduces according to five different models of self-incompatibility.

Author
------
Tara Furstenau  
Biodesign Institute  
Center for Evolutionary Medicine and Informatics  
Arizona State University  
[Website](http://tfursten.github.io)  

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
5. Compile  
  ```
  make
  ```

Dependencies
-------------
The Boost c++ Library is required for compilation and usage.
* Foreach
* Program Options

Run
----
Usage:
```
./si config.txt
./si --help
Allowed Options:

General Options:
  --help                Produce help message

Configuration:
  -x [ --maxX ] arg (=100)              Set X dimension
  -y [ --maxY ] arg (=100)              Set Y dimension
  --landscape arg (=torus)              Set landscape boundary (torus or 
                                        rectangle)
  -g [ --generations ] arg (=10)        Set number of Generations to run after 
                                        burn-in
  -p [ --pollen ] arg (=10)             Set number of pollen produced per 
                                        individual
  -o [ --ovule ] arg (=10)              Set number of ovules per individual
  -n [ --markers ] arg (=3)             Set number of markers
  -u [ --smut ] arg (=1.0000000000000001e-05)
                                        Set S locus mutation rate
  -m [ --mmut ] arg (=1.0000000000000001e-05)
                                        Set marker mutation rate
  --pdel arg (=1)                       Set deleterious selection coefficient
  --dmut arg (=0.0001)                  Set deleterious mutation rate for 
                                        unlinked locus
  -d [ --distribution ] arg (=disk)     Set Dispersal Distribution
  -q [ --sigmaP ] arg (=2)              Set dispersal parameter for pollen
  -r [ --sigmaS ] arg (=2)              Set dispersal parameter for seed
  -b [ --burn ] arg (=0)                Set Burn-in Period
  -t [ --sample ] arg (=1)              Sample every n generations after 
                                        burn-in
  -f [ --output_file ] arg (=data)      Output File Name
  --seed arg (=0)                       Set PRNG seed
  -s [ --si ] arg (=nsi)                Set self-incompatibility system
  --pparam arg (=0)                     Extra Parameter for pollen dispersal
  --sparam arg (=0)                     Extra Parameter for seed dispersal
  --fast arg (=1)                       Use fast dispersal when available





```
