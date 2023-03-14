# ESS-color


ESS-Color is a bioinformatics tool for constructing compressed representation of sets of k-mer sets.


## Requirements

- Linux operating system (64 bit)
- GCC >= 4.8 or a C++11 capable compiler
- Snakemake
- Git
- CMake 3.12+
- [ESSCompress](https://github.com/medvedevgroup/ESSCompress)
- [joinCounts](https://github.com/Transipedia/dekupl-joinCounts)
- [KMC](https://github.com/refresh-bio/KMC)


## Quick start

First, install all the pre-requisites and make sure the executables are in your `PATH`. Then, install additional executables from source:

    git clone https://github.com/medvedevgroup/ESSColor.git
    cd ESSColor
    bash compile.sh
    
You can move/copy ALL the executables in `ess-color/bin` to the bin directory that is already in your PATH. For instance, considering `/usr/bin` is already in PATH, you need to run the command `mv ess-color/bin/* /usr/bin` to move all executables for ESS-Color software. An alternative to moving/copying executables is adding the location of `ess-color/bin` to your PATH.


After compiling, set up the `config.yaml`. Change to the directory and use [TODO: add details of config file]

    snakemake -s Snakefile.smk -j 8
    



## Command lines

### Color matrix generation

Generates a color matrix from a KMC database list.

    genmatrix [OPTION...]
    -c, --count-list arg  [Mandatory] Path to KMC database files. One line per database
    -d, --debug-verif     Debug flag to verify if the output coresponds to the input (Time consuming).
    -o, --outmatrix arg   [Mandatory] Path to the output color matrix
    -s, --strout          String output (default is binary)

